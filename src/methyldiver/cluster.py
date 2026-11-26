"""Select, fit and apply frozen regional methylation-pattern medoids."""
from collections import Counter
from contextlib import ExitStack
import csv
from dataclasses import dataclass
import hashlib
import json
from math import fsum
from pathlib import Path
from tempfile import TemporaryDirectory

import pysam

from . import __version__, full_read
from .coordinates import cpg_positions
from .modbam import (Region, _index_path, _validate_parameters, decode_calls, eligible_reads,
                     modification_config, read_regions_bed, require_modification)

SCHEMA = "methyldiver.cluster/2"
ALGORITHM = "weighted_pam_normalized_hamming"
LEGACY = ("methyldiver.cluster/1", "weighted_pam_two_normalized_hamming")
SETTINGS = dict(modification="5mC", call_threshold=0.8, min_mapq=20, min_base_quality=10,
                min_complete_reads=5, min_cluster_reads=2, max_reads_per_region=100000)
STATUSES = {"fitted", "no_coverage", "fewer_than_three_cpgs", "insufficient_complete_reads",
            "one_distinct_pattern", "insufficient_distinct_patterns", "insufficient_cluster_support"}
REASONS = ("all_missing", "partial", "no_model", "tie")


def add_cluster_command(commands):
    group = commands.add_parser("cluster", help="Select cluster counts, fit regional medoids or assign reads.")
    modes = group.add_subparsers(dest="mode", required=True, title="Commands")
    descriptions = {"select": "Compare candidate cluster counts using training reads.",
                    "fit": "Learn medoids independently at each BED locus.",
                    "apply": "Assign reads to frozen regional medoids."}
    for mode, description in descriptions.items():
        parser = modes.add_parser(mode, allow_abbrev=False, help=description)
        parser.add_argument("--input", type=Path, required=True, help="Indexed coordinate-sorted modBAM.")
        parser.add_argument("--reference", type=Path, required=True, help="Matching indexed FASTA.")
        parser.add_argument("--output", type=Path, required=True, help="New output directory.")
        parser.add_argument("--sample-label", help="Sample identifier [input filename without .bam].")
        if mode == "apply":
            parser.add_argument("--model", type=Path, required=True, help="model.json from cluster fit.")
            continue
        parser.add_argument("--full-read", action="store_true",
                            help="Use the whole aligned molecule, including flanking CpGs [off].")
        parser.add_argument("--min-shared-cpgs", type=int, metavar="N",
                            help="Minimum jointly confident reference CpGs per full-read comparison [3].")
        if mode == "fit":
            parser.add_argument("--clusters", default="2", metavar="K|auto",
                                help="Positive cluster count or auto for per-locus selection [2].")
        parser.add_argument("--max-clusters", type=int, metavar="K",
                            help="Largest candidate for selection, including the one-cluster baseline [6].")
        parser.add_argument("--regions-bed", type=Path, required=True, help="Nonoverlapping BED3/BED4 loci.")
        parser.add_argument("--modification", choices=("5mC", "5hmC"), help="Selected CpG modification [5mC].")
        for name, help_text in (
            ("call_threshold", "Binary-call confidence (0.5, 1]"),
            ("min_mapq", "Minimum mapping quality (0–254)"),
            ("min_base_quality", "Minimum base quality (0–254)"),
            ("min_complete_reads", "Minimum training reads per locus (complete in regional mode, usable in full-read mode)"),
            ("min_cluster_reads", "Minimum uniquely assigned training molecules per medoid"),
            ("max_reads_per_region", "Maximum alignments fetched per locus")):
            flags = ["--" + name.replace("_", "-")]
            if name == "min_complete_reads":
                flags.append("--min-training-reads")
            parser.add_argument(*flags, type=type(SETTINGS[name]),
                                help=help_text+" [%(default)s].")
        parser.set_defaults(**SETTINGS)


def fit_medoids(counts, k=2):
    """Deterministic abundance-weighted PAM on observed binary strings."""
    patterns = sorted(counts)
    if not 1 <= k <= len(patterns):
        raise ValueError("Cluster count must be between one and the number of distinct patterns.")
    values = {p: int(p, 2) for p in patterns}
    first = min(patterns, key=lambda p: (-counts[p], p))
    centers = [first]
    nearest = {p: (v ^ values[first]).bit_count() for p, v in values.items()}
    while len(centers) < k:
        seed = min(patterns, key=lambda p: (-nearest[p], -counts[p], p))
        centers.append(seed)
        for p, v in values.items():
            nearest[p] = min(nearest[p], (v ^ values[seed]).bit_count())
    current = tuple(sorted(centers))
    score = sum(counts[p] * nearest[p] for p in patterns)
    weights = [counts[p] for p in patterns]
    while score:
        best_score, best = score, current
        for removed in current:
            kept = [m for m in current if m != removed]
            # Cost of keeping all but the removed center, reused across candidates.
            remaining = [min(((v ^ values[m]).bit_count() for m in kept), default=len(first)+1)
                         for v in values.values()]
            for candidate in patterns:
                if candidate in current:
                    continue
                x = values[candidate]
                cost = sum(n * min(d, (v ^ x).bit_count())
                           for v, n, d in zip(values.values(), weights, remaining))
                if cost <= best_score and cost < score:
                    trial = tuple(sorted([*kept, candidate]))
                    if (cost, trial) < (best_score, best):
                        best_score, best = cost, trial
        if best_score == score:
            break
        score, current = best_score, best
    return list(current), score


def assign(pattern, medoids):
    empty = (None,) * (len(medoids) + 1)
    if not pattern or set(pattern) == {"."}:
        return "unassigned", "all_missing", *empty
    if "." in pattern:
        return "unassigned", "partial", *empty
    if not medoids:
        return "unassigned", "no_model", *empty
    x = int(pattern, 2)
    mismatches = [(x ^ int(m, 2)).bit_count() for m in medoids]
    distances = [v / len(pattern) for v in mismatches]
    ordered = sorted(mismatches)
    margin = (ordered[1] - ordered[0]) / len(pattern) if len(medoids) > 1 else None
    if margin == 0:
        return "unassigned", "tie", *distances, margin
    return f"C{mismatches.index(ordered[0])+1}", "", *distances, margin


@dataclass
class Candidate:
    k: int
    status: str
    medoids: tuple = ()
    cost: float | None = None
    support: tuple = ()
    ties: int | None = None
    score: float | None = None


def silhouette(groups, total):
    """Mean over complete molecules, assigning zero to ties and singletons."""
    # Summing Hamming distances requires only the per-CpG counts in each group.
    n_sites = len(next(iter(groups[0])))
    sizes = [sum(g.values()) for g in groups]
    ones = [[sum(n for p, n in g.items() if p[j] == "1") for j in range(n_sites)] for g in groups]
    contributions = []
    for own, group in enumerate(groups):
        if sizes[own] == 1:
            continue
        for pattern, count in sorted(group.items()):
            sums = [sum(size-v if bit == "1" else v for bit, v in zip(pattern, column))
                    for size, column in zip(sizes, ones)]
            a = sums[own] / (sizes[own]-1)
            b = min(s / size for i, (s, size) in enumerate(zip(sums, sizes)) if i != own)
            contributions.append(count * (b-a) / max(a,b) if max(a,b) else 0.0)
    return fsum(contributions) / total


def evaluate(counts, n_reads, n_cpgs, settings, k, *, score=False):
    complete = sum(counts.values())
    if n_cpgs < 3:
        return Candidate(k, "fewer_than_three_cpgs")
    if not n_reads:
        return Candidate(k, "no_coverage")
    if complete < settings["min_complete_reads"]:
        return Candidate(k, "insufficient_complete_reads")
    if len(counts) < k:
        return Candidate(k, "one_distinct_pattern" if len(counts) == 1 else "insufficient_distinct_patterns")
    medoids, cost = fit_medoids(counts, k)
    groups = [Counter() for _ in medoids]
    ties = 0
    for pattern, n in counts.items():
        cluster, reason, *_ = assign(pattern, medoids)
        if reason:
            ties += n
        else:
            groups[int(cluster[1:])-1][pattern] = n
    support = tuple(sum(g.values()) for g in groups)
    status = "fitted" if min(support) >= settings["min_cluster_reads"] else "insufficient_cluster_support"
    value = None
    if score and status == "fitted" and k > 1:
        value = silhouette(groups, complete)
    return Candidate(k, status, tuple(medoids), cost, support, ties, value)


def select_count(candidates):
    valid = (c for c in candidates if c.status == "fitted")
    # A single group is the no-split baseline. A split must improve on zero.
    return max(valid, key=lambda c: (round(c.score or 0, 12), -c.k), default=None)


def fit_region(counts, n_reads, n_cpgs, settings, k=2):
    result = evaluate(counts, n_reads, n_cpgs, settings, k)
    return (result.status, list(result.medoids), result.cost) if result.status == "fitted" else (result.status, [], None)


def validate_settings(s):
    if not isinstance(s, dict) or s.keys() != SETTINGS.keys():
        raise ValueError("Invalid clustering model settings.")
    modification_config(s["modification"])
    _validate_parameters(3, s["call_threshold"], s["min_mapq"], s["min_base_quality"],
                         s["min_complete_reads"], 1, s["max_reads_per_region"], 1)
    if type(s["min_cluster_reads"]) is not int or s["min_cluster_reads"] < 1:
        raise ValueError("min_cluster_reads must be a positive integer.")


def check_regions(regions):
    if not regions or len({r.region_id for r in regions}) != len(regions):
        raise ValueError("Supply nonempty regions with unique region IDs.")
    ordered = sorted(regions, key=lambda r: (r.chrom, r.start))
    if any(a.chrom == b.chrom and a.end > b.start for a, b in zip(ordered, ordered[1:])):
        raise ValueError("Clustering regions overlap.")


def load_model(path):
    model = json.loads(path.read_text())
    try:
        identity = model["schema"], model["algorithm"]
        full = identity == full_read.IDENTITY
        expected_assignment = full_read.ASSIGNMENT if full else "complete_only_untied"
        if (identity not in (LEGACY, (SCHEMA, ALGORITHM), full_read.IDENTITY) or
                model["assignment"] != expected_assignment or not isinstance(model["regions"], list)):
            raise ValueError("Unsupported clustering model schema, algorithm or assignment mode.")
        if full and (type(model["min_shared_cpgs"]) is not int or model["min_shared_cpgs"] < 1):
            raise ValueError("Invalid full-read minimum shared CpG count.")
        if identity == LEGACY:
            model["clusters"] = model["max_clusters"] = 2
        k, limit = model["clusters"], model["max_clusters"]
        if type(limit) is not int or limit < 1 or (k != "auto" and (type(k) is not int or k != limit)):
            raise ValueError("Invalid cluster count in model.")
        validate_settings(model["settings"])
        regions = []
        for record in model["regions"]:
            region = Region(record["chrom"], record["start"], record["end"], record["region_id"])
            if full:
                if (not region.name or type(record["contig_length"]) is not int
                        or record["contig_length"] < region.end or record["fit_status"] not in full_read.STATUSES
                        or not isinstance(record["medoids"], list)):
                    raise ValueError("Invalid full-read region in model.")
                profiles = full_read.load_profiles(record, model["min_shared_cpgs"])
                if record["fit_status"] == "fitted":
                    if not 1 <= len(profiles) <= limit or (k != "auto" and len(profiles) != k):
                        raise ValueError("Invalid full-read medoid count.")
                elif profiles:
                    raise ValueError("Failed region unexpectedly contains medoids.")
                regions.append(region)
                continue
            anchors, medoids = record["cpg_anchors"], record["medoids"]
            if (not region.name or not isinstance(anchors, list) or
                    any(type(p) is not int or not region.start <= p < region.end-1 for p in anchors) or
                    anchors != sorted(set(anchors)) or record["fit_status"] not in STATUSES or
                    not isinstance(medoids, list)):
                raise ValueError("Invalid region or CpG anchors in model.")
            if record["fit_status"] == "fitted":
                if (len(anchors) < 3 or not 1 <= len(medoids) <= limit or
                        (k != "auto" and len(medoids) != k) or medoids != sorted(set(medoids)) or
                        any(not isinstance(m, str) or len(m) != len(anchors) or set(m)-set("01") for m in medoids)):
                    raise ValueError("Invalid fitted medoids in model.")
            elif medoids:
                raise ValueError("Failed region unexpectedly contains medoids.")
            regions.append(region)
        check_regions(regions)
    except (KeyError, TypeError, AttributeError) as exc:
        raise ValueError(f"Malformed clustering model: {exc}") from exc
    return model, regions


def tsv(stack, directory, name, columns):
    handle = stack.enter_context((directory/name).open("w", newline="", encoding="utf-8"))
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(columns.split())
    return writer


def analyze_regions(args, model, regions, fasta, bam, dest, stack, label):
    fitting = args.mode == "fit"
    training = args.mode != "apply"
    settings = model["settings"]
    requested, limit = model["clusters"], model["max_clusters"]
    selecting = training and requested == "auto"
    if selecting:
        selection = tsv(stack, dest, "cluster_selection.tsv",
            "sample region_id k eligible_reads complete_reads distinct_patterns status "
            "cluster_support tied_reads mismatch_cost adjusted_silhouette recommended")
    if args.mode != "select":
        cluster_names = [f"C{i}" for i in range(1, limit+1)]
        reads_out = tsv(stack, dest, "assignments.tsv",
            "sample region_id read_group read_id calls observed_cpgs model_cpgs cluster reason "
            + " ".join("distance_"+c for c in cluster_names) + " margin")
        summary = tsv(stack, dest, "regions.tsv",
            "sample region_id fit_status model_cpgs fetched excluded_flags excluded_mapq eligible "
            "complete partial all_missing distinct_complete_patterns used_to_fit "
            + " ".join("assigned_"+c for c in cluster_names) + " unassigned ties no_model")
        patterns = tsv(stack, dest, "patterns.tsv", "sample region_id pattern complete_reads")
    total_qc = Counter()
    for i, r in enumerate(regions):
        if r.chrom not in bam.references or r.chrom not in fasta.references:
            raise ValueError(f"{r.region_id}: contig missing from BAM or FASTA.")
        length = fasta.get_reference_length(r.chrom)
        if length != bam.get_reference_length(r.chrom) or r.end > length:
            raise ValueError(f"{r.region_id}: incompatible reference length or boundaries.")
        seq = fasta.fetch(r.chrom, r.start, r.end).upper()
        anchors = cpg_positions(seq, r.start)
        digest = hashlib.sha256(seq.encode("ascii")).hexdigest()
        if not training:
            record = model["regions"][i]
            if anchors != record["cpg_anchors"] or digest != record.get("reference_sha256") or length != record.get("contig_length"):
                raise ValueError(f"{r.region_id}: model/reference mismatch in anchors, sequence or contig length.")
        anchor_set = set(anchors)
        qc, counts = Counter(), Counter()
        rows = []
        for read in eligible_reads(bam, r, settings, qc):
            calls, _ = decode_calls(read, anchor_set, settings["call_threshold"],
                                    settings["min_base_quality"], qc, settings["modification"])
            pattern = "".join("." if calls.get(a) is None else str(calls[a]) for a in anchors)
            if args.mode != "select":
                rows.append((read.get_tag("RG") if read.has_tag("RG") else "", read.query_name, pattern))
            if pattern and "." not in pattern:
                counts[pattern] += 1
        total_qc.update(qc)
        if selecting:
            candidates = [evaluate(counts, qc["reads_eligible"], len(anchors), settings, k, score=True)
                          for k in range(1, limit+1)]
            winner = select_count(candidates)
            for c in candidates:
                selection.writerow([label, r.region_id, c.k, qc["reads_eligible"], sum(counts.values()),
                    len(counts), c.status, ",".join(map(str, c.support)), c.ties, c.cost,
                    c.score, int(c is winner)])
            if not fitting:
                continue
            chosen = winner or candidates[0]
            status = chosen.status
            medoids = list(chosen.medoids) if winner else []
            cost = chosen.cost if winner else None
        elif fitting:
            status, medoids, cost = fit_region(counts, len(rows), len(anchors), settings, requested)
        if fitting:
            model["regions"].append(dict(region_id=r.region_id, chrom=r.chrom, start=r.start,
                end=r.end, contig_length=length, reference_sha256=digest, cpg_anchors=anchors,
                medoids=medoids, fit_status=status, training_mismatch_cost=cost))
        else:
            status, medoids = record["fit_status"], record["medoids"]
        outcomes = Counter()
        for rg, name, pattern in sorted(rows):
            cluster, reason, *distances = assign(pattern, medoids)
            outcomes[reason or cluster] += 1
            reads_out.writerow([label, r.region_id, rg, name, pattern, len(anchors)-pattern.count("."),
                                len(anchors), cluster, reason, *distances[:-1],
                                *([None] * (limit-len(medoids))), distances[-1]])
        order = [format(n, "03b") for n in range(8)] if len(anchors) == 3 else sorted(counts)
        for pattern in order:
            patterns.writerow([label, r.region_id, pattern, counts[pattern]])
        complete = sum(counts.values())
        used = complete if fitting and status in ("fitted", "insufficient_cluster_support") else 0
        summary.writerow([label, r.region_id, status, len(anchors), qc["reads_fetched"],
            qc["reads_excluded_flags"], qc["reads_excluded_mapq"], len(rows), complete,
            outcomes["partial"], outcomes["all_missing"], len(counts), used, *(outcomes[c] for c in cluster_names),
            sum(outcomes[k] for k in REASONS), outcomes["tie"], outcomes["no_model"]])
    require_modification({label: total_qc}, settings["modification"])


def run_cluster(args):
    fitting = args.mode == "fit"
    training = args.mode != "apply"
    if training:
        requested = getattr(args, "clusters", "auto")
        if requested == "auto":
            limit = 6 if args.max_clusters is None else args.max_clusters
        else:
            if args.max_clusters is not None:
                raise ValueError("--max-clusters is for cluster select or fit --clusters auto.")
            if not requested.isascii() or not requested.isdecimal():
                raise ValueError("--clusters requires a positive integer or auto.")
            requested = limit = int(requested)
        if limit < 1:
            raise ValueError("The cluster count must be positive.")
        settings = {name: getattr(args, name) for name in SETTINGS}
        validate_settings(settings)
        regions = read_regions_bed(args.regions_bed)
        check_regions(regions)
        model = dict(schema=SCHEMA, software_version=__version__, algorithm=ALGORITHM,
                     assignment="complete_only_untied", settings=settings, regions=[],
                     clusters=requested, max_clusters=limit)
        if args.full_read:
            minimum = 3 if args.min_shared_cpgs is None else args.min_shared_cpgs
            if minimum < 1:
                raise ValueError("--min-shared-cpgs must be positive.")
            model.update(schema=full_read.IDENTITY[0], algorithm=full_read.IDENTITY[1],
                         assignment=full_read.ASSIGNMENT, min_shared_cpgs=minimum)
        elif args.min_shared_cpgs is not None:
            raise ValueError("--min-shared-cpgs requires --full-read.")
    else:
        model, regions = load_model(args.model)
        requested, limit = model["clusters"], model["max_clusters"]
    selecting = training and requested == "auto"
    label = args.sample_label if args.sample_label is not None else args.input.stem
    if not label or any(c in label for c in "\t\n\r"):
        raise ValueError("Sample label must be nonempty and contain no tabs or newlines.")
    if args.output.exists():
        raise ValueError(f"Output already exists: {args.output}")
    if not Path(str(args.reference)+".fai").is_file():
        raise ValueError("Missing FASTA index (.fai).")
    index = _index_path(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=args.output.parent) as tmp:
        dest = Path(tmp)/"results"
        dest.mkdir()
        with ExitStack() as stack:
            bam = stack.enter_context(pysam.AlignmentFile(str(args.input), "rb", index_filename=str(index)))
            fasta = stack.enter_context(pysam.FastaFile(str(args.reference)))
            if not bam.is_bam or bam.header.to_dict().get("HD", {}).get("SO") != "coordinate":
                raise ValueError("Input must be a coordinate-sorted BAM with @HD SO:coordinate.")
            if model["schema"] == full_read.IDENTITY[0]:
                full_read.analyze_full_reads(args, model, regions, fasta, bam, dest, stack, label)
            else:
                analyze_regions(args, model, regions, fasta, bam, dest, stack, label)
            if fitting:
                model["training_sample"] = label
                if selecting:
                    model["selection_method"] = ("shared_cpg_silhouette_abstentions_zero"
                        if args.full_read else "complete_read_weighted_silhouette_ties_zero")
                (dest/"model.json").write_text(json.dumps(model, indent=2)+"\n")
        dest.rename(args.output)
    print(f"methyldiver: cluster {args.mode} complete | {len(regions)} loci | {args.output}")
