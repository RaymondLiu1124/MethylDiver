"""Clustering molecules selected by a BED interval using sparse calls along the alignment."""
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from math import fsum


IDENTITY = ("methyldiver.cluster/full-read/1", "shared_cpg_pam")
ASSIGNMENT = "unique_nearest_with_overlap_to_all_medoids"
STATUSES = {"fitted", "no_coverage", "insufficient_usable_reads", "insufficient_distinct_patterns",
            "insufficient_comparable_medoids", "insufficient_cluster_support"}


class ReferenceCpGs:
    """Check the two reference bases around a projected anchor, including block boundaries."""

    def __init__(self, fasta, chrom):
        self.length = fasta.get_reference_length(chrom)
        self.block_size = 65536
        width, length = self.block_size, self.length

        @lru_cache(maxsize=32)
        def block(index):
            start = index * width
            return fasta.fetch(chrom, start, min(start+width+1, length)).upper()
        self.block = block

    def __contains__(self, anchor):
        if not 0 <= anchor < self.length-1:
            return False
        index, offset = divmod(anchor, self.block_size)
        return self.block(index)[offset:offset+2] == "CG"


def contig_digest(fasta, chrom, length):
    digest = sha256()
    for start in range(0, length, 1048576):
        digest.update(fasta.fetch(chrom, start, min(start+1048576, length)).upper().encode("ascii"))
    return digest.hexdigest()


def unpack(profile):
    return [v//2 for v in profile], "".join(str(v%2) for v in profile)


def load_profiles(record, minimum):
    profiles = []
    for medoid in record["medoids"]:
        anchors, calls = medoid["cpg_anchors"], medoid["calls"]
        if (not isinstance(anchors, list) or not isinstance(calls, str) or len(calls) != len(anchors)
                or len(calls) < minimum or set(calls)-set("01")
                or any(type(a) is not int or not 0 <= a < record["contig_length"]-1 for a in anchors)
                or anchors != sorted(set(anchors))):
            raise ValueError("Invalid full-read medoid coordinates or calls.")
        profiles.append(tuple(2*a+int(v) for a,v in zip(anchors,calls)))
    if profiles != sorted(set(profiles)):
        raise ValueError("Full-read medoids must be distinct and in canonical order.")
    return profiles


class Comparisons:
    """Pair-specific Hamming distances over observed coordinates, with bounded caching."""

    def __init__(self, profiles, minimum):
        positions = sorted({v//2 for p in profiles for v in p})
        indices = {p:i for i,p in enumerate(positions)}
        self.minimum = minimum
        masks = {}
        for profile in profiles:
            observed = modified = 0
            for value in profile:
                bit = 1 << indices[value//2]
                observed |= bit
                if value%2:
                    modified |= bit
            masks[profile] = (observed, modified)

        @lru_cache(maxsize=65536)
        def compare(a, b):
            mask_a, ones_a = masks[a]
            mask_b, ones_b = masks[b]
            shared = mask_a & mask_b
            overlap = shared.bit_count()
            distance = ((ones_a ^ ones_b) & shared).bit_count()/overlap if overlap >= minimum else None
            return distance, overlap
        self.compare = compare

    def assignment(self, p, centers):
        pairs = [self.compare(p, m) for m in centers]
        distances = [d for d, n in pairs]
        margin = None
        if not p:
            reason = "all_missing"
        elif len(p) < self.minimum:
            reason = "too_few_calls"
        elif not centers:
            reason = "no_model"
        elif None in distances:
            reason = "insufficient_overlap"
        else:
            nearest = sorted(distances)
            if len(nearest) > 1:
                margin = nearest[1]-nearest[0]
            if margin != 0:
                return distances.index(nearest[0]), "", pairs, margin
            reason = "tie"
        return None, reason, pairs, margin


def choose_medoids(counts, k, comparisons):
    """Prioritize comparable training coverage, then minimize summed normalized distance."""
    profiles = sorted(counts)
    distance = comparisons.compare

    def objective(centers):
        missing = 0
        costs = []
        for p in profiles:
            ds = [distance(p,m)[0] for m in centers]
            if None in ds:
                missing += counts[p]
            else:
                costs.append(counts[p]*min(ds))
        return missing, round(fsum(costs), 12)

    def distinct(candidate, kept):
        return all(distance(candidate,m)[0] not in (None, 0) for m in kept)

    centers = ()
    while len(centers) < k:
        options = (p for p in profiles if distinct(p, centers))
        # Among equally good seeds prefer well-supported, well-covered reads.
        seed = min(options, key=lambda p: (objective((*centers,p)), -counts[p], -len(p), p), default=None)
        if seed is None:
            return (), None
        centers = tuple(sorted((*centers,seed)))
    score = objective(centers)
    while True:
        best = score, centers
        for removed in centers:
            kept = tuple(p for p in centers if p != removed)
            for p in profiles:
                if p in centers or not distinct(p, kept):
                    continue
                trial = tuple(sorted((*kept,p)))
                trial_score = objective(trial)
                if trial_score < score:
                    best = min(best, (trial_score, trial))
        if best[0] == score:
            return centers, score[1]
        score, centers = best


def separation(groups, total, comparisons):
    """Mean silhouette over usable molecules. Missing group evidence and abstentions score zero."""
    contributions = []
    for own, group in enumerate(groups):
        for p,n in group.items():
            means = []
            for i, other in enumerate(groups):
                size = 0
                sums = []
                for q, weight in other.items():
                    d, _ = comparisons.compare(p,q)
                    if d is None:
                        continue
                    weight -= int(i == own and p == q)
                    size += weight
                    sums.append(d*weight)
                if not size:
                    break
                means.append(fsum(sums)/size)
            if len(means) == len(groups):
                a = means[own]
                b = min(v for i,v in enumerate(means) if i != own)
                if max(a,b):
                    contributions.append(n*(b-a)/max(a,b))
    return fsum(contributions)/total


def evaluate_profiles(counts, eligible, settings, k, comparisons, score):
    from .cluster import Candidate
    if not eligible:
        return Candidate(k, "no_coverage")
    total = sum(counts.values())
    if total < settings["min_complete_reads"]:
        return Candidate(k, "insufficient_usable_reads")
    if len(counts) < k:
        return Candidate(k, "insufficient_distinct_patterns")
    centers, cost = choose_medoids(counts, k, comparisons)
    if not centers:
        return Candidate(k, "insufficient_comparable_medoids")
    groups = [Counter() for _ in centers]
    ties = 0
    for p,n in counts.items():
        index, reason, _, _ = comparisons.assignment(p, centers)
        if index is not None:
            groups[index][p] = n
        elif reason == "tie":
            ties += n
    support = tuple(sum(g.values()) for g in groups)
    status = "fitted" if min(support) >= settings["min_cluster_reads"] else "insufficient_cluster_support"
    value = separation(groups,total,comparisons) if score and k>1 and status=="fitted" else None
    return Candidate(k, status, centers, cost, support, ties, value)


def analyze_full_reads(args, model, regions, fasta, bam, dest, stack, label):
    from .cluster import select_count, tsv
    from .modbam import decode_calls, eligible_reads, require_modification
    training = args.mode != "apply"
    selecting = training and model["clusters"] == "auto"
    settings, minimum, limit = model["settings"], model["min_shared_cpgs"], model["max_clusters"]
    clusters = [f"C{i}" for i in range(1,limit+1)]
    if selecting:
        selection = tsv(stack, dest, "cluster_selection.tsv",
            "sample region_id k eligible_reads usable_reads distinct_profiles status cluster_support "
            "tied_reads insufficient_overlap_reads normalized_distance_cost adjusted_silhouette recommended")
    if args.mode != "select":
        reads = tsv(stack, dest, "assignments.tsv",
            "sample region_id read_group read_id cpg_anchors calls observed_cpgs cluster reason "
            + " ".join("distance_"+c+" shared_cpgs_"+c for c in clusters)+" margin")
        summary = tsv(stack, dest, "regions.tsv",
            "sample region_id fit_status fetched excluded_flags excluded_mapq excluded_gap eligible usable "
            + " ".join("assigned_"+c for c in clusters)
            + " unassigned all_missing too_few_calls insufficient_overlap tie no_model")
    digests = {}
    total_qc = Counter()
    for i,r in enumerate(regions):
        if r.chrom not in fasta.references or r.chrom not in bam.references:
            raise ValueError(f"{r.region_id}: contig missing from BAM or FASTA.")
        length = fasta.get_reference_length(r.chrom)
        if r.end > length or bam.get_reference_length(r.chrom) != length:
            raise ValueError(f"{r.region_id}: incompatible reference length or boundaries.")
        reference = ReferenceCpGs(fasta,r.chrom)
        if args.mode != "select" and r.chrom not in digests:
            digests[r.chrom] = contig_digest(fasta,r.chrom,length)
        centers = []
        if not training:
            record = model["regions"][i]
            if length != record["contig_length"] or digests[r.chrom] != record["reference_sha256"]:
                raise ValueError(f"{r.region_id}: full-read model/reference mismatch.")
            centers = load_profiles(record, minimum)
            if any(v//2 not in reference for m in centers for v in m):
                raise ValueError(f"{r.region_id}: medoid anchor is not a reference CpG.")
        qc, counts = Counter(), Counter()
        rows = []
        for read in eligible_reads(bam, r, settings, qc):
            if not read.get_overlap(r.start,r.end):
                qc["excluded_gap"] += 1
                qc["reads_eligible"] -= 1
                continue
            calls, _ = decode_calls(read, reference, settings["call_threshold"], settings["min_base_quality"],
                                    qc, settings["modification"])
            profile = tuple(sorted(2*a+v for a,v in calls.items() if v is not None))
            counts[profile] += 1
            if args.mode != "select":
                rows.append((read.get_tag("RG") if read.has_tag("RG") else "", read.query_name, profile))
        total_qc.update(qc)
        usable = {p:n for p,n in counts.items() if len(p)>=minimum}
        n_usable = sum(usable.values())
        comparisons = Comparisons(set(counts).union(centers), minimum)
        if training:
            ks = range(1,limit+1) if selecting else [model["clusters"]]
            candidates = [evaluate_profiles(usable,qc["reads_eligible"],settings,k,comparisons,selecting) for k in ks]
            selected = select_count(candidates)
            if selecting:
                for c in candidates:
                    lost = n_usable - sum(c.support) - c.ties if c.support else None
                    selection.writerow([label,r.region_id,c.k,qc["reads_eligible"],n_usable,len(usable),c.status,
                        ",".join(map(str,c.support)),c.ties,lost,c.cost,c.score,int(c is selected)])
            if args.mode == "select":
                continue
            status = (selected or candidates[0]).status
            centers = selected.medoids if selected else ()
            stored = []
            for center in centers:
                anchors, calls = unpack(center)
                stored.append(dict(cpg_anchors=anchors, calls=calls))
            model["regions"].append(dict(region_id=r.region_id, chrom=r.chrom, start=r.start, end=r.end,
                contig_length=length, reference_sha256=digests[r.chrom], fit_status=status, medoids=stored))
        else:
            status = record["fit_status"]
        outcomes = Counter()
        for rg,name,p in sorted(rows):
            index, reason, pairs, margin = comparisons.assignment(p,centers)
            cluster = clusters[index] if index is not None else "unassigned"
            outcomes[reason or cluster] += 1
            anchors, calls = unpack(p)
            fields = [v for pair in pairs for v in pair] + [None]*(2*(limit-len(centers)))
            reads.writerow([label,r.region_id,rg,name,",".join(map(str,anchors)),calls,len(p),
                            cluster,reason,*fields,margin])
        reasons = ["all_missing","too_few_calls","insufficient_overlap","tie","no_model"]
        summary.writerow([label,r.region_id,status,qc["reads_fetched"],qc["reads_excluded_flags"],
            qc["reads_excluded_mapq"],qc["excluded_gap"],qc["reads_eligible"],n_usable,
            *(outcomes[c] for c in clusters),sum(outcomes[s] for s in reasons),*(outcomes[s] for s in reasons)])
    require_modification({label:total_qc},settings["modification"])
