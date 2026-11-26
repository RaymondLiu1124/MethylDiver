from collections import Counter
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam

from methyldiver.cluster import SETTINGS, assign, evaluate, fit_medoids, fit_region, select_count
from methyldiver.modbam import Region, analyze_single_modbam
from test_modbam import HEADER, SEQUENCE, make_read


class TestCluster(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.ref = self.root/"ref.fa"
        self.ref.write_text(">chrToy\n"+SEQUENCE+"\n")
        pysam.faidx(str(self.ref))
        self.bed = self.root/"loci.bed"
        self.bed.write_text("chrToy\t2\t12\ttoy\n")
        self.reads = [make_read("a"+str(i), "001", reverse=i%2 == 1) for i in range(3)]
        self.reads += [make_read("b"+str(i), "100", opposite=i%2 == 1) for i in range(3)]
        self.reads += [make_read("tie", "000"), make_read("partial", ml=[255, 127, 0]),
                       make_read("missing", no_tags=True)]
        self.bam = self.save_bam("train", self.reads)

    def save_bam(self, name, reads):
        path = self.root/(name+".bam")
        with pysam.AlignmentFile(str(path), "wb", header=HEADER) as bam:
            for read in reads:
                bam.write(read)
        pysam.index(str(path))
        return path

    def cli(self, mode, out, bam=None, model=None, ref=None, ok=True, options=()):
        args = ["--model", model] if mode == "apply" else ["--regions-bed", self.bed]
        result = subprocess.run([sys.executable, "-m", "methyldiver", "cluster", mode,
            "--input", bam or self.bam, "--reference", ref or self.ref, "--sample-label", "test",
            "--output", self.root/out, *args, *options], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0 if ok else 2, result.stderr)
        if not ok:
            self.assertFalse((self.root/out).exists())
        return result

    def table(self, output, name):
        with (self.root/output/name).open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_weighted_objective_and_statuses(self):
        counts = Counter({"000":1,"001":1,"011":1,"110":10,"111":10})
        self.assertEqual(fit_medoids(counts), (["110","111"],5))
        for counts, n, k, status in [
                ({},0,3,"no_coverage"), ({},7,3,"insufficient_complete_reads"),
                ({"11":8},8,2,"fewer_than_three_cpgs"), ({"111":6},6,3,"one_distinct_pattern"),
                ({"000":1,"111":9},10,3,"insufficient_cluster_support")]:
            self.assertEqual(fit_region(Counter(counts), n, k, SETTINGS)[0], status)
        self.assertEqual(assign("0000", ["0011","1100"]), ("unassigned","tie",0.5,0.5,0))
        self.assertEqual(assign("00.1", ["0011","1100"]), ("unassigned","partial",None,None,None))

    def test_fit_apply_counts_freezing_and_reference(self):
        self.cli("fit", "fit")
        model_file = self.root/"fit/model.json"
        model = json.loads(model_file.read_text())
        locus = model["regions"][0]
        self.assertEqual(locus["cpg_anchors"], [2,6,10])
        self.assertEqual(locus["medoids"], ["001","100"])
        counts = [int(row["complete_reads"]) for row in self.table("fit","patterns.tsv")]
        old = analyze_single_modbam(self.bam, self.ref, [Region("chrToy", 2, 12)])
        self.assertEqual(counts, old["units"][0]["sample"]["counts"])
        self.cli("apply", "same", model=model_file)
        self.assertEqual(self.table("fit","assignments.tsv"), self.table("same","assignments.tsv"))
        rows = self.table("same","assignments.tsv")
        self.assertEqual(Counter(r["cluster"] for r in rows), {"C1":3,"C2":3,"unassigned":3})
        self.assertEqual(Counter(r["reason"] for r in rows if r["reason"]), {"tie":1,"partial":1,"all_missing":1})
        self.assertEqual([float(r["distance_"+r["cluster"]]) for r in rows if not r["reason"]], [0]*6)
        summary, = self.table("same","regions.tsv")
        self.assertEqual(int(summary["eligible"]), len(rows))
        self.assertEqual(int(summary["unassigned"]),3)
        for i,r in enumerate(self.reads):
            r.set_tag("HP", 1+i%2)
        shuffled = self.save_bam("shuffled", reversed(self.reads))
        self.cli("fit", "shuffled_fit", bam=shuffled)
        for name in ("model.json","assignments.tsv","patterns.tsv","regions.tsv"):
            self.assertEqual((self.root/"fit"/name).read_bytes(), (self.root/"shuffled_fit"/name).read_bytes())
        heldout = self.save_bam("heldout", [make_read(str(i), "100") for i in range(7)])
        self.cli("apply", "heldout_result", bam=heldout, model=model_file)
        self.assertEqual([r["cluster"] for r in self.table("heldout_result","assignments.tsv")], ["C2"]*7)
        ref = self.root/"other.fa"
        ref.write_text(">chrToy\n"+SEQUENCE.replace("CG","TG",1)+"\n")
        pysam.faidx(str(ref))
        self.assertIn("mismatch",self.cli("apply","changed",model=model_file,ref=ref,ok=False).stderr)
        model["regions"][0]["medoids"] = ["00","111"]
        corrupt = self.root/"corrupt.json"
        corrupt.write_text(json.dumps(model))
        self.assertIn("medoids",self.cli("apply","invalid",model=corrupt,ok=False).stderr)

    def test_four_cpgs_are_unique_features_and_alignment_filters_reconcile(self):
        seq = SEQUENCE + "CGAA"
        header = {"HD":{"VN":"1.6","SO":"coordinate"}, "SQ":[{"SN":"chrToy","LN":len(seq)}]}
        ref = self.root/"four.fa"
        ref.write_text(">chrToy\n"+seq+"\n")
        pysam.faidx(str(ref))
        reads = [make_read(str(i), "0011" if i%2 else "1100", sequence=seq, header=header)
                 for i in range(8)]
        reads[0].mapping_quality = 19
        reads[1].flag |= 2048
        self.bed.write_text(f"chrToy\t0\t{len(seq)}\tfour\n")
        bam = self.root/"four.bam"
        with pysam.AlignmentFile(str(bam),"wb",header=header) as out:
            for read in reads:
                out.write(read)
        pysam.index(str(bam))
        self.cli("fit","four",bam=bam,ref=ref)
        record, = json.loads((self.root/"four/model.json").read_text())["regions"]
        self.assertEqual(record["cpg_anchors"], [2,6,10,14])
        self.assertEqual(record["medoids"], ["0011","1100"])
        summary, = self.table("four","regions.tsv")
        self.assertEqual([summary[k] for k in ("fetched","excluded_mapq","excluded_flags","eligible","complete")],
                         ["8","1","1","6","6"])

    def test_uncovered_fit_retains_the_locus_on_apply(self):
        empty = self.save_bam("empty", [])
        self.cli("fit", "empty_fit", bam=empty)
        model = self.root/"empty_fit/model.json"
        self.assertEqual(json.loads(model.read_text())["regions"][0]["fit_status"], "no_coverage")
        self.cli("apply", "untrained", model=model)
        summary, = self.table("untrained", "regions.tsv")
        self.assertEqual([summary[k] for k in ("eligible","no_model","partial","all_missing","unassigned")],
                         ["9","7","1","1","9"])

    def test_boundary_excludes_half_dyad(self):
        self.bed.write_text("chrToy\t2\t11\tshort\n")
        self.cli("fit","short")
        record = json.loads((self.root/"short/model.json").read_text())["regions"][0]
        self.assertEqual(record["cpg_anchors"], [2,6])
        self.assertEqual(record["fit_status"], "fewer_than_three_cpgs")

    def test_multicluster_selection_and_application(self):
        reads = [make_read(f"{pattern}-{i}", pattern) for pattern in ("000", "011", "101")
                 for i in range(3)]
        bam = self.save_bam("three", reads)
        options = ["--max-clusters", "4"]
        self.cli("select", "scores", bam=bam, options=options)
        self.assertEqual([p.name for p in (self.root/"scores").iterdir()], ["cluster_selection.tsv"])
        candidates = self.table("scores", "cluster_selection.tsv")
        self.assertEqual([c["k"] for c in candidates if c["recommended"] == "1"], ["3"])
        self.assertEqual(candidates[2]["adjusted_silhouette"], "1.0")
        self.assertEqual(candidates[3]["status"], "insufficient_distinct_patterns")
        self.cli("fit", "auto", bam=bam, options=["--clusters", "auto", *options])
        self.assertEqual(candidates, self.table("auto", "cluster_selection.tsv"))
        model_path = self.root/"auto/model.json"
        model = json.loads(model_path.read_text())
        self.assertEqual(model["regions"][0]["medoids"], ["000", "011", "101"])
        self.cli("fit", "fixed", bam=bam, options=["--clusters", "3"])
        self.assertEqual(json.loads((self.root/"fixed/model.json").read_text())["regions"], model["regions"])
        changed = self.save_bam("changed", [make_read("z", "101"), make_read("t", "001"),
                                            make_read("p", ml=[255, 127, 0])])
        self.cli("apply", "predictions", bam=changed, model=model_path)
        by_id = {r["read_id"]: r for r in self.table("predictions", "assignments.tsv")}
        self.assertEqual(by_id["z"]["cluster"], "C3")
        self.assertEqual(by_id["z"]["distance_C3"], "0.0")
        self.assertEqual(by_id["z"]["distance_C4"], "")
        self.assertAlmostEqual(float(by_id["z"]["margin"]), 2/3)
        self.assertEqual(by_id["t"]["reason"], "tie")
        self.assertEqual(by_id["p"]["distance_C3"], "")

    def test_baseline_and_unsupported_selection(self):
        counts = Counter({"111": 6})
        baseline = evaluate(counts, 6, 3, SETTINGS, 1, score=True)
        self.assertIsNone(baseline.score)
        self.assertEqual(select_count([baseline, evaluate(counts, 6, 3, SETTINGS, 2)]).k, 1)
        self.assertEqual(assign("100", ["111"]), ("C1", "", 2/3, None))
        self.assertIsNone(select_count([evaluate(Counter(), 0, 3, SETTINGS, 1)]))
        for args in (["--clusters", "0"], ["--clusters", "two"], ["--max-clusters", "4"],
                     ["--clusters", "auto", "--max-clusters", "0"]):
            self.cli("fit", "invalid", options=args, ok=False)

    def test_silhouette_weighting_ties_and_singletons(self):
        # Distinct groups have three 001 and three 100 reads. 000 ties.
        result = evaluate(Counter({"001":3, "100":3, "000":1}), 7, 3, SETTINGS, 2, score=True)
        self.assertEqual(result.ties, 1)
        self.assertAlmostEqual(result.score, 6/7)
        # 001 and 000 form one 4-read cluster; 111 forms a singleton.
        # Within-cluster mean for 000 is 1/3, between-cluster mean is 3.
        # For 001 those means are 1 and 2. Singleton contributes zero.
        result = evaluate(Counter({"000":3, "001":1, "111":1}), 5, 3,
                          dict(SETTINGS, min_cluster_reads=1), 2, score=True)
        self.assertAlmostEqual(result.score, (3*(1-1/9) + .5)/5)

    def test_schema_one_application(self):
        self.cli("fit", "fit")
        model = json.loads((self.root/"fit/model.json").read_text())
        model.update(schema="methyldiver.cluster/1", algorithm="weighted_pam_two_normalized_hamming")
        del model["clusters"], model["max_clusters"]
        path = self.root/"legacy.json"
        path.write_text(json.dumps(model))
        self.cli("apply", "legacy", model=path)
        self.assertEqual(self.table("legacy","assignments.tsv"), self.table("fit","assignments.tsv"))

    def test_full_read_uses_aligned_flanks_and_selects_by_bed(self):
        # Only CpG 6 is inside the BED. All molecules call it 0, while both flanks differ.
        self.bed.write_text("chrToy\t6\t8\tgate\n")
        training = [make_read(str(i), "000" if i<3 else "101", reverse=i%2==1) for i in range(6)]
        skipped = make_read("gap", sequence=SEQUENCE[:6]+SEQUENCE[8:], cigar="6M2N6M",
                            mm="C+m?,0,0;", ml=[0,255])
        outside = make_read("outside", sequence=SEQUENCE[:4], mm="C+m?,0;", ml=[255])
        bam = self.save_bam("flanks", [*training, skipped, outside])
        self.cli("fit","regional",bam=bam)
        self.assertEqual(self.table("regional","regions.tsv")[0]["fit_status"], "fewer_than_three_cpgs")
        opts = ["--full-read","--min-shared-cpgs","2"]
        self.cli("fit","full",bam=bam,options=opts)
        record = json.loads((self.root/"full/model.json").read_text())["regions"][0]
        self.assertEqual(record["medoids"], [dict(cpg_anchors=[2,6,10],calls=p) for p in ("000","101")])
        summary, = self.table("full","regions.tsv")
        self.assertEqual([summary[key] for key in ["eligible","excluded_gap","assigned_C1","assigned_C2"]],
                         ["6","1","3","3"])
        model = self.root/"full/model.json"
        self.cli("apply","same_full",bam=bam,model=model)
        self.assertEqual(self.table("full","assignments.tsv"),self.table("same_full","assignments.tsv"))
        # Internal uncertain call, a soft-clipped C and a non-CpG mismatch supply no feature.
        clipped = make_read("clipped",sequence="C"+"C"+SEQUENCE[1:],cigar="1S14M",
                            mm="C+m?,0,0,0,0,0;",ml=[255,255,255,127,255])
        partial = self.save_bam("partial_full", [clipped, make_read("low",ml=[127,255,127])])
        self.cli("apply","partial_full_run",bam=partial,model=model)
        by_id = {r["read_id"]:r for r in self.table("partial_full_run","assignments.tsv")}
        self.assertEqual([by_id["clipped"][k] for k in ("cpg_anchors","calls","cluster","shared_cpgs_C1")],
                         ["2,10","11","C2","2"])
        self.assertEqual(by_id["low"]["reason"],"too_few_calls")
        self.assertEqual(by_id["low"]["distance_C2"],"")
        self.cli("select","full_scores",bam=bam,options=opts+["--max-clusters","3"])
        self.cli("fit","full_auto",bam=bam,options=opts+["--clusters","auto","--max-clusters","3"])
        self.assertEqual(self.table("full_auto","cluster_selection.tsv"),self.table("full_scores","cluster_selection.tsv"))
        self.assertEqual([c["k"] for c in self.table("full_scores","cluster_selection.tsv") if c["recommended"]=="1"],["2"])
        ref = self.root/"changed_flank.fa"
        ref.write_text(">chrToy\nT"+SEQUENCE[1:]+"\n")
        pysam.faidx(str(ref))
        self.assertIn("mismatch",self.cli("apply","bad_flank",model=model,ref=ref,ok=False).stderr)

    def test_full_read_distance_uses_shared_coordinates_and_requires_comparable_centers(self):
        from methyldiver.full_read import Comparisons, choose_medoids
        # Packed (BED0 anchor * 2 + state). a and b disagree at one of two shared sites.
        a,b,c = (0,2,4),(3,4,7),(8,10,12)
        data = Comparisons([a,b,c],2)
        self.assertEqual(data.compare(a,b),(.5,2))
        self.assertEqual(data.compare(a,c),(None,0))
        index,reason,pairs,margin = data.assignment(b,[a,c])
        self.assertEqual((index,reason,pairs,margin),(None,"insufficient_overlap",[(.5,2),(None,0)],None))
        self.assertEqual(choose_medoids({a:3,c:3},2,data),((),None))
        self.assertEqual(data.assignment(b,[a,b])[0],1)
        self.cli("fit","zero_shared",options=["--full-read","--min-shared-cpgs","0"],ok=False)
        self.cli("fit","wrong_scope",options=["--min-shared-cpgs","5"],ok=False)

    def test_full_read_heldout_span_and_novel_sites(self):
        self.bed.write_text("chrToy\t6\t7\tgate\n")
        train = self.save_bam("left_flank", [make_read(str(i), "00" if i<3 else "11",
                                         sequence=SEQUENCE[:9]) for i in range(6)])
        self.cli("fit","left_model",bam=train,options=["--full-read","--min-shared-cpgs","2"])
        heldout = self.save_bam("right_flank",[make_read("long","111"),
            make_read("short","11",sequence=SEQUENCE[4:],start=4)])
        self.cli("apply","novel",bam=heldout,model=self.root/"left_model/model.json")
        long, short = self.table("novel","assignments.tsv")
        self.assertEqual([long[k] for k in ("cpg_anchors","observed_cpgs","shared_cpgs_C1","cluster")],
                         ["2,6,10","3","2","C2"])
        self.assertEqual([short[k] for k in ("cpg_anchors","observed_cpgs","shared_cpgs_C2","reason")],
                         ["6,10","2","1","insufficient_overlap"])
        self.assertEqual(short["distance_C2"], "")
