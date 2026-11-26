# Regional read clustering

Use `cluster fit` to learn representative methylation patterns at each predefined locus. Use
`cluster apply` to assign reads from another sequencing run to the same patterns.
All clustering commands start from the aligned modBAM. The BED selects the loci independently
of any previous ED analysis.

```bash
methyldiver cluster fit --input training.bam --reference genome.fa \
  --regions-bed loci.bed --sample-label training --output cluster_training

methyldiver cluster apply --input heldout.bam --reference genome.fa \
  --model cluster_training/model.json --sample-label heldout --output cluster_heldout
```

Supply a coordinate-sorted BAM with an index, its indexed reference FASTA and a new
output directory. BED3 and BED4 use zero-based, half-open intervals. Regions must be
nonoverlapping and have unique IDs. BED4 column four supplies the ID. In BED3 the
coordinate string serves as the ID.

In the default regional mode, each fully contained reference CpG contributes one feature, ordered by the forward
reference cytosine coordinate. A locus must contain at least three CpGs. Each
eligible primary molecule contributes one row per locus. Calls use the existing
alignment, strand and modification rules. The defaults are 5mC, confidence 0.8,
MAPQ 20 and base quality 10. Choose a single modification at fit time. Apply
inherits the settings and checks each local reference sequence digest, contig
length and CpG list.

## Full-read comparison

Use `--full-read` with `cluster fit` or `cluster select` to use the BED as a
read selector. An eligible primary read enters a locus when at least one aligned
query base overlaps the interval. An interval may contain zero CpGs. A deletion
or reference skip spanning the entire interval contributes to `excluded_gap`.
For selected reads, the existing modification caller collects confident calls
across the whole primary alignment, including both flanks. The output uses the
forward reference cytosine coordinate for either strand of each CpG dyad.

```bash
methyldiver cluster fit --full-read --input training.bam --reference genome.fa \
  --regions-bed loci.bed --output full_read_fit

methyldiver cluster apply --input heldout.bam --reference genome.fa \
  --model full_read_fit/model.json --output full_read_test
```

Each distance is the number of discordant calls divided by the number of reference
CpGs confidently called on both molecules. For example, one disagreement among
four shared CpGs gives a distance of 0.25, regardless of how many other CpGs
either molecule covers. The default minimum is three shared CpGs. Set
`--min-shared-cpgs N` at training time to adjust it. Apply inherits the
threshold, feature mode and medoids. It checks the sequence digest of each
analyzed contig, including the flanks.

A training read is usable when it contains at least that many confident calls.
The minimum training depth (`--min-training-reads`, default five) and minimum
untied support per cluster (`--min-cluster-reads`, default two) still apply.
Medoids are observed molecules. Candidate centers must be mutually comparable
and differ at one or more shared CpGs. Greedy initialization followed by medoid
swaps first maximizes the number of training reads comparable to every center,
then minimizes their summed nearest distances. Each distance uses its own
pair of reads' shared sites. A read needs a defined distance to every medoid
and a unique nearest medoid to receive a cluster assignment.

Use `cluster select --full-read` or `cluster fit --full-read --clusters auto`
to compare cluster counts. For the adjusted silhouette, each group mean uses
the available pairwise distances. Reads with an undefined group mean, singleton
groups and abstentions contribute zero to the average over all usable reads.

In full-read mode, `assignments.tsv` contains the exact BED0 `cpg_anchors`
and the corresponding binary `calls` for each molecule. A call outside every
training medoid is still recorded on a held-out read. Each `distance_C1`,
`distance_C2`, etc. is paired with its denominator (`shared_cpgs_C1`,
`shared_cpgs_C2`, etc.). A blank distance indicates insufficient overlap or
an unavailable medoid. `regions.tsv` reports assignment counts and reasons for
abstention. Fit adds the frozen sparse profiles in `model.json`. Selection
writes `cluster_selection.tsv` with support, ties, insufficient-overlap counts
and the sum of normalized distances among comparable reads.

## Choosing a cluster count

The default is two clusters per locus. For a known design, set a positive integer
with `cluster fit --clusters 3`. To examine the training evidence first, run:

```bash
methyldiver cluster select --input training.bam --reference genome.fa \
  --regions-bed loci.bed --max-clusters 6 --output count_check
```

Open `count_check/cluster_selection.tsv`. Each row gives one locus and candidate
count, the number of complete reads, uniquely assigned reads per cluster, ties,
total mismatch cost, an adjusted silhouette score and a recommendation flag.
`select` writes only this table. To learn the recommended number independently at
each locus and save the model in one pass, use:

```bash
methyldiver cluster fit --input training.bam --reference genome.fa \
  --regions-bed loci.bed --clusters auto --max-clusters 6 --output cluster_auto
```

The search considers 1 through `--max-clusters` (default 6) on the same complete
training molecules. Candidates must pass the depth, distinct-pattern and
per-cluster support requirements. For each supported split, MethylDiver computes
the average silhouette on uniquely assigned reads using Hamming distances and
each pattern's read count. For each read it compares mean distance to other
reads in its own cluster with the lowest mean distance to another cluster.
Singletons contribute zero. Ties are excluded from those cluster means and
contribute zero to the final average over **all complete training reads**.
[Silhouette definition](https://doi.org/10.1016/0377-0427(87)90125-7)

The highest supported adjusted score wins. One cluster provides a zero baseline
and its silhouette field is blank. Equal scores to 12 decimal places favor the
smaller count. If all candidates fail, every recommendation flag is zero.
Read the support counts alongside the score, particularly when a small cluster
accounts for a slight improvement. The recommendation describes separation
among the candidate patterns observed in this training sample. A biological
interpretation for the groups requires independent evidence.

## Regional fitting and assignment

Fit uses molecules with a confident call at every model CpG. It requires five
complete reads (`--min-complete-reads`) and two untied training reads per medoid
(`--min-cluster-reads`) by default. Use compact intervals because completeness
usually falls as the number of CpGs increases. The input cap is 100,000 fetched
alignments per locus (`--max-reads-per-region`).

The distance between two complete reads is the fraction of CpGs with opposite
states. Each medoid is an observed pattern. PAM medoid swaps minimize the
read-count-weighted sum of nearest distances until no single swap improves it.
The first seed is the most frequent pattern, each subsequent seed maximizes its distance
to the nearest existing seed. Frequency and lexical order resolve initialization ties.
Equal-cost improving swaps use the lexically smallest pair. This procedure
finds a local optimum. Sorted medoids receive stable IDs C1, C2 and so on.

A complete read is assigned to its unique nearest medoid. A tie for nearest
medoid remains unassigned and has margin zero. A partial or all-missing read has
blank distances. The margin is the absolute difference between the two smallest distances. A one-cluster fit has a blank margin.
Parent-of-origin interpretation requires independent genetic evidence.

## Regional output files

| File | Contents |
| --- | --- |
| `cluster_selection.tsv` (select or automatic fit) | Per-candidate scores, support, fit status and the recommended count at each locus. |
| `model.json` (fit) | Versioned schema, settings, reference checks, fit status, CpG anchors and frozen medoids for every locus. |
| `assignments.tsv` | One row per eligible molecule per locus with read group, read ID, ordered call string, observed and model CpG counts, cluster, reason, distances to the medoids and margin. |
| `regions.tsv` | Training fit status, sample-specific alignment exclusions, eligible and complete read counts, support for each cluster, and unassigned counts. |
| `patterns.tsv` | Counts of complete ordered patterns. Three-CpG loci include all eight patterns, including zeros. Longer loci list observed patterns. |

The assignment denominator is all eligible primary molecules overlapping the locus.
Eligible equals complete + partial + all-missing and also equals the sum of all assigned cluster counts plus unassigned.
Unassigned equals ties + partial + all-missing + complete reads without a fitted
model. The fit status always refers to the training run. For example,
`no_coverage` means zero eligible training molecules. Compare fetched and
excluded counts to distinguish zero overlapping alignments from filtering.
