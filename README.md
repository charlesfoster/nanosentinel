# nanosentinel

Monitor a live Oxford Nanopore run from a samplesheet (in the same format as h2seq's, so one sheet can serve both):
reads are mapped as MinKNOW writes them, coverage is shown on a local dashboard, and the time each sample takes to
reach a target depth is recorded. Optionally concatenates each sample's reads for analysis afterwards.

- **Dashboard** on a local port: per-sample progress to target, genome coverage with amplicons and genes, reference
  composition, amplicon-by-sample depth (dropouts show as hatched cells), progress curves and read lengths. A
  self-contained `dashboard.html` snapshot is written when you stop, so you keep a record of the run.
- Sequencing-time **time to target**: when each sample first reached `--target_depth` over `--target_fraction` of the
  genome (default 20x over 90%), measured from the reads' own timestamps. You are notified when every non-control sample
  has reached it; monitoring continues until you stop it.
- **Protocols**: prepare a reference panel once with `prep-protocol` (align every reference to a coordinate reference,
  place genes, regions and primers on each, QC the result), then monitor with `--protocol`. Samples of different
  genotypes are drawn on one shared axis, so genes, regions and amplicons line up. A ready-made **HCV protocol** is
  bundled (`--protocol hcv`).
- Handles **several similar references** in the panel: a read that fits sibling genotypes equally well counts towards
  all of them, so coverage is not split between near-identical references.
- Waits automatically if MinKNOW has not created the run's `fastq_pass` folder yet.
- Quiet terminal: one status line and a few event messages. Everything is also logged to `nanosentinel.log`.

## Install

Needs [pixi](https://pixi.sh). The environment (Python, minimap2, MAFFT, numpy) installs itself on first run, and again if
`pixi.lock` changes. Works on Linux and macOS.

```bash
git clone <this repo> ~/Programs/nanosentinel
mkdir -p ~/.local/bin
ln -s ~/Programs/nanosentinel/bin/nanosentinel ~/.local/bin/nanosentinel
# make sure ~/.local/bin is on your PATH (add to ~/.bashrc / ~/.zshrc if it is not):
export PATH="$HOME/.local/bin:$PATH"
nanosentinel --version
```

## Samplesheet

The h2seq samplesheet, plus columns for monitoring. Only `sample` and `barcode` are required.

```csv
sample,long_reads,barcode,run_name,neg_control
patient1,,BC01,2026-11-09_run,
patient2,,BC02,2026-11-09_run,
NTC,,BC12,2026-11-09_run,true
```

| column | meaning |
| --- | --- |
| `sample` | as in h2seq (no spaces). Repeated sample IDs are fine (one row per barcode / run). |
| `barcode` | `BC01`, `barcode01`, `1` ... all mean `barcode01`. |
| `run_name` | Optional. A folder name under the MinKNOW data directory; also selects rows if the sheet spans several runs. |
| `neg_control` | Optional. `true` for negative controls: they are excluded from "all samples at target" and a control reaching target depth raises a contamination warning. |
| `long_reads` | Leave empty while sequencing; `--concatenate_reads` fills it in. |

## Monitor a run

```bash
nanosentinel monitor --input samplesheet.csv --run_name 2026-11-09_run --protocol hcv
```

or, with ordinary files and no preparation (quick mode: only the coordinate reference is annotated, so other
references are plotted without genes or amplicons):

```bash
nanosentinel monitor --input samplesheet.csv --run_name 2026-11-09_run \
    --possible_references references.fasta --primer_bed primers.bed
```

Opens the dashboard at <http://localhost:6010>. If the run has not started yet it waits for `fastq_pass` to appear
(Ctrl+C to cancel). Press **Enter** or **Ctrl+C** to stop.

In quick mode the inputs are ordinary files: `--possible_references` is a multi-FASTA panel the reads are mapped to (include one entry per
genotype or subtype you might see), `--reference_fasta` picks the reference the primer scheme and genes are drawn
against (default: the first panel entry), `--primer_bed` gives amplicons (ARTIC-style `<amplicon>_LEFT` /
`<amplicon>_RIGHT` names on that reference) and `--gff` adds gene labels.

Outputs in `nanosentinel_<run_name>/` (reports are refreshed about once a minute and written again when you stop):

| file | contents |
| --- | --- |
| `dashboard.html` | Self-contained snapshot of the dashboard at the end of monitoring (opens straight from disk). |
| `time_to_target.tsv` | When each sample first reached the target (sequencing time), with reads at that point and the reference it was reached on. |
| `coverage_timeseries.tsv` | Progress over time: one line per sample per annotation file with reads, mapped bases and the best fraction at target depth. |
| `sample_summary.tsv` | Per sample: reads, % mapped, dominant reference and its mean depth, and read length statistics (mean, median, N50, min, max). |
| `reference_summary.tsv` | Per sample per reference with coverage: primary read share, mean depth and fraction at target depth. Shows genotype composition and mixed samples. |
| `annotations/` | One small TSV per FASTQ chunk (per-read alignments). Reused if you restart with the same `--outdir`. |
| `reference/` | The reference panel and its minimap2 index. |
| `nanosentinel.log` | The tool's log. |

In `sample_summary.tsv`, `dominant_reference` is the reference with the most reads, while `best_fraction` is the highest
fraction at target depth over any reference; they can differ when a sample's reads split between similar references.

Useful options (see `--help` for all):

- `--denominator amplicons` measures the fraction over each reference's amplicons instead of its whole length (a
  reference with no amplicons cannot reach target in this mode).
- `--target_region NAME` (with `--protocol`) measures the target over a named region, for example
  `--target_region E1E2 --target_depth 20 --target_fraction 0.9` for "20x over 90% of E1E2". The region is found on each
  sample's own reference through its coordinate map and is held to its full length on the coordinate axis, so a reference
  that lacks part of it (a partial sequence, a deletion) cannot reach target by covering what it does have. Time to target,
  the dashboard bar and the notification all refer to the region; `time_to_target.tsv` records it as `region:E1E2`.
- `--target_depth`, `--target_fraction`, `--min_identity`, `--webhook URL`, `--port`, `--no_browser`.
- `--host 0.0.0.0` to view the dashboard from another machine on the network (there is no authentication).
- `--run_dir` points at an exact MinKNOW run folder (use it if a run name matches more than one folder).
- `--dry_run` checks the inputs and builds the reference index without monitoring.
- `--clear_annotations` discards earlier annotations and maps every read again.

## Concatenate reads for analysis

At the end of monitoring:

```bash
nanosentinel monitor ... --concatenate_reads
```

or for a finished run, without monitoring:

```bash
nanosentinel monitor --input samplesheet.csv --run_name 2026-11-09_run --skip_monitor --concatenate_reads
```

This writes `reads/<sample>_<run>_<barcode>.fastq.gz` and `samplesheet_with_reads.csv` (this run's rows, `long_reads`
filled in) ready to pass to an analysis pipeline. Gzipped FASTQs are joined without recompression; `--gzip` compresses
plain-text input.

## Protocols

A protocol is a reference panel prepared once and reused for every run. It takes a minute or two to build, loads in well
under a second, and refuses to be used if its files have changed since it was built.

```bash
nanosentinel prep-protocol --possible_references panel.fasta --coordinate_ref REF_ID \
    --gff genes.gff --region E1E2=915-2579 --primer_fasta primers.fa --amplicon_pattern '^(Seg[0-9]+)_' \
    --outdir protocols/my-protocol
nanosentinel monitor --input samplesheet.csv --run_name run1 --protocol protocols/my-protocol
```

`prep-protocol` aligns every reference to the **coordinate reference** with MAFFT (about a minute for 238 HCV genomes;
or give `--panel_alignment` to skip it). That gives each reference a *coordinate map*: for each of its bases, the
position of the aligned base on the coordinate reference. Depth measured on any reference is then drawn on the
coordinate reference's axis, so a 1a, a 3a and a 2b sample show the same genes, the same E1E2 region and the same
amplicon columns. Positions on the dashboard are therefore coordinate-reference numbers (H77 numbers for HCV).

| input | meaning |
| --- | --- |
| `--possible_references` | The panel reads are mapped to. Keep every reference you want considered. |
| `--coordinate_ref` / `--reference_fasta` | The reference everything is drawn on (default: first panel entry). |
| `--gff` | Gene (or CDS) features on the coordinate reference. |
| `--region NAME=START-END` | Named regions on the coordinate reference, 1-based inclusive. Repeatable. |
| `--primer_bed` | Amplicons on the coordinate reference (ARTIC-style `<amplicon>_LEFT` / `_RIGHT` names); carried onto every reference through the alignment. |
| `--primer_fasta` | Primer sequences. Their sites are **searched for on every reference** (see below). |
| `--group_pattern` | Regular expression naming a reference's sibling group (default: the text before the first `_`, which is the subtype for `1a_AF009606.1`). |

**Finding primers.** Each primer is searched for on every reference, allowing IUPAC codes (`R`, `Y`, `S`, ...). A site
counts if it has at most `--max_mismatch` mismatches (default 5), at most `--max_3prime_mismatch` of them in the last
5 bases of the primer's 3' end (default 1), and lies near where that primer belongs, meaning within 150 bases (on the
coordinate axis) of the position at least three references agree on. Without that last rule, genotype-specific primers
would pick up chance matches elsewhere. An amplicon's span on a reference runs from its best-matching LEFT primer to its
best-matching RIGHT primer, so genotype-specific primer variants resolve themselves: a reference gets an amplicon only
if both ends have a site. The mismatch limits are a judgement about how permissive annealing is; the sensitivity is
easy to see by rebuilding with different values, and `primer_sites.tsv` records the mismatches of every site.

**Outputs** (all in the protocol folder):

| file | contents |
| --- | --- |
| `protocol.json` | Name, date, nanosentinel and MAFFT versions, genes, regions, QC summary, and the hashes `monitor` checks. |
| `panel.fasta`, `coordinate_maps.npz` | The panel and the coordinate maps. |
| `amplicons.tsv` | Each amplicon's span on each reference (1-based, inclusive) and the primers and mismatches behind it. |
| `primer_sites.tsv` | Every primer site found, with mismatches. |
| `features.tsv` | Genes and regions lifted onto each reference, with a status: `ok`, `odd_span` (lifted span far from the expected length), `partial` or `absent`. |
| `qc.txt` | The QC report, also printed by `prep-protocol`: lift status per gene and region, amplicons per reference and the genotype groups missing each, poorly aligned references, primers that could not be placed. |

In the dashboard a sample is shown on the best-covered reference of its **sibling group**. The group is chosen after the
first 200 mapped reads ("taste") and kept ("lock") unless another group leads by 1.5x and 100 reads, so the display does
not flip between near-identical references that share reads. The group is reported as `reference_group` in
`sample_summary.tsv`. Time to target is still "reached on any reference".

### The HCV protocol

`protocols/hcv/` (use it as `--protocol hcv`) is built by `protocols/hcv/build.sh`.

- **Panel:** the 238 references of h2seq (`assets/reference_data/hcv_references.fasta`); every one is kept. This is a
  judgement, not a measured win: in a test with unseen 1a, 1b, 3g and 4k genomes, the full panel and one representative
  per subtype covered about the same (within ±3 percentage points, better for some genomes and worse for others) and
  both assigned the right subtype every time. The reasons to keep everything are that alignment and lifting are one-off
  costs, nothing has to be chosen as "the" representative, the panel stays identical to the one h2seq uses, and the
  diversity is there if you later want finer calls. The main argument for representatives, near-identical siblings
  splitting reads and making the displayed reference unstable, is handled by the group lock above.
- **Coordinate reference:** H77, `1a_AF009606.1`, so positions are H77 numbers (E1E2 = 915-2579).
- **Genes:** the mature proteins (Core ... NS5B) from the standard polyprotein cleavage sites, checked against the GenBank
  translation (`genes.gff`).
- **Primers:** `H2Seq_HCVtiling_primerseq.fa`, 40 primers in 10 segments (`Seg1` ... `Seg10`); the amplicons are named
  `Seg1` ... `Seg10`. Some primers are genotype-specific, so a genotype the scheme was not designed for has fewer
  segments; the heatmap shows those as grey "no primer site" cells, and `qc.txt` lists which groups lack which segment.
  Segment positions on H77 run from 281 to 9,304.

## Known limitations

- Coverage comes from minimap2 alignment spans (deletions count as covered). It is a monitoring estimate; the depth from
  your final analysis is the arbiter.
- Reads are mapped to the panel, not to one chosen reference per sample. A read credited to a sibling reference counts
  towards that reference's coverage, so for similar genotypes `best_reference` can be any of the near-equal siblings.
- A protocol does not move depth between references: a sample is drawn from the depth on its one displayed reference,
  and positions that reference lacks (insertions or deletions against the coordinate reference) show as gaps.
- Reads are mapped to the whole panel, so a read that fits several references equally well is credited to each (see above).
- Read quality is not reported (only lengths); that would need the FASTQs to be read.
- Run time is measured from the first read of any listed barcode, not from the start of the run.
- The dashboard has no authentication; keep it on `127.0.0.1` (the default) unless you are on a trusted network.

## Development

```bash
pixi run -e dev pytest
bash protocols/hcv/build.sh                           # rebuild the bundled HCV protocol (about 80 seconds)
python dev/simulate_run.py --panel panel.fasta --outdir OUT_DIR \
    --barcode barcode01=3a_X76918.1:1500 --barcode barcode02=1a_AF009606.1:1800   # simulated reads, one genotype per barcode
```

### Trying it without a sequencer

**Want to see it work first?** `examples/hcv-e1e2-amplicons/` is a small simulated run (1.5 MB; two 1a and two 3a
whole-genome amplicon samples and a negative control) with the commands to monitor it, aiming for 20x over E1E2, and
what you should see.

`dev/simulate_run.py` makes reads from panel references (`REF[+REF2]:COUNT`). Add `@START-END` to draw reads only from a
window, or `@amplicons` (with `--protocol`) for tiling-amplicon reads that start and end at the reference's primer sites,
in uneven amounts per amplicon. `--divergence 0.03` gives each barcode its own genome, 3% away from its panel reference,
as a real sample would be. `dev/replay_run.py` plays any finished run back as if it were arriving
now, so you can watch the monitor work:

```bash
# terminal 1: start the monitor first; it waits for the folder to appear
nanosentinel monitor --input samplesheet.csv --run_dir replay/fastq_pass --protocol hcv
# terminal 2: one hour of sequencing in a minute
python dev/replay_run.py /path/to/finished/fastq_pass replay/fastq_pass --speed 60
```

- Chunks are released in the order the reads were produced, using the reads' own `start_time=` stamps (not file
  modification times, which a copy loses). A chunk appears when its last read was produced, as MinKNOW would write it.
- Each chunk is copied under a hidden name and then renamed, so the monitor never sees a half-copied file.
- `--speed 0` copies everything at once; `--dry_run` prints the schedule.
- Time to target is computed from the read timestamps, so it does not depend on the replay speed (checked: the same
  times at 120x as when the whole run is already on disk).
- To look at a finished run all at once, no replay is needed: point `--run_dir` at it.

To validate the tool against your analysis pipeline: replay a real run, run `nanosentinel monitor ... --concatenate_reads`,
put the concatenated reads through the pipeline, and compare its depth with the "covered at >=Nx" reported here.
