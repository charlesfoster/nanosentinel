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
- Handles **several similar references** in the panel: a read that fits sibling genotypes equally well counts towards
  all of them, so coverage is not split between near-identical references.
- Waits automatically if MinKNOW has not created the run's `fastq_pass` folder yet.
- Quiet terminal: one status line and a few event messages. Everything is also logged to `nanosentinel.log`.

## Install

Needs [pixi](https://pixi.sh). The environment (Python, minimap2, numpy) installs itself on first run, and again if
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
nanosentinel --input samplesheet.csv --run_name 2026-11-09_run \
    --possible_references references.fasta --primer_bed primers.bed
```

Opens the dashboard at <http://localhost:6010>. If the run has not started yet it waits for `fastq_pass` to appear
(Ctrl+C to cancel). Press **Enter** or **Ctrl+C** to stop.

Inputs are ordinary files: `--possible_references` is a multi-FASTA panel the reads are mapped to (include one entry per
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

- `--denominator amplicons` measures the fraction over the primer-scheme span instead of the whole genome.
- `--target_depth`, `--target_fraction`, `--min_identity`, `--webhook URL`, `--port`, `--no_browser`.
- `--host 0.0.0.0` to view the dashboard from another machine on the network (there is no authentication).
- `--run_dir` points at an exact MinKNOW run folder (use it if a run name matches more than one folder).
- `--dry_run` checks the inputs and builds the reference index without monitoring.
- `--clear_annotations` discards earlier annotations and maps every read again.

## Concatenate reads for analysis

At the end of monitoring:

```bash
nanosentinel ... --concatenate_reads
```

or for a finished run, without monitoring:

```bash
nanosentinel --input samplesheet.csv --run_name 2026-11-09_run --skip_monitor --concatenate_reads
```

This writes `reads/<sample>_<run>_<barcode>.fastq.gz` and `samplesheet_with_reads.csv` (this run's rows, `long_reads`
filled in) ready to pass to an analysis pipeline. Gzipped FASTQs are joined without recompression; `--gzip` compresses
plain-text input.

## Known limitations

- Coverage comes from minimap2 alignment spans (deletions count as covered). It is a monitoring estimate; the depth from
  your final analysis is the arbiter.
- Reads are mapped to the panel, not to one chosen reference per sample. A read credited to a sibling reference counts
  towards that reference's coverage, so for similar genotypes `best_reference` can be any of the near-equal siblings.
- Read quality is not reported (only lengths); that would need the FASTQs to be read.
- Run time is measured from the first read of any listed barcode, not from the start of the run.
- The dashboard has no authentication; keep it on `127.0.0.1` (the default) unless you are on a trusted network.

## Development

```bash
pixi run -e dev pytest
python dev/make_fake_run.py reads.fastq.gz OUT_DIR   # a fake MinKNOW run to try it on (use --minknow_data OUT_DIR --run_name exp1)
```
