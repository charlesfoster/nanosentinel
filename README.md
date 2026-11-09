# nanosentinel

Monitor a live Oxford Nanopore run with [RAMPART](https://github.com/artic-network/rampart), using a samplesheet in
the same format as h2seq's (so one sheet can serve both), and record how long each sample takes to reach a target depth.

- Generates the RAMPART protocol (references, genome, primer scheme, coverage thresholds) from plain FASTA / BED / GFF
  files, so there are no hand-maintained protocol folders. It is independent of h2seq; you can simply reuse the same
  reference and primer files.
- Quiet terminal: one status line, a few event messages. RAMPART's own output goes to `rampart.log`.
- Writes `time_to_target.tsv`: when each sample first reached `--target_depth` over `--target_fraction` of the
  genome (default 20x over 90%), measured in sequencing time from the reads' own timestamps.
- Notifies you when every non-control sample has reached target; monitoring continues until you stop it.
- Waits automatically if MinKNOW has not created the run's `fastq_pass` folder yet.
- Records progress over time and writes end-of-run summaries (read counts, read length statistics, per-reference composition).
- Optionally concatenates each sample's reads and writes an h2seq-ready samplesheet.

## Install

Needs [pixi](https://pixi.sh). The environment (RAMPART, minimap2, Snakemake, ...) installs itself on first run.

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

Opens RAMPART at <http://localhost:6010>. If the run has not started yet it waits for `fastq_pass` to appear (Ctrl+C to
cancel). Press **Enter** or **Ctrl+C** to stop.

Outputs in `nanosentinel_<run_name>/` (reports are refreshed about once a minute and written again when you stop):

| file | contents |
| --- | --- |
| `time_to_target.tsv` | When each sample first reached the target (sequencing time), with reads at that point and the reference it was reached on. |
| `coverage_timeseries.tsv` | Progress over time: one line per sample per annotation file with reads, mapped bases and the best fraction at target depth. Plot `elapsed_min` against `best_fraction`. |
| `sample_summary.tsv` | Per sample: reads, % mapped, dominant reference and its mean depth, and read length statistics (mean, median, N50, min, max) for all reads. |
| `reference_summary.tsv` | Per sample per reference that received reads: read share, estimated mean depth and fraction at target depth. Shows genotype composition and mixed samples. |
| `nanosentinel.log`, `rampart.log` | The tool's log and RAMPART's (very chatty) output. |
| `annotations/`, `rampart_protocol/` | RAMPART's per-file annotation CSVs and the protocol generated for this run. |

In `sample_summary.tsv`, `dominant_reference` is the reference with the most reads, while `best_fraction` is the highest
fraction at target depth over any reference; they can differ when a sample's reads split between similar references.

Useful options (see `--help` for all):

- `--reference_fasta` coordinate reference the coverage plot is drawn against (default: first panel entry); `--gff` adds genes.
- `--denominator amplicons` measures the fraction over the primer-scheme span instead of the whole genome.
- `--target_depth`, `--target_fraction`, `--webhook URL`, `--rampart_port`, `--no_browser`.
- `--run_dir` points at an exact MinKNOW run folder (use it if a run name matches more than one folder).
- `--dry_run` checks the inputs and writes the protocol without launching anything.

Restarting with the same `--outdir` reuses earlier annotations; add `--clear_annotations` to start afresh.

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
filled in) ready for `nextflow run charlesfoster/h2seq --input samplesheet_with_reads.csv`. Gzipped FASTQs are joined without
recompression; `--gzip` compresses plain-text input.

## Known limitations

- Coverage is estimated from RAMPART's per-read start/end coordinates (no CIGAR, no MAPQ), so it is slightly
  optimistic; h2seq's depth is the final word. Depth is tracked per sample against each reference in the panel and the
  best one is reported, so a sample whose reads split across near-identical references is under-counted.
- Read quality is not reported: RAMPART's annotations hold read lengths but not quality scores.
- RAMPART ignores reads under any hidden (dot-prefixed) folder, so nanosentinel refuses to start if the `fastq_pass` path contains one.
- Run time is measured from the first read of any listed barcode, not from the start of the run.
- RAMPART 1.2.0 ignores `--annotatedPath`; nanosentinel works around it by running RAMPART with the output directory as its working directory.

## Development

```bash
pixi run -e dev pytest
```
