# Roadmap

Planned work, roughly in the order we intend to do it. Items 1 and 3 are done; the rest is not started.

Sequencing always writes `fastq(.gz)` and `pod5`; BAM output is out of scope and `pod5` files are ignored.

## 1. Replay a run (also our test harness) (done)

`dev/replay_run.py SRC_fastq_pass DEST_fastq_pass [--speed N] [--dry_run]`; see the README ("Trying it without a
sequencer"). Differences from the original sketch, on purpose:

- A chunk is released when its **last** read was produced (that is when MinKNOW would have finished writing it), with the
  first read as the tie-break; the order still comes from the reads' `start_time=`, never from file modification times.
- Chunks are copied under a hidden name and renamed, so a half-copied chunk is never visible.
- A chunk without timestamps is placed after the previous chunk of its barcode, with a warning.

Still to do: replay a real MinKNOW run (item 2) and run the validation recipe in the README against the analysis
pipeline. `dev/simulate_run.py` (with `@START-END` windows) covers synthetic runs.

## 2. Check on real data (partly done)

Done: a finished real run on linux, monitored with the `hcv` protocol and an E1E2 target (real file layout and read
headers, end to end).

Still to check:

- A **live** MinKNOW run, as opposed to a finished one: chunks arriving while the run is in progress, and a run that
  has not started when the monitor does. (`dev/replay_run.py` can stand in until then, using the real dataset.)
- Annotation throughput at real chunk sizes and output rates, on the sequencing machine (and a large panel).
- A restart mid-run reuses annotations and resumes cleanly.
- Compare the "covered at >=Nx" it reports against the depth from the analysis pipeline for the same samples (the
  validation recipe in the README).

## 3. Subcommands, protocol prep and per-reference annotation (done)

`prep-protocol` and `monitor` are subcommands; a protocol (`protocols/<name>/`, `--protocol <path or bundled name>`)
holds the aligned, annotated panel. The README's "Protocols" section describes it. Where it differs from the original
sketch, on purpose:

- **Coordinate maps, not lifted features.** The panel is aligned to the coordinate reference with `mafft --add` (no
  `--keeplength`, so no base is lost and every reference base maps exactly). Each reference gets a map to coordinate
  positions; depth is projected onto the coordinate axis, so every sample is drawn with the same genes, regions and
  amplicon columns. `features.tsv` (genes and regions lifted onto each reference, with a status) is a QC output, not
  something `monitor` reads, so there is no hand-correcting of a lift; supply a better `--panel_alignment` instead.
- **Primers by search.** `--primer_fasta` places each primer on every reference (IUPAC-aware, mismatch limits, 3' end
  rule, constrained to the primer's consensus position); an amplicon spans its best-matching LEFT to RIGHT primer.
  `--primer_bed` (on the coordinate reference) is also supported and carried through the maps.
- **Taste then lock** is by sibling group (for example subtype), not by reference.
- **The whole panel is kept** (the HCV protocol has all 238); a one-per-subtype panel gave similar coverage.
- A protocol is checked by hash on load, and `--protocol` takes a bundled name (`hcv`) as well as a path, so a small
  version of the "named protocol library" (item 6) exists already.

Still open from the sketch:

- Regions are drawn and lifted, and one region can be the target (`--target_region`); several targets at once is item 4.
- A reference with a poor alignment is left off the shared axis and plotted on its own axis without annotations; there
  is no per-protocol tolerance setting for that.
- Primer mismatch limits are defaults to be tuned against the real scheme (`--max_mismatch`, `--max_3prime_mismatch`).

## 4. Unified targets (replaces the single depth target)

Depends on item 3, because regions must exist on every reference.

Today there is one target: `--target_depth` over `--target_fraction` of the genome, or of one named region
(`--target_region`) or of the amplicons (`--denominator amplicons`). Replace it with a list of named targets, each
tracked, reported and notified separately. The single-region case is done and is the building block: the tracker already
measures a target over a per-reference mask with a fixed denominator.

- `--targets targets.tsv` with columns `name`, `region`, `depth`, `fraction`, `reads`. The existing flags remain as a
  shortcut for a single whole-genome target.
- **Read-count targets:** track progress towards a target number of reads, with percentages accordingly. Count mapped
  reads by default (decide whether to count only reads inside the target region).
- **Region targets:** several at once (for example E1E2 and the whole genome), each with its own result.
- Each target gets its own column and bar in the table, its own time-to-target, and its own notification.

## 5. Decision support

- **Time-to-target estimate** per sample (for example "about 35 min at the current rate"), extrapolated from the
  coverage curve.
- **Plateau and stall warnings:** coverage no longer rising despite more reads (typically amplicon dropout); a barcode
  whose read rate has dried up. Add a "dropped amplicons" count per sample to the table.
- **Run-level QC panel:** yield over time (reads per hour), barcode balance (to catch unequal pooling), and the fraction
  of reads in `unclassified`, which is currently ignored.

## 6. Quality of life

- `run_summary.json` per run, machine-readable (for a LIMS or database).
- Optional filtering at concatenation: minimum and maximum length (and quality). Byte-concatenation cannot filter, so
  this means decompressing; keep it opt-in.
- Per-sample protocol or reference (a samplesheet column), for runs that mix pathogens.
- A named protocol library (`--protocol hcv-e1e2` resolving under a standard directory).
- Optional access token when the dashboard is served with `--host 0.0.0.0`; an SSH-tunnel recipe in the README.
- A small CI workflow that runs the tests and a tiny replay.

## 7. Optional RAMPART viewer (last)

- `monitor --use_rampart` runs RAMPART as an alternative viewer; the built-in dashboard stays the default.
- RAMPART, Snakemake, Biopython and Node return as an **optional pixi environment**, so the default install stays light.
- A small adapter turns RAMPART's annotation CSVs into the same per-read table the tracker reads.
- Left last: the most maintenance for the least new capability. The earlier RAMPART-based implementation is in the
  first commit of this repository.

## Open questions

- With `monitor` as a subcommand, where does concatenating the reads of a finished run live: keep
  `monitor --skip_monitor --concatenate_reads`, or give it its own subcommand (for example `concat`)?
- Should time-to-target be measured on the locked reference (cleaner, but depends on the early choice) or stay as "reached
  on any reference"?

## Not planned

- Sample-swap detection (genotype is unknown before sequencing).
- A warning for negative-control read counts during the run (only a control reaching target depth warns).
- BAM input and adaptive sampling (no barcodes).
