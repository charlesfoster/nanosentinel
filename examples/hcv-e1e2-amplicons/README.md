# Example: HCV whole-genome amplicons, aiming for 20x over E1E2

A small simulated MinKNOW run (about 1.5 MB, 6,020 reads) to try nanosentinel without a sequencer. It mimics
tiling-amplicon sequencing of the HCV genome with the bundled `hcv` protocol:

| barcode | sample | what it is |
| --- | --- | --- |
| barcode01 | `1a_patient_A` | subtype 1a (near `1a_EF407457.1`), 1,500 reads |
| barcode02 | `1a_patient_B` | subtype 1a (near `1a_M62321.1`), 1,500 reads |
| barcode03 | `3a_patient_A` | subtype 3a (near `3a_D17763.1`), 1,500 reads |
| barcode04 | `3a_patient_B` | subtype 3a (near `3a_JN714194.1`), 1,500 reads |
| barcode05 | `NTC` | negative control, 20 stray reads |

The reads are **simulated, not real**. Each read starts and ends at an amplicon's primer sites (Seg1 ... Seg10 of the
lab's tiling scheme), amplicons are drawn in uneven amounts as in a real pool, errors are about 5%, and each sample's
genome is 3% away from its nearest panel reference, so nothing is an exact copy of a panel sequence. Base qualities are
flat (one value per read) to keep the files small; nanosentinel does not read them.

The goal: **20x over 90% of E1E2** (H77 positions 915-2579), which is covered by Seg1, Seg2 and Seg3.

## Run it

From the nanosentinel folder (use `./bin/nanosentinel` if it is not on your PATH):

```bash
nanosentinel monitor \
    --input examples/hcv-e1e2-amplicons/samplesheet.csv \
    --run_dir examples/hcv-e1e2-amplicons/fastq_pass \
    --protocol hcv --target_region E1E2 --target_depth 20 --target_fraction 0.9 \
    --outdir examples/hcv-e1e2-amplicons/out
```

The dashboard opens at <http://localhost:6010>. The whole run is already on disk, so the monitor reads it in one go
(a few seconds); press **Enter** or **Ctrl+C** when it says all samples have reached the target. Then look in
`examples/hcv-e1e2-amplicons/out/`: `dashboard.html` (open it straight from disk), `time_to_target.tsv`,
`sample_summary.tsv` and the other reports.

### Watch it arrive as if live

Start the monitor on a folder that does not exist yet, then replay the run into it (two terminals):

```bash
# terminal 1
nanosentinel monitor \
    --input examples/hcv-e1e2-amplicons/samplesheet.csv \
    --run_dir examples/hcv-e1e2-amplicons/replay/fastq_pass \
    --protocol hcv --target_region E1E2 --target_depth 20 --target_fraction 0.9 \
    --outdir examples/hcv-e1e2-amplicons/out

# terminal 2: 1 hour 40 minutes of sequencing in about 100 seconds
python3 dev/replay_run.py examples/hcv-e1e2-amplicons/fastq_pass examples/hcv-e1e2-amplicons/replay/fastq_pass --speed 60
```

### Compare with the whole-genome target

Leave out `--target_region E1E2` (and use a different `--outdir`) to measure 20x over 90% of the whole genome instead.

## What you should see

Time to target, from the read timestamps (so it is the same whether you replay or not):

| sample | E1E2 target | reads needed | whole-genome target | reads needed |
| --- | --- | --- | --- | --- |
| `1a_patient_A` | 00:36:52 | 585 | 00:51:40 | 766 |
| `1a_patient_B` | 00:53:32 | 803 | 00:50:20 | 758 |
| `3a_patient_A` | 00:13:46 | 181 | 00:42:28 | 627 |
| `3a_patient_B` | 00:12:35 | 202 | 00:24:37 | 382 |
| `NTC` | control, never reaches it | | | |

(A different minimap2 version can move these by a little.)

Things worth noticing:

- **One genome, one axis.** The 1a and 3a samples are mapped to their own references but drawn on the same H77 axis, so the
  genes, the E1E2 bracket and the amplicon bands line up for everyone.
- **A region target can be stricter than the genome.** E1E2 is only 1,665 bases and 90% of it leaves just 167 bases of
  slack, while Seg1, Seg2 and Seg3 each cover about a quarter of it that no other amplicon covers. In `1a_patient_B` Seg3
  gets few reads, so E1E2 is reached *after* the whole genome. In the amplicon heatmap, the Seg1-Seg3 cells are the ones
  to watch when E1E2 is the aim.
- **The negative control stays flat** (about 2x), and is flagged as a control, not a sample.

## Regenerate the data

`bash examples/hcv-e1e2-amplicons/make_data.sh` rewrites `fastq_pass/` (needs only `python3`; the result is deterministic).
The command inside is a good template for other scenarios: see `dev/simulate_run.py --help`.
