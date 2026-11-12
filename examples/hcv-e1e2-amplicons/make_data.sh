#!/usr/bin/env bash
# Regenerates fastq_pass/ (the data in this folder). Run from anywhere:  bash examples/hcv-e1e2-amplicons/make_data.sh
# Needs only python3. The result is deterministic (seeded), so you get the files that ship here.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$here/../.."
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
python3 "$root/dev/simulate_run.py" --panel "$root/protocols/hcv/panel.fasta" --protocol "$root/protocols/hcv" \
  --outdir "$tmp" --divergence 0.03 --flat_quality --chunk 250 \
  --barcode barcode01=1a_EF407457.1:1500@amplicons \
  --barcode barcode02=1a_M62321.1:1500@amplicons \
  --barcode barcode03=3a_D17763.1:1500@amplicons \
  --barcode barcode04=3a_JN714194.1:1500@amplicons \
  --barcode barcode05=3a_D17763.1:20@amplicons
rm -rf "$here/fastq_pass"
cp -R "$tmp"/exp1/sample1/*/fastq_pass "$here/fastq_pass"
echo "Wrote $here/fastq_pass"
