#!/usr/bin/env bash
# Rebuilds this protocol. Run from the nanosentinel folder:  bash protocols/hcv/build.sh [path/to/hcv_references.fasta]
#
# Inputs
#   panel     h2seq/assets/reference_data/hcv_references.fasta (238 references; the id's prefix is the subtype)
#   primers   H2Seq_HCVtiling_primerseq.fa (10 segments, some primers genotype-specific; copied here as primers.fa)
#   genes     genes.gff (mature proteins of H77, from the polyprotein cleavage sites)
# The coordinate reference is H77 (1a_AF009606.1), so positions are H77 numbers. E1E2 is H77 915-2579.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
panel="${1:-$HOME/Programs/h2seq/assets/reference_data/hcv_references.fasta}"
"$here/../../bin/nanosentinel" prep-protocol \
  --possible_references "$panel" \
  --coordinate_ref 1a_AF009606.1 \
  --gff "$here/genes.gff" \
  --region E1E2=915-2579 \
  --primer_fasta "$here/primers.fa" \
  --amplicon_pattern '^(Seg[0-9]+)_' \
  --name hcv \
  --outdir "$here"
