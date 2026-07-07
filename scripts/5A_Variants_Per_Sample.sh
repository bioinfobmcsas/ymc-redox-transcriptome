#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="$(pwd)"

PASS_VCF_DIR="$BASE_DIR/mutect2/pass_tlod_vcfs"
OUT_DIR="$BASE_DIR/mutect2/merged_sites"
INPUTS_DIR="$BASE_DIR/inputs"

RECURRENT_TSV="$OUT_DIR/all_pass_snps_recurrent.tsv"
OUT_TSV="$OUT_DIR/all_pass_snps_per_sample.tsv"
INPUTS_TSV="$INPUTS_DIR/all_pass_snps_per_sample.tsv"

mkdir -p "$OUT_DIR" "$INPUTS_DIR"

tmp_sample_variants="$(mktemp)"
tmp_recurrent="$(mktemp)"

trap 'rm -f "$tmp_sample_variants" "$tmp_recurrent"' EXIT

if [[ ! -s "$RECURRENT_TSV" ]]; then
    echo "ERROR: recurrent SNP table not found: $RECURRENT_TSV" >&2
    exit 1
fi

awk 'BEGIN { FS = OFS = "\t" } NR > 1 { print $1 }' "$RECURRENT_TSV" > "$tmp_recurrent"

shopt -s nullglob
vcfs=("$PASS_VCF_DIR"/WRS*.pass.tlod.snps.vcf.gz)

if [[ ${#vcfs[@]} -eq 0 ]]; then
    echo "ERROR: no WRS*.pass.tlod.snps.vcf.gz files found in $PASS_VCF_DIR" >&2
    exit 1
fi

for vcf in "${vcfs[@]}"; do
    sample="$(basename "$vcf" .pass.tlod.snps.vcf.gz)"
    echo "Processing $sample" >&2

    gzip -cd "$vcf" | awk -v sample="$sample" '
        BEGIN { OFS = "\t" }
        /^#/ { next }
        {
            variant_id = $1"_"$2"_"$4"_"$5
            print sample, variant_id, $1, $2, $4, $5
        }
    '
done > "$tmp_sample_variants"

{
    printf "sample\tvariant_id\tCHROM\tPOS\tREF\tALT\n"

    awk '
    BEGIN { FS = OFS = "\t" }
    NR == FNR {
        recurrent[$1] = 1
        next
    }
    ($2 in recurrent) {
        print
    }
    ' "$tmp_recurrent" "$tmp_sample_variants"
} > "$OUT_TSV"

cp "$OUT_TSV" "$INPUTS_TSV"

echo "Done: $OUT_TSV"
echo "Rows: $(( $(wc -l < "$OUT_TSV") - 1 ))"
echo "Analysis input: $INPUTS_TSV"
