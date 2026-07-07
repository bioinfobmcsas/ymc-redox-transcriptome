#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="$(pwd)"

FILTERED_DIR="$BASE_DIR/mutect2/filtered_vcfs"
OUT_DIR="$BASE_DIR/mutect2/merged_sites"
PASS_VCF_DIR="$BASE_DIR/mutect2/pass_tlod_vcfs"

MIN_SAMPLES=1
MIN_TLOD=6

mkdir -p "$OUT_DIR" "$PASS_VCF_DIR"

TMP="$OUT_DIR/all_pass_snps_tlod.tmp.tsv"
TMP_BODY="$OUT_DIR/all_pass_snps_recurrent.body.tsv"
OUT_TSV="$OUT_DIR/all_pass_snps_recurrent.tsv"
OUT_BED="$OUT_DIR/all_pass_snps_recurrent.bed"

: > "$TMP"
: > "$TMP_BODY"

echo "=== Filtering VCFs (PASS + SNP + TLOD>$MIN_TLOD) ==="

shopt -s nullglob
vcfs=("$FILTERED_DIR"/*.mutect2.filtered.vcf.gz)

if [[ ${#vcfs[@]} -eq 0 ]]; then
    echo "ERROR: no *.mutect2.filtered.vcf.gz files found in $FILTERED_DIR" >&2
    exit 1
fi

for vcf in "${vcfs[@]}"; do
    sample="$(basename "$vcf" .mutect2.filtered.vcf.gz)"

    # фильтр только WRS
    if [[ ! "$sample" =~ ^WRS ]]; then
        continue
    fi

    out_vcf="$PASS_VCF_DIR/${sample}.pass.tlod.snps.vcf.gz"

    echo "Processing $sample"

    bcftools view \
        -i "FILTER='PASS' && INFO/TLOD>${MIN_TLOD}" \
        -m2 -M2 \
        -v snps \
        -Oz -o "$out_vcf" \
        "$vcf"

    bcftools index -t "$out_vcf"

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$out_vcf" >> "$TMP"
done

echo "=== Building recurrent sites (>= $MIN_SAMPLES samples) ==="

awk -v min="$MIN_SAMPLES" '
BEGIN { OFS="\t" }
{
    key = $1"\t"$2"\t"$3"\t"$4
    count[key]++
}
END {
    for (k in count) {
        if (count[k] >= min) {
            split(k, a, "\t")
            variant_id = a[1]"_"a[2]"_"a[3]"_"a[4]
            print variant_id, a[1], a[2], a[3], a[4], count[k]
        }
    }
}
' "$TMP" | sort -k2,2 -k3,3n > "$TMP_BODY"

echo -e "variant_id\tCHROM\tPOS\tREF\tALT\tn_samples_pass" > "$OUT_TSV"
cat "$TMP_BODY" >> "$OUT_TSV"

awk 'BEGIN { OFS = "\t" } NR > 1 { print $2, $3 - 1, $3 }' "$OUT_TSV" > "$OUT_BED"

rm -f "$TMP" "$TMP_BODY"

echo "=== Done ==="
echo "Filtered VCFs: $PASS_VCF_DIR"
echo "Final variants: $OUT_TSV"
echo "BED: $OUT_BED"
