#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="$(pwd)"

PASS_VCF_DIR="$BASE_DIR/mutect2/pass_tlod_vcfs"
OUT_DIR="$BASE_DIR/mutect2/merged_sites"
INPUTS_DIR="$BASE_DIR/inputs"

OUT_TSV="$OUT_DIR/all_pass_snps_presence_absence_matrix.tsv"
INPUTS_MATRIX="$INPUTS_DIR/all_pass_snps_presence_absence_matrix.tsv"

mkdir -p "$OUT_DIR" "$INPUTS_DIR"

tmp_sample_variants="$(mktemp)"
tmp_samples="$(mktemp)"

trap 'rm -f "$tmp_sample_variants" "$tmp_samples"' EXIT

shopt -s nullglob
vcfs=("$PASS_VCF_DIR"/WRS*.pass.tlod.snps.vcf.gz)

if [[ ${#vcfs[@]} -eq 0 ]]; then
    echo "ERROR: no WRS*.pass.tlod.snps.vcf.gz files found in $PASS_VCF_DIR" >&2
    exit 1
fi

for vcf in "${vcfs[@]}"; do
    sample="$(basename "$vcf" .pass.tlod.snps.vcf.gz)"
    echo "$sample" >> "$tmp_samples"

    echo "Processing $sample" >&2

    gzip -cd "$vcf" | awk -v sample="$sample" '
        BEGIN { OFS="\t" }
        /^#/ { next }
        {
            variant_id = $1"_"$2"_"$4"_"$5
            print sample, variant_id, $1, $2, $4, $5
        }
    '
done > "$tmp_sample_variants"

awk '
BEGIN {
    FS = OFS = "\t"
}

FILENAME == ARGV[1] {
    samples[++n_samples] = $1
    next
}

FILENAME == ARGV[2] {
    sample = $1
    variant_id = $2

    if (!(variant_id in seen_variant)) {
        seen_variant[variant_id] = 1
        variants[++n_variants] = variant_id

        chrom[variant_id] = $3
        pos[variant_id]   = $4
        ref[variant_id]   = $5
        alt[variant_id]   = $6
    }

    present[variant_id, sample] = 1
    next
}

END {
    printf "CHROM\tPOS\tREF\tALT\tvariant_id"
    for (i = 1; i <= n_samples; i++) {
        printf "\t%s", samples[i]
    }
    printf "\n"

    for (v = 1; v <= n_variants; v++) {
        variant_id = variants[v]

        printf "%s\t%s\t%s\t%s\t%s", \
            chrom[variant_id], \
            pos[variant_id], \
            ref[variant_id], \
            alt[variant_id], \
            variant_id

        for (i = 1; i <= n_samples; i++) {
            sample = samples[i]
            printf "\t%d", ((variant_id, sample) in present ? 1 : 0)
        }

        printf "\n"
    }
}
' "$tmp_samples" "$tmp_sample_variants" > "$OUT_TSV"

cp "$OUT_TSV" "$INPUTS_MATRIX"

echo "Done: $OUT_TSV"
echo "Variants: $(( $(wc -l < "$OUT_TSV") - 1 ))"
echo "Samples: ${#vcfs[@]}"
echo "Analysis matrix: $INPUTS_MATRIX"
