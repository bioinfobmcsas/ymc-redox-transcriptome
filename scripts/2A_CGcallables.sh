#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="$(pwd)"


REF="$BASE_DIR/S288C_Reference/S288C_reference_sequence_R64-5-1_20240529_v2.fa"
BED_DIR="$BASE_DIR/callable_loci/bed"
TMP_DIR="$BASE_DIR/callable_loci/tmp_cg_count"
OUT="$BASE_DIR/callable_loci/callable_base_counts.tsv"
INPUTS_DIR="$BASE_DIR/inputs"
INPUTS_OUT="$INPUTS_DIR/callable_base_counts.tsv"

mkdir -p "$TMP_DIR" "$INPUTS_DIR"

if [[ ! -f "$REF" ]]; then
    echo "ERROR: reference fasta not found: $REF" >&2
    exit 1
fi

if [[ ! -d "$BED_DIR" ]]; then
    echo "ERROR: BED directory not found: $BED_DIR" >&2
    exit 1
fi

shopt -s nullglob
beds=("$BED_DIR"/*.callable.bed)

if [[ ${#beds[@]} -eq 0 ]]; then
    echo "ERROR: no .callable.bed files found in $BED_DIR" >&2
    exit 1
fi

printf "sample\tA\tC\tG\tT\tN\tCG_total\tcallable_total\n" > "$OUT"

for bed in "${beds[@]}"; do
    sample="$(basename "$bed" .callable.bed)"
    pass_bed="$TMP_DIR/${sample}.PASS.bed"
    pass_fa="$TMP_DIR/${sample}.PASS.fa"

    echo "Processing $sample..."

    awk '$4=="CALLABLE"' "$bed" > "$pass_bed"

    if [[ ! -s "$pass_bed" ]]; then
        printf "%s\t0\t0\t0\t0\t0\t0\t0\n" "$sample" >> "$OUT"
        continue
    fi

    bedtools getfasta \
        -fi "$REF" \
        -bed "$pass_bed" \
        -fo "$pass_fa" >/dev/null 2>&1

    counts=$(grep -v '^>' "$pass_fa" | tr -d '\n' | awk '
    BEGIN {a=0; c=0; g=0; t=0; n=0}
    {
        seq=toupper($0)
        for(i=1;i<=length(seq);i++){
            b=substr(seq,i,1)
            if(b=="A") a++
            else if(b=="C") c++
            else if(b=="G") g++
            else if(b=="T") t++
            else n++
        }
    }
    END {
        callable_total = a + c + g + t + n
        cg_total = c + g
        printf "%d\t%d\t%d\t%d\t%d\t%d\t%d\n", a, c, g, t, n, cg_total, callable_total
    }')

    printf "%s\t%s\n" "$sample" "$counts" >> "$OUT"
done

cp "$OUT" "$INPUTS_OUT"

echo "Done."
echo "Output: $OUT"
echo "Analysis input: $INPUTS_OUT"
