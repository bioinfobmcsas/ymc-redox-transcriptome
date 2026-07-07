#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="$(pwd)"

REF="$BASE_DIR/S288C_Reference/SacCer3.fa"
VAR_BAM_DIR="$BASE_DIR/variant_bams"
CALLABLE_DIR="$BASE_DIR/callable_loci"
BED_DIR="$CALLABLE_DIR/bed"
SUMMARY_DIR="$CALLABLE_DIR/summary"
TMP_DIR="$BASE_DIR/tmp"
LOG_DIR="$BASE_DIR/logs_callableloci"

JOBS=9
JAVA_MEM="8g"

MIN_DEPTH=10
MAX_DEPTH=1000000
MIN_BASE_QUAL=10
MIN_MAPQ=10
MAX_LOW_MAPQ=1
MAX_FRAC_LOW_MAPQ=0.1

mkdir -p "$BED_DIR" "$SUMMARY_DIR" "$TMP_DIR" "$LOG_DIR"

run_callableloci_one() {
    set -euo pipefail

    local bam="$1"
    local sample
    local out_bed
    local out_summary
    local sample_tmp

    sample="$(basename "$bam" .splitncigar.bam)"
    out_bed="$BED_DIR/${sample}.callable.bed"
    out_summary="$SUMMARY_DIR/${sample}.callable.summary.tsv"
    sample_tmp="$TMP_DIR/${sample}_callableloci"

    mkdir -p "$sample_tmp"

    echo "=== Processing $sample ==="

    gatk --java-options "-Xmx${JAVA_MEM} -Djava.io.tmpdir=${sample_tmp}" \
        CallableLoci \
        -R "$REF" \
        -I "$bam" \
        -O "$out_bed" \
        --summary "$out_summary" \
        --min-depth "$MIN_DEPTH" \
        --max-depth "$MAX_DEPTH" \
        --min-base-quality "$MIN_BASE_QUAL" \
        --min-mapping-quality "$MIN_MAPQ" \
        --max-low-mapq "$MAX_LOW_MAPQ" \
        --max-fraction-of-reads-with-low-mapq "$MAX_FRAC_LOW_MAPQ"

    [[ -s "$out_bed" ]] || {
        echo "ERROR: CallableLoci BED failed for $sample" >&2
        exit 1
    }

    [[ -s "$out_summary" ]] || {
        echo "ERROR: CallableLoci summary failed for $sample" >&2
        exit 1
    }

    rm -rf "$sample_tmp"

    echo "=== Done $sample ==="
}

export REF BED_DIR SUMMARY_DIR TMP_DIR LOG_DIR JAVA_MEM
export MIN_DEPTH MAX_DEPTH MIN_BASE_QUAL MIN_MAPQ MAX_LOW_MAPQ MAX_FRAC_LOW_MAPQ
export -f run_callableloci_one

find "$VAR_BAM_DIR" -maxdepth 1 -type f -name "*.splitncigar.bam" -print0 | \
    xargs -0 -n 1 -P "$JOBS" bash -euo pipefail -c '
        bam="$1"
        sample="$(basename "$bam" .splitncigar.bam)"
        run_callableloci_one "$bam" > "$LOG_DIR/${sample}.log" 2>&1
    ' _

echo "=== Building combined summary table ==="

COMBINED="$CALLABLE_DIR/callable_summary_all_samples.tsv"

{
    printf "sample\tstate\tnBases\n"

    for f in "$SUMMARY_DIR"/*.callable.summary.tsv; do
        sample="$(basename "$f" .callable.summary.tsv)"

        awk -v sample="$sample" '
            BEGIN { OFS="\t" }
            NR == 1 { next }
            NF >= 2 { print sample, $1, $2 }
        ' "$f"
    done
} > "$COMBINED"

echo "=== Done ==="
echo "BED intervals:    $BED_DIR"
echo "Per-sample stats: $SUMMARY_DIR"
echo "Combined table:   $COMBINED"
