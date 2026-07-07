#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="$(pwd)"

REF="$BASE_DIR/S288C_Reference/SacCer3.fa"

VAR_BAM_DIR="$BASE_DIR/variant_bams"

CALLABLE_DIR="$BASE_DIR/callable_loci"
CALLABLE_BED_DIR="$CALLABLE_DIR/bed"

MPILEUP_DIR="$BASE_DIR/mpileup_callable_loci"
CALLABLE_ONLY_BED_DIR="$MPILEUP_DIR/callable_only_bed"
LOG_DIR="$BASE_DIR/logs_mpileup_callable_loci"

JOBS=9

# match CallableLoci settings
MIN_BASE_QUAL=10
MIN_MAPQ=10
MAX_DEPTH=1000000

mkdir -p "$MPILEUP_DIR" "$CALLABLE_ONLY_BED_DIR" "$LOG_DIR"

echo "=== Checking inputs ==="

[[ -s "$REF" ]] || {
    echo "ERROR: reference not found: $REF" >&2
    exit 1
}

[[ -d "$VAR_BAM_DIR" ]] || {
    echo "ERROR: BAM dir not found: $VAR_BAM_DIR" >&2
    exit 1
}

[[ -d "$CALLABLE_BED_DIR" ]] || {
    echo "ERROR: CallableLoci BED dir not found: $CALLABLE_BED_DIR" >&2
    exit 1
}

command -v samtools >/dev/null 2>&1 || {
    echo "ERROR: samtools not found in PATH" >&2
    exit 1
}

command -v bgzip >/dev/null 2>&1 || {
    echo "ERROR: bgzip not found in PATH" >&2
    exit 1
}

run_mpileup_callable_one() {
    set -euo pipefail

    local bam="$1"
    local sample
    local callable_bed
    local callable_only_bed
    local out

    sample="$(basename "$bam" .splitncigar.bam)"

    callable_bed="$CALLABLE_BED_DIR/${sample}.callable.bed"
    callable_only_bed="$CALLABLE_ONLY_BED_DIR/${sample}.CALLABLE_ONLY.bed"
    out="$MPILEUP_DIR/${sample}.callable.mpileup.gz"

    echo "=== Processing $sample ==="

    [[ -s "$bam" ]] || {
        echo "ERROR: BAM not found or empty: $bam" >&2
        exit 1
    }

    [[ -s "$callable_bed" ]] || {
        echo "ERROR: CallableLoci BED not found for $sample: $callable_bed" >&2
        exit 1
    }

    # GATK CallableLoci BED usually has state in column 4.
    # Keep only CALLABLE intervals.
    awk 'BEGIN{OFS="\t"}
        $0 !~ /^#/ && NF >= 4 && $4 == "CALLABLE" {
            print $1, $2, $3
        }
    ' "$callable_bed" > "$callable_only_bed"

    [[ -s "$callable_only_bed" ]] || {
        echo "ERROR: no CALLABLE intervals extracted for $sample from $callable_bed" >&2
        echo "First lines of BED:" >&2
        head "$callable_bed" >&2
        exit 1
    }

    echo "Callable intervals:"
    wc -l "$callable_only_bed"

    echo "Running mpileup..."

    samtools mpileup \
        -aa \
        -f "$REF" \
        -l "$callable_only_bed" \
        -q "$MIN_MAPQ" \
        -Q "$MIN_BASE_QUAL" \
        --ff 0x900 \
        -d "$MAX_DEPTH" \
        "$bam" \
    | bgzip -c > "$out"

    [[ -s "$out" ]] || {
        echo "ERROR: mpileup failed for $sample" >&2
        exit 1
    }

    echo "=== Done $sample ==="
    echo "Output: $out"
}

export REF VAR_BAM_DIR CALLABLE_BED_DIR MPILEUP_DIR CALLABLE_ONLY_BED_DIR
export MIN_BASE_QUAL MIN_MAPQ MAX_DEPTH
export -f run_mpileup_callable_one

find "$VAR_BAM_DIR" -maxdepth 1 -type f -name "WRS*.splitncigar.bam" -print0 | \
xargs -0 -n 1 -P "$JOBS" bash -euo pipefail -c '
        bam="$1"
        sample="$(basename "$bam" .splitncigar.bam)"
        run_mpileup_callable_one "$bam" > "'"$LOG_DIR"'/${sample}.mpileup_callable.log" 2>&1
    ' _

echo "=== Done ==="
echo "Callable-only BEDs: $CALLABLE_ONLY_BED_DIR"
echo "mpileup outputs:    $MPILEUP_DIR"
echo "logs:               $LOG_DIR"
