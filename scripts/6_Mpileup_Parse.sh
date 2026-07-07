#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="$(pwd)"

MPILEUP_DIR="$BASE_DIR/mpileup_callable_loci"
OUT_DIR="$BASE_DIR/inputs/loci_counts"
PARSER="$BASE_DIR/scripts/mpileup_to_ref_alt_csv.py"

JOBS=6

mkdir -p "$OUT_DIR"

[[ -d "$MPILEUP_DIR" ]] || { echo "ERROR: no dir $MPILEUP_DIR" >&2; exit 1; }
[[ -f "$PARSER" ]] || { echo "ERROR: no parser $PARSER" >&2; exit 1; }

process_one() {
    set -euo pipefail

    local mpileup="$1"
    local sample

    sample="$(basename "$mpileup" .callable.mpileup.gz)"

    echo "Processing $sample"

    gzip -dc "$mpileup" | python3 "$PARSER" | \
        gzip > "$OUT_DIR/${sample}.ref_alt_counts.csv.gz"

    echo "Done $sample"
}

export OUT_DIR PARSER
export -f process_one

find "$MPILEUP_DIR" \
    -maxdepth 1 \
    -type f \
    -name "WRS*.callable.mpileup.gz" | \
    sort | \
    xargs -n 1 -P "$JOBS" bash -euo pipefail -c '
        process_one "$1"
    ' _

echo "=== DONE ==="
echo "Output: $OUT_DIR"
