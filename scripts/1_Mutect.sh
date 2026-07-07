#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="$(pwd)"

REF="$BASE_DIR/S288C_Reference/SacCer3.fa"

IN_DIR="$BASE_DIR/all_bams"
VAR_BAM_DIR="$BASE_DIR/variant_bams"

MUTECT_DIR="$BASE_DIR/mutect2"
RAW_DIR="$MUTECT_DIR/raw_vcfs"
FILTERED_DIR="$MUTECT_DIR/filtered_vcfs"

METRICS_DIR="$BASE_DIR/metrics"
TMP_DIR="$BASE_DIR/tmp"
LOG_DIR="$BASE_DIR/logs_mutect2_preprocessing"

JOBS=9
SAMTOOLS_THREADS=1
JAVA_MEM="8g"

mkdir -p \
    "$VAR_BAM_DIR" \
    "$RAW_DIR" \
    "$FILTERED_DIR" \
    "$METRICS_DIR" \
    "$TMP_DIR" \
    "$LOG_DIR"

[[ -s "$REF" ]] || {
    echo "ERROR: reference FASTA not found: $REF" >&2
    exit 1
}

if [[ ! -s "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi

REF_DICT="${REF%.*}.dict"

if [[ ! -s "$REF_DICT" ]]; then
    gatk CreateSequenceDictionary \
        -R "$REF" \
        -O "$REF_DICT"
fi

process_bam_mutect2() {
    set -euo pipefail

    local bam="$1"
    local sample
    local rg_bam
    local markdup_bam
    local split_bam
    local metrics
    local raw_vcf
    local filtered_vcf
    local stats
    local sample_tmp

    sample="$(basename "$bam" .bam)"

    rg_bam="$VAR_BAM_DIR/${sample}.rg.bam"
    markdup_bam="$VAR_BAM_DIR/${sample}.markdup.bam"
    split_bam="$VAR_BAM_DIR/${sample}.splitncigar.bam"

    metrics="$METRICS_DIR/${sample}.markdup.metrics.txt"

    raw_vcf="$RAW_DIR/${sample}.mutect2.vcf.gz"
    filtered_vcf="$FILTERED_DIR/${sample}.mutect2.filtered.vcf.gz"
    stats="${raw_vcf}.stats"

    sample_tmp="$TMP_DIR/${sample}_mutect2"
    mkdir -p "$sample_tmp"

    echo "=== Processing $sample ==="

    echo "[1/5] AddOrReplaceReadGroups"
    gatk --java-options "-Xmx${JAVA_MEM} -Djava.io.tmpdir=${sample_tmp}" \
        AddOrReplaceReadGroups \
        -I "$bam" \
        -O "$rg_bam" \
        -RGID "$sample" \
        -RGLB "$sample" \
        -RGPL ILLUMINA \
        -RGPU "$sample" \
        -RGSM "$sample"

    [[ -s "$rg_bam" ]] || {
        echo "ERROR: AddOrReplaceReadGroups failed for $sample" >&2
        exit 1
    }

    samtools index -@ "$SAMTOOLS_THREADS" "$rg_bam"

    echo "[2/5] MarkDuplicates"
    gatk --java-options "-Xmx${JAVA_MEM} -Djava.io.tmpdir=${sample_tmp}" \
        MarkDuplicates \
        -I "$rg_bam" \
        -O "$markdup_bam" \
        -M "$metrics" \
        --CREATE_INDEX false

    [[ -s "$markdup_bam" ]] || {
        echo "ERROR: MarkDuplicates failed for $sample" >&2
        exit 1
    }

    samtools index -@ "$SAMTOOLS_THREADS" "$markdup_bam"

    echo "[3/5] SplitNCigarReads"
    gatk --java-options "-Xmx${JAVA_MEM} -Djava.io.tmpdir=${sample_tmp}" \
        SplitNCigarReads \
        -R "$REF" \
        -I "$markdup_bam" \
        -O "$split_bam"

    [[ -s "$split_bam" ]] || {
        echo "ERROR: SplitNCigarReads failed for $sample" >&2
        exit 1
    }

    samtools index -@ "$SAMTOOLS_THREADS" "$split_bam"

    echo "[4/5] Mutect2"
    gatk --java-options "-Xmx${JAVA_MEM} -Djava.io.tmpdir=${sample_tmp}" \
        Mutect2 \
        -R "$REF" \
        -I "$split_bam" \
        -O "$raw_vcf" \
        --dont-use-soft-clipped-bases true

    [[ -s "$raw_vcf" ]] || {
        echo "ERROR: Mutect2 failed for $sample" >&2
        exit 1
    }

    [[ -s "$stats" ]] || {
        echo "ERROR: missing Mutect2 stats file for $sample" >&2
        exit 1
    }

    echo "[5/5] FilterMutectCalls"
    gatk --java-options "-Xmx${JAVA_MEM} -Djava.io.tmpdir=${sample_tmp}" \
        FilterMutectCalls \
        -R "$REF" \
        -V "$raw_vcf" \
        --stats "$stats" \
        -O "$filtered_vcf"

    [[ -s "$filtered_vcf" ]] || {
        echo "ERROR: FilterMutectCalls failed for $sample" >&2
        exit 1
    }

    rm -rf "$sample_tmp"

    echo "=== Done $sample ==="
}

export REF
export VAR_BAM_DIR RAW_DIR FILTERED_DIR METRICS_DIR TMP_DIR LOG_DIR
export SAMTOOLS_THREADS JAVA_MEM
export -f process_bam_mutect2

find "$IN_DIR" -maxdepth 1 -type f -name "*.bam" -print0 | \
    xargs -0 -n 1 -P "$JOBS" bash -euo pipefail -c '
        bam="$1"
        sample="$(basename "$bam" .bam)"
        process_bam_mutect2 "$bam" > "$LOG_DIR/${sample}.log" 2>&1
    ' _

echo "=== Building union of PASS sites across samples ==="

PASS_VCFS=()

while IFS= read -r -d '' f; do
    PASS_VCFS+=("$f")
done < <(find "$FILTERED_DIR" -maxdepth 1 -type f -name "*.mutect2.filtered.vcf.gz" -print0 | sort -z)

if [[ ${#PASS_VCFS[@]} -eq 0 ]]; then
    echo "ERROR: no filtered VCF files found" >&2
    exit 1
fi

tmp_sites="$MUTECT_DIR/all_pass_sites.tmp.tsv"

: > "$tmp_sites"

for vcf in "${PASS_VCFS[@]}"; do
    bcftools view -f PASS "$vcf" | \
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' >> "$tmp_sites"
done

sort -u "$tmp_sites" > "$MUTECT_DIR/all_pass_sites.tsv"
rm -f "$tmp_sites"

cut -f1,2 "$MUTECT_DIR/all_pass_sites.tsv" | \
    sort -u > "$MUTECT_DIR/all_pass_sites.positions.tsv"

echo "=== Done ==="
echo "Variant BAMs:        $VAR_BAM_DIR"
echo "Raw Mutect2 VCFs:    $RAW_DIR"
echo "Filtered VCFs:       $FILTERED_DIR"
echo "Full site list:      $MUTECT_DIR/all_pass_sites.tsv"
echo "Positions only:      $MUTECT_DIR/all_pass_sites.positions.tsv"
