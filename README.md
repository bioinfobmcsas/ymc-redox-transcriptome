# ymc-redox-transcriptome
Analysis code for whole-cell RNA-seq data, including redox-phase transcriptome, RNA-seq mismatch, and variant-rate analyses across the yeast metabolic cycle.

## Installation

Create the conda environment from the provided `environment.yml` file:

```bash
conda env create -f environment.yml
conda activate ceres-rnaseq
```

## Usage

After activating the environment, run the figure-reproduction scripts in `scripts_analysis/` or execute the preprocessing workflows in `scripts/`.

Some large intermediate inputs, including `inputs/loci_counts/`, are expected to come from `dataset.zip` distributed via GitHub Releases rather than from the repository itself.

Expression and GO-analysis inputs used for Figures 2-5 are provided in `inputs/`: `wrs3_s1_counts.txt`, `deseq2_LvsH_wrs3_s1_genes.xlsx`, and `sgd.gaf.gz`.

`scripts_analysis/deseq2_wrs3_gene.r` runs the paired Low-DO versus High-DO DESeq2 analysis from `wrs3_s1_counts.txt` and writes `deseq2_LvsH_wrs3_s1_genes.xlsx` with GAF-compatible gene symbols.

## Pipeline Outputs Used As Inputs

The preprocessing scripts in `scripts/` produce the intermediate files consumed by the figure-reproduction scripts:

| Step | Script | Main pipeline output | Analysis-ready input |
| --- | --- | --- | --- |
| Callable base counts | `scripts/2A_CGcallables.sh` | `callable_loci/callable_base_counts.tsv` | `inputs/callable_base_counts.tsv` |
| PASS SNV filtering | `scripts/3_SNPs.sh` | `mutect2/pass_tlod_vcfs/*.pass.tlod.snps.vcf.gz`, `mutect2/merged_sites/all_pass_snps_recurrent.tsv`, `mutect2/merged_sites/all_pass_snps_recurrent.bed` | Used by the next workflow step |
| Variant presence matrix | `scripts/5_Variants_Presence.sh` | `mutect2/merged_sites/all_pass_snps_presence_absence_matrix.tsv` | `inputs/all_pass_snps_presence_absence_matrix.tsv` |
| Per-sample variant table | `scripts/5A_Variants_Per_Sample.sh` | `mutect2/merged_sites/all_pass_snps_per_sample.tsv` | `inputs/all_pass_snps_per_sample.tsv` |
| Callable-locus ref/alt counts | `scripts/6_Mpileup_Parse.sh` | `inputs/loci_counts/*.ref_alt_counts.csv.gz` | Used by `scripts/callable_postprocessing.R` |
| Locus-level mismatch summaries | `scripts/callable_postprocessing.R` | `inputs/loci_counts_summed_by_site_cond.csv.gz`, `inputs/loci_delta_by_site.csv.gz` | `inputs/loci_delta.csv`, `inputs/loci_per_sample_alt_fraction.csv` |

# Figure Reproduction Map

This repository reproduces manuscript Figures 2-10 from the scripts/notebooks below. Run them from the project root or from the `scripts_analysis/` directory; they resolve `inputs/` in either case.

| Figure | Script / Notebook | Main output |
| --- | --- | --- |
| 2 | `pca_expression.R` | `PCA_featureCounts_Low_High.pdf` |
| 3 | `volcano_dge.py` | `Volcano_deseqLvsH.pdf` |
| 4 | `GSEA.R` | `GO_BP_dotplot_top15_up_down.pdf` |
| 5 | `GSEA.R` | `GO_MF_dotplot_top15_up_down.pdf` |
| 6 | `mismatches_delta.R` | `mean_delta_by_6_muttypes_ordered.pdf` |
| 7 | `statistics_mismatches_mutations.R` | `betabinomial_estimated_MR_and_fold_change.pdf` |
| 8 | `statistics_mismatches_mutations.R` | `nbinom_estimated_VR_and_fold_change.pdf` |
| 9 | `recurrence_calculation.R` | `high_low_recurrence.pdf` |
| 10 | `shared_sites.R` | `shared_site_subtraction.pdf` |
