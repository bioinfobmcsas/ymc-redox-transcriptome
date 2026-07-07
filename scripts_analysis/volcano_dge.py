#!/usr/bin/env python
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


project_root = Path(".").resolve() if Path("inputs").is_dir() else Path("..").resolve()
input_file = project_root / "inputs" / "deseq2_LvsH_wrs3_s1_genes.xlsx"
output_file = project_root / "Volcano_deseqLvsH.pdf"

padj_cut = 0.05
lfc_col = "log2FoldChange"
padj_col = "padj"
ymax = 40.0

plt.rcParams["font.family"] = "Times New Roman"

df = pd.read_excel(input_file)

required_cols = {lfc_col, padj_col}
missing_cols = required_cols.difference(df.columns)

if missing_cols:
    missing = ", ".join(sorted(missing_cols))
    raise ValueError(f"Missing required columns in {input_file}: {missing}")

df = df.replace([np.inf, -np.inf], np.nan)
df = df.dropna(subset=[lfc_col, padj_col])
df = df[df[padj_col] >= 0].copy()

min_positive = np.nextafter(0, 1)
df["neg_log10_padj"] = -np.log10(df[padj_col].clip(lower=min_positive))
df["neg_log10_padj"] = df["neg_log10_padj"].clip(upper=ymax)

sig = df[padj_col] <= padj_cut
pos = df[lfc_col] > 1
neg = df[lfc_col] < -1
sig_no_fc = sig & ~pos & ~neg

fig, ax = plt.subplots(figsize=(8, 6))

not_sig = ax.scatter(
    df.loc[~sig, lfc_col],
    df.loc[~sig, "neg_log10_padj"],
    c="lightgrey",
    alpha=1.0,
    edgecolors="none",
)

up = ax.scatter(
    df.loc[sig & pos, lfc_col],
    df.loc[sig & pos, "neg_log10_padj"],
    c="red",
    alpha=0.8,
    edgecolors="none",
)

sig_mid = ax.scatter(
    df.loc[sig_no_fc, lfc_col],
    df.loc[sig_no_fc, "neg_log10_padj"],
    c="#C8A23A",
    alpha=0.9,
    edgecolors="none",
)

down = ax.scatter(
    df.loc[sig & neg, lfc_col],
    df.loc[sig & neg, "neg_log10_padj"],
    c="blue",
    alpha=0.8,
    edgecolors="none",
)

ax.axhline(-np.log10(padj_cut), color="grey", linestyle="--", linewidth=1)
ax.axvline(0, color="black", linestyle=":", linewidth=1)
ax.set_ylim(0, ymax)
ax.set_xlabel(r"$\mathregular{\log_{2}}$ fold change (LowDO / HighDO)", fontsize=14)
ax.set_ylabel(r"$\mathregular{-\log_{10}}$ adjusted p-value", fontsize=14)
ax.tick_params(axis="both", which="major", labelsize=14)
ax.legend(
    (not_sig, up, sig_mid, down),
    (
        f"padj > {padj_cut:.2g}",
        f"padj <= {padj_cut:.2g}, log2FC > 1",
        f"padj <= {padj_cut:.2g}, |log2FC| <= 1",
        f"padj <= {padj_cut:.2g}, log2FC < -1",
    ),
    loc="upper right",
    shadow=True,
    fontsize=14,
)

fig.tight_layout()
fig.savefig(output_file)
plt.close(fig)

print(f"Volcano plot saved to {output_file}")
