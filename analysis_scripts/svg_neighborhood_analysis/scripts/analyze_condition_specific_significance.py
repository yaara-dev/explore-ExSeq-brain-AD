#!/usr/bin/env python3
"""
Condition-Specific Consistent CELINA Significance
=================================================

Finds genes that are reproducibly CELINA-significant within one genotype
(WT or 5xFAD) across biological replicates for a given cell type, but NOT
reproducibly significant in the other genotype.

A gene is "consistently significant" in a condition/cell-type when it is
significant (adjusted p < --alpha) in **all** applicable samples of that
genotype, and there are at least --min-replicates applicable samples
(e.g. 3 of 3 for Xenium with three replicates per condition).

Classification per (gene, cell_type):
  - WT_exclusive     : consistent in WT, not consistent in 5xFAD
  - 5xFAD_exclusive  : consistent in 5xFAD, not consistent in WT
  - both_consistent  : consistent in both
  - neither_consistent: consistent in neither

All outputs carry explicit n_significant/n_applicable replicate counts so the
consistency call is fully transparent (e.g. WT "3 of 3", 5x "1 of 4").

This is distinct from the Table S10 co-clustering consistency analysis:
here we ask whether individual genes pass the CELINA significance threshold
consistently, not whether they co-cluster.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Canonical sample -> genotype mapping
# (matches analyze_clustering_patterns.identify_sample_conditions)
SAMPLE_CONDITIONS = {
    # Stellar / ExSeq
    "fem2_5x_F5_B": "5xFAD",
    "fem3_5x_E7_A_left": "5xFAD",
    "fem4_5x_F8_A_R": "5xFAD",
    "fem5_5x_left": "5xFAD",
    "fem6_5x_13C6_right": "5xFAD",
    "fem7_5x_14C7": "5xFAD",
    "fem2_WT_F3_B_left": "WT",
    "fem3_WTE1_B": "WT",
    "fem4_WT_F11": "WT",
    "fem5_WT_left": "WT",
    "fem6_WT_11C9_left": "WT",
    "fem7_WT_12C8": "WT",
    # Xenium
    "WT_6_section1": "WT",
    "WT_6_section2": "WT",
    "WT_5": "WT",
    "5xFAD_6": "5xFAD",
    "5xFAD_5_section1": "5xFAD",
    "5xFAD_5_section2": "5xFAD",
}


def load_sample_pvalues(input_dir: Path, pvalues_file: str) -> dict[str, pd.DataFrame]:
    """Load the adjusted p-value matrix for every known stellar sample.

    Returns {sample_name: DataFrame} where the DataFrame is indexed by gene and
    columns are cell types (adjusted p-values, NaN where not tested).
    """
    sample_data: dict[str, pd.DataFrame] = {}
    for sample in SAMPLE_CONDITIONS:
        path = input_dir / sample / pvalues_file
        if not path.exists():
            print(f"  WARNING: missing p-value file for {sample}: {path}")
            continue
        df = pd.read_csv(path)
        if "gene" not in df.columns:
            print(f"  WARNING: {path} has no 'gene' column, skipping")
            continue
        df = df.set_index("gene")
        # Coerce all cell-type columns to numeric (blank -> NaN)
        df = df.apply(pd.to_numeric, errors="coerce")
        sample_data[sample] = df
        print(f"  Loaded {sample} ({SAMPLE_CONDITIONS[sample]}): "
              f"{df.shape[0]} genes x {df.shape[1]} cell types")
    return sample_data


def samples_by_condition(sample_data: dict[str, pd.DataFrame], condition: str) -> list[str]:
    return [s for s in sample_data if SAMPLE_CONDITIONS[s] == condition]


def collect_genes_and_celltypes(sample_data: dict[str, pd.DataFrame]) -> tuple[list[str], list[str]]:
    genes: set[str] = set()
    cell_types: set[str] = set()
    for df in sample_data.values():
        genes.update(df.index.astype(str))
        cell_types.update(df.columns.astype(str))
    return sorted(genes), sorted(cell_types)


def evaluate_condition(
    gene: str,
    cell_type: str,
    condition_samples: list[str],
    sample_data: dict[str, pd.DataFrame],
    alpha: float,
) -> dict:
    """Return replicate stats for one gene/cell_type within one condition.

    'applicable' = sample where the cell type column exists and the gene has a
    non-null adjusted p-value (i.e. the gene was actually tested there).
    """
    applicable_samples: list[str] = []
    significant_samples: list[str] = []
    q_values: list[float] = []

    for sample in condition_samples:
        df = sample_data[sample]
        if cell_type not in df.columns or gene not in df.index:
            continue
        q = df.at[gene, cell_type]
        if pd.isna(q):
            continue
        applicable_samples.append(sample)
        q_values.append(float(q))
        if q < alpha:
            significant_samples.append(sample)

    return {
        "applicable_samples": applicable_samples,
        "significant_samples": significant_samples,
        "n_applicable": len(applicable_samples),
        "n_significant": len(significant_samples),
        "min_q": min(q_values) if q_values else np.nan,
        "max_q": max(q_values) if q_values else np.nan,
    }


def format_consistency(n_significant: int, n_applicable: int) -> str:
    # Use "of" instead of "/" so spreadsheet tools don't parse values as dates.
    return f"{n_significant} of {n_applicable}"


def is_consistent(
    n_significant: int,
    n_applicable: int,
    min_replicates: int,
    require_all_applicable: bool = True,
) -> bool:
    if n_applicable < min_replicates:
        return False
    if require_all_applicable:
        return n_significant == n_applicable
    return n_significant >= min_replicates


def classify(wt_consistent: bool, x5_consistent: bool) -> str:
    if wt_consistent and x5_consistent:
        return "both_consistent"
    if wt_consistent:
        return "WT_exclusive"
    if x5_consistent:
        return "5xFAD_exclusive"
    return "neither_consistent"


def build_gene_table(
    sample_data: dict[str, pd.DataFrame],
    alpha: float,
    min_replicates: int,
    require_all_applicable: bool = True,
) -> pd.DataFrame:
    wt_samples = samples_by_condition(sample_data, "WT")
    x5_samples = samples_by_condition(sample_data, "5xFAD")
    genes, cell_types = collect_genes_and_celltypes(sample_data)

    rows = []
    for cell_type in cell_types:
        for gene in genes:
            wt = evaluate_condition(gene, cell_type, wt_samples, sample_data, alpha)
            x5 = evaluate_condition(gene, cell_type, x5_samples, sample_data, alpha)

            # Only keep (gene, cell_type) pairs that were tested somewhere
            if wt["n_applicable"] == 0 and x5["n_applicable"] == 0:
                continue

            wt_consistent = is_consistent(
                wt["n_significant"], wt["n_applicable"], min_replicates, require_all_applicable
            )
            x5_consistent = is_consistent(
                x5["n_significant"], x5["n_applicable"], min_replicates, require_all_applicable
            )

            rows.append({
                "cell_type": cell_type,
                "gene": gene,
                "classification": classify(wt_consistent, x5_consistent),
                "wt_consistent": wt_consistent,
                "x5_consistent": x5_consistent,
                "wt_consistency": format_consistency(wt["n_significant"], wt["n_applicable"]),
                "x5_consistency": format_consistency(x5["n_significant"], x5["n_applicable"]),
                "n_wt_significant": wt["n_significant"],
                "n_wt_applicable": wt["n_applicable"],
                "n_5x_significant": x5["n_significant"],
                "n_5x_applicable": x5["n_applicable"],
                "wt_significant_samples": ", ".join(wt["significant_samples"]),
                "x5_significant_samples": ", ".join(x5["significant_samples"]),
                "wt_applicable_samples": ", ".join(wt["applicable_samples"]),
                "x5_applicable_samples": ", ".join(x5["applicable_samples"]),
                "min_q_wt": wt["min_q"],
                "max_q_wt": wt["max_q"],
                "min_q_5x": x5["min_q"],
                "max_q_5x": x5["max_q"],
            })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["cell_type", "classification", "gene"]).reset_index(drop=True)
    return df


def build_summary_table(
    gene_table: pd.DataFrame,
    sample_data: dict[str, pd.DataFrame],
    min_replicates: int,
) -> pd.DataFrame:
    wt_samples = samples_by_condition(sample_data, "WT")
    x5_samples = samples_by_condition(sample_data, "5xFAD")

    def n_samples_with_celltype(samples, cell_type):
        return sum(
            1 for s in samples
            if cell_type in sample_data[s].columns
            and sample_data[s][cell_type].notna().any()
        )

    rows = []
    for cell_type, group in gene_table.groupby("cell_type"):
        n_wt = n_samples_with_celltype(wt_samples, cell_type)
        n_5x = n_samples_with_celltype(x5_samples, cell_type)
        counts = group["classification"].value_counts()
        rows.append({
            "cell_type": cell_type,
            "n_wt_applicable_samples": n_wt,
            "n_5x_applicable_samples": n_5x,
            "wt_consistency_eligible": n_wt >= min_replicates,
            "x5_consistency_eligible": n_5x >= min_replicates,
            "WT_exclusive": int(counts.get("WT_exclusive", 0)),
            "5xFAD_exclusive": int(counts.get("5xFAD_exclusive", 0)),
            "both_consistent": int(counts.get("both_consistent", 0)),
            "neither_consistent": int(counts.get("neither_consistent", 0)),
            "n_tested_pairs": len(group),
        })
    return pd.DataFrame(rows).sort_values("cell_type").reset_index(drop=True)


def build_exclusive_pivot(gene_table: pd.DataFrame, classification: str, consistency_col: str) -> pd.DataFrame:
    subset = gene_table[gene_table["classification"] == classification]
    rows = []
    for cell_type, group in subset.groupby("cell_type"):
        group = group.sort_values("gene")
        genes = group["gene"].tolist()
        details = [f"{g} ({c})" for g, c in zip(group["gene"], group[consistency_col])]
        rows.append({
            "cell_type": cell_type,
            "n_genes": len(genes),
            "genes": ", ".join(genes),
            "genes_with_consistency": "; ".join(details),
        })
    if not rows:
        return pd.DataFrame(columns=[
            "cell_type", "n_genes", "genes", "genes_with_consistency",
        ])
    return pd.DataFrame(rows).sort_values("cell_type").reset_index(drop=True)


def write_markdown_report(
    path: Path,
    gene_table: pd.DataFrame,
    summary: pd.DataFrame,
    alpha: float,
    min_replicates: int,
    pvalues_file: str,
    require_all_applicable: bool = True,
) -> None:
    lines = []
    lines.append("# Condition-Specific Consistent CELINA Significance")
    lines.append("")
    consistency_rule = (
        f"significant in **all** applicable samples (≥{min_replicates} required), "
        f"e.g. 3 of 3 for three replicates"
        if require_all_applicable
        else f"significant in at least **{min_replicates} replicates (samples)** "
             f"(not necessarily every applicable sample)"
    )
    lines.append(
        f"Genes that are reproducibly CELINA-significant (adjusted p < {alpha}) "
        f"across biological replicates within one genotype for a cell type, but "
        f"not the other. A gene is 'consistent' in a condition when it is "
        f"{consistency_rule}."
    )
    lines.append("")
    lines.append(f"**P-value source:** `{pvalues_file}`  ")
    lines.append(f"**Significance threshold:** adjusted p < {alpha}  ")
    lines.append(f"**Minimum applicable samples:** {min_replicates}")
    lines.append(
        f"**Require all applicable significant:** {require_all_applicable}"
    )
    lines.append("")
    lines.append("Consistency strings are shown as `n_significant of n_applicable`.")
    lines.append("")

    lines.append("## Summary by cell type")
    lines.append("")
    lines.append("| Cell type | WT samples | 5xFAD samples | WT eligible | 5xFAD eligible | WT-exclusive | 5xFAD-exclusive | Both | Neither |")
    lines.append("|-----------|-----------|---------------|-------------|----------------|--------------|-----------------|------|---------|")
    for _, r in summary.iterrows():
        lines.append(
            f"| {r['cell_type']} | {r['n_wt_applicable_samples']} | {r['n_5x_applicable_samples']} "
            f"| {r['wt_consistency_eligible']} | {r['x5_consistency_eligible']} "
            f"| {r['WT_exclusive']} | {r['5xFAD_exclusive']} | {r['both_consistent']} | {r['neither_consistent']} |"
        )
    lines.append("")

    for cell_type, group in gene_table.groupby("cell_type"):
        srow = summary[summary["cell_type"] == cell_type].iloc[0]
        lines.append(f"## {cell_type}")
        lines.append("")
        lines.append(
            f"Applicable samples: {srow['n_wt_applicable_samples']} WT, "
            f"{srow['n_5x_applicable_samples']} 5xFAD "
            f"(min {min_replicates} required for a consistency call). "
            f"WT eligible: {srow['wt_consistency_eligible']}, "
            f"5xFAD eligible: {srow['x5_consistency_eligible']}."
        )
        lines.append("")

        wt_only = group[group["classification"] == "WT_exclusive"].sort_values("gene")
        x5_only = group[group["classification"] == "5xFAD_exclusive"].sort_values("gene")
        both = group[group["classification"] == "both_consistent"].sort_values("gene")

        lines.append(f"**WT-exclusive ({len(wt_only)}):**")
        if len(wt_only):
            lines.append("")
            for _, r in wt_only.iterrows():
                lines.append(f"- {r['gene']} — WT {r['wt_consistency']}, 5x {r['x5_consistency']}")
        else:
            lines.append(" none")
        lines.append("")

        lines.append(f"**5xFAD-exclusive ({len(x5_only)}):**")
        if len(x5_only):
            lines.append("")
            for _, r in x5_only.iterrows():
                lines.append(f"- {r['gene']} — WT {r['wt_consistency']}, 5x {r['x5_consistency']}")
        else:
            lines.append(" none")
        lines.append("")

        lines.append(f"**Consistent in both ({len(both)}):**")
        if len(both):
            lines.append("")
            for _, r in both.iterrows():
                lines.append(f"- {r['gene']} — WT {r['wt_consistency']}, 5x {r['x5_consistency']}")
        else:
            lines.append(" none")
        lines.append("")

    path.write_text("\n".join(lines))


def plot_counts_by_celltype(summary: pd.DataFrame, path: Path) -> None:
    if summary.empty:
        return
    categories = ["WT_exclusive", "5xFAD_exclusive", "both_consistent", "neither_consistent"]
    cell_types = summary["cell_type"].tolist()
    x = np.arange(len(cell_types))

    fig, ax = plt.subplots(figsize=(max(8, len(cell_types) * 1.4), 6))
    bottom = np.zeros(len(cell_types))
    for cat in categories:
        vals = summary[cat].to_numpy()
        ax.bar(x, vals, bottom=bottom, label=cat)
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels(cell_types, rotation=45, ha="right")
    ax.set_ylabel("Number of (gene, cell type) pairs")
    ax.set_title("Condition-specific CELINA significance by cell type")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", default="data/stellar",
                        help="Directory containing per-sample stellar folders")
    parser.add_argument("--pvalues-file", default="all_p_values_FDR_global_bonferroni.csv",
                        help="Adjusted p-value filename inside each sample folder")
    parser.add_argument("--alpha", type=float, default=0.01,
                        help="Significance threshold on adjusted p-values")
    parser.add_argument("--min-replicates", type=int, default=3,
                        help="Minimum applicable samples per condition for a "
                             "consistency call (default: 3)")
    parser.add_argument(
        "--allow-partial-replicates",
        action="store_true",
        help="Count as consistent when significant in at least --min-replicates "
             "samples (old behavior), instead of all applicable samples",
    )
    parser.add_argument("--output-dir",
                        default="condition_specific_significance_analysis/results",
                        help="Output directory for tables, report, and plots")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    plots_dir = output_dir / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    print("Loading per-sample adjusted p-values...")
    sample_data = load_sample_pvalues(input_dir, args.pvalues_file)
    if not sample_data:
        raise SystemExit("No sample p-value files loaded; aborting.")

    n_wt = len(samples_by_condition(sample_data, "WT"))
    n_5x = len(samples_by_condition(sample_data, "5xFAD"))
    print(f"\nLoaded {len(sample_data)} samples ({n_wt} WT, {n_5x} 5xFAD)")
    require_all = not args.allow_partial_replicates
    print(f"Alpha = {args.alpha}, min applicable = {args.min_replicates}, "
          f"require all applicable significant = {require_all}\n")

    gene_table = build_gene_table(
        sample_data, args.alpha, args.min_replicates, require_all
    )
    if gene_table.empty:
        raise SystemExit("No tested (gene, cell_type) pairs found; aborting.")

    summary = build_summary_table(gene_table, sample_data, args.min_replicates)
    wt_pivot = build_exclusive_pivot(gene_table, "WT_exclusive", "wt_consistency")
    x5_pivot = build_exclusive_pivot(gene_table, "5xFAD_exclusive", "x5_consistency")

    genes_csv = output_dir / "condition_specific_significant_genes.csv"
    summary_csv = output_dir / "condition_specific_significance_summary.csv"
    wt_csv = output_dir / "WT_exclusive_genes_by_celltype.csv"
    x5_csv = output_dir / "5xFAD_exclusive_genes_by_celltype.csv"
    report_md = output_dir / "CONDITION_SPECIFIC_CELINA_SIGNIFICANCE.md"
    counts_png = plots_dir / "counts_by_celltype.png"

    gene_table.to_csv(genes_csv, index=False)
    summary.to_csv(summary_csv, index=False)
    wt_pivot.to_csv(wt_csv, index=False)
    x5_pivot.to_csv(x5_csv, index=False)
    write_markdown_report(report_md, gene_table, summary, args.alpha,
                          args.min_replicates, args.pvalues_file, require_all)
    plot_counts_by_celltype(summary, counts_png)

    print("Outputs written:")
    for p in [genes_csv, summary_csv, wt_csv, x5_csv, report_md, counts_png]:
        print(f"  {p}")

    total_wt = int(summary["WT_exclusive"].sum())
    total_5x = int(summary["5xFAD_exclusive"].sum())
    total_both = int(summary["both_consistent"].sum())
    print(f"\nTotals across cell types: "
          f"WT-exclusive={total_wt}, 5xFAD-exclusive={total_5x}, both={total_both}")


if __name__ == "__main__":
    main()
