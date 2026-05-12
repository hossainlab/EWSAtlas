"""Publication-ready plotting utilities."""

from __future__ import annotations

import anndata as ad
import matplotlib.pyplot as plt
import scanpy as sc


def plot_qc_violin(
    adata: ad.AnnData,
    groupby: str = "sample_id",
    save: str | None = None,
) -> None:
    """Violin plots of QC metrics per sample."""
    metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    fig, axes = plt.subplots(1, len(metrics), figsize=(6 * len(metrics), 5))

    for ax, metric in zip(axes, metrics):
        sc.pl.violin(adata, metric, groupby=groupby, ax=ax, show=False,
                     rotation=45, stripplot=False)
        ax.set_title(metric)

    plt.tight_layout()
    if save:
        fig.savefig(save, dpi=300, bbox_inches="tight")
    plt.show()


def plot_umap_grid(
    adata: ad.AnnData,
    color: list[str],
    ncols: int = 3,
    save: str | None = None,
) -> None:
    """UMAP grid colored by multiple metadata fields."""
    sc.pl.umap(
        adata,
        color=color,
        ncols=ncols,
        show=save is None,
        save=save,
    )
