"""MAD-based per-sample QC filtering.

Scientific basis: Lun et al. (2016) Genome Biology; Heumos et al. (2023) Nature Methods.
Applying MAD per sample (not per dataset) avoids batch-driven false positives.
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def _mad(arr: np.ndarray) -> float:
    median = np.median(arr)
    return float(np.median(np.abs(arr - median)))


def is_outlier(
    series: pd.Series,
    n_mads: float = 5.0,
    log_transform: bool = True,
    tail: str = "both",
) -> pd.Series:
    """Flag cells more than ``n_mads`` MADs from the per-sample median.

    Parameters
    ----------
    series
        QC metric values for one sample.
    n_mads
        MAD multiplier. 5 is recommended for scRNA-seq (Lun et al. 2016).
    log_transform
        Log1p-transform before computing MAD (appropriate for count-based metrics).
    tail
        ``"both"`` | ``"upper"`` | ``"lower"`` — which tail(s) to flag.
    """
    vals = np.log1p(series.values.astype(float)) if log_transform else series.values.astype(float)
    median = np.median(vals)
    mad = _mad(vals)

    deviation = vals - median
    outlier = np.zeros(len(vals), dtype=bool)
    if tail in ("both", "lower"):
        outlier |= deviation < -n_mads * mad
    if tail in ("both", "upper"):
        outlier |= deviation > n_mads * mad
    return pd.Series(outlier, index=series.index)


def compute_qc_metrics(adata: ad.AnnData, mt_prefix: str = "MT-") -> ad.AnnData:
    """Annotate cells with standard QC metrics in-place.

    Adds to ``adata.obs``: ``n_genes_by_counts``, ``total_counts``, ``pct_counts_mt``.
    """
    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    return adata


def flag_outliers(
    adata: ad.AnnData,
    sample_key: str = "sample_id",
    n_mads: float = 5.0,
) -> ad.AnnData:
    """Flag QC outliers per sample using MAD-based detection.

    Adds ``adata.obs["outlier"]`` (bool) and per-metric flags.
    All computed bounds are stored in ``adata.uns["qc_bounds"]``.
    """
    obs = adata.obs.copy()
    obs["outlier_n_genes_low"] = False
    obs["outlier_n_genes_high"] = False
    obs["outlier_counts_low"] = False
    obs["outlier_counts_high"] = False
    obs["outlier_mt_high"] = False

    qc_bounds = {}

    for sample, idx in obs.groupby(sample_key).groups.items():
        sub = obs.loc[idx]

        # n_genes: flag both tails (too low = empty; too high = doublet signal)
        out_genes_low = is_outlier(sub["n_genes_by_counts"], n_mads=n_mads, tail="lower")
        out_genes_high = is_outlier(sub["n_genes_by_counts"], n_mads=n_mads, tail="upper")

        # total_counts: both tails
        out_counts_low = is_outlier(sub["total_counts"], n_mads=n_mads, tail="lower")
        out_counts_high = is_outlier(sub["total_counts"], n_mads=n_mads, tail="upper")

        # pct_mt: upper tail only (dying cells)
        out_mt = is_outlier(sub["pct_counts_mt"], n_mads=n_mads, log_transform=False, tail="upper")

        obs.loc[idx, "outlier_n_genes_low"] = out_genes_low.values
        obs.loc[idx, "outlier_n_genes_high"] = out_genes_high.values
        obs.loc[idx, "outlier_counts_low"] = out_counts_low.values
        obs.loc[idx, "outlier_counts_high"] = out_counts_high.values
        obs.loc[idx, "outlier_mt_high"] = out_mt.values

        qc_bounds[sample] = {
            "n_cells_before": len(idx),
            "n_genes_median": float(np.median(sub["n_genes_by_counts"])),
            "n_genes_mad": float(_mad(np.log1p(sub["n_genes_by_counts"].values.astype(float)))),
            "total_counts_median": float(np.median(sub["total_counts"])),
            "pct_mt_median": float(np.median(sub["pct_counts_mt"])),
        }

    obs["outlier"] = (
        obs["outlier_n_genes_low"]
        | obs["outlier_n_genes_high"]
        | obs["outlier_counts_low"]
        | obs["outlier_counts_high"]
        | obs["outlier_mt_high"]
    )

    adata.obs = obs
    adata.uns["qc_bounds"] = qc_bounds
    return adata


def apply_qc_filters(
    adata: ad.AnnData,
    min_cells_per_gene: int = 3,
    report_path: str | Path | None = None,
) -> ad.AnnData:
    """Remove outlier cells and lowly-expressed genes.

    Expects ``flag_outliers`` to have been run first (adds ``adata.obs["outlier"]``).
    Optionally writes per-sample QC report to CSV.
    """
    n_before = adata.n_obs

    # Cell filter
    adata = adata[~adata.obs["outlier"]].copy()

    # Gene filter: expressed in at least min_cells_per_gene
    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)

    n_after = adata.n_obs
    print(f"QC: retained {n_after:,} / {n_before:,} cells "
          f"({n_before - n_after:,} removed, {100 * n_after / n_before:.1f}%)")

    if report_path is not None:
        bounds = adata.uns.get("qc_bounds", {})
        rows = []
        for sample, stats in bounds.items():
            rows.append({"sample": sample, **stats})
        pd.DataFrame(rows).to_csv(report_path, index=False)
        print(f"QC report saved → {report_path}")

    return adata
