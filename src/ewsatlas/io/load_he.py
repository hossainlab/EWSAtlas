"""Loader for He et al. 2025 — 10x Chromium MTX format."""

from pathlib import Path

import anndata as ad
import scanpy as sc


_SAMPLE_METADATA = {
    "TN.1": {"treatment": "naive", "patient_id": "He_TN1"},
    "TN.2": {"treatment": "naive", "patient_id": "He_TN2"},
    "TN.3": {"treatment": "naive", "patient_id": "He_TN3"},
    "NAC.1": {"treatment": "neoadjuvant", "patient_id": "He_NAC1"},
    "NAC.3": {"treatment": "neoadjuvant", "patient_id": "He_NAC3"},
    "RC.1": {"treatment": "relapsed", "patient_id": "He_RC1"},
    "RC.2": {"treatment": "relapsed", "patient_id": "He_RC2"},
    "RC.3": {"treatment": "relapsed", "patient_id": "He_RC3"},
}


def load_he2025(data_dir: str | Path) -> ad.AnnData:
    """Load all He et al. 2025 samples from MTX format into a single AnnData.

    Parameters
    ----------
    data_dir
        Path to ``data/raw/He2025/`` containing per-sample subdirectories.

    Returns
    -------
    Concatenated AnnData with standardized obs fields.
    """
    data_dir = Path(data_dir)
    adatas = []

    for sample_name, meta in _SAMPLE_METADATA.items():
        sample_dir = data_dir / sample_name
        if not sample_dir.exists():
            continue

        # genes.tsv is uncompressed 2-col format (ENSEMBL_ID \t SYMBOL).
        # var_names="gene_symbols" uses column index 1 (the symbol).
        adata = sc.read_10x_mtx(
            sample_dir,
            var_names="gene_symbols",
            make_unique=True,
            cache=False,
            gex_only=False,
        )
        adata.obs_names = [f"{sample_name}_{bc}" for bc in adata.obs_names]
        adata.obs["sample_id"] = sample_name
        adata.obs["patient_id"] = meta["patient_id"]
        adata.obs["treatment"] = meta["treatment"]
        adata.obs["dataset"] = "He2025"
        adata.obs["platform"] = "10x_Chromium"
        adata.obs["tissue"] = "primary_tumor"

        adatas.append(adata)

    combined = ad.concat(adatas, join="outer", fill_value=0)
    combined.var_names_make_unique()
    combined.layers["counts"] = combined.X.copy()
    return combined
