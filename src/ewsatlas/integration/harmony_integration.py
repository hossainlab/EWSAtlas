"""Harmony integration — fast batch correction for sanity checking.

Corrects PCA embedding per patient_id. Use for rapid UMAP inspection
before committing to the longer scVI run.
"""

import anndata as ad
import scanpy as sc


def run_harmony(
    adata: ad.AnnData,
    batch_key: str = "patient_id",
    n_hvgs: int = 3000,
    n_pcs: int = 50,
    random_state: int = 0,
) -> ad.AnnData:
    """Normalize, select HVGs, run PCA, then Harmony batch correction.

    Parameters
    ----------
    adata
        AnnData with raw counts in ``adata.layers["counts"]``.
        ``adata.X`` will be overwritten with normalized values.
    batch_key
        obs column for batch correction.
    n_hvgs
        Number of highly variable genes.
    n_pcs
        Number of PCs before Harmony.
    random_state
        Reproducibility seed.

    Returns
    -------
    AnnData with ``adata.obsm["X_pca_harmony"]`` and ``adata.obsm["X_umap"]``.
    """
    try:
        import harmonypy as hm
    except ImportError as e:
        raise ImportError("harmonypy not installed. Run: uv add harmonypy") from e

    adata = adata.copy()
    adata.X = adata.layers["counts"].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_hvgs,
        batch_key=batch_key,
        flavor="seurat_v3",
    )
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()

    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=n_pcs, random_state=random_state)

    ho = hm.run_harmony(
        adata_hvg.obsm["X_pca"],
        adata_hvg.obs,
        batch_key,
        random_state=random_state,
    )
    adata.obsm["X_pca_harmony"] = ho.Z_corr.T

    sc.pp.neighbors(adata, use_rep="X_pca_harmony", random_state=random_state)
    sc.tl.umap(adata, random_state=random_state)
    sc.tl.leiden(adata, resolution=0.5, random_state=random_state)

    return adata
