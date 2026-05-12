"""scVI integration — probabilistic batch correction for final atlas.

Uses negative binomial likelihood to model count data.
Batch key: patient_id (not platform) to preserve biological signal
while correcting patient-specific technical variation.
"""

import anndata as ad
import scanpy as sc


def run_scvi(
    adata: ad.AnnData,
    batch_key: str = "patient_id",
    n_hvgs: int = 3000,
    n_latent: int = 30,
    n_layers: int = 2,
    n_hidden: int = 128,
    max_epochs: int = 400,
    early_stopping: bool = True,
    random_state: int = 0,
) -> ad.AnnData:
    """Select HVGs, train scVI model, embed into latent space, compute UMAP.

    Parameters
    ----------
    adata
        AnnData with raw counts in ``adata.layers["counts"]``.
    batch_key
        obs column used as batch covariate in scVI.
    n_hvgs
        Number of highly variable genes.
    n_latent
        Dimensionality of scVI latent space.
    n_layers
        Number of layers in encoder/decoder.
    n_hidden
        Hidden layer width.
    max_epochs
        Maximum training epochs.
    early_stopping
        Stop early if validation loss stagnates.
    random_state
        Reproducibility seed.

    Returns
    -------
    AnnData with:
    - ``adata.obsm["X_scVI"]``: latent representation
    - ``adata.obsm["X_umap"]``: UMAP on scVI neighbors
    - ``adata.obs["leiden"]``: Leiden clusters
    - ``adata.uns["scvi_model_params"]``: training config
    """
    try:
        import scvi
    except ImportError as e:
        raise ImportError("scvi-tools not installed. Run: uv add scvi-tools") from e

    scvi.settings.seed = random_state

    adata = adata.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_hvgs,
        batch_key=batch_key,
        flavor="seurat_v3",
    )

    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    adata_hvg.layers["counts"] = adata.layers["counts"][
        :, adata.var["highly_variable"]
    ]

    scvi.model.SCVI.setup_anndata(
        adata_hvg,
        layer="counts",
        batch_key=batch_key,
    )

    model = scvi.model.SCVI(
        adata_hvg,
        n_latent=n_latent,
        n_layers=n_layers,
        n_hidden=n_hidden,
        gene_likelihood="nb",
    )
    model.train(
        max_epochs=max_epochs,
        early_stopping=early_stopping,
        plan_kwargs={"lr": 1e-3},
    )

    adata_hvg.obsm["X_scVI"] = model.get_latent_representation()

    sc.pp.neighbors(adata_hvg, use_rep="X_scVI", random_state=random_state)
    sc.tl.umap(adata_hvg, random_state=random_state)
    sc.tl.leiden(adata_hvg, resolution=0.5, random_state=random_state)

    # Transfer embeddings back to full-gene AnnData
    adata.obsm["X_scVI"] = adata_hvg.obsm["X_scVI"]
    adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
    adata.obs["leiden"] = adata_hvg.obs["leiden"]
    adata.obsp = adata_hvg.obsp
    adata.uns["neighbors"] = adata_hvg.uns["neighbors"]
    adata.uns["scvi_model_params"] = {
        "batch_key": batch_key,
        "n_latent": n_latent,
        "n_layers": n_layers,
        "n_hidden": n_hidden,
        "gene_likelihood": "nb",
        "max_epochs": max_epochs,
        "n_hvgs": n_hvgs,
        "random_state": random_state,
    }

    return adata, model
