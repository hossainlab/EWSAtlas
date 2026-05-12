"""Doublet detection — must be run per sample before merging.

Primary: scrublet (pure Python, no R dependency).
Adds adata.obs["doublet_score"] and adata.obs["predicted_doublet"].

Reference: Wolock et al. (2019) Cell Systems; Germain et al. (2022) for scDblFinder benchmarks.
"""

import anndata as ad
import numpy as np
import scipy.sparse as sp


def detect_doublets(
    adata: ad.AnnData,
    sample_key: str = "sample_id",
    expected_doublet_rate: float = 0.06,
    random_state: int = 0,
) -> ad.AnnData:
    """Run scrublet doublet detection per sample.

    Parameters
    ----------
    adata
        AnnData with raw counts in ``adata.layers["counts"]``.
    sample_key
        Column in ``adata.obs`` identifying samples.
    expected_doublet_rate
        Expected multiplet rate. 10x default ~6% for ~8k cells/sample.
    random_state
        Reproducibility seed.

    Returns
    -------
    AnnData with ``adata.obs["doublet_score"]`` and ``adata.obs["predicted_doublet"]``.
    """
    try:
        import scrublet as scr
    except ImportError as e:
        raise ImportError(
            "scrublet not installed. Run: uv add scrublet"
        ) from e

    scores = np.zeros(adata.n_obs, dtype=float)
    calls = np.zeros(adata.n_obs, dtype=bool)

    for sample, idx in adata.obs.groupby(sample_key).groups.items():
        counts = adata[idx].layers["counts"]
        if sp.issparse(counts):
            counts = counts.toarray()

        scrub = scr.Scrublet(
            counts,
            expected_doublet_rate=expected_doublet_rate,
            random_state=random_state,
        )
        ds, pd_ = scrub.scrub_doublets(verbose=False)
        scores[adata.obs_names.get_indexer(idx)] = ds
        calls[adata.obs_names.get_indexer(idx)] = pd_
        n_doublets = pd_.sum()
        print(f"  {sample}: {n_doublets} predicted doublets / {len(idx)} cells "
              f"({100 * n_doublets / len(idx):.1f}%)")

    adata.obs["doublet_score"] = scores
    adata.obs["predicted_doublet"] = calls
    return adata
