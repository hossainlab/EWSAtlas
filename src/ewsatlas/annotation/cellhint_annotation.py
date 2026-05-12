"""Cell type annotation using CellHint + EwS-specific marker scoring.

CellHint: automated hierarchical cell type annotation.
EwS markers: NKX2-2, CD99, VCAN, PRKCB, EWSR1 for tumor cell identification.
EWSR1-FLI1 activity: Franzetti 2020 high/low signatures for heterogeneity axis.
"""

from __future__ import annotations

import anndata as ad
import scanpy as sc


_EWS_MARKERS = ["NKX2-2", "CD99", "VCAN", "PRKCB", "EWSR1"]

# Franzetti et al. 2020 EWSR1-FLI1 activity signatures
# High activity = EWSR1-FLI1-high (proliferative state)
# Low activity = EWSR1-FLI1-low (migratory/mesenchymal state)
_FLI1_HIGH_SIGNATURE = [
    "NKX2-2", "CCND1", "EZH2", "CAV1", "GLG1", "STEAP1", "STEAP2",
    "NPY1R", "FCGRT", "SLC3A2",
]
_FLI1_LOW_SIGNATURE = [
    "SERPINE1", "TGFB1", "VIM", "FN1", "SPARC", "COL1A1", "COL3A1",
    "IGFBP3", "ACTA2", "CDH2",
]


def annotate_celltypes(
    adata: ad.AnnData,
    reference: str = "hlca",
    min_confidence: float = 0.5,
) -> ad.AnnData:
    """Run CellHint cell type annotation.

    Parameters
    ----------
    adata
        AnnData with normalized log1p counts in ``adata.X``.
    reference
        CellHint reference atlas (e.g. ``"hlca"`` for Human Lung Cell Atlas).
    min_confidence
        Minimum confidence threshold for accepting CellHint labels.

    Returns
    -------
    AnnData with ``adata.obs["cellhint_label"]`` and ``adata.obs["cellhint_confidence"]``.
    """
    try:
        import cellhint
    except ImportError as e:
        raise ImportError("cellhint not installed. Run: uv add cellhint") from e

    result = cellhint.tl.predict(adata, reference=reference)
    adata.obs["cellhint_label"] = result["predicted_label"]
    adata.obs["cellhint_confidence"] = result["confidence"]
    adata.obs["cellhint_label_filtered"] = adata.obs.apply(
        lambda r: r["cellhint_label"] if r["cellhint_confidence"] >= min_confidence else "uncertain",
        axis=1,
    )
    return adata


def score_ews_markers(adata: ad.AnnData) -> ad.AnnData:
    """Score EwS tumor markers and EWSR1-FLI1 activity on each cell.

    Adds to ``adata.obs``:
    - ``score_ews_markers``: composite EwS marker score
    - ``score_fli1_high``: EWSR1-FLI1 high-activity signature score
    - ``score_fli1_low``: EWSR1-FLI1 low-activity signature score
    - ``fli1_activity``: continuous activity index (high - low)
    """
    present_markers = [g for g in _EWS_MARKERS if g in adata.var_names]
    if present_markers:
        sc.tl.score_genes(adata, present_markers, score_name="score_ews_markers")
    else:
        print("Warning: no EwS markers found in var_names")

    present_high = [g for g in _FLI1_HIGH_SIGNATURE if g in adata.var_names]
    present_low = [g for g in _FLI1_LOW_SIGNATURE if g in adata.var_names]

    if present_high:
        sc.tl.score_genes(adata, present_high, score_name="score_fli1_high")
    if present_low:
        sc.tl.score_genes(adata, present_low, score_name="score_fli1_low")

    if present_high and present_low:
        adata.obs["fli1_activity"] = (
            adata.obs["score_fli1_high"] - adata.obs["score_fli1_low"]
        )

    return adata
