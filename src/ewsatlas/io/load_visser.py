"""Loader for Visser et al. 2023 — CEL-Seq2 transcript count format.

Each file: genes × wells (384-well plate), tab-separated, gzip-compressed.
Gene IDs: ``ENSG00000000003__TSPAN6`` (ENSEMBL__SYMBOL).
"""

import gzip
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


# Treatment assignments from Visser et al. 2023 (Nature Communications)
# TM codes correspond to patient tumor biopsies.
# Samples not listed here default to treatment="unknown" via fallback.
_PATIENT_MAP = {
    "TM338": {"patient_id": "Visser_TM338", "treatment": "naive"},
    "TM339": {"patient_id": "Visser_TM339", "treatment": "naive"},
    "TM344": {"patient_id": "Visser_TM344", "treatment": "naive"},
    "TM348": {"patient_id": "Visser_TM348", "treatment": "naive"},
    "TM416": {"patient_id": "Visser_TM416", "treatment": "naive"},
    "TM417": {"patient_id": "Visser_TM417", "treatment": "neoadjuvant"},
    "TM424": {"patient_id": "Visser_TM424", "treatment": "neoadjuvant"},
    "TM425": {"patient_id": "Visser_TM425", "treatment": "neoadjuvant"},
    "TM505": {"patient_id": "Visser_TM505", "treatment": "relapsed"},
    "TM506": {"patient_id": "Visser_TM506", "treatment": "relapsed"},
    "TM507": {"patient_id": "Visser_TM507", "treatment": "relapsed"},
    "TM508": {"patient_id": "Visser_TM508", "treatment": "relapsed"},
    "TM547": {"patient_id": "Visser_TM547", "treatment": "unknown"},
    "TM548": {"patient_id": "Visser_TM548", "treatment": "unknown"},
    "TM549": {"patient_id": "Visser_TM549", "treatment": "unknown"},
    "TM552": {"patient_id": "Visser_TM552", "treatment": "unknown"},
    "TM564": {"patient_id": "Visser_TM564", "treatment": "unknown"},
    "TM570": {"patient_id": "Visser_TM570", "treatment": "unknown"},
    "TM572": {"patient_id": "Visser_TM572", "treatment": "unknown"},
    "TM574": {"patient_id": "Visser_TM574", "treatment": "unknown"},
    "TM707": {"patient_id": "Visser_TM707", "treatment": "unknown"},
    "TM709": {"patient_id": "Visser_TM709", "treatment": "unknown"},
    "TM712": {"patient_id": "Visser_TM712", "treatment": "unknown"},
    "TM734": {"patient_id": "Visser_TM734", "treatment": "unknown"},
    "TM736": {"patient_id": "Visser_TM736", "treatment": "unknown"},
    "TM737": {"patient_id": "Visser_TM737", "treatment": "unknown"},
    "TM739": {"patient_id": "Visser_TM739", "treatment": "unknown"},
    "TM768": {"patient_id": "Visser_TM768", "treatment": "unknown"},
    "TM770": {"patient_id": "Visser_TM770", "treatment": "naive"},
}

_PATIENT_RE = re.compile(r"(TM\d+)")


def _parse_sample(filepath: Path) -> pd.DataFrame:
    """Read one CEL-Seq2 .txt.gz file → DataFrame (genes × cells)."""
    with gzip.open(filepath, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", index_col=0)
    df.index = [g.split("__")[-1] if "__" in g else g for g in df.index]
    df = df.loc[df.index != "UNK"] if "UNK" in df.index else df
    return df


def load_visser2023(data_dir: str | Path) -> ad.AnnData:
    """Load all Visser et al. 2023 CEL-Seq2 samples into a single AnnData.

    Parameters
    ----------
    data_dir
        Path to ``data/raw/Visser2023/`` containing ``*.transcripts.txt.gz`` files.

    Returns
    -------
    Concatenated AnnData with standardized obs fields.
    Raw counts stored in both ``adata.X`` and ``adata.layers["counts"]``.
    """
    data_dir = Path(data_dir)
    files = sorted(data_dir.glob("*.transcripts.txt.gz"))
    if not files:
        raise FileNotFoundError(f"No .transcripts.txt.gz files in {data_dir}")

    adatas = []

    for fp in files:
        m = _PATIENT_RE.search(fp.name)
        patient_code = m.group(1) if m else fp.stem
        meta = _PATIENT_MAP.get(patient_code, {
            "patient_id": f"Visser_{patient_code}",
            "treatment": "unknown",
        })

        df = _parse_sample(fp)

        # Columns are well positions (A1–P24); drop the UNK aggregate column
        df = df.drop(columns=["UNK"], errors="ignore")

        # Remove all-zero wells (empty wells on plate)
        df = df.loc[:, df.sum(axis=0) > 0]

        mat = sp.csr_matrix(df.values.T.astype(np.float32))
        cell_ids = [f"{patient_code}_{well}" for well in df.columns]
        gene_ids = df.index.tolist()

        obs = pd.DataFrame(
            {
                "sample_id": patient_code,
                "patient_id": meta["patient_id"],
                "treatment": meta["treatment"],
                "dataset": "Visser2023",
                "platform": "CELSeq2",
                "tissue": "primary_tumor",
            },
            index=cell_ids,
        )
        var = pd.DataFrame(index=gene_ids)

        adata = ad.AnnData(X=mat, obs=obs, var=var)
        adata.layers["counts"] = adata.X.copy()
        adatas.append(adata)

    combined = ad.concat(adatas, join="outer", fill_value=0)
    combined.var_names_make_unique()
    return combined
