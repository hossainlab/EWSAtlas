"""Loader for Goodspeed et al. 2025 — 10x Chromium scRNA-seq.

STATUS: BLOCKED — GEO accession GSE176190 contains Affymetrix .CEL.gz files,
not 10x scRNA-seq matrices. The paper describes 10x Chromium (~6,228 cells,
7 treatment-naive patients, primary tumors + CTCs).

Resolution steps:
1. Check GSE176190 for separate scRNA-seq sub-accessions (GSM with
   matrix.mtx / barcodes.tsv / features.tsv)
2. If absent, contact authors (Goodspeed lab) for data deposit or transfer

Once 10x matrices are available, place per-sample directories under
data/raw/Goodspeed2025/<sample_id>/ and update _SAMPLE_METADATA below.
"""

from pathlib import Path

import anndata as ad


_SAMPLE_METADATA: dict[str, dict] = {
    # Populate once data is available.
    # Example entry:
    # "GS_PT1": {"patient_id": "Goodspeed_PT1", "treatment": "naive", "tissue": "primary_tumor"},
    # "GS_CTC1": {"patient_id": "Goodspeed_PT1", "treatment": "naive", "tissue": "CTC"},
}


def load_goodspeed2025(data_dir: str | Path) -> ad.AnnData:
    """Load Goodspeed et al. 2025 10x Chromium samples.

    Parameters
    ----------
    data_dir
        Path to ``data/raw/Goodspeed2025/`` with per-sample MTX directories.

    Raises
    ------
    NotImplementedError
        Until scRNA-seq matrices are obtained from GEO or authors.
    """
    raise NotImplementedError(
        "Goodspeed 2025 data not yet available in scRNA-seq format. "
        "GSE176190 currently contains Affymetrix .CEL.gz files. "
        "See module docstring for resolution steps."
    )
