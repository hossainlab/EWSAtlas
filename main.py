"""EWSAtlas — entry point for the full harmonized atlas pipeline."""

from pathlib import Path

import yaml


def load_config(path: str = "configs/params.yaml") -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def main() -> None:
    cfg = load_config()
    raw = Path(cfg["paths"]["raw_data"])
    processed = Path(cfg["paths"]["processed_data"])
    integrated = Path(cfg["paths"]["integrated_data"])
    processed.mkdir(parents=True, exist_ok=True)
    integrated.mkdir(parents=True, exist_ok=True)

    from ewsatlas.io import load_he2025, load_visser2023
    from ewsatlas.qc import (
        apply_qc_filters,
        compute_qc_metrics,
        detect_doublets,
        flag_outliers,
    )
    from ewsatlas.integration import run_harmony, run_scvi
    from ewsatlas.annotation import annotate_celltypes, score_ews_markers

    qc_cfg = cfg["qc"]
    int_cfg = cfg["integration"]

    # ── Phase 1: Load ────────────────────────────────────────────────────────
    print("Loading He et al. 2025 ...")
    he = load_he2025(raw / "He2025")
    he.write_h5ad(processed / "he2025_raw.h5ad")

    print("Loading Visser et al. 2023 ...")
    visser = load_visser2023(raw / "Visser2023")
    visser.write_h5ad(processed / "visser2023_raw.h5ad")

    # ── Phase 2: QC ──────────────────────────────────────────────────────────
    datasets = {"he2025": he, "visser2023": visser}
    qc_passed = {}

    for name, adata in datasets.items():
        print(f"\nQC: {name}")
        adata = compute_qc_metrics(adata, mt_prefix=qc_cfg["mt_gene_prefix"])
        adata = detect_doublets(adata)
        adata = adata[~adata.obs["predicted_doublet"]].copy()
        adata = flag_outliers(adata, n_mads=qc_cfg["mad_n_mads"])
        adata = apply_qc_filters(
            adata,
            min_cells_per_gene=qc_cfg["min_cells_per_gene"],
            report_path=processed / f"{name}_qc_report.csv",
        )
        adata.write_h5ad(processed / f"{name}_qc.h5ad")
        qc_passed[name] = adata

    # ── Phase 3: Integration ──────────────────────────────────────────────────
    import anndata as ad

    print("\nConcatenating datasets ...")
    combined = ad.concat(
        list(qc_passed.values()),
        join="inner",
        label="dataset",
        keys=list(qc_passed.keys()),
    )
    combined.layers["counts"] = combined.X.copy()

    print("Running Harmony (sanity check) ...")
    adata_harmony = run_harmony(
        combined,
        batch_key=int_cfg["batch_key"],
        n_hvgs=int_cfg["n_hvgs"],
    )
    adata_harmony.write_h5ad(integrated / "ewsatlas_harmony.h5ad")

    print("Running scVI (final integration) ...")
    adata_scvi, model = run_scvi(
        combined,
        batch_key=int_cfg["batch_key"],
        n_hvgs=int_cfg["n_hvgs"],
        n_latent=int_cfg["scvi_n_latent"],
        n_layers=int_cfg["scvi_n_layers"],
        max_epochs=int_cfg["scvi_max_epochs"],
        early_stopping=int_cfg["scvi_early_stopping"],
    )
    model.save(str(integrated / "scvi_model"), overwrite=True)
    adata_scvi.write_h5ad(integrated / "ewsatlas_scvi.h5ad")

    # ── Phase 4: Annotation ───────────────────────────────────────────────────
    print("Annotating cell types ...")
    adata_scvi = score_ews_markers(adata_scvi)
    adata_scvi.write_h5ad(integrated / "ewsatlas_annotated.h5ad")
    print("Done. Atlas saved →", integrated / "ewsatlas_annotated.h5ad")


if __name__ == "__main__":
    main()
