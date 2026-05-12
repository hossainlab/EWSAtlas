import pandas as pd
from pathlib import Path

def create_sample_metadata():
    samples = []

    # ── Visser et al. 2023 ──────────────────────────────────────────────────
    # Based on GEO GSE243347 and Cancer Res Commun 2023
    visser_mapping = {
        "TM338": {"patient_id": "Visser_ES024", "treatment": "post-treatment", "tissue": "primary_tumor", "sample_type": "resection", "location": "pelvis"},
        "TM339": {"patient_id": "Visser_ES024", "treatment": "post-treatment", "tissue": "primary_tumor", "sample_type": "resection", "location": "pelvis"},
        "TM344": {"patient_id": "Visser_ES016", "treatment": "post-treatment", "tissue": "primary_tumor", "sample_type": "resection", "location": "lung"},
        "TM348": {"patient_id": "Visser_ES016", "treatment": "post-treatment", "tissue": "primary_tumor", "sample_type": "resection", "location": "lung"},
        "TM416": {"patient_id": "Visser_ES025", "treatment": "pre-treatment", "tissue": "primary_tumor", "sample_type": "biopsy", "location": "tibia"},
        "TM417": {"patient_id": "Visser_ES025", "treatment": "pre-treatment", "tissue": "primary_tumor", "sample_type": "biopsy", "location": "tibia"},
        "TM424": {"patient_id": "Visser_ES016", "treatment": "pre-treatment", "tissue": "primary_tumor", "sample_type": "biopsy", "location": "lung"},
        "TM425": {"patient_id": "Visser_ES016", "treatment": "pre-treatment", "tissue": "metastasis", "sample_type": "biopsy", "location": "lung"},
        "TM505": {"patient_id": "Visser_ES027", "treatment": "pre-treatment", "tissue": "primary_tumor", "sample_type": "biopsy", "location": "femur"},
        "TM506": {"patient_id": "Visser_ES027", "treatment": "pre-treatment", "tissue": "primary_tumor", "sample_type": "biopsy", "location": "femur"},
        "TM547": {"patient_id": "Visser_ES030", "treatment": "post-treatment", "tissue": "primary_tumor", "sample_type": "resection", "location": "pelvis"},
        "TM548": {"patient_id": "Visser_ES030", "treatment": "post-treatment", "tissue": "primary_tumor", "sample_type": "resection", "location": "pelvis"},
        "TM549": {"patient_id": "Visser_ES030", "treatment": "post-treatment", "tissue": "primary_tumor", "sample_type": "resection", "location": "pelvis"},
        "TM552": {"patient_id": "Visser_ES006", "treatment": "post-treatment", "tissue": "metastasis", "sample_type": "resection", "location": "lung nodule"},
        "TM564": {"patient_id": "Visser_ES036", "treatment": "pre-treatment", "tissue": "primary_tumor", "sample_type": "biopsy", "location": "rib"},
        "TM570": {"patient_id": "Visser_ES036", "treatment": "pre-treatment", "tissue": "metastasis", "sample_type": "biopsy", "location": "lung"},
        "TM572": {"patient_id": "Visser_ES006", "treatment": "post-treatment", "tissue": "metastasis", "sample_type": "resection", "location": "lung nodule"},
        "TM574": {"patient_id": "Visser_ES039", "treatment": "pre-treatment", "tissue": "primary_tumor", "sample_type": "biopsy", "location": "fibula"},
        "TM707": {"patient_id": "Visser_ES042", "treatment": "post-treatment", "tissue": "primary_tumor", "sample_type": "resection", "location": "pelvis"},
        "TM709": {"patient_id": "Visser_ES042", "treatment": "post-treatment", "tissue": "primary_tumor", "sample_type": "resection", "location": "pelvis"},
        "TM737": {"patient_id": "Visser_ES042", "treatment": "post-treatment", "tissue": "primary_tumor", "sample_type": "resection", "location": "pelvis"},
        "TM739": {"patient_id": "Visser_ES042", "treatment": "post-treatment", "tissue": "primary_tumor", "sample_type": "resection", "location": "pelvis"},
        "TM712": {"patient_id": "Visser_ES010", "treatment": "post-treatment", "tissue": "primary_tumor", "sample_type": "resection", "location": "tibia"},
        "TM734": {"patient_id": "Visser_ES162", "treatment": "post-treatment", "tissue": "primary_tumor", "sample_type": "resection", "location": "pelvis"},
        "TM736": {"patient_id": "Visser_ES039_ES030_Mixed", "treatment": "mixed", "tissue": "mixed", "sample_type": "mixed", "location": "mixed"},
        "TM768": {"patient_id": "Visser_ES048", "treatment": "post-treatment", "tissue": "metastasis", "sample_type": "resection", "location": "lung"},
        "TM770": {"patient_id": "Visser_ES039", "treatment": "post-treatment", "tissue": "primary_tumor", "sample_type": "resection", "location": "fibula"},
    }

    for sample_id, meta in visser_mapping.items():
        samples.append({
            "sample_id": sample_id,
            "patient_id": meta["patient_id"],
            "dataset": "Visser2023",
            "treatment": meta["treatment"],
            "tissue": meta["tissue"],
            "sample_type": meta["sample_type"],
            "location": meta["location"],
            "platform": "CELSeq2",
            "fusion": "EWS-FLI1"
        })

    # ── He et al. 2025 ──────────────────────────────────────────────────────
    # Based on Cell Comm Sig 2025 and sample IDs
    he_mapping = {
        "TN.1": {"patient_id": "He_TN1", "treatment": "pre-treatment"},
        "TN.2": {"patient_id": "He_TN2", "treatment": "pre-treatment"},
        "TN.3": {"patient_id": "He_TN3", "treatment": "pre-treatment"},
        "NAC.1": {"patient_id": "He_NAC1", "treatment": "post-treatment"},
        "NAC.3": {"patient_id": "He_NAC3", "treatment": "post-treatment"},
        "RC.1": {"patient_id": "He_RC1", "treatment": "relapsed"},
        "RC.2": {"patient_id": "He_RC2", "treatment": "relapsed"},
        "RC.3": {"patient_id": "He_RC3", "treatment": "relapsed"},
    }

    for sample_id, meta in he_mapping.items():
        samples.append({
            "sample_id": sample_id,
            "patient_id": meta["patient_id"],
            "dataset": "He2025",
            "treatment": meta["treatment"],
            "tissue": "primary_tumor",
            "sample_type": "unknown",
            "location": "unknown",
            "platform": "10x_Chromium",
            "fusion": "unknown"
        })

    # ── Goodspeed et al. 2025 (GSE277083) ──────────────────────────────────
    # Note: Currently missing from raw/ folder (folder contains cell line data instead)
    # Adding for completeness of metadata
    goodspeed_patients = ["001", "002", "005", "006", "009", "010", "012", "038", "051", "061", "066"]
    for p in goodspeed_patients:
        # Primary Tumor
        samples.append({
            "sample_id": f"Goodspeed_{p}_PT",
            "patient_id": f"Goodspeed_{p}",
            "dataset": "Goodspeed2025",
            "treatment": "pre-treatment",
            "tissue": "primary_tumor",
            "sample_type": "biopsy",
            "location": "unknown",
            "platform": "10x_Chromium",
            "fusion": "unknown"
        })
        # CTCs (where available)
        samples.append({
            "sample_id": f"Goodspeed_{p}_CTC",
            "patient_id": f"Goodspeed_{p}",
            "dataset": "Goodspeed2025",
            "treatment": "pre-treatment",
            "tissue": "peripheral_blood",
            "sample_type": "CTC",
            "location": "unknown",
            "platform": "10x_Chromium",
            "fusion": "unknown"
        })

    df = pd.DataFrame(samples)
    df.to_csv("data/metadata/sample_metadata.csv", index=False)
    print(f"Sample metadata created with {len(df)} entries.")

def create_patient_metadata():
    patients = []

    # ── Visser et al. 2023 ──────────────────────────────────────────────────
    visser_patients = [
        {"patient_id": "Visser_ES006", "age": 13.7, "sex": "M", "disease_status": "metastatic"},
        {"patient_id": "Visser_ES016", "age": 12.0, "sex": "M", "disease_status": "metastatic"},
        {"patient_id": "Visser_ES024", "age": 14.0, "sex": "M", "disease_status": "localized"},
        {"patient_id": "Visser_ES025", "age": 15.0, "sex": "M", "disease_status": "localized"},
        {"patient_id": "Visser_ES027", "age": 17.0, "sex": "M", "disease_status": "localized"},
        {"patient_id": "Visser_ES030", "age": 16.0, "sex": "M", "disease_status": "metastatic"},
        {"patient_id": "Visser_ES036", "age": 11.0, "sex": "M", "disease_status": "metastatic"},
        {"patient_id": "Visser_ES039", "age": 13.0, "sex": "F", "disease_status": "localized"},
        {"patient_id": "Visser_ES042", "age": 10.0, "sex": "M", "disease_status": "localized"},
        {"patient_id": "Visser_ES010", "age": 6.0, "sex": "M", "disease_status": "localized"},
        {"patient_id": "Visser_ES162", "age": 14.0, "sex": "M", "disease_status": "localized"},
        {"patient_id": "Visser_ES048", "age": 15.0, "sex": "M", "disease_status": "metastatic"},
    ]
    for p in visser_patients:
        p.update({"dataset": "Visser2023", "fusion": "EWS-FLI1"})
        patients.append(p)

    # ── He et al. 2025 ──────────────────────────────────────────────────────
    he_patients = ["He_TN1", "He_TN2", "He_TN3", "He_NAC1", "He_NAC3", "He_RC1", "He_RC2", "He_RC3"]
    for p in he_patients:
        patients.append({
            "patient_id": p,
            "age": None,
            "sex": "unknown",
            "disease_status": "unknown",
            "dataset": "He2025",
            "fusion": "unknown"
        })

    # ── Goodspeed et al. 2025 ──────────────────────────────────────────────
    goodspeed_patients = ["001", "002", "005", "006", "009", "010", "012", "038", "051", "061", "066"]
    for p in goodspeed_patients:
        patients.append({
            "patient_id": f"Goodspeed_{p}",
            "age": None,
            "sex": "unknown",
            "disease_status": "unknown",
            "dataset": "Goodspeed2025",
            "fusion": "unknown"
        })

    df = pd.DataFrame(patients)
    df.to_csv("data/metadata/patient_metadata.csv", index=False)
    print(f"Patient metadata created with {len(df)} entries.")

if __name__ == "__main__":
    create_sample_metadata()
    create_patient_metadata()
