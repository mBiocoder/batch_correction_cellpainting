import scanpy as sc
from bcc.scib_metrics import evaluate, metrics_preparation


adata_unintegrated = sc.read_h5ad(
    f"batch_correction_cellpainting/data/jump/jump_target2_spherized.h5ad"
)

adata_integrated = sc.read_h5ad(f"batch_correction_cellpainting/data/jump/scanvi.h5ad")

adata_unintegrated = adata_unintegrated[
    adata_unintegrated.obs.Metadata_Source == "source_9"
].copy()
adata_unintegrated = metrics_preparation(adata_unintegrated, "Metadata_JCP2022")

adata_integrated = adata_integrated[
    adata_integrated.obs.Metadata_Source == "source_9"
].copy()
adata_integrated = metrics_preparation(adata_integrated, "Metadata_JCP2022", True)

df = evaluate(
    adata_unintegrated,
    adata_integrated,
    "Metadata_JCP2022",
    "Metadata_Plate",
    "Unintegrated",
    compute_unintegrated=False,
)
print("Computation successful!")
