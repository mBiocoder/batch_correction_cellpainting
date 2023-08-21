import scanpy as sc
from bcc.scib_metrics import evaluate_low_level, evaluate, metrics_preparation
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "-data",
    type=str,
    default="batch_correction_cellpainting/data/jump/",
    help="Path to the data directory",
)
args = parser.parse_args()

adata_unintegrated = sc.read_h5ad(f"{args.data}jump_target2_spherized.h5ad")

# metrics computation low level
adata_integrated = sc.read_h5ad(f"{args.data}harmony_low.h5ad")
df = evaluate_low_level(
    adata_unintegrated,
    adata_integrated,
    "Metadata_JCP2022",
    "Metadata_Plate",
    "Metadata_Source",
    "Harmony",
)
df.to_csv(f"{args.data}scib_harmony_low.csv")

adata_integrated = sc.read_h5ad(f"{args.data}scanorama_low.h5ad")
df = evaluate_low_level(
    adata_unintegrated,
    adata_integrated,
    "Metadata_JCP2022",
    "Metadata_Plate",
    "Metadata_Source",
    "Scanorama",
    compute_unintegrated=False,
)
df.to_csv(f"{args.data}scib_scanorama_low.csv")

adata_integrated = sc.read_h5ad(f"{args.data}scgen_low.h5ad")
df = evaluate_low_level(
    adata_unintegrated,
    adata_integrated,
    "Metadata_JCP2022",
    "Metadata_Plate",
    "Metadata_Source",
    "scGen",
    compute_unintegrated=False,
    integrated_obsm_key="corrected_latent",
)
df.to_csv(f"{args.data}scib_scgen_low.csv")

adata_integrated = sc.read_h5ad(f"{args.data}scanvi_low.h5ad")
df = evaluate_low_level(
    adata_unintegrated,
    adata_integrated,
    "Metadata_JCP2022",
    "Metadata_Plate",
    "Metadata_Source",
    "scanVI",
    compute_unintegrated=False,
    integrated_obsm_key="X_emb",
)
df.to_csv(f"{args.data}scib_scanvi_low.csv")


# metrics computation high level
adata_unintegrated = metrics_preparation(adata_unintegrated, "Metadata_JCP2022")

adata_integrated = sc.read_h5ad(f"{args.data}harmony_high.h5ad")
adata_integrated = metrics_preparation(adata_integrated, "Metadata_JCP2022", True)
df = evaluate(
    adata_unintegrated,
    adata_integrated,
    "Metadata_JCP2022",
    "Metadata_Source",
    "Harmony",
)
df.to_csv(f"{args.data}scib_harmony_high.csv")

adata_integrated = sc.read_h5ad(f"{args.data}scanorama_high.h5ad")
adata_integrated = metrics_preparation(
    adata_integrated, "Metadata_JCP2022", True, integrated_obsm_key="X_emb"
)
df = evaluate(
    adata_unintegrated,
    adata_integrated,
    "Metadata_JCP2022",
    "Metadata_Source",
    "Scanorama",
    compute_unintegrated=False,
)
df.to_csv(f"{args.data}scib_scanorama_high.csv")
