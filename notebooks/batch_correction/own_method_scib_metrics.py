import scanpy as sc
from bcc.scib_metrics import evaluate, metrics_preparation
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "-data",
    type=str,
    default="batch_correction_cellpainting/data/jump/",
    help="Path to the data directory",
)
parser.add_argument(
    "--compute_kbet",
    action="store_true",
    default=False,
    help="When set, the kBET scores are computed",
)
args = parser.parse_args()


### JUMP data
adata_unintegrated = sc.read_h5ad(f"{args.data}jump_target2_spherized.h5ad")
adata_unintegrated = metrics_preparation(adata_unintegrated, "Metadata_JCP2022")

# own method on JUMP, scanvi setting, correct for source
adata_integrated = sc.read_h5ad(f"{args.data}jump_sourcecorrected_scanviVAE.h5ad")
adata_integrated = metrics_preparation(adata_integrated, "Metadata_JCP2022", True)
df = evaluate(
    adata_unintegrated,
    adata_integrated,
    "Metadata_JCP2022",
    "Metadata_Source",
    "Own_scanvi",
    compute_kbet=args.compute_kbet,
)
df.to_csv(f"{args.data}scib_own_scanvi_jump_source.csv")

# own method on JUMP, scanvi setting, correct for plate
adata_integrated = sc.read_h5ad(f"{args.data}jump_platecorrected_scanviVAE.h5ad")
adata_integrated = metrics_preparation(adata_integrated, "Metadata_JCP2022", True)
df = evaluate(
    adata_unintegrated,
    adata_integrated,
    "Metadata_JCP2022",
    "Metadata_Plate",
    "Own_scanvi",
    compute_kbet=args.compute_kbet,
)
df.to_csv(f"{args.data}scib_own_scanvi_jump_plate.csv")

# own method on JUMP, scvi setting, correct for plate
adata_integrated = sc.read_h5ad(f"{args.data}jump_platecorrected_scviVAE.h5ad")
adata_integrated.obsm["X_emb"] = adata_integrated.X
adata_integrated = metrics_preparation(adata_integrated, "Metadata_JCP2022", True)
df = evaluate(
    adata_unintegrated,
    adata_integrated,
    "Metadata_JCP2022",
    "Metadata_Plate",
    "Own_scvi",
    compute_kbet=args.compute_kbet,
)
df.to_csv(f"{args.data}scib_own_scvi_jump_plate.csv")

# own method on JUMP, scvi setting, correct for source
adata_integrated = sc.read_h5ad(f"{args.data}jump_sourcecorrected_scviVAE.h5ad")
adata_integrated.obsm["X_emb"] = adata_integrated.X
adata_integrated = metrics_preparation(adata_integrated, "Metadata_JCP2022", True)
df = evaluate(
    adata_unintegrated,
    adata_integrated,
    "Metadata_JCP2022",
    "Metadata_Source",
    "Own_scvi",
    compute_kbet=args.compute_kbet,
)
df.to_csv(f"{args.data}scib_own_scvi_jump_source.csv")


### rxrx1
adata_unintegrated = sc.read_h5ad(f"{args.data}rxrx1_integrated_CVAE.h5ad")
adata_unintegrated = metrics_preparation(adata_unintegrated, "sirna")

# own method on rxrx1, scvi setting, correct for experiment
adata_integrated = sc.read_h5ad(f"{args.data}rxrx1_integrated_CVAE.h5ad")
adata_integrated = metrics_preparation(adata_integrated, "sirna", True)
df = evaluate(
    adata_unintegrated,
    adata_integrated,
    "sirna",
    "experiment",
    "Own_scvi",
    compute_kbet=args.compute_kbet,
)
df.to_csv(f"{args.data}scib_own_CVAE_rxrx1.csv")

# own method on rxrx1, scanvi setting, correct for experiment
adata_integrated = sc.read_h5ad(f"{args.data}rxrx1_integrated_scanviVAE.h5ad")
adata_integrated = metrics_preparation(adata_integrated, "sirna", True)
df = evaluate(
    adata_unintegrated,
    adata_integrated,
    "sirna",
    "experiment",
    "Own_scanvi",
    compute_kbet=args.compute_kbet,
)
df.to_csv(f"{args.data}scib_own_scanviVAE_rxrx1.csv")


### rxrx19b
adata_unintegrated = sc.read_h5ad(f"{args.data}rxrx19b_unintegrated.h5ad")
adata_unintegrated.obs.plate = adata_unintegrated.obs.plate.astype("category")
adata_unintegrated = metrics_preparation(adata_unintegrated, "disease_condition")

# own method on rxrx19b, scvi setting, correct for plate
adata_integrated = sc.read_h5ad(f"{args.data}rxrx19b_platecorrected_scviVAE.h5ad")
adata_integrated.obs.plate = adata_integrated.obs.plate.astype("category")
adata_integrated.obsm["X_emb"] = adata_integrated.X
adata_integrated = metrics_preparation(adata_integrated, "disease_condition", True)
df = evaluate(
    adata_unintegrated,
    adata_integrated,
    "disease_condition",
    "plate",
    "Own_scvi",
    compute_kbet=args.compute_kbet,
)
df.to_csv(f"{args.data}scib_own_scvi_rxrx19b.csv")

# own method on rxrx19b, scanvi setting, correct for plate
adata_integrated = sc.read_h5ad(f"{args.data}rxrx19b_platecorrected_scanviVAE.h5ad")
adata_integrated.obs.plate = adata_integrated.obs.plate.astype("category")
adata_integrated.obsm["X_emb"] = adata_integrated.X
adata_integrated = metrics_preparation(adata_integrated, "disease_condition", True)
df = evaluate(
    adata_unintegrated,
    adata_integrated,
    "disease_condition",
    "plate",
    "Own_scanvi",
    compute_kbet=args.compute_kbet,
)
df.to_csv(f"{args.data}scib_own_scanvi_rxrx19b.csv")
