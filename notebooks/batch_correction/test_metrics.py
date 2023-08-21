import scanpy as sc
from bcc.scib_metrics import evaluate_low_level
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "-unintegrated",
    type=str,
    required=True,
    help="Path to the unintegrated JUMP file (spherized)",
)
parser.add_argument(
    "-integrated",
    type=str,
    required=True,
    help="Path to the low level Harmony integration",
)
args = parser.parse_args()

unintegrated = sc.read_h5ad(args.unintegrated)
integrated = sc.read_h5ad(args.integrated)

# for testing: use only one source
unintegrated = unintegrated[unintegrated.obs.Metadata_Source == "source_9"].copy()
integrated = integrated[integrated.obs.Metadata_Source == "source_9"].copy()

df = evaluate_low_level(
    unintegrated,
    integrated,
    "Metadata_JCP2022",
    "Metadata_Plate",
    "Metadata_Source",
    "Harmony",
)
print("Finished metrics computation, everything worked :)")
