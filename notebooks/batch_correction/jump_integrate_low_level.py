import scanpy as sc
from bcc.batch_correction import (
    scanvi_integration,
    scvi_integration,
)
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    type=str,
    required=True,
    help="Path to the unintegrated JUMP file (spherized)",
)
parser.add_argument(
    "-o",
    type=str,
    required=True,
    help="Path to where the integrated .h5ad file should be saved",
)
parser.add_argument(
    "-method",
    type=str,
    default="scanvi",
    choices=["scanvi", "scvi"],
    help="Integration method that should be used (scgen, scanvi or scvi)",
)
args = parser.parse_args()


adata = sc.read_h5ad(args.i)

if args.method == "scanvi":
    adata_per_source = scanvi_integration(
        adata,
        batch="Metadata_Plate",
        labels="Metadata_JCP2022",
        hierarchical="Metadata_Source",
    )
else:
    adata_per_source = scvi_integration(
        adata,
        batch="Metadata_Plate",
        hierarchical="Metadata_Source",
    )

# save output
adata_per_source.write_h5ad(args.o)
