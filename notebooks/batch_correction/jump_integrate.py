import scanpy as sc
from bcc.batch_correction import (
    scgen_integration,
    scanvi_integration,
    scvi_integration,
)
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    type=str,
    required=True,
    help="Path to the spherized JUMP file",
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
    default="scgen",
    choices=["scgen", "scanvi", "scvi"],
    help="Integration method that should be used (scgen, scanvi or scvi)",
)
args = parser.parse_args()


adata = sc.read_h5ad(args.i)

if args.method == "scgen":
    adata_overall = scgen_integration(
        adata, batch="Metadata_Source", labels="Metadata_JCP2022"
    )
elif args.method == "scanvi":
    adata_overall = scanvi_integration(
        adata, batch="Metadata_Source", labels="Metadata_JCP2022"
    )
else:
    adata_overall = scvi_integration(adata, batch="Metadata_Source")

# save output
adata_overall.write_h5ad(args.o)
