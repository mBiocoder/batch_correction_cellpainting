import scanpy as sc
import scib
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i", type=str, help="Path to the input h5ad data file to integrate"
)
parser.add_argument(
    "-o", type=str, help="Path to the h5ad output file (integrated data set)"
)
parser.add_argument(
    "-batch", type=str, default="week", help="Batch key (techincal variation to remove)"
)
parser.add_argument(
    "-label", type=str, default="moa", help="Label (biological variation to retain)"
)
parser.add_argument(
    "-max_epochs",
    type=int,
    default=400,
    help="Maximum number of epochs for the integration",
)

args = parser.parse_args()

# read input file
adata = sc.read_h5ad(args.i)

# integration with scanVI
scanvi = adata.copy()
scanvi.layers["counts"] = scanvi.X
scib.ig.scanvi(scanvi, batch=args.batch, labels=args.label, max_epochs=args.max_epochs)

# save result
scanvi.write_h5ad(args.o)
