import scanpy as sc
import scib
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, help="Path to input file (h5ad)")
parser.add_argument("-out", type=str, help="Path to output file (h5ad)")
parser.add_argument(
    "-emb",
    action="store_true",
    help="Whether input file contains integrated data or not",
)
parser.add_argument(
    "-dmso",
    action="store_true",
    help="Whether to subset the input file to the DMSO samples",
)

args = parser.parse_args()

# read adata object
adata = sc.read_h5ad(args.i)

print(adata)

# subset to DMSO if wanted
if args.dmso:
    print("Subset to only DMSO samples")
    adata = adata[adata.obs.treatment == "DMSO"].copy()
    print(adata)

# perform processing
print("### start processing ###")
if not args.emb:
    print("start PCA")
    sc.pp.pca(adata)
    print("start neighbor computation")
    sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca")
else:
    print("start neighbor computation")
    sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_emb")
print("start clustering")
scib.me.cluster_optimal_resolution(adata, cluster_key="cluster", label_key="moa")
print("### processing done ###")

# save result
print(adata)
print(adata.obs.cluster)
adata.write_h5ad(args.out)
