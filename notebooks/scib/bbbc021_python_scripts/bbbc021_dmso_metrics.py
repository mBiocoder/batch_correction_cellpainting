import scanpy as sc
import scib
import argparse
import pandas as pd
import rpy2.robjects.packages as rpackages

kbet = rpackages.importr("kBET")

parser = argparse.ArgumentParser()
parser.add_argument(
    "-unintegrated", type=str, help="Path to the unintegrated h5ad file"
)
parser.add_argument(
    "-integrated",
    type=str,
    help="Path to the integrated h5ad file",
)
parser.add_argument(
    "-out",
    type=str,
    help="Path to output file",
)
parser.add_argument(
    "-method_name",
    type=str,
    default="Integtated",
    help="Name of the integration method, used as column name of the output table",
)

args = parser.parse_args()

# read files
unintegrated = sc.read_h5ad(args.unintegrated)
integrated = sc.read_h5ad(args.integrated)

# subset to DMSO samples
unintegrated = unintegrated[unintegrated.obs.treatment == "DMSO"]
integrated = integrated[integrated.obs.treatment == "DMSO"]

# compute scIB batch correction metrics
print("### Start metrics computation ###")

print("Graph connectivity...")
graph_conn_pre = scib.me.graph_connectivity(unintegrated, label_key="moa")
graph_conn_post = scib.me.graph_connectivity(integrated, label_key="moa")

print("iLISI...")
ilisi_pre = scib.me.ilisi_graph(unintegrated, batch_key="week", type_="full")
ilisi_post = scib.me.ilisi_graph(
    integrated, batch_key="week", type_="embed", use_rep="X_emb"
)

print("kBET...")
kbet_pre = scib.me.kBET(
    unintegrated, batch_key="week", label_key="moa", type_="full", embed="X_pca"
)
kbet_post = scib.me.kBET(
    integrated, batch_key="week", label_key="moa", type_="embed", embed="X_emb"
)

print("PC regression...")
pcr_pre = pd.NA
pcr_post = scib.me.pcr_comparison(
    unintegrated, integrated, covariate="week", embed="X_emb"
)

print("Batch ASW...")
basw_pre = scib.me.silhouette_batch(
    unintegrated, batch_key="week", label_key="moa", embed="X_pca"
)
basw_post = scib.me.silhouette_batch(
    integrated, batch_key="week", label_key="moa", embed="X_emb"
)

print("### Metrics computation done ###")


# create result data table
results = pd.DataFrame(
    {
        "Metric": [
            "Graph_connectivity",
            "iLISI",
            "kBET",
            "PCR_comparison",
            "Batch_ASW",
        ],
        "Unintegrated": [
            graph_conn_pre,
            ilisi_pre,
            kbet_pre,
            pcr_pre,
            basw_pre,
        ],
        args.method_name: [
            graph_conn_post,
            ilisi_post,
            kbet_post,
            pcr_post,
            basw_post,
        ],
    }
)

# save result
results.to_csv(args.out, index=False)
