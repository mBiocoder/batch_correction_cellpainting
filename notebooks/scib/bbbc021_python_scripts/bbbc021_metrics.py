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

# compute scIB metrics
print("### Start metrics computation ###")

# biological conservation
print("ARI...")
ari_pre = scib.me.ari(unintegrated, cluster_key="cluster", label_key="moa")
ari_post = scib.me.ari(integrated, cluster_key="cluster", label_key="moa")

print("cLISI...")
clisi_pre = scib.me.clisi_graph(unintegrated, label_key="moa", type_="full")
clisi_post = scib.me.clisi_graph(
    integrated, label_key="moa", type_="embed", use_rep="X_emb"
)

print("Isolated Labels ASW...")
il_asw_pre = scib.me.isolated_labels_asw(
    unintegrated, label_key="moa", batch_key="week", embed="X_pca"
)
il_asw_post = scib.me.isolated_labels_asw(
    integrated, label_key="moa", batch_key="week", embed="X_emb"
)

print("Isolated Labels F1...")
il_f1_pre = scib.me.isolated_labels_f1(
    unintegrated, batch_key="week", label_key="moa", embed=None
)
il_f1_post = scib.me.isolated_labels_f1(
    integrated, batch_key="week", label_key="moa", embed=None
)

print("NMI...")
nmi_pre = scib.me.nmi(unintegrated, cluster_key="cluster", label_key="moa")
nmi_post = scib.me.nmi(integrated, cluster_key="cluster", label_key="moa")

print("Silhouette...")
silhouette_pre = scib.me.silhouette(unintegrated, label_key="moa", embed="X_pca")
silhouette_post = scib.me.silhouette(integrated, label_key="moa", embed="X_emb")

# batch correction
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
            "ARI",
            "cLISI",
            "Isolated_labels_ASW",
            "Isolated_labels_F1",
            "NMI",
            "Silhouette",
            "Graph_connectivity",
            "iLISI",
            "kBET",
            "PCR_comparison",
            "Batch_ASW",
        ],
        "Unintegrated": [
            ari_pre,
            clisi_pre,
            il_asw_pre,
            il_f1_pre,
            nmi_pre,
            silhouette_pre,
            graph_conn_pre,
            ilisi_pre,
            kbet_pre,
            pcr_pre,
            basw_pre,
        ],
        args.method_name: [
            ari_post,
            clisi_post,
            il_asw_post,
            il_f1_post,
            nmi_post,
            silhouette_post,
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
