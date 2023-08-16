import scanpy as sc
import scib
import pandas as pd
import rpy2.robjects.packages as rpackages
from anndata import AnnData

kbet = rpackages.importr("kBET")


def evaluate_methods(
    adata_unintegrated: AnnData,
    adatas_integrated: list(AnnData),
    method_names: list(str),
    label_key: str,
    batch_key: str,
) -> pd.DataFrame:
    """Evaluate different integration methods by computing the scIB metrics for each

    Results are returned in one combined DataFrame

    Args:
        adata_unintegrated: AnnData object containing unintegrated data
        adatas_integrated:  List containing the integrated AnnData objects
                            of the different integration methods to evaluate
        method_names:       List containing the names of the integration
                            methods used for creating `adatas_integrated`.
                            Should not contain duplicates!
        label_key:          Column name of the biological labels
        batch_key:          Column name of the batch labels

    Returns:
        DataFrame:          scIB metrics for each method saved in a DataFrame
    """
    # prepare AnnData objects
    adata_unintegrated = _metrics_preparation(adata_unintegrated, label_key)
    adatas_integrated = [
        _metrics_preparation(adata, label_key, integrated=True)
        for adata in adatas_integrated
    ]

    # compute metrics
    result_frames = []
    for i, adata in enumerate(adatas_integrated):
        if i == 0:
            result_frames.append(
                _evaluate(
                    adata_unintegrated, adata, label_key, batch_key, method_names[i]
                )
            )
        else:
            result_frames.append(
                _evaluate(
                    adata_unintegrated,
                    adata,
                    label_key,
                    batch_key,
                    method_names[i],
                    compute_unintegrated=False,
                )
            )

    # combine resulting DataFrames
    result_frames = [df.set_index("Metric") for df in result_frames]
    return pd.concat(result_frames, axis=1)


def _evaluate(
    adata_unintegrated: AnnData,
    adata_integrated: AnnData,
    label_key: str,
    batch_key: str,
    integration_method: str,
    compute_unintegrated: bool = True,
) -> pd.DataFrame:
    """Compute all scIB metrics

    Compute all scIB metrics for unintegrated and integrated data
    and return the results in form of a pandas DataFrame.
    When specified, metrics are only computed for the integrated data.

    Args:
        adata_unintegrated:     AnnData object containing unintegrated data
                                Should already be prepared!
        adata_integrated:       AnnData object containing integrated data
                                Should already be prepared!
        label_key:              Column name of the biological labels
        batch_key:              Column name of the batch labels
        integration_method:     Name of the used integration method
        compute_unintegrated:   If False, skip metrics computation for
                                unintegrated data

    Returns:
        DataFrame:          All scIB metrics saved in a DataFrame
    """
    # compute all metrics for unintegrated and integrated data
    if compute_unintegrated:
        unintegrated_results = [
            _compute_ari(adata_unintegrated, label_key),
            _compute_clisi(adata_unintegrated, label_key),
            _compute_il_asw(adata_unintegrated, label_key, batch_key),
            _compute_il_f1(adata_unintegrated, label_key, batch_key),
            _compute_nmi(adata_unintegrated, label_key),
            _compute_silhouette(adata_unintegrated, label_key),
            _compute_graph_connectivity(adata_unintegrated, label_key),
            _compute_ilisi(adata_unintegrated, batch_key),
            _compute_kbet(adata_unintegrated, label_key, batch_key),
            _compute_pcr(adata_unintegrated, batch_key),
            _compute_batch_asw(adata_unintegrated, label_key, batch_key),
        ]
    integrated_results = [
        _compute_ari(adata_unintegrated, label_key, integrated=adata_integrated),
        _compute_clisi(adata_unintegrated, label_key, adata_integrated),
        _compute_il_asw(adata_unintegrated, label_key, batch_key, adata_integrated),
        _compute_il_f1(adata_unintegrated, label_key, batch_key, adata_integrated),
        _compute_nmi(adata_unintegrated, label_key, integrated=adata_integrated),
        _compute_silhouette(adata_unintegrated, label_key, adata_integrated),
        _compute_graph_connectivity(adata_unintegrated, label_key, adata_integrated),
        _compute_ilisi(adata_unintegrated, batch_key, adata_integrated),
        _compute_kbet(adata_unintegrated, label_key, batch_key, adata_integrated),
        _compute_pcr(adata_unintegrated, batch_key, adata_integrated),
        _compute_batch_asw(adata_unintegrated, label_key, batch_key, adata_integrated),
    ]

    # return result as pandas DataFrame
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
            integration_method: integrated_results,
        }
    )
    if compute_unintegrated:
        results["Unintegrated"] = unintegrated_results
    return results


def _metrics_preparation(
    adata: AnnData,
    label_key: str,
    integrated: bool = False,
    cluster_key: str = "cluster",
) -> AnnData:
    """Prepare an AnnData object for scIB metrics computation

    Neighbor computation and clustering is performed which is
    necessary for computing the scIB metrics.
    In case the input AnnData contains unintegrated data, PCA
    is performed first.

    Args:
        adata:          AnnData object containing the data
        label_key:      Column name of the biological labels
        integrated:     Whether the AnnData contains integrated
                        data or not
        cluster_key:    Key for storing the cluster labels

    Returns:
        AnnData:        Batch corrected AnnData object
    """
    if not integrated:
        sc.pp.pca(adata)
        sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca")
    else:
        sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_emb")
    scib.me.cluster_optimal_resolution(
        adata, cluster_key=cluster_key, label_key=label_key
    )
    return adata


def _compute_ari(
    unintegrated: AnnData,
    label_key: str,
    cluster_key: str = "cluster",
    integrated: AnnData = None,
) -> float:
    """Compute biological conservation metric ARI

    ARI = Adjusted Rand Index.
    The ARI compares the overlap of two clusterings, in this case
    the biological labels and the Louvain clustering.
    Can be applied to all integration output types.

    Args:
        unintegrated:   AnnData object with unintegrated data.
                        This is used for metric computation when
                        no integrated AnnData is given
        label_key:      Column name of the biological labels
        cluster_key:    Key that stores the cluster labels
        integrated:     AnnData object with integrated data.
                        When None, the unintegrated AnnData is
                        used for metric computation

    Returns:
        float:        ARI score in [0, 1]
    """
    if not isinstance(integrated, AnnData):
        return scib.me.ari(unintegrated, cluster_key=cluster_key, label_key=label_key)
    else:
        return scib.me.ari(integrated, cluster_key=cluster_key, label_key=label_key)


def _compute_clisi(
    unintegrated: AnnData,
    label_key: str,
    integrated: AnnData = None,
) -> float:
    """Compute biological conservation metric cLISI

    Quantifies how many items can be drawn from a
    neighbor list before one batch is observed twice.
    Can be applied to all integration output types.

    Args:
        unintegrated:   AnnData object with unintegrated data.
                        This is used for metric computation when
                        no integrated AnnData is given
        label_key:      Column name of the biological labels
        integrated:     AnnData object with integrated data.
                        When None, the unintegrated AnnData is
                        used for metric computation

    Returns:
        float:        cLISI score in [0, 1]
    """
    if not isinstance(integrated, AnnData):
        return scib.me.clisi_graph(unintegrated, label_key=label_key, type_="full")
    else:
        return scib.me.clisi_graph(
            integrated, label_key=label_key, type_="embed", use_rep="X_emb"
        )


def _compute_il_asw(
    unintegrated: AnnData,
    label_key: str,
    batch_key: str,
    integrated: AnnData = None,
) -> float:
    """Compute biological conservation metric isolated labels ASW

    Quantify how well isolated labels are distinguished from all
    other labels using the average-width silhouette score (ASW).
    TODO: output types?

    Args:
        unintegrated:   AnnData object with unintegrated data.
                        This is used for metric computation when
                        no integrated AnnData is given
        label_key:      Column name of the biological labels
        batch_key:      Column name of the batch labels
        integrated:     AnnData object with integrated data.
                        When None, the unintegrated AnnData is
                        used for metric computation

    Returns:
        float:        isolated labels ASW score in [0, 1]
    """
    if not isinstance(integrated, AnnData):
        return scib.me.isolated_labels_asw(
            unintegrated, label_key=label_key, batch_key=batch_key, embed="X_pca"
        )
    else:
        return scib.me.isolated_labels_asw(
            integrated, label_key=label_key, batch_key=batch_key, embed="X_emb"
        )


def _compute_il_f1(
    unintegrated: AnnData,
    label_key: str,
    batch_key: str,
    integrated: AnnData = None,
) -> float:
    """Compute biological conservation metric isolated labels F1

    Quanify how well isolated labels are distinguished from other
    labels by data-driven clustering.
    The F1 score is used to evaluate clustering w.r.t. ground truth labels.
    TODO: output types?

    Args:
        unintegrated:   AnnData object with unintegrated data.
                        This is used for metric computation when
                        no integrated AnnData is given
        label_key:      Column name of the biological labels
        batch_key:      Column name of the batch labels
        integrated:     AnnData object with integrated data.
                        When None, the unintegrated AnnData is
                        used for metric computation

    Returns:
        float:        isolated labels F1 score in [0, 1]
    """
    if not isinstance(integrated, AnnData):
        return scib.me.isolated_labels_f1(
            unintegrated, batch_key=batch_key, label_key=label_key, embed=None
        )
    else:
        return scib.me.isolated_labels_f1(
            integrated, batch_key=batch_key, label_key=label_key, embed=None
        )


def _compute_nmi(
    unintegrated: AnnData,
    label_key: str,
    cluster_key: str = "cluster",
    integrated: AnnData = None,
) -> float:
    """Compute biological conservation metric NMI

    NMI = normalized mutual information.
    Quantifies the overlap of the biological labels
    and Louvain clusters.
    Can be applied to all integration output types.

    Args:
        unintegrated:   AnnData object with unintegrated data.
                        This is used for metric computation when
                        no integrated AnnData is given
        label_key:      Column name of the biological labels
        cluster_key:    Key that stores the cluster labels
        integrated:     AnnData object with integrated data.
                        When None, the unintegrated AnnData is
                        used for metric computation

    Returns:
        float:        NMI score in [0, 1]
    """
    if not isinstance(integrated, AnnData):
        return scib.me.nmi(unintegrated, cluster_key=cluster_key, label_key=label_key)
    else:
        return scib.me.nmi(integrated, cluster_key=cluster_key, label_key=label_key)


def _compute_silhouette(
    unintegrated: AnnData,
    label_key: str,
    integrated: AnnData = None,
) -> float:
    """Compute biological conservation metric Silhouette

    Average Silhouette width.
    Only for feature and embedding integration outputs.

    Args:
        unintegrated:   AnnData object with unintegrated data.
                        This is used for metric computation when
                        no integrated AnnData is given
        label_key:      Column name of the biological labels
        integrated:     AnnData object with integrated data.
                        When None, the unintegrated AnnData is
                        used for metric computation

    Returns:
        float:        Silhouette score in [0, 1]
    """
    if not isinstance(integrated, AnnData):
        return scib.me.silhouette(unintegrated, label_key=label_key, embed="X_pca")
    else:
        return scib.me.silhouette(integrated, label_key=label_key, embed="X_emb")


def _compute_graph_connectivity(
    unintegrated: AnnData,
    label_key: str,
    integrated: AnnData = None,
) -> float:
    """Compute batch correction metric Graph Connectivity

    Assess whether the kNN graph directly connects all items
    with the same biological label.
    Can be used with all integration outputs.

    Args:
        unintegrated:   AnnData object with unintegrated data.
                        This is used for metric computation when
                        no integrated AnnData is given
        label_key:      Column name of the biological labels
        integrated:     AnnData object with integrated data.
                        When None, the unintegrated AnnData is
                        used for metric computation

    Returns:
        float:        graph connectivity score in [0, 1]
    """
    if not isinstance(integrated, AnnData):
        return scib.me.graph_connectivity(unintegrated, label_key=label_key)
    else:
        return scib.me.graph_connectivity(integrated, label_key=label_key)


def _compute_ilisi(
    unintegrated: AnnData,
    batch_key: str,
    integrated: AnnData = None,
) -> float:
    """Compute batch correction metric iLISI

    Quantify number of items that can be drawn from a neighbor list
    before one batch is observed twice.
    Can be used with all integration outputs.

    Args:
        unintegrated:   AnnData object with unintegrated data.
                        This is used for metric computation when
                        no integrated AnnData is given
        batch_key:      Column name of the batch labels
        integrated:     AnnData object with integrated data.
                        When None, the unintegrated AnnData is
                        used for metric computation

    Returns:
        float:        iLISI score in [0, 1]
    """
    if not isinstance(integrated, AnnData):
        return scib.me.ilisi_graph(unintegrated, batch_key=batch_key, type_="full")
    else:
        return scib.me.ilisi_graph(
            integrated, batch_key=batch_key, type_="embed", use_rep="X_emb"
        )


def _compute_kbet(
    unintegrated: AnnData,
    label_key: str,
    batch_key: str,
    integrated: AnnData = None,
) -> float:
    """Compute batch correction metric kBET

    Determine whether the label composition of a k nearest
    neighborhood of an item is similar to the expected (global)
    label composition.
    Can be used with all integration outputs.

    Args:
        unintegrated:   AnnData object with unintegrated data.
                        This is used for metric computation when
                        no integrated AnnData is given
        label_key:      Column name of the biological labels
        batch_key:      Column name of the batch labels
        integrated:     AnnData object with integrated data.
                        When None, the unintegrated AnnData is
                        used for metric computation

    Returns:
        float:        kBET score in [0, 1]
    """
    if not isinstance(integrated, AnnData):
        return scib.me.kBET(
            unintegrated,
            batch_key=batch_key,
            label_key=label_key,
            type_="full",
            embed="X_pca",
        )
    else:
        return scib.me.kBET(
            integrated,
            batch_key=batch_key,
            label_key=label_key,
            type_="embed",
            embed="X_emb",
        )


def _compute_pcr(
    unintegrated: AnnData,
    covariate: str,
    integrated: AnnData = None,
) -> float:
    """Compute batch correction metric PCR

    Principal component regression score.
    Compares variance contribution of variable batch before
    and after integration, so: Cannot be computed before integration!
    Can be used with feature & embedding integration outputs.

    Args:
        unintegrated:   AnnData object with unintegrated data.
                        This is used for metric computation when
                        no integrated AnnData is given
        covariate:      Column name to regress against (batch key)
        integrated:     AnnData object with integrated data.
                        When None, the unintegrated AnnData is
                        used for metric computation

    Returns:
        float:        PCR score in [0, 1]
    """
    if not isinstance(integrated, AnnData):
        return pd.NA
    else:
        return scib.me.pcr_comparison(
            unintegrated, integrated, covariate=covariate, embed="X_emb"
        )


def _compute_batch_asw(
    unintegrated: AnnData,
    label_key: str,
    batch_key: str,
    integrated: AnnData = None,
) -> float:
    """Compute batch correction metric batch ASW

    Batch ASW = batch average silhouette width
    Average ASW over each biological label. Silhouette measures how
    similar an object is to its own cluster compared to other clusters
    Can be used with feature & embedding integration outputs.

    Args:
        unintegrated:   AnnData object with unintegrated data.
                        This is used for metric computation when
                        no integrated AnnData is given
        label_key:      Column name of the biological labels
        batch_key:      Column name of the batch labels
        integrated:     AnnData object with integrated data.
                        When None, the unintegrated AnnData is
                        used for metric computation

    Returns:
        float:        batch ASW score in [0, 1]
    """
    if not isinstance(integrated, AnnData):
        return scib.me.silhouette_batch(
            unintegrated, batch_key=batch_key, label_key=label_key, embed="X_pca"
        )
    else:
        return scib.me.silhouette_batch(
            integrated, batch_key=batch_key, label_key=label_key, embed="X_emb"
        )
