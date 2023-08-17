import scanpy as sc
import scib
import anndata as ad
from anndata import AnnData
from sklearn.preprocessing import MinMaxScaler
import pertpy as pt


def reset_corrected_anndata(adata: AnnData, obsm_key: str = "X_emb") -> AnnData:
    """Prepare corrected AnnData object for next level of batch correction

    After batch correction, this method creates a new AnnData object from
    a corrected one where the corrected features (usually in X_emb) are
    saved in X. This is helpful if you want to iterate batch correction
    methods on different levels to remove nested batch effects.

    Args:
        adata:          AnnData object containing the corrected data
        obsm_key:       Key under which the corrected data is stored
                        in obsm of `adata`

    Returns:
        AnnData:        AnnData object with X_emb saved in X
    """
    new_adata = ad.AnnData(adata.obsm[obsm_key].copy())
    new_adata.obs = adata.obs.copy()
    return new_adata


def scanorama_integration(
    adata: AnnData,
    batch: str,
    hierarchical: str = None,
) -> AnnData:
    """Function to perfrom integration with Scanorama

    Computes corrected feature matrix and embedding representation
    of corrected data.
    When this method is called on a higher level (i.e. use output
    of previous run as input), one should call `reset_corrected_anndata`
    on the input before. Else, the original data from adata.X is used
    for the next level and not the corrected values.

    Args:
        adata:          AnnData object containing the data
        batch:          Name of the column in adata.obs which contains
                        the batch labels
        hierarchical:   Doesn't have an influence on the batch correction
                        when None. Else: correction of a "lower level"
                        of nested batch effects. The batch correction is
                        then performed separately for each group that is
                        defined by the column hierarchical in adata.obs

    Returns:
        AnnData:        Batch corrected AnnData object
    """
    data = adata.copy()
    data.layers["counts"] = data.X

    if isinstance(hierarchical, str):
        # do batch correction on lower level (iterate over high-level batches)
        results = []
        for group in data.obs[hierarchical].unique():
            adata_part = data[data.obs[hierarchical] == group]
            adata_part = scib.ig.scanorama(adata_part, batch=batch)
            results.append(adata_part)
        data = ad.concat(results, merge="same")
        return data
    else:
        # do batch correction over all observations at once
        data = scib.ig.scanorama(data, batch=batch)
        return data


def scgen_integration(
    adata: AnnData,
    batch: str,
    labels: str,
    max_epochs: int = 100,
    hierarchical: str = None,
) -> AnnData:
    """Function to perfrom integration with scGen

    Computes corrected feature matrix.
    When this method is called on a higher level (i.e. use output
    of previous run as input), one should call `reset_corrected_anndata`
    on the input before. Else, the original data from adata.X is used
    for the next level and not the corrected values.

    Args:
        adata:          AnnData object containing the data
        batch:          Name of the column in adata.obs which contains
                        the batch labels
        labels:         Name of the column in adata.obs which contains
                        the biological labels
        max_epochs:     Max number of epochs for training the model
        hierarchical:   Doesn't have an influence on the batch correction
                        when None. Else: correction of a "lower level"
                        of nested batch effects. The batch correction is
                        then performed separately for each group that is
                        defined by the column hierarchical in adata.obs

    Returns:
        AnnData:        Batch corrected AnnData object
    """
    # import scgen as sg

    data = adata.copy()
    data.layers["counts"] = data.X

    if isinstance(hierarchical, str):
        # do batch correction on lower level (iterate over high-level batches)
        results = []
        for group in data.obs[hierarchical].unique():
            adata_part = data[data.obs[hierarchical] == group].copy()
            pt.tl.SCGEN.setup_anndata(adata_part, batch_key=batch, labels_key=labels)
            model = pt.tl.SCGEN(adata_part)
            model.train(
                max_epochs=max_epochs,
                batch_size=32,
                early_stopping=True,
                early_stopping_patience=25,
            )
            adata_part = model.batch_removal()
            results.append(adata_part)
        data = ad.concat(results, merge="same")
        return data
    else:
        # do batch correction over all observations at once
        pt.tl.SCGEN.setup_anndata(data, batch_key=batch, labels_key=labels)
        model = pt.tl.SCGEN(data)
        model.train(
            max_epochs=max_epochs,
            batch_size=32,
            early_stopping=True,
            early_stopping_patience=25,
        )
        data = model.batch_removal()
        return data


def harmony_integration(
    adata: AnnData,
    batch: str,
    hierarchical: str = None,
) -> AnnData:
    """Function to perfrom integration with Harmony

    Computes corrected feature matrix and embedding representation.
    When this method is called on a higher level (i.e. use output
    of previous run as input), one should call `reset_corrected_anndata`
    on the input before. Else, the original data from adata.X is used
    for the next level and not the corrected values.

    Args:
        adata:          AnnData object containing the data
        batch:          Name of the column in adata.obs which contains
                        the batch labels
        hierarchical:   Doesn't have an influence on the batch correction
                        when None. Else: correction of a "lower level"
                        of nested batch effects. The batch correction is
                        then performed separately for each group that is
                        defined by the column hierarchical in adata.obs

    Returns:
        AnnData:        Batch corrected AnnData object
    """
    data = adata.copy()
    data.layers["counts"] = data.X

    if isinstance(hierarchical, str):
        # do batch correction on lower level (iterate over high-level batches)
        results = []
        for group in data.obs[hierarchical].unique():
            adata_part = data[data.obs[hierarchical] == group]
            adata_part = scib.ig.harmony(adata_part, batch=batch)
            results.append(adata_part)
        data = ad.concat(results, merge="same")
        return data
    else:
        # do batch correction over all observations at once
        data = scib.ig.harmony(data, batch=batch)
        return data


def scanvi_integration(
    adata: AnnData,
    batch: str,
    labels: str,
    max_epochs: int = 400,
    hierarchical: str = None,
) -> AnnData:
    """Function to perfrom integration with scanVI

    Computes corrected feature matrix and embedding representation.
    When this method is called on a higher level (i.e. use output
    of previous run as input), one should call `reset_corrected_anndata`
    on the input before. Else, the original data from adata.X is used
    for the next level and not the corrected values.

    Args:
        adata:          AnnData object containing the data
        batch:          Name of the column in adata.obs which contains
                        the batch labels
        labels:         Name of the column in adata.obs which contains
                        the biological labels
        max_epochs:     Max number of epochs for training the model
        hierarchical:   Doesn't have an influence on the batch correction
                        when None. Else: correction of a "lower level"
                        of nested batch effects. The batch correction is
                        then performed separately for each group that is
                        defined by the column hierarchical in adata.obs

    Returns:
        AnnData:        Batch corrected AnnData object
    """
    scanvi = adata.copy()
    # min-max scaling because else you get exploding gradients
    scaler = MinMaxScaler()
    scanvi.X = scaler.fit_transform(scanvi.X)
    scanvi.layers["counts"] = scanvi.X

    if isinstance(hierarchical, str):
        # do batch correction on lower level (iterate over high-level batches)
        scanvi_results = []
        for group in scanvi.obs[hierarchical].unique():
            adata_part = scanvi[scanvi.obs[hierarchical] == group]
            adata_part = scib.ig.scanvi(
                adata_part, batch=batch, labels=labels, max_epochs=max_epochs
            )
            scanvi_results.append(adata_part)
        scanvi = ad.concat(scanvi_results, merge="same")
        return scanvi
    else:
        # do batch correction over all observations at once
        scanvi = scib.ig.scanvi(
            scanvi, batch=batch, labels=labels, max_epochs=max_epochs
        )
        return scanvi


def scvi_integration(
    adata: AnnData,
    batch: str,
    max_epochs: int = 400,
    hierarchical: str = None,
) -> AnnData:
    """Function to perfrom integration with scVI

    Computes corrected feature matrix and embedding representation.
    When this method is called on a higher level (i.e. use output
    of previous run as input), one should call `reset_corrected_anndata`
    on the input before. Else, the original data from adata.X is used
    for the next level and not the corrected values.

    Args:
        adata:          AnnData object containing the data
        batch:          Name of the column in adata.obs which contains
                        the batch labels
        max_epochs:     Max number of epochs for training the model
        hierarchical:   Doesn't have an influence on the batch correction
                        when None. Else: correction of a "lower level"
                        of nested batch effects. The batch correction is
                        then performed separately for each group that is
                        defined by the column hierarchical in adata.obs

    Returns:
        AnnData:        Batch corrected AnnData object
    """
    scvi = adata.copy()
    # min-max scaling because else you get exploding gradients
    scaler = MinMaxScaler()
    scvi.X = scaler.fit_transform(scvi.X)
    scvi.layers["counts"] = scvi.X

    if isinstance(hierarchical, str):
        # do batch correction on lower level (iterate over high-level batches)
        scvi_results = []
        for group in scvi.obs[hierarchical].unique():
            adata_part = scvi[scvi.obs[hierarchical] == group]
            adata_part = scib.ig.scvi(adata_part, batch=batch, max_epochs=max_epochs)
            scvi_results.append(adata_part)
        scvi = ad.concat(scvi_results, merge="same")
        return scvi
    else:
        # do batch correction over all observations at once
        scvi = scib.ig.scvi(scvi, batch=batch, max_epochs=max_epochs)
        return scvi
