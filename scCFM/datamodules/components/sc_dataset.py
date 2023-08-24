from typing import Optional
import numpy as np
import scanpy as sc
from sklearn import preprocessing
import torch
import torch.nn as nn
import joblib

def adata_dataset(path: str, 
                  x_layer: str = None,
                  batch_key: str = "source",
                  label_key: str = "compound",
                  meta_key: str = "welltype", # another metadata column that should be extracted for further analysis
                  path_le: str = "label_encoder.joblib",
                  path_le_classes: str = "label_encoder_classes.joblib", 
                  use_pca: bool = False,
                  n_dimensions: Optional[int] = None
                  ):
    """Loads AnnData from a path

    Args:
        path (str): the path to the AnnData object
        x_layer (str, optional): the layer in AnnData to use 
        time_key (str, optional): the key in adata.obs representing the path. Defaults to "experimental_time".
        path_le (str, optional): path to save the label_encoder to 
        use_pca (bool, optional): whether to use pca or gene counts. Defaults to False.
        n_dimensions (Optional[int], optional): number of dimensions to keep in PCA. Defaults to None.
    """
    
    # Read AnnData 
    adata = sc.read_h5ad(path)

        
    # 1. Define X 
    if use_pca:
        X = adata.obsm["X_pca"][:, n_dimensions]
    else:
        if x_layer in adata.layers:
            X = adata.layers[x_layer]
        elif x_layer in adata.obsm:
            X = adata.obsm[x_layer]
        else:
            X = adata.X
            
    # 2. Define batch 
    #convert batch keys to numbers 
    le = preprocessing.LabelEncoder()
    targets = le.fit_transform(adata.obs[batch_key])
    joblib.dump(le, path_le) # save label encoder for later use 
    # one-hot encode the classes
    onehot = nn.functional.one_hot(torch.as_tensor(targets))
    
    #3. Extract further metadata
    meta = list(adata.obs[meta_key])

    #4. Extract biological label as integer
    le_label = preprocessing.LabelEncoder()
    classes = le_label.fit_transform(adata.obs[label_key])
    joblib.dump(le_label, path_le_classes) 
    


    return X, onehot, torch.tensor(classes), meta

def load_dataset(path: str, 
                 x_layer: str = None, 
                 batch_key: str = "source",
                 label_key: str = "compound",
                 meta_key: str = "welltype", 
                 path_le: str = "label_encoder.joblib",
                 path_le_classes: str = "label_encoder_classes.joblib", 
                 use_pca: bool = False, 
                 n_dimensions: int = None
                 ):
    """Wrapper around adata_dataset function to implement controls
    """
    if path.endswith("h5ad"):
        return adata_dataset(path=path,
                             x_layer=x_layer,
                             batch_key = batch_key,
                             label_key=label_key,
                             meta_key = meta_key, 
                             path_le = path_le, 
                             path_le_classes = path_le_classes,
                             use_pca=use_pca, 
                             n_dimensions=n_dimensions)
    else:
        raise NotImplementedError()
    