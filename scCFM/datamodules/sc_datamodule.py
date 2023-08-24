import os
import numpy as np
from scCFM.scCFM.datamodules.components.sc_dataset import load_dataset
from functools import partial
from typing import Optional, List
import torch
from torch.utils.data import DataLoader, random_split, Dataset
from pytorch_lightning import LightningDataModule
from lightning.pytorch.utilities.combined_loader import CombinedLoader
import warnings
warnings.filterwarnings("ignore")


class CellDataset(Dataset):
    def __init__(self, data, batch, label, meta):
        self.data = data
        self.batch = batch
        self.label = label
        self.meta = meta

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        X = torch.tensor(self.data[idx])
        batch = torch.tensor(self.batch[idx])
        return {"X": X.to(torch.float32), "batch_key": batch.to(torch.float32), "label_key": self.label[idx],"meta": self.meta[idx]}

    
    
    
class scDataModule(LightningDataModule):
    pass_to_model = True

    def __init__(
        self,
        path: str,
        x_layer: str,
        batch_key: str,
        label_key: str,
        meta_key: os.strerror, 
        path_le: str,
        path_le_classes: str, 
        use_pca: bool, 
        n_dimensions: Optional[int] = None, 
        train_val_test_split: List = [1.0, 0.0, 0.0],
        batch_size: Optional[int] = 64,
        num_workers: int = 0
    ):
        super().__init__()
                
        assert os.path.exists(path), "The data path does not exist"
        
        # Collect dataset 
        self.data, self.batch, self.label, self.meta = load_dataset(path=path, 
                                            x_layer=x_layer,
                                            batch_key = batch_key,
                                            label_key=label_key,
                                            meta_key = meta_key, 
                                            path_le = path_le,
                                            path_le_classes = path_le_classes,
                                            use_pca=use_pca, 
                                            n_dimensions=n_dimensions)
        
        self.in_dim = self.data.shape[1] 
        self.n_cond = self.batch.shape[1] # number of batches/conditions 
        self.n_labels = max(self.label)+1 # number of classes for clf
        self.train_val_test_split = train_val_test_split
        self.batch_size = batch_size
        self.num_workers = num_workers
        
        # Associate time to the process 
        #self.idx2cond = {idx: val for idx, val in enumerate(np.unique(self.cond))}
        
        # Create dataset 
        self.dataset = CellDataset(self.data, self.batch, self.label, self.meta)
        
        # Create a list of subset data
        self.split()

    def split(self):
        self.split_data = random_split(self.dataset, self.train_val_test_split)

    def train_dataloader(self):
        return DataLoader(
            dataset=self.split_data[0],
            batch_size=self.batch_size,
            num_workers=self.num_workers,
            drop_last=True,
            shuffle=True
        )

    def val_dataloader(self):
        return DataLoader(
            dataset=self.split_data[1],
            batch_size=self.batch_size,
            num_workers=self.num_workers,
            drop_last=True,
            shuffle=False
        )
    
    def test_dataloader(self):
        return DataLoader(
            dataset=self.split_data[2],
            batch_size=self.batch_size,
            num_workers=self.num_workers,
            drop_last=True,
            shuffle=False
        )
        