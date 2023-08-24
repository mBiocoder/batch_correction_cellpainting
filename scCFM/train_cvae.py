import sys
import pytorch_lightning as pl
import torch

sys.path.insert(0, "../")
#from paths import EXPERIMENT_FOLDER

from scCFM.scCFM.datamodules.sc_datamodule import scDataModule
from scCFM.scCFM.models.base.vae import CVAE
from scCFM.scCFM.models.base.vae import scANVI_CVAE

from pytorch_lightning import Trainer
from pytorch_lightning.callbacks import ModelCheckpoint, EarlyStopping
from pytorch_lightning.loggers import WandbLogger


# Training function
class Solver:
    def __init__(self,config):

        self.task_name = config["training"]["task_name"]
        # Fix seed for reproducibility
        self.seed = config["training"]["seed"]
        torch.manual_seed(self.seed)
        if self.seed:
            pl.seed_everything(self.seed, workers=True)

        # Initialize datamodule
        self.init_datamodule(config["datamodule"])
         # Initialize model
        self.init_model(config["model"])
         # Initialize callbacks
        self.init_checkpoint_callback(config["model_checkpoint"])
        # Initialize callbacks
        self.init_early_stopping_callback(config["early_stopping"])
        # Initialize logger
        self.init_logger(config["logger"])
         # Initialize the lightning trainer
        self.init_trainer(config["trainer"])


    def init_datamodule(self,config_datamodule):
        self.datamodule = scDataModule(path=config_datamodule["path"],
                                       x_layer=config_datamodule["x_layer"],
                                       batch_key=config_datamodule["batch_key"],
                                       label_key=config_datamodule["label_key"],
                                       meta_key=config_datamodule["meta_key"],
                                       path_le = config_datamodule["path_le"],
                                       path_le_classes = config_datamodule["path_le_classes"],
                                       use_pca=config_datamodule["use_pca"],
                                       n_dimensions=config_datamodule["n_dimensions"],
                                       train_val_test_split=config_datamodule["train_val_test_split"],
                                       batch_size=config_datamodule["batch_size"],
                                       num_workers=config_datamodule["num_workers"])

    def init_model(self, config_model):
        self.model = CVAE(in_dim=self.datamodule.in_dim,
                         hidden_dims=config_model["hidden_dims"],
                         batch_norm=config_model["batch_norm"],
                         dropout=config_model["dropout"],
                         dropout_p=config_model["dropout_p"],
                         n_epochs=config_model["n_epochs"],
                         likelihood=config_model["likelihood"],
                         per_cell_var=config_model["per_cell_var"],
                         conditional_ae=config_model["conditional_ae"], 
                         n_cond=self.datamodule.n_cond,
                         kl_warmup_fraction=config_model["kl_warmup_fraction"],
                         kl_weight=config_model["kl_weight"])

    def init_checkpoint_callback(self,config_checkpoint):
        self.model_ckpt_callbacks = ModelCheckpoint(filename=config_checkpoint["filename"],
                                                    monitor=config_checkpoint["monitor"],
                                                    mode=config_checkpoint["mode"],
                                                    save_last=config_checkpoint["save_last"],
                                                    auto_insert_metric_name=config_checkpoint["auto_insert_metric_name"])

    def init_early_stopping_callback(self,config_es):
        self.early_stopping_callbacks = EarlyStopping(monitor=config_es["monitor"],
                                                      patience=config_es["patience"],
                                                      mode=config_es["mode"],
                                                      min_delta=config_es["min_delta"],
                                                      verbose=config_es["verbose"],
                                                      strict=config_es["strict"],
                                                      check_finite=config_es["check_finite"],
                                                      stopping_threshold=config_es["stopping_threshold"],
                                                      divergence_threshold=config_es["divergence_threshold"],
                                                      check_on_train_epoch_end=config_es["check_on_train_epoch_end"]
                                                      )

    def init_logger(self,config_logger):
        self.logger = WandbLogger(offline=config_logger["offline"],
                                  name=config_logger["name"],
                                  id=config_logger["id"],
                                  project=config_logger["project"],
                                  log_model=config_logger["log_model"],
                                  prefix=config_logger["prefix"],
                                  group=config_logger["group"],
                                  tags=config_logger["tags"],
                                  job_type=config_logger["job_type"])

    def init_trainer(self,config_trainer):
        self.trainer = Trainer(callbacks=[self.model_ckpt_callbacks],#, self.early_stopping_callbacks],
                               logger=self.logger,
                               max_epochs=config_trainer["max_epochs"],
                               accelerator=config_trainer["accelerator"],
                               devices=config_trainer["devices"],
                               log_every_n_steps=config_trainer["log_every_n_steps"],
                               limit_val_batches=0,
                               num_sanity_val_steps=0)

    def train(self):
        # Fit the model
        self.trainer.fit(model=self.model,
                         train_dataloaders=self.datamodule.train_dataloader())#,
                         #val_dataloaders=self.datamodule.val_dataloader())

        train_metrics = self.trainer.callback_metrics
        return train_metrics




class Solver_scANVI_CVAE:
    def __init__(self,config):

        self.task_name = config["training"]["task_name"]
        # Fix seed for reproducibility
        self.seed = config["training"]["seed"]
        torch.manual_seed(self.seed)
        if self.seed:
            pl.seed_everything(self.seed, workers=True)

        # Initialize datamodule
        self.init_datamodule(config["datamodule"])
         # Initialize model
        self.init_model(config["model"])
         # Initialize callbacks
        self.init_checkpoint_callback(config["model_checkpoint"])
        # Initialize callbacks
        self.init_early_stopping_callback(config["early_stopping"])
        # Initialize logger
        self.init_logger(config["logger"])
         # Initialize the lightning trainer
        self.init_trainer(config["trainer"])


    def init_datamodule(self,config_datamodule):
        self.datamodule = scDataModule(path=config_datamodule["path"],
                                       x_layer=config_datamodule["x_layer"],
                                       batch_key=config_datamodule["batch_key"],
                                       label_key=config_datamodule["label_key"],
                                       meta_key=config_datamodule["meta_key"],
                                       path_le = config_datamodule["path_le"],
                                       path_le_classes = config_datamodule["path_le_classes"],
                                       use_pca=config_datamodule["use_pca"],
                                       n_dimensions=config_datamodule["n_dimensions"],
                                       train_val_test_split=config_datamodule["train_val_test_split"],
                                       batch_size=config_datamodule["batch_size"],
                                       num_workers=config_datamodule["num_workers"])

    def init_model(self, config_model):
        self.model = scANVI_CVAE(in_dim=self.datamodule.in_dim,
                         hidden_dims=config_model["hidden_dims"],
                         hidden_dims_clf=config_model["hidden_dims_clf"],
                         batch_norm=config_model["batch_norm"],
                         dropout=config_model["dropout"],
                         dropout_p=config_model["dropout_p"],
                         n_epochs=config_model["n_epochs"],
                         likelihood=config_model["likelihood"],
                         per_cell_var=config_model["per_cell_var"],
                         conditional_ae=config_model["conditional_ae"], 
                         n_cond=self.datamodule.n_cond,
                         n_labels=self.datamodule.n_labels,
                         kl_warmup_fraction=config_model["kl_warmup_fraction"],
                         kl_weight=config_model["kl_weight"],
                         ce_weight=config_model["ce_weight"])

    def init_checkpoint_callback(self,config_checkpoint):
        self.model_ckpt_callbacks = ModelCheckpoint(filename=config_checkpoint["filename"],
                                                    monitor=config_checkpoint["monitor"],
                                                    mode=config_checkpoint["mode"],
                                                    save_last=config_checkpoint["save_last"],
                                                    auto_insert_metric_name=config_checkpoint["auto_insert_metric_name"])

    def init_early_stopping_callback(self,config_es):
        self.early_stopping_callbacks = EarlyStopping(monitor=config_es["monitor"],
                                                      patience=config_es["patience"],
                                                      mode=config_es["mode"],
                                                      min_delta=config_es["min_delta"],
                                                      verbose=config_es["verbose"],
                                                      strict=config_es["strict"],
                                                      check_finite=config_es["check_finite"],
                                                      stopping_threshold=config_es["stopping_threshold"],
                                                      divergence_threshold=config_es["divergence_threshold"],
                                                      check_on_train_epoch_end=config_es["check_on_train_epoch_end"]
                                                      )

    def init_logger(self,config_logger):
        self.logger = WandbLogger(offline=config_logger["offline"],
                                  name=config_logger["name"],
                                  id=config_logger["id"],
                                  project=config_logger["project"],
                                  log_model=config_logger["log_model"],
                                  prefix=config_logger["prefix"],
                                  group=config_logger["group"],
                                  tags=config_logger["tags"],
                                  job_type=config_logger["job_type"])

    def init_trainer(self,config_trainer):
        self.trainer = Trainer(callbacks=[self.model_ckpt_callbacks],#, self.early_stopping_callbacks],
                               logger=self.logger,
                               max_epochs=config_trainer["max_epochs"],
                               accelerator=config_trainer["accelerator"],
                               devices=config_trainer["devices"],
                               log_every_n_steps=config_trainer["log_every_n_steps"],
                               limit_val_batches=0,
                               num_sanity_val_steps=0)

    def train(self):
        # Fit the model
        self.trainer.fit(model=self.model,
                         train_dataloaders=self.datamodule.train_dataloader())#,
                         #val_dataloaders=self.datamodule.val_dataloader())

        train_metrics = self.trainer.callback_metrics
        return train_metrics



