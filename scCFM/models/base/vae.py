import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.optim import Adam
from torch.distributions import Normal, kl_divergence
import pytorch_lightning as pl

from scCFM.scCFM.models.base.mlp import MLP
from scCFM.scCFM.models.base.mlp import Classifier


class BaseAutoencoder(pl.LightningModule):
    def __init__(
        self,
        in_dim,
        hidden_dims,
        batch_norm,
        dropout,
        dropout_p,
        per_cell_var, #boolean
        conditional_ae, #boolean
        n_cond, # number of conditions
        activation=torch.nn.ReLU,
        likelihood="gaussian"
    ):
        super(BaseAutoencoder, self).__init__()

        # Attributes
        self.in_dim = in_dim 
        self.hidden_dims = hidden_dims
        self.batch_norm = batch_norm
        self.activation = activation
        self.dropout = dropout
        self.dropout_p = dropout_p
        self.likelihood = likelihood
        self.latent_dim = hidden_dims[-1]
        self.per_cell_var = per_cell_var
        self.conditional_ae = conditional_ae
        self.n_cond = n_cond 
        
        if conditional_ae: 
            # Encoder
            self.encoder_layers = MLP(
                hidden_dims=[in_dim + n_cond, *hidden_dims[:-1]], # changed in_dim 
                batch_norm=batch_norm,
                dropout=dropout,
                dropout_p=dropout_p,
                activation=activation,
            )
            # Decoder
            h = hidden_dims[::-1]
            h[0] = h[0] + n_cond
            self.decoder_layers = MLP(
                hidden_dims = h,  # make sure first layer has dim latent_dim+n_cond
                batch_norm=batch_norm,
                dropout=dropout,
                dropout_p=dropout_p,
                activation=activation,
            )
        else: 
            # Encoder
            self.encoder_layers = MLP(
                hidden_dims=[in_dim, *hidden_dims[:-1]],
                batch_norm=batch_norm,
                dropout=dropout,
                dropout_p=dropout_p,
                activation=activation,
            )
            # Decoder
            self.decoder_layers = MLP(
                hidden_dims=[*hidden_dims[::-1]], 
                batch_norm=batch_norm,
                dropout=dropout,
                dropout_p=dropout_p,
                activation=activation,
            )
           
        
        if likelihood == "gaussian":
            if per_cell_var:
                self.log_sigma = torch.nn.Linear(hidden_dims[0], self.in_dim)
                self.decoder_mu = torch.nn.Linear(hidden_dims[0], self.in_dim)

            else:
                self.log_sigma = torch.nn.Parameter(torch.randn(self.in_dim))
                self.decoder_mu = torch.nn.Linear(hidden_dims[0], self.in_dim)

        else:
            raise NotImplementedError

    def encode(self, x):
        pass

    def decode(self, z):
        h = self.decoder_layers(z)
        
        if self.likelihood == "gaussian":
            if self.per_cell_var:   
                log_sigma = self.log_sigma(h)
                mu = self.decoder_mu(h)
                return dict(mu=mu,
                            log_sigma=log_sigma)
            else: 
                mu = self.decoder_mu(h)
                return dict(mu=mu)
        else:
            raise NotImplementedError

            
    def forward(self, batch):
        pass

    def reconstruction_loss(self, x, decoder_output):
        
        if self.likelihood == "gaussian":
            if self.per_cell_var:
                mu = decoder_output["mu"]
                sigma = decoder_output["log_sigma"]
                distr = Normal(loc=mu, scale=torch.exp(sigma)) #scale – standard deviation
                recon_loss = -distr.log_prob(x).sum(-1)
            else:
                mu = decoder_output["mu"]
                log_sigma_reshaped = self.log_sigma.repeat(mu.shape[0], 1)
                distr = Normal(loc=mu, scale=torch.exp(log_sigma_reshaped))
                recon_loss = -distr.log_prob(x).sum(-1)

        else:
            raise NotImplementedError
        return recon_loss
    
    def configure_optimizers(self):
        return Adam(self.parameters(), lr=1e-4)


class scANVI_CVAE(BaseAutoencoder):
  
    def __init__(
        self,
        in_dim,
        hidden_dims,
        hidden_dims_clf, # added 
        batch_norm,
        dropout,
        dropout_p,
        per_cell_var, 
        conditional_ae, 
        n_cond, 
        n_labels, # added, number of eg compounds, classes for the classifier
        n_epochs: int,
        activation=torch.nn.ReLU,
        likelihood:str="gaussian",
        kl_warmup_fraction: float=0.5,
        kl_weight: float=None,
        ce_weight: float=1 #added
    ):
        super(scANVI_CVAE, self).__init__(
            in_dim,
            hidden_dims,
            batch_norm,
            dropout,
            dropout_p,
            per_cell_var, 
            conditional_ae,
            n_cond,
            activation,
            likelihood
        )
       
        # Latent space
        self.mu = nn.Linear(hidden_dims[-2], self.latent_dim)
        self.logvar = nn.Linear(hidden_dims[-2], self.latent_dim)
        # classifier
        self.classifier = Classifier(
                hidden_dims=[self.latent_dim, *hidden_dims_clf, n_labels],  
                batch_norm=batch_norm,
                dropout=dropout,
                dropout_p=dropout_p,
                activation=activation)
        # cross entropy
        self.cross_entropy = nn.CrossEntropyLoss() 

        if kl_weight==None:
            self.kl_weight = 0
            self.kl_weight_increase = 1/(kl_warmup_fraction*n_epochs)
            self.anneal_kl = True
        else:
            self.kl_weight = kl_weight
            self.anneal_kl = False

        self.ce_weight = ce_weight
        
    def encode(self, x):
        h = self.encoder_layers(x)
        mu, logvar = self.mu(h), self.logvar(h)
        z = self.reparameterize(mu, logvar) 
        return dict(z=z,
                    mu=mu,
                    logvar=logvar)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        z = mu + eps * std
        return z

    def forward(self, batch):
        x = batch["X"]
        cond = batch["batch_key"]
        if x.shape[0] > 1: #concatenate
          x_cond = torch.cat((x, cond), dim = 1) 
        else:
          x_cond = torch.cat((x, cond), dim = 0) 
        z, mu, logvar = self.encode(x_cond).values()
        # add classifier
        clf_output = self.classifier(z)

        if z.shape[0] > 1: #concatenate
          z_cond = torch.cat((z, cond), dim = 1)
        else: 
          z_cond = torch.cat((z, cond), dim = 0)

        return self.decode(z_cond), mu, logvar, clf_output  #also return output of classifier 

        
    def step(self, batch, prefix):
        x = batch["X"]
        decoder_output, mu, logvar, clf_output = self.forward(batch) # also output of classifier returned 
        recon_loss = self.reconstruction_loss(x, decoder_output)
        kl_div = self.kl_divergence(mu, logvar)
        kl_weight = min([self.kl_weight, 1])
        # add cross entropy 
        ce_loss = self.cross_entropy(clf_output, batch["label_key"])

        loss = torch.mean(recon_loss + kl_weight * kl_div + self.ce_weight*ce_loss) # add ce
        dict_losses = {f"{prefix}/loss": loss,
                       f"{prefix}/kl": kl_div.mean(),
                       f"{prefix}/lik": recon_loss.mean(),
                       f"{prefix}/ce": ce_loss}
        self.log_dict(dict_losses, prog_bar=True)
        if prefix == "train":
            return loss
    
    def kl_divergence(self, mu, logvar):
        p =  Normal(mu, torch.sqrt(torch.exp(0.5 * logvar))) 
        q = Normal(torch.zeros_like(mu), torch.ones_like(mu))
        kl_div = kl_divergence(p, q).sum(dim=-1)
        return kl_div
    
    def training_step(self, batch, batch_idx):
        loss = self.step(batch, "train")
        return loss

    def validation_step(self, batch, batch_idx):
        self.step(batch, "val")
    
    def on_train_epoch_end(self):
        if self.anneal_kl:
            self.kl_weight += self.kl_weight_increase
        else:
            self.kl_weight += 0




class CVAE(BaseAutoencoder):
  
    def __init__(
        self,
        in_dim,
        hidden_dims,
        batch_norm,
        dropout,
        dropout_p,
        per_cell_var, #boolean
        conditional_ae, #boolean
        n_cond, # number of conditions
        n_epochs: int,
        activation=torch.nn.ReLU,
        likelihood:str="gaussian",
        kl_warmup_fraction: float=0.5,
        kl_weight: float=None
    ):
        super(CVAE, self).__init__(
            in_dim,
            hidden_dims,
            batch_norm,
            dropout,
            dropout_p,
            per_cell_var, 
            conditional_ae,
            n_cond,
            activation,
            likelihood
        )
       
        # Latent space
        self.mu = nn.Linear(hidden_dims[-2], self.latent_dim)
        self.logvar = nn.Linear(hidden_dims[-2], self.latent_dim)
        if kl_weight==None:
            self.kl_weight = 0
            self.kl_weight_increase = 1/(kl_warmup_fraction*n_epochs)
            self.anneal_kl = True
        else:
            self.kl_weight = kl_weight
            self.anneal_kl = False
        
    def encode(self, x):
        h = self.encoder_layers(x)
        mu, logvar = self.mu(h), self.logvar(h)
        z = self.reparameterize(mu, logvar) 
        return dict(z=z,
                    mu=mu,
                    logvar=logvar)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        z = mu + eps * std
        return z

    def forward(self, batch):
        x = batch["X"]
        cond = batch["batch_key"]
        if x.shape[0] > 1: #concatenate
          x_cond = torch.cat((x, cond), dim = 1) 
        else:
          x_cond = torch.cat((x, cond), dim = 0) 
        z, mu, logvar = self.encode(x_cond).values()
        if z.shape[0] > 1: #concatenate
          z_cond = torch.cat((z, cond), dim = 1)
        else: 
          z_cond = torch.cat((z, cond), dim = 0)
        return self.decode(z_cond), mu, logvar

        
    def step(self, batch, prefix):
        x = batch["X"]
        decoder_output, mu, logvar = self.forward(batch)
        recon_loss = self.reconstruction_loss(x, decoder_output)
        kl_div = self.kl_divergence(mu, logvar)
        kl_weight = min([self.kl_weight, 1])

        loss = torch.mean(recon_loss + kl_weight * kl_div)
        dict_losses = {f"{prefix}/loss": loss,
                       f"{prefix}/kl": kl_div.mean(),
                       f"{prefix}/lik": recon_loss.mean()}
        self.log_dict(dict_losses, prog_bar=True)
        if prefix == "train":
            return loss
    
    def kl_divergence(self, mu, logvar):
        p =  Normal(mu, torch.sqrt(torch.exp(0.5 * logvar))) 
        q = Normal(torch.zeros_like(mu), torch.ones_like(mu))
        kl_div = kl_divergence(p, q).sum(dim=-1)
        return kl_div
    
    def training_step(self, batch, batch_idx):
        loss = self.step(batch, "train")
        return loss

    def validation_step(self, batch, batch_idx):
        self.step(batch, "val")
    
    def on_train_epoch_end(self):
        if self.anneal_kl:
            self.kl_weight += self.kl_weight_increase
        else:
            self.kl_weight += 0
          
  
    
    
    
    
