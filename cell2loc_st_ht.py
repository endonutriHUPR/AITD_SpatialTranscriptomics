import numpy as np #matrices
import scanpy as sc #models
import matplotlib.pyplot as plt #images
import matplotlib as mpl # images
import cell2location 
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs
import seaborn as sns

### Single cell

adata_ref = sc.read_h5ad(
    f'together_filt_rawCounts_singlecellHT.h5ad'
)
adata_ref.obs["Manual_annotation"]=adata_ref.obs["Manual_annotation"].astype('category')



from cell2location.utils.filtering import filter_genes

selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=2)

# filter the object
adata_ref = adata_ref[:, selected].copy()


# create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)

# view anndata_setup as a sanity check
mod.view_anndata_setup()


mod.train(max_epochs=250,
    #batch_size=2500,
    train_size=1,
    lr=0.002)
    
    

mod.plot_history(20)



adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
)

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
adata_file


mod.plot_QC()


### ST samples

## HT1
HT1=sc.read_visium(path="/HT1/outs/",
                    count_file="filtered_feature_bc_matrix.h5") 
HT1.var_names_make_unique()
HT1.obs['sample'] = list(HT1.uns['spatial'].keys())[0]
HT1.var['SYMBOL'] = HT1.var_names

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(HT1.var_names, inf_aver.index)
HT1 = HT1[:, intersect].copy()
inf_aver_HT1 = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=HT1)

# create and train the model
mod_HT1_ncells30_epoch20000_8000genes = cell2location.models.Cell2location(
    HT1, cell_state_df=inf_aver_HT1
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=30,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)
mod_HT1_ncells30_epoch20000_8000genes.view_anndata_setup()

mod_HT1_ncells30_epoch20000_8000genes.train(max_epochs=20000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1
         )

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod_HT1_ncells30_epoch20000_8000genes.plot_history(1000)
plt.legend(labels=['full data training']);

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
HT1 = mod_HT1_ncells30_epoch20000_8000genes.export_posterior(
    HT1, sample_kwargs={'num_samples': 1000, 'batch_size': mod_HT1_ncells30_epoch20000_8000genes.adata.n_obs}
)

# Save model


mod_HT1_ncells30_epoch20000_8000genes.save(f"{run_name}", overwrite=True)

#mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# Save anndata object with results
adata_file = f"{run_name}/sp_HT2.h5ad"
HT1.write(adata_file)
adata_file

## HT2

HT2=sc.read_visium(path="/HT2/outs/",
                    count_file="filtered_feature_bc_matrix.h5") 
HT2.var_names_make_unique()
HT2.obs['sample'] = list(HT2.uns['spatial'].keys())[0]
HT2.var['SYMBOL'] = HT2.var_names

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(HT2.var_names, inf_aver.index)
HT2 = HT2[:, intersect].copy()
inf_aver_HT2 = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=HT2)

# create and train the model
mod_HT2_ncells30_epoch20000_8000genes = cell2location.models.Cell2location(
    HT2, cell_state_df=inf_aver_HT2
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=30,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)
mod_HT2_ncells30_epoch20000_8000genes.view_anndata_setup()

mod_HT2_ncells30_epoch20000_8000genes.train(max_epochs=20000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1
         )

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod_HT2_ncells30_epoch20000_8000genes.plot_history(1000)
plt.legend(labels=['full data training']);

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
HT2 = mod_HT2_ncells30_epoch20000_8000genes.export_posterior(
    HT2, sample_kwargs={'num_samples': 1000, 'batch_size': mod_HT2_ncells30_epoch20000_8000genes.adata.n_obs}
)

# Save model


mod_HT2_ncells30_epoch20000_8000genes.save(f"{run_name}", overwrite=True)

#mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# Save anndata object with results
adata_file = f"{run_name}/sp_HT2.h5ad"
HT2.write(adata_file)
adata_file

## HT3

HT3=sc.read_visium(path="/HT3/outs/",
                    count_file="filtered_feature_bc_matrix.h5") 
HT3.var_names_make_unique()
HT3.obs['sample'] = list(HT2.uns['spatial'].keys())[0]
HT3.var['SYMBOL'] = HT2.var_names

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(HT3.var_names, inf_aver.index)
HT3 = HT3[:, intersect].copy()
inf_aver_HT3 = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=HT3)

# create and train the model
mod_HT3_ncells30_epoch20000_8000genes = cell2location.models.Cell2location(
    HT3, cell_state_df=inf_aver_HT3
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=30,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)
mod_HT3_ncells30_epoch20000_8000genes.view_anndata_setup()

mod_HT3_ncells30_epoch20000_8000genes.train(max_epochs=20000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1
         )

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod_HT3_ncells30_epoch20000_8000genes.plot_history(1000)
plt.legend(labels=['full data training']);

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
HT3 = mod_HT3_ncells30_epoch20000_8000genes.export_posterior(
    HT3, sample_kwargs={'num_samples': 1000, 'batch_size': mod_HT3_ncells30_epoch20000_8000genes.adata.n_obs}
)

# Save model


mod_HT3_ncells30_epoch20000_8000genes.save(f"{run_name}", overwrite=True)

#mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# Save anndata object with results
adata_file = f"{run_name}/sp_HT2.h5ad"
HT3.write(adata_file)
adata_file
