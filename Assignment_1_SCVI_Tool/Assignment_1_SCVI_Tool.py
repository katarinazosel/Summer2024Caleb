#!/usr/bin/env python
# coding: utf-8

# ## Project Configuration

# In[1]:


from scvi_colab import install
install(run_outside_colab=True)


# In[137]:


import os
import tempfile

import anndata
import muon
import numpy as np
import pooch
import scanpy as sc
import scvi 
import seaborn as sns
import torch
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt


# In[179]:


sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")

#Change to Ouputs when done
save_dir = 'C:/Users/chpar/OneDrive/Documents/GitHub/Summer2024Caleb/Assignment_1_SCVI_Tool/Saved_Model'

#Inputs:
BATCHES_CSV_PATH = "Test_Paths.csv"
N_TOP_GENES = 1200
OUTPUT_PATH = "Outputs"


#Variable equivalencies between tutorial and actual
#Tutorial Name:          Actual Name:
#"cell_source"              "batch"


# In[4]:


scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)


# ## Loading Data from CSV Paths

# In[322]:


df = pd.read_csv(BATCHES_CSV_PATH)
batches = df.to_dict(orient="index")

# Get dataframe labels from CSV file
csv_labels = df.columns.to_list()

batches


# In[323]:


print(batches.items())
print (csv_labels)


# In[324]:


adatas = {}

for batch in batches.items():
    path = batch[1][csv_labels[0]]
    batch_label = batch[1][csv_labels[1]]

    batch_adata = sc.read_10x_mtx(path)
    batch_adata.var.add_suffix(suffix=batch_label) 
    adatas[batch_label] = batch_adata

adata = ad.concat(adatas, label="batch")
adata.obs_names_make_unique()

adata


# In[325]:


print(adata.obs.sample(n=5))


# In[326]:


adata


# ## Loading Data from Tutorial

# In[71]:


adata = scvi.data.heart_cell_atlas_subsampled(save_path=save_dir)
adata


# ## Preprocessing Data

# In[327]:


# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")


# In[328]:


sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)


# In[329]:


print("# cells, # genes before filtering:", adata.shape)

sc.pp.filter_genes(adata, min_counts=3)
sc.pp.filter_cells(adata, min_counts=3)

print("# cells, # genes after filtering:", adata.shape)


# In[330]:


#Save counts information before normalizing
adata.layers["counts"] = adata.X.copy()


# In[331]:


#Normalization of data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

#store normalized values in .raw
adata.raw = adata


# In[336]:


adata.obs_keys


# In[332]:


# Identify highly variable ggenes
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=N_TOP_GENES,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key="batch",
)


# In[333]:


sc.pl.highly_variable_genes(adata)


# ## Running SCVI Model

# In[337]:


scvi.model.SCVI.setup_anndata(
    adata, 
    layer="counts", 
    categorical_covariate_keys=["batch"],
    continuous_covariate_keys=["pct_counts_mt", "pct_counts_ribo"],
    batch_key="batch")

model = scvi.model.SCVI(adata)


# In[338]:


model


# In[339]:


model.train()


# ## Save/Load Model

# In[340]:


model_dir = save_dir 
#model.save(save_dir, overwrite=True)   


# In[341]:


#Run to load model from save
model_dir = save_dir 
model = scvi.model.SCVI.load(model_dir, adata=adata)


# In[342]:


model


# ## Obtaining Outputs

# In[343]:


SCVI_LATENT_KEY = "X_scVI"

latent = model.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent
latent.shape


# In[344]:


SCVI_NORMALIZED_KEY = "scvi_normalized"
adata.layers[SCVI_NORMALIZED_KEY] = model.get_normalized_expression(library_size=10e4)


# #### Visualization without Batch Correction

# In[414]:


sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=20)
sc.tl.umap(adata, min_dist=0.3)

sc.tl.leiden(adata, resolution=0.3)


# In[370]:


np.random.seed(0)
random_indices = np.random.permutation(list(range(adata.shape[0])))

sc.pl.umap(
    adata[random_indices, :],
    color=["leiden", "batch"],
    ncols=2,
    frameon=False,
)


# In[371]:


sc.tl.paga(adata, groups="leiden") # change groups to the name of the clusters

sc.pl.paga(adata)


# In[372]:


adata.obs_keys


# #### Visualization with batch correction

# In[430]:


#Use scVI latent space for UMAP
sc.pp.neighbors(adata, n_pcs=10, n_neighbors=20) #, use_rep=SCVI_LATENT_KEY
sc.tl.umap(adata, min_dist=0.3)


# In[431]:


SCVI_CLUSTERS_KEY = "leiden_scVI"
sc.tl.leiden(adata, key_added=SCVI_CLUSTERS_KEY, resolution=0.05)


# In[432]:


sc.pl.umap(
    adata,
    color=[SCVI_CLUSTERS_KEY],
    frameon=False,
)
sc.pl.umap(
    adata,
    color=["batch"],
    ncols=2,
    frameon=False,
)


# In[448]:


sc.tl.paga(adata, groups="leiden") # change groups to the name of the clusters

sc.pl.paga(adata)


# In[449]:


sc.tl.draw_graph(adata, init_pos='paga', random_state=123, layout='fa')
sc.pl.draw_graph(adata, color=["leiden"], legend_loc='on data') # plotting


# In[354]:


adata.var_keys


# In[385]:


markers = {"Alpha":["GCG", "LOXL4", "GC", "CRYBA2","TTR","IRX2","TM4SF4","ARX","MAFB", "PGR"], #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5092539/
           "Beta":["INS", "IAPP","MAFA","NPTX2","DLK1","TGFBR3","SYT13","SMAD9","CDKN1C","TFCP2L1","SIX3","MNX1","MBP5","PIR", "G6PC2","ADCYAP1","PDX1"], #"ADCYAP1","PDX1",#ADCYAP1#https://www.nature.com/artiMAFBcles/s41587-022-01219-z#MOESM10 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5092539/
           "New Beta":[  "SPD", "INS", "MMP-2","CK-19"],#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2842599/
            "Delta":["SST", "RBP4", "HHEX", "PCSK1", "PRG4" "BHLHE41", "LEPR"],#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5092539/
            "Acinar":["GP2", "CTRC", "AMY2A", "CPA1", "SYCN"], #https://www.nature.com/articles/s41598-019-40481-1/tables/1
            "Cancer stem cells":[ "OCT4","MYC","KLF4","SOX2"],
            "Ductal":["KRT19","SOX9"], #https://www.nature.com/articles/s41598-019-40481-1/tables/1
            "Endothelial":["CD93","ACE","CD143","C1qR1","CD31","CD34"], #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7749106/
            "Macrophages":["CD206","CD68","CD163"],
            "PP": ["PPY","ETV1","ARX","PAX6","SERTM1""CARTPT","SLITRK6"], #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5092539/
            "Pancreatic progenitor": ["CPA1","NEUROG3","NKX6.1","PTF1A","SOX9","C-MYC","SOX4","NEUROD1"], #https://www.nature.com/articles/s41467-017-00561-0
            "Additional factors": ["SOX2","OCT4","NANOG","PDX1"] 
            }


# In[457]:


for cell_type, genes in markers.items():
        print(f"marker gene of {cell_type}")
        for gene in genes:
            if gene in adata.var_names:
                sc.pl.umap(adata, color=gene, use_raw=True, ncols=2,save=f'{gene}_expression_umap.png')
                plt.show()
os.rename("figures", "gene_expression")


# In[450]:


adata.write("./8human_fetal_sander_lynn_scVI_batch_corrected.h5ad")


# In[451]:


adata = sc.read_h5ad("./8human_fetal_sander_lynn_scVI_batch_corrected.h5ad")


# In[405]:




