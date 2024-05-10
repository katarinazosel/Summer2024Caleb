## Differential Expression with pyDESeq2  

link: https://www.youtube.com/watch?v=wIvxFEMQVwg

Prepping Counts table
 - Set index to ensemble ID(first col)
 - Remove rows with only 0s
 - Transpose table

Create metadata dataframe if it isn't already in a separate file

Create deseqdataset object (adata object) using counts, metadata, and conditions/batches

Run deseq2 function on object, adds rows and cols for statistical values

Run DeseqStats and then .summary() to do the actual statistical calculations

Get the results as a data frame using .results_df 

Convert ensemble IDs to gene symbols by mapping a dictionary to the index column in a new column of the data frame (sanbomics has a tool to help)

Filter genes with low expression using baseMean parameter(threshold is somewhat arbitrary but video recommends bM < 10 should be removed)

To find differentially expressed genes filter using padj < 0.05 and log2Fold change > 0.5

#### Analysis

###### PCA

Uses scanpy, pass deseqdataset object as argument scanpy pca function, then plot(Dot size must be >= 200 because there are few dots)

###### GSEA

Link: https://www.youtube.com/watch?v=yOQcrUMCALw

Uses gseapy, Uses filtered dataframe of differentially expressed genes either output from deseq or imported file

Add rank column which is -log10(padj) * log2FoldChange
 - padj, and log2FoldChange are columns in the dataframe

Sort dataframe by rank

Create new dataframe using only gene column and rank column

If a gene list is unavailable, built in lists are available. 
To make a custom list you can just make a dictionary of {term : list of genes}

Call the gsea function. gp.prerank(rank, geneset, rng_seed)

.results are a dictionary. Loop over to pull out values and place in a dataframe
    - Pulled out fdr, es, and nes values

plot a term using the gseaplot() function

###### Heatmap

Uses numpy and seaborn
Uses the deseqdataset object

Add a layer of the log1P values of the normed_counts layer

To make heatmap need lop1P values, var names(cols) and obs names(rows). These can be put in a dataframe for graphing(dataframe must be transposed so that var_names are the rows, and obs_names are the columns)

Plot using seaborn clustermap

###### Volcano Plot

Using Sanbomics package and the result of the deseq2 process

run the volcano() function passing the results and the symbol that we want to highlight

## Single-cell pseudotime and gene regulatory analysis with CellOracle

Allows the simulation of gene perturbations using scRNA-seq data

Link:https://www.youtube.com/watch?v=in5_keFSArU
Documentation: https://morris-lab.github.io/CellOracle.documentation/

Uses scanpy, numpy, pandas, matplotlib, and celloracle packages

Do preprocessing/filtering

Normalize data and calculate log1P

Calculate highly variable genes

Do PCA, neighbors, then cluster

Do a umap, annotate cells

Initialize pseudotime object, pass adata, dimensional reduction data name, clustering data name

Create a lineage dictionary of the cells

Set and plot the lineage

identify and find the root cells

set the root cell in the pseudotime object

calculate the diffusiton map for the adat

calculate pseudotime lineage

pseudotime values are stored in the dataframe

#### In Silico Perturbation Analysis
###### Set up and calculate GRN data

Need a base_GRN data, follow base_GRN input data tutorial in documentation if using scATAC-seq or Bulk ATAC-seq data
Built in base_GRN datasets are available too

initialize oracle object, import adata with raw counts, and base_GRN

Run a PCA and select the pcs to use
- Autoselect function is in oracle documentation
- k should be 0.025 * number of cells (according to documentation)

Compute knn imputations using the values given above

Calculate GRN for each cell group

Filter out insignificant links (p= 0.001) and keep top 10k(default number)

Compute network scores and add them to dataframe

GRN data is used for perturbation analysis

###### Perturbation Analysis

Add GRN data to oracle data and fit GRN for simulation

Simulate shift by setting target gene expression between 0 and 2 times the max expression

Using previously calculated pseudotime to calculate developmental flow

add pseudotime observation column

create a gradient object using the oracle object and pseudotime

To compare, create an oracle devlopment object, pass the gradient and perturbation simulation data
Then calucalte the inner product

Inner product
- 1 if initial and perturbed vectors point in same direction
- -1 if initial and perturbed vectors point in opposite directions


































