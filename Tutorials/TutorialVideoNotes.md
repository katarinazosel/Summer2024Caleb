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

To find differentially expressed genes filter using p < 0.05 and log2Fold change > 0.5







