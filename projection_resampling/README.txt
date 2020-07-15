Projection Resampling script
----------------------------

This python script can be used to project insertion counts from TnSeq
across multiple samples onto principle components determined by PCA or
varimax dimensions.  All input and output files are in tab-separated
format, which can be opened as spreadsheets in Excel.

The main input (counts.txt) is a spreadsheet (tab-separated format) of
the *normalized* counts at each TA site (rows) across individual datasets
(columns).  It is critical that the counts are properly normalized
across datasets (to account for differences in total number of reads).
It is the user's responsibility to do this before running this script.

The metadata is a spreadsheet indicating which datasets are associated
with which condition (i.e. as replicates).

The loadings.txt file incorporates the loadings of the conditions onto
PCs or varimax dimensions (determined independently, e.g. using R).

The "prot_table" is customized format for the genes in the genome (annotation).
  see https://transit.readthedocs.io/en/stable/transit_running.html#prot-tables-annotations

The last file (gene_list.txt) is for selecting genes; the first column
is the ORF id, and the last column is an adjusted p-value.  The
results will be printed out only for genes with Padj<0.05.  One way of
using this is to perform ANOVA on all the genes in the genome, to
identify the subset of genes exhibiting significant variability, and
put the adjusted p-value in the final column.  Alternatively, you can
just provide a target list of genes of interest (with a 0 in the final
column).

Usage and Example
-----------------

> python projection_resampling_Abx.py <counts.txt> <samples_metadata.txt> <genome.prot_table> <loadings.txt> <gene_list.txt> 

Using the example files provided, the command would be:

> python projection_resampling_Abx.py multiTnSeq_Abx_3b_counts_TTR.txt Abx_3b_samples_metadata_all.txt H37Rv3.prot_table multiTnSeq_Abx_3b_varimax_loadings.txt michelles_list.txt > projection_resampling_out.txt

The output file can be opened as a spreadsheet.  It shows the original
mean counts in each of the conditions, the mean counts projected onto
the PCs (or varimax dimensions), the log-fold-changes (LFCs) of the
projected means for each dimension (computed relative to the pooled
mean count across the other PCs), and the adjusted p-values based on
resampling of each PC against the all the others), indicating the
significance of association of each gene with each PC (or varimax
dimension).
