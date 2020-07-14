Resampling Statistics for Evaluating Conditional Fitness (Treatment vs. Control)

Perl program 'AllStats_tb_v9.3.3mm3.pl' performs resampling statistics comparing Tn-seq libraries grown under different conditions, i.e., treatment vs. control. Input files are tab-delimited text containing insertion counts (column 2) at each genomic TA coordinate (column 1). Program was designed for replicate libraries for each condition, although a single library for either is permitted. Program reads the input from a tab-delimited text configuration file in which count file names are listed in column 1 with either a zero (0) or one (1) in column 2, zero to indicate control and one indicating treatment.

The program also requires an annotation file corresponding to the genome build from which the TA coordinates of the count files were obtained. The file provided here corresponds to M. tuberculosis H37Rv build NC_018143.1 (Broad Institute), which can be used to reproduce the results reported in Bellerose et al. (2020) Distinct bacterial pathways influence the efficacy of antibiotics against Mycobacterium tuberculosis. mSystems vol:pp.

Prerequisites (in PATH):
	Perl (e.g. perl v5.18.1) with the following modules installed:
		List::Util
		List::MoreUtils
		Statistics::R
	R (e.g. R v3.5.1) with libraries
		DESeq
		NOISeq

usage: perl AllStats_tb_v9.3.3mm3.pl <notations.table2> <basename for output> <config_file.txt>

Output directory contains the following files:

	basename_AllStats_v9.3.3_out.txt - table of normalized insertion counts and resampling statistics summarized by ORF	basename_counts_raw.txt - table of raw insertion counts by ORF	basename_debug.txt - for debugging use only	basename_ecdf_plots.pdf - empirical cumulative distribution (ECDF) plots of insertion counts by library	basename_ecdf.R - R script to generate ECDF plots	basename_missed_TAs.txt - list of TA sites with no insertions in any library	basename_normalized.igv - IGV browser tracks of normalized insertion counts	basename_normalized.txt - table of normalized insertion counts (by TA site)	basename_orf_counts.txt - table of normalized insertion counts (by ORF)
	basename_qq_plots - directory containing quantile-quantile plots of normalized counts for each library vs. a standardized reference library 	basename_R - directory of R scripts that can be executed to produce histograms of resampling results for each ORF	basename_ratio_histogram.pdf - density plot of log2(ratios), treatment/control
	basename_hist.R - R script to generate basename_ratio_histogram.pdf	basename_raw.igv - IGV browser tracks of raw counts (not normalized)	basename_volcano.pdf - volcano plot of resampling results	basename.R - R script to produce basename_volcano.pdf

Sample files are included here and can be processed using command line:

perl AllStats_tb_v9.3.3mm3.pl NC_018143_noTerm2.ptt.table2 day0 pretreat_config.txt

NOTE: The resampling algorithm implemented here has been reproduced in a user friendly suite of programs at https://orca1.tamu.edu/essentiality/transit described in DeJesus, M.A., Ambadipudi, C., Baker, R., Sassetti, C., and Ioerger, T.R. (2015). TRANSIT - a Software Tool for Himar1 TnSeq Analysis. PLOS Computational Biology, 11(10):e1004401.
