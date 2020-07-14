Data Processing Pipeline for Counting Tn-seq Insertions

Prerequisites:

	Executables in PATH:
		samtools v0.0.19
		bwa v0.7.12
		perl (e.g. perl v5.18.1)

	Helper programs in working directory:
		barcode_use_v1.0.pl
		get_raw_counts_piped.pl
		edit_fastqs_v5.0.pl
		barcodes4bwa_v2.0.pl

	Perl Module: List::MoreUtils

	bwa indexes in subdirectory 'bwa_indexes_tb'

Processing executable is TnSeq_pipe7.0.local.sh.

Call: TnSeq_pipe7.0.local.sh <reads1.fastq|reads1.fastq.gz> <reads2.fastq|reads2.fastq.gz> <name for output directory>"

Full pathnames to the reads files are required! If gzipped, read file .fq.gz is allowed. Requires read lengths ≥75 bases!

To test with included sample files:

./TnSeq_pipe7.0.local.sh ../../sample_files/L8_CGATCATG_L008_R1_001.fastq.gz ../../sample_files/L8_CGATCATG_L008_R2_001.fastq.gz testout








