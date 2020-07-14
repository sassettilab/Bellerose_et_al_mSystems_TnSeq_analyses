#!/bin/bash

# TnSeq_pipe6.0.sh is TnSeq_pipe5.0.3.sh modified to use barcodes4bwa_v2.0
# No changes otherwise.
#
# requires edit_fastqs_v4.1.pl to account for CASAVA 1.8+ read IDs
# REQUIRES FULL PATH NAME TO INPUT FILES (v5.0.2)
# BWA indexes in local subdirectioy --> ./bwa_indexes_tb
#
# As for TnSeq_pipe5.0.3.sh
#	requires barcode_use_v1.0.pl, get_raw_counts_piped.pl
#	in working directory along with edit_fastqs_v4.1.pl and barcodes4bwa_v2.0

# !!! local version calls only 2 cpu's (recommend 8 if working on cluster).

# Required in PATH:
#	perl v5.18.1 (or later) and module List::MoreUtils
#	samtools v0.0.19
#	bwa v0.7.12

if [ $# -eq 0 ]
then
	echo "usage: TnSeq_pipe6.0.sh <reads1.fastq|reads1.fastq.gz> <reads2.fastq|reads2.fastq.gz><basename for output>"
	exit 1 # exit shell script
fi

FILETAG=$RANDOM
filename="$1"
fileext=${filename##*.}

mkdir -p output_v6.0_$3
cd output_v6.0_$3
cp ../bwa_indexes_tb/* ./

echo '##############################################################'
echo 'TnSeq_pipe6.0.'
echo `date -R`
echo "input files: $1, $2"
echo "name for output: $3"
echo '##############################################################'

if [ ${fileext} = "gz" ]
then
    echo "### expanding $1..."
	zcat $1 > $FILETAG\_1.fastq
	echo "### expanding $2..."
	zcat $2 > $FILETAG\_2.fastq
else
	echo '### input files are not gzipped'
	ln -s $1 $FILETAG\_1.fastq
	ln -s $2 $FILETAG\_2.fastq
fi

echo '### filtering valid reads and stripping adapter sequences...'
perl ../edit_fastqs_v4.1.pl $FILETAG\_1.fastq $FILETAG\_2.fastq $3
echo '### aligning...'
bwa aln -q 5 -l 32 -k 2 -o 1 -t 2 NC_018143.1.fa $3\_reads.txt > $3\_reads.sai
bwa aln -q 5 -l 32 -k 2 -o 1 -t 2 NC_018143.1.fa $3\_mates.txt > $3\_mates.sai
bwa sampe -P NC_018143.1.fa $3\_reads.sai $3\_mates.sai $3\_reads.txt $3\_mates.txt > $3\_sampe.sam

echo '### counting hits...'
perl ../barcodes4bwa_v2.0.pl $3\_sampe.sam $3
echo '### analyzing barcode use...'
perl ../barcode_use_v1.0.pl $3\_details.txt
echo '### writing files of raw counts...'
perl ../get_raw_counts_piped.pl $3\_summary.txt 0 > $3\_raw_breaks.txt
perl ../get_raw_counts_piped.pl $3\_summary.txt 1 > $3\_raw_barcodes.txt
perl ../get_raw_counts_piped.pl $3\_summary.txt 2 > $3\_raw_templates.txt
echo '### compressing bad reads files...'
gzip $3\_bad_reads.txt
gzip $3\_bad_mates.txt
echo "### converting ${3}_sampe.sam to BAM..."
samtools view -Sb $3\_sampe.sam > $3\_sampe.bam
echo  "### run completed `date -R`"

# clean up
rm $FILETAG*
rm NC_018143.1.*
rm $3\_reads.txt
rm $3\_mates.txt
rm *.sai
rm $3\_sampe.sam

