# program processes paired-end qseq files for subsequent Bowtie alignment:
# !!! exactly the same as edit_fastqs_v3.4.pl except code is added to handle
#	CASAVA 1.8+ read IDs
#
# version 4.1 changed read ID naming to avoid readID mismatch errors in bwa v0.7.12
#	27Feb2015
#
use strict;
die "usage: edit_fastqs_v4.1.pl <read_fastq.txt> <mate_fastq.txt> <basename for output files>\n"
	unless (@ARGV == 3);
my $outfile = @ARGV[2];
open READS, "@ARGV[0]"
	or die "Can't open @ARGV[0] ($!)\n";
open MATES, "@ARGV[1]"
	or die "Can't open @ARGV[1] ($!)\n";
open OUT1, ">$outfile\_reads.txt";
open OUT2, ">$outfile\_mates.txt";
open OUT3, ">$outfile\_bad_reads.txt";
open OUT4, ">$outfile\_bad_mates.txt";
my ($count, $junk_reads, $junk_mates, $good_reads) = (0,0,0,0);
my $mate_clip;
my @read;
my @mate;
my $primary;
while (my $seqID = <READS>) {
	@read = ();
	@mate = ();
	$count++;
	unless ($seqID =~ /^@/) {
		die "Something wrong in read_fastq file! No \"@\" (read $count).\n";
	}
	push @read, $seqID; #retain CR
	my $seq = <READS>;
	chomp $seq;
	my $index = index($seq, 'TGTTA');
	if ($index < 0 || $index > 26) { # TGTTA not in adaptor sequence or not found at all
		$junk_reads++;
		print OUT3 ">$count\n$seq\n";
		<READS>; #throw away next two lines in READS and all four MATES lines
		<READS>;
		foreach (1..4) {<MATES>;}
		next;
	}
	my $clip = substr($seq, $index + 3);
	push @read, $clip;
	<READS>; #skip next line
	my $qual = <READS>;
	chomp $qual;
	push @read, substr($qual, $index + 3);
	#
	#now process mates file
	#
	$seqID = <MATES>;
	unless ($seqID =~ /^@/) {
		die "Something wrong in mate_fastq file! No \"@\" (read $count).\n";
	}
	push @mate, $seqID;
	$seq = <MATES>;
	chomp $seq;
	if ($seq =~ /^[GATCN]{21,26}GTG([GATC]{2}A[GATC]{2}A[GATC]{3})TGG/) { # $1 = barcode
		$mate_clip = index($seq, $1) + 21;
	} else { # no barcode found
#		if ($seq =~ /TG([GATC]{2}A[GATC]{2}A[GATC]{3})TGGTCGTGGTAT/) {print "$seq\n";}
		$junk_mates++;
		print OUT4 ">$count\n$seq\n";
		<MATES>;
		<MATES>;
		next;
	}
	my $clip = substr($seq, $mate_clip, 50);
	if (length($clip) < 50) { # read too short; this should not happen
		$junk_mates++;
		print OUT4 ">$count\n$seq\n";
		<MATES>;
		<MATES>;
		next;
	}
	my $barcode = $1;
	push @mate, $clip;
	<MATES>; #skip next line
	my $qual = <MATES>;
	chomp $qual;
	push @mate, substr($qual, $mate_clip, 50);
	#append barcode to read ID
	if ($read[0] =~ s/#.*1/_$barcode/) { # CASAVA 1.7 = @HWUSI-EAS100R:6:73:941:1973#NNNNNN/1
		$mate[0] =~ s/#.*\d/_$barcode/;
	} else { # CASAVA 1.8 = @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
		chomp $read[0];
		$read[0] =~ /^(.*)\s/;
		$read[0] = "$1_$barcode\n";
		$mate[0] = "$1_$barcode\n";
	}
	#now write results
	print OUT1 "$read[0]$read[1]\n+\n$read[2]\n";
	print OUT2 "$mate[0]$mate[1]\n+\n$mate[2]\n";
	$good_reads++;
}
close READS;
close MATES;
close OUT1;
close OUT2;
close OUT3;
close OUT4;
print "Processed $count alignments\n";
printf "%d (%.2f%%) failed 'TGTTA' filter on reads\n", $junk_reads, $junk_reads/$count*100;
printf "%d (%.2f%%) failed mate filters\n", $junk_mates, $junk_mates/$count*100;
printf "\t%d (%.2f%%) passed primary filter\n", $good_reads, $good_reads/$count * 100;
print "$good_reads alignments written to $outfile\_reads.txt, $outfile\_mates.txt\n";
