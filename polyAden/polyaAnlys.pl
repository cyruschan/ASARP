#!/usr/bin/perl
use warnings;
use strict;
BEGIN { push @INC, "../" }
require "anlysUtilities.pl"; #sub's for compSamples
require "mirnaParser.pl";


# set autoflush for error and output
select(STDERR);
$| = 1;
select(STDOUT);
$| = 1;

if(@ARGV < 4){

  print <<EOT;

USAGE: perl $0 asarp.gene.prediction trx.annotation.file sam.folder out.folder

This pipeline generates the allele-specifc FASTA files in the output folder, which are used to
analyze the allele-specific miRNA targeting of the alternative polyadenylation (alternative 
termination), i.e. ASAT, cases

asarp.gene.prediction	The ASARP .gene.prediction result file
trx.annotation.file	The gene transcript annotation file (also used in ASARP config as xiaofile)
sam.folder		The folder containing all aligned reads in SAM format; the folder should
			contain files chr*.sam (all suffixes are assumed to be .sam)
out.folder		The output folder (without the trailing "/")
			fasta file will be the AT.chromosome.gene.fa where chromosome and gene are
			extracted from asarp.gene.prediction. Each fa contains two sequences, one 
			for the reference allele and the other for the alternative allele.
EOT
  exit;
}
my ($asarp, $xiaoF, $u87Fd, $output) = @ARGV;

#my $u87Fd = "/home/cyruschan/asarp_perl/U87.siRNA/U87.siRNA.merged"; # siRNA
# get all ASARP results
my ($asarpGeneRef) = getAsarpAll($asarp);
my $spType = "AT";
my %at = %{$asarpGeneRef->{$spType}};

my $transRef = readTranscriptFile($xiaoF);
my $altRef = getGeneAltTransEnds($transRef); #get alternative initiation/termination (AI/AT) events from transcripts

my ($faRef) = genAlleleLocation(\%at, 25, $altRef, $u87Fd);
print "Output fasta files to folder: $output\n";
my %fa = %$faRef;
for my $key (keys %fa){
  my ($chr, $gene) = split(';', $key);
  open(FP, ">", "$output/$spType.$chr.$gene.fa") or die "ERROR: cannot output: $output/$spType.$chr.$gene.fa\n";
  print FP $fa{$key};
  close(FP);
}

