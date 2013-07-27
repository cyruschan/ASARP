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

#USAGE: perl $0 asarp.gene.prediction trx.annotation.file miRNA_prediction out.folder
  print <<EOT;

USAGE: perl $0 asarp.gene.prediction trx.annotation.file out.folder

This pipeline analyzes the allele-specific alternative polyadenylation
(alternative termination), i.e. ASAT, cases.
EOT
  exit;
}
#my ($asarp, $xiaoF, $miFile, $output) = @ARGV;
my ($asarp, $xiaoF, $output) = @ARGV;
# get all ASARP results
my ($asarpGeneRef) = getAsarpAll($asarp);
my $spType = "AT";
my %at = %{$asarpGeneRef->{$spType}};

my $transRef = readTranscriptFile($xiaoF);
my $altRef = getGeneAltTransEnds($transRef); #get alternative initiation/termination (AI/AT) events from transcripts

my ($faRef) = genAlleleLocation(\%at, 20, $altRef);
print "Output fasta files to folder: $output\n";
my %fa = %$faRef;
for my $key (keys %fa){
  my ($chr, $gene) = split(';', $key);
  open(FP, ">", "$output/ASAT.$chr.$gene.fa") or die "ERROR: cannot output: $output/ASAT.$chr.$gene.fa\n";
  print FP $fa{$key};
  close(FP);
}

###########################################################################
# sub-routines

sub checkMiRanda{
  ####### overlap with miRanda
  my ($atRef, $miFile) = @_;
  my %at = %$atRef;
  my ($miHitRef, $miInfoRef) = miRandaResParser(\%at, $miFile);

  my %hits = %$miHitRef;
  my %miInfo = %$miInfoRef;

  for my $key (keys %hits){
  
    #print "$key\n"; #chr;gene
    #get the genome location information
    my @mirs = split (/\t/, $hits{$key});
    for(@mirs){
      my ($miAc, $miSt, $miEnd, $miStrand) = split(':', $_);

  
      # just the raw part in gene.prediction:
      # e.g. AI:5+(0.1362),AS:^RI(0.2010)^ASS(0.2010);44683815 na T>C 47:89;+
      my @snvs = split(/\t/, $at{$key}); #dummyInfo;snpInfo;$strandInfo
      for(@snvs){
        my ($dummyInfo, $snpInfo, $strandInfo) = split(';', $_);
        # get 3' UTR region
        my ($pos, $id, $alleles, $count) = split(' ', $snpInfo);
        if($pos <= $miEnd && $pos >= $miSt){ # overlaps
          print "HIT: $miAc: $miSt-$miEnd, $miStrand\n$snpInfo;$dummyInfo\n$miInfo{$miAc}\n\n";
        }
      }
    }
  }
}
