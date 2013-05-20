#!/usr/bin/perl
use warnings;
use strict;
require "randUtilities.pl";

if(@ARGV < 4){
  print <<EOT;
  
USAGE: perl $0 ase_result rand_results indices output_p [p_cutoff]

NOTE: evaluation of p-values for the individual ASARP SNV results based on randomization for certain indices

ARGUMENTS:
ase_result		ASARP AllASE gene prediction output, in particular
			.ase.prediction is expected
rand_results		the output folder and file prefix for the randomized AllASE gene result files
			ASARP pipeline consistent suffix is assumed for result files, i.e. .ase.prediction
			rand_results will be used as the common prefix for these files:
			e.g. "result/rand" for all "result.rand.*.ase.prediction" files
indices			The number of randomized SNV files, or the index range 
			for the randomized batch
			e.g. "20" means 20 randomized SNV lists indexed by 0 to 19
			e.g. "20:29" means 10 randomized SNV lists indexed by 20 to 29
			i.e. rand_results.20.gene.prediction, 
			     rand_results.21.gene.prediction, ...
			indice range >= N is needed to get enough power
			for p-value resolution ~1/N (e.g. N=1000 for p~0.001)
output_p		output ASARP SNV results with p-values associated

OPTIONAL:
p_cutoff		p-value cutoff (<=) for the output results output_p, 
			if not input, all p values are output
EOT
  exit;
}


my ($input, $randF, $freq, $outputF, $pCutoff) = @ARGV;

if(!defined($pCutoff)){ $pCutoff = 1; }
my ($from, $to) = getIndex($freq);
my $N = $to-$from+1; #the total number for simulations

# get the real ASARP AllASE results
my ($aseHsRef) = getAseResult($input);

# get the randomized ASARP AllAse results
# the difference is that, rather than using
# $randAseHsRef as the list (which may contain ASE genes beyond $aseHsRef),
# $aseHsRef is always used as the actual AllASE gene list

my %randAse = ();
for(my $i=$from; $i<=$to; $i++){
  my $randFileName = "$randF.$i.ase.prediction";
  my ($randAseHsRef) = getAseResult($randFileName);
  my %oneRandHs = %$randAseHsRef;
  for(keys %oneRandHs){
    if($oneRandHs{$_} != 1){
      die "ERROR: each gene is expected to appear at most once in AllASE: $_ appears $oneRandHs{$_} times\n";
    }
    $randAse{$_} += 1;
  }
}

my $toOutputP = "";
my %realAse = %$aseHsRef;
my ($sumP, $pCnt) = (0, 0);
$pCnt = keys %realAse; # total number of genes
print "There are in total $pCnt AllASE genes to be evaluated\n";

for(keys %realAse){
  my $p = 0;
  if(defined($randAse{$_})){
    $p = $randAse{$_}/$N; # averaged
    $sumP += $p;
  }
  if($p <= $pCutoff){
    my ($chr, $gene) = split(';', $_);
    $toOutputP .= "$chr $gene $p\n";
  }
}
print $toOutputP;
print "OVerall ASE FDR estimated: ";
print "%.4f\n", $sumP/$pCnt;
print "Output P-Values of ASE Genes to $outputF\n";
open(OP, ">", $outputF) or die "ERROR: cannot open output: $outputF\n";
print OP $toOutputP;
close(OP);
print "Finished\n";
