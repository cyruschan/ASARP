#!/usr/bin/perl
use warnings;
use strict;
require "randUtilities.pl";

if(@ARGV < 4){
  print <<EOT;
  
USAGE: perl $0 asarp_result org_snvs rand_folder indices [seed]

NOTE: generation of random SNVs for the ASARP SNV results for certain indices

ARGUMENTS:
asarp_result		ASARP prediction output, in particular
			.gene.prediction is expected
org_snvs		file of the original SNV list for the ASARP pipeline
rand_folder		the output folder for the randomized SNV files
indices			The number of randomized SNV files, or the index range 
			for the randomized batch
			e.g. "20" means 20 randomized SNV lists indexed by 0 to 19
			e.g. "20:29" means 10 randomized SNV lists indexed by 20 to 29
			i.e. rand_folder/rand.SNV.20, rand_folder/rand.SNV.21, ...
OPTIONAL:
seed			fixed random seed for debugging

EOT
  exit;
}

die "This script is obsolete and replaced by more reaonsable randAllSnvs.pl\n";

my ($input, $orgSnvs, $outFolder, $freq, $seed) = @ARGV;

if(defined($seed)){
  print "\nIMPORTANT: setting seed for srand will NOT generate RANDOM results.\n";
  print "           this seed setting is used only for debugging.\n";
  srand($seed);
}

my ($from, $to) = getIndex($freq);

# steps
# 1. get the SNVs from xx.gene.prediction for randomization
# X2. get the SNVs from xx.ASE.prediction for exclusion (no need to be considered again)
#    in this case, is it the SNVs contained by those ASE genes that are to be exclueded?
#    however, this may remove control SNVs that are used by other genes as well.
#    therefore, it is better to keep all even the ASE SNVs (at the cost of longer comp' time)
# 3. retain all other SNVs in the original SNV list to be control candidates

# procedure for step 1
my ($snvsRef, $sumsRef, $snvHsRef) = getAsarpResult($input);
my @snvs = @$snvsRef;
my @sums = @$sumsRef;
my %snvHs = %$snvHsRef;

# the result will be output to a backup file as well: $outFolder/asarp.snv.bk
open(BFP, ">", "$outFolder/asarp.snv.bk") or die "ERROR: cannot open backup file to write for $input: $outFolder/asarp.snv.bk\n";
for(my $i=0; $i<@snvs; $i++){
  print BFP "$snvs[$i]\t$sums[$i]\n";
}
close(BFP);

# procedure for step 3
my $orgOutput = ();
my ($cnt, $kept, $skip) = (0, 0, 0);
open(SP, $orgSnvs) or die "ERROR: cannot open original SNV file $orgSnvs to read\n";
while(<SP>){
  $cnt++;
  chomp $_;
  my ($chr, $pos, $al, $id) = split(' ', $_); 
  my $check = join(" ", $chr, $pos, $al, $id);
  if(defined($snvHs{$check})){
    # this is in the ASARP result
    print "Skip $_\n";
    $skip++;
  }else{
    $orgOutput .= $_."\n";
    $kept++;
  }
}
close(SP);

print "\ntotal original $kept SNvs retained: skip $skip from $cnt original SNVs\n";

# the result will be output to a backup file as well: $outFolder/org_other.snv.bk
open(BSP, ">", "$outFolder/org_other.snv.bk") or die"ERROR: cannot open backup file to write for $input: $outFolder/org_other.snv.bk\n";
print BSP "$orgOutput";
close(BSP);

# generate random reads and put them to the outFolder
for(my $i=$from; $i<=$to; $i++){
  my $toOutput = "";
  for(my $j=0; $j<@snvs; $j++){
    my ($read1, $read2) = (0, 0);
    for(my $k=0; $k<$sums[$j]; $k++){
      # generate samples from a binomial distribution
      if(rand() < 0.5){
        $read1++;
      }else{
        $read2++;
      }
    }
    $toOutput .= "$snvs[$j] $read1:$read2:0\n"; #random reads
  }
  #print "\nRandomized SNV:\n$toOutput\n";
  open(RFP, ">", "$outFolder/rand.snv.$i") or print "WARNING: fail to open random SNV file for index $i\n";
  print RFP $toOutput;
  print RFP $orgOutput;
  close(RFP);
}

