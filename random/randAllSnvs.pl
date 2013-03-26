#!/usr/bin/perl
use warnings;
use strict;
require "randUtilities.pl";

if(@ARGV < 3){
  print <<EOT;
  
USAGE: perl $0 org_snvs rand_folder indices [seed]

NOTE: generation of random SNVs for the ASARP SNV results for certain indices

ARGUMENTS:
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


my ($orgSnvs, $outFolder, $freq, $seed) = @ARGV;

if(defined($seed)){
  print "\nIMPORTANT: setting seed for srand will NOT generate RANDOM results.\n";
  print "           this seed setting is used only for debugging.\n";
  srand($seed);
}

my ($from, $to) = getIndex($freq);

# get the SNV list of input snvs
my @snvs = ();
my @sums = (); # to store the read sums
my $cnt = 0;
open(SP, $orgSnvs) or die "ERROR: cannot open original SNV file $orgSnvs to read\n";
while(<SP>){
  $cnt++;
  chomp $_;
  my ($chr, $pos, $al, $id, $reads) = split(' ', $_); 
  my $info = join(" ", $chr, $pos, $al, $id);
  my ($r1, $r2, $wnt) = split(':', $reads);
  $r1 += $r2;
  push @snvs, $info;
  push @sums, $r1;
}
close(SP);
print "\ntotal $cnt original SNVs\n";

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
  close(RFP);
}

