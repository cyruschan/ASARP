#!/usr/bin/perl -w
use strict;

require "fileParser.pl"; #sub's for input annotation files
require "snpParser.pl"; #sub's for snps
use MyConstants qw( $CHRNUM $supportedList $supportedTags );

# set autoflush for error and output
select(STDERR);
$| = 1;
select(STDOUT);
$| = 1;

if(@ARGV < 2){
  print <<EOT;
  
USAGE: perl $0 output_file config_file [optional: parameter_file]

Get the ASE and Powerful SNVs from the input SNV file (in config_file)
The criteria (powerful and ASE, i.e. p-value, cutoffs) are contained in the parameter_file
There will be two result files: 
"output_file.ase" with the ASE SNVs (subset of powerful SNVs), and 
"output_file.pwr" with the powerful SNVs

EOT
  exit;
}

# input arguments: $outputFile--output, $configs--configuration file for input files, $params--configuration file for parameters
my ($outputFile, $configs, $params) = getArgs(@ARGV); 
my ($snpF, $bedF, $rnaseqF, $xiaoF, $splicingF, $estF) = getRefFileConfig($configs); # input annotation/event files
my ($POWCUTOFF, $SNVPCUTOFF, $FDRCUTOFF, $ASARPPCUTOFF, $NEVCUTOFFLOWER, $NEVCUTOFFUPPER, $ALRATIOCUTOFF) = getParameters($params); # parameters

my ($snpRef, $pRef) = initSnp($snpF, $POWCUTOFF);
#print "SNV List:\n";
#printListByKey($snpRef, 'powSnps');
# suggested, get the Chi-Squared Test p-value cutoff from FDR ($FDRCUTOFF)
if(defined($FDRCUTOFF)){
  print "Calculating the Chi-Squared Test p-value cutoff for FDR <= $FDRCUTOFF...\n";
  if(defined($SNVPCUTOFF)){
    print "NOTE: user-provided p-value in config: $SNVPCUTOFF is ignored.\n";
  }
  $SNVPCUTOFF = fdrControl($pRef, $FDRCUTOFF);
  print "Chi-Squared Test p-value cutoff: $SNVPCUTOFF\n";
}

#basic statistics for powerful SNVs but not genes
my $powSnvCnt = 0;
my $aseSnvCnt = 0;

open(AP, ">", "$outputFile.ase") or die "ERROR: cannot open $outputFile.ase for ASE SNVs\n";
open(PP, ">", "$outputFile.pwr") or die "ERROR: cannot open $outputFile.pwr for Powerful SNVs\n";

for(my $i=1; $i<=$CHRNUM; $i++){
  my ($snpChrRef) = getListByKeyChr($snpRef, 'powSnps', $i);
  my %powSnps = %$snpChrRef;
  my $powCntChr = keys %powSnps;  
  my $aseCntChr = 0;
  
  #basic statistics
  $powSnvCnt += $powCntChr; # keys %powSnps;
  for(keys %powSnps){
    my @allSnpInfo = split(';', $powSnps{$_}); #separate by ;, if there are multiple snps at the same position
    for(@allSnpInfo){
      my ($p, $pos, $alleles, $snpId, $refAl, $altAl) = getSnpInfo($_);
      #get allelic ratios
      my $r = $refAl/($refAl+$altAl);
      my $toOutput = join(" ", formatChr($i), $pos, $alleles, $snpId, $r, "$refAl:$altAl:0", $p);
      if($p <= $SNVPCUTOFF){
        $aseCntChr += 1; #each SNV **location** added once
	print AP $toOutput."\n"; #ase
      }
      print PP $toOutput."\n"; #powerful, which is a superset of ase
      last; # I assume there can be multiple SNVs at the same position, but usually ppl don't
      #therefore, a "last" is shot to get only the first one in the position
    }
  }
  if($powCntChr >0){  printChr($i); print" has $powCntChr powerful SNVs and $aseCntChr SNVs\n";  }
  $aseSnvCnt += $aseCntChr;
}

close(AP);
close(PP);
print "\nThere are in total $aseSnvCnt ASE SNVs out of $powSnvCnt powerful SNVs\n";
print "Done\n";

