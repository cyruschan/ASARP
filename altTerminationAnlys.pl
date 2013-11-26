#!/usr/bin/perl
use warnings;
use strict;
require "anlysUtilities.pl"; #sub's for compSamples

# set autoflush for error and output
select(STDERR);
$| = 1;
select(STDOUT);
$| = 1;

# modify specificAseAnlys.pl to customize 3'UTR analysis
if(@ARGV < 9){

  print <<EOT;

USAGE: perl $0 prefix1 config1 param1 asarp1 prefix2 config2 param2 asarp2 output [asarp_path]

This pipeline combines several ASARP scripts to analyze allele-specific alternative
polyadenylation / 3' termination (AT) of the 3'UTR ASE SNVs which can be miRNA
targets. The analysis can be performed between two samples (datasets), e.g. two 
poly(A)+ samples from different celluar compartments: cytocol (C.A+) and nucleus 
(N.A+).

The first input will be the focused case, e.g. set C.A+ to be the first input and 
N.A+ second to focus on C.A+ specific 3'UTR ASE SNVs (excluding N.A+)

In the input, 1, 2 indicate the corresponding in/output for samples (datasets) 
1 and 2, respectively. 

prefix		the intermediate output prefix to distinguish
		severval files with: .ase .ase_distri.list, etc.
		They are generated as the intermediate output for each sample
		NOTE: make sure the prefix names are **DIFFERENT** between 
		prefix1 and prefix2
config		the .config file which is required for ASARP pipeline entry
param		the .param file, optional in ASARP pipeline, is required in
		this pipeline
asarp		the ASARP output result (main summary file) path predicted.
		the corresponding .ase.prediction and .gene.prediction will
		be used (have to be present with the main summary's folder)
output		the final output file of this analysis pipeline

OPTIONAL:
asarp_path	the -I path for the ASARP pipeline folder if it is run outside
		the original folder. Because this pipeline executes several
		scripts in ASARP: aseSnvs.pl, snp_distri.pl, etc., the ASARP
		path is needed for them.

EOT

exit;
}

my ($intron1, $config1, $param1, $result1, $intron2, $config2, $param2, $result2, $output, $iPath) = @ARGV;
if(!defined($iPath)){
  $iPath = ".";
}
if(substr($iPath, -1, 1) ne "/" ){
  $iPath .= "/"; #add the path slash
}

if($intron1 eq $intron2){
  die "ERROR: prefixes 1 and 2 must be different: now $intron1 and $intron2\n";
}


my $specificType = "3'UTR";
# focus on the specific SNVs only and compare them with the 2 sample results
print "Analyzing $specificType ASE SNVs\n";


print "\n[[I]]. Obtain all $specificType ASE SNVs from sample 1 (case)\n";
my ($iAseGenesRef1, $iAseSnvsRef1, $iAsarpGenesRef1, $iAsarpSnvsRef1) = specificAsePipeline($intron1, $config1, $param1, $result1, $specificType, $iPath);

print "\n[[II]]. Obtain all $specificType ASE SNVs from sample 2 (background)\n";
my ($iAseGenesRef2, $iAseSnvsRef2, $iAsarpGenesRef2, $iAsarpSnvsRef2) = specificAsePipeline($intron2, $config2, $param2, $result2, $specificType, $iPath);

# now you can intersect what ever you want
# no need to do ASE for ASAT cases
# print "\n[[III]]. Intersect ASE Genes with $specificType ASE SNVs\n";
my $asarpType = 'AT';
my $asarpName = "Alternative Termination/Polyadenlylation";
print "\n[[III]]. Compare $asarpType ($asarpName) genes with $specificType SNVs\n";
# may use part of samplesAnlys.pl to get only common RI results
my $riGeneRef1 =$iAsarpGenesRef1->{$asarpType};
my $riGeneRef2 =$iAsarpGenesRef2->{$asarpType};

my ($iRiRef1, $iRiComRef, $iRiRef2) = intersectHashes($riGeneRef1, $riGeneRef2);

my $noRiGenes1 = keys %$iRiRef1;
my $noRiGenes2 = keys %$iRiRef2;
my $noRiComGenes = keys %$iRiComRef;
print "$specificType ASE SNVs common in genes with $asarpType ($asarpName):\n";
print "SAMPLE1_ONLY\tSAMPLES_COM\tSAMPLE2_ONLY\n";
print "$noRiGenes1\t$noRiComGenes\t$noRiGenes2\n";
outputDetailsByType($output, "$specificType.$asarpType.txt", $iRiRef1, $iRiRef2, $iRiComRef);

print "We should focus on $asarpName in C.A+ only and common\n";

# check whether the ASE SNVs in N.A+ are powerful but not significant
