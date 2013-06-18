#!/usr/bin/perl
use warnings;
use strict;
require "anlysUtilities.pl"; #sub's for compSamples

# set autoflush for error and output
select(STDERR);
$| = 1;
select(STDOUT);
$| = 1;

if(@ARGV < 10){

  print <<EOT;

USAGE: perl $0 specific_type specific_output1 config1 param1 result1 specific_output2 config2 param2 result2 output [asarp_path]

This pipeline combines several perl scripts to analyze the specific ASE SNVs
from two samples (datasets), especially in the context where the two samples 
represent different celluar compartments, e.g. cytocol and nucleuse.

In the input, 1, 2 indicate the corresponding in/output for samples (datasets) 
1 and 2, respectively. 

specific_type	specific type of ASE SNVs to be anlayzed, e.g. 3'UTR or 5'UTR

specific_output	the specific ASE snv intermediate output prefix
		severval files with suffixes: .ase .ase_distri.list, etc.
		will be generated as the intermediate output for one sample
		NOTE: make sure the file names are **DIFFERENT** between 
		specific_output1 and specific_output2
config		the .config file which is required for ASARP pipeline entry
param		the .param file, optional in the ASARP pipeline, is required
		in this pipeline
result		the result (main output summary file) of ASARP
		the corresponding .ase.prediction and .gene.prediction will
		be needed (assumed to be present with the main file)
output		the finaly output file of this analysis pipeline

OPTIONAL:
asarp_path	the -I path for the ASARP pipeline folder if it is run outside
		the original folder. Because this pipeline executes several
		scripts in ASARP: aseSnvs.pl, snp_distri.pl, etc., the ASARP
		path is needed for them.

EOT

exit;
}

my ($specificType, $intron1, $config1, $param1, $result1, $intron2, $config2, $param2, $result2, $output, $iPath) = @ARGV;
if(!defined($iPath)){
  $iPath = ".";
}
if(substr($iPath, -1, 1) ne "/" ){
  $iPath .= "/"; #add the path slash
}
#my $specificType = 'INTRON';
# focus on the specific SNVs only and compare them with the 2 sample results
print "Analyzing $specificType ASE SNVs\n";


print "\n[[I]]. Obtain all $specificType ASE SNVs from sample 1\n";
my ($iAseGenesRef1, $iAseSnvsRef1, $iAsarpGenesRef1, $iAsarpSnvsRef1) = specificAsePipeline($intron1, $config1, $param1, $result1, $specificType, $iPath);

print "\n[[II]]. Obtain all $specificType ASE SNVs from sample 2\n";
my ($iAseGenesRef2, $iAseSnvsRef2, $iAsarpGenesRef2, $iAsarpSnvsRef2) = specificAsePipeline($intron2, $config2, $param2, $result2, $specificType, $iPath);

# now you can intersect what ever you want
print "\n[[III]]. Intersect ASE Genes with $specificType ASE SNVs\n";

print "ASE Gene-level\n";
my ($iAseGNo1, $iAseGNo2) = hashRefNo($iAseGenesRef1, $iAseGenesRef2);
print "SAMPLE1: $iAseGNo1 SAMPLE2: $iAseGNo2\n";
my ($iAse1, $iCom, $iAse2) = intersectHashes($iAseGenesRef1, $iAseGenesRef2);
print "ASE genes containing $specificType ASE SNVs:\n";
my ($noAG1, $noComAG, $noAG2) = hashRefNo($iAse1, $iCom, $iAse2);
print "SAMPLE1_ONLY\tSAMPLES_COM\tSAMPLE2_ONLY\n";
print "$noAG1\t$noComAG\t$noAG2\n";  
outputDetailsByType($output, "$specificType.ASE.txt", $iAse1, $iAse2, $iCom);

print "ASE SNV-level\n";
my ($iAseNo1, $iAseNo2) = hashRefNo($iAseSnvsRef1, $iAseSnvsRef2);
print "SAMPLE1: $iAseNo1 SAMPLE2: $iAseNo2\n";
my ($iAseS1, $iComS, $iAseS2) = intersectHashes($iAseSnvsRef1, $iAseSnvsRef2);
print "$specificType ASE SNVs in ASE genes:\n";
my ($noAS1, $noComAS, $noAS2) = hashRefNo($iAseS1, $iComS, $iAseS2);
print "SAMPLE1_ONLY\tSAMPLES_COM\tSAMPLE2_ONLY\n";
print "$noAS1\t$noComAS\t$noAS2\n";  
outputDetailsByType($output, "$specificType.ASE_snvs.txt", $iAseS1, $iAseS2, $iComS);

# what to intersect then?
my $asarpType = 'RI';
my $asarpName = "Retained Intron";
print "\n[[IV]]. Compare $asarpType ($asarpName) genes with $specificType SNVs\n";
# may use part of samplesAnlys.pl to get only common RI results
my $riGeneRef1 =$iAsarpGenesRef1->{$asarpType};
my $riGeneRef2 =$iAsarpGenesRef2->{$asarpType};

my ($riGNo1, $riGNo2) = hashRefNo($riGeneRef1, $riGeneRef2);
print "SAMPLE1: $riGNo1 SAMPLE2: $riGNo2\n";
my ($iRiRef1, $iRiComRef, $iRiRef2) = intersectHashes($riGeneRef1, $riGeneRef2);
my ($noRiGenes1, $noRiComGenes, $noRiGenes2) = hashRefNo($iRiRef1, $iRiComRef, $iRiRef2);
print "$specificType ASE SNVs common in genes with $asarpType ($asarpName):\n";
print "SAMPLE1_ONLY\tSAMPLES_COM\tSAMPLE2_ONLY\n";
print "$noRiGenes1\t$noRiComGenes\t$noRiGenes2\n";
outputDetailsByType($output, "$specificType.$asarpType.txt", $iRiRef1, $iRiRef2, $iRiComRef);


