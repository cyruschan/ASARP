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
my $asarpType = 'RI';
my $asarpName = "Retained Intron";
print "Analyzing $specificType ASE SNVs\n";
my @steps = (
"[[I]]. Obtain all $specificType ASE SNVs from sample 1",
"[[II]]. Obtain all $specificType ASE SNVs from sample 2",
"[[III]]. Intersect ASE Genes with $specificType ASE SNVs",
"[[IV]]. Compare $asarpType genes with $specificType SNVs"
);

open(my $fp, ">", $output) or die "ERROR: cannot open output file: $output\n";

print "\n$steps[0]\n";
print $fp "$steps[0]\n"; # to output file
my ($iAseGenesRef1, $iAseSnvsRef1, $iAsarpGenesRef1, $iAsarpSnvsRef1, $aseGenesRef1, $exclGenesRef1) = specificAsePipeline($intron1, $config1, $param1, $result1, $specificType, $iPath);

print "\n$steps[1]\n";
print $fp "$steps[1]\n"; # to output file
my ($iAseGenesRef2, $iAseSnvsRef2, $iAsarpGenesRef2, $iAsarpSnvsRef2, $aseGenesRef2, $exclGenesRef2) = specificAsePipeline($intron2, $config2, $param2, $result2, $specificType, $iPath);

# now you can intersect what ever you want
print "\n$steps[2]\n";
print $fp "$steps[2]\n"; # to output file
my $outAse = ""; 
$outAse .= "ASE Gene-level\n";
my ($iAseGNo1, $iAseGNo2) = hashRefNo($iAseGenesRef1, $iAseGenesRef2);
$outAse .= "SAMPLE1: $iAseGNo1 SAMPLE2: $iAseGNo2\n";
my ($iAse1, $iCom, $iAse2) = intersectHashes($iAseGenesRef1, $iAseGenesRef2);
$outAse .= "Intersect: ASE genes containing $specificType ASE SNVs:\n";
my ($noAG1, $noComAG, $noAG2) = hashRefNo($iAse1, $iCom, $iAse2);
$outAse .= "SAMPLE1_ONLY\tSAMPLES_COM\tSAMPLE2_ONLY\n$noAG1\t$noComAG\t$noAG2\n"; 
outputDetailsByType($output, "$specificType.ASE.txt", $iAse1, $iAse2, $iCom);

print "SAMPLE1 $specificType ASE SNVs in SAMPLE2 ASE genes (not in SAMPLE1 ASE genes)\n";
my ($ex1, $ex1Ase2, $ase2) = intersectHashes($exclGenesRef1, $aseGenesRef2);
my ($noEx1Ase2) = hasRefNo($ex1Ase2); 
print "There are $noEx1Ase2\n";
my %hs = %$ex1Ase2;
for(keys %hs){
  print "$_\n$hs{$_}\n";
}

print "SAMPLE2 $specificType ASE SNVs in SAMPLE1 ASE genes (not in SAMPLE2 ASE genes)\n";
my ($ex2, $ex2Ase1, $ase1) = intersectHashes($exclGenesRef2, $aseGenesRef1);
my ($noEx2Ase1) = hasRefNo($ex2Ase1); 
print "There are $noEx2Ase1\n";
my %hs2 = %$ex2Ase1;
for(keys %hs2){
  print "$_\n$hs2{$_}\n";
}

$outAse .= "ASE SNV-level\n";
my ($iAseNo1, $iAseNo2) = hashRefNo($iAseSnvsRef1, $iAseSnvsRef2);
$outAse .= "SAMPLE1: $iAseNo1 SAMPLE2: $iAseNo2\n";
my ($iAseS1, $iComS, $iAseS2) = intersectHashes($iAseSnvsRef1, $iAseSnvsRef2);
$outAse .=  "Intersect: $specificType ASE SNVs in ASE genes:\n";
my ($noAS1, $noComAS, $noAS2) = hashRefNo($iAseS1, $iComS, $iAseS2);
$outAse .= "SAMPLE1_ONLY\tSAMPLES_COM\tSAMPLE2_ONLY\n$noAS1\t$noComAS\t$noAS2\n";  
outputDetailsByType($output, "$specificType.ASE_snvs.txt", $iAseS1, $iAseS2, $iComS);
print "$outAse\n";
print $fp "$outAse\n";

# what to intersect then?
print "\n$steps[3]\n";
print $fp "$steps[3]\n";
# may use part of samplesAnlys.pl to get only common RI results
my $riGeneRef1 =$iAsarpGenesRef1->{$asarpType};
my $riGeneRef2 =$iAsarpGenesRef2->{$asarpType};

my $outAsarp = "";
$outAsarp .= "$specificType SNVs in $asarpType Gene-level\n";
my ($riGNo1, $riGNo2) = hashRefNo($riGeneRef1, $riGeneRef2);
$outAsarp .= "SAMPLE1: $riGNo1 SAMPLE2: $riGNo2\n";
my ($iRiRef1, $iRiComRef, $iRiRef2) = intersectHashes($riGeneRef1, $riGeneRef2);
$outAsarp .= "Intersect: $asarpType genes\n";
my ($noRiGenes1, $noRiComGenes, $noRiGenes2) = hashRefNo($iRiRef1, $iRiComRef, $iRiRef2);
$outAsarp .= "SAMPLE1_ONLY\tSAMPLES_COM\tSAMPLE2_ONLY\n$noRiGenes1\t$noRiComGenes\t$noRiGenes2\n";
outputDetailsByType($output, "$specificType.$asarpType.txt", $iRiRef1, $iRiRef2, $iRiComRef);

print "$outAsarp\n";
print $fp "$outAsarp\n";

# End of the program
print $fp "\nALL STEPS FINISHED\n";
close($fp);
print "\nALL STEPS FINISHED\n";

