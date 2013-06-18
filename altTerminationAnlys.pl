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

USAGE: perl $0 3'UTR_output1 config1 param1 result1 3'UTR_output2 config2 param2 result2 output [asarp_path]

This pipeline combines several perl scripts to analyze allele-specific alternative
polyadenylation / 3' termination (AT) of the 3'UTR ASE SNVs which can be miRNA
targets. The analysis can be performed between two samples (datasets), e.g. two 
samples represent different celluar compartments: cytocol and nucleus.

In the input, 1, 2 indicate the corresponding in/output for samples (datasets) 
1 and 2, respectively. 

3'UTR_output	the 3'UTR ASE snv intermediate output prefix
		severval files with suffixes: .ase .ase_distri.list, etc.
		will be generated as the intermediate output for one sample
		NOTE: make sure the file names are **DIFFERENT** between 
		3'UTR_output1 and 3'UTR_output2
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
my $specificType = "3'UTR";
# focus on the specific SNVs only and compare them with the 2 sample results
print "Analyzing $specificType ASE SNVs\n";


print "\n[[I]]. Obtain all $specificType ASE SNVs from sample 1\n";
my ($iAseGenesRef1, $iAseSnvsRef1, $iAsarpGenesRef1, $iAsarpSnvsRef1) = specificAsePipeline($intron1, $config1, $param1, $result1, $specificType, $iPath);

print "\n[[II]]. Obtain all $specificType ASE SNVs from sample 2\n";
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
my %iRiGenes1 = %$iRiRef1; 
my %iRiGenes2 = %$iRiRef2; 
my %iRiComGenes = %$iRiComRef; 

my $noRiGenes1 = keys %iRiGenes1;
my $noRiGenes2 = keys %iRiGenes2;
my $noRiComGenes = keys %iRiComGenes;
print "$specificType ASE SNVs common in genes with $asarpType ($asarpName):\n";
print "SAMPLE1_ONLY\tSAMPLES_COM\tSAMPLE2_ONLY\n";
print "$noRiGenes1\t$noRiComGenes\t$noRiGenes2\n";
outputDetailsByType($output, "$specificType.$asarpType.txt", $iRiRef1, $iRiRef2, $iRiComRef);

print "We should focus on $asarpName in C.A+ only and common\n";


