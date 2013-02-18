#!/usr/bin/perl -w
use strict;

require "fileParser.pl"; #sub's for input annotation files
require "snpParser.pl"; #sub's for snps
require "bedHandler.pl"; #for bed.sam
my ($outputFile, $configs, $params) = getArgs(@ARGV); 
# all the annotations, events, reference files
my ($snpF, $bedF, $rnaseqF, $xiaoF, $splicingF, $estF) = getRefFileConfig($configs); 
# all the p-values cutoffs, thresholds
my ($POWCUTOFF, $SNVPCUTOFF, $ASARPPCUTOFF, $NEVCUTOFFLOWER, $NEVCUTOFFUPPER, $ALRATIOCUTOFF) = getParameters($params);

#my ($bedRef) = readBedByChr($bedF, 5);
#simpleTestReads($bedF, 5);

my $transRef = readTranscriptFile($xiaoF);
#printListByKey($transRef, 'trans');
my ($genesRef, $geneNamesRef) = getGeneIndex($transRef);
my $altRef = getGeneAltTransEnds($transRef);
#printAltEnds($altRef);

# read all events from annotations and optionally rna-seq and est event files
my $allEventsListRef = readAllEvents($splicingF, $rnaseqF, $estF, $transRef, $geneNamesRef);
my $splicingRef = compileGeneSplicingEvents($genesRef, values %$allEventsListRef);

my $snpRef = initSnp($snpF, $POWCUTOFF);
#print "SNV List:\n";
#printListByKey($snpRef, 'powSnps');
#printListByKey($snpRef, 'snps');

my $geneSnpRef = setGeneSnps($snpRef, $transRef);
#print "Significant Snvs: \n";
#printGetGeneSnpsResults($geneSnpRef,'gPowSnps', $snpRef,'powSnps', 1); #$SNVPCUTOFF);
#print "Ordinary Snvs: \n";
#printGetGeneSnpsResults($geneSnpRef,'gSnps', $snpRef,'snps', 1);


my ($snpEventsRef) = setSnpEvents($geneSnpRef, $altRef, $splicingRef); #match snps with events
#print "Pow Alt: \n";
#printSnpEventsResultsByType($snpEventsRef,'powSnpAlt'); 
#print "Snp Alt: \n";
#printSnpEventsResultsByType($snpEventsRef,'snpAlt'); 
#print "Pow Sp: \n";
#printSnpEventsResultsByType($snpEventsRef,'powSnpSp'); 
#print "Ord Sp: \n";
#printSnpEventsResultsByType($snpEventsRef,'snpSp'); 


print "\n\nCalculating NEV\n";
my ($snpsNevRef) = filterSnpEventsWithNev($snpRef, $geneSnpRef, $snpEventsRef, $bedF, $allEventsListRef, $NEVCUTOFFLOWER, $NEVCUTOFFUPPER); 
#print "Pow NEV Alt: \n";
#printSnpEventsResultsByType($snpsNevRef,'nevPowSnpAlt'); 
#print "NEV Alt: \n";
#printSnpEventsResultsByType($snpsNevRef,'nevSnpAlt'); 
#print "\n\n";
#print "Pow NEV Sp: \n";
#printSnpEventsResultsByType($snpsNevRef,'nevPowSnpSp'); 
#print "NEV Sp: \n";
#printSnpEventsResultsByType($snpsNevRef,'nevSnpSp'); 

print "processing ASE's\n";
my ($allAsarpsRef) = processASEWithNev($snpRef, $geneSnpRef, $snpsNevRef, $SNVPCUTOFF, $ASARPPCUTOFF, $ALRATIOCUTOFF);

print "\n";
my $outputASE = 'test.ase.predicted.txt';
my $outputSnv = 'test.snv.predicted.txt';
my $outputGene ='test.gene.predicted.txt';
outputRawASARP($allAsarpsRef, 'ASEgene', $outputASE);
outputRawASARP($allAsarpsRef, 'ASARPgene', $outputGene);
outputRawASARP($allAsarpsRef, 'ASARPsnp', $outputSnv);

my $allNarOutput = formatOutputVerNAR($allAsarpsRef);
if(defined($outputFile)){
  my $isOpen = open(my $fp, ">", $outputFile);
  if(!$isOpen){
    print "Warning: cannot open file: $outputFile to write the predicted results; will print on screen instead\n";
    print $allNarOutput;
  }else{
    print $fp $allNarOutput;
    close($fp);
  }
}else{
  print $allNarOutput;
}
