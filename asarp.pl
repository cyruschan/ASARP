#!/usr/bin/perl -w
use strict;

require "fileParser.pl"; #sub's for input annotation files
require "snpParser.pl"; #sub's for snps
require "bedHandler.pl"; #for bed.sam

# global variables as potential parameters here
our $POWCUTOFF = 20; #threshold of powerful snvs
our $SNVPCUTOFF = 0.014;  #p-value cutoff for the Chi-squared test
our $ASARPPCUTOFF = 0.05; #p-value cutoff for the Fisher'exact test
our $NEVCUTOFF = 0.8; #normalized expression value (NEV) cutoff for alternative regions/exons
our $ALRATIOCUTOFF = 0.2; #allelic ratio difference cutoff for AS SNVs

#my $defaultConfig='test.config';
my $defaultConfig='default.config';
my ($snpF, $bedF, $genomeF, $rnaseqF, $xiaoF, $splicingF, $estF) = getRefFileConfig($defaultConfig); 

#my ($bedRef) = readBedByChr($bedF, $genomeF, 5);
#simpleTestReads($bedF, $genomeF, 5);

my $transRef = readTranscriptFile($xiaoF);
#printListByKey($transRef, 'trans');
my $genesRef = getGeneIndex($transRef);
my $altRef = getGeneAltTransEnds($transRef);
#printAltEnds($altRef);

# read all events from annotations and optionally rna-seq and est event files
my $allEventsListRef = readAllEvents($splicingF, $rnaseqF, $estF, $transRef);
my $splicingRef = compileGeneSplicingEvents($genesRef, values %$allEventsListRef);

my $snpRef = initSnp($snpF, $POWCUTOFF);
#print "SNV List:\n";
#printListByKey($snpRef, 'powSnps');
#printListByKey($snpRef, 'snps');

my $geneSnpRef = setGeneSnps($snpRef, $transRef);
#print "Significant Snvs: \n";
#printGetGeneSnpsResults($geneSnpRef,'gPowSnps', $snpRef,'powSnps', $PCUTOFF);
#print "Ordinary Snvs: \n";
#printGetGeneSnpsResults($geneSnpRef,'gSnps', $snpRef,'snps', 1);

my ($snpEventsRef) = setSnpEvents($geneSnpRef, $altRef, $splicingRef); #match snps with events
#print "Pow Alt: \n";
#printSnpEventsResultsByType($snpEventsRef,'powSnpAlt'); 
#print "\n";

#print "Pow Sp: \n";
#printSnpEventsResultsByType($snpEventsRef,'powSnpSp'); 
#print "Ord Sp: \n";
#printSnpEventsResultsByType($snpEventsRef,'snpSp'); 

print "\n\nCalculating NEV\n";
my ($snpsNevRef) = filterSnpEventsWithNev($snpRef, $geneSnpRef, $snpEventsRef, $bedF, $genomeF, $allEventsListRef, $NEVCUTOFF); 
print "Pow NEV Alt: \n";
printSnpEventsResultsByType($snpsNevRef,'nevPowSnpAlt'); 
print "\n";
print "NEV Alt: \n";
printSnpEventsResultsByType($snpsNevRef,'nevSnpAlt'); 
print "\n\n";
print "Pow NEV Sp: \n";
printSnpEventsResultsByType($snpsNevRef,'nevPowSnpSp'); 
print "\n";
print "NEV Sp: \n";
printSnpEventsResultsByType($snpsNevRef,'nevSnpSp'); 
print "\n";


print "processing ASE's\n";
my ($allAsarpsRef) = processASEWithNev($snpRef, $geneSnpRef, $snpsNevRef, $SNVPCUTOFF, $ASARPPCUTOFF, $NEVCUTOFF, $ALRATIOCUTOFF);

print "\n";
my $outputSnv = 'snv.predicted.txt';
my $outputGene ='gene.predicted.txt';
outputRawASARP($allAsarpsRef, 'ASEgene');
outputRawASARP($allAsarpsRef, 'ASARPgene'); #, $outputGene);
outputRawASARP($allAsarpsRef, 'ASARPsnp'); #, $outputSnv);

