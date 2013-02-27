#!/usr/bin/perl -w
use strict;

require "fileParser.pl"; #sub's for input annotation files
require "snpParser.pl"; #sub's for snps
require "bedHandler.pl"; #for bed.sam

if(@ARGV < 2){
  print "USAGE: perl $0 output_file config_file [optional: parameter_file]\n";
  exit;
}

# input arguments: $outputFile--output, $configs--configuration file for input files, $params--configuration file for parameters
my ($outputFile, $configs, $params) = getArgs(@ARGV); 
my ($snpF, $bedF, $rnaseqF, $xiaoF, $splicingF, $estF) = getRefFileConfig($configs); # input annotation/event files
my ($POWCUTOFF, $SNVPCUTOFF, $ASARPPCUTOFF, $NEVCUTOFFLOWER, $NEVCUTOFFUPPER, $ALRATIOCUTOFF) = getParameters($params); # parameters

# read the transcript annotation file
my $transRef = readTranscriptFile($xiaoF);
#printListByKey($transRef, 'trans'); #utility sub: show transcripts (key: trans)
my ($genesRef, $geneNamesRef) = getGeneIndex($transRef); #get indices of gene transcript starts and gene names
my $altRef = getGeneAltTransEnds($transRef); #get alternative initiation/termination (AI/AT) events from transcripts
#printAltEnds($altRef); #utility sub: show AI/AT events derived from transcripts

# read all annotations, optionally rna-seq and est, splicing events and compile them
my $allEventsListRef = readAllEvents($splicingF, $rnaseqF, $estF, $transRef, $geneNamesRef);
my $splicingRef = compileGeneSplicingEvents($genesRef, values %$allEventsListRef); #compile events from different sources

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
my $outputASE = $outputFile.'.ase.prediction';
my $outputSnv = $outputFile.'.snv.prediction';
my $outputGene = $outputFile.'.gene.prediction';
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

=head1 NAME

asarp.pl -- The main application script, i.e. the entry program, of the ASARP pipeline.

=head1 SYNOPSIS

Look at the source: F<../asarp.pl> and it is self-explanatory. There are basically 3 steps:

1. parse the input files and compile alternative mRNA processing events. see L<fileParser>

2. get the SNVs and match them with the events. see L<snpParser>

3. process ASARP (including ASE) patterns and output the formatted results. see source and L<snpParser>

=head1 DESCRIPTION

The methodology in detail is explained in the paper: 

Li G, Bahn JH, Lee JH, Peng G, Chen Z, Nelson SF, Xiao X. Identification of allele-specific alternative mRNA processing via transcriptome sequencing, Nucleic Acids Research, 2012, 40(13), e104

and its Supplementary Materials

G<img/Intro.png>

See http://nar.oxfordjournals.org/content/40/13/e104 for more details of the paper.


=head1 SEE ALSO

L<fileParser>, L<snpParser>, L<MyConstants>

=head1 COPYRIGHT

This pipeline is free software; you can redistribute it and/or modify it given that the related works and authors are cited and acknowledged.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.

=head1 AUTHOR

Cyrus Tak-Ming CHAN

Xiao Lab, Department of Integrative Biology & Physiology, UCLA

=cut
