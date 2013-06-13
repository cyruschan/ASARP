#!/usr/bin/perl
use warnings;
use strict;
require "anlysUtilities.pl"; #sub's for compSamples
use MyConstants qw( $asarpTags );
our @asarpChecks = split(' ', $asarpTags);

# set autoflush for error and output
select(STDERR);
$| = 1;
select(STDOUT);
$| = 1;

if(@ARGV<4){
  print <<EOT;

USAGE: perl $0 sample1 sample2 type output

This script gets the ASE/ASARP SNV cases from two samples
(datasets) and compare their common and uncommon cases at 
both gene and SNV levels.

For example, sample1 and sample2 can be tumor and normal
ASARP results, and one can find the common ASE genes shared
in the two samples. Another example can be ASARP results of
different cell compartments (e.g. cytosol VS nucleus), and
one can analyze the potential kinetics of ASE/ASARP in the 
two snapshots.

ARGUMENTS:
sample1/2	the main result file path of sample1/2
		in particular, the auxiliary predcition files
		e.g. .gene.prediction, .ase.prediction, etc.
		will be used, so make sure they are with sample1/2
type		the type of interest, namely "ASE" or "ASARP"
		(no quotes when input)

output		the output file path specified. Depending on the
		type chosen, output.ase (for "ASE") OR
		output.asarp (for "ASARP") AND
		output.as/ai/ass...
		will be output accordingly, specifying the overlap
		results of a particular type in ASE or ASARP

EOT
  exit;
}

my ($sample1, $sample2, $type, $output) = @ARGV;
#print " ($sample1, $sample2, $type, $output)\n";

$type = uc $type;
if($type ne "ASE" && $type ne "ASARP"){
  die "ERROR: type not supported: $type. ASE or ASARP expected\n";
}
print "Analyzing $type results\nSAMPLE1: $sample1\tSAMPLE2: $sample2\n";
open(my $FPO, ">", $output) or die "ERROR: cannot open output: $output\n";
print $FPO "samples\t$sample1\t$sample2\tcommon\n";
print $FPO "$type\tSAMPLE1_ONLY\tSAMPLE2_ONLY\tSAMPLES12_COMMON\n";

if($type eq "ASE"){
  my $resultType = 'ase.prediction';
  my ($ref1, $snvRef1) = getAseAll("$sample1.$resultType");
  my ($ref2, $snvRef2) = getAseAll("$sample2.$resultType");

  #get their intron/exon categories?
  # no, this is done in intronAseAnlys.pl specifically for introns

  # get intersection
  my ($s1Ref, $comRef, $s2Ref) = intersectHashes($ref1, $ref2);
  
  my %s1Ex = %$s1Ref; #sample 1 exclusive
  my %s2Ex = %$s2Ref; #sample 2 exclusive
  my %sCom = %$comRef; #common
  # brief: print this type to a one liner
  my $no1 = keys %s1Ex;
  my $no2 = keys %s2Ex;
  my $noCom = keys %sCom;
  my $oneLiner = join("\t", ("ASE", $no1, $no2, $noCom));
  print $FPO "$oneLiner\n"; 
  print "SAMPLE1_ONLY: $no1\tSAMPLE2_ONLY: $no2\tSAMPLES12_COMMON: $noCom\n";
    
  # details: output the overlapping details by type
  outputDetailsByType($output, "ASE", \%s1Ex, \%s2Ex, \%sCom);

}else{ # ASARP
  my $resultType = 'gene.prediction';
  my ($ref1) = getAsarpAll("$sample1.$resultType");
  my ($ref2) = getAsarpAll("$sample2.$resultType");

  # show results (debug)
  #for(@asarpChecks){
  #  showAsarpByType($ref1, $_);
  #  showAsarpByType($ref2, $_);
  #}
  # compare this two
  for(@asarpChecks){
    my $asType = $_;
    print "Comparing $asType\n"; 
    my ($s1Ref, $comRef, $s2Ref) = intersectHashes($ref1->{$_}, $ref2->{$_});
  
    my %s1Ex = %$s1Ref; #sample 1 exclusive
    my %s2Ex = %$s2Ref; #sample 2 exclusive
    my %sCom = %$comRef; #common
    # brief: print this type to a one liner
    my $no1 = keys %s1Ex;
    my $no2 = keys %s2Ex;
    my $noCom = keys %sCom;
    my $oneLiner = join("\t", ($asType, $no1, $no2, $noCom));
    print $FPO "$oneLiner\n"; 
    
    print "SAMPLE1_ONLY: $no1\tSAMPLE2_ONLY: $no2\tSAMPLES12_COMMON: $noCom\n";
    # details: output the overlapping details by type
    outputDetailsByType($output, $asType, \%s1Ex, \%s2Ex, \%sCom);
  }
    
}
close($FPO);

