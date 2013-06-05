#!/usr/bin/perl -w
use strict;

require "fileParser.pl"; #sub's for input annotation files
require "snpParser.pl"; #sub's for snps
require "bedHandler.pl"; #for bed.sam

# set autoflush for error and output
select(STDERR);
$| = 1;
select(STDOUT);
$| = 1;

if(@ARGV < 2){
  print "USAGE: perl $0 output_plots config_file [optional: parameter_file]\n";
  print "To plot and analyze ASE SNV p-values\n";
  exit;
}

# input arguments: $outputFile--output, $configs--configuration file for input files, $params--configuration file for parameters
my ($outputFile, $configs, $params) = getArgs(@ARGV); 
#strand-specific setting is more related to file configs (data dependent)
my ($snpF, $bedF, $rnaseqF, $xiaoF, $splicingF, $estF, $STRANDFLAG) = getRefFileConfig($configs); # input annotation/event files
my ($POWCUTOFF, $SNVPCUTOFF, $FDRCUTOFF, $ASARPPCUTOFF, $NEVCUTOFFLOWER, $NEVCUTOFFUPPER, $ALRATIOCUTOFF) = getParameters($params); # parameters


# if strand-specific flag is set, need to get two separate SNV lists (+ and - respectively)
# the extra Rc (reference complement) references are for - strand if $STRANDFLAG is set
my ($snpRef, $snpRcRef, $pRef, $pRcRef) = (undef, undef, undef, undef);

if($STRANDFLAG){
  ($snpRef, $pRef) = initSnp($snpF, $POWCUTOFF, '+');
  ($snpRcRef, $pRcRef) = initSnp($snpF, $POWCUTOFF, '-');
  #print "SNV List: +\n";  printListByKey($snpRef, 'powSnps');  printListByKey($snpRef, 'snps');
  #print "SNV List: -\n";  printListByKey($snpRcRef, 'powSnps');  printListByKey($snpRcRef, 'snps');
  
  my @pList = @$pRef;
  my @pRcList = @$pRcRef;
  #merge two lists
  my @joinList = (@pList, @pRcList);
  $pRef = \@joinList; # renew the reference to be the full list
}else{
  ($snpRef, $pRef) = initSnp($snpF, $POWCUTOFF);
  #print "SNV List:\n";
  #printListByKey($snpRef, 'powSnps');  #printListByKey($snpRef, 'snps');
}

# suggested, get the Chi-Squared Test p-value cutoff from FDR ($FDRCUTOFF)
my ($finalP, $modiP, $byP) = (-1, -1, -1);
if(defined($FDRCUTOFF)){
  print "Calculating the Chi-Squared Test p-value cutoff for FDR <= $FDRCUTOFF...\n";
  if(defined($SNVPCUTOFF)){
    print "NOTE: user-provided p-value in config: $SNVPCUTOFF is ignored.\n";
  }
  ($finalP, $modiP, $byP) = fdrControl($pRef, $FDRCUTOFF, 1); #1--verbose
  $SNVPCUTOFF = $finalP;
  print "Chi-Squared Test adjusted p-value cutoff: $SNVPCUTOFF\n\n";
}else{
  print "FDR control NOT used. User set Chi-Squared Test p-value cutoff: $SNVPCUTOFF\n";
}

#plot p-values
plotPvalues($pRef, $FDRCUTOFF, $modiP, $byP);

######################################################################################
## sub-routines
sub plotPvalues
{
  my ($pRef, $cutoff, $modifiedP, $byP) = @_;
  my @pList = @$pRef;
  my $pListSize = @pList; #size
  @pList = sort{$a<=>$b}@pList; #sorted
  my $pSize = @pList;
  # Create a communication bridge with R and start R
  my $R = Statistics::R->new();

  my $rVar = "plist";
  passListToR($R, $pRef, $rVar);
  # load "fdrtool" library and p-values
  $R->run('library("fdrtool")'); 
  
  $R->run("fdr = fdrtool($rVar, statistic=\"pvalue\")");
  $R->run("qpos <- max(which(fdr\$qval <= $cutoff))");
  my $qval = $R->get('fdr$qval[qpos]');
  $R->run("lpos <- max(which(fdr\$lfdr <= $cutoff))");
  my $lfdr = $R->get('fdr$lfdr[lpos]');
  my $pFdr = $R->get("$rVar\[qpos\]");
  my $plfdr = $R->get("$rVar\[lpos\]");
  
  print "Fdr q-value: $qval with p-value: $pFdr\n\n";
  print "Local fdr  : $lfdr with p-value: $plfdr\n\n";
  
  # do the plot
  $R->run("png(filename=\"$outputFile\")"); #, width=3.25,height=3.25,units=\"in\", res = 1200)");
  $R->run("plot($rVar)"); #\[1:qpos+5\])");
  $R->run("abline(h=$pFdr,col=3,lty=2)"); # Fdr (BH method)
  $R->run("abline(h=$plfdr,col=4,lty=3)"); # local FDR
  if($byP >0){
    $R->run("abline(h=$byP,col=5,lty=4)"); # BY method's p
  }
  if($modifiedP >0){
    $R->run("abline(h=$modifiedP,col=2,lty=1)"); # modified FDR
  }
  $R->run('title(main="p-values")');
  $R->run("dev.off()");
  
  $R->stop;
}
