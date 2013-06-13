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
my ($finalP, $modiP, $bhP, $byP) = (-1, -1, -1, -1);
if(defined($FDRCUTOFF)){
  print "Calculating the Chi-Squared Test p-value cutoff for FDR <= $FDRCUTOFF...\n";
  if(defined($SNVPCUTOFF)){
    print "NOTE: user-provided p-value in config: $SNVPCUTOFF is ignored.\n";
  }
  ($finalP, $modiP, $bhP, $byP) = fdrControl($pRef, $FDRCUTOFF, 1); #1--verbose
  $SNVPCUTOFF = $finalP;
  print "Chi-Squared Test adjusted p-value cutoff: $SNVPCUTOFF\n\n";
}else{
  print "FDR control NOT used. User set Chi-Squared Test p-value cutoff: $SNVPCUTOFF\n";
}

#plot p-values
plotPvalues($pRef, $FDRCUTOFF, $modiP, $bhP, $byP);

######################################################################################
## sub-routines
sub plotPvalues
{
  my ($pRef, $cutoff, $modifiedP, $bhP, $byP) = @_;
  print "Plotting different FDR methods at threshold: $cutoff\n";
  print "modified Fdr: $modifiedP\nBH Fdr: $bhP\nBY Fdr: $byP\n";
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
  my $qpos = $R->get('qpos');
  my $qval = $R->get('fdr$qval[qpos]');
  $R->run("lpos <- max(which(fdr\$lfdr <= $cutoff))");
  my $lpos = $R->get('lpos');
  my $lfdr = $R->get('fdr$lfdr[lpos]');
  my $pFdr = $R->get("$rVar\[qpos\]");
  my $plfdr = $R->get("$rVar\[lpos\]");
 
  #getPrintInR($R, $rVar);
  #getPrintInR($R, 'fdr$qval');
  #getPrintInR($R, 'fdr$lfdr');

  print "fdrtool Fdr: $pFdr\n";
  print "fdrtool fdr: $plfdr\n";
  #print "fdrtool Fdr: $qval (p-value: $pFdr) pos: $qpos\n";
  #print "fdrtool fdr: $lfdr (p-value: $plfdr) pos: $lpos\n";
 
  plotInR($R, $outputFile, $rVar, $pFdr, $modifiedP, $byP, $plfdr, $bhP);
  my $upper = $qpos;
  if($upper+1 < @pList){ $upper += 1; } #leave some margin
  plotInR($R, "$outputFile.zoom.pdf", $rVar, $pFdr, $modifiedP, $byP, $plfdr, $bhP, $upper);
  
  simplePlotInR($R, "$outputFile.upper.pdf", $rVar, 0.8); # the uppder corner for modiFdr
  
  $R->stop;
}

sub getPrintInR
{
  my ($R, $var) = @_;
  my $ref = $R->get($var);
  my @list = @$ref;
  print "$var\n@list\n";
}


sub plotInR
{
  my ($R, $outputFile, $rVar, $pFdr, $modifiedP, $byP, $plfdr, $bhP, $upper) = @_;
  # do the plot
  #$R->run("png(filename=\"$outputFile\", width=5,height=5,units=\"in\", res = 600)");
  $R->run("pdf(file=\"$outputFile\", width=7,height=3.5)");
  $R->run("par(mar=c(7.2, 3.8, 0.2, 0.2))");

  my $settings = "xlab=\"indices (sorted)\", ylab=\"p-values\", cex.lab=0.8, cex=0.1"; #pch=\".\""; #, lwd=1";
  if(defined($upper)){
    $R->run("plot($rVar\[1:$upper\], $settings)"); 
    $upper = $R->get("$rVar\[$upper\]");
    $upper = $upper*0.9; # adjust the margin
    #print "Upper is $upper\n";
  }else{
    $upper = 0.8; # for legend positioning
    $R->run("plot($rVar, $settings)"); #\[1:qpos+5\])");
  }
  $R->run("abline(a=0, b=1/length($rVar), col=7, lty=1)"); # diagonal line
  $R->run("abline(h=$pFdr,col=2,lty=2)"); # Fdr (fdrtool)
  if($modifiedP >=0){
    $R->run("abline(h=$modifiedP,col=3,lty=3)"); # modified FDR
  }
  if($bhP >=0){
    $R->run("abline(h=$bhP,col=4,lty=4)"); # BH p.adjust method's p
  }
  if($byP >=0){
    $R->run("abline(h=$byP,col=5,lty=5)"); # BY method's p
  }
  $R->run("abline(h=$plfdr,col=6,lty=6)"); # local fdr (fdrtool)
  
  #$R->run('title(main="p-values")');
  $R->run('legend(1, '.$upper.', c("fdrtool(Fdr)","modified(Fdr)","BH(Fdr)","BY(Fdr)", "fdrtool(fdr)", "diagonal"), 
  cex=0.8, bg="white", col=c(2,3,4,5,6,7), lty=c(2,3,4,5,6,1))');

  $R->run("dev.off()");
  
}

sub simplePlotInR
{
  my ($R, $outputFile, $rVar, $lower) = @_;
  # do the plot
  #$R->run("png(filename=\"$outputFile\", width=5,height=5,units=\"in\", res = 600)");
  $R->run("pdf(file=\"$outputFile\", width=3.5,height=3.5)");
  $R->run("par(mar=c(4.2, 3.8, 0.2, 0.2))");
  
  my $settings = "xlab=\"indices (sorted)\", ylab=\"p-values\", cex.lab=0.8, cex=0.1"; #pch=\".\""; #, lwd=1";
  $R->run("lowerPos <- min(which($rVar >= $lower))");
  $R->run("plot($rVar\[lowerPos:length($rVar)\], lowerPos:length($rVar), $settings)");
  $R->run("abline(a=0, b=length($rVar), col=7, lty=1)"); # diagonal line
  $R->run("dev.off()");
}
