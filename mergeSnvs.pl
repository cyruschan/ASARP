#!/usr/bin/perl
use warnings;
use strict;
use Statistics::R; #interact with R
require 'snpParser.pl';
require "readUtilities.pl";

if(@ARGV < 4){
  print <<EOT;

USAGE perl $0 prefix suffix mono=x output [isStrandSp]
merge SNVs of different chromosomes into one SNV list
the equivalent can be achieved using cat and grep commands
for strand-specific cases '+\$' and '\\-\$' can get all SNVs from the
+ and - strands respectively
prefix	the folder path and prefix string for all chr* SNV files
	e.g. "/home/N.A+/" for all chr* SNVs in that folder
suffix	the suffix after chr[1..22/X/Y/M]
mono=	options (no space): 
	"mono=0": allow mono-allelic SNVs; 
	"mono=x": discard SNVs with any allele < x reads; 
	useful when genotype is unknown and dbSNPs are used
	where mono-allelic SNVs would be called
output	the output file name for the merged SNVs
	note that output.plus and output.minus are by-product
	outputs for + and - strands for strand-specific cases
	make sure the output and chr* SNV file names do not conflict
OPTIONAL
isStrandSp	set 1 if the SNV files are strand-specific
		default: non-strand-specific

EOT
  exit;
}

my ($prefix, $suffix, $mono, $output, $isStranded) = @ARGV;
$mono = lc $mono;
my $isMono = 0;
if($mono =~ /mono=(\d+)/){
  $isMono = $1; # matched
}else{
  die "ERROR: mono=x expected: $mono input by user\n";
}

if($isMono >= 0){
  if($isMono == 0){
    print "Note: Allow mono-allelic SNVs (suppose heterozyous SNVs are known from genomic sequencing)\n";
  }else{
    print "Note: Discard SNVs with any allele (reference or alternative) < $isMono reads (in the case SNVs are called from RNA-Seq)\n";
  }
}else{
  die "ERROR: mono= should be non-negative ($mono input)\n";
}
my ($rat, $cov, $refcnt, $altcnt) = mergeSnvs($prefix, $suffix, $isMono, $output, $isStranded);
  
my $n = @$rat;
my $m = @$cov;
print "$n ratios; $m coverages;\ntotal reference reads: $refcnt; alternative reads: $altcnt;\noutput $output\n";

scatterPlot($rat, $cov, "$output.scatt.pdf");
histogramPlot($rat, $refcnt, $altcnt, "$output.hist.pdf");



#########################################################################
###  sub-routines

sub scatterPlot
{
  my ($rat, $cov, $output) = @_;
  my $R = Statistics::R->new();
  #plot(x,y, main="PDF Scatterplot Example", col=rgb(0,100,0,50,maxColorValue=255), pch=16)
  my ($rVarRatio, $rVarCov) = ('alratio', 'alcov'); # R variable names
  passListToR($R, $rat, $rVarRatio, 'nosort'); 
  passListToR($R, $cov, $rVarCov, 'nosort');
  print "Generating scatter plot to $output... ";
  $R->run("xnum <- length($rVarRatio)");
  my $xnum = $R->get('xnum');
  my $cmd = "plot($rVarRatio, $rVarCov, main=\"Scatterplot\", xlab=\"allelic ratio\", ylab=\"read coverage\", sub=\"$xnum SNVs in total\", col=rgb(0,0,100,50,maxColorValue=255), pch=16)";
  plotToPdf($R, $cmd, $output);
  $R->stop;
  print "Done\n";
}

sub histogramPlot
{
  my ($rat, $refcnt, $altcnt, $output) = @_;
  my $R = Statistics::R->new();
  my $rVarRatio = 'alratio';
  passListToR($R, $rat, $rVarRatio, 'nosort');

  # perform bionomial test to check the statistics:
  #
  $R->run("xmed <-median($rVarRatio)");
  my $xMed = $R->get('xmed');
  $R->run("xmean <-mean($rVarRatio)");
  my $xMean = $R->get('xmean'); 
 
  # get the counts
  #$R->run("scnt <- length(which($rVarRatio<0.5))");
  #$R->run("mcnt <- length(which($rVarRatio==0.5))");
  #$R->run("lcnt <- length(which($rVarRatio>0.5))");
  #my $scnt = $R->get('scnt');
  #my $mcnt = $R->get('mcnt');
  #my $lcnt = $R->get('lcnt');
  # split the .5 cnt
  #$R->run("mcnt2 <- as.integer(mcnt/2)");
  $R->run("pv = binom.test($refcnt, $refcnt+$altcnt, 0.5)\$p.value");
  my $pVal = $R->get('pv');
  print "median:\t$xMed\tmedian\t$xMean\tp-value:\t$pVal\n";

  $xMean = sprintf("%.2f", $xMean);
  $xMed = sprintf("%.2f", $xMed);
  #$pVal = sprintf("%.4f", $pVal);

  print "Generating histogram to $output... ";
  my $cmd = "hist($rVarRatio, breaks=50, main=\"Histogram\", xlab=\"allelic ratio\", sub=\"median: $xMed mean: $xMean p-value: $pVal\", ylab=\"frequency\")";
  plotToPdf($R, $cmd, $output);
  $R->stop;
  print "Done\n";
}

sub plotToPdf{
  my ($R, $cmd, $file) = @_; #pass the R object
  print "Plot to $file\n";
  # figure settings
  $R->run("pdf(file=\"$file\", height=6, width=6)");
  my $setMar = "par(mar=c(5.2, 4.8, 1.2, 0.8), cex.lab=1.6, cex.axis=1.8, lwd=4)";
  $R->run($setMar);
  #set line weight
  #print "$cmd\n";
  $R->run("$cmd");
  $R->run("dev.off()"); # end of plotting
  #boxplot(x,y, names=c(\"$caseLab\", \"$bgLab\"), ylim = c(ymin, ymax), ylab=\"log($yLab)\", $boxLwd)");
}


__END__

=head1 NAME

mergeSnvs.pl -- merge SNVs of different chromosomes into one SNV list
the equivalent can be achieved using cat and grep commands
for strand-specific cases '+$' and '\-$' can get all SNVs from the
+ and - strands respectively

After the full pre-processing, one can run ASARP for each chromosome's
SNVs, or merge all chromosomes' SNVs into a single SNV list using this
script (equivalent to using Unix cat and/or grep).

=head1 SYNOPSIS

USAGE:

 perl $0 prefix suffix mono=0/1 output [isStrandSp]

 prefix		the folder path and prefix string for all chr* SNV files
 		e.g. "/home/N.A+/" for all chr* SNVs in that folder
 suffix		the suffix after chr[1..22/X/Y/M]
 mono=		2 options: 
		"mono=1": allow mono-allelic SNVs; 
		"mono=0": discard mono-allelic SNVs; 
		useful when genotype is unknown and dbSNPs are used
 output		the output file name for the merged SNVs
 		note that output.plus and output.minus are by-product
 		outputs for + and - strands for strand-specific cases
 		make sure the output and chr* SNV file names do not conflict

 OPTIONAL
 isStrandSp	set 1 if the SNV files are strand-specific


=head1 SEE ALSO

L<rmDup>, L<mergeSam>, L<procReads>, L<procReadsJ>, L<asarp>

=head1 COPYRIGHT

This pipeline is free software; you can redistribute it and/or modify it given that the related works and authors are cited and acknowledged.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.

=head1 AUTHOR

Cyrus Tak-Ming CHAN

Xiao Lab, Department of Integrative Biology & Physiology, UCLA

=cut
