#!/usr/bin/perl -w
use strict;

require "fileParser.pl"; #sub's for input annotation files
require "snpParser.pl"; #sub's for snps
use MyConstants qw( $CHRNUM $supportedList $supportedTags );

# set autoflush for error and output
select(STDERR);
$| = 1;
select(STDOUT);
$| = 1;

if(@ARGV < 2){
  print <<EOT;
  
USAGE: perl $0 output_file config_file [optional: parameter_file]

To get the ASE and Powerful SNVs from the input SNV file (in config_file) along with their Chi-Squared Goodness of Fit Test p-values
The config_file (same for the ASARP pipeline) is used to obtain the input SNV file path. Other files/paths are read through and ignored.
The criteria (powerful and ASE, i.e. p-value, cutoffs) are contained in the parameter_file (same for the ASARP pipelien). Other parameters are read through and ignored.

There will be two result files: 
"output_file.ase" with the ASE SNVs (a subset of the powerful SNVs), and 
"output_file.pwr" with the powerful SNVs

EOT
  exit;
}

# input arguments: $outputFile--output, $configs--configuration file for input files, $params--configuration file for parameters
my ($outputFile, $configs, $params) = getArgs(@ARGV); 
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
if(defined($FDRCUTOFF)){
  print "Calculating the Chi-Squared Test p-value cutoff for FDR <= $FDRCUTOFF...\n";
  if(defined($SNVPCUTOFF)){
    print "NOTE: user-provided p-value in config: $SNVPCUTOFF is ignored.\n";
  }
  $SNVPCUTOFF = fdrControl($pRef, $FDRCUTOFF, 1); #1--verbose
  print "Chi-Squared Test p-value cutoff: $SNVPCUTOFF\n\n";
}

my ($AP, $PP) = (undef, undef); #file handles
open($AP, ">", "$outputFile.ase") or die "ERROR: cannot open $outputFile.ase for ASE SNVs\n";
open($PP, ">", "$outputFile.pwr") or die "ERROR: cannot open $outputFile.pwr for Powerful SNVs\n";

if($STRANDFLAG){
  print "+ strand only powerful and ASE SNVs\n";
  my ($ase, $pwr) = outputAsePwrSnvsOneStrand($snpRef, $SNVPCUTOFF, $AP, $PP);
  print "- strand only powerful and ASE SNVs\n";
  my ($aseRc, $pwrRc) = outputAsePwrSnvsOneStrand($snpRcRef, $SNVPCUTOFF, $AP, $PP);
  $ase+=$aseRc;
  $pwr+=$pwrRc;
  print "\nCombining both + and - strands, there are in total $ase ASE SNVs out of $pwr powerful SNVs\n";
}else{

  outputAsePwrSnvsOneStrand($snpRef, $SNVPCUTOFF, $AP, $PP);
}

close($AP);
close($PP);
print "Done\n";

sub outputAsePwrSnvsOneStrand{
	my ($snpRef, $SNVPCUTOFF, $AP, $PP) = @_;
	#basic statistics for powerful SNVs but not genes
	my $powSnvCnt = 0;
	my $aseSnvCnt = 0;

	for(my $i=1; $i<=$CHRNUM; $i++){
	  my ($snpChrRef) = getListByKeyChr($snpRef, 'powSnps', $i);
	  my %powSnps = %$snpChrRef;
	  my $powCntChr = keys %powSnps;  
	  my $aseCntChr = 0;
	  
	  #basic statistics
	  $powSnvCnt += $powCntChr; # keys %powSnps;
	  for(keys %powSnps){
	    my @allSnpInfo = split(';', $powSnps{$_}); #separate by ;, if there are multiple snps at the same position
	    for(@allSnpInfo){
	      my ($p, $pos, $alleles, $snpId, $refAl, $altAl) = getSnpInfo($_);
	      #get allelic ratios
	      my $r = $refAl/($refAl+$altAl);
	      my $toOutput = join(" ", formatChr($i), $pos, $alleles, $snpId, $r, "$refAl:$altAl:0", $p);
	      if($p <= $SNVPCUTOFF){
		$aseCntChr += 1; #each SNV **location** added once
		print $AP $toOutput."\n"; #ase
	      }
	      print $PP $toOutput."\n"; #powerful, which is a superset of ase
	      last; # I assume there can be multiple SNVs at the same position, but usually ppl don't
	      #therefore, a "last" is shot to get only the first one in the position
	    }
	  }
	  if($powCntChr >0){  printChr($i); print" has $powCntChr powerful SNVs and $aseCntChr SNVs\n";  }
	  $aseSnvCnt += $aseCntChr;
	}
	print "\nThere are in total $aseSnvCnt ASE SNVs out of $powSnvCnt powerful SNVs\n";
	return ($aseSnvCnt, $powSnvCnt);
}

=head1 NAME

aseSnvs.pl -- To get the ASE and Powerful SNVs from the input SNV file (in config_file) along with their Chi-Squared Goodness of Fit Test p-values. It is an introductory application script to get familiar with the ASARP pipeline (L<asarp>). 

=head1 SYNOPSIS

  perl aseSnvs.pl output_file config_file [optional: parameter_file]

The C<config_file> (same in the ASARP pipeline) is used to obtain the input SNV file path. Other files/paths are read through and ignored.
The criteria (powerful and ASE, i.e. p-value, cutoffs) are contained in the C<parameter_file> (same in the ASARP pipelien). Other parameters are read through and ignored.

=head2 INPUTS

The application counts the powerful SNVs and ASE SNVs. It accepts the same configuration and parameter file formats as the full ASARP pipeline. See below and L<asarp> for more details.

C<config_file> is the input configuration file which contains all the input file keys and their paths. The format is <key>tab<path>. Line starting with # are comments. Example: F<../default.config>

C<parameter_file> is the parameter configuration file which contains all the thresholds and cutoffs, e.g. p-value cuttoffs and bounds for absolute allelic ratio difference. The format of each line is <parameter>tab<value>. Lines starting with # are comments. It is optional and the default is: F<../default.param>

See the format descriptions in L<snpParser> for more details about the input SNV list content.

=head2 OUTPUTS

There will be two result files: 
"output_file.ase" with the ASE SNVs (a subset of the powerful SNVs), and 
"output_file.pwr" with the powerful SNVs

Format (space dilimited): chromosome position ref>alt dbSNP_ID allelic_ratio ref#:alt#:0 p-value
ref indicates the reference allele, and alt the alternative allele. ref# and alt# are their respective counts (wrong nucleotides are not used so 0 is output in ref#:alt#:0). allelic_ratio = ref#/(ref#+alt#). p-value is the Chi-Squared Test p-value. Examples:

 chr1 68591173 A>G rs1046835 0.651898734177215 103:55:0 0.0001341704
 chr5 96111371 G>A rs13160562 0.640350877192982 73:41:0 0.00272584
 chr5 96215303 C>T na 0.96969696969697 96:3:0 9.029544e-21

=head1 SEE ALSO

L<asarp>, L<snp_distri>, L<fileParser>, L<snpParser>, L<MyConstants>

=head1 COPYRIGHT

This pipeline is free software; you can redistribute it and/or modify it given that the related works and authors are cited and acknowledged.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.

=head1 AUTHOR

Cyrus Tak-Ming CHAN

Xiao Lab, Department of Integrative Biology & Physiology, UCLA

=cut
