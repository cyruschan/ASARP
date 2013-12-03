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

if(@ARGV<8){
  print <<EOT;

USAGE: perl $0 sample1 snv1 sample2 snv2 type powerful p-value output

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
snv1/2		the RNA SNV files (allelic read counts) for the
		two samples

type		the type of interest, namely "ASE" or "ASARP"
		(no quotes when input)

powerful	the cutoff (e.g. 20) for powerful SNVs

p-value		the cutoff (e.g. 0.001) for further ASE SNV filtering

output		the output file path specified. Depending on the
		type chosen, output.ase (for "ASE") OR
		output.asarp (for "ASARP") AND
		output.as/ai/ass...
		will be output accordingly, specifying the overlap
		results of a particular type in ASE or ASARP

EOT
  exit;
}

my ($sample1, $snvFile1, $sample2, $snvFile2, $type, $powerful, $pcutoff, $output) = @ARGV;
#print " ($sample1, $sample2, $type, $output)\n";

$type = uc $type;
if($type ne "ASE" && $type ne "ASARP"){
  die "ERROR: type not supported: $type. ASE or ASARP expected\n";
}
print "Analyzing $type results\nSAMPLE1: $sample1\tSAMPLE2: $sample2\n";

print "Powerful cutoff: $powerful; P-value cutoff: $pcutoff\n";

open(my $FPO, ">", $output) or die "ERROR: cannot open output: $output\n";
print $FPO "samples\t$sample1\t$sample2\tcommon\n";
print $FPO "$type\tSAMPLE1_ONLY\tSAMPLE2_ONLY\tSAMPLES12_COMMON\n";

if($type eq "ASE"){
  my $resultType = 'ase.prediction';
  my ($ref1, $snvRef1, $stFlag1) = getAseAll("$sample1.$resultType");
  my ($ref2, $snvRef2, $stFlag2) = getAseAll("$sample2.$resultType");

  if($stFlag1 != $stFlag2){
    print STDERR "Strand-specific settings: SAMPLE_1: $stFlag1; SAMPLE_2: $stFlag2\n";
    die "ERROR: the two samples are not with the same (non) strand-specific settings\n";
  }

  #get their intron/exon categories?
  # no, this is done in intronAseAnlys.pl specifically for introns

  # get intersection
  my ($g1Ref, $gcomRef, $g2Ref) = intersectHashes($ref1, $ref2);
  
  # get intersection of ASE SNVs
  my ($s1Ref, $comRef, $s2Ref) = intersectHashes($snvRef1, $snvRef2);
  
  # brief: print this type to a one liner
  my $no1 = keys %$s1Ref;
  my $no2 = keys %$s2Ref;
  my $noCom = keys %$comRef;
  my $oneLiner = join("\t", ("ASE", $no1, $no2, $noCom));
  print $FPO "$oneLiner\n"; 
  print "SAMPLE1_ONLY: $no1\tSAMPLE2_ONLY: $no2\tSAMPLES12_COMMON: $noCom\n";
    
  # further analyze the Sample_1 specific cases
  #my $snvFile2 = ''; #input by user
  #my $powerful = 20; #input by user
  my ($caseSnv, $powCnt1, $nonCnt1, $unCnt1) = getPowerfulControlSnvs($snvFile1, $powerful, $s1Ref, $stFlag1);
  my ($contrlSnv, $powCnt2, $nonCnt2, $unCnt2) = getPowerfulControlSnvs($snvFile2, $powerful, $s1Ref, $stFlag2);
  print "Case: Powerful: $powCnt1\tNon-powerful: $nonCnt1\n"; #\tAbsent: $unCnt1\n";
  print "Control: Powerful: $powCnt2\tNon-powerful: $nonCnt2\n"; #\tAbsent: $unCnt2\n";
  
 
  my %select = ();
  my $cnt = 0;
  for(keys %$contrlSnv){
    my ($p1, $c1) = split(/\t/, $caseSnv->{$_});
    my ($p2, $c2) = split(/\t/, $contrlSnv->{$_});
    if($p2 <= $pcutoff){
      next;
    }
    $cnt += 1;
    #print "$_\n";
    my ($ctrlP, $ctrlCnt, $ctrlLine) = split(/\t/, $contrlSnv->{$_});
    my (@sar) = split(/ /, $ctrlLine);
    my $ctrlAlCnt = $sar[-1];
    if($stFlag2){
      $ctrlAlCnt = $sar[-2];
    }
    if(!($ctrlAlCnt =~ /\d+:\d+:\d+/)){
      die "ERROR: parsing the allelic count of SAMPLE2: $ctrlAlCnt in $contrlSnv->{$_}\n";
    }
    #
    #$select{$_} = "$caseSnv->{$_}\n$contrlSnv->{$_}\n";
    #$select{$_} = "$caseSnv->{$_}\t(ctrl: $ctrlP;$ctrlCnt;ctrlAlCnt)\n";
    $select{$_} = "$caseSnv->{$_}\t(ctrl: $ctrlP;$ctrlAlCnt)\n";
    #print "$caseSnv->{$_}\n";
    #print "$contrlSnv->{$_}\n";
    #print "\n";
  }


  # gene level results
  my $geneCnt = 0;
  my $content = '';
  for my $gene (keys %$g1Ref){
    my @allSnvs = split(/\t/, $g1Ref->{$gene});
    my $outSnv = '';
    for(@allSnvs){
      if(defined($select{$_})){
        $outSnv .= "$select{$_}";
      }
    }
    if($outSnv ne ''){
      $geneCnt += 1;
      $outSnv = "$gene\n$outSnv";
      #print "$outSnv\n";
      $content .= "$outSnv\n";
    }
  }

  my $geneCntS1 = keys %$g1Ref;
  my $header = '';
  $header .= "#SAMPLE1 Specific $type genes: $geneCnt out of $geneCntS1 $type genes\n";
  $header .= "#SAMPLE1 Specific $type SNVs: $cnt out of $powCnt1 in-$type-gene SNVs ($powCnt2 powerful SAMPLE2 controls)\n";

  print $header;
  $content = $header.$content;

  # SAMPLE1 Specific details
  outputSampleSp($output, "ASE", $content);


  # details: output the overlapping details by type
  outputDetailsByType($output, "ASE", $s1Ref, $s2Ref, $comRef);

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
  
    # brief: print this type to a one liner
    my $no1 = keys %$s1Ref;
    my $no2 = keys %$s2Ref;
    my $noCom = keys %$comRef;
    my $oneLiner = join("\t", ($asType, $no1, $no2, $noCom));
    print $FPO "$oneLiner\n"; 
    
    print "SAMPLE1_ONLY: $no1\tSAMPLE2_ONLY: $no2\tSAMPLES12_COMMON: $noCom\n";
    # details: output the overlapping details by type
    outputDetailsByType($output, $asType, $s1Ref, $s2Ref, $comRef);
  }
    
}
close($FPO);

###################################################################################
#  sub-routines

sub outputSampleSp{

  my ($out, $ty, $content) = @_;
  my $outp = "$out.$ty.sample1sp.txt";
  open(WR, ">", $outp) or die "ERROR: cannot write to $outp\n";
  print WR $content;
  close(WR);
}

sub getPowerfulControlSnvs{

  my ($snvFile, $powerful, $sRef, $flag) = @_;
  my $DILIMIT = ' '; # space separated
  my %ctrl = (); #control SNVs

  my ($pCnt, $nCnt, $uCnt) = (0, 0, 0); # powerful and non-powerful

  # Create a communication bridge with R and start R
  my $R = Statistics::R->new();
  open(FP, $snvFile) or die "ERROR: cannot open SNV file: $snvFile\n";
  while(<FP>){
    chomp;
    my $line = $_;
    my ($chr, $pos, $al, $id, $cnt, $strand) = split(/$DILIMIT/, $line);
    my $snvKey = "$chr;$pos";
    if($flag){
      if(defined($strand) && ($strand eq '+' || $strand eq '-')){
        $snvKey .=";$strand";
      }else{
        die "ERROR: cannot parse strand: $strand in $_\n";
      }
    }

    my $pValue = 1; #not performed

    if(defined($sRef->{$snvKey})){
      my ($rc, $ac, $wc) = split(':', $cnt);
      my $cnt = $rc+$ac;
      if($cnt >= $powerful){ 
        $R->set('x', [$rc, $ac]);
        $R->run('p = chisq.test(x)$p.value'); #default expected dist is c(0.5, 0.5)
        $pValue = $R->get('p');
        #$ctrl{$snvKey} = $cnt; 
        $ctrl{$snvKey} = "$pValue\t$cnt\t$line"; 
        $pCnt += 1; 
      }
      else{ $nCnt += 1; }
      
    }else{
      $uCnt += 1;
    }
  }
  close(FP);
  

  $R->stop;

  return (\%ctrl, $pCnt, $nCnt, $uCnt);

}




