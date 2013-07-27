#!/usr/bin/perl

use strict;
use warnings;
require "mirnaParser.pl"; # needs rc

if(@ARGV < 4){

  print <<EOT;

USAGE: perl $0 miRanda_result extracted_output short_list comp_list [mirSVR_results]

EOT
  exit;
}

my ($in, $out, $short, $comp, $mirSVR) = @ARGV;

open(FP, $in) or die "ERROR: cannot onpen miRanda result: $in\n";
my $result = "";
my $shortlist = "";

#my $miHeader = "Seq1,Seq2,Tot Score,Tot Energy,Max Score,Max Energy,Strand,Len1,Len2,Positions";
my $miHeader = "Seq1\tSeq2\tScore\tEnergy\tSeq1_rg\tSeq2_rg\tmatch_len\tp1\tp2";
my $miHeader2 = "Scores for this hit:";

my $n = 0;
my %ans = ();

while(<FP>){
  if(index($_, "Energy:") != -1){
    my $hit = 0;
    my $isScore = 0;
    my $temp = "$_";
    while(<FP>){
      if($isScore){
        # check if the result overlaps the snv position
	my ($miS, $mS, $Sc, $enrg, $miRg, $mRg, $aliLen) = split(/\t/, $_);
	my $miRNA = substr($miS, 1); # remove '>'
        my ($mFrom, $mTo) = split(' ', $mRg);
	my ($chr, $gene, $pos, $id, $als, $cnts, $at, $inPos, $al, $cnt) = split(';', $mS);

        my $flag = filterWithMirSVR($chr, $gene, $pos, $miRNA, $mirSVR);
	if($flag){
	  #print " ($miS, $mS, $Sc, $enrg, $miRg, $mRg: $mFrom, $mTo, $aliLen)\n ($chr, $gene, $pos, $id, $als, $cnts, $at, $inPos, $al, $cnt)\n";
	  if($inPos <= $mTo && $inPos >= $mFrom){
            $shortlist .= $_;
	    $n += 1; 
	    my $key = "$chr;$gene;$pos;$id;$als;$cnts;$at";
	    my $subkey = "$al;$cnt";
	    my $val = "$miRNA\t$Sc\t$enrg"; #\t$miRg\t$mRg\t$aliLen
	    if(!defined($ans{$key})){
	      my %hs = ($subkey =>"$val\n"); 
	      $ans{$key} = \%hs;
	    }else{
	      if(!defined($ans{$key}->{$subkey})){
	        $ans{$key}->{$subkey} = "$val\n";
	      }else{
	        $ans{$key}->{$subkey} .= "$val\n";
	      }
	    }
	    #print "$_\n";
	  }
	}
	$isScore = 0;
      }
      if($_ eq "Complete\n"){
        $temp .= "$_\n";
        last;
      }
      $temp .= $_; # keep concatenating
      if($_ eq "$miHeader2\n"){
        $hit = 1; # there is a hit inside
	$isScore = 1;
      }
    }
    if($hit){
      $result .= $temp;
    }
  }
}
close(FP);
print "In total $n targets extracted\n";

open(OP, ">", $out) or die "ERROR: cannot write to extract file $out\n";
print OP $result;
close(OP);

open(SP, ">", $short) or die "ERROR: cannot write to shortlist file $short\n";
print SP "$miHeader\n";
print SP $shortlist;
close(SP);

open(CP, ">", $comp) or die "ERROR: cannot write to comparison file $comp\n";
# compare ref>alt miRanda target prediction scores
for my $key (keys %ans){
  print CP "$key\n";
  my %hs = %{$ans{$key}};
  my $m = keys %hs;
  for my $al (keys %hs){
    print CP "$al\n$hs{$al}";
  }
  print CP "\n";
}
close(CP);

#compareEnergy(\%ans);
my $kdFile = "/home/cyruschan/asarp_perl/U87.DroshaKD/Dro.snvs/U87.DroshaKD.snv.lst";
compareScoreEnergy(\%ans, $kdFile);



##################################################################################
###########		sub-routines		########
#
#
#
sub filterWithMirSVR{
  my ($chr, $gene, $pos, $miRNA, $mirSVR) = @_;
  if(!defined($mirSVR)){
    return 1;
  }
  my $check = open(MFP, $mirSVR);
  if(!$check){
    print STDERR "$mirSVR not readable! ignore filtering\n";
    return 1;
  }
  close(MFP);

  my $flag = 0;
  # command:
  #my $cmd = "grep -w $gene $mirSVR | cut -f 2,4,7,14,15,16,17,18,19";
  my $cmd = "grep -w $gene $mirSVR | cut -f 2,4,14 | grep '$miRNA'"; # $mir, $gene, $loc
  my $res = qx($cmd);
  my @selected = split(/\n/, $res);
  for my $mtch (@selected){
    my ($mir, $miGene, $loc) = split(/\t/, $mtch);
    #[hg19:3:38524951-38524971:+]
    if($loc =~ /\[hg19:(X|Y|M|\d+):([,\-\d]+):[\+|\-]\]/){
      my ($miChr, $rg) = ($1, $2);
      my @range = split(',', $rg);
      for(@range){
        my ($s, $e) = split('-', $_);
        #print "$mtch\n$miChr, $s, $e\n";
        #if($mir eq $miRNA && $pos >= $s && $pos <= $e){
        if($pos >= $s && $pos <= $e){ # do not require the same miRNA and mir
          print "MATCH!!! $mir $pos\n";
          $flag = 1; last; # no need to search other 
	}
      }
    }else{ die "ERROR: cannot parse $mir $miGene: $loc\n"; }
  }

  return $flag;
}

sub compareScoreEnergy{
  my ($ansRef, $kdSnv) = @_;
  my %ans = %$ansRef;
  # compare ref>alt miRanda target prediction scores
      
  # Create a communication bridge with R and start R
  my $R = Statistics::R->new();
  my ($sig, $insig) = (0, 0); #significant and insignificant counts
  my $TH = 0.05; # p-value threshold

  my ($scY, $scN, $enY, $enN, $seY, $seN) = (0, 0, 0, 0, 0, 0);
  for my $key (keys %ans){
    my %hs = %{$ans{$key}};
    my $m = keys %hs; # should be 1 or 2
    
    my %ens = ();
    my %scs = ();
    #print "$key\n";
    my ($chr, $gene, $pos, $id, $als, $cnts, $at) = split(';', $key);
    my $std = '+';
    my $kd = qx(grep $pos $kdSnv);
    #print "$kd\n";
    my ($ref, $alt) = split('>', $als);
    my ($refCnt, $altCnt) = split(':', $cnts);
    if(index($at, 'AT:3-')!= -1){
      $ref = rc($ref);
      $alt = rc($alt); #reverse complement
      $std = '-';
    }
    my @kds = split(/\n/, $kd);
    my $kdStr = '';
    for(@kds){
      if($_ ne '' && substr($_, -1, 1) eq $std){ # check string
        $kdStr = $_; last;
      }
    }
    if($kdStr ne ''){
      # perform parsing and chi-squared test
      my ($kdChr, $kdPos, $kdAls, $kdId, $kdCnts, $kdStd) = split(' ', $kdStr);
      my ($kdRefCnt, $kdAltCnt) = split(':', $kdCnts);
      # chi-squared
      $R->run("M <- as.table(rbind(c($refCnt, $altCnt), c($kdRefCnt, $kdAltCnt)))");
      $R->run('p = fisher.test(M)$p.value'); #default expected dist is c(0.5, 0.5)
      my $pValue = $R->get('p');
      if($pValue <= $TH){
        $sig += 1;
	print "$key\n$kdStr\n";
        print "M <- as.table(rbind(c($refCnt, $altCnt), c($kdRefCnt, $kdAltCnt)))\np-value: $pValue\n\n";
      }else{
        $insig += 1;
	next;
      }
    }#else{ print "$key not found in Drosha KD\n"; }
    else{ next; }
    $scs{"$ref;$refCnt"} = 0; #score
    $scs{"$alt;$altCnt"} = 0; #score
    $ens{"$ref;$refCnt"} = 0; #energy
    $ens{"$alt;$altCnt"} = 0; #energy
    my $details = ''; 
    #my $key = "$chr;$gene;$pos;$id;$als;$cnts;$at";
    #my $subkey = "$al;$cnt";
    #my $val = "$miRNA\t$Sc\t$enrg"; #\t$miRg\t$mRg\t$aliLen
    for my $al (keys %hs){ # sub-key
      #print "$al\n$hs{$al}";
      #score
      my $alSc = 0; # score, larger the better
      my $alMirSc = '';
      #energy
      my $alEn = 0; # energy, smaller the better
      my $alMir = '';
      my @vals = split(/\n/, $hs{$al});
      for(@vals){
        my ($miR, $sc, $enrg) = split(/\t/, $_);
	if($sc > $alSc){
	  $alSc = $sc;
	  $alMirSc = $miR;
	}
	if($enrg < $alEn){
	  $alEn = $enrg;
	  $alMir = $miR;
	}
      }
      $scs{$al} = $alSc;
      $ens{$al} = $alEn;
      
      $details .= "$al: max score $alMirSc: $alSc\n";
      $details .= "$al: min energy $alMir: $alEn\n";
    }
    # comparison:
    if($refCnt > $altCnt){ # reference > alt
      # compare score
      if($scs{"$ref;$refCnt"} < $scs{"$alt;$altCnt"}){ # Y: larger count smaller score
        $scY += 1;
      }elsif($scs{"$ref;$refCnt"} > $scs{"$alt;$altCnt"}){ # N: larger count larger score
        $scN += 1;
      }
      # compare energy
      if($ens{"$ref;$refCnt"} > $ens{"$alt;$altCnt"}){ # Y: larger count higher energy
        $enY += 1;
      }elsif($ens{"$ref;$refCnt"} < $ens{"$alt;$altCnt"}){ # N: larger count lower energy
        $enN += 1;
      }
      # compare score and energy
      if($scs{"$ref;$refCnt"} < $scs{"$alt;$altCnt"} && $ens{"$ref;$refCnt"} > $ens{"$alt;$altCnt"}){ # Y: larger count smaller score, higher energy
        $seY += 1;
	print "$key\n$kd\n$details\n";

      }elsif($scs{"$ref;$refCnt"} > $scs{"$alt;$altCnt"} && $ens{"$ref;$refCnt"} < $ens{"$alt;$altCnt"}){ # N: larger count larger score, smaller energy
        $seN += 1;
      }
    }elsif($refCnt < $altCnt){
      # compare score
      if($scs{"$ref;$refCnt"} > $scs{"$alt;$altCnt"}){ # Y: larger count smaller score
        $scY += 1;
      }elsif($scs{"$ref;$refCnt"} < $scs{"$alt;$altCnt"}){ # N: larger count larger score
        $scN += 1;
      }
      # compare energy
      if($ens{"$ref;$refCnt"} < $ens{"$alt;$altCnt"}){ # Y: larger count higher energy
        $enY += 1;
      }elsif($ens{"$ref;$refCnt"} > $ens{"$alt;$altCnt"}){ # N: larger count lower energy
        $enN += 1;
      }
      # compare score and energy
      if($scs{"$ref;$refCnt"} > $scs{"$alt;$altCnt"} && $ens{"$ref;$refCnt"} < $ens{"$alt;$altCnt"}){ # Y: larger count smaller score, higher energy
        $seY += 1;
	print "$key\n$kd\n$details\n";
      }elsif($scs{"$ref;$refCnt"} < $scs{"$alt;$altCnt"} && $ens{"$ref;$refCnt"} > $ens{"$alt;$altCnt"}){ # N: larger count larger score, smaller energy
        $seN += 1;
      }
    }

  }
  $R->stop;
  print "Significant $TH: $sig\nInsignificant $TH: $insig\n";
  print "Y (higher score, fewer reads): $scY\nN (higher score, more reads): $scN\n";
  print "Y (lower energy, fewer reads): $enY\nN (lower energy, more reads): $enN\n";
  print "Y (hi sc  lo en, fewer reads): $seY\nN (hi sc  lo en, more reads): $seN\n";
}

