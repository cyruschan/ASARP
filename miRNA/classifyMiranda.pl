#!/usr/bin/perl
use strict;
use warnings;

if(@ARGV < 4 ){
  print <<EOT;

USAGE: perl $0 input_summary score_th energy_th output_summary

EOT
  exit;
}
my ($input, $SCORE, $ENERGY, $output) = @ARGV;

my $res = readSummary($input, $SCORE, $ENERGY);
outputResults($output, $res);


sub readSummary
{
  my ($in, $SCORETH, $ENERGYTH) = @_;
  my %hits = ();
  my $content = '';
  open(FP, $in) or die "ERROR: cannot open $in\n";
  print "Reading miRNA target info from $in\n";
  my ($gCnt, $sCnt, $tCnt) = (0, 0, 0);
  while(<FP>){
    chomp;
    my ($mir, $snv, $sc, $en, $mirRg, $mrnaRg, $len, $mirP, $mrnaP) = split(/\t/, $_);
    my ($isRef, $gKey, $sKey, $off, $al, $alCnt) = parseSnvInfo($snv);
    my ($start, $end) = split(' ', $mrnaRg);
    if($off < $start || $off > $end){ next; } # skip out-of-range cases
    my $val = "$isRef;$mir;$sc;$en"; #isRef, mir, score and engergy
    if(!defined($hits{$gKey})){ # a new case
      my %ss = ();
      $ss{$sKey} = $val;
      $hits{$gKey} = \%ss; # hash reference
      $gCnt += 1;
      $sCnt += 1; # snv count!
    }else{
      my $sRef = $hits{$gKey};
      if(!defined($sRef->{$sKey})){
        $sRef->{$sKey} = $val;
	$sCnt += 1;
      }else{
        $sRef->{$sKey} .= "\t$val";
      }
    }
    $tCnt += 1;
    #print "Add $gKey, $sKey with value: $val\n";
    #print "$hits{$gKey}->{$sKey}\n";
    #exit;
  }
  close(FP);
  print "$gCnt genes with $sCnt SNVs and $tCnt initial mirna targets\n";
  print "Filtering allele-specific miRNA targets\n";

  my ($dsrpt, $crt, $swtch, $shrd, $miTargets) = (0, 0, 0, 0, 0);
  my ($refYes, $altYes) = (0, 0);
  for my $gk (keys %hits){
    print "$gk\n";
    for my $sk (keys %{$hits{$gk}}){
      print "$sk\n"; # for one SNV
      #"$pos;$id;$all;$cnt"
      my ($spos, $sid, $sall, $scnt) = split(';', $sk);
      my ($refNo, $altNo) = split (':', $scnt);

      my ($miDs, $miCr, $miSh) = (0, 0, 0); # micro RNA disrupt, creation, shared
      my @sHits = split(/\t/, $hits{$gk}->{$sk});
      my %refHs = ();
      my %altHs = ();
      my ($minRefEn, $minAltEn) = (0, 0);
      for(@sHits){ # for all miRNA of one SNV
        my ($isRef, $mir, $sc, $en) = split(';', $_);
	if($isRef){
	  $refHs{$mir} = "$sc;$en";
	  if($en < $minRefEn){
	    $minRefEn = $en;
	  }
	}else{
	  $altHs{$mir} = "$sc;$en";
	  if($en < $minAltEn){
	    $minAltEn = $en;
	  }
	}
      }
      # clear all non-allele-specific targets
      # remove shared
      my ($refRef, $altRef, $tShrd, $minShrdEn) = filterShared(\%refHs, \%altHs);
      $miSh += $tShrd;

      $minShrdEn = 0; # ignore minimal shared energy
      $refRef = filterHits($refRef, $ENERGYTH, $SCORETH, $minShrdEn);
      $altRef = filterHits($altRef, $ENERGYTH, $SCORETH, $minShrdEn);
      %refHs = %$refRef;
      %altHs = %$altRef;

      #print "SHARED\n";
      #print "DISRUPT (ref only)\n";
      $miDs += keys %refHs;
      for(keys %refHs){
        my $msg = "DISRUPT:REF $_\t$refHs{$_}\n";
	$content .=$msg;
	print $msg;
      }
      #print "CREATE (alt only)\n";
      $miCr += keys %altHs;
      for(keys %altHs){
        my $msg = "CREATE:ALT $_\t$altHs{$_}\n";
	$content .= $msg;
      }
      if($miDs && $miCr){ # both disrupt and create at the same SNV: that's switch
        $swtch += 1;
      }else{
        if($miDs){	$dsrpt += 1;	}
        if($miCr){	$crt += 1;	}
	if(!$miCr && !$miDs){
	  $shrd += 0;
	}
	# REF should have fewer reads
	if($miDs && $refNo < $altNo){ #disrupt
	  #print "Yes!!\n";
	  $refYes += 1;
	}
	# ALT should have more reads
	if($miCr && $refNo > $altNo){ # creation
	  #print "Yes!!\n";
	  $altYes += 1;
	}
      }
      $miTargets += $miSh+$miCr+$miDs; # total mir Targets

    } # end of one SNV
  }
  my $ovrMsg = "SNV: disrupted: $dsrpt ($refYes)\tcreated: $crt ($altYes)\tswitched: $swtch\tshared: $shrd\n";
  $ovrMsg .= "gene: $gCnt\tsnv: $sCnt\tmiRNA: $miTargets\nTarget ratios: gene: ".($miTargets/$gCnt)." snv: ".($miTargets/$sCnt)."\n";
  print $ovrMsg;
  $content .= $ovrMsg;
  return $content;
}

sub outputResults{
  my ($out, $cont) = @_;
  open(WR, ">", $out) or die "ERROR: cannot output to $out\n";
  print WR $cont;
  close(WR);

}

sub filterShared{

   my ($refHs, $altHs) = @_;
   my $shrd = 0;
   my $minEn = 0; # minimal energy in the shared targets

      for(keys %$refHs){
        if(defined($altHs->{$_})){
	  print "SHARED: $_\t$refHs->{$_}\t$altHs->{$_}\n";
	  $shrd += 1;
	  my ($scr, $enr) = split(';', $refHs->{$_});
	  my ($sca, $ena) = split(';', $altHs->{$_});
	  if(abs($ena) > abs($enr)){ $enr = $ena; }
	  if(abs($enr) > abs($minEn)){
	    $minEn = $enr; # the lowest energy (with maximal abs)
	  }

	  delete $refHs->{$_};
	  delete $altHs->{$_};
	}
      }
   print "MIN_SHARED_ENERGY: $minEn\n";
   return ($refHs, $altHs, $shrd, $minEn);
}      

sub parseSnvInfo{

  my ($snv) = @_;
  my ($chr, $gene, $pos, $id, $all, $cnt, $type, $off, $al, $alCnt) = split(';', $snv);
  my ($ref, $alt) = split('>', $all);
  my ($rcnt,$acnt) = split(':', $cnt);
  # check the strand and determine whether it is ref or alt
  my $actal = $al; # actual allele
  if($type =~ /ASE:3([\+|-])/){
    my $strand = $1;
    #print "strand: $strand\n";
    $actal = rc($al);
    #exit;
  }else{
    die "ERROR: not ASE:3+/-\n";
  }
  my $isRef = 0;
  if($actal eq $ref){
    $isRef = 1;
  }
  return ($isRef, "$chr;$gene", "$pos;$id;$all;$cnt", $off, $al, $alCnt);
}

sub rc
{
  my ($str) = @_;
  $str = uc $str;
  # Make a new copy of the DNA (see why we saved the original?)
  my $revcom = reverse $str;
  #
  # # See the text for a discussion of tr///
  $revcom =~ tr/ACGT/TGCA/;
  return $revcom;
}


sub filterHits
{
  my ($refHs, $ENERGY, $SCORE, $SHRDEN) = @_;

      for(keys %$refHs){
        my ($temps, $tempe) = split(';', $refHs->{$_});
	if(abs($tempe) < abs($ENERGY) || $temps < $SCORE || abs($tempe) < abs($SHRDEN)){
	  #print "REMOVE $_ $refHs->{$_}\n";
	  delete $refHs->{$_}; # not used
	}
      }
      return $refHs;
}
