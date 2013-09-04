#!/usr/bin/perl

use strict;
use warnings;
require "mirnaParser.pl"; # needs rc

if(@ARGV < 6){

  print <<EOT;

USAGE: perl $0 ASAT_fasta miRanda_result output ext_miRanda_out energy_th percent_th [mirSVR_results]

ASAT_fasta	the merged fasta file containing all fasta sequences with the reference and alternative alleles from the ASAT results
		The fasta identifier information is important and used in the script
		format: >chr;gene;snv_cordinate;snv_id;ref>alt;ref#:alt#;ASAT_info;snv_pos_in_seq;ref[alt]; ref#[alt#]
		Example:
		>chr10;NRP1;33467312;rs10827206;C>T;370:232;AT:3-(0.6122);26;G;370
		GCTCACAAAGAATAAGCCTGCCTTAGGGCTGGCAACATCTAAGCCTCTAAC
		>chr10;NRP1;33467312;rs10827206;C>T;370:232;AT:3-(0.6122);26;A;232;
		GCTCACAAAGAATAAGACTGCCTTAAGGCTGGCAACATCTAAGCCTCTAAC
		
		Note: the two-allele (ref and alt) sequences should go in pairs; the sequence should be the reverse 
		complement if the SNV is on - strand. The last ref[alt] indicates the allele version of the sequence,
		and is the reverse on - strand (i.e. C>T;370:232 converted to G;370 and A;232 respectively)

miRanda_result	the combined miRanda results for all ASAT cases predicted from ASAT_fasta

output		the output result file
ext_miRanda_out	the compact extracted results from miRanda_result, with 2 different suffixes of ext_miRanda_out
		.extlst (extracted list containing miRanda miRNA targets only)
		.complst (comparison list containing alignment scores and energies of all miRanda predictions)

energy_th	the (negated) energy cutoff for allele-specific miRNA target comparisons, 
		e.g. 20 means only miRNA targets with energies <= -20 will be compared
percent_th	the percentage cutoff for average normalized energy changes to be considered and output
		e.g. 0.1 and 1.0 mean 10% and 100% average normalized energy changes respectively

EOT
#		.shrtlst (extracted list containing miRanda miRNA targets only)
  exit;
}

my ($asatin, $in, $out, $miRandaExt, $ENRGTH, $MFETH, $mirSVR) = @ARGV;

my $miRnaOut = $miRandaExt.".extlst";
#my $short = $miRandaExt.".shrtlst";
my $comp = $miRandaExt.".complst";

if($ENRGTH >= 0){
  die "ERROR: minimum free energy (MFE) threshold should be <0 (now $ENRGTH)\n";
}
if($MFETH <= 0){
  die "ERROR: absolute minimum free energy (MFE) change percentage threshold change should be >0, e.g. 0.5 for 50% (now $MFETH)\n";
}

# get all allele-specific version fasta information
my ($ansRef, $k) = getAlleleFasta($asatin);

my %ans = %$ansRef;

# get all miRanda results
my $n = 0;
open(FP, $in) or die "ERROR: cannot onpen miRanda result: $in\n";
my $result = "";
#my $shortlist = "";

#my $miHeader = "Seq1,Seq2,Tot Score,Tot Energy,Max Score,Max Energy,Strand,Len1,Len2,Positions";
my $miHeader = "Seq1\tSeq2\tScore\tEnergy\tSeq1_rg\tSeq2_rg\tmatch_len\tp1\tp2";
my $miHeader2 = "Scores for this hit:";

# get all miRanda_result
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
            #$shortlist .= $_;
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

open(OP, ">", $miRnaOut) or die "ERROR: cannot write to extract file $miRnaOut\n";
print OP $result;
close(OP);

#open(SP, ">", $short) or die "ERROR: cannot write to shortlist file $short\n";
#print SP "$miHeader\n";
#print SP $shortlist;
#close(SP);

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
#my $kdFile = "/home/cyruschan/asarp_perl/U87.DroshaKD/Dro.snvs/U87.DroshaKD.snv.lst";
my $kdFile = "/home/cyruschan/asarp_perl/gm12878/N.A+/snv.ss.S/N.A+.snv.lst";
my ($outInfo) = compareScoreEnergy(\%ans, $kdFile, $ENRGTH, $MFETH);

open(SP, ">", $out) or die "ERROR: cannot output results to $out\n";
print SP $outInfo;
close(SP);

print "FINISHED\n";


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
  my ($ansRef, $kdSnv, $ENRGTH, $MFETH) = @_;
  my %ans = %$ansRef;
  my $outStr = "Y/N\tchr\tgene\t+/-\tsnv\tid\tASAT\tavgMFE\talleles\tmiRNA(maxMFE)\n";
  # compare ref>alt miRanda target prediction scores
      
  # Create a communication bridge with R and start R
  my $R = Statistics::R->new();
  my ($sig, $insig) = (0, 0); #significant and insignificant counts
  my $TH = 1; #0.05; #0.05; #0.05; #0.05; # p-value threshold

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
    if(index($at, 'AT:3-')!= -1 || index($at, 'ASE:3-')!= -1){
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

      # direction check:
      #if(($refCnt-$altCnt)*($kdRefCnt-$kdAltCnt) < 0){
      # p-value check
      if($pValue <= $TH){
        $sig += 1;
	print "$key\n$kdStr\n";
        print "Test ($refCnt, $altCnt), ($kdRefCnt, $kdAltCnt)\np-value: $pValue\n";
      }else{
        $insig += 1;
	#next;
      }
    }else{ print "$key not found in Drosha KD\n\n";  } #next; }
    #else{ next; }
    $scs{"$ref;$refCnt"} = 0; #score
    $scs{"$alt;$altCnt"} = 0; #score
    $ens{"$ref;$refCnt"} = 0; #energy
    $ens{"$alt;$altCnt"} = 0; #energy
    my $details = ''; 
   
    # store the absolute min energy, absolute max score, and the list of hits (no need to distinguish ref and alt)
    my %hits = (); # keys are the two alleles, values are references to miRNA hashes of the two alleles
    for my $al (keys %hs){ # sub-key
      #print "$al\n$hs{$al}";
      #results
      my %mir = ();
      #score
      my $alSc = 0; # score, larger the better
      my $alMirSc = '';
      #energy
      my $alEn = 0; # energy, smaller the better
      my $alMir = '';
      my @vals = split(/\n/, $hs{$al});
      for(@vals){
        my ($miR, $sc, $enrg) = split(/\t/, $_);
        $mir{$miR} = $enrg;
	if($sc > $alSc){ # not used
	  $alSc = $sc;
	  $alMirSc = $miR;
	  # associate the engery with the best alignment score
	  #$alEn = $enrg;
	  #$alMir = $miR;

	}
	if($enrg < $alEn){
	  $alEn = $enrg;
	  $alMir = $miR;
	}
      }
      $scs{$al} = $alSc;
      $ens{$al} = $alEn;
      $hits{$al} = \%mir;
      
      $details .= "$al: max score $alMirSc: $alSc\n";
      #$details .= "$al: aso energy $alMir: $alEn\n";
      $details .= "$al: min energy $alMir: $alEn\n";
    }
    # find the maximal change of MFE (minimum free energy)
    # must distinguish ref and alt (always ref VS alt, otherwise the validation would be invalid)
    my %refhs = ();
    my %alths = ();
    if(defined($hits{"$ref;$refCnt"})){
      %refhs = %{$hits{"$ref;$refCnt"}};
    }
    if(defined($hits{"$alt;$altCnt"})){
      %alths = %{$hits{"$alt;$altCnt"}};
    }
    my $maxMFE = 0;
    my $maxMir = '';
    my ($avgMFE, $avgCnt, $bgMFE, $bgCnt) = (0, 0, 0, 0); #average MEF
    #my ($refmir, $altmir) = ('', '');
    # fill in the highest possilbe minimal energy
    for(keys %refhs){
      if(!defined($alths{$_})){ $alths{$_} = -1; } 
    }
    for(keys %alths){
      if(!defined($refhs{$_})){ $refhs{$_} = -1; }
    }
    for(keys %refhs){
      # # now alths{$_} must be defined

      # normalized diff:
      my $norm = abs($refhs{$_});
      if(abs($alths{$_}) < $norm){ $norm = abs($alths{$_}); } # normalization over the smaller one
      my $diff = ($refhs{$_} - $alths{$_})/$norm; #normalized diff
      # a new filter: energy scores cannot be too low for both alleles
      #print "$_ diff: $diff = ref $refhs{$_} - alt $alths{$_}\n";
      if(abs($diff) > abs($maxMFE)  && ($refhs{$_} <= $ENRGTH || $alths{$_} <= $ENRGTH))
      {
        $maxMFE = $diff; $maxMir = $_;	
      } # only when they are <= threshold
      if(abs($diff) >= $MFETH  && ($refhs{$_} <= $ENRGTH || $alths{$_} <= $ENRGTH))
      {
        $avgMFE += $diff; $avgCnt += 1; 
      }
      # all are considered as background
      if(abs($diff) >= $MFETH)
      {
        $bgMFE += $diff; $bgCnt += 1;
      }
    }
    if($avgCnt){ $avgMFE /= $avgCnt;  }
    if($bgCnt){ $bgMFE /= $bgCnt;  }

    #### debug
    #$avgMFE = $maxMFE;

    if(defined($refhs{$maxMir})){ # it indicates that $alths{$maxMir} is also defined
      print "maxMFE: $maxMir $maxMFE: $ref;$refCnt:$refhs{$maxMir} -- $alt;$altCnt:$alths{$maxMir}\n";
    }else{ print "No maxMFE cases for $ref, $alt\n"; }

=cut
    if($refCnt > $altCnt){ # reference > alt
      if($maxMFE >0){ $enY += 1; print "Yes\n\n"; }
      elsif($maxMFE <0){ $enN += 1; print "No\n\n";}
    }elsif($refCnt < $altCnt){
      if($maxMFE <0){ $enY += 1; print "Yes\n\n";}
      elsif($maxMFE >0){ $enN += 1; print "No\n\n";}
    }
=cut

    my $case = join("\t", $chr, $gene, $std, $pos, $id, $at, $avgMFE, "$ref;$refCnt>$alt;$altCnt", $maxMir)."\n";
    if($refCnt > $altCnt){ # reference > alt
      if($avgMFE >=$MFETH){ 
        print $case;
	$outStr .= "Y\t$case";
        $enY += 1; #print "Yes\n\n"; 
      }
      elsif($avgMFE <= -$MFETH){ 
	$outStr .= "N\t$case";
        $enN += 1; 
      }# print "No\n\n";}
    }elsif($refCnt < $altCnt){
      if($avgMFE <= -$MFETH){ 
        #print "avgMFE: $avgMFE: $ref;$refCnt -- $alt;$altCnt\n";
        print $case;
	$outStr .= "Y\t$case";
        $enY += 1; #print "Yes\n\n";
      }
      elsif($avgMFE >= $MFETH){
	$outStr .= "N\t$case";
        $enN += 1; 
      }# print "No\n\n";}
    }
    next;

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
	#print "$key\n$kd\n$details\n";

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
	#print "$key\n$kd\n$details\n";
      }elsif($scs{"$ref;$refCnt"} < $scs{"$alt;$altCnt"} && $ens{"$ref;$refCnt"} > $ens{"$alt;$altCnt"}){ # N: larger count larger score, smaller energy
        $seN += 1;
      }
    }
    print "$key\n$kd\n$details\n";

  }
  $R->stop;
  print "Energy (MFE) threshold: -$ENRGTH; Change percentage threshold: $MFETH\n"; 
  print "Significant $TH: $sig\nInsignificant $TH: $insig\n";
  #print "Y (higher score, fewer reads): $scY\nN (higher score, more reads): $scN\n";
  print "Y (lower energy, fewer reads): $enY\nN (lower energy, more reads): $enN\n";
  #print "Y (hi sc  lo en, fewer reads): $seY\nN (hi sc  lo en, more reads): $seN\n";
  #
  $outStr .= "\nY (lower energy, fewer reads): $enY\nN (lower energy, more reads): $enN\n";
  return $outStr;
}

