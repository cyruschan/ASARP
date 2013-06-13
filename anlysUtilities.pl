#!/usr/bin/perl
use warnings;
use strict;

require "fileParser.pl"; #sub's for input annotation files
require "snpParser.pl"; #sub's for snps
use MyConstants qw( $CHRNUM $asarpTags );
our @asarpChecks = split(' ', $asarpTags);

################################################################################3
#####  all sub-routines in compSamples.pl ###########


# this sub-routine is slightly different than getAsarp in random/randUtilities.pl
# for a different purpose
sub getAsarpAll{

  my ($asarp) = @_;
 
  #get all predictions for processing
  open(FP, $asarp) or die "cannot open $asarp\n";
  my @pred = <FP>;
  chomp @pred;
  close(FP);
  my $checkText = "ASARP gene level";
  if($pred[0] ne $checkText){
    die "ERROR: expecting $asarp as gene.prediction (with header: $checkText)\n";
  }
  
  my %results = ();
  for(@asarpChecks){
    $results{$_} = getAsarpByType(undef, \@pred, $_);
  }
  return \%results;
}

sub getAsarpByType{

  my ($asasRef, $predRef, $selectType) = @_;
  my @pred = @$predRef;
  my %asas = ();
  if(defined($asasRef)){ %asas = %$asasRef; }
  # type can be AS, AI, AT, ASS, SE RI
  $selectType = checkValidAsarpType($selectType);


  for(my $i = 1; $i < @pred; $i++){
    if($pred[$i] =~ /^chr/){ #starting of a new gene
      my ($chr, $gene) = split('\t', $pred[$i]);
      $i++; #go to the next line

      #get the snps for this gene
      while($i<@pred && $pred[$i] ne "" && !($pred[$i]=~/^chr/)){
        #get SNPs
        my ($dummyInfo, $snpInfo, $strandInfo) = split(';', $pred[$i]);
        if(index($dummyInfo, $selectType.':') != -1 || (index($dummyInfo, 'AS:') != -1 && index($dummyInfo, '^'.$selectType.'(') != -1)){ #1st: AI/AS/AT; 2nd: RI/SE/ASS
	  #get snpInfo    #my($pos, $id, $al, $reads) = split(' ', $snpInfo);
	  my $key = "$chr;$gene"; #gene is already strand specific
          if(defined($asas{$key})){
	    $asas{$key} .="\t$pred[$i]";
	  }else{
	    $asas{$key} ="$pred[$i]";
	  }
	}
        $i++;
      }
    }
  }
  return \%asas;
}

sub showAsarpByType
{
  my ($asasRef, $selectType) = @_;
  my %asas = ();
  if(defined($asasRef)){ %asas = %$asasRef; }
  # type can be AS, AI, AT, ASS, SE RI
  $selectType = checkValidAsarpType($selectType);
  if(defined($asas{$selectType})){
    print "$selectType gene-level results\n";
    my %hs = %{$asas{$selectType}}; 
    for(keys %hs){
      print "$_\n";
      my @instances = split(/\t/, $hs{$_});
      for(@instances){
        print "$_\n";
      }
    }
    print "\n";
  }
}


sub outputDetailsByType{
  my ($output, $type, $s1ExRef, $s2ExRef, $sComRef) = @_;

  open(SFP, ">", "$output.$type") or die "ERROR: cannot open output: $output.$type\n";
  print SFP "SAMPLE1_ONLY_GENES:$type\n";
  print SFP formatDetailsCore($s1ExRef);

  print SFP "SAMPLE2_ONLY_GENES:$type\n";
  print SFP formatDetailsCore($s2ExRef);
  
  print SFP "SAMPLES12_COMMON_GENES:$type\n";
  print SFP formatDetailsCore($sComRef);
  
  close(SFP);
}

# want to output the cases in the correct chromosome order
sub formatDetailsCore
{
  my ($sRef) = @_;
  my $toPrint = "";
  my %sHs = %$sRef;
  
  my @chrs = (); 
  for(my $i=0; $i<=$CHRNUM; $i++){
    $chrs[$i] = "";
  }
  my @keys = keys %sHs;
  # get all chr indices
  for(@keys){
    my ($chr, $gene) = split(';', $_);
    my $chrNo = getChrID($chr);
    $chrs[$chrNo] .= "$_\t";
  }
  
  for(my $i=1; $i<=$CHRNUM; $i++){
    if($chrs[$i] ne ""){
      my @allKeys = split(/\t/, $chrs[$i]);
      for(@allKeys){
        my ($chr, $gene) = split (';', $_);
	$toPrint .= "$chr\t$gene\n$sHs{$_}\n";
      }
    }
  }
  return $toPrint;
}

# this sub-routine is modified from the one in random/randUtilities.pl which is used in random/evalAsePValues.pl
# get ase.prediction from ASARP results

# the difference is that, in this case, all the reported SNVs are also desired for more detailed comparisons

# get ase.prediction from ASARP results
# input: 	gene.prediction result file
# output:	$aseHsRef: 
#		reference to a hash containing information to 
#		distinguish all ASE genes and ASE gene counts
sub getAseAll{
 
  my ($input) = @_;

  open(FP, $input) or die "ERROR: cannot open $input to read\n";
  my @pred = <FP>;
  chomp @pred;
  close(FP);

  my $checkText = "ASE gene level (all powerful SNVs are ASEs)";
  if($pred[0] ne $checkText){
    die "ERROR: expecting $input to be ase.prediction (with header: $checkText)\n";
  }

  my %snvHs = (); #hash for all the ASE SNV identifiers
  my %aseHs = (); #hash for all the ASE gene identifiers
  for(my $i = 1; $i < @pred; $i++){
    if($pred[$i] =~ /^chr/){ #starting of a new gene
      my ($chr, $gene) = split('\t', $pred[$i]);
      my $key = "$chr;$gene";
      if(!defined($aseHs{$key})){
        $aseHs{$key} = 1;
      }else{
        $aseHs{$key} += 1; #redundant here, but necessary in counting random cases
      }
      $i++; #go to next line!
      my @allKeySnvs = (); # new: store all SNVs of that ASE gene
      # need to get the SNVs as well (snv-level FDR is needed as well)
      while($i<@pred && $pred[$i] ne "" && !($pred[$i]=~/^chr/)){
        #get SNPs
        #sample: rs324419,1.620014e-06,T>C,46871986,0:23;+
        my ($snpInfo, $strandInfo) = split(';', $pred[$i]);
        my($id, $p, $al, $pos, $reads) = split(',', $snpInfo);
        my ($r1, $r2) = split(':', $reads);
        my $info = join(" ", $chr, $pos, $al, $id);
        my $keySnv = "$chr;$pos"; #chr, gene and reads: X:Y are also needed
        
	if(defined($strandInfo) && ($strandInfo eq '+' || $strandInfo eq '-')){
          $info .= " $strandInfo"; # add strand information
	  $keySnv = "$keySnv;$strandInfo";
        }
	push @allKeySnvs, $keySnv; # a new snv added
        
	#print "Adding $keySnv with $info\n";     
        $snvHs{$keySnv} = $pred[$i]; #information should be added only once
        $i++;
      }

      # update the ASE gene entry (different from randomization)
      $aseHs{$key} = join("\t", @allKeySnvs);
    }
  }

  return (\%aseHs, \%snvHs);
}

sub intersectHashes
{
  my ($ref1, $ref2) = @_;

  #simply do the comparison
  my %s1Hs = %$ref1;
  my %s2Hs = %$ref2;
  
  my %s1Ex = (); #sample 1 exclusive
  my %s2Ex = (); #sample 2 exclusive
  my %sCom = (); #common

  for(keys %s1Hs){
      if(defined($s2Hs{$_})){
        if(!defined($sCom{$_})){
	  $sCom{$_} = $s1Hs{$_}."\n".$s2Hs{$_}; 
	}
      }else{
        $s1Ex{$_} = $s1Hs{$_}; # exclusive
      }
  }
  for(keys %s2Hs){
      if(defined($s1Hs{$_})){ #should be done
      }else{
        $s2Ex{$_} = $s2Hs{$_}; # exclusive
      }
  }
  return (\%s1Ex, \%sCom, \%s2Ex);
}
