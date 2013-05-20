#!/usr/bin/perl
use warnings;
use strict;

# get gene.prediction from ASARP results
# input: 	gene.prediction result file
# output:	$snvsRef: reference to an array containing all SNVs	
#		$sumRef: reference to an array containing all the corresponding SNV read sums
#		$snvHsRef: reference to a hash containing information to distinguish all SNVs
sub getAsarpResult{
 
  my ($input) = @_;

  open(FP, $input) or die "ERROR: cannot open $input to read\n";
  my @pred = <FP>;
  chomp @pred;
  close(FP);

  my $checkText = "ASARP gene level";
  if($pred[0] ne $checkText){
    die "ERROR: expecting $input as gene.prediction (with header: $checkText)\n";
  }

  my %snvHs = (); #hash for all the SNV identifiers
  my @snvs = ();
  my @sums =  ();
  for(my $i = 1; $i < @pred; $i++){
    if($pred[$i] =~ /^chr/){ #starting of a new gene
      my ($chr, $gene) = split('\t', $pred[$i]);
      $i++; #go to the next line

      #get the snps for this gene
      while($i<@pred && $pred[$i] ne "" && !($pred[$i]=~/^chr/)){
        #get SNPs
        my ($dummyInfo, $snpInfo, $strandInfo) = split(';', $pred[$i]);
        my($pos, $id, $al, $reads) = split(' ', $snpInfo);
        my ($r1, $r2) = split(':', $reads);
        my $info = join(" ", $chr, $pos, $al, $id);
	if(!defined($strandInfo) && ($strandInfo eq '+' || $strandInfo eq '-')){
	  $info .= " $strandInfo"; # add strand information
	}

	$dummyInfo .= " $chr $gene $reads"; #chr, gene and reads: X:Y are also needed
        if(!defined($snvHs{$info})){
          $snvHs{$info} = $dummyInfo; #add dummyInfo which will be used for evaluating random results
          push @snvs,$info; 
          $r1 += $r2;
          push @sums,$r1;
        }else{
	  $snvHs{$info} .= "\t".$dummyInfo;
	}
        $i++;
      }
    }
  }

  return (\@snvs, \@sums, \%snvHs);
}

sub getIndex{

  my ($freq) = @_;
  my ($from, $to) = (0, 0);
  
  if($freq =~ /\:/){
    ($from, $to) = split(/\:/, $freq);
  }else{
    $to = $freq-1;
  }
  print "Randomized SNV list index: $from to $to\n";
  return ($from, $to);
}

1;

