#!/usr/bin/perl
use warnings;
use strict;

# all sub-routine utilities here
# for processing reads and pipe-ups

# mask for certain read positions (with low qualities)
# input:	a position list (one-based start and end, inclusive)
# output:	an array with the masked positions marked 1, and others 0's
# NOTE: the position list is a comma separated list, where each element is in a range form: "a:b"
# indicating read position range from a to b (one-based, inclusive) is discarded
# For the special case of "a:a", one can simply write "a" in the comma list as the element
# For example, 1,4:6,100 is equivalent to 1:1,4:6,100:100
# **IMPORTANT NOTE**: the return mask array is for INTERNAL use, and the positions are converted to **zero-based** ones.
sub readPosMask
{
  my ($discardPos) = @_;
  my @discard = ();

  #print "discard: $discardPos\n";
  print "Evaluating syntax of discarded read positions (one-based, range start-end inclusive)...\n";
  my $previousTo = 0;
  my @poss = split(',', $discardPos);
  for(@poss){
    my ($from, $to, $extra) = split(':', $_);
    if(defined($extra)){
      print "\nERROR: each range only allows one ':': $_\n";
      exit;
    }
    if(!defined($to)){
      $to = $from;
    }
    if(!($from=~/\d+/) || !($to=~/\d+/)){
      print "\nERROR: cannot pass $from or $to as integers\n";
      exit;
    }
    if($previousTo > $from || $from > $to){
      print "\nERROR: discarded read positions/ranges must be sorted: ..:$previousTo,$from:$to\n";
      exit;
    }
    print "$from:$to,";
    for(my $i=$from-1; $i<=$to-1; $i++){ #zero-based internally
      $discard[$i] = 1; #fill-in discarded positions
    }
    $previousTo = $to;
  }
  if($previousTo>0){
    print "\n";
  }
  for(my $i = 0; $i< $previousTo; $i++){
    if(!defined($discard[$i])){	$discard[$i] = 0; }
  }

  return @discard;

}

1;
