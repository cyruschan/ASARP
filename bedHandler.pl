#!/usr/bin/perl/ -w
use strict;
#use IO::Handle;

use MyConstants qw( $CHRNUM $TENKB $supportedList );

require "fileParser.pl"; # need to use binarySearch

## This file handles the bed files parsing, reading and data retrieval.
## To manage a small memory footprint, a binned approach is used.
## Therefore, it is strongly recommended to use the getter subroutines
## e.g.getReadSlice, getEffReadSumLength 
## to retrieve certain read counts and their sum of interest


##################################################################################
### the following sub's are about bed reading and processing

# sub routine to read in the bed files
#	input:	bed file name for certain chromosome, the path to human genome (hg19)
#	output: the reference to the bed arrays (index: genomic coordinates, elements: read counts)
sub readBedByChr
{
  #my ($bedFiles, $genomePath) = @_;
  #my @files = glob($bedFiles);  
  #for my $file (@files){
  #  open(my $fh, "<", $file) or die "Cannot open bed file: $file\n";
  #  close($fh);
  #}
  my ($bedFolder, $genomePath, $chr) = @_;
  if(!($genomePath =~ /[\\|\/]$/)){
    $genomePath = $genomePath."/";  }
  my $chrFile = $genomePath."chr".$chr.".fa";

  open(my $cFh, "<", $chrFile) or die "Cannot open chromosome file: $chrFile\n";
  #get all the data without newlines
  my @seqs = <$cFh>;
  close($cFh);
  my $header = $seqs[0]; #\n included
  my $noNewlines = @seqs -1; #first newline in $header
  if(!$seqs[-1] =~ /\n$/){ $noNewlines-=1;  }
  my $fileSize = -s $chrFile; #get file size
  my $charSize = $fileSize - length($header) - $noNewlines; 
  
  my @beds = (); #use a binned function to save memory
  
  my @reads = ();
  my @reads_idx = ();

  if(!($bedFolder =~ /[\\|\/]$/)){
    $bedFolder .= "/";
  }
  my $bedFile = $bedFolder."chr".$chr.".sam.bed";

  open(my $bFh, "<", $bedFile) or die "Cannot open bed file: $bedFile\n";
  my $dummy = <$bFh>; chomp $dummy;
  print "Reading bed file: $bedFile...\t";
  #print "Chromosome size: $charSize\n";
  if(!$dummy =~ /chr$chr/){  die "Inconsistent chr $chr with $dummy in $bedFile\n"; }
  $dummy = <$bFh>; #get rid of the first 2 header lines
  my $count = 0;
  while(<$bFh>){
    chomp;
    $count++;
    #if($count%($TENKB*100)==0){ print $count." ";  }
    my ($ch, $start, $end, $n) = split(' ', $_);
    #$start; #already zero-based?
    $end-=1;#1-based, converted back to zero based
    my $sBin = int($start/$TENKB);
    my $eBin = int($end/$TENKB);

    my $sOff = $start-$sBin*$TENKB;
    my $eOff = $end-$eBin*$TENKB;
    if($sOff >= $TENKB || $eOff >= $TENKB){
      die "Error in calculating bin pos $sOff-$eOff for $start-$end\n";
    }

    for(my $i=$sBin; $i<=$eBin; $i++){
      my $startIn = $sOff;
      if($i>$sBin){ $startIn = 0; }
      my $endIn = $eOff;
      if($i<$eBin){ $endIn = $TENKB-1; }
      $beds[$i] .= "$startIn,$endIn,$n\t";
    }
  }
  close($bFh);
  
  my $binSize = @beds;
  #print "Bins: $binSize\n";

#  print "Putting counts in positions\n";
#  for(my $i=0; $i<$binSize; $i++){
#    if(defined($beds[$i])){
#      my @posArray = (0)x$TENKB;
#      my @pos = split('\t', $beds[$i]);
#      foreach(@pos){
#        my ($sOff, $eOff, $n) = split(',', $_);
#	for($sOff..$eOff){
#	  if($posArray[$_] >0 && $n>0){ die "Location $_ already filled in ".$posArray[$_]."\n"; }
#	  $posArray[$_]=$n; #they won't overlap, otherwise there is error
#	}
#      }
#      $beds[$i]=''; #clear all first
#      for(0..$TENKB-1){
#        $beds[$i] = $beds[$i].$posArray[$_].",";
#      }
#    }
#  }

  for(my $i=0; $i<$binSize; $i++){
    if(defined($beds[$i])){
      my %posHash = ();
      my @pos = split('\t', $beds[$i]);
      foreach(@pos){
        my ($sOff, $eOff, $n) = split(',', $_);
	if(defined($posHash{$sOff})){
	  die "ERROR: no beds are assumed to have the same starting position:\n$posHash{$sOff} and $_\n";
	}
	$posHash{$sOff} = "$eOff,$n";
      }
      $beds[$i] = ''; #clear all first to save memory

      $reads[$i] = \%posHash;
      $reads_idx[$i] = [sort {$a<=>$b} keys %posHash];
    }
  }
  my %allReads = (
    'reads'=>\@reads,
    'reads_idx' => \@reads_idx,
  );

  print "done\n";
  return \%allReads;
}


# sub routine to get the counts for a certain slice (sub-array) from the beds array
#the input positions ($from, $to) are assumed to be 1-based (incl)
sub getReadSlice
{
  my ($bedRef, $from, $to) = @_;
  my @beds = @{$bedRef->{'reads'}};
  my @bedsIdx = @{$bedRef->{'reads_idx'}};
  
  #print "from: $from to $to\n";
  if($from > $to){ die "ERROR: wrong slice range: [$from, $to]\n"; }
  if($to - $from + 1 > $TENKB){
    die "ERROR: read count slice [$from, $to] longer than $TENKB is not yet supported. You need to modify the source to implement that\n";
  }
  $from -= 1; $to -= 1; #the input positions are assumed to be 1-based (incl), now internally converted to 0-based
  my $sBin = int($from/$TENKB);
  my $eBin = int($to/$TENKB);

  my $sOff = $from-$sBin*$TENKB;
  my $eOff = $to-$eBin*$TENKB;
  my @readArray = ();

  #print "bins: $sBin, $eBin, offs: $sOff, $eOff\n";

  if($eBin > $sBin){ #if the item is across bins
    if($eOff+$TENKB-$sOff+1 != $to-$from+1){
      die "Wrong calculation: $eOff+$TENKB-$sOff+1 != $to-$from+1\n";
    }

    my @firstHalf = getSliceInBin($bedsIdx[$sBin], $beds[$sBin], $sOff, $TENKB-1);
    my @secondHalf = getSliceInBin($bedsIdx[$eBin], $beds[$eBin], 0, $eOff);
    @readArray = (@firstHalf, @secondHalf);

  }else{
    @readArray = getSliceInBin($bedsIdx[$sBin], $beds[$sBin], $sOff, $eOff); 
  }
  if(@readArray != $to-$from+1){
    die "$from to $to: Wrong @readArray\n";
  }
  return @readArray;
}

sub getSliceInBin{
  
  my ($idxRef, $ref, $from, $to) = @_;
  my @idx = @$idxRef;
  my %reads = %$ref;
  my @res = ();
  for(my $i = 0; $i<$to-$from+1; $i++){ $res[$i] = 0;  }
  my ($loc, $unMatchFlag) = binarySearch($idxRef, $from, 0, @idx-1, 'left');
  # need to shift to the previous one if $from is < $idx[$loc] but > $idx[$loc-1] (assumed > the end position of $idx[$loc-2], if any)
  if($unMatchFlag == 1 && $loc>0){ $loc -= 1; }
 
  for(my $i = $loc; $i < @idx; $i++){
   my $s = $idx[$i];
   my ($e, $n) = split(',', $reads{$s}); 

   #get the intersect
   if($s > $to){ last; } #end already
   if($e < $from){ next; } #see the next one
   
   print "$s\t$e\t$n\t";
   
   if($e > $to){ $e = $to; }
   if($s < $from){ $s = $from; }
   for(my $j = $s; $j<=$e; $j++){
     if($res[$j-$from] >0){ die "Error overlap @res at pos $j-$from (the reads from bed.sam shouldn't overlap)\n";  }
     $res[$j-$from]=$n;
   }

  }
  return @res;
}

# sub routine to get the effective sum and effective length
#the input positions ($from, $to) are assumed to be 1-based (incl)
sub getEffReadSumLength
{
  my ($bedRef, $from, $to) = @_;
  my @readArray = getReadSlice($bedRef, $from, $to);
  my ($readSum, $readLen) = (0, 0);
  foreach(@readArray){ #i.e. @readArray
    if($_>0){
      $readSum += $_;
      ++$readLen;
    }
  }
  #print "reads: $readSum, len: $readLen\n";
  return ($readSum, $readLen);
}

############################################
### minor auxiliary (mainly for debug prints)

sub printBedReads{

  # read bed by chromosome
  #my $bedRef = readBedByChr($bedF, $genomeF, 10);
  my ($bedRef) = @_;
  my @beds = @{$bedRef->{'reads'}};
  my @bedsIdx = @{$bedRef->{'reads_idx'}};
  
  print "No of bins: ", scalar @beds, "\n";

  my $cnt=0;
  for(my $i=0; $i<@beds; $i++){
    if(defined($beds[$i])){
      my @idx = @{$bedsIdx[$i]};
      print "$i: idx: @idx\t"; print "\n";
      my $c = 0;
      my %hash = %{$beds[$i]};
      for(keys %hash){
        print "$_-$hash{$_}\t";
        $c++;
	if($c>20){ last; }
      }
      print "\n";
      $cnt++;
    }
    if($cnt > 20){ last; }
  }
  print "\n";

}

sub simpleTestReads{
  my ($bedF, $genomeF) = @_;
  my ($bedRef) = readBedByChr($bedF, $genomeF, 10);

  print "\ngetReadSlice\n";
  my ($c, $l) = getEffReadSumLength($bedRef, 179987, 180064);
  print "$c, $l\n";

  ($c, $l) = getEffReadSumLength($bedRef, 174710, 174750);
  print "$c, $l\n";

  ($c, $l) = getEffReadSumLength($bedRef, 135387008, 135387050);
  print "$c, $l\n";
}

1;
