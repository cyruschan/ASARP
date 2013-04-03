#!/usr/bin/perl
use strict;
use warnings;

if(@ARGV < 3){
  print <<EOT;

USAGE: perl $0 snv_line_no rand_file_prefix rand_times

Check the sample distribution of read counts of a certain SNV, provided its random folder and file prefix
the default output will be binom.snv_line_no.dist

snv_line_no		from which line number the SNV is to be investigated
rand_file_prefix	the file prefix, e.g. "test/rand.snv." for 
			test/rand.snv.1, test/rand.snv.2, ... test/rand.snv.1000 (assume rand_times is 1000)
rand_times		the largest index of the rand_file suffix, e.g. 1000 means the randomized snvs will
			have 1, 2, 3, ..., 1000 as the suffixes
EOT
}


my ($specNo, $files, $FREQ) = @ARGV;

if(!defined($specNo) || $specNo <=0){ 
  die "need to provide a line number for a particular SNV\n"; 
}
#my $files = "test/rand.snv.";
#my $N = 1000000; #just check the first 10
my @ratios = ();
#my $FREQ = 1000;
my $ttlReads = 0;
#my $specNo = 1;# just pick a number to investigate

my $j = 0;
for(my $i=1; $i<=$FREQ; $i++){
  open(FP, $files.$i) or die "ERROR: cannot open $files $i\n";
  my $lineNo = 0;
  while(<FP>){
    $lineNo++;
    if($lineNo < $specNo){
      next;
    }elsif($lineNo>$specNo){
      last;
    }
    chomp;
    my ($chr, $pos, $al, $id, $reads) = split(' ', $_);
    my ($r1, $r2) = split(':', $reads);
    #print "$r1, $r2\n";
    if($ttlReads == 0){
      $ttlReads = $r1+$r2;
    }else{
      if($ttlReads != $r1+$r2){
        die "ERROR: inconsistent total read count in rand $i case with $r1+$r2 != $ttlReads\n";
      }
    }
    #exit;
    $ratios[$j] = $r1;
    $j++;
    #if($j>=$N){ last; }
  }
  close(FP);
}
if($j!= $FREQ){ print STDERR "problem: $j != $FREQ\n"; }

my $outFile = "binom.$specNo.$ttlReads.dist";
open(OP, ">", $outFile) or die "ERROR: cannot open $outFile to output\n"; 
for(@ratios){
  print OP "$_\n";
}
close(OP);
