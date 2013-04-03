#!/usr/bin/perl
use strict;
use warnings;

if(@ARGV < 2){
  print <<EOT;

USAGE: perl $0 rand_file_prefix rand_times

Check the overall averaged read count ratios of all SNVs across the rand_times ranomized cases prefixed with rand_file_prefix
The ratios are directly output so you need to redirect the output to a file if necessary.
The read count ratio for one SNV is calculated as 

r = sum(allel1/(allel1+allel2))/rand_times
sum is taken over 1,2, ..., rand_times

rand_file_prefix	the file prefix, e.g. "test/rand.snv." for 
			test/rand.snv.1, test/rand.snv.2, ... test/rand.snv.1000 (assume rand_times is 1000)
rand_times		the largest index of the rand_file suffix, e.g. 1000 means the randomized snvs will
			have 1, 2, 3, ..., 1000 as the suffixes
EOT

exit;
}


my ($files, $FREQ) = @ARGV;
#my $files = "test/rand.snv.";
#my $N = 1000000; #just check the first 10
my @ratios = ();
#my $FREQ = 1000;
my $SNVNO = 0;

for(my $i=1; $i<=$FREQ; $i++){
  #print "$i\n"; next;
  open(FP, $files.$i) or die "cannot open $files $i\n";
  my $j = 0;
  while(<FP>){
    chomp;
    my ($chr, $pos, $al, $id, $reads) = split(' ', $_);
    my ($r1, $r2) = split(':', $reads);
    #print "$r1, $r2\n";
    #exit;
    if($r1+$r2>=1){
      $ratios[$j] += $r1/($r1+$r2);
      $j++;
    }
    #if($j>=10){ last; }
  }
  if($i==1){
    $SNVNO = $j;
  }elsif($j != $SNVNO){
    die "Wrong number of SNVs ($SNVNO VS $j)\n";
  }
  close(FP);
}
print "There are $SNVNO non-trivial SNVs\n";

for(@ratios){
  printf "%.4f\n", $_/$FREQ;
}
