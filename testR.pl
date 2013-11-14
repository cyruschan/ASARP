#! /usr/bin/perl/ -w
use strict;
use Statistics::R; #interact with R


print "A simple test if Statistics::R is working\n";

# Create a communication bridge with R and start R
my $R = Statistics::R->new();

my ($tAllel1, $tAllel2, $cAllel1, $cAllel2) = (1, 9, 11, 3);

#use the R object to make a Fisher's exact test
$R->set('x', [$tAllel1, $tAllel2, $cAllel1, $cAllel2]);
my $xx = $R->get('x');
print "x: @$xx\n";

$R->run('xm = matrix(data = x, nrow=2)');

$xx = $R->get('xm');
#print "xm: @$xx\n";

print "Perform a fisher's exact test\n";

$R->run('p = fisher.test(xm)$p.value');
my $testFisher = $R->get('p');
print "fisher test result: $testFisher\n";
#my $pValue = @$testFisher[-1];   

if(defined($testFisher) && $testFisher <= 1 && $testFisher >= 0){
  print "Test SUCCESSFUL\n";
}else{
  print "Test FAILED (test result = $testFisher)\n";
}


$R->stop;
