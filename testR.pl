#! /usr/bin/perl/ -w
use strict;
use Statistics::R; #interact with R

# Create a communication bridge with R and start R
my $R = Statistics::R->new();

my ($tAllel1, $tAllel2, $cAllel1, $cAllel2) = (1, 9, 11, 3);

#use the R object to make a Fisher's exact test
$R->set('x', [$tAllel1, $tAllel2, $cAllel1, $cAllel2]);
my $xx = $R->get('x');
print "x: @$xx\n";

$R->run('xm = matrix(data = x, nrow=2)');

$xx = $R->get('xm');
print "xm: @$xx\n";

$R->run('p = fisher.test(xm)$p.value');
my $testFisher = $R->get('p');
print "fisher test result: $testFisher\n";
#my $pValue = @$testFisher[-1];   



$R->stop;
