#!/usr/bin/perl
use warnings;
use strict;
require "randUtilities.pl";

if(@ARGV < 4){
  print <<EOT;
  
USAGE: perl $0 asarp_result rand_results indices output_p

NOTE: evaluation of p-values for the individual ASARP SNV results based on randomization for certain indices

ARGUMENTS:
asarp_result		ASARP prediction output, in particular
			.gene.prediction is expected
rand_results		the output folder and file prefix for the randomized ASARP SNV result files
			ASARP consistent suffix is assumed for result files, i.e. .gene.prediction
			e.g. "result/rand" for all "result.rand.*.gene.prediction" files
indices			The number of randomized SNV files, or the index range 
			for the randomized batch
			e.g. "20" means 20 randomized SNV lists indexed by 0 to 19
			e.g. "20:29" means 10 randomized SNV lists indexed by 20 to 29
			i.e. rand_results.20.gene.prediction, 
			     rand_results.21.gene.prediction, ...
			indice range >= N is needed to get enough power
			for p-value resolution ~1/N (e.g. N=1000 for p~0.001)
output_p		output ASARP SNV results with p-values associated

EOT
  exit;
}


my ($input, $randF, $freq, $outputF) = @ARGV;

my ($from, $to) = getIndex($freq);
my $N = $to-$from+1; #the total number for simulations

# get the real ASARP SNV results
my @snvTypes = ();
my ($snvsRef, $sumsRef, $snvHsRef) = getAsarpResult($input);
my $asarpSnvCntRef = getAsarpTypeCounts(\@snvTypes, $snvsRef, $snvHsRef);
#showSnvTypeCnt($snvsRef, $asarpSnvCntRef);

# get the randomized ASARP SNV results which correspond to @snvs
# the difference is that, rather than using
# $randSnvsRef as the list (which may contain SNVs beyond $snvsRef),
# $snvsRef is always used as the SNV list

my @randTypes = ();
my $randomSnvCntRef = \@randTypes; # accumulates all the AI/AS/AT/COMP counts for the randomized SNVs
for(my $i=$from; $i<=$to; $i++){
  my $randFileName = "$randF.$i.gene.prediction";
  my ($randSnvsRef, $randSumsRef, $randSnvHsRef) = getAsarpResult($randFileName);
  $randomSnvCntRef = getAsarpTypeCounts($randomSnvCntRef, $snvsRef, $randSnvHsRef); #use $snvsRef as the list
}
#showSnvTypeCnt($snvsRef, $randomSnvCntRef);

# get p-values
my $pRef = computePValues($snvsRef, $asarpSnvCntRef, $randomSnvCntRef, $N);

my %sumP = ();
my %pCnt = ();
for(qw(AI AS AT COMP)){
  $sumP{$_} = 0;
  $pCnt{$_} = 0;
}
my %pHash = %$pRef;
for(keys %pHash){
  #print "$_:\n";
  my @info = split('\t', $pHash{$_});
  for(@info){ 
    my ($p, $chr, $pos, $al, $id, $type) = split(' ', $_);
    $sumP{$type} += $p;
    $pCnt{$type} += 1;
    #print "$_\n"; 
  }
  #print "\n";
}

print "Overall ASARP FDR estimated:\n";
for(keys %pCnt){
  if($pCnt{$_}>0){
    printf "%s: %.4f\n", $_, $sumP{$_}/$pCnt{$_};
  }
}
print "\n";

# result printed to file
print "Output P-values of ASARP SNVs to $outputF\n";
open(OP, ">", $outputF) or die "ERROR: cannot open output: $outputF\n";
for(keys %pHash){
  print OP "$_:\n";
  my @info = split('\t', $pHash{$_});
  for(@info){ 
    print OP "$_\n"; 
  }
  print OP "\n";
}
close(OP);
print "Finished\n";

# compute p-values for each SNV at each gene
sub computePValues
{
  my ($snvsRef, $asarpRef, $randomRef, $N) = @_;
  my %geneResults = ();

  my @snvs = @$snvsRef;
  my @asarp = @$asarpRef;
  my @random = @$randomRef;
  for(my $i=0; $i<@snvs; $i++){
    my %aHs = %{$asarp[$i]};
    my %rHs = ();
    if(defined($random[$i])){
      %rHs = %{$random[$i]};
    }
    for(keys %aHs){
      my $p = 0;
      if(defined($rHs{$_})){ # non-zero p-value
        $p = $rHs{$_}/$N;
      }
      my ($gene, $type) = split('\t', $_);
      $geneResults{$gene} .= "$p $snvs[$i] $type\t";
    }
  }

  return \%geneResults;
}

#### sub-routines ####
# get all ASARP type counts for SNVs (both real and randomized ASARP ones)

sub getAsarpTypeCounts
{
  my ($snvTypesRef, $snvsRef, $snvHsRef) = @_;
  my @snvTypes = @$snvTypesRef;
  my @snvs = @$snvsRef;
  my %snvHs = %$snvHsRef;

  my @types = qw(AI AS AT);
  for(my $i=0; $i<@snvs; $i++){
    my %hs = ();
    if(defined($snvTypes[$i])){
      %hs = %{$snvTypes[$i]}; 
    }
    if(defined($snvHs{$snvs[$i]})){
      my @geneInfo = split('\t', $snvHs{$snvs[$i]});
      #if(@geneInfo > 1){  print "geneInfo\n@geneInfo\n"; }
      for my $info (@geneInfo){
          my $typeCnt = 0; #type count specific to this SNV and this Gene
	  my ($anno, $chr, $gene) = split(' ', $info);
          for my $type (@types){
            # unique info: gene + type
	    if($anno =~ /$type\:/){ # this type
              $hs{"$gene\t$type"} += 1;
	      $typeCnt++;
            }
	  }
          if($typeCnt > 1){ #COMP
            $hs{"$gene\tCOMP"} += 1; 
          }
      }
    }
    $snvTypes[$i] = \%hs; #get reference
  }
  return \@snvTypes;
}

# show the SNVs and the corresponding type counts gather
sub showSnvTypeCnt{
  my ($snvsRef, $asarpSnvCntRef) = @_;
  my @snvs = @$snvsRef;
  my @snvsCnts = @$asarpSnvCntRef;
  for(my $i = 0; $i< @snvs; $i++){
    print "$snvs[$i]:\n";
    my %hs = %{$snvsCnts[$i]};
    for(keys %hs){
      print "$_ ($hs{$_})\n";
    }
  }
  print "\n";
}

