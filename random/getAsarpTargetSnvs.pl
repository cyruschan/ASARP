#!/usr/bin/perl
use warnings;
use strict;
require "randUtilities.pl";

if(@ARGV < 4){
  print <<EOT;
  
USAGE: perl $0 ctrl_snvs rand_infolder rand_outfolder indices [seed]

NOTE: generation of random SNVs (only target and control SNVs output by the ASARP pipeline) for certain indices

ARGUMENTS:
ctrl_snvs		file of controlSNV output of the ASARP pipeline
			with the extension: controlSNV.prediction 
rand_infolder		the input folder for the randomized SNV files (already done with randAllSnvs.pl)
			all rand.SNV.* files will be read and compared with the list from ctrl_snvs
rand_outfolder		the output folder for the randomized SNV files
indices			The number of randomized SNV files, or the index range 
			for the randomized batch
			e.g. "20" means 20 randomized SNV lists indexed by 0 to 19
			e.g. "20:29" means 10 randomized SNV lists indexed by 20 to 29
			i.e. rand_folder/rand.SNV.20, rand_folder/rand.SNV.21, ...
OPTIONAL:
seed			fixed random seed for debugging

EOT
  exit;
}



my ($ctrlSnvs, $inFolder, $outFolder, $freq, $seed) = @ARGV;

if(defined($seed)){
  print "\nIMPORTANT: setting seed for srand will NOT generate RANDOM results.\n";
  print "           this seed setting is used only for debugging.\n";
  srand($seed);
}

if($inFolder eq $outFolder){
  die "ERROR: input random folder cannot be the same as the output one: $inFolder\n The input rand.snv will all be overwritten!\n";
}

my ($from, $to) = getIndex($freq);

# get the SNV list of input snvs
my %snvs = ();
open(FP, $ctrlSnvs) or die "ERROR: cannot open controlSNV.prediction file: $ctrlSnvs\n";
while(<FP>){
  chomp;
  if($_ =~/;/){
    my ($types, $p, $trgt, $ctrl) = split(';', $_);
    #print "$trgt $ctrl\n";
    $snvs{$trgt} = 1;
    $snvs{$ctrl} = 1;
  }
}
close(FP);

my $snvNo = keys %snvs;
print "There are $snvNo for the ASARP target and control SNVs\n";

# get all rand_inFolder SNVs
# assume all the rand_inFolder SNV lists are the same (differ only in the allele reads)
my $ff = 0; # the flag indicating whether the SNV filter has been done
my @filters = ();
for(my $i=$from; $i<=$to; $i++){
  my $file = "$inFolder/rand.snv.$i";
  my $toOutput = "";
  open(IP, $file) or die "ERROR: cannot open $file\n";
  if(!$ff){ #need to read at least one file to get the filtered list;
    my $j = 0;
    while(<IP>){
      # no need to chomp because we don't care the last $read
      # in this case, we need to parse the input file
      my ($chr, $pos, $al, $id, $read) = split(' ', $_);
      if(defined($snvs{$pos})){
        $toOutput .= $_; 
        push @filters, $j; #filters has only the line number
      }
      $j++;
    }
    $ff = 1; #the filters list has been done

  }else{
    my ($j, $fi) = (0, 0);
    while($fi<@filters){
      my $line = <IP>; #read line
      if($j == $filters[$fi]){ #filters is an ordered subset of the SNV lines
        # we can actually double-check if this line is correct
	my ($chr, $pos, $al, $id, $read) = split(' ', $line);
	if(!defined($snvs{$pos})){
	  die "ERROR: in line $fi+1 SNV: $pos is not included in the target/control ASARP SNV list (probably input random files are not consistent in their SNV lists)\n";
	}
	$toOutput .= $line; 
	$fi++;
      }
      $j++;
    }
  
  }
  close(IP);

  # output the file:
  my $outFile = "$outFolder/rand.snv.$i";
  open(OP, ">", $outFile) or die "ERROR: cannot write $outFile\n";
  print OP $toOutput;
  close(OP);
}

print "Finished\n";
