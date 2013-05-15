#!/usr/bin/perl
use strict;
use warnings;

require "fileParser.pl"; #sub's for input annotation files
require "snpParser.pl"; #sub's for snps

if(@ARGV<4){
  print <<EOT;

USAGE: perl $0 asarp_file db outfolder annotations

This script streamlines significant ASARP SNV cases (i.e.
keeping only the target-control ASARP SNV pairs with the 
lowest p-values) and provides a compact view for UCSC 
genome browser visualization

It generates pgSNP tracks and the associated custom .aux files to specify
the appropriate ranges for visualization and to provide detailed
information of ASARP SNVs in the pgSNP tracks.
pgSNP information can be found at: 
http://main.genome-browser.bx.psu.edu/FAQ/FAQformat.html#format10

ARGUMENTS:
asarp_file	the ASARP result file, serving as the prefix for accessing
		"asarp_file.gene.prediction" and 
		"asarp_file.controlSNV.prediction"
		of the ASARP pipeline
db		db (species code) for UCSC genome browser, e.g. hg19 or mm9
outfolder	output folder for the pgSnp and .aux files
annotations	transcript annotation file, which is also required for
		the ASARP pipeline (e.g. hg19.merged...txt)

EOT
  exit;
}

my ($input, $db, $outFolder, $xiaoF) = @ARGV;

open(FP, "$input.controlSNV.prediction") or die "ERROR: cannot open $input.controlSNV.prediction to read\n";
my @pred = <FP>;
chomp @pred;
close(FP);
my $checkText = "ASARP control SNVs";
if($pred[0] ne $checkText){
  die "ERROR: expecting $input with controlSNV.prediction (with header: $checkText)\n";
}

open(GP, "$input.gene.prediction") or die "ERROR: cannot open $input.gene.prediction to read\n";
my @predGenes = <GP>;
chomp @predGenes;
close(GP);
$checkText = "ASARP gene level";
if($predGenes[0] ne $checkText){
  die "ERROR: expecting $input with gene.prediction (with header: $checkText)\n";
}

print "Getting the target SNV information from gene.prediction...\n";
my %pgSnp = ();
my %geneTrgts = (); #just store the target information
for(my $i = 1; $i < @predGenes; $i++){
  if($predGenes[$i] =~ /^chr/){ #starting of a new gene
    my ($chr, $gene) = split('\t', $predGenes[$i]);
    $i++; #go to the next line
    # get the target SNVs
    while($i<@predGenes && $predGenes[$i] ne "" && !($predGenes[$i]=~/^chr/)){
      #get SNPs
      my ($catInfo, $snpInfo) = split(';', $predGenes[$i]);
      my($pos, $id, $al, $reads) = split(' ', $snpInfo);
      $geneTrgts{"$gene;$pos"} = $snpInfo;
      
      my ($al1, $al2) = split('>', $al);
      my ($r1, $r2) = split(':', $reads);
      my $snpTrack = join("\t", $chr, $pos-1, $pos, "$al1/$al2", 2, "$r1,$r2", "0,0"); #extra #, $gene, $catInfo);
      $pgSnp{$chr} .= $snpTrack."\n";
      
      $i++;
    }
  }
}
# finished (just want the SNP ID, allelels, etc.)

my %pgAux = ();
my %chrStart = ();
my %chrEnd = ();

# read the transcript annotation file
my $transRef = readTranscriptFile($xiaoF);
#printListByKey($transRef, 'trans'); #utility sub: show transcripts (key: trans)
my ($genesIdxRef, $geneStartsRef, $geneEndsRef) = getGeneIndex($transRef); #get indices of gene transcript starts and gene names
my @geneStarts = @$geneStartsRef;
my @geneEnds = @$geneEndsRef;

print "Finding gene and flanking regions for the SNVs...\n";
my $curGene = '';
my %exs = (); #the exons for curGene
my @sortedExs = (); #the sorted keys

for(my $i = 1; $i < @pred; $i++){
  if($pred[$i] =~ /^chr/){ #starting of a new gene
    my ($chr, $gene) = split('\t', $pred[$i]);
    $i++; #go to the next line
      
    #try to get information for the Aux file to determine the range to enquiry
    my $chrId = getChrID($chr);
    my ($chrTransRef, $chrTransIdxRef) = getListByKeyChr($transRef, 'trans', $chrId);
    my %chrTrans = %$chrTransRef;
    my @chrTransIdx = @$chrTransIdxRef;
     
    if($gene ne $curGene){ #not yet parsed
      $curGene = $gene; #update current gene
      #get all the information for the current gene
      my $geneS = $geneStarts[$chrId]->{$gene};
      my $geneE = $geneEnds[$chrId]->{$gene};
     
      print "searching for $gene start:end $geneS:$geneE\n";
      my ($loc, $unMatchFlag) = binarySearch($chrTransIdxRef, $geneS, 0, @chrTransIdx-1, 'left');  
      # it must match otherwise the TransIdx is wrong
      if($unMatchFlag){
        die "ERROR: $gene start: $geneS must match some transcripts in $chr\n";
      }
      print "$gene starts in idx[$loc]\n"; #: $chrTransIdx[$loc] of \n"; #@chrTransIdx\n";
      
      #trans structure from fileParser.pl
      #  $trans[$chrId]{$txStart}.=$txEnd.";".$cdsStart.";".$cdsEnd.";".$exonStarts.";".$exonEnds.";".$ID.";".uc($geneName).";".$isCoding.";".$strand."\t";
      my $ti = $loc;
      %exs =(); #re-initialize
      @sortedExs = ();
      while($ti < @chrTransIdx && $chrTransIdx[$ti] <= $geneE){
        if($chrTrans{$chrTransIdx[$ti]} =~ /;$gene;/){
          my @allEnds = split('\t', $chrTrans{$chrTransIdx[$ti]});	
	  for(@allEnds){
	    my ($txEnd, $cdsS, $cdsE, $exonSs, $exonEs, $txId, $txGene) = split(';', $_);
	    if($txGene ne $gene){ next; }
	    my @exonSSet = split(',', $exonSs);
	    my @exonESet = split(',', $exonEs);
	    for(my $j = 0; $j < @exonSSet; $j++){
	      if(!defined($exs{$exonSSet[$j]}) || $exs{$exonSSet[$j]} < $exonESet[$j]){
	        $exs{$exonSSet[$j]} = $exonESet[$j];
	      }
	    }
	  }
	}
	$ti++;
      }
      @sortedExs = sort {$a <=> $b} keys %exs;
      #print "$gene exs:\n@sortedExs\n";
    }
    
    # now gene transcript info ready, go to get aux
    my %coms = (); # compact results, for the same target SNV, keep only the lowest p-value
    my %bestP = (); #the best p-value for each target SNV
    while($i<@pred && $pred[$i] ne "" && !($pred[$i]=~/^chr/)){
      #get target-control pairs
      my ($catInfo, $pVal, $trgtPos, $ctrlPos) = split(';', $pred[$i]);
      if(defined($bestP{$trgtPos}) && $bestP{$trgtPos} <= $pVal){ #no need to look at this case
	next; # no need to see this any more
      }
      #update info
      $bestP{$trgtPos} = $pVal;

      #use @sortedEx and %exs to get the left and right flanking regions
      my ($rangeL, $rangeR) = getBoundary($trgtPos, \%exs, \@sortedExs, 1);
      my ($ctrlRangeL, $ctrlRangeR) = getBoundary($ctrlPos, \%exs, \@sortedExs, 0); # no need to get flanking regions as NEV is not required for ctrlPos
      if($ctrlRangeL < $rangeL){ $rangeL = $ctrlRangeL; }
      if($ctrlRangeR > $rangeR){ $rangeR = $ctrlRangeR; }

      # now we need %geneTrgts
      my $snpInfo = $geneTrgts{"$gene;$trgtPos"};
      $coms{"$gene;$trgtPos"} = join("\t", $chr, $rangeL, $rangeR, "$gene $snpInfo (vs $ctrlPos p: $pVal) type: $catInfo");
      
      if(!defined($chrStart{$chr})){  $chrStart{$chr} = $trgtPos;  }
      elsif($chrStart{$chr} > $trgtPos){	$chrStart{$chr} = $trgtPos;	}
      
      if(!defined($chrEnd{$chr})){  $chrEnd{$chr} = $trgtPos;  }
      elsif($chrEnd{$chr} < $trgtPos){	$chrEnd{$chr} = $trgtPos;	}
      
      $i++;
    }
    for(keys %coms){
      $pgAux{$chr} .= $coms{$_}."\n";
    }
  }
}
print "\nOutput pgSnp and aux to $outFolder\n";
for(keys %pgSnp){
  my $chr = $_;
  print "$chr\n";
  open(OP, ">", "$outFolder/$chr.pgSnp") or die "ERROR: cannot open $outFolder/$chr.pgSnp for writing\n";
  print OP "track type=pgSnp visibility=3 db=$db name=\"SNV_$chr\" description=\"ASARP SNVs in $chr\"\n";
  print OP "browser position $chr:$chrStart{$chr}-$chrEnd{$chr}\n";
  print OP $pgSnp{$chr};
  close(OP);
  
  open(AP, ">", "$outFolder/$chr.pgSnp.aux") or die "ERROR: cannot open $outFolder/$chr.pgSnp.aux for writing\n";
  print AP $pgAux{$chr};
  close(AP);
}

sub getBoundary{

   my ($pos, $exsHsRef, $srtExsArrRef, $jump) = @_;
   # jump is the flag to specify whether you want to get only the boundaries
   # surrounding the position (0) or 
   # the left and right boundaries of the previous and next exons (1)
   my @sortedExs = @$srtExsArrRef;
   my %exs = %$exsHsRef;

      my ($exIdx, $unMatch) = binarySearch(\@sortedExs, $pos, 0, @sortedExs-1, 'left'); 
      my ($rangeL, $rangeR) = ($sortedExs[0], $exs{$sortedExs[-1]}); #initialize $rangeL and $rangeR to the possible min and max (in case no matches)
      my $lIdx = 0;
      if($exIdx > 0){  $lIdx = $exIdx - 1; #the previous exon start
      }
      $rangeL = $sortedExs[$lIdx]; # lIdx;
      my $jExIdx = $exIdx;
      if(!$unMatch){ # it means it matches 
        $jExIdx++; #next one
      }
      if($jExIdx <@sortedExs){
        $rangeR = $exs{$sortedExs[$jExIdx]};
      }
      #print "$chr $gene $pos range: $rangeL to $rangeR\n";
      if($jump){
        return ($rangeL, $rangeR);
      }

      # jump == 0 need to get the accurate range of the exon start and end (if it is on an exon):
      if(!$unMatch){
        return ($sortedExs[$exIdx], $exs{$sortedExs[$exIdx]}); # it's directly on an exon
      }
      return ($sortedExs[$lIdx], $exs{$sortedExs[$lIdx]});
}      
