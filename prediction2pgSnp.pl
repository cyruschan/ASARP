#!/usr/bin/perl
use strict;
use warnings;

require "fileParser.pl"; #sub's for input annotation files
require "snpParser.pl"; #sub's for snps

if(@ARGV<4){
  print <<EOT;

USAGE: perl $0 prediction db outfolder annotations

This script generates the pgSNP tracks which can be uploaded to
USCS genome browser, most likely with the bedgraph tracks, for 
visual investigation of the ASARP cases.
associated with each pgSNP there will be an .aux file to specify
the appropriate range for visualization and provide more detailed
information of ASARP SNVs in the pgSNP tracks.
pgSNP information can be found at: ....

ARGUMENTS:
prediction	gene.prediction file of the ASARP pipeline
db		db (species code) for UCSC genome browser, e.g. hg19 or mm9
outfolder	output folder for the pgSNP files
annotations	transcript annotation file, which is also required for
		the ASARP pipeline (e.g. hg19.merged...txt)

EOT
  exit;
}

my ($input, $db, $outFolder, $xiaoF) = @ARGV;

open(FP, $input) or die "ERROR: cannot open $input to read\n";
my @pred = <FP>;
chomp @pred;
close(FP);

my %pgSnp = ();
my %pgAux = ();
my %chrStart = ();
my %chrEnd = ();


my $checkText = "ASARP gene level";
if($pred[0] ne $checkText){
  die "ERROR: expecting $input as gene.prediction (with header: $checkText)\n";
}

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
    
    # now gene transcript info ready, go to get pgSnp and aux
    while($i<@pred && $pred[$i] ne "" && !($pred[$i]=~/^chr/)){
      #get SNPs
      my ($catInfo, $snpInfo) = split(';', $pred[$i]);
      my($pos, $id, $al, $reads) = split(' ', $snpInfo);
      my ($al1, $al2) = split('>', $al);
      my ($r1, $r2) = split(':', $reads);
      my $snpTrack = join("\t", $chr, $pos-1, $pos, "$al1/$al2", 2, "$r1,$r2", "0,0"); #extra #, $gene, $catInfo);
      $pgSnp{$chr} .= $snpTrack."\n";

      #use @sortedEx and %exs to get the left and right flanking regions
      my ($exIdx, $unMatch) = binarySearch(\@sortedExs, $pos, 0, @sortedExs-1, 'left'); 
      my ($rangeL, $rangeR) = ($sortedExs[0], $exs{$sortedExs[-1]});
      if($exIdx > 0){ 
	$rangeL = $sortedExs[$exIdx-1];
      }
      if(!$unMatch){
        $exIdx++; #next one
      }
      if($exIdx <@sortedExs){
        $rangeR = $exs{$sortedExs[$exIdx]};
      }
      #print "$chr $gene $pos range: $rangeL to $rangeR\n";
      
      $pgAux{$chr} .= join("\t", $chr, $rangeL, $rangeR, "$gene $snpInfo type: $catInfo")."\n";
      
      if(!defined($chrStart{$chr})){  $chrStart{$chr} = $pos;  }
      elsif($chrStart{$chr} > $pos){	$chrStart{$chr} = $pos;	}
      
      if(!defined($chrEnd{$chr})){  $chrEnd{$chr} = $pos;  }
      elsif($chrEnd{$chr} < $pos){	$chrEnd{$chr} = $pos;	}
      
      $i++;
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
