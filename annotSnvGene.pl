#!/usr/bin/perl
use warnings;
use strict;

require "fileParser.pl"; #sub's for input annotation files
require "snpParser.pl"; #sub's for snps
use MyConstants qw( $CHRNUM $supportedList $supportedTags );

if(@ARGV < 2){
  #die "Need snp_file, hg19_file\n";
  print <<EOT;
USAGE: perl $0 input_snv_file gene_annotation_file

This script is to fetch gene region annotation for SNVs in detail. 
The results will be organized by genes. For each gene, the SNVs 
will be output, along with all of their hosting transcript IDs, 
their region annotation type, and the neighboring exons.

ARGUMENTS:

input_snv_file		SNV input file in the following format
			chr position [any additional fields] strand (examples below)
			chr10 97624709 +
			chr10 97625099 G>A rs4512761 267:186:0 +
			(strand should be the last attribute)

gene_annotation_file	the gene transcript annotation file used in the ASARP
			config, i.e. "xiaofile". 
			Example: data/hg19.merged.to.ensg.all.tx.03.18.2011.txt 

NOTE: currently only strand-specific SNV files are supported
      for non strand specific input, user has to add the strand
      attribute: +/- at the end of each SNV accordingly (e.g. the 
      strand of the gene containing the SNV)
      All the fields, if any, between the position and the strand are ignored
      e.g. only the 1st, 2nd and 6th attributes from the following SNV is used 
      chr10 97624709 A>C rs2275759 111:137:0 +

      The results will be output to screen and they can be redirected
      to an output file via adding ">> output_file"

EOT
  exit;
}
my ($snpFile, $xiaoF, $output) = @ARGV;

# Demonstration of using the fileParser
# read the transcript annotation file
my ($transRef) = readTranscriptFile($xiaoF);
my ($transGeneRef) = readTranscriptFileByGene($xiaoF);
#printListByKey($transRef, 'trans'); #utility sub: show transcripts (key: trans)

my ($snpRef) = initFakeSnpSS($snpFile, '+');
my ($snpRcRef) = initFakeSnpSS($snpFile, '-');

# check annotations with SNVs
my ($geneSnpRef) = setGeneSnps($snpRef, $transRef, '+');
my ($geneSnpRcRef) = setGeneSnps($snpRcRef, $transRef, '-');
print "ALL + STRAND CASES\n";
printIntronGeneSnpsResults($geneSnpRef,'gPowSnps', $snpRef,'powSnps', $transGeneRef, '+');

print "ALL - STRAND CASES\n";
printIntronGeneSnpsResults($geneSnpRcRef,'gPowSnps', $snpRcRef,'powSnps', $transGeneRef, '-');



######################################################################################

# a fake version just to get all gene names for the SNVs
sub initFakeSnpSS{
  #snps and powSnps
  my @pList = ();
  my @snps = (); my @powSnps = ();
  for(my $i=0; $i<=$CHRNUM; $i++){
    push @snps, {};
    push @powSnps, {};
  }
  
  
  my ($snpFile, $strandType) = @_;
  open(my $fh, "<", $snpFile) or die "Cannot open snp file: $snpFile for reading.\n";
  print "Reading from $snpFile\n";

  # Create a communication bridge with R and start R
  #my $R = Statistics::R->new();
  
  my $count = 0;
  while(<$fh>){
    $count++;
    if(!$count%1000){
      #print $count." ";
    }
    chomp;
    my ($chrRaw, $pos, @remain)=split(/ /, $_);
    if(!defined($remain[-1])){
      die "ERROR: missing strand at the last attribute in: $_\n";
    }
    my $strandInLine = $remain[-1]; #last attribute
    #print "$_\n";
    #strand-specific handling
    if(defined($strandType)){
      if(!defined($strandInLine)){
        die "ERROR: SNV data must contain strand (+/-) at the end when strand-specific flag is set\n";
      }else{
        if($strandType ne $strandInLine){ #handle only SNVs with the specified $strandType
	  next;
	}
      }
    }else{ #setting is non-strand specific
      if(defined($strandInLine)){ # there will be errors if strand specific data are handled in a non-strand specific way
        die "ERROR: strand-specific SNV data are not handled when strand=specific flag is unset\n";
      }
    }

    my $chrID = getChrID($chrRaw); #auxiliary from fileparser.pl
    #check numeric
    if(!($chrID=~/^\d+$/)){
      #print "$chrID\n"; 
      next; 
    }
    #do the R Chi-squared Test
    my $infoKept = "$chrRaw $pos";
    $powSnps[$chrID]{$pos} = $infoKept;	

  }
  close($fh);
  #$R->stop;

  #create array indices
  my @snp_idx = (); my @powSnps_idx = ();
  for(my $i=1; $i<=$CHRNUM; $i++){
     $snp_idx[$i] = [sort {$a<=>$b} keys %{$snps[$i]}];
     $powSnps_idx[$i] = [sort {$a<=>$b} keys %{$powSnps[$i]}];
  }

  my %snpList = (
    'snps' => \@snps,
    'powSnps' => \@powSnps,
    'snps_idx' => \@snp_idx,
    'powSnps_idx' => \@powSnps_idx,
  );

  return (\%snpList);
}

# a fake one to print out the intron things only
# print out the gene snps results for certain type (powerful or ordinary snps)
sub printIntronGeneSnpsResults
{
  my ($geneSnpRef, $geneSnpKey, $snpRef, $snpKey, $transGeneRef, $strand) = @_;
  #print "Gene level SNP ($geneSnpKey) VS transcript results\n";

  for (my $i=1; $i<=$CHRNUM; $i++){
    my ($geneRef) = getListByKeyChr($geneSnpRef, $geneSnpKey, $i); 
    my $geneChrRef = getChrGeneSnpsSorted($geneSnpRef, $geneSnpKey, $i);
    my ($snpInfoRef) = getListByKeyChr($snpRef, $snpKey, $i);
    my %snps = %$snpInfoRef; #snps information
    my @genes = @$geneChrRef;
    my %genesPrinted = (); #to store genes printed already in order not to double-print
    if(@genes){
      my $chr = formatChr($i); #print "\n";
      #print "@genes\t";
      foreach(@genes){
        my @allGenes = split('\t', $_);
        foreach(@allGenes){
	  #get the snp information
	  my @geneMatches = split('\t', $geneRef->{$_});
	  for(@geneMatches){
	    my ($snpPos, $matchInfo, $geneName, $txStart, $id, $regStart, $regEnd) = split(';', $_);
	    if(!defined($genesPrinted{$geneName})){
	      print "$chr\t$geneName\n";
	      $genesPrinted{$geneName} = 1;
	    }else{ $genesPrinted{$geneName} += 1; }
	    #print $snps{$snpPos}."\n";
            my ($prevRef, $nextRef) = findNextExons($chr, $snpPos, $geneName, $id, $transGeneRef, $strand); 
	    my $prevE = join(";", keys %$prevRef); 
	    my $nextE = join(";", keys %$nextRef); 
	    my @allSnpInfo = split(';', $snps{$snpPos}); #separate by ;, if there are multiple snps at the same position
            foreach(@allSnpInfo){
	      my $toPrint = $matchInfo."\t".$geneName."\t".$id.":".$regStart."-".$regEnd."\t$prevE\t$nextE\n";
	      print $_."\t".$toPrint;
	    }
	  }
        }
      }
      print "\n";
    }
  }

}


# get annotations arranged in a chr, gene manner
sub readTranscriptFileByGene
{

  #array of hashes
  my @trans=();
  for(my $i=0; $i<=$CHRNUM; $i++){
    push @trans, {}; #empty initialization
  }

  #keep track of discarded chromosomes
  my %discardedChrs = ();

  my ($fileName)=@_;
  open(my $fh, "<", $fileName) or die "Cannot open merged transcrimptome file: $fileName\n";
  print "Reading the transcriptome file ", $fileName, "\n"; 
  
  my $count = 0;
  while(<$fh>){
    $count++;
    if(!($count%10000)){
      #print $count, "\t";
      #STDOUT->autoflush(1);#need to flush if STDOUT not attached to terminal
    }
    chomp;
    my ($ID, $chrRaw, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd, 
    $exonCount, $exonStarts, $exonEnds, $geneName, $cdsStartStat, $cdsEndStat)
    = split(/\t/, $_);# ID, chr, strand, txStart, txEnd, 
    #cdsstart, cdsend, exoncount, exonstarts, 
    # exonends, genename, cdsstartstat,cdsendstat
    my $chrID = getChrID($chrRaw);
    if(!($chrID=~/^\d+$/)){
      if(!defined($discardedChrs{$chrID})){
        $discardedChrs{$chrID}=1;
      }else{ $discardedChrs{$chrID}+=1;  }
      #print "$chrID\n"; 
      #print "Not valid chr ID: $chrID\n"; 
      next;
    }
    $ID = uc $ID; # to upper
    #change 0-based to 1-based for all starts
    $txStart+=1;
    $cdsStart+=1;
    my @allExonStarts = split(/,/, $exonStarts);
    for(my $i=0; $i< @allExonStarts; $i++){
      $allExonStarts[$i]+=1;
    }
    $exonStarts = join(",", @allExonStarts);
    #but no need to add 1 for ends
    $exonEnds = substr($exonEnds, 0, -1); #because there is an additional "," in the original file

    my $isCoding = 1;
    #if($geneName=~/noncoding/){  $isCoding = 0; }
    if($cdsStart > $cdsEnd){ $isCoding = 0; }
    #if($isCoding==0){   print "Non-coding: $ID, $geneName\n"; }
    # currently skip calculating the 5UTR 3UTR and exon length here
    $geneName = uc $geneName; 
    $trans[$chrID]{$geneName}.=$txStart.";".$txEnd.";".$cdsStart.";".$cdsEnd.";".$exonStarts.";".$exonEnds.";".$ID.";".$isCoding.";".$strand."\t";
  }
  close($fh);

  #set the indices (such that event starts are sorted) for chrs
  my @trans_idx = ();

  for(my $i=1; $i<=$CHRNUM; $i++){
    #printChr($i); 
    $trans_idx[$i] = [sort keys %{$trans[$i]}]; #now it's sorted by gene name.
  }
  #debug print
  #my $chrToPrint = 4;
  #printChr($chrToPrint); print "\n";
  #foreach (@{$trans_idx[$chrToPrint]}){ print $_-1,"\n"; }
  
  #printDiscardedChrs(\%discardedChrs);  
  
  my %transList = (
    'trans' => \@trans,
    'trans_idx' => \@trans_idx,
    'type' => 'transcripts (arranged by genes); not the one arranged by transcript starts in ASARP',
  );
  return \%transList;
}

# find the next exons given a position
sub findNextExons{
  my ($chr, $pos, $gene, $transId, $transGeneRef, $strand) = @_;
  my $chrId = getChrID($chr);

  my %prevE = (); #all previous exons
  my %nextE = (); #all next exons
  
  my ($chrTransRef, $chrTransIdxRef) = getListByKeyChr($transGeneRef, 'trans', $chrId);
  my %trans = %$chrTransRef; #chromosome specific

  if(defined($trans{$gene})){
    my @transcripts = split('\t', $trans{$gene});
    #$trans[$chrID]{$geneName}.=$txStart.";".$txEnd.";".$cdsStart.";".$cdsEnd.";".$exonStarts.";".$exonEnds.";".$ID.";".$isCoding.";".$strand."\t"
    for(@transcripts){
      my ($txS, $txE, $cdsS, $cdsE, $exSs, $exEs, $id, $isCode, $strandTx) = split(';', $_);
      if($transId ne $id){ next; } # only check per-ID transcripts

      if(defined($strand) && $strandTx ne $strand){
        next;
	#die "Strand inconsistent: $chr $pos $strand while $gene has transcript $id on $strandTx\n";
      }
      my @exStarts = split(',', $exSs);
      my @exEnds = split(',', $exEs);
      my ($loc, $unMatchFlag) = binarySearch(\@exStarts, $pos, 0, @exStarts-1, 'left');
      #print "\n index $loc ($unMatchFlag) of pos $pos  from @exStarts\n";
      
      if(!$unMatchFlag){ #matches, that's not possible
        die "$chr $pos should be an INTRONIC SNV, should not match any exons\n";
      }else{
        if($loc > 0 && $loc < @exStarts){
	  my $key = "$exStarts[$loc-1],$exEnds[$loc-1]";
	  if(defined($prevE{$key})){
	   $prevE{$key} .= "\t$id";
	  }else{
	   $prevE{$key} = "$id";
	  }
	}
	if($loc < @exStarts){
	  my $key = "$exStarts[$loc],$exEnds[$loc]";
	  if(defined($nextE{$key})){
	   $nextE{$key} .= "\t$id";
	  }else{
	   $nextE{$key} = "$id";
	  }
	}
      }
    }
  }else{
    my $dieMsg = "$chr, $pos, $gene";
    if(defined($strand)){
      $dieMsg .=", $strand";
    }
    $dieMsg .= " not found\n";
    die $dieMsg;
  }

  return (\%prevE, \%nextE);
}


__END__

=head1 NAME

annotSnvGene.pl -- to fetch gene region annotation for SNVs in detail. 
The results will be organized by genes. For each gene, the SNVs 
will be output, along with all of their hosting transcript IDs, 
their region annotation type, and the neighboring exons.


=head1 SYNOPSIS

USAGE: 
 
 perl $0 input_snv_file gene_annotation_file

ARGUMENTS:

input_snv_file		SNV input file in the following format
			chr position [any additional fields] strand (examples below)
			chr10 97624709 +
			chr10 97625099 G>A rs4512761 267:186:0 +
			(strand should be the last attribute)

gene_annotation_file	the gene transcript annotation file used in the ASARP
			config, i.e. "xiaofile". 
			Example: data/hg19.merged.to.ensg.all.tx.03.18.2011.txt 

NOTE:

currently only strand-specific SNV files are supported
for non strand specific input, user has to add the strand
attribute: +/- at the end of each SNV accordingly (e.g. the 
strand of the gene containing the SNV)
All the fields, if any, between the position and the strand are ignored
e.g. only the 1st, 2nd and 6th attributes from the following SNV is used 
chr10 97624709 A>C rs2275759 111:137:0 +

The results will be output to screen and they can be redirected
to an output file via adding ">> output_file"

=head1 DESCRIPTIONS

It will first print out all + strand cases arranged by chromosomes and genes, and then 
all - strand cases. The output format will be 
 
 chr position region_type:[subtype:]strand gene transcript_ID:SNV_hosting_region[ left_exon right_exon]

 *_exon format: exon_start,exon_end (not applicable for UTR regions)
 
The output will look like this:

 ALL + STRAND CASES
 chr1    AK4
 chr1 65694819   exon:3'UTR:+    AK4     NM_203464:65691746-65697828
 chr1 65694819   exon:3'UTR:+    AK4     NM_001005353:65691746-65697828  
 chr1    AK3L1
 chr1 65694819   exon:3'UTR:+    AK3L1   ENST00000327299:65691746-65697819
 chr1 65695853   exon:3'UTR:+    AK3L1   ENST00000327299:65691746-65697819

 chr10   ENTPD1
 chr10 97624709  intron:+        ENTPD1  ENST00000453258:97624619-97625933       97624481,97624618       97625934,97629452
 chr10 97624709  intron:+        ENTPD1  UC001KLI.3:97624619-97625933    97624481,97624618       97625934,97637022
 chr10 97624709  intron:+        ENTPD1  NM_001098175:97624619-97625933  97624481,97624618       97625934,97637023
 ...
 ALL - STRAND CASES
 ...
 chr5 96130836   exon:normal:-   ERAP1   UC003KMN.2:96130745-96130865    96130745,96130865       96132878,96133012
 chr5 96139595   exon:normal:-   ERAP1   NM_016442:96139106-96139646     96139106,96139646       96143563,96143892
 chr5 96139595   exon:normal:-   ERAP1   UC003KML.2:96139106-96139646    96139106,96139646       96143563,96143892

Note that the results will be lengthy so it is better to analyze a small subset of selected SNVs using this script.

=head1 SEE ALSO

L<Overview>, L<asarp>

=head1 COPYRIGHT

This pipeline is free software; you can redistribute it and/or modify it given that the related works and authors are cited and acknowledged.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.

=head1 AUTHOR

Cyrus Tak-Ming CHAN

Xiao Lab, Department of Integrative Biology & Physiology, UCLA

=cut
