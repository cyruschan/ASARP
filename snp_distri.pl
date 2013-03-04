#!/usr/bin/perl -w
use strict;
# a simple demonstration with the ASARP pipeline:
# in this application, the transcripts and SNVs are read in and the SNV distribution in 5'/3' UTRs, exons and introns is reported

# NOTE: one SNV can match several transcripts

use MyConstants qw( $CHRNUM $supportedList $supportedTags );
require "fileParser.pl"; #sub's for input annotation files
require "snpParser.pl"; #sub's for snps

# input arguments: $outputFile--output, $snpF--SNV file path, $xiaoF--transcript annotation file
# optional: $POWCUTOFF--powerful SNV cutoff, default: 20
my ($outputFile, $snpF, $xiaoF, $POWCUTOFF) = @ARGV;
if(@ARGV < 3){ # !defined($outputFile) || !defined($snpF) || !defined($xiaoF)){
  print "USAGE: perl $0 output_file snv_file transcript_file [powerful_snv_cutoff]\n";
  exit;
}
if(!defined($POWCUTOFF)){
  $POWCUTOFF = 20;
  print "Default powerful_snv_cutoff used: $POWCUTOFF\n";
}

my $fp = undef;
open($fp, ">", $outputFile) or die "ERROR: Cannot open $outputFile to write\n";
print $fp "Type\tExon\tIntron\t5'UTR\t3'UTR\tComplex\tIn-gene\tIntergenic\tTotal\n";


# Demonstration of using the fileParser
# read the transcript annotation file
my $transRef = readTranscriptFile($xiaoF);
#printListByKey($transRef, 'trans'); #utility sub: show transcripts (key: trans)

# Demonstration of using the snpParser
my $snpRef = initSnp($snpF, $POWCUTOFF);
#print "SNV List:\n";
#printListByKey($snpRef, 'powSnps');
#printListByKey($snpRef, 'snps');

my $geneSnpRef = setGeneSnps($snpRef, $transRef);
#print "Significant Snvs: \n";
#printGetGeneSnpsResults($geneSnpRef,'gPowSnps', $snpRef,'powSnps', 1); #$SNVPCUTOFF);
#print "Ordinary Snvs: \n";
#printGetGeneSnpsResults($geneSnpRef,'gSnps', $snpRef,'snps', 1);

print "Calculating powerful SNV distribution...\n";
my @part1 = getGeneSnpsDistri($geneSnpRef,'gPowSnps', $snpRef,'powSnps'); 
print $fp "Powerful(>=$POWCUTOFF)\t".join("\t",@part1)."\n";
printSnvDist(@part1);

print "Calculating non-powerful SNV distribution...\n";
my @part2 = getGeneSnpsDistri($geneSnpRef,'gSnps', $snpRef,'snps'); 
print $fp "Non-powerful(<$POWCUTOFF)\t".join("\t",@part2)."\n";
printSnvDist(@part2);

my @total = ();
for(my $i = 0; $i < @part1; $i++){
  $total[$i] = $part1[$i] + $part2[$i];
}
#the elements contained in the total array:
#($tEx, $tIn, $t5UTR, $t3UTR, $tCompPosNo, $tGeneSnpPosNo, $tIntergSnp, $allSnpNo) = @total;
print "\nOverall SNV distribution:\n";
printSnvDist(@total);
print $fp "Overall\t".join("\t",@total)."\n";

close($fp);

# print out the gene snps results for certain type (powerful or ordinary snps)
sub getGeneSnpsDistri
{
  my ($geneSnpRef, $geneSnpKey, $snpRef, $snpKey) = @_;

  my @dist = (); #snp distribution results
  for (my $i=0; $i<=$CHRNUM; $i++){
    push @dist, {};
  }
  my $allSnpNo = 0; #all SNV positions, including those without matching any gene transcripts

  for (my $i=1; $i<=$CHRNUM; $i++){
    my ($geneRef) = getListByKeyChr($geneSnpRef, $geneSnpKey, $i); 
    my $geneChrRef = getChrGeneSnpsSorted($geneSnpRef, $geneSnpKey, $i);
    my ($snpInfoRef) = getListByKeyChr($snpRef, $snpKey, $i);
    my %snps = %$snpInfoRef; #snps information
    $allSnpNo += keys %snps;
    my @genes = @$geneChrRef; 

    my %genesPrinted = (); #to store genes printed already in order not to double-print
    if(@genes){
      #print "@genes\t";
      foreach(@genes){
        my @allGenes = split('\t', $_);
        foreach(@allGenes){
	  #get the snp information
	  my @geneMatches = split('\t', $geneRef->{$_});
	  for(@geneMatches){
	    my ($snpPos, $matchInfo, $geneName, $txStart, $id, $regStart, $regEnd) = split(';', $_);
	    #use bit to represent $matchInfo types: exon, intron, 5'UTR, 3'UTR
            my($inEx, $inIn, $in5UTR, $in3UTR) = (0,0,0,0);

	    if($matchInfo =~ /exon/ && !($matchInfo =~ /UTR/)){ #UTR also has the exon: key word
	      $inEx = 1;
	    }
	    if($matchInfo =~ /intron/){
	      $inIn = 1;
	    }
	    if($matchInfo =~ /5'UTR/){
	      $in5UTR = 1;
	    }
	    if($matchInfo =~ /3'UTR/){
	      $in3UTR = 1;
	    }
            #print "($inEx, $inIn, $in5UTR, $in3UTR)\n";
	    #exit;
	    #for the snp position, get all the match info for it
	    if(!defined($dist[$i]{$snpPos})){
	      $dist[$i]{$snpPos} = join(',', $inEx, $inIn, $in5UTR, $in3UTR); 
	    }else{
	      my ($oEx, $oIn, $o5UTR, $o3UTR) = split(',', $dist[$i]{$snpPos});
	      if(!$oEx){ $oEx = $inEx;	}
	      if(!$oIn){ $oIn = $inIn;	}
	      if(!$o5UTR){ $o5UTR = $in5UTR;	}
	      if(!$o3UTR){ $o3UTR = $in3UTR;	}
	      $dist[$i]{$snpPos} = join(',', $oEx, $oIn, $o5UTR, $o3UTR);
	    }
	  }
        }
      }
    }
  }

  # clean-up for absolute concept of introns: only when a SNV never overlaps with any transcripts
  for(my $i=1; $i<=$CHRNUM; $i++){
    for(keys %{$dist[$i]}){
      my ($inEx, $inIn, $in5UTR, $in3UTR) = split(',', $dist[$i]{$_});
      if($inEx || $in5UTR || $in3UTR){ #once SNV matches some exon, 5'UTR, or 3'UTR, it cannot be considered as in an intron
        $inIn = 0;
      }
      $dist[$i]{$_} = join(',', $inEx, $inIn, $in5UTR, $in3UTR);
    }
  }

  #collect and sumamrize the results:
  my $tGeneSnpPosNo = 0;
  my $tCompPosNo = 0;
  my ($tEx, $tIn, $t5UTR, $t3UTR) = (0, 0, 0, 0);
  for(my $i=1; $i<=$CHRNUM; $i++){
    my $chr = formatChr($i);
    my %hs = %{$dist[$i]};
    for(keys %hs){
      $tGeneSnpPosNo++;
      my ($inEx, $inIn, $in5UTR, $in3UTR) = split(',', $hs{$_});
      $tEx += $inEx;
      $tIn += $inIn;
      $t5UTR += $in5UTR;
      $t3UTR += $in3UTR;
      if($inEx + $inIn + $in5UTR + $in3UTR > 1){
        $tCompPosNo ++;
      }
    }
  }
  # No of intergenic SNV positions
  my $tIntergSnp = $allSnpNo - $tGeneSnpPosNo;

  return ($tEx, $tIn, $t5UTR, $t3UTR, $tCompPosNo, $tGeneSnpPosNo, $tIntergSnp, $allSnpNo);
}

# auxiliary sub-routine to print out the SNV position distribution:
# input:	an array containing position counts of SNVs (in the following order): 
#		in Exon, in Intron, in 5' UTR, in 3' UTR, of complex types, in gene regions, in intergenic regions, and, of the total
# output:	printouts of the percentages (over the position count of total SNVs)
sub printSnvDist{
  my ($tEx, $tIn, $t5UTR, $t3UTR, $tCompPosNo, $tGeneSnpPosNo, $tIntergSnp, $allSnpNo) = @_;
  print "SNV Position Distribution:\n";
  if($allSnpNo>0){
    print "Exon: $tEx "; printf "(%.2f%%)\n", ($tEx/$allSnpNo*100); 
    print "Intron: $tIn "; printf "(%.2f%%)\n", ($tIn/$allSnpNo*100); 
    print "5' UTR: $t5UTR ";  printf "(%.2f%%)\n", ($t5UTR/$allSnpNo*100); 
    print "3' UTR: $t3UTR ";  printf "(%.2f%%)\n", ($t3UTR/$allSnpNo*100); 
    print "Complex (in-gene): $tCompPosNo ";  printf "(%.2f%%)\n", ($tCompPosNo/$allSnpNo*100); 
    print "In-gene total: $tGeneSnpPosNo "; printf "(%.2f%%)\n", ($tGeneSnpPosNo/$allSnpNo*100);
    print "Intergenic total: $tIntergSnp ";  printf "(%.2f%%)\n", ($tIntergSnp/$allSnpNo*100); 
  }else{
    print "WARNING: no SNVs are found\n";
  }
  print "Total SNV Positions: $allSnpNo\n";
  print "\n";

}

=head1 NAME

snp_distri.pl -- A simple introductory application script to get familiar with the ASARP pipeline (L<asarp>) 

The application calculates the SNV (position) distribution in exons, introns, etc. according to annotations. A SNV is disarded (not contributing to the total SNV number) if it is covered by < 2 RNA-Seq reads.

=head1 SYNOPSIS

  perl snp_distri.pl output_file snv_file transcript_file [powerful_snv_cutoff]

C<snv_file>: a SNV list, see the format description in L<snpParser>

C<powerful_snv_cutoff>: an optional cutoff to categorize SNVs into powerful (>= C<powerful_snv_cutoff>) and non-powerful types. Default: 20

C<transcript_file>: Transcript and gene annotation file, see the format description in L<fileParser>

=over 6

=item An example file is ../data/hg19.merged.to.ensg.all.tx.03.18.2011.txt, 
which was created by merging ensembl Refseq, UCSC knowngene, Gencode
gene, and Vegagene. 

=back

C<output_file>: Tab dilimited counts of SNV positions in different gene regions. Headers are included as illustrated in the following example:

	Type    Exon    Intron  5'UTR   3'UTR   Complex In-gene Intergenic      Total
	Powerful(>=20)  3788    972     586     4097    1464    7937    342     8279
	Non-powerful(<20)       3846    32143   1246    3679    1323    39566   6580    46146
	Overall 7634    33115   1832    7776    2787    47503   6922    54425

=head1 DESCRIPTION

The application is SNV position (ref genome location) oriented. In other words, if (rarely, and not seen in our data) one position contains multiple SNVs, it will be still considered as one SNV (position). A SNV may overalp multiple transcripts of multiple genes. The rule to determine the categories of a SNV is as follows:

If the SNV overlaps certain transcript exon blocks, it is considered as in (coding) exons, regardless the times it overlaps.

If the SNV overlaps certain 5'/3' UTRs, it is considered as in 5'/3' UTRs (non-coding), regardless the times it overlaps.

Only when a SNV never overlaps any exons nor 5'/3' UTRs, and its genome location is within certain transcript span, it is considered as in introns, regardless of the times it overlaps. As a result, intron SNVs are exclusive to exons, 5'/3' UTRs.

A SNV can be categorized into multiple categories, denoted as complex. Therefore, complex is the union of SNVs with any combinations of Exon and 5'/3' UTR types. In-gene SNVs are the union of all these categories. Therefore, the sum of exon and 5'+3' UTR SNV counts will be larger than the total in-gene SNV count if complex SNVs exist. In-gene and Intergenic SNVs are exclusive to each other.

For any of the above cases, a SNV is considered in-gene. If a SNV is not in-gene, it is considered in the intergenic regions.

The application also outputs percentage (over the total SNV positions) sumamries for powerful, non-powerful and overall SNV distributions. Sample output:

	Calculating powerful SNV distribution...
	SNV Position Distribution:
	Exon: 3788 (45.75%)
	Intron: 972 (11.74%)
	5' UTR: 586 (7.08%)
	3' UTR: 4097 (49.49%)
	Complex (in-gene): 1464 (17.68%)
	In-gene total: 7937 (95.87%)
	Intergenic total: 342 (4.13%)
	Total SNV Positions: 8279

	Calculating non-powerful SNV distribution...
	SNV Position Distribution:
	Exon: 3846 (8.33%)
	Intron: 32143 (69.66%)
	5' UTR: 1246 (2.70%)
	3' UTR: 3679 (7.97%)
	Complex (in-gene): 1323 (2.87%)
	In-gene total: 39566 (85.74%)
	Intergenic total: 6580 (14.26%)
	Total SNV Positions: 46146
	...

=head1 SEE ALSO

L<asarp>, L<fileParser>, L<snpParser>, L<MyConstants>

=head1 COPYRIGHT

This pipeline is free software; you can redistribute it and/or modify it given that the related works and authors are cited and acknowledged.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.

=head1 AUTHOR

Cyrus Tak-Ming CHAN

Xiao Lab, Department of Integrative Biology & Physiology, UCLA

=cut
