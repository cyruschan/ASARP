#!/usr/bin/perl
use warnings;
use strict;
# a simple demonstration with the ASARP pipeline:
# in this application, the transcripts and SNVs are read in and the SNV distribution in 5'/3' UTRs, exons and introns is reported

# NOTE: one SNV can match several transcripts

use MyConstants qw( $CHRNUM $supportedList $supportedTags );
require "fileParser.pl"; #sub's for input annotation files
require "snpParser.pl"; #sub's for snps

# input arguments: $outputFile--output, $snpF--SNV file path, $xiaoF--transcript annotation file
# optional: $POWCUTOFF--powerful SNV cutoff, default: 20
my ($outputFile, $snpF, $xiaoF, $strandFlag, $POWCUTOFF, $snvPwrOut, $snvOrdOut, $isDetailed) = @ARGV;
if(@ARGV < 4){ # !defined($outputFile) || !defined($snpF) || !defined($xiaoF)){
  print <<EOT; 
  
USAGE: perl $0 output_file snv_file transcript_file strand_flag [powerful_snv_cutoff pwr_snv_details ordinary_snv_details is_detailed]

ARGUMENTS:
strand_flag	0--non-strand-specific; 1/2--strand-specific (not distinguishable at this step).
		If set, there will be 2 snv_details files for each category (pwr or ordinary), 
		with .plus and .minus indicating the distribtuion in the + and - strands
OPTIONAL:
The optional arguments must be input in order.
pwr_snv_details	
ordinary_snv_details 
		These two are the output files for the detailed SNV categories of 
		powerful and non-powerful (ordinary) SNVs respectively. 
is_detailed	If is_detailed is set 1 after both *_snv_details, detailed information of the 
		SNVs will be output: chr;geneName;snvPos;type
EOT
  exit;
}
if($strandFlag == 1 || $strandFlag == 2){
  print "NOTE: strand-specific output is enabled. There will be *.plus and *.minus for each snv detail output file if it is specified\n";
}elsif($strandFlag == 0){
  print "NOTE: Non-strand-specific output is set.\n"; 
}else{
  die "ERROR: strand_flag is expected to be 0 or 1/2. But $strandFlag is input\n";
}

if(!defined($POWCUTOFF)){
  $POWCUTOFF = 20;
  print "Default powerful_snv_cutoff used: $POWCUTOFF\n";
}

my $fp = undef;
open($fp, ">", $outputFile) or die "ERROR: Cannot open $outputFile to write\n";

# Demonstration of using the fileParser
# read the transcript annotation file
my ($transRef) = readTranscriptFile($xiaoF);
my $transGeneRef = undef;
if(defined($isDetailed) && $isDetailed eq "1"){ #one is the only setting
  print "NOTE: detailed SNV information (gene, type, transcript ID, neighboring exons) will be output as pwr_snv_details.lst ordinary_snv_details.lst \n";
  ($transGeneRef) = readTranscriptFileByGene($xiaoF);
}
#printListByKey($transRef, 'trans'); #utility sub: show transcripts (key: trans)

# need to make it strand specific:
if($strandFlag){

  my ($snpRef, $pRef) = initSnp($snpF, $POWCUTOFF, '+');
  my ($snpRcRef, $pRcRef) = initSnp($snpF, $POWCUTOFF, '-');
  my ($pwrOut, $ordOut, $pwrOutRc, $ordOutRc) = (undef, undef, undef, undef);
  if(defined($snvPwrOut)){
    $pwrOut = $snvPwrOut.".plus";
    $pwrOutRc = $snvPwrOut.".minus";
    print "Two powerful SNV detail files: $pwrOut and $pwrOutRc will be output for plus and minus strands\n";
  }
  if(defined($snvOrdOut)){
    $ordOut = $snvOrdOut.".plus";
    $ordOutRc = $snvOrdOut.".minus";
    print "Two ordinary SNV detail files: $ordOut and $ordOutRc will be output for plus and minus strands\n";
  }
  print "+ strand only SNP distribution:\n";
  print $fp "Type (+ only)\tExon\tIntron\t5'UTR\t3'UTR\tComplex\tIn-gene\tIntergenic\tTotal\n";
  my @total = snpDistOneStrand($snpRef, $transRef, $fp, '+', $pwrOut, $ordOut, $transGeneRef);
  print "- strand only SNP distribution:\n";
  print $fp "Type (- only)\tExon\tIntron\t5'UTR\t3'UTR\tComplex\tIn-gene\tIntergenic\tTotal\n";
  my @totalRc = snpDistOneStrand($snpRcRef, $transRef, $fp, '-', $pwrOutRc, $ordOutRc, $transGeneRef);

  print "\nOverall SNV distribution (+ strand only):\n";
  printSnvDist(@total);
  print $fp "Overall (+)\t".join("\t",@total)."\n";
  
  print "\nOverall SNV distribution (- strand only):\n";
  printSnvDist(@totalRc);
  print $fp "Overall (-)\t".join("\t",@totalRc)."\n";
  
  close($fp);
}else{
  #non-strand-specific version
  
  print $fp "Type\tExon\tIntron\t5'UTR\t3'UTR\tComplex\tIn-gene\tIntergenic\tTotal\n";
  # Demonstration of using the snpParser
  my ($snpRef) = initSnp($snpF, $POWCUTOFF);
  #print "SNV List:\n";
  #printListByKey($snpRef, 'powSnps');
  #printListByKey($snpRef, 'snps');
  my @total = snpDistOneStrand($snpRef, $transRef, $fp, undef, $snvPwrOut, $snvOrdOut, $transGeneRef);
  #the elements contained in the total array:
  #($tEx, $tIn, $t5UTR, $t3UTR, $tCompPosNo, $tGeneSnpPosNo, $tIntergSnp, $allSnpNo) = @total;
  print "\nOverall SNV distribution:\n";
  printSnvDist(@total);
  print $fp "Overall\t".join("\t",@total)."\n";
  close($fp);
}


sub snpDistOneStrand{
  # $fp is the file handle being open, bad (need to finish it quick though)
  my ($snpRef, $transRef, $fp, $strandInfo, $snvPwrOut, $snvOrdOut, $transGeneRef) = @_;

  my ($geneSnpRef) = setGeneSnps($snpRef, $transRef, $strandInfo);
  #print "Significant Snvs: \n";
  #printGetGeneSnpsResults($geneSnpRef,'gPowSnps', $snpRef,'powSnps', 1); #$SNVPCUTOFF);
  #print "Ordinary Snvs: \n";
  #printGetGeneSnpsResults($geneSnpRef,'gSnps', $snpRef,'snps', 1);

  print "Calculating powerful SNV distribution...\n";
  my @part1 = getGeneSnpsDistri($geneSnpRef,'gPowSnps', $snpRef,'powSnps', $strandInfo, $snvPwrOut, $transGeneRef); 
  print $fp "Powerful(>=$POWCUTOFF)\t".join("\t",@part1)."\n";
  printSnvDist(@part1);

  print "Calculating non-powerful SNV distribution...\n";
  my @part2 = getGeneSnpsDistri($geneSnpRef,'gSnps', $snpRef,'snps', $strandInfo, $snvOrdOut, $transGeneRef); 
  print $fp "Non-powerful(<$POWCUTOFF)\t".join("\t",@part2)."\n";
  printSnvDist(@part2);

  my @total = ();
  for(my $i = 0; $i < @part1; $i++){
    $total[$i] = $part1[$i] + $part2[$i];
  }
  return @total;
}

# print out the gene snps results for certain type (powerful or ordinary snps)
sub getGeneSnpsDistri
{
  my ($geneSnpRef, $geneSnpKey, $snpRef, $snpKey, $strandInfo, $detailedOutput, $transGeneRef) = @_;
  my $snvDetails = ""; #super-detailed SNV information

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
      my $chr = formatChr($i); #print "\n";
      #print "@genes\t";
      foreach(@genes){
        my @allGenes = split('\t', $_);
        foreach(@allGenes){
	  #get the snp information
	  my @geneMatches = split('\t', $geneRef->{$_});
	  for(@geneMatches){
	    my ($snpPos, $matchInfo, $geneName, $txStart, $id, $regStart, $regEnd) = split(';', $_);
	    
	    # To get even more details of the SNVs
	    if(defined($transGeneRef)){ # super-detailed information 
	      my $key = join(";", $chr, $geneName, $snpPos, uc $matchInfo);
	      if(!defined($genesPrinted{$key})){
	        #$snvDetails .= "$chr\t$geneName\n";
	        $genesPrinted{$key} = 1;
		$snvDetails .= "$key\n";
	      }else{ $genesPrinted{$key} += 1; }
	      #print $snps{$snpPos}."\n";
              #my ($prevRef, $nextRef) = findNextExons($chr, $snpPos, $geneName, $transGeneRef, $strandInfo); 
	      #my $prevE = join(";", keys %$prevRef); 
	      #my $nextE = join(";", keys %$nextRef); 
	      #my @allSnpInfo = split(';', $snps{$snpPos}); #separate by ;, if there are multiple snps at the same position
              #foreach(@allSnpInfo){
	      #  my $toPrint = $matchInfo."\t".$geneName."\t".$id.":".$regStart."-".$regEnd."\t$prevE\t$nextE\n";
	      #  $snvDetails .=  $snpPos."\t".$toPrint;
	      #}
	    }
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

  #add detailed output file
  my $fp = undef;
  if(defined($detailedOutput)){
    if(!defined($strandInfo)){
      open($fp, ">", $detailedOutput) or die "ERROR: Cannot open $detailedOutput\n";
      print $fp "chr\tpos\tcategory\n";
    }else{
      open($fp, ">", $detailedOutput) or die "ERROR: Cannot open $detailedOutput\n";
      print $fp "chr\tpos\tcategory\tstrand\n"; #additional strand
    }
  }

  if(defined($detailedOutput) && defined($transGeneRef)){
    open(my $dp, ">", "$detailedOutput.lst") or die "ERROR: Cannot open $detailedOutput for SNV details\n";
    print $dp $snvDetails;
    close($dp);
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

      if(defined($detailedOutput)){
        my $category = "";
	if($inEx+$in5UTR+$in3UTR == 0){
	  $category = "INTRON;";
	  #if($inIn){
	  #  $category = "INTRON";
	  #}else{
	  #  $category = "INTERGENIC"; #no such cases in in-gene SNVs
	  #}
	}else{
	  if($inEx){ $category .="EXON;"; }
	  if($in5UTR){ $category .="5'UTR;"; }
	  if($in3UTR){ $category .="3'UTR;"; }
	}
        if(!defined($strandInfo)){
	  print $fp "$chr\t$_\t$category\n";
	}else{
	  print $fp "$chr\t$_\t$category\t$strandInfo\n";
	}
      }
    }
  }
  close($fp) if defined($detailedOutput);

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


# find the next exons given an Intronic position
sub findNextExons{
  my ($chr, $pos, $gene, $transGeneRef, $strand) = @_;
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

=head1 NAME

snp_distri.pl -- To calculate the SNV (position) distribution in exons, UTRs, introns, etc., according to transcript annotations. A SNV is disarded (not contributing to the total SNV number) if it is covered by < 2 RNA-Seq reads.

It also serves as an introductory application script to get familiar with the ASARP pipeline (L<asarp>) 

=head1 SYNOPSIS

  perl snp_distri.pl output_file snv_file transcript_file [powerful_snv_cutoff pwr_snv_details ordinary_snv_details [is_detailed]]

C<snv_file>: a SNV list, see the format description in L<snpParser>

C<powerful_snv_cutoff>: an optional cutoff to categorize SNVs into powerful (>= C<powerful_snv_cutoff>) and non-powerful types. Default: 20

C<transcript_file>: Transcript and gene annotation file, see the format description in L<fileParser>

The optional arguments must be input in order.

C<pwr_snv_details> and C<ordinary_snv_details> are the output files for the detailed SNV categories of powerful and non-powerful (ordinary) SNVs respectively.

C<is_detailed>:	default: 0; when it is set 1 after both *_snv_details, detailed information of the SNVs will be output.
		The detail format: chr;geneName;snvPos;type

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

Only when a SNV never overlaps any exons nor 5'/3' UTRs, and its genome location is within certain transcript span, it is considered as in introns, regardless of the times it overlaps. As a result, intron SNVs are exclusive to exons and 5'/3' UTRs.

A SNV can be categorized into multiple categories, denoted as complex. Therefore, complex is the union of SNVs with any combinations of Exon and 5'/3' UTR types. In-gene SNVs are the union of all the categories: exon, 5' UTR, 3' UTR, intron. Therefore, the sum of intron, exon, 5' and 3' UTR SNV counts will be larger than the total in-gene SNV count if complex SNVs exist. In-gene and Intergenic SNVs are exclusive to each other.

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


If C<pwr_snv_details> and C<ordinary_snv_details> are input, you will have two addtional output files providing detailed categories of all individual B<in-gene> SNVs. E.g.

 chr	pos	category
 chr1	68586620	INTRON
 chr1	68591173	3'UTR;
 chr1	68590177	INTRON
 chr1	68608003	INTRON
 chr1	68591253	3'UTR;
 chr1	68591405	3'UTR;
 chr1	68624878	EXON;
 chr1	68589935	INTRON
 chr1	68585708	INTRON

=head1 SEE ALSO

L<asarp>, L<fileParser>, L<snpParser>, L<MyConstants>

=head1 COPYRIGHT

This pipeline is free software; you can redistribute it and/or modify it given that the related works and authors are cited and acknowledged.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.

=head1 AUTHOR

Cyrus Tak-Ming CHAN

Xiao Lab, Department of Integrative Biology & Physiology, UCLA

=cut
