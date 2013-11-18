#!/usr/bin/perl
use warnings;
use strict;

#BEGIN {push @INC, '.'}
require "fileParser.pl"; # only to use binarySearch
require "readUtilities.pl";
# set autoflush for error and output
select(STDERR);
$| = 1;
select(STDOUT);
$| = 1;

our $INTRVL = 100000; #interval to output processed counts
my $samType = "a standard SAM file with
	11 attributes. You can also use samtools and bedtools; or
	to use procReadsJ.pl on the special jsam files introduced
	by Dr. Jae-Hyung Lee for RNA-editin and allele specific 
	expression (ASE) studies";
#my $samType = "Dr. Jae-Hyung Lee's
#	20-attribute SAM file output format, used in RNA-editing
#	or allele specific expression (ASE) studies";
if(@ARGV < 6){
  printUsage($samType);
}

#####################################################################
# get parameters
my ($chrToCheck, $samFile, $snvFile, $outputSnvs, $outputBedgraph, $pairEnded, $strandFlag, $title, $discardPos) = @ARGV;
($title, $strandFlag) = parseArg($title, $strandFlag);

# masking
my @discard = ();
if(defined($discardPos)){
  print STDERR "WARNING: masking not tested in this version; ignored\n";
  @discard = readPosMask($discardPos);
}

#####################################################################
# read DNA SNV list
print "Processing input SNV list file of $chrToCheck: $snvFile...";
my $dnaSnvList = readSnvList($snvFile, $chrToCheck);
my ($dnaSnvsRef, $dnaSnvsIdxRef) = procDnaSnv(\$dnaSnvList, $chrToCheck);


#####################################################################
# parsing input sam file
#my ($samRef) = readSamFileJ($samFile);
my ($samRef) = readSamFile($samFile);
my ($blocksRef, $snvRef, $blocksRcRef, $snvRcRef) = parseSamReads($chrToCheck, $samRef, $pairEnded, $strandFlag); 

#####################################################################
# bedgraph handling
#
my ($bedRef, $bedIdxRef) = procBedgraph($chrToCheck, $blocksRef);
if($strandFlag){ #need to handle the minus information as well
  my ($bedRcRef, $bedRcIdxRef) = procBedgraph($chrToCheck, $blocksRcRef);
  outputBedgraph($chrToCheck, $outputBedgraph, $title, $strandFlag, $bedRef, $bedIdxRef, $bedRcRef, $bedRcIdxRef);
}else{
  outputBedgraph($chrToCheck, $outputBedgraph, $title, $strandFlag, $bedRef, $bedIdxRef); # non-strand-specific
}
#####################################################################
# get all SNV reads from the string "snv", bedgraph is not used anymore

# quick termination case
my @dnaSnvs_idx = @$dnaSnvsIdxRef;
if(@dnaSnvs_idx == 0){
  #speed up as no need to check the others
  print "Skipped procesing SNVs with RNA reads as there are no DNA SNVs. Finished with 0 SNVs\n";
  exit;
}

# output SNVs
my $ttlSnvs = 0;
open(SP, ">", $outputSnvs) or die "ERROR: cannot open $outputSnvs to output candidate SNVs\n";
if($strandFlag){ #need to handle the minus information as well
  my ($rnaSnvs, $dCnt) = getSnvReads($chrToCheck, $snvRef, $dnaSnvsRef, $dnaSnvsIdxRef, '+');
  my ($rnaSnvsRc, $dRcCnt) = getSnvReads($chrToCheck, $snvRcRef, $dnaSnvsRef, $dnaSnvsIdxRef, '-');
  print SP $rnaSnvs;
  print SP $rnaSnvsRc;
  $ttlSnvs += $dCnt + $dRcCnt;
}else{
  my ($rnaSnvs, $dCnt) = getSnvReads($chrToCheck, $snvRef, $dnaSnvsRef, $dnaSnvsIdxRef);
  print SP $rnaSnvs;
  $ttlSnvs = $dCnt;
}
close(SP);

print "\nFinished with $ttlSnvs SNVs\n";


########################################################################################################################################
##########################################################################
#		sub-routines (standard SAM file)
##########################################################################

## read sam file and cut most of the useless attributes out
sub readSamFile{

  my ($samFile) = @_;
  my @sam = ();

  my $cnt = 0;
  print "Processing SAM file: $samFile and counting SNV reads...";
  open(FP, "<", $samFile) or die "ERROR: Can't open $samFile";
  my $hdFlag = 0; #header flag
  while(<FP>){
    chomp $_;
    if($_ =~ /^\@/){
      if($hdFlag ==0){
        print "NOTE: Header lines of the sam file will be ignored\n";
      }
      $hdFlag += 1;
      next;
    }
    my @attr = split('\t', $_);
    $sam[$cnt]  = join("\t", $attr[0], $attr[1], $attr[2], $attr[3], $attr[5], $attr[9]); # id strand chr start cigar and read: good enough
    ++$cnt;
    if($cnt%($INTRVL*10) == 0){ print "$cnt...";  }
  }
  close (FP);
  my $unit = 'lines';
  print " $cnt $unit. Done.\n";
  if($hdFlag){
    print "$hdFlag header lines ignored. Done.\n";
  }

  return \@sam;
}

sub parseSamReads
{
  my ($chrToCheck, $samRef, $pairEnded, $strandFlag) = @_;
  print "Processing the reads to get blocks and raw SNV statistics...\n";
  
  my @sam = @$samRef;
  my $N = @sam; # no. of reads
  my $O = 0; #overlap count

  my $snv = ''; # to get the snv candidates
  my $blocks = ''; # all the blocks
  #strand-specific: -
  my $snvRc = '';
  my $blocksRc = '';

  for(my $cnt = 0; $cnt < $N; $cnt ++){
    if($cnt%($INTRVL) == 0){ print "$cnt...";  }
    #if($cnt > $INTRVL/100*4){ last; }
    my @attr = split('\t', $sam[$cnt]); # pair1
    my $strand = getStrandInRead($attr[1], $strandFlag);

    my %snvInPair = (); my $sipRef = \%snvInPair;
    # standard sam format;	alternative: # jsam: JH's sam file  ($sipRef, my $block) = addSnvInPairJ($sipRef, \@attr);   #my $block = $attr[13];
    ($sipRef, my $block) = addSnvInPair($sipRef, \@attr, $dnaSnvsIdxRef, \@discard); 
  
    # paired-end RNA-Seq cases; strand-specific is also handled here
    my @attr2 = ();
    if($pairEnded){ #get pair2
      ++$cnt;
      if(!defined($sam[$cnt])){ # pair2
       die "ERROR: missing pair 2 in $chrToCheck at line $cnt\n";
      }
      my @attr2 = split('\t', $sam[$cnt]); # pair2
      checkPairId($attr[0], $attr2[0], $cnt);
      checkPairChr($attr[2], $attr2[2], $cnt);
      #standard sam format;	alternative: #jsam: #($sipRef, my $block2) = addSnvInPairJ($sipRef, \@attr2);   #my $block2 = $attr2[13];
      ($sipRef, my $block2) = addSnvInPair($sipRef, \@attr2, $dnaSnvsIdxRef, \@discard);
   
      # if pair2 is sense ($strandFlag == 2) and strand eq '+', pair1 position is larger than pair2!
      # if pair1 is sense ($strandFlag == 1) and strand eq '-', pair1 position is larger than pair2!
      my $isOverlap = 0;
      if(($strandFlag == 2 && $strand eq '+') || ($strandFlag == 1 && $strand eq '-') || ($strandFlag ==0 && $strand eq '-')){
         #non-strand specific: - here means antisense
         ($block, $isOverlap) = mergeBlockInPair($block2, $block, $attr[3]); #p2 p1
      }else{ # incl. $strandFlag == 0
         ($block, $isOverlap) = mergeBlockInPair($block, $block2, $attr2[3]); # p1 p2
      }
      $O += $isOverlap;
    }
    #flaten the snvs in pair
    my $snvToAdd = "";
    my %sip = %$sipRef;
    for(keys %sip){
      $snvToAdd .= "$_\t$sip{$_}\n"; # location and allele
    }

    # get all the (masked) blocks and snvs for future processing
    if(!$strandFlag || $strand eq '+'){ #non-strand specific or +; following pair1
      $snv .= $snvToAdd;
      if($blocks eq ''){
        $blocks = $block;
      }else{
        $blocks .= ",".$block;
      }
    }else{ # - strand
      $snvRc .= $snvToAdd; #strand-specific: -; following pair1!
      if($blocksRc eq ''){
        $blocksRc = $block;
      }
      else{
        $blocksRc .= ",".$block;
      }
    }
  }

  print "Done.\n";
  print "OVERLAP_READS\t$O\t".sprintf("%d\t%.3f\n", $N/2, $O/($N/2));
  #print "$snv\n$snvRc\n";

  return (\$blocks, \$snv, \$blocksRc, \$snvRc);
}


sub getSnvReads{
  my ($chr, $snvRef, $dnaSnvsRef, $dnaSnvsIdxRef, $strand) = @_;
  my $snv = $$snvRef; #just a string
  my $dCnt = 0;
  my $rnaSnvs = ""; # final output with reference allele, alt allele and wrong nt
  if(defined($strand)){
    $strand = " $strand"; # for the output format
  }else{ $strand = "";	} # empty string
  print "Processing".$strand." SNVs to get RNA read counts...\n";
  my %rawSnvsA = ();
  my %rawSnvsC = ();
  my %rawSnvsG = ();
  my %rawSnvsT = ();
  my %rawSnvsN = (); # N also appears

  my %dnaSnvs = %$dnaSnvsRef;
  my @dnaSnvs_idx = @$dnaSnvsIdxRef;
  
  my @snvs = split(/\n/, $snv); # flat snv string into array
  for(@snvs){
    # get the raw statistics
    if($_ ne ''){
      my ($loc, $al) = split(/\t/, $_);
      if($al eq 'A'){
        $rawSnvsA{$loc} += 1;
      }elsif($al eq 'C'){
        $rawSnvsC{$loc} += 1;
      }elsif($al eq 'G'){
        $rawSnvsG{$loc} += 1;
      }elsif($al eq 'T'){
        $rawSnvsT{$loc} += 1;
      }elsif($al eq 'N'){
        $rawSnvsN{$loc} += 1;
      }else{
        print STDERR "WARNING: unrecognized allele: $al at $loc\n";
      }
    }
  }

  my %rawSnvs = (
  'A' => \%rawSnvsA,
  'C' => \%rawSnvsC,
  'G' => \%rawSnvsG,
  'T' => \%rawSnvsT,
  'N' => \%rawSnvsN,
  );
  my @alpha = qw(A C G T N);
  for my $loc (@dnaSnvs_idx){
    # get the reference and alternative for further categorization
    #print "$dnaSnvs{$loc}\n";
    my ($als, $id) = split(" ", $dnaSnvs{$loc});
    my ($ref, $alt) = split('>', $als);
    my $refCnt = getRawSnvAl(\%rawSnvs, $loc, $ref); 
    my $altCnt = getRawSnvAl(\%rawSnvs, $loc, $alt);
    if($refCnt || $altCnt){ # either of them should be != 0
      my $wrgCnt = 0;
      for(@alpha){
        if($_ ne $ref && $_ ne $alt){
          $wrgCnt += getRawSnvAl(\%rawSnvs, $loc, $_);
	}
      }
      # only when there is no wt, or wt < alt
      if(!$wrgCnt || $wrgCnt < $altCnt){
        $rnaSnvs .= "$chr $loc $dnaSnvs{$loc} $refCnt:$altCnt:$wrgCnt".$strand."\n"; # strand automatically handled
	#print "$rnaSnvs\n";
	$dCnt += 1;
      }else{
	print "WARNING: DISCARD: $chrToCheck:$loc alt $alt rna count <= other nts: $altCnt <= $wrgCnt\n"; 
      }
    }
    
  }
  return ($rnaSnvs, $dCnt);
}

sub getRawSnvAl
{
  my ($ref, $loc, $al) = @_;
  my %hs = %$ref;
  if(defined($hs{$al})){
    my %als = %{$hs{$al}};
    if(defined($als{$loc})){
      return $als{$loc};
    }
  }
  return 0;
}

# this is for the general SAM files (not using extra attributes)
# add SNV (attributes of the original SAM file) to the read pair list
sub addSnvInPair{
  my ($ref, $attr, $dnaSnvsIdx, $discardRef) = @_;
  my %snv = %$ref;
  #my %dnaSnvs = %$dnaSnvRef;
  my $n = @{$dnaSnvsIdx};
  my @discard = ();
  if(defined($discardRef)){	 @discard = @$discardRef;	}
  my $start = $$attr[3];
  my $cigar = $$attr[4];
  my $block = parseCigar($start, $cigar);
  if(!defined($block)){ # not supported CIGAR, directly return
    return ($ref, '');
  }
  #print "$cigar, $block | $$attr[11]\n";
  my $snvToAdd = "";

  my @blocks = split(',', $block);
  my $rePos = 0; # relative position
  for(@blocks){
    my ($s, $e) = split(':', $_);
#=pod    
    my ($lBound, $lUnMatch) = binarySearch($dnaSnvsIdx, $s, 0, $n-1, 'left');
    for(my $j = $lBound; $j < $n && $$dnaSnvsIdx[$j]<= $e; $j++){ #snv is more sparse!
      my $pos = $$dnaSnvsIdx[$j];
      my $readPos = $rePos + $pos-$s; # 0-based
      if(!defined($discard[$readPos]) || !$discard[$readPos]){ # not discarded
	  my $al = substr($$attr[5], $readPos, 1);
	  if(defined($snv{$pos})){ # existing
	    if($snv{$pos} ne $al){
	      print "WARNING: INCONSISTENT SNV at $pos: $al diff from $snv{$pos}\n";
	      delete $snv{$pos};
	    }
	  }else{
	    $snv{$pos} = $al; #\t$qual\n"; # add new
	  }
      }

    }
    $rePos += $e-$s+1; # next block
#=cut

  }
  #handling of discarded position
  if(@discard > 0){ # need to handle discarded positions
      my $readLen = length($$attr[5]);
      $block = maskBlock($block, $readLen, $discardRef);
  }
  return (\%snv, $block);
}

__END__

=head1 NAME

procReads.pl -- Processing a duplicate-removed SAM file (L<rmDup>) of a chromosome to generate the chromosome specific SNV list and the bedgraph file. The output files are used as input files for the ASARP pipeline.

The new procReads.pl supports strand-specific paired-end RNA-Seq data (i.e. the SAM file strand information is reliable) for more accurate results.

=head1 SYNOPSIS

This is part of the full pre-processing:

=over 6

1. rmDup (removing PCR duplicates for SAM files (including Dr. JH Lee's SAM format); samtools/bedtools can be used for standard SAM files)

2. mergeSam (merging SAM files if there are independent duplicates)

3. B<procReads> (processing SAM files to get SNV read counts and generate bedgraph files) 

=back

USAGE:

 perl procReads.pl input_sam_file input_snvs output_snvs output_bedgraph is_paired_end

NOTE:

the read processing script is for a standard SAM file with
11 attributes. You can also use samtools and bedtools; or
to use procReadsJ.pl on the special jsam files introduced
by Dr. Jae-Hyung Lee for RNA-editin and allele specific 
expression (ASE) studies

There are some assumptions and requirements for the input SAM/JSAM files.
See L<Files> for more details.

ARGUMENTS:

 chr			chromosome to be investigated (correspond to the input_sam_file)
 input_sam_file		SAM file input after duplicate removal (use rmDup.pl)
 intput_snvs		input SNV list (without read counts)
 output_snvs		output SNV candidates with read counts
 output_bedgraph	output bedgraph file, see below for the details:
			http://genome.ucsc.edu/goldenPath/help/bedgraph.html
 is_paired_end		0: single-end; 1: paired-end
			For paired-end reads, all reads should be paired up, 
			where pair-1 should be always followed by pair-2 in the next line.


OPTIONAL [strongly recommended to be input]:

 is_strand_sp		0: non-strand specific (or no input); 
			1: strand-specific with pair 1 **sense**;
			2: strand-specific with pair 1 **anti-sense**
			Be careful with the strand-specific setting as it will give totally opposite 
			strand information if wrongly set.

			The strand-specific option is used for strand-specific RNA-Seq data.
			When set, specialized bedgraph files will be output (output_bedgraph)
			where there is a 5th extra attribute specifying the strand: + or -
			besides the standard ones: http://genome.ucsc.edu/goldenPath/help/bedgraph.html
			One can use grep and cut to get +/- strand only bedgraphs.

OPTIONAL [if input, must be input in order following is_strand_sp]:

 bedgraph_title		a short title for the output bedgraph files (will be put in description of the header line)
                        if there are spaces in between it should be quoted
			e.g. "nbt.editing reads: distinct after dup removal"
			if not input, "default" will be used as the short title


=head1 DESCRIPTION

=head2 INPUT

C<input_sam_file> should contain only 1 chromosome, and it should be in SAM/JSAM format. If it is in JSAM format (check out http://www.ncbi.nlm.nih.gov/pubmed/21960545 for more details), L<procReadsJ> can also be used.

The SNP list should be in a format like this:

 chr1 20129 C>T rs12354148 2:0:0
 chr1 118617 T>C na 1:0:0
 chr1 237763 G>A rs79665216 1:0:0
 chr1 565508 G>A rs9283150 0:1:0

Each line is space separated, with

 chromosome
 location 
 ref_allele>alt_allele 
 dbSnp_id 
 [read_counts ref:alt:others]

Only the first 4 fields will be parsed so
C<read_counts> are not needed and ignored which are from DNA genomic sequencing.

To avoid unnecessary computational time to read SNVs of other chrosomes than the one in C<inptu_sam_file>, it is suggested to keep SNVs of the same chromosome in one seperate file.

=head2 OUTPUT

By default, the strand specific flag C<is_strand_sp> is unset, and the C<inptu_sam_file> strand information is considered unreliable and not used. The output files are standard bedgraph files http://genome.ucsc.edu/goldenPath/help/bedgraph.html with space as the dilimiter, and SNV files should follow the file format described in L<Files>:

Example SNV output lists look like this:

 chr10 1046712 G>A rs2306409 50:39:0
 chr10 1054444 A>G rs11253567 2:0:0
 chr10 1055866 A>T rs4880751 1:2:1
 chr10 1055949 G>A rs12355506 7:2:0
 chr10 1055968 G>A rs72478237 6:5:0
 chr10 1060218 G>A rs3207775 42:37:0

When C<is_strand_sp> is set, the program outputs bedgraph and SNV files with the extra last strand attributes. Output SNV examples would look like:

 chr10 1046712 G>A rs2306409 30:23:0 +
 chr10 1055866 A>T rs4880751 1:2:0 +
 chr10 1055949 G>A rs12355506 2:0:0 +
 chr10 1055968 G>A rs72478237 2:2:0 +
 chr10 1060218 G>A rs3207775 27:22:0 +
 ...
 chr10 1046712 G>A rs2306409 20:16:0 -
 chr10 1054444 A>G rs11253567 2:0:0 -
 chr10 1055949 G>A rs12355506 5:2:0 -
 chr10 1055968 G>A rs72478237 4:3:0 -
 chr10 1060218 G>A rs3207775 15:15:0 -

As a result, one SNV may appear twice if it has valid ref:alt:wrnt read counts on both + and - strands. In the example above, one can have more accurate information for the RNA read counts, especially when there are genes on the opposite strands.

The output bedgraph lines would look like:

 chr10 181481 181482 7 +
 chr10 181482 181483 9 +
 chr10 181483 181499 10 +
 ...
 chr10 181479 181482 5 -
 chr10 181482 181483 8 -
 chr10 181483 181499 9 -
 ...

Note that all + lines are output before any - lines output. While the 5th strand attribute is not specified in the bedgraph standard, 
two additional bedgraph files are output, with suffixes .plus and .minus, to provide the + and - only standard bedgraph tracks respectively.

=head1 SEE ALSO

L<rmDup>, L<mergeSam>, L<asarp>

=head1 COPYRIGHT

This pipeline is free software; you can redistribute it and/or modify it given that the related works and authors are cited and acknowledged.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.

=head1 AUTHOR

Cyrus Tak-Ming CHAN

Xiao Lab, Department of Integrative Biology & Physiology, UCLA

=cut
