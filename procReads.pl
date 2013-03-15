#!/usr/bin/perl
use warnings;
use strict;

require "readUtilities.pl";
# set autoflush for error and output
select(STDERR);
$| = 1;
select(STDOUT);
$| = 1;

if(@ARGV < 5){
  print <<EOT;
USAGE: perl $0 input_sam_file input_snvs output_snvs output_bedgraph is_paired_end [discarded_read_pos]

NOTE: 	the read processing script is for Dr. Jae-Hyung Lee's
	20-attribute SAM file output format, used in RNA-editing
	or allele specific expression (ASE) studies

input_sam_file 		SAM file input after duplicate removal (use rmDup.pl)
intput_snvs		input SNV list (without read counts)
output_snvs		output SNV candidates with read counts
output_bedgraph		output bedgraph file, see below for the details:
			http://genome.ucsc.edu/goldenPath/help/bedgraph.html
is_paired_end		0: single-end; 1: paired-end
			For paired-end reads, all reads should be paired up, 
			where pair-1 should be always followed by pair-2 in the next line.

OPTIONAL:
discarded_read_pos	masked-out (low-quality) read positions in calculating 
			the max read quality scores, 
			in 1-based, inclusive, interval (a:b,c:d,... no space) format:
			e.g. 1:1,61:70 will discard the 1st, 61st-70th read positions.
			NOTE: the remaining reads will still contain the positions.

EOT
  exit;
}

my ($samFile, $snvFile, $outputSnvs, $outputBedgraph, $pairEnded, $discardPos) = @ARGV;

our $INTRVL = 100000; #interval to output processed counts

my @discard = ();
if(defined($discardPos)){
  @discard = readPosMask($discardPos);
}

print "Processing SAM file: $samFile...";
open(FP, "<", $samFile) or die "ERROR: Can't open $samFile";
my %snv = (); # to get the snv candidates
my %blocks = (); # all the blocks
my %pileups = (); # to get the (pre-) pileups of the pairs

my $cnt = 0;
while(<FP>){
  $cnt++;
  if($cnt%$INTRVL == 0){ print "$cnt...";  }
  my $pair1 = $_;
  chomp $pair1;
  my @attr = split('\t', $pair1);
  # we need no. of matches (attr 14, incl. known snps) to further process
  my ($id, $strand, $chr, $start, $block, $mismatches) = ($attr[0], $attr[1], $attr[2], $attr[3], $attr[11], $attr[13]); #read quality rather than mapping quality is used
  
  if($mismatches > 0){ #we can record SNV candidates
    # attr 16-20
    $snv{$chr} .= "$attr[15]\t$attr[16]\t$attr[17]\t$attr[18]\t$attr[19]\n";
  }
  if(@discard > 0){
    my $readLen = length($attr[9]); #read length
    $block = maskBlock($block, $readLen, \@discard); 
  }
  if($pairEnded){ #get pair2
    my $pair2 = <FP>;
    if(defined($pair2)){
      chomp $pair2;
      my @attr2 = split('\t', $pair2);
      my ($id2, $strand2, $chr2, $start2, $block2, $mismatches2) = ($attr2[0], $attr2[1], $attr2[2], $attr2[3], $attr2[11], $attr2[13]);
      if($chr ne $chr2){  #not checking ids #($id ne $id2 || $chr ne $chr2){
	print "ERROR: reads have to be in pairs in the same chromosome: $id in $chr different from $id2 in $chr2 at line $cnt\n";
	#print "ERROR: reads have to be in pairs: $id in $chr different from $id2 in $chr2 at line $cnt\n";
	exit;
      }
      if($mismatches2 > 0){
        $snv{$chr} .= "$attr2[15]\t$attr2[16]\t$attr2[17]\t$attr2[18]\t$attr2[19]\n";
      }
      if(@discard > 0){
        my $readLen2 = length($attr[9]);
        $block .= ",".maskBlock($block2, $readLen2, \@discard);
      }else{
        $block .= ",".$block2;
      }
    }else{
      print "ERROR: missing pair 2 in $chr at line $cnt\n"; exit;
    }
  }

  # get all the (masked) blocks
  if(!defined($blocks{$chr})){
    $blocks{$chr} = $block;
  }
  else{
    $blocks{$chr} .= ",".$block;
  }
}
close (FP);
my $unit = 'lines';
if($pairEnded){	$unit = 'pairs'; }
print " $cnt $unit. Done.\n";

my @allChrsInSam = keys %blocks;
if(@allChrsInSam > 1){
  die "ERROR: the SAM file should contain reads from only one chromosome: @allChrsInSam\n";
}

print "Processing input SNV list file: $snvFile...";
my %snvList = ();
my $dSnvCnt = 0;
my $dSnvChrCnt = 0;
open(SNP, "<", $snvFile) or die "ERROR: Can't open $samFile";
while(<SNP>){
  chomp;
  $dSnvCnt++;
  if($dSnvCnt%$INTRVL == 0){ print "$dSnvCnt...";  }
  my ($chr, $pos, $alleles, $id) = split(' ', $_);
  if($chr ne $allChrsInSam[0]){
    next; #ignore all other chromosomes to save time
  }
  my $str = join(" ", $pos, $alleles, $id);
  if(defined($snvList{$chr})){
    $snvList{$chr} .= "\t".$str;
  }else{
    $snvList{$chr} .= $str;
  }
  $dSnvChrCnt++;
}
close(SNP);
print " $dSnvCnt SNVs ($dSnvChrCnt kept for $allChrsInSam[0]). Done.\n";
#print "SNVs in chromosomes\n";
#for(keys %snvList){
#  my @chrSnvs = split('\t', $snvList{$_});
#  print "$_: ".(scalar @chrSnvs)."\n";
#}

#############################################################
# output SNVs
open(SP, ">", $outputSnvs) or die "ERROR: cannot open $outputSnvs to output candidate SNVs\n";
# output bedgraph
open(BP, ">", $outputBedgraph) or die "ERROR: cannot open $outputBedgraph to output bedgraph\n";
my $ttlSnvs = 0;
# sort all the block starts and make pile-ups
for my $chr (keys %blocks){
  print "Processing $chr mapped blocks...";
  my ($blkRef, $blkIdxRef) = procBlocks($blocks{$chr});
  my %thisChrBlks = %$blkRef;
  my @blkStarts = @$blkIdxRef;
  my $blkNo = @blkStarts;
  print " $blkNo block starts. Done.\n";

  print "Processing $chr bedgraph...";
  my @bedgraph = (); #bedgraph content in increasing order
  my @bedgraph_idx = (); #sorted index array of the bedgraph

  my $left4Next = ""; #those parts running into the next block
  for(my $i=0; $i<$blkNo; $i++){
    my $thisStart = $blkStarts[$i];
    my $nextStart = -1;
    if($i<$blkNo-1){ #if not the last one, there still a next start
      $nextStart = $blkStarts[$i+1]; 
    }
    # get the previous left-overs
    if($left4Next ne ""){
      $thisChrBlks{$thisStart} .= "\t".$left4Next;
      # $left4Next can be cleared up to here
      $left4Next = "";
    }
    my @thisEnds = split("\t", $thisChrBlks{$thisStart});
    @thisEnds = sort{$a<=>$b} @thisEnds; #sort all the ends

    my $pileupNo = @thisEnds;
    my $s = $thisStart-1; #make it zero-based to be consistent with bedgraph
    for(my $j = 0; $j < $pileupNo; $j++){
      #$nextStart == -1 means the last block, and i need to handle it without left-overs
      if($nextStart != -1 && $thisEnds[$j] >= $nextStart){ # i won't do it
        # dump left-overs to the next guy whenever possible
	if($left4Next eq ""){
	  $left4Next = $thisEnds[$j];
	}else{
          $left4Next .= "\t".$thisEnds[$j];
	}
	$thisEnds[$j] = $nextStart-1; # i will handle up to here
      } 

      #ok, it's my responsibiity to handle it
      if($s == $thisEnds[$j]){
        next; #i.e. there are two blks with the same ends, and it's already counted
      }
      #if the first pileup is an extension to the previous one: prv end = this start, counts equal
      if($j == 0 && @bedgraph>0){
        #extendable from the previous
	my ($pChr, $pS, $pE, $pCnt) = split(' ', $bedgraph[-1]); #chr start end count
	if($pE == $s && $pCnt == $pileupNo){
	  pop @bedgraph; #need to update the previous one
	  pop @bedgraph_idx; #idex needs the pop as well
	  $s = $pS; #extend the start (use the previous one)
	}
      }
      #the start needs to be zero-based
      push(@bedgraph, "$chr $s $thisEnds[$j] ".($pileupNo-$j));
      push(@bedgraph_idx, $s); #the start as the index
      $s = $thisEnds[$j]; # the next start converted to zero-based: +1 then -1
    }
  }
 
  ####################################################################
  # output bedgraph for this chromosome
  # get the track opt first
  my ($dummyChr, $dummyS, $bedLastPos) = split(" ", $bedgraph[-1]);
  my $trackRange = "$chr:".($bedgraph_idx[0]+1).":".$bedLastPos; #start converted back to 1-based for chr range in UCSC
  my $trackOpt = "track type=bedGraph name=\"reads:distinct in $trackRange\" description=\"nbt.editing reads: distinct after dup removal in $trackRange\" visibility=full autoScale=on gridDefault=on graphType=points yLineOnOff=on yLineMark=0 smoothingWindow=off alwaysZero=on\n";
  print BP $trackOpt;
  my $bedCnt = 0;
  for(@bedgraph){
    print BP $_."\n";
    $bedCnt++;
    if($bedCnt%$INTRVL == 0){ print "$bedCnt...";  }
  }
  print " $bedCnt lines. Done.\n";
  ####################################################################
  
  print "Processing $chr genomic SNVs...";
  my @dnaSnvVals = split('\t', $snvList{$chr});
  my %dnaSnvs = ();
  for(@dnaSnvVals){
    my ($pos, $als, $id) = split(' ', $_);
    $dnaSnvs{$pos} = $als." ".$id;
  }
  my @dnaSnvs_idx = sort{$a <=> $b} keys %dnaSnvs;
  my $dnaSnvNo = keys %dnaSnvs;
  print " $dnaSnvNo SNVs. Done.\n";
  ####################################################################
  
  print "Processing $chr candidate SNVs (mismatches)...";
  # here we can forget about the blocks and focus on the pileups (bedgraph)
  my ($snvRef) = procSnv($snv{$chr}, \@discard); 
  my %snvs = %$snvRef;
  my $rSnvNo = keys %snvs;
  print " $rSnvNo SNV locations. Done.\n";

  print "Processing $chr SNVs with read counts...\n";
  my ($bi, $si) = (0, 0);#b for bedgraph_idx, s for dnaSnv_idx, 
  my $dCnt = 0; #$dCnt for matched dnaSnvs_idx
  while($bi<@bedgraph_idx && $si<@dnaSnvs_idx){
    # print "bed idx $bi: $bedgraph_idx[$bi]\nsnv idx: $si: $dnaSnvs_idx[$si]\n";
    if($bedgraph_idx[$bi]+1 <= $dnaSnvs_idx[$si]){
      $bi++;
      if($bi<@bedgraph_idx){
        next;
      } #if $bi is the last of bedgraph, must match
    }
    
    # no match for this dna snv
    if($bi-1 < 0 || $bedgraph_idx[$bi-1]+1 > $dnaSnvs_idx[$si]){
      print "DISCARD: $chr:$dnaSnvs_idx[$si] matches no RNA-Seq reads\n";
      $si++;
      next;
    }
    
    # a match
    my ($chrBed, $startBed, $endBed, $cntBed) = split(' ', $bedgraph[$bi-1]);
    #print "SNV $si ($dnaSnvs_idx[$si]) matches bedgraph $bi-1: $bedgraph[$bi-1]\n";
    # further filter with RNA SNVs (mismatch counts)
    if(defined($snvs{$dnaSnvs_idx[$si]})){
      my ($dAls, $dId) = split(' ', $dnaSnvs{$dnaSnvs_idx[$si]});
      my ($dRefAl, $dAltAl) = split('>', $dAls);
    
      my %snvHs = %{$snvs{$dnaSnvs_idx[$si]}};
      my $altAl = '';
      my $maxAltCnt = 0;
    
      my $misCnt = 0;
      my $misType = 0;
      for(keys %snvHs){
        if(length($_) == 1){ #only mismatches
          $misType++;
	  $misCnt += $snvHs{$_};
	  if($snvHs{$_} > $maxAltCnt){
	    $maxAltCnt = $snvHs{$_};
	    $altAl = $_;
	  }
        }
      }
      #checking the genomic Snv with the RNA Snv
      if($dRefAl ne $snvHs{"ref"}){
        print "DISCARD: $chr:$dnaSnvs_idx[$si] ref inconsistent: dna: $dRefAl, rna: ".$snvHs{"ref"}."\n";
      }
      else{
        my $refCnt = $cntBed - $misCnt;
	my $altCnt = 0;
	if(defined($snvHs{$dAltAl})){  
	  $altCnt = $snvHs{$dAltAl};  
	}
        my $wrongCnt = $misCnt - $altCnt;
        if($altCnt == $maxAltCnt && $altCnt>$wrongCnt){ 
          my $snvStr = join(" ", $chr, $dnaSnvs_idx[$si], $dAls, $dId, $refCnt.":".$maxAltCnt.":".$wrongCnt);
          print SP $snvStr, "\n";
	  $dCnt++;
        }else{
	  print "DISCARD: $chr:$dnaSnvs_idx[$si] alt $dAltAl rna count <= other alts: $altCnt <= $wrongCnt\n"; 
	}
      }
    }else{
      # fully match, i.e. only the reference SNV appears here
      # $refCnt is just $cntBed here
      my $snvStr = join(" ", $chr, $dnaSnvs_idx[$si], $dnaSnvs{$dnaSnvs_idx[$si]}, $cntBed.":0:0");
      print SP $snvStr, "\n";
      $dCnt++;
    }
    $si++;
    if($dCnt>0 && $dCnt%$INTRVL==0){
      print "$dCnt...\n";
    }
  }
  print " $dCnt matched SNVs. Done.\n";
  $ttlSnvs += $dCnt;
}
close(SP);
close(BP);
print "\nFinished with $cnt read $unit and $ttlSnvs SNVs\n";

##########################################################################
#		sub-routines
##########################################################################

# simply process the block (which may contain sub-blocks) into the blocks hash
# input:	the block to be merged
# output:	reference to the hash of blocks
# output:	reference to the index array of the block starts
sub procBlocks{
  my ($block) = @_; 
  my %allBlocks = ();
  my @blks = split(',', $block);
  for(@blks){
    my ($blkS, $blkE) = split(':', $_);
    if(defined($allBlocks{$blkS})){
      $allBlocks{$blkS} .= "\t".$blkE;
    }else{
      $allBlocks{$blkS} = $blkE;
    }
  }
  my @blkStarts = sort{$a<=>$b} keys %allBlocks;
  return (\%allBlocks, \@blkStarts);
}

# mask (i.e. eliminate some parts of) the blocks according to the read position mask
# Assumption: the block lengths should add up to the read length (i.e. no in-dels)
# input:	a block in its original string form in the sam file ("a:b,c:d...")
# input:	the reference to the position mask
# output:	masked block in its string form ("a:b,c:d...")
sub maskBlock{
  my ($block, $readLen, $maskRef) = @_;
  my @mask = @$maskRef;
  my $ttlBlkLen = 0; #for checking in-dels
  #print "maskBlock\norg: $block\n";
  my @blks = split(',', $block);
  my @newBlks = ();
  my $readOffset = 0;
  for(@blks){
    my ($s, $e) = split(':', $_);
    my $len = $e-$s+1;
    #print "analyzing $s:$e, len $len\n";
    $ttlBlkLen += $len;
    my $ts = $s; #temp start for broken blocks
    my $i = 0;
    for($i = 0; $i+$readOffset < @mask && $i<$len; $i++){
      if($mask[$i+$readOffset]){ # break point appears
        my $break = $s+$i;
        #print "break at $break (pos $i+$readOffset)\n"; 
        #check if it is a non-empty sub-block
        if($ts <= $break-1){
          push(@newBlks, $ts.":".($break-1));
        }
        #set the next s, so the skipped position is exactly $break
        $ts = $break+1;
      }
    }
    # have to handle the left-overs
    if($ts<=$e){ #there are left-overs
      push(@newBlks, $ts.":".$e);
    }
    $readOffset += $len; #this block is over
  }
  #print "new: ", join(',', @newBlks), "\n";
  #print "Read len: $ttlBlkLen\n";
  if($ttlBlkLen != $readLen){ die "ERROR: wrong total block length: $ttlBlkLen (read length $readLen)\n"; }
  return join(',', @newBlks);
}

# process SNV candidates collected from Jae-Hyung's SAM attributes
# input:	the string containing all raw information of SNVs
# input:	the reference to the position mask
# output:	the reference to the SNV hash, in which SNVs are added
#   column 16: all mismatch genomic coordinates, separated by space " "
#   column 17: all mismatch reference sequences, separated by space " " (sequence on + strand of genome)
#   column 18: all mismatch read sequences, separated by space " " (sequence on + strand of genome)
#   column 19: all mismatch read qualities, speparated by space " "
#   column 20: all mismatch positions in the read (relative to the read sequence which is already converted to + strand of genome sequence), seprated by space " "
sub procSnv{
  my ($allSnvs, $maskRef) = @_;

  my @SnvArray = split('\n', $allSnvs);
  my %hs = ();
  my @mask = @$maskRef;
  
  for(@SnvArray){
    my ($coord, $refAl, $misAl, $qual, $pos) = split('\t', $_);
    #print "($coord, $refAl, $misAl, $qual, $pos)\n";
    
    my @coords = split(' ', $coord);
    my @refAls = split(' ', $refAl);
    my @misAls = split(' ', $misAl);
    my @quals = split(' ', $qual);
    my @poss = split(' ', $pos);

    my $no = @coords;
    for(my $i = 0; $i < $no; $i++){
      my $qScr = ord($quals[$i]) - 33;
      my $loc = $coords[$i];

      if($poss[$i]-1 <@mask && $mask[$poss[$i]-1]){ #sorry, 0-based for CSers, internally
        print "Skip mismatch [ref: $refAls[$i] alt: $misAls[$i]] at $loc as position $poss[$i] of the read is masked\n";
        next; # if the position is in the mask, skip it (note mask is 0-based)
      }
      if(!defined($hs{$loc})){
        my %snv =
        (
          "ref" => $refAls[$i],
	  $misAls[$i] => 1,
	);
	  #"$misAls[$i]_q" => $qScr,
	  #"$misAls[$i]_p" => $poss[$i],
        #);
        $hs{$loc} = \%snv; 
      }else{
        my %snv = %{$hs{$loc}};
     
        if($snv{"ref"} ne $refAls[$i]){
          print "ERROR: different reference allele at $loc: ".$snv{"ref"}."VS $refAls[$i]\n"; exit;
        }
        if(!defined($snv{$misAls[$i]})){
          $snv{$misAls[$i]} = 1;
          #$snv{"$misAls[$i]_q"} = $qScr;
          #$snv{"$misAls[$i]_p"} = $poss[$i];
        }else{
          $snv{$misAls[$i]} += 1;
          #$snv{"$misAls[$i]_q"} .= "\t".$qScr;
          #$snv{"$misAls[$i]_p"} .="\t".$poss[$i];
        }
        $hs{$loc} = \%snv;
      }
    }
  }
  #my @snvs_idx = sort{$a<=>$b} keys %hs;
  
  return \%hs; #, \@snvs_idx);
}

__END__

=head1 NAME

procReads.pl -- Processing a duplicate-removed SAM file (L<rmDup>) of a chromosome (Dr. JH Lee's format) to generate the chromosome specific SNV list and the bedgraph file. The output files are used as input files for the ASARP pipeline.

=head1 SYNOPSIS

This is part of the full pre-processing:

=over 6

1. rmDup (removing PCR duplicates for SAM files in Dr. JH Lee's format)

2. mergeSam (merging SAM files if there are independent duplicates)

3. B<procReads> (processing SAM files to get SNV read counts and generate bedgraph files) 

=back

USAGE:

 perl procReads.pl input_sam_file input_snvs output_snvs output_bedgraph is_paired_end [discarded_read_pos]

NOTE:

the read processing script is for Dr. Jae-Hyung Lee's
20-attribute SAM file output format, used in RNA-editing
or allele specific expression (ASE) studies

ARGUMENTS:

 input_sam_file		SAM file input after duplicate removal (use rmDup.pl)
 intput_snvs		input SNV list (without read counts)
 output_snvs		output SNV candidates with read counts
 output_bedgraph	output bedgraph file, see below for the details:
			http://genome.ucsc.edu/goldenPath/help/bedgraph.html
 is_paired_end		0: single-end; 1: paired-end
			For paired-end reads, all reads should be paired up, 
			where pair-1 should be always followed by pair-2 in the next line.

OPTIONAL:

 discarded_read_pos	masked-out (low-quality) read positions in calculating 
			the max read quality scores, 
			in 1-based, inclusive, interval (a:b,c:d,... no space) format:
			e.g. 1:1,61:70 will discard the 1st, 61st-70th read positions.
			NOTE: the remaining reads will still contain the positions.

=head1 DESCRIPTION

C<input_sam_file> should contain only 1 chromosome, and it should be in Dr. Jae-Hyung Lee's SAM format (check out http://www.ncbi.nlm.nih.gov/pubmed/21960545 for more details)

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

=head1 SEE ALSO

L<rmDup>, L<mergeSam>, L<asarp>

=head1 COPYRIGHT

This pipeline is free software; you can redistribute it and/or modify it given that the related works and authors are cited and acknowledged.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.

=head1 AUTHOR

Cyrus Tak-Ming CHAN

Xiao Lab, Department of Integrative Biology & Physiology, UCLA

=cut
