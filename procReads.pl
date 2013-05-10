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
USAGE: perl $0 input_sam_file input_snvs output_snvs output_bedgraph is_paired_end [is_strand_sp bedgraph_title discarded_read_pos]

NOTE: 	the read processing script is for Dr. Jae-Hyung Lee's
	20-attribute SAM file output format, used in RNA-editing
	or allele specific expression (ASE) studies

input_sam_file 		SAM file input after duplicate removal (use rmDup.pl, and optionally mergeSam.pl)
intput_snvs		input SNV list (without read counts)
output_snvs		output SNV candidates with read counts
output_bedgraph		output bedgraph file, see below for the details:
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

discarded_read_pos	masked-out (low-quality) read positions in calculating 
			the max read quality scores, 
			in 1-based, inclusive, interval (a:b,c:d,... no space) format:
			e.g. 1:1,61:70 will discard the 1st, 61st-70th read positions.
			NOTE: the remaining reads will still contain the positions.

EOT
  exit;
}

my ($samFile, $snvFile, $outputSnvs, $outputBedgraph, $pairEnded, $strandFlag, $title, $discardPos) = @ARGV;
# handling empty title
if(!defined($title)){
  $title = "default"; #default
}
# handling undefined $strandFlag for backward compatibility
if(!defined($strandFlag)){
  print "WARNING: is_strand_sp should be input to specify whether the RNA-Seq data (input_sam_file) are strand specific.\nSet to be 0: non-strand specific for backward compatibility\n";
  $strandFlag = 0;
}elsif($strandFlag == 0){
  print "NOTE: non-strand-specific (default)\n";
}elsif($strandFlag == 1){
  print "IMPORTANT NOTE: pair 1 is sense in the strand-specific setting\n";
}elsif($strandFlag == 2){
  print "IMPORTANT NOTE: pair 1 is anti-sense in the strand-specific setting\n";
}else{
  die "ERROR: have to set a specific strand-specific flag: 0/1/2; unkonwn flag set: $strandFlag\n";
}

our $INTRVL = 100000; #interval to output processed counts

my @discard = ();
if(defined($discardPos)){
  @discard = readPosMask($discardPos);
}

print "Processing SAM file: $samFile...";
open(FP, "<", $samFile) or die "ERROR: Can't open $samFile";
my %snv = (); # to get the snv candidates
my %blocks = (); # all the blocks

#strand-specific
my %snvRc = ();
my %blocksRc = ();

my $cnt = 0;
while(<FP>){
  $cnt++;
  if($cnt%$INTRVL == 0){ print "$cnt...";  }
  my $pair1 = $_;
  chomp $pair1;
  my @attr = split('\t', $pair1);
  # we need no. of matches (attr 14, incl. known snps) to further process
  my ($id, $strand, $chr, $start, $block, $mismatches) = ($attr[0], $attr[1], $attr[2], $attr[3], $attr[11], $attr[13]); #read quality rather than mapping quality is used
  my $snvToAdd = "";
  # handle the strand info
  if($strand & 16){ # i.e. 0x10 according to SAM format: http://samtools.sourceforge.net/SAM1.pdf
    $strand = '-';
    if(defined($strandFlag) && $strandFlag == 2){ # pair 1 is anti-sense
      $strand = '+'; #flip pair 1's strand
    }
  }else{
    $strand = '+'; 
    if(defined($strandFlag) && $strandFlag == 2){ # pair 1 is anti-sense
      $strand = '-'; #flip pair 1's strand
    }
  }
  if($mismatches > 0){ #we can record SNV candidates
    # attr 16-20
    $snvToAdd .= "$attr[15]\t$attr[16]\t$attr[17]\t$attr[18]\t$attr[19]\n";
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
      if($chr ne $chr2){ 
	die "ERROR: reads have to be in pairs in the same chromosome: $id in $chr different from $id2 in $chr2 at line $cnt\n";
      }
      # add back the check for id and id2: they should be identical or differ only with the last 2 characters, typically: /1 and /2
      if($id ne $id2){
        my $idLen = length($id)-2; # ignoring the last 2 characters
        my $id1p = substr($id, 0, $idLen);
	my $id2p = substr($id2, 0, $idLen);
	if($id1p ne $id2p){
	  die "ERROR: reads have to be in pairs: $id in $chr differ from $id2 in $chr2 (more than the last 2 characters) at line $cnt\n";
	}
      
      }
      if($mismatches2 > 0){
        $snvToAdd .= "$attr2[15]\t$attr2[16]\t$attr2[17]\t$attr2[18]\t$attr2[19]\n";
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
  if(!$strandFlag || $strand eq '+'){ #non-strand specific or +; following pair1
    $snv{$chr} .= $snvToAdd;
    if(!defined($blocks{$chr})){
      $blocks{$chr} = $block;
    }
    else{
      $blocks{$chr} .= ",".$block;
    }
  }else{ # - strand
    $snvRc{$chr} .= $snvToAdd; #strand-specific: -; following pair1!
    if(!defined($blocksRc{$chr})){
      $blocksRc{$chr} = $block;
    }
    else{
      $blocksRc{$chr} .= ",".$block;
    }
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
open(SNP, "<", $snvFile) or die "ERROR: Can't open $snvFile\n";
while(<SNP>){
  chomp;
  $dSnvCnt++;
  if($dSnvCnt%$INTRVL == 0){ print "$dSnvCnt...";  }
  my ($chr, $pos, $alleles, $id) = split(' ', $_);
  if(@allChrsInSam == 0){ # no sam file
    last; # go to the last step
  }
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
if(@allChrsInSam>0){
  print " $dSnvCnt SNVs ($dSnvChrCnt kept for $allChrsInSam[0]). Done.\n";
}else{
  print "WARNING: no data in the SAM file: $samFile\n";
}

#print "SNVs in chromosomes\n";
#for(keys %snvList){
#  my @chrSnvs = split('\t', $snvList{$_});
#  print "$_: ".(scalar @chrSnvs)."\n";
#}

#############################################################
my $ttlSnvs = 0;
# sort all the block starts and make pile-ups
for my $chr (keys %blocks){
  print "Processing $chr mapped blocks and bedgraph...\n";
  my ($blkRef, $blkIdxRef) = procBlocks($blocks{$chr});
  my ($bedRef, $bedIdxRef) = procBeds($blkRef, $blkIdxRef, $chr);
  my @bedgraph = @$bedRef;
  my @bedgraph_idx = @$bedIdxRef;

  my @bedgraphRc = ();
  my @bedgraphRc_idx = ();
  if($strandFlag){ #need to handle the minus information as well
    print "\nstrand-specific set: additional - blocks and bedgraph are processed...\n";
    my ($blkRcRef, $blkRcIdxRef) = procBlocks($blocksRc{$chr});
    my ($bedRcRef, $bedIdxRcRef) = procBeds($blkRcRef, $blkRcIdxRef, $chr);
    @bedgraphRc = @$bedRcRef;
    @bedgraphRc_idx = @$bedIdxRcRef;
    
  }

  ####################################################################
  # output bedgraph
  open(BP, ">", $outputBedgraph) or die "ERROR: cannot open $outputBedgraph to output bedgraph\n";
  # output bedgraph for this chromosome
  my $addDescr = "$title; ";
  # get the track opt first
  my ($dummyChr, $dummyS, $bedLastPos) = split(" ", $bedgraph[-1]);
  my $bedFirstPos = $bedgraph_idx[0];
  if($strandFlag){
    $addDescr .= "strand-specific: the extra field is strand";
    my ($dummyChr2, $dummyS2, $bedLastPos2) = split(" ", $bedgraphRc[-1]);
    if(defined($bedLastPos2) && $bedLastPos2 > $bedLastPos ){
      $bedLastPos = $bedLastPos2;
    }
    if($bedgraphRc_idx[0] < $bedFirstPos){
      $bedFirstPos = $bedgraphRc_idx[0];
    }
  }
  my $trackRange = "$chr:".($bedFirstPos+1).":".$bedLastPos; #start converted back to 1-based for chr range in UCSC
  my $trackOpt = "track type=bedGraph name=\"reads_$chr\" description=\"$trackRange $addDescr\" visibility=full autoScale=on gridDefault=on graphType=bar yLineOnOff=on yLineMark=0 smoothingWindow=off alwaysZero=on\n";
  print BP $trackOpt;
  my $bedCnt = 0;
  print "Outputting bedgraph file...";
  if($strandFlag){
    print "\nstrand-specific set: additional - lines are output in the bedgraph...\n";
    for(@bedgraph){
      print BP $_." +\n"; #strand added
      $bedCnt++;
      if($bedCnt%$INTRVL == 0){ print "$bedCnt...";  }
    }
    for(@bedgraphRc){
      print BP $_." -\n"; #strand added
      $bedCnt++;
      if($bedCnt%$INTRVL == 0){ print "$bedCnt...";  }
    }

  }else{ # non-strand specific
    for(@bedgraph){
      print BP $_."\n";
      $bedCnt++;
      if($bedCnt%$INTRVL == 0){ print "$bedCnt...";  }
    }
  }
  close(BP);
  print " $bedCnt lines. Done.\n";

  if($strandFlag){
    print "Outputting standard 4-attribute bedgraph files ($outputBedgraph.plus and .minus), one for each strand...";
    open(PP, ">", "$outputBedgraph.plus") or die "ERROR: cannot open $outputBedgraph.plus to output bedgraph\n";
    my ($bpp1, $bpp2) = (0, 1); 
    if(@bedgraph_idx >0){
      $bpp1 = $bedgraph_idx[0];
      my ($dummyChr, $dummyS, $bedLastPos) = split(" ", $bedgraph[-1]);
      $bpp2 = $bedLastPos;
    } 
    my $plusRange = "$chr:".($bpp1+1).":".$bpp2; #start converted back to 1-based for chr range in UCSC
    print PP "track type=bedGraph name=\"reads_".$chr."_plus\" description=\"$plusRange $title; + strand only\" visibility=full autoScale=on gridDefault=on graphType=bar yLineOnOff=on yLineMark=0 smoothingWindow=off alwaysZero=on\n";
    for(@bedgraph){
      print PP $_."\n"; #strand added
    }
    close(PP);
    
    
    open(MP, ">", "$outputBedgraph.minus") or die "ERROR: cannot open $outputBedgraph.minus to output bedgraph\n";
    my ($bmp1, $bmp2) = (0, 1); 
    if(@bedgraphRc_idx >0){
      $bmp1 = $bedgraphRc_idx[0];
      my ($dummyChr2, $dummyS2, $bedLastPos2) = split(" ", $bedgraphRc[-1]);
      $bmp2 = $bedLastPos2;
    } 
    my $minusRange = "$chr:".($bmp1+1).":".$bmp2; #start converted back to 1-based for chr range in UCSC
    print MP "track type=bedGraph name=\"reads_".$chr."_minus\" description=\"$minusRange $title; - strand only\" visibility=full autoScale=on gridDefault=on graphType=bar yLineOnOff=on yLineMark=0 smoothingWindow=off alwaysZero=on\n";
    for(@bedgraphRc){
      print MP $_."\n"; #strand added
    }
    close(MP);
    print " ".(scalar @bedgraph)." + lines and ".(scalar @bedgraphRc)." - lines. Done. \n";

  }
  ####################################################################
  
  print "Processing $chr genomic SNVs...";
  my @dnaSnvVals = ();
  if(defined($snvList{$chr})){
    @dnaSnvVals = split('\t', $snvList{$chr});
  }
  my %dnaSnvs = ();
  for(@dnaSnvVals){
    my ($pos, $als, $id) = split(' ', $_);
    $dnaSnvs{$pos} = $als." ".$id;
  }
  my @dnaSnvs_idx = sort{$a <=> $b} keys %dnaSnvs;
  my $dnaSnvNo = keys %dnaSnvs;
  print " $dnaSnvNo SNVs. Done.\n";
  ####################################################################
  
  print "Processing $chr candidate SNVs (mismatches)...\n";
  if(@dnaSnvs_idx == 0){
    #speed up as no need to check the others
    print "Skipped the rest as no matched SNVs. Done\n";
    next;
  }

  # for strand-specific version, need to process them separately
  
  # output SNVs
  open(SP, ">", $outputSnvs) or die "ERROR: cannot open $outputSnvs to output candidate SNVs\n";
  my ($snvStr, $snvRcStr) = ('', '');
  my ($dCnt, $dRcCnt) = (0, 0);
  ($snvStr, $dCnt) = getSnvReads($chr, \$snv{$chr}, \%dnaSnvs, \@dnaSnvs_idx, \@bedgraph, \@bedgraph_idx, \@discard);
  if($strandFlag){
    print "\nstrand-specific set: additional - SNVs are output in the snv file...\n";
    ($snvRcStr, $dRcCnt) = getSnvReads($chr, \$snvRc{$chr}, \%dnaSnvs, \@dnaSnvs_idx, \@bedgraphRc, \@bedgraphRc_idx, \@discard);
    # need to merge the SNVs from forward and backward strands
    my @snvsAll = split(/\n/, $snvStr);
    my @snvsAllRc = split(/\n/, $snvRcStr);

    for(@snvsAll){
      print SP "$_ +\n";
    }
    for(@snvsAllRc){
      print SP "$_ -\n";
    }
  }else{
    print SP $snvStr; #no strand-specific info at the end
  
  }
  close(SP);
  $ttlSnvs += $dCnt + $dRcCnt;
}
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

# process bedgraph (non-strand specific)
# input:	the references to the blocks and their indices
# output:	the references to the begraph and its index 
sub procBeds{

  my ($blkRef, $blkIdxRef, $chr) = @_;

  my %thisChrBlks = %$blkRef;
  my @blkStarts = @$blkIdxRef;
  my $blkNo = @blkStarts;
  print " $blkNo block starts. Done.\n";

  #print "Processing $chr bedgraph...";
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

  return (\@bedgraph, \@bedgraph_idx);

}

# process SNV candidates collected from Jae-Hyung's SAM attributes
# input:	the reference to the string containing all raw information of SNVs
# input:	the reference to the position mask
# output:	the reference to the SNV hash, in which SNVs are added
#   column 16: all mismatch genomic coordinates, separated by space " "
#   column 17: all mismatch reference sequences, separated by space " " (sequence on + strand of genome)
#   column 18: all mismatch read sequences, separated by space " " (sequence on + strand of genome)
#   column 19: all mismatch read qualities, speparated by space " "
#   column 20: all mismatch positions in the read (relative to the read sequence which is already converted to + strand of genome sequence), seprated by space " "
sub procSnv{
  my ($allSnvsRef, $maskRef) = @_;
  my $allSnvs = $$allSnvsRef;
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
        if(!defined($snv{"ref"})){
	  print "WARNING: no ref allele for $loc; check sam line containing: $coord\t$refAl\t$misAl\t$qual\t$pos\n";
	  next; 
	}
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

# process the SNVs and get their read counts
# input: references to RNA SNVs, DNA SNVs, bedgraph [discard list] and their indices
# output: the string containing all SNVs with reads and the count

sub getSnvReads{

  my ($chr, $snvChrRef, $dnaSnvsRef, $dnaSnvsIdxRef, $bedRef, $bedIdxRef, $discardRef) = @_;

  my %dnaSnvs = %$dnaSnvsRef;
  my @dnaSnvs_idx = @$dnaSnvsIdxRef;

  my @bedgraph = @$bedRef;
  my @bedgraph_idx = @$bedIdxRef;

  my $snvStr = ""; #result to be return

  # here we can forget about the blocks and focus on the pileups (bedgraph)
  my ($snvRef) = procSnv($snvChrRef, $discardRef); 
  my %snvs = %$snvRef;
  my $rSnvNo = keys %snvs;
  print " $rSnvNo SNV locations. Done.\n";

  print "Processing $chr SNVs with read counts...\n";
  my ($bi, $si) = (0, 0);#b for bedgraph_idx, s for dnaSnv_idx, 
  my $dCnt = 0; #$dCnt for matched dnaSnvs_idx
  while($bi<@bedgraph_idx && $si<@dnaSnvs_idx){
    #print "bed idx $bi: ".($bedgraph_idx[$bi]+1)."\nsnv idx: $si: $dnaSnvs_idx[$si]\n";
    if($bedgraph_idx[$bi]+1 > $dnaSnvs_idx[$si]){ # the current SNV is smaller, need to check next one
      $si++;
      next;
    }
    # now the SNV potentially overlaps with the bedgraph range
    # i.e. $bedgraph_idx[$bi]+1 <= $dnaSnvs_idx[$si]
    my ($chrBed, $startBed, $endBed, $cntBed) = split(' ', $bedgraph[$bi]);
    if($endBed < $dnaSnvs_idx[$si]){ # the SNV is larger than the whole area
      $bi++;
      next;
    }
    # now $bedgraph_idx[$bi]+1 <= $dnaSnvs_idx[$si] and $endBed >= $dnaSnvs_idx[$si] so it is a match
    # a match
    #print "SNV $si ($dnaSnvs_idx[$si]) matches bedgraph $bi: $bedgraph[$bi]+1 = $startBed to $endBed\n";
    #if(!($dnaSnvs_idx[$si] > $startBed && $dnaSnvs_idx[$si] <= $endBed)){
    #  print "Wrong: !($dnaSnvs_idx[$si] > $startBed && $dnaSnvs_idx[$si] < $endBed)\n";
    #  exit;
    #}
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
          $snvStr .= join(" ", $chr, $dnaSnvs_idx[$si], $dAls, $dId, $refCnt.":".$maxAltCnt.":".$wrongCnt)."\n";
	  $dCnt++;
        }else{
	  print "DISCARD: $chr:$dnaSnvs_idx[$si] alt $dAltAl rna count <= other alts: $altCnt <= $wrongCnt\n"; 
	}
      }
    }else{
      # fully match, i.e. only the reference SNV appears here
      # $refCnt is just $cntBed here
      $snvStr .= join(" ", $chr, $dnaSnvs_idx[$si], $dnaSnvs{$dnaSnvs_idx[$si]}, $cntBed.":0:0")."\n";
      $dCnt++;
    }
    $si++;
    if($dCnt>0 && $dCnt%$INTRVL==0){
      print "$dCnt...\n";
    }
  }
  print " $dCnt matched SNVs. Done.\n";

  return ($snvStr, $dCnt);

}

__END__

=head1 NAME

procReads.pl -- Processing a duplicate-removed SAM file (L<rmDup>) of a chromosome (Dr. JH Lee's format) to generate the chromosome specific SNV list and the bedgraph file. The output files are used as input files for the ASARP pipeline.

The new procReads.pl supports strand-specific RNA-Seq data (i.e. the SAM file strand information is reliable) for more accurate results.

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

 discarded_read_pos	masked-out (low-quality) read positions in calculating 
			the max read quality scores, 
			in 1-based, inclusive, interval (a:b,c:d,... no space) format:
			e.g. 1:1,61:70 will discard the 1st, 61st-70th read positions.
			NOTE: the remaining reads will still contain the positions.

=head1 DESCRIPTION

=head2 INPUT

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

=head2 OUTPUT

By default, the strand specific flag C<is_strand_sp> is unset, and the C<inptu_sam_file> strand information is considered unreliable and not used. The output files are standard bedgraph files http://genome.ucsc.edu/goldenPath/help/bedgraph.html with space as the dilimiter, and SNV files with the following format:

Each line is space separated, with

 chromosome
 location 
 ref_allele>alt_allele 
 dbSnp_id 
 read_counts ref:alt:wrnt

Note that C<read_counts> are RNA read counts obtained from the SAM (a.k.a the bedgraph) file. C<ref> indicates the read count of the reference allele, C<alt> the alternative allele, C<wrnt> (wrong nt) indicates other alleles neither ref nor alt. It is required that C<alt> > C<wrnt>, otherwise that SNV is discarded (dicarded on a particular strand if strand-specific option is on). Output SNV examples would look like:

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
