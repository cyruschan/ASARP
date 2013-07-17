#!/usr/bin/perl
use warnings;
use strict;

require "fileParser.pl"; # only to use binarySearch
require "readUtilities.pl";
# set autoflush for error and output
select(STDERR);
$| = 1;
select(STDOUT);
$| = 1;

our $INTRVL = 100000; #interval to output processed counts
my $samType = "a standard SAM file with
	10 attributes. You are recommended to use more efficient
	tools such as samtools and bedtools for this task; or
	to use procReadsJ.pl on the special jsam file introduced
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
  print "Skipped procesing SNVs with RNA reads as there are no DNA SNVs. Done\n";
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
  while(<FP>){
    chomp $_;
    my @attr = split('\t', $_);
    $sam[$cnt]  = join("\t", $attr[0], $attr[1], $attr[2], $attr[3], $attr[5], $attr[9]); # id strand chr start cigar and read: good enough
    ++$cnt;
    if($cnt%($INTRVL*10) == 0){ print "$cnt...";  }
  }
  close (FP);
  my $unit = 'lines';
  print " $cnt $unit. Done.\n";

  return \@sam;
}

sub parseSamReads
{
  my ($chrToCheck, $samRef, $pairEnded, $strandFlag) = @_;
  print "Processing the reads to get blocks and raw SNV statistics...\n";
  
  my @sam = @$samRef;
  my $N = @sam; # no. of reads

  my $snv = ''; # to get the snv candidates
  my $blocks = ''; # all the blocks
  #strand-specific: -
  my $snvRc = '';
  my $blocksRc = '';

  for(my $cnt = 0; $cnt < $N; $cnt ++){
    if($cnt%($INTRVL/10) == 0){ print "$cnt...";  }
    #if($cnt > $INTRVL/10*4){ last; }
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
      if(($strandFlag == 2 && $strand eq '+') || ($strandFlag == 1 && $strand eq '-') || ($strandFlag ==0 && $strand eq '-')){
         #non-strand specific: - here means antisense
         $block = mergeBlockInPair($block2, $block, $attr[3]); #p2 p1
      }else{ # incl. $strandFlag == 0
         $block = mergeBlockInPair($block, $block2, $attr2[3]); # p1 p2
      }
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
      }else{
        die "ERROR: unrecognized allele: $al at $loc\n";
      }
    }
  }

  my %rawSnvs = (
  'A' => \%rawSnvsA,
  'C' => \%rawSnvsC,
  'G' => \%rawSnvsG,
  'T' => \%rawSnvsT,
  );
  my @alpha = qw(A C G T);
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
	print "WARNING: DISCARD SNV $dnaSnvs{$loc} as alt count <= other nts: $altCnt <= $wrgCnt\n"; 
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

# parse the CIGAR string to get block information
sub parseCigar{
  my ($start, $cigar) = @_;
  my $position = $start;
  my $block = '';
  while ($cigar !~ /^$/){
    # handle only the matched parts
    if ($cigar =~ /^([0-9]+[MIDSN])/){
      my $cigar_part = $1;
      if ($cigar_part =~ /(\d+)M/){
        if($block ne ''){ #already some blocks
	  $block .= ",";
	}
	$block .= "$position:";
        $position += $1;
	my $end = $position - 1;
	$block .= "$end";
      } elsif ($cigar_part =~ /(\d+)N/){
        $position += $1; # skipped positions
      } else { die "ERROR: CIGAR not supported: $cigar\n"; }
      #} elsif ($cigar_part =~ /(\d+)I/){
      #} elsif ($cigar_part =~ /(\d+)D/){
      #} elsif ($cigar_part =~ /(\d+)S/){
      $cigar =~ s/$cigar_part//;
    }else{
      die "ERROR: CIGAR not supported: $cigar\n";
    }
  }
  return $block;
}

# this is for the general SAM files (not using extra attributes)
# add SNV (attributes of the original SAM file) to the read pair list
sub addSnvInPair{
  my ($ref, $ref2, $dnaSnvIdxRef, $discardRef) = @_;
  my %snv = %$ref;
  my @attr = @$ref2;
  my @dnaSnvsIdx = @$dnaSnvIdxRef;
  #my %dnaSnvs = %$dnaSnvRef;
  my $n = @dnaSnvsIdx;
  my @discard = ();
  if(defined($discardRef)){	 @discard = @$discardRef;	}
  my $start = $attr[3];
  my $cigar = $attr[4];
  my $block = parseCigar($start, $cigar);
  #print "$cigar, $block | $attr[11]\n";
  my $snvToAdd = "";

  my @blocks = split(',', $block);
  my $rePos = 0; # relative position
  for(@blocks){
    my ($s, $e) = split(':', $_);
#=pod    
    my ($lBound, $lUnMatch) = binarySearch($dnaSnvsIdxRef, $s, 0, $n-1, 'left');
    for(my $j = $lBound; $dnaSnvsIdx[$j]<= $e; $j++){ #snv is more sparse!
      my $pos = $dnaSnvsIdx[$j];
      my $readPos = $rePos + $pos-$s; # 0-based
      if(!defined($discard[$readPos]) || !$discard[$readPos]){ # not discarded
	  my $al = substr($attr[5], $readPos, 1);
	  if(defined($snv{$pos})){ # existing
	    if($snv{$pos} ne $al){
	      print "WARNING: INCONSISTENT SNV at $pos: $snv{$pos} $al\n";
	    }
	  }else{
	    $snv{$pos} = $al; #\t$qual\n"; # add new
	  }
      }

    }
    $rePos += $e-$s+1; # next block
#=cut

=pod
    for my $pos ($s .. $e){
      if(defined($dnaSnvs{$pos})){ # a snv position
        #++$c;
        # check if there is any problem with the SNV first
	if(!defined($discard[$rePos]) || !$discard[$rePos]){ # not discarded
	  my $al = substr($attr[5], $rePos, 1);
	  if(defined($snv{$pos})){ # existing
	    if($snv{$pos} ne $al){
	      print "WARNING: INCONSISTENT SNV at $_: $snv{$_} $al\n";
	    }
	  }else{
	    $snv{$pos} = $al; #\t$qual\n"; # add new
	    print "SNV $pos: $al\n";
	  }
	}
      }
      ++$rePos;
    }
=cut

  }
  #handling of discarded position
  if(@discard > 0){ # need to handle discarded positions
      my $readLen = length($attr[5]);
      $block = maskBlock($block, $readLen, $discardRef);
  }
  return (\%snv, $block);
}

########################################################################################################################################
##########################################################################
#		sub-routines (JSAM version)
##########################################################################


## read jsam file and cut most of the useless attributes out
sub readSamFileJ{

  my ($samFile) = @_;
  my @sam = ();

  my $cnt = 0;
  print "Processing SAM file: $samFile and counting SNV reads...";
  open(FP, "<", $samFile) or die "ERROR: Can't open $samFile";
  while(<FP>){
    #if($cnt > $INTRVL/100*4){ exit; }
    chomp $_;
    my @attr = split('\t', $_);
    $sam[$cnt]  = join("\t", $attr[0], $attr[1], $attr[2], $attr[3], $attr[11], $attr[13], $attr[15], $attr[16], $attr[17]); # id strand chr start blocks mismatch_cnt loc(s) ref_allele(s) alt_allele(s): good enough
    ++$cnt;
    if($cnt%($INTRVL*10) == 0){ print "$cnt...";  }
  }
  close (FP);
  my $unit = 'lines';
  print " $cnt $unit. Done.\n";

  return \@sam;
}


# this is specifically for JH's SAM file (with extra attributes)
# add JH's SNV (attribute 13, 15-19) to the read pair list
sub addSnvInPairJ{
  my ($ref, $ref2) = @_;
  my %snv = %$ref;
  my @attr = @$ref2;
  
  # we need no. of matches (attr 14, incl. known snps) to further process
  #my ($id, $strand, $chr, $start, $block, $mismatches) = ($attr[0], $attr[1], $attr[2], $attr[3], $attr[11], $attr[13]); #read quality rather than mapping quality is used
  if($attr[5]){ # need to consider match
    my @locs = split(' ', $attr[6]); #genome locations
    #my @refs = split(' ', $attr[7]); #reference alleles
    my @alts = split(' ', $attr[8]); #alternative alleles
    # all not used

    for(my $j=0; $j<@locs; $j++){
      if(defined($snv{$locs[$j]})){
        if($alts[$j] ne $snv{$locs[$j]}){
	  print "WARNING: INCONSISTENT SNV at $locs[$j]: $alts[$j] $snv{$locs[$j]}\n";
	  delete $snv{$locs[$j]}; # inconsistent case
	}
      }else{
        $snv{$locs[$j]} = $alts[$j];
      }
    }
  }
  return (\%snv, $attr[4]); # the block
}

  
#######################################################################################################################################
# process SNV candidates collected from Jae-Hyung's SAM attributes
# input:	the reference to the string containing all raw information of SNVs
# input:	the reference to the position mask
# output:	the reference to the SNV hash, in which SNVs are added
#   column 16: all mismatch genomic coordinates, separated by space " "
#   column 17: all mismatch reference sequences, separated by space " " (sequence on + strand of genome)
#   column 18: all mismatch read sequences, separated by space " " (sequence on + strand of genome)
#   column 19: all mismatch read qualities, speparated by space " "
#   column 20: all mismatch positions in the read (relative to the read sequence which is already converted to + strand of genome sequence), seprated by space " "
sub procSnvJ{
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
    #my @quals = split(' ', $qual);
    #my @poss = split(' ', $pos);

    my $no = @coords;
    for(my $i = 0; $i < $no; $i++){
      #my $qScr = ord($quals[$i]) - 33;
      my $loc = $coords[$i];

      #if($poss[$i]-1 <@mask && $mask[$poss[$i]-1]){ #sorry, 0-based for CSers, internally
      #  print "Skip mismatch [ref: $refAls[$i] alt: $misAls[$i]] at $loc as position $poss[$i] of the read is masked\n";
      #  next; # if the position is in the mask, skip it (note mask is 0-based)
      #}
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
# jsam version: JH's sam format. Extra fields are included 
# (no reference allele information directly available, need to use bedgraph to infer)

sub getSnvReadsJ{

  my ($chr, $snvChrRef, $dnaSnvsRef, $dnaSnvsIdxRef, $bedRef, $bedIdxRef, $discardRef) = @_;

  my %dnaSnvs = %$dnaSnvsRef;
  my @dnaSnvs_idx = @$dnaSnvsIdxRef;

  my @bedgraph = @$bedRef;
  my @bedgraph_idx = @$bedIdxRef;

  my $snvStr = ""; #result to be return

  # here we can forget about the blocks and focus on the pileups (bedgraph)
  my ($snvRef) = procSnvJ($snvChrRef, $discardRef); 
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
