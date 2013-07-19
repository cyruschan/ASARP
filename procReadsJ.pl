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
my $samType = "Dr. Jae-Hyung Lee's
	20-attribute SAM file output format, used in RNA-editing
	or allele specific expression (ASE) studies. If you use
	standard 10-attribute SAM files, you can use procReads.pl
	provided in this pipeline, samtools or bedtools";
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
my ($samRef) = readSamFileJ($samFile);
my ($blocksRef, $snvRef, $delRef, $blocksRcRef, $snvRcRef, $delRcRef) = parseSamReadsJ($chrToCheck, $samRef, $pairEnded, $strandFlag); 

#####################################################################
# bedgraph handling
#
my ($bedRef, $bedIdxRef) = procBedgraph($chrToCheck, $blocksRef);
my ($bedRcRef, $bedRcIdxRef) = (undef, undef);
if($strandFlag){ #need to handle the minus information as well
  ($bedRcRef, $bedRcIdxRef) = procBedgraph($chrToCheck, $blocksRcRef);
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
  my ($rnaSnvs, $dCnt) = getSnvReadsJ($chrToCheck, $snvRef, $delRef, $dnaSnvsRef, $dnaSnvsIdxRef, $bedRef, $bedIdxRef, '+');
  my ($rnaSnvsRc, $dRcCnt) = getSnvReadsJ($chrToCheck, $snvRcRef, $delRcRef, $dnaSnvsRef, $dnaSnvsIdxRef, $bedRcRef, $bedRcIdxRef, '-');
  print SP $rnaSnvs;
  print SP $rnaSnvsRc;
  $ttlSnvs += $dCnt + $dRcCnt;
}else{
  my ($rnaSnvs, $dCnt) = getSnvReadsJ($chrToCheck, $snvRef, $delRef, $dnaSnvsRef, $dnaSnvsIdxRef, $bedRef, $bedIdxRef);
  print SP $rnaSnvs;
  $ttlSnvs = $dCnt;
}
close(SP);

print "\nFinished with $ttlSnvs SNVs\n";

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
    $sam[$cnt]  = join("\t", $attr[0], $attr[1], $attr[2], $attr[3], length($attr[9]), $attr[11], $attr[13], $attr[15], $attr[16], $attr[17], $attr[19]); # id strand chr start read-len blocks mismatch_cnt loc(s) ref_allele(s) alt_allele(s) relative_pos: good enough
    ++$cnt;
    if($cnt%($INTRVL*10) == 0){ print "$cnt...";  }
  }
  close (FP);
  my $unit = 'lines';
  print " $cnt $unit. Done.\n";

  return \@sam;
}


# this is specifically for JH's SAM file (with extra attributes)
# add JH's SNVs to the read pair list: new indices used
# 0	1    2	  3	4 	5	6	    7   	8 		9	10
# id strand chr start read-len blocks mismatch_cnt loc(s) ref_allele(s) alt_allele(s) relative_pos: good enough
sub addSnvInPairJ{
  my ($ref, $attr, $discardRef) = @_;
  my %snv = %$ref;
  my @discard = ();
  my $del = '';
  if(defined($discardRef)){	 @discard = @$discardRef;	}
  
  if($$attr[6]){ # need to consider match
    my @locs = split(' ', $$attr[7]); #genome locations
    my @refs = split(' ', $$attr[8]); #reference alleles
    my @alts = split(' ', $$attr[9]); #alternative alleles
    my @poss = split(' ', $$attr[10]); #alternative alleles
    #not all used

    for(my $j=0; $j<@locs; $j++){
      my $val = "$alts[$j]\t$refs[$j]";
      if(defined($snv{$locs[$j]})){
        if($val ne $snv{$locs[$j]}){
	  print "WARNING: INCONSISTENT SNV at $locs[$j]: $val diff from $snv{$locs[$j]}\n";
	  delete $snv{$locs[$j]}; # inconsistent case
	  if($del ne ''){
	    $del .= "\t";
	  }
	  $del .= "$locs[$j]";
	}
      }else{
        if(@discard > 0 && $discard[$poss[$j]-1]){ next; }
        $snv{$locs[$j]} = $val;
      }
    }
  }
  #handling of discarded position
  if(@discard > 0){ # need to handle discarded positions
      $$attr[5] = maskBlock($$attr[5], $$attr[4], $discardRef);
  }
  return (\%snv, $$attr[5], $del); # the block
}


sub parseSamReadsJ
{
  my ($chrToCheck, $samRef, $pairEnded, $strandFlag) = @_;
  print "Processing the reads to get blocks and raw SNV statistics...\n";
  
  my @sam = @$samRef;
  my $N = @sam; # no. of reads

  my $snv = ''; # to get the snv candidates
  my $blocks = ''; # all the blocks
  my $dels = '';
  #strand-specific: -
  my $snvRc = '';
  my $blocksRc = '';
  my $delsRc = '';

  for(my $cnt = 0; $cnt < $N; $cnt ++){
    if($cnt%($INTRVL) == 0){ print "$cnt...";  }
    #if($cnt > $INTRVL/10*4){ last; }
    my @attr = split('\t', $sam[$cnt]); # pair1
    my $strand = getStrandInRead($attr[1], $strandFlag);

    my %snvInPair = (); my $sipRef = \%snvInPair;
    # jsam: JH's sam file  
    ($sipRef, my $block) = addSnvInPairJ($sipRef, \@attr);   #my $block = $attr[13];
  
    # paired-end RNA-Seq cases; strand-specific is also handled here
    my @attr2 = ();
    my $del = "";
    if($pairEnded){ #get pair2
      ++$cnt;
      if(!defined($sam[$cnt])){ # pair2
       die "ERROR: missing pair 2 in $chrToCheck at line $cnt\n";
      }
      my @attr2 = split('\t', $sam[$cnt]); # pair2
      checkPairId($attr[0], $attr2[0], $cnt);
      checkPairChr($attr[2], $attr2[2], $cnt);
      #jsam: $del only happen when pair reads overlap and SNV is inconsistent 
      ($sipRef, my $block2, $del) = addSnvInPairJ($sipRef, \@attr2);   #my $block2 = $attr2[13];
   
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
      if($dels eq ''){
        $dels .= $del;
      }else{
        $dels .= "\t$del";
      }
      if($blocks eq ''){
        $blocks = $block;
      }else{
        $blocks .= ",".$block;
      }
    }else{ # - strand
      $snvRc .= $snvToAdd; #strand-specific: -; following pair1!
      if($delsRc eq ''){
        $delsRc .= $del;
      }else{
        $delsRc .= "\t$del";
      }
      if($blocksRc eq ''){
        $blocksRc = $block;
      }
      else{
        $blocksRc .= ",".$block;
      }
    }
  }
  print "Done.\n";
  #print "$snv\n$snvRc\n";

  return (\$blocks, \$snv, \$dels, \$blocksRc, \$snvRc, \$delsRc);
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
  my @mask = ();
  if(defined($maskRef)){
    @mask = @$maskRef;
  }
  
  for(@SnvArray){
      my ($loc, $al, $refAl) = split(/\t/, $_);
      if(!defined($hs{$loc})){
        my %snv =
        (
          "ref" => $refAl, # no ref needed in the simplified version
	  $al => 1,
	);
        $hs{$loc} = \%snv; 
      }else{
        my %snv = %{$hs{$loc}};
        if(!defined($snv{"ref"})){
	  print "WARNING: no ref allele for $loc; check sam line containing: $loc\t$al\t$refAl\n";
	  next; 
	}
        if($snv{"ref"} ne $refAl){
          print "ERROR: different reference allele at $loc: ".$snv{"ref"}."VS $refAl\n"; exit;
        }
        if(!defined($snv{$al})){
          $snv{$al} = 1;
        }else{
          $snv{$al} += 1;
        }
        $hs{$loc} = \%snv;
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

  my ($chr, $snvChrRef, $delRef, $dnaSnvsRef, $dnaSnvsIdxRef, $bedRef, $bedIdxRef, $strand, $discardRef) = @_;

  my %dels = ();
  if($$delRef ne ''){
    my @array = split(/\t/, $$delRef);
    for(@array){ $dels{$_} += 1; }
  }

  if(!defined($strand)){  $strand = "";  }
  else{
    $strand = " $strand";
  }
  my %dnaSnvs = %$dnaSnvsRef;
  my @dnaSnvs_idx = @$dnaSnvsIdxRef;

  my @bedgraph = @$bedRef;
  my @bedgraph_idx = @$bedIdxRef;

  my $snvStr = ""; #result to be return

  # here we can forget about the blocks and focus on the pileups (bedgraph)
  my ($snvRef) = procSnvJ($snvChrRef, $discardRef); 
  my %snvs = %$snvRef;
  my $rSnvNo = keys %snvs;
  print " $rSnvNo SNV locations with mismatches. Done.\n";

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
	if(defined($dels{$dnaSnvs_idx[$si]}) && $dels{$dnaSnvs_idx[$si]}){
	  $refCnt -= $dels{$dnaSnvs_idx[$si]};
	  #print "NOTE: adjust reference count by $dels{$dnaSnvs_idx[$si]} overlapping reads with inconsistent SNVs recorded at $dnaSnvs_idx[$si]\n";
	}
	if($refCnt < 0){
	  print STDERR "ERROR: there are complicating overlapping read pairs at $dnaSnvs_idx[$si], making the total read count and reference count ($refCnt) inaccurate; check the original SAM file for overlapping read pairs and inconsistent SNVs\n";
	  exit;
	  $refCnt = 0;
	}
	my $altCnt = 0;
	if(defined($snvHs{$dAltAl})){  
	  $altCnt = $snvHs{$dAltAl};  
	}
        my $wrongCnt = $misCnt - $altCnt;
        if($altCnt == $maxAltCnt && $altCnt>$wrongCnt){ 
          $snvStr .= join(" ", $chr, $dnaSnvs_idx[$si], $dAls, $dId, $refCnt.":".$maxAltCnt.":".$wrongCnt).$strand."\n";
	  $dCnt++;
        }else{
	  print "DISCARD: $chr:$dnaSnvs_idx[$si] alt $dAltAl rna count <= other alts: $altCnt <= $wrongCnt\n"; 
	}
      }
    }else{
      # fully match, i.e. only the reference SNV appears here
      # $refCnt is just $cntBed here
      $snvStr .= join(" ", $chr, $dnaSnvs_idx[$si], $dnaSnvs{$dnaSnvs_idx[$si]}, $cntBed.":0:0").$strand."\n";
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

procReadsJ.pl -- Processing a duplicate-removed JSAM file (L<rmDup>) of a chromosome (Dr. JH Lee's format) to generate the chromosome specific SNV list and the bedgraph file. The output files are used as input files for the ASARP pipeline. See L<procReads> for the version on standard SAM files.

=head1 SYNOPSIS

This is part of the full pre-processing:

=over 6

1. rmDup (removing PCR duplicates for SAM files in Dr. JH Lee's format)

2. mergeSam (merging SAM files if there are independent duplicates)

3. B<procReads> (processing SAM files to get SNV read counts and generate bedgraph files) 

=back

USAGE:

 perl procReadsJ.pl chr input_sam_file input_snvs output_snvs output_bedgraph is_paired_end [discarded_read_pos]

NOTE:

the read processing script is for Dr. Jae-Hyung Lee's
20-attribute SAM file output format, used in RNA-editing
or allele specific expression (ASE) studies. If you use
standard 10-attribute SAM files, you can use procReads.pl
provided in this pipeline, samtools or bedtools

The complete documentation with the same arguments and usage can be found in L<procReads>
The only difference between procReadsJ.pl and procReads.pl is that they work on jsam and 
standard sam (jsam also treated as standard sam) files respectively.

=head1 SEE ALSO

L<rmDup>, L<mergeSam>, L<asarp>

=head1 COPYRIGHT

This pipeline is free software; you can redistribute it and/or modify it given that the related works and authors are cited and acknowledged.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.

=head1 AUTHOR

Cyrus Tak-Ming CHAN

Xiao Lab, Department of Integrative Biology & Physiology, UCLA

=cut
