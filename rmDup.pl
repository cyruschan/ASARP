#!/usr/bin/perl
use warnings;
use strict;

require "readUtilities.pl";
# set autoflush for error and output
select(STDERR);
$| = 1;
select(STDOUT);
$| = 1;

if(@ARGV < 3){
  print <<EOT;
USAGE: perl $0 input_sam_file output_sam_file is_paired_end [discarded_read_pos]
NOTE: 	the duplicate removal script is for standard SAM and Dr. Jae-Hyung Lee's
	20-attribute SAM file output formats, used in RNA-editing
	or allele specific expression (ASE) studies

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

my ($samFile, $output, $pairEnded, $discardPos) = @ARGV;

my @discard = ();
if(defined($discardPos)){
  @discard = readPosMask($discardPos);
}

open(FP, "<", $samFile) or die "ERROR: Can't open $samFile";
my %readDups = (); # to get them pair-up
my %reads = (); # to get the lines for the pair
print "Processing SAM file...";
my $cnt = 0;
my $dup = 0;
while(<FP>){
  $cnt++;
  my $pair1 = $_;
  chomp $pair1;
  my @attr = split('\t', $pair1);
  my ($id, $strand, $chr, $start, $quality) = ($attr[0], $attr[1], $attr[2], $attr[3], $attr[10]); #read quality rather than mapping quality is used
  my $key;
  if(!defined($attr[11])){ #standard SAM
   $key = parseCigar($start, $attr[5]);
  }else{
   $key = $attr[11]; #attribute 12 is the region blocks separated by ','
  }
  my $value = $pair1;

  if($pairEnded){ #get pair2
    my $pair2 = <FP>;
    $cnt++;
    if(defined($pair2)){
      chomp $pair2;
      my @attr2 = split('\t', $pair2);
      my ($id2, $strand2, $chr2, $start2, $quality2) = ($attr2[0], $attr2[1], $attr2[2], $attr2[3], $attr2[10]);
      checkPairId($id, $id2, $cnt);
      checkPairChr($chr, $chr2, $cnt);
      
      # they should be combined into a pair
      if(!defined($attr2[11])){ #standard SAM
        $key .= ",".parseCigar($start2, $attr2[5]);
      }else{
        $key .= ",".$attr2[11]; #attribute 12 is the region blocks separated by ','
      }
      $value .= "\n".$pair2;
      $quality .= $quality2;
    }else{
      print "ERROR: missing pair 2 in $chr at line $cnt\n"; exit;
    }
  }
  if($cnt%100000==0){ print "$cnt..."; }
  #calculate score for the read (pair)
  my $score = getReadScore($quality, $pairEnded, \@discard);
  #print "$id, $strand, $chr, $key: scr: $score\n";

  if(!defined($reads{$chr})){
    $reads{$chr} = {};
    $reads{$chr}->{$key} = $score."\n".$value;
  }elsif(!defined($reads{$chr}->{$key})){
    $reads{$chr}->{$key} = $score."\n".$value;
  }else{
    $dup += 1;
    #compare existing to get the one with better score
    my $oScore = 0;
    if($reads{$chr}->{$key} =~ /^(.+)\n/){
      $oScore = $1;
      #print "old score: $oScore\n";
    }else{
      print "ERROR: parsing ".$reads{$chr}->{$key}."\n";
      exit;
    }
    if($score > $oScore){ #replace
      #print "replace: using $id with $key\n";
      #print $reads{$chr}->{$key}. "\nreplaced by\n".$value."\n";
      $reads{$chr}->{$key} = $score."\n".$value;
    }#elsif($score == $oScore){
    #  print "tie: $id with $key\n";
    #}else{
    #  print "remove: $id with $key\n";
    #}

  }

}
close(FP);
print "$cnt lines\n";
print "Duplicates in total: $dup\n";
my $total = $cnt;
my $unit = 'reads';
if($pairEnded){
  $unit = 'pairs';
  $total = $cnt/2;
}
print "Total $unit $total\n";
if($total){
printf "Duplicate Percentage: %.2f %%\n", ($dup/$total*100);
}else{ print "Duplicate Percentage: N/A\n"; }

# write to the output sam file
open (OP, ">", $output) or die "ERROR: Cannot open $output for writing\n";
for(keys %reads){
  my %hs = %{$reads{$_}};
  for(keys %hs){
    my @values = split("\n", $hs{$_});
    for(my $j = 1; $j<@values; $j++){
      print OP $values[$j]."\n";
    }
  }
}
close(OP);

print "DONE\n";

sub getReadScore{
  my ($qual, $pairEnded, $maskRef) = @_;
  my @mask = @{$maskRef};
  my @quals = split("", $qual);
  my $readLen = @quals;
  if($pairEnded){
    $readLen /= 2; #assume both read pairs have the same read length
  }
  my $scr = 0;
  #print "raw: @quals\nscr: ";
  my $offset = 0;
  my $len = 0;
  for(my $i=0; $i<$readLen; $i++){
    #skip all the positions that are masked
    if($i<@mask && $mask[$i]){
      next;
    }
    #Sanger/Illumina 1.8+ encoding: Phred+33
    #http://en.wikipedia.org/wiki/FASTQ_format
    $scr += ord($quals[$offset+$i])-33;
    $len += 1;
  }

  if($pairEnded){
    # for paired-end: pair 2
    $offset = $readLen;
    for(my $i=0; $i<$readLen; $i++){
      #skip all the positions that are masked
      if($i<@mask && $mask[$i]){
        next;
      }
      #Sanger/Illumina 1.8+ encoding: Phred+33
      #http://en.wikipedia.org/wiki/FASTQ_format
      $scr += ord($quals[$offset+$i])-33;
      $len += 1;
    }
  }

  if($len>0){
    $scr = $scr/$len;
  }else{
    $scr = 0;
  }
  return $scr;
}

__END__

=head1 NAME

rmDup.pl -- Removing duplicates in a SAM file of a chromosome (Dr. JH Lee's format), where the extra 12th attribute (mapped read blocks) are used to identify distinct reads (read pairs). Reads (read pairs) are considered as duplicates only if all of their mapped read blocks have the same coordinates. Only the read (pair) with the highest read quality will be kept. The output SAM file can be used as the input file for merging of multiple independent replicates (L<mergeSam>), or read processing to generate SNV and bedgraph files (L<procReads>).

=head1 SYNOPSIS

This is part of the full pre-processing:

=over 6

1. B<rmDup> (removing PCR duplicates for SAM files (including Dr. JH Lee's SAM format); samtools/bedtools can be used for standard SAM files)

2. mergeSam (merging SAM files if there are independent duplicates)

3. procReads (processing SAM files to get SNV read counts and generate bedgraph files) 

=back


USAGE: 

 perl rmDup.pl input_sam_file output_sam_file is_paired_end [discarded_read_pos]

NOTE:

the duplicate removal script is for standard SAM and Dr. Jae-Hyung Lee's
20-attribute SAM file output formats, used in RNA-editing
or allele specific expression (ASE) studies

ARGUMENTS:

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

C<input_sam_file> should contain only 1 chromosome, and it should be in standard SAM format or Dr. Jae-Hyung Lee's SAM format (check out http://www.ncbi.nlm.nih.gov/pubmed/21960545 for more details)

=head1 SEE ALSO

L<mergeSam>, L<procReads>, L<asarp>

=head1 COPYRIGHT

This pipeline is free software; you can redistribute it and/or modify it given that the related works and authors are cited and acknowledged.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.

=head1 AUTHOR

Cyrus Tak-Ming CHAN

Xiao Lab, Department of Integrative Biology & Physiology, UCLA

=cut
