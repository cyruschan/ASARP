#!/usr/bin/perl
use warnings;
use strict;
BEGIN { push @INC, "../" }
require 'fileParser.pl';
require 'snpParser.pl';
require 'readUtilities.pl';


# a work-around for U87 as the file format is not jsam
# this sub-routine **hard-codes** the U87 sam file location
# each time we retrieve the reads containing the current snv location, with the desired alle
sub getU87RawReads
{
  my ($chr, $pos, $al, $lrds, $rrds) = @_;
  my %cans = (); # all string candidates
  #my $u87 = "/home/cyruschan/asarp_perl/U87.DroshaKD/Dro.merged"; # DroshaKD
  # CAUTION: different files with different read lengths: 60M for WT u87; 90M for DroshaKD and siRNA
  #only work for U87 (Gang version sam data)
  my $u87 = "/home/gxxiao/apps/data/u87/cofactor2/chr.one"; # WT u87
  my $cmd = "grep $pos:$al $u87/$chr.one.sam | grep '60M' | cut -f 4,6,8 | sort -u"; # for U87 Gang version only
  print "$cmd\n";
  my $res = qx($cmd);
  print "$res\n";
  my @reads = split('\n', $res);
  for(@reads){
    my ($range, $cigar, $read) = split('\t', $_);
    if(index($range, ',')!=-1){
      #print "Skip junction: $range\t$read\n";
      next;
    }
    my ($s, $e) = split(':', $range);
    if(($pos-$lrds >= $s) && $pos+$rrds <= $e){ # in range
      my $len = $lrds+1+$rrds;
      my $str = substr($read, $pos-$lrds-$s, $len);
      $cans{$str} += 1;
    }
  }

  my ($maxK, $maxV) = (undef, 0);
  for(keys %cans){
    if($cans{$_} > $maxV){
      $maxK = $_;
      $maxV = $cans{$_};
    }
    #print "$_\t$cans{$_}\n";
  }
  if(!defined($maxK)){
     print "need genomefetch:\n";
     my ($left, $right) = ($pos-$lrds, $pos+$rrds);
     $maxK = qx(~/GenomeFetch.py -o hg19 -c $chr -f $left -t $right -s 1); # always on + strand w.r.t ref genome
     chomp $maxK;
     #print "need to replace alt\n";
	#  for(my $j = 0; $j<$lrds; $j++){
	#    print " ";
	#  } print "$al\n";
     #print "$maxK\n";
     $maxK = substr($maxK, 0, $lrds).$al.substr($maxK, $lrds+1, $rrds);
     #print "$maxK\n";
  }
  return ($maxK, $maxV);
}

sub getU87RawReads2 # using the formal way
{
  my ($chr, $pos, $al, $lrds, $rrds, $readsRef) = @_;
  my %cans = (); # all string candidates
  
  my $res = $readsRef->{$chr}->{$pos}->{$al};
  my @reads = split('\t', $res);
  my $n = @reads;
  print "@reads\n$n reads\n"; 
  for(@reads){
    my ($rPos, $read) = split(';', $_);
    if(($rPos-$lrds > 0) && $rPos+$rrds < length($read) ){ # in range
      my $len = $lrds+1+$rrds;
      my $str = substr($read, $rPos-$lrds, $len);
      print "Extracted $str from $read in which rela pos is $rPos for allele: $pos $al\n"; #exit;
      $cans{$str} += 1;
    }
  }

  my ($maxK, $maxV) = (undef, 0);
  for(keys %cans){
    if($cans{$_} > $maxV){
      $maxK = $_;
      $maxV = $cans{$_};
    }
    #print "$_\t$cans{$_}\n";
  }
  if(!defined($maxK)){
     print "need genomefetch:\n";
     my ($left, $right) = ($pos-$lrds, $pos+$rrds);
     $maxK = qx(~/GenomeFetch.py -o hg19 -c $chr -f $left -t $right -s 1); # always on + strand w.r.t ref genome
     chomp $maxK;
     #print "need to replace alt\n";
	#  for(my $j = 0; $j<$lrds; $j++){
	#    print " ";
	#  } print "$al\n";
     #print "$maxK\n";
     $maxK = substr($maxK, 0, $lrds).$al.substr($maxK, $lrds+1, $rrds);
     #print "$maxK\n";
  }
  return ($maxK, $maxV);
}


sub getRawReadLines{
  my ($gListRef, $folder, $suffix, $rds) = @_; # the folder path and suffix for the sam files
  my %genes = %$gListRef;
  my %chrSnvs = ();
  # get the raw read list for all SNVs needed
  for my $key (keys %genes){
    my ($chr, $gene) = split(';', $key);
    my %hs = ();
    if(defined($chrSnvs{$chr})){
      %hs = %{$chrSnvs{$chr}};  
    }
    my @snvs = split(/\t/, $genes{$key}); #dummyInfo;snpInfo;$strandInfo
    for my $snvLine (@snvs){
      my ($dummyInfo, $snpInfo, $strandInfo) = split(';', $snvLine);
      my ($pos, $id, $alleles, $count) = split(' ', $snpInfo);
      my ($ref, $alt) = split('>', $alleles);
      $hs{$pos}->{$ref} = '';
      $hs{$pos}->{$alt} = '';
    }
    $chrSnvs{$chr} = \%hs;
  }
  
  #print "all chrs ", keys %chrSnvs, "\n"; exit;
  print "Reading sam files ($folder/*.$suffix) and get raw reads\n";
  # get the SAM reads
  for my $chrName (keys %chrSnvs){
    my @poss = sort{$a <=>$b} keys %{$chrSnvs{$chrName}}; # all poss
    my ($minP, $maxP) = ($poss[0], $poss[-1]); # filter out useless positions
    print "$chrName range: ($minP, $maxP)\n@poss\n"; #exit;
    my $samFile = "$folder/$chrName.$suffix";
    open(FP, "<", $samFile) or die "ERROR: Can't open $samFile";
    while(<FP>){
      chomp $_;
      my @attr = split('\t', $_); # standard sam
      if($attr[3] <= $maxP && $attr[3] >= $minP && $attr[5] =~ /^\d+M$/){ # in range; not splicing
        # further check: get the relative position and 
        my $block = parseCigar($attr[3], $attr[5]);
	my ($s, $e) = split(':', $block);

	for my $pos (@poss){
	  my ($ref, $alt) = keys %{$chrSnvs{$chrName}->{$pos}};
	  if($pos-$rds >= $s && $pos+$rds <= $e){ #overlaps and in range
	    my $rPos = $pos - $s; # 0 based
	    #check the allele
	    my $nt = substr($attr[9], $rPos, 1); # the allele at that location
	    if($nt eq $ref){
	      $chrSnvs{$chrName}->{$pos}->{$ref} .= "$rPos;$attr[9]\t"; # add pos and read
	      #print "$pos $ref: $rPos;$attr[9]\n";
	    }
	    if($nt eq $alt){
	      $chrSnvs{$chrName}->{$pos}->{$alt} .= "$rPos;$attr[9]\t"; # add pos and read
	      #print "$pos $alt: $rPos;$attr[9]\n";
	    }
	  }
	}
      }
    }
    close(FP);
  }
  return \%chrSnvs;
}
# generate allele-specific fasta file
# input: 
# ASAT list
# radius (for both left and right)
#
sub genAlleleLocation
{
  my ($gListRef, $rds, $altRef) = @_;
  my %genes = %$gListRef;
  my %fasta = ();

  my $u87Fd = "/home/cyruschan/asarp_perl/U87.siRNA/U87.siRNA.merged"; # siRNA
  my $samSuf = "sam";
  my ($readsRef) = getRawReadLines($gListRef, $u87Fd, $samSuf, $rds);

  for my $key (keys %genes){
    print "$key\n";
    my ($chr, $gene) = split(';', $key);
    my $chrID = getChrID($chr);
    my ($utrRef, $strd) = getUtrByChrGene($altRef, $chrID, $gene);
    my ($luRef) = getLongestUtrs($utrRef, $strd);
    my %lUtrs = %$luRef;
    my @utrKeys = keys %lUtrs;
    @utrKeys = sort { $a <=> $b } @utrKeys; 
    my @snvs = split(/\t/, $genes{$key}); #dummyInfo;snpInfo;$strandInfo
    for my $snvLine (@snvs){
      #print "$snvLine\n";
      my ($dummyInfo, $snpInfo, $strandInfo) = split(';', $snvLine);
      # get 3' UTR region
      my ($pos, $id, $alleles, $count) = split(' ', $snpInfo);
     
      my ($lRds, $rRds, $maxL) = (0, 0, 0);
      # get the max range from the UTRs
      for(@utrKeys){
        if(($strd eq '+' && $_ <= $pos && $lUtrs{$_} >= $pos) || ($strd eq '-' && $_ >= $pos && $lUtrs{$_} <= $pos)){
	  
	  my ($left, $right) = ($_, $lUtrs{$_});
	  if($left > $right){ 
	    ($left, $right) = ($lUtrs{$_}, $_);
	  }
	  my ($lDis, $rDis) = ($pos - $left, $right - $pos);
	  my $l = $lDis+1+$rDis;
	  if($l > $maxL){
	    ($lRds, $rRds, $maxL) = ($lDis, $rDis, $l);
	  }
	
	}     
      }
      print "maxmal range for $pos: left $lRds right $rRds with len: $maxL\n";
      
      my ($ref, $alt) = split('>', $alleles);
      my ($refCnt, $altCnt) = split(':', $count);
      # restrict the range:
      if($lRds > $rds){ $lRds = $rds; }
      if($rRds > $rds){ $rRds = $rds; }
      print "final UTR range for $pos: left $lRds right $rRds\n";
      $lRds = $rds; $rRds = $rds; # use just the range set
      print "updated range (not UTR constraint) for $pos: left $lRds right $rRds\n";

      #my ($dna, $refReads) = getU87RawReads($chr, $pos, $ref, $lRds, $rRds);
      #print "Ref: $dna\t$refReads\n";
      #my ($altDna, $altReads) = getU87RawReads($chr, $pos, $alt, $lRds, $rRds);
      #print "Alt: $altDna\t$altReads\n";
      
      my ($dna, $refReads) = getU87RawReads2($chr, $pos, $ref, $lRds, $rRds, $readsRef);
      print "Ref: $dna\t$refReads\n";
      my ($altDna, $altReads) = getU87RawReads2($chr, $pos, $alt, $lRds, $rRds, $readsRef);
      print "Alt: $altDna\t$altReads\n";
      
      my $inPos = $lRds+1;
      if($strd eq '-'){
        $ref = rc($ref);
	$dna = rc($dna);
	$alt = rc($alt);
	$altDna = rc($altDna);
	$inPos = $rRds+1;
      }
     
      my $header = "$chr;$gene;$pos;$id;$alleles;$count;$dummyInfo;$inPos"; #location (1-based)
      my $newCase = ">$header;$ref;$refCnt\n$dna\n>$header;$alt;$altCnt;\n$altDna\n";
      if(!defined($fasta{$key})){
	$fasta{$key} = $newCase;
      }else{
	$fasta{$key} .= $newCase;
      }
	  #for(my $j = 0; $j<$rds; $j++){
	  #  print " ";
	  #} print $ref;
	  #print "$trimRefDna\n$trimAltDna\n";
	  #for(my $j = 0; $j<$rds; $j++){
	  #  print " ";
	  #}  print $alt;
    }
  }
  return \%fasta;
}

sub rc
{
  my ($str) = @_;
  $str = uc $str;
  # Make a new copy of the DNA (see why we saved the original?)
  my $revcom = reverse $str;
  #
  # # See the text for a discussion of tr///
  $revcom =~ tr/ACGT/TGCA/;
  return $revcom;
}


sub getShortestUtrs{
  my ($ref, $strand)  =@_;
  my %utr = %$ref;
  my %lu = ();
  for my $uS (keys %utr){
    my $le = -1;
    my $len = 99999999;
    my @ends = split(';', $utr{$uS});
    for my $ed (@ends){
      my $temp = $ed - $uS;
      if($strand eq '-'){
        $temp = 0 - $temp;
      }
      if($temp < 0){ die "ERROR: $temp <0 for UTR starting at $uS: $utr{$uS}\n"; }
      if($temp < $len){
        $len = $temp;
	$le = $ed;
      }
    }
    #print "shortest UTR: $uS to $le on $strand\n";
    $lu{$uS} = $le;
  }
  return (\%lu, $strand);
}

sub getLongestUtrs{
  my ($ref, $strand)  =@_;
  my %utr = %$ref;
  my %lu = ();
  for my $uS (keys %utr){
    my $le = -1;
    my $len = 0;
    my @ends = split(';', $utr{$uS});
    for my $ed (@ends){
      my $temp = $ed - $uS;
      if($strand eq '-'){
        $temp = 0 - $temp;
      }
      if($temp < 0){ die "ERROR: $temp <0 for UTR starting at $uS: $utr{$uS}\n"; }
      if($temp > $len){
        $len = $temp;
	$le = $ed;
      }
    }
    #print "longest UTR: $uS to $le on $strand\n";
    $lu{$uS} = $le;
  }
  return (\%lu, $strand);
}

sub getUtrByChrGene{
  my ($altEventRef, $i, $gene) = @_;
  my %altsGene = %{getAltEndsListByChr($altEventRef, $i)};
  my %alts = %{$altsGene{$gene}};
  if(defined($alts{'3+'})){ # && keys %{$alts{'3+'}} > 0){
    my %hs = %{$alts{'3+'}};
    my $no = keys %hs;
    if($no >0){
      return ($alts{'3+'}, '+'); # UTR on +
    }
  }
  return ($alts{'3-'}, '-'); # UTR on -
  
}

# procedure to parse miRanda results
#format:
#mirbase_acc    mirna_name      gene_id gene_symbol     transcript_id   ext_transcript_id       mirna_alignment alignment       gene_alignment  mirna_start     mirna_end       gene_start      gene_end        genome_coordinates      conservation    align_score     seed_cat        energy  mirsvr_score
#
# MIMAT0000062	hsa-let-7a	5270	SERPINE2	uc002vnu.2	NM_006216	uuGAUAUGUUGGAUGAU-GGAGu	  | :|: ||:|| ||| |||| 	aaCGGUGAAAUCU-CUAGCCUCu	2	21	495516	[hg19:2:224840068-224840089:-]	0.5684	122	0	-14.73	-0.7269
# get the results and arrange them in a gene manner
sub miRandaResParser{
  my ($gListRef, $miFile) = @_;
  my %genes = %$gListRef;
  my %hits = ();
  my %miRanda = (); #store results (gene name as key)

  print "Reading miRanda results: $miFile\n";
  open(FP, $miFile) or die ("ERROR: cannot open $miFile\n");
  <FP>; # ignore first line
  while(<FP>){
    chomp;
    #my ($mi_acc, $mi_na, $g_id, $g_sym, $tx_id, $extx_id, $mi_ali, $ali, $g_ali, $mi_st, $mi_end, $g_st, $g_end, $genom, $cons, $al_sc, $seed, $enrg, $misvr_sc) = split('\t', $_);
    my ($mi_acc, $mi_na, $g_id, $g_sym, $tx_id, $extx_id, $mi_ali, $ali, $g_ali, $mi_st, $mi_end, $g_st, $g_end, $genom, $cons, $al_sc, $seed, $enrg, $misvr_sc) = split('\t', $_);
    $genom = substr($genom, 1, length($genom)-2); # remove [ and ]
    my ($hg, $chrId, $startend, $strand) = split(':', $genom);
    my $key = "chr$chrId;$g_sym";
    if(defined($genes{$key})){ # in our list
      my @locations = split(',', $startend);
      for(@locations){
        my ($start, $end) = split('-', $_);

        if(!defined($hits{$key})){
          $hits{$key} = "$mi_acc:$start:$end:$strand";
        }else{
          $hits{$key} .= "\t$mi_acc:$start:$end:$strand";
        }
      }
      #info:
      $miRanda{$mi_acc} = $_; # the whole result
      #print "($mi_acc, $mi_na, $g_id, $g_sym, $tx_id, $extx_id, $mi_ali, $ali, $g_ali, $genom)\n";
    }
  }

  close(FP);

  return (\%hits, \%miRanda);

}
1;
