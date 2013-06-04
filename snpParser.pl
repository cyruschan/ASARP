#! /usr/bin/perl/ -w
use strict;
use Statistics::R; #interact with R
require 'fileParser.pl';
require 'bedHandler.pl';

use MyConstants qw( $CHRNUM $supportedList $supportedTags );

# the sub routines to process snps
# input
#	$snpFile	the snp file name (path)
#	$powCount	the threshold for powerful (frequent) snps
# optional input
#	$strandType	'+' or '-' if strand-specific flag is set
#			when input, only SNVs ending with the 
#			specified $strandType are handled.
# output
#	\%snpList	the hash reference containing both
#			powerful and non-powerful snps as well as
#			their indices (_idx).
#	\@pList		all the p-values for powerful SNVs, used to
#			get the FDR using R functions
sub initSnp{
  #snps and powSnps
  my @pList = ();
  my @snps = (); my @powSnps = ();
  for(my $i=0; $i<=$CHRNUM; $i++){
    push @snps, {};
    push @powSnps, {};
  }
  
  
  my ($snpFile, $powCount, $strandType) = @_;
  open(my $fh, "<", $snpFile) or die "Cannot open snp file: $snpFile for reading.\n";
  print "Reading from $snpFile\n";

  # Create a communication bridge with R and start R
  my $R = Statistics::R->new();
  
  my $count = 0;
  while(<$fh>){
    $count++;
    if(!$count%1000){
      #print $count." ";
    }
    chomp;
    my ($chrRaw, $pos, $alleles, $snpName, $reads, $strandInLine)=split(/ /, $_);
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
    my ($refAl, $altAl, $wrngAl) = split(/:/, $reads);

    my $chrID = getChrID($chrRaw); #auxiliary from fileparser.pl
    #check numeric
    if(!($chrID=~/^\d+$/)){
      #print "$chrID\n"; 
      next; 
    }
    #do the R Chi-squared Test
    my $totalAl = $refAl+$altAl;
    my $infoKept = $pos."\t".$alleles."\t".$snpName."\t".$refAl."\t".$altAl.";";
    #print "Ref: $refAl, Alt: $altAl, Total: $totalAl\n";
    if($totalAl>=$powCount){ #powerful snp
      
      #use the R object to make a Chi-squared goodness-of-fit test
      $R->set('x', [$refAl, $altAl]);
      $R->run('p = chisq.test(x)$p.value'); #default expected dist is c(0.5, 0.5)
      my $pValue = $R->get('p');
      push(@pList, $pValue);

      #print "Chr $chrID: Ref: $refAl, Alt: $altAl, Total: $totalAl\n";
      #print "Chi-squared test\t$pValue\n";
      
      #bionomial test (exact):
      #$R->run("p2 = binom.test($refAl, $totalAl, 0.5)");
      #my $testBin = $R->get('p2');
      #print "Binomial test\t\t";
      #my $resBin = join(' ', @$testBin);
      #if($resBin=~/p-value\s*[<|=]\s*(\S+)/){
      #  print "$1\n";
      #}

      if(defined($powSnps[$chrID]{$pos})){
	  print "Warning: multiple powerful SNVs at the same position: $pos with ".$powSnps[$chrID]{$pos}."\n";
      }	  
      $powSnps[$chrID]{$pos} .= $pValue."\t".$infoKept;	
    }elsif($totalAl>=2){ #non-trivial snps? following Gang's version
      if(defined($snps[$chrID]{$pos})){
        print "Warning: multiple SNVs at the same position: $pos with ".$snps[$chrID]{$pos}."\n";
      }
      $snps[$chrID]{$pos} .= "1\t".$infoKept;  #add 1 as the fake p-value for ordinary snps
    }

  }
  close($fh);
  $R->stop;

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

  return (\%snpList, \@pList);
}


# Reimplementation of snp2asdensity.snp2exon
# To get all snp information located by genes (see snpVSTrans for details)
# "Locating the positions of the SNV in the transcriptome"
# input
#	$snpRef		reference to the SNP list
#	$transRef	reference to the trasncript list
# optional:
#	$STRAND		the strand information: +/- or undef
# output
#	reference to the integrated gene snp list with snp info and gene locations
sub setGeneSnps{
  my ($snpRef, $transRef, $STRAND) = @_;
  my ($powSnps, $powSnps_idx, $ordSnps, $ordSnps_idx) = (undef, undef, undef, undef);
  if(!defined($STRAND)){ #non-strand specific, as usual
    ($powSnps, $powSnps_idx) = snpVsTrans($snpRef, $transRef, 'powSnps'); #powerful snps
    ($ordSnps, $ordSnps_idx) = snpVsTrans($snpRef, $transRef, 'snps'); #non-trivial snps
  }else{
    ($powSnps, $powSnps_idx) = snpVsTrans($snpRef, $transRef, 'powSnps', $STRAND); #powerful snps
    ($ordSnps, $ordSnps_idx) = snpVsTrans($snpRef, $transRef, 'snps', $STRAND); #non-trivial snps
  }
  my %geneSnps = (
   'gSnps' => $ordSnps,
   'gSnps_idx' => $ordSnps_idx,
   'gPowSnps' => $powSnps,
   'gPowSnps_idx' => $powSnps_idx,
  );

  return \%geneSnps;
}


# The corresponding core part to snp2asdensity.snp2exon
# A wrapper call setGeneSnps will call this twice using $snpTypeKey 'snps' and 'powSnps'
# to get all snp information located by genes
# "Locating the positions of the SNV in the transcriptome"
# input
#	$snpRef		reference to the SNP list
#	$transRef	reference to the trasncript list
#	$snpTypeKey	SNP type to be tested from the SNP list ('snps' or 'powSnps')
# optional
#	$setStrand	the strand specified for the SNV input
# output
#	references to the gene snp info and gene locations

sub snpVsTrans{
  my ($snpRef, $transRef, $snpTypeKey, $setStrand) = @_;

  my @geneSnps = ();
  my @geneLocations = ();
  for(my $i=0; $i<=$CHRNUM; $i++){
    push @geneSnps, {};
    push @geneLocations, {};
  }

  for(my $i=1; $i<=$CHRNUM; $i++){ #for each chromosome
    my %genes = %{$geneSnps[$i]};
    my %index = %{$geneLocations[$i]};
    my %geneMinTxStarts = (); #auxiliary to store minimal txStart for each gene


    my ($chrSnps, $chrSnps_idx) = getListByKeyChr($snpRef, $snpTypeKey, $i);
    my ($chrTrans, $chrTrans_idx) = getListByKeyChr($transRef, 'trans', $i);
    my %snps = %$chrSnps; my @snps_idx = @$chrSnps_idx; # for this chr only
    my %trans = %$chrTrans; my @trans_idx = @$chrTrans_idx; # this chr only

    if(@snps_idx==0 || @trans_idx==0){
      next;
    }
    #print "Processing "; printChr($i); print "\t";
    #print "SNVs ($snpTypeKey): ", scalar @snps_idx, "; ";
    #print "Transcripts: ", scalar @trans_idx, "\n";
    
    my ($si, $ti) = (0, 0);
    while($si<@snps_idx && $ti<@trans_idx){
      my $sPos = $snps_idx[$si];
      #print "$si VS $ti: $sPos VS $trans_idx[$ti]\n";
      if($sPos < $trans_idx[$ti]){
        #print "\$si++\n";
        $si+=1; next;
      }else{ #now the SNP is larger than some transcript start
        my $newTi = $ti; #store the old ti position
	# we have to check all $trans_idx[$newTi] until $sPos is not longer larger
        while($newTi<@trans_idx && $sPos >= $trans_idx[$newTi]){
	  my $tPos = $trans_idx[$newTi];
	  my @tSet = split('\t', $trans{$tPos});
	  my $maxTxEnd = -1; #to store the largest transcript end of @tSet
	  foreach(@tSet){
	    my ($txEnd, $cdsStart, $cdsEnd, $exonStarts, $exonEnds, $id, $gene, $isCoding, $txStrand) = split(';', $_);
	    #strand-specific handling
	    if(defined($setStrand) && $setStrand ne $txStrand){
	      next;
	    
	    }
	    if($txEnd > $maxTxEnd){	$maxTxEnd = $txEnd;	}#always store the max transcript end
	    if($sPos<=$txEnd){ # there is a hit
	       #$sPos<= $txEnd means the next $sPos may still match txEnd in @tSet, can't increase $ti.
	       # further checking
	       my @exss = split(',', $exonStarts);
	       my @exes = split(',', $exonEnds);
	       my $exNo = @exss; 
	       if($exNo != @exes){ die "exon no.s inconsistent: $exNo VS ".(scalar @exes);  }
	       
	       my $snpInfoToAdd = '';
	       for(my $j=0; $j<$exNo; $j++){
	         my $type = '';
	         if($sPos >= $exss[$j] && $sPos <= $exes[$j]){ #in exon
	           #5'UTR or 3'UTR check
	           if($j==0){
		     if($txStrand eq '+' && $sPos <$cdsStart){ #5'UTR
		       $type = '5\'UTR:+';
		     }elsif($txStrand eq '-' && $sPos <$cdsStart){ #3'UTR reverse strand
		       $type = '3\'UTR:-';
		     }else{
		       if($txStrand eq '+'){  $type = 'first:+'; }
		       else{ $type = 'last:-';  }
		     }
		   }elsif($j == $exNo-1){ #last exon start/first exon end: hv to know strand
		     if($txStrand eq '+' && $sPos >$cdsEnd){ #3'UTR
		       $type = '3\'UTR:+';
		     }
		     elsif($txStrand eq '-' && $sPos >$cdsEnd){ #5'UTR on reverse strand
		       $type = '5\'UTR:-';
		     }else{ 
		       if($txStrand eq '-'){  $type = 'first:-'; }
		       else{ $type = 'last:+';  }
		     }
		   }else{
		     $type = 'normal:'.$txStrand;
		   }
		   $snpInfoToAdd = 'exon:'.$type.';'.$gene.';'.$tPos.";".$id.';'.$exss[$j].';'.$exes[$j];
		   last;
	         }
	       
	         if($j >0 && $sPos > $exes[$j-1] && $sPos < $exss[$j]){ # in intron
	           $snpInfoToAdd = 'intron:'.$txStrand.';'.$gene.';'.$tPos.";".$id.';'.($exes[$j-1]+1).';'.($exss[$j]-1);
		   last;
	         }
               }
	       # add to the geneSnps list:
	       if(!defined($geneMinTxStarts{$gene}) || $tPos < $geneMinTxStarts{$gene}){
	         $geneMinTxStarts{$gene} = $tPos; #initial or minimal position
	       }
	       $genes{$gene} .= $sPos.";".$snpInfoToAdd."\t";
	       #print $sPos,";".$snpInfoToAdd."\n";
	       # now just store all genes at location $tPos, it will be amended by %geneMinTxStarts
	       my $geneStub = $gene."\t";
	       if(!defined($index{$tPos}) || !($index{$tPos}=~/$geneStub/)){
	         $index{$tPos} .= $gene."\t"; #gene at this exon start
	       }
             }
	  } #end of foreach(@tSet)
          if($sPos > $maxTxEnd && $newTi == $ti){ #no longer need to check this
	    $ti += 1;
	  }
	  $newTi++; #hv to check all transcript starts that are < $sPos
	} #end while
	$si += 1; #$si has gone through all $newTi (i.e. all transcript starts <= $sPos)
      }
    }
    #sort out the locations for indices
    #foreach(keys %geneMinTxStarts){
    #  print "$_: ", $geneMinTxStarts{$_}, "\n";
    #}

    foreach my $pos (keys %index){
      my $newGeneList = '';
      my @genesStartAt = split('\t', $index{$pos});
      
      foreach(@genesStartAt){
        if(defined($geneMinTxStarts{$_}) && $geneMinTxStarts{$_}==$pos){ #this indexed position is really the minimal txStart
	  $newGeneList .= $_."\t"; #re-construct the list
	}#else{ print "$_ not starting minimally at $pos but ".$geneMinTxStarts{$_}."\n"; }
      }
      if($newGeneList eq ''){ #actually all the previous genes stored here are not minimal txStart
        #delete this position in hash
	delete $index{$pos};
      }else{ #update it with the newGeneList
        $index{$pos} = $newGeneList;
      }
    } #now the index is done
 
    $geneSnps[$i] = \%genes; #gene-based arrangement
    $geneLocations[$i] = \%index; #each gene index is its first exon start position in the chromosome
  } #end of for each chromosome
  
  return (\@geneSnps, \@geneLocations);
}

################ auxiliary subroutines #######################
# get an array of genes sorted according their locations
# input: 
#	$listRef	the reference to the gene SNP list
#	$key		the key for snps: 'gSnps' or 'gPowSnps'
#	$chr		the chromosome number of interest
# output:
#	\@sortedGenes	an array containing gene names, sorted in their locations

sub getChrGeneSnpsSorted{
  my ($listRef, $key, $chr) = @_;
  my ($gListRef, $gLocRef) = getListByKeyChr($listRef, $key, $chr); 

  my %geneLocs = %{$gLocRef};
  my @sortedGenes = ();
  for(sort {$a<=>$b} keys %geneLocs){
    push @sortedGenes, $geneLocs{$_};
  }
  return \@sortedGenes;
}

# get the snp info in a list
# input: 
#	$snpInfo	the information stored for a particular individual SNP
# output:
#	@		the SNP information as an array
sub getSnpInfo{
  my ($snpInfo) = @_;
  
  #my $infoKept = $pos."\t".$alleles."\t".$snpName."\t".$refAl."\t".$altAl.";";
  #$powSnps[$chrID]{$pos} .= $pValue."\t".$infoKept;
  return split('\t', $snpInfo);
}

#################################################
sub processASEWithNev
{
  my ($snpRef, $geneSnpRef, $snpEventsNevRef, $snvPValueCutoff, $asarpPValueCutoff, $alleleRatioCutoff) = @_;
  my %ss = %$snpEventsNevRef;

  #basic statistics for powerful SNVs and genes
  my ($aseSnvCnt, $powSnvCnt, $non0AseGeneCnt, $powGeneCnt) = (0, 0, 0, 0);
  my ($aseSnvStr, $powSnvStr) = ('', '');

  my ($powAltRef, $snpAltRef, $powSpRef, $snpSpRef) = ($ss{'nevPowSnpAlt'}, $ss{'nevSnpAlt'}, $ss{'nevPowSnpSp'}, $ss{'nevSnpSp'});
  
  my @allPowAlts = @$powAltRef;
  my @allSnpAlts = @$snpAltRef;
  my @allPowSps = @$powSpRef; 
  my @allSnpSps = @$snpSpRef;
  
  # init the results
  my @aseGenes = ();
  my @asarpGenes = ();
  my @asarpSnps = ();
  my @asarpControls = ();
  for(my $i=0; $i<=$CHRNUM; $i++){
    push @aseGenes, {};
    push @asarpGenes, {};
    push @asarpSnps, {};
    push @asarpControls, {};
  }

  # Create a communication bridge with R and start R
  my $R = Statistics::R->new();
  
  for(my $i=1; $i<=$CHRNUM; $i++){
     #init
     my %aseGeneHash = ();
     my %asarpGeneHash = ();
     my %asarpGeneControls = (); # just used to store all the control SNVs we have
     my %outTabu = (); #not double-outputing the target SNVs
     my %asarpSnpHash = ();
     # new p-value scheme
     my %pValueSnpHash = (); # to store the counts (BF Correction) or the p-value lists (FDR control)

     my ($powGeneSnpChrRef) = getListByKeyChr($geneSnpRef, 'gPowSnps', $i);
     my %powGenes = %$powGeneSnpChrRef;
     my ($geneSnpChrRef) = getListByKeyChr($geneSnpRef, 'gSnps', $i);
     my %snpGenes = %$geneSnpChrRef;
     
     my ($snpChrRef) = getListByKeyChr($snpRef, 'powSnps', $i);
     my %powSnps = %$snpChrRef;
     my ($ordSnpChrRef) = getListByKeyChr($snpRef, 'snps', $i);
     my %ordSnps = %$ordSnpChrRef;


     my %powAlts = (); my %snpAlts = ();
     my %powSps = (); my %snpSps = ();

     # get ready for ASARP (AS and 5'/3' UTR alternation)
     if(defined($allPowAlts[$i])){  %powAlts = %{$allPowAlts[$i]};  }
     if(defined($allSnpAlts[$i])){  %snpAlts = %{$allSnpAlts[$i]};  }

     if(defined($allPowSps[$i])){  %powSps = %{$allPowSps[$i]};  }
     if(defined($allSnpSps[$i])){  %snpSps = %{$allSnpSps[$i]};  }
   
     if(keys %powGenes >0){  printChr($i); print"\n";  }
     #basic statistics
     $powGeneCnt += keys %powGenes;
     $powSnvCnt += keys %powSnps;
     for(keys %powSnps){
       my $snpPos = $_;
       $powSnvStr .= "$i^$snpPos\t"; #use the internal chr id to lable it to save space
       my @allSnpInfo = split(';', $powSnps{$_}); #separate by ;, if there are multiple snps at the same position
       for(@allSnpInfo){
         my ($p, $pos, $alleles, $snpId) = getSnpInfo($_);
	 if($p <= $snvPValueCutoff){
	   $aseSnvCnt += 1; #each SNV **location** added once
	   $aseSnvStr .="$i^$snpPos\t";
	   last;
	 }
       }
     }
     #gene level ASE

     # Stage 1: check gene level ASE's\n
     # Step 1.1: Genes with >= 2 Snps, and >= 1 powSnp
     for(keys %powGenes){ #genes with >= 1 powSnp
       my $gene = $_;
       #print "Gene $gene\n";
       my $snpGroupRef = groupGeneSnps($powGenes{$gene});
       my %snpGroup = %$snpGroupRef;
       #print "grouped snp:\t";
       #for(keys %snpGroup){ 
       #  print "SNP $_: "; # $snpGroup{$_}\n";
       #   my %snpHs = %{$snpGroup{$_}};
       #   for(keys %snpHs){
       #     print "$_: $snpHs{$_}\t";
       #   }
       #  print "\n";
       #} print "\n";


       my %aseList = (); #to store all the ASE snps not to be processed again
       
       my $aseCount = 0; # count of ASE snps
       my $aseInfo = '';
       for(keys %snpGroup){
         my @allSnpInfo = split(';', $powSnps{$_}); #separate by ;, if there are multiple snps at the same position
	 foreach(@allSnpInfo){
	   my ($p, $pos, $alleles, $snpId, $refCnt, $altCnt) = getSnpInfo($_);
           if($p <= $snvPValueCutoff){
	      $aseList{$pos} = 1; 
	      $aseInfo .= "$snpId,$p,$alleles,$pos,$refCnt:$altCnt\t";
	      $aseCount += 1;
	   }
	 }
       }
       # basic statistics
       if($aseCount){
	 $non0AseGeneCnt += 1;
       }

       # Step 1.2: all powSnps are ASEs and there are >=2 ASEs
       if($aseCount == keys %snpGroup && $aseCount >= 2){
         #print "Gene: $gene is a gene with all $aseCount ASE's: $aseInfo\n";
	 $aseGeneHash{$gene} = $aseInfo;

       }
       else{  
	 
	 #print "Gene: $gene is not gene-level ASE ($aseCount out of ".scalar(keys %snpGroup).")\n";
         #stage 2: ASARP: including Alternative splicing and 5'/3' Alt init/term
         #Step 2.1 get target (all snps passing NEV filter, incl. non-powerful snps) and control (non-ASE powerful snps) SNPs
	 #besids %snpGroup, we also need %ordSnpGroup
	 my $targetNev = 0; #need to get NEV to double check results

	 my %ordSnpGroup = ();
	 if(defined($snpGenes{$gene})){
           my $ordSnpGroupRef = groupGeneSnps($snpGenes{$gene});
           %ordSnpGroup = %$ordSnpGroupRef;
	 }
	 my %allSnpGroup = (%snpGroup, %ordSnpGroup); #merge all the snp groups
	 
	 for my $trgtPos (keys %allSnpGroup){ #each key is a snp
	   my ($targetFlagIT, $targetFlagS) = (0, 0); #set if it satisfies the target SNV condition
	   my ($altInit, $altTerm, $altSpInfo) = ('', '', '');
	   #print "SNV: $trgtPos\n"; 
	   # Step 2.1.1. check if this snp is with any events, i.e. in any alternatively spliced regions (5'/3' or AS)
	   # check if it is a target for AI/AT: 
	   if(defined($powAlts{$gene})){
	     ($altInit, $altTerm, $targetFlagIT) = getTargetInitTermInfo($powAlts{$gene}, $trgtPos);
	   }
	   # the first $altInit ne '' is imposed to save unnecessary time as the snp
	   # is either in %powAlts or %snpAlts
	   if($altInit eq '' && $altTerm eq '' && defined($snpAlts{$gene})){
	     ($altInit, $altTerm, $targetFlagIT) = getTargetInitTermInfo($snpAlts{$gene}, $trgtPos);
	   }

	   # check if it is a target for AS: for Splicing the format is different
	   if(defined($powSps{$gene})){
	     ($altSpInfo, $targetFlagS) = getTargetSplicingInfo($powSps{$gene}, $trgtPos); 
	   }
	   # the first $altSpInfo ne '' is imposed to save unnecessary time as the snp
	   # is either in %powSps or %snpSps
	   if($altSpInfo eq '' && defined($snpSps{$gene})){
	     ($altSpInfo, $targetFlagS) = getTargetSplicingInfo($snpSps{$gene}, $trgtPos);
	   }
           
	   if($targetFlagIT || $targetFlagS){
	     #Step 2.1.2 try to locate the control reference SNV (only from non-ASE powerful SNVs)
             for my $ctrlPos (keys %snpGroup){ #have to be powerful
	       if($ctrlPos == $trgtPos || defined($aseList{$ctrlPos})){ 
	         next; #cannot be the same pos, cannot be ASE SNP
	       }
	       # Step 2.1.3 make sure that $trgtPos and $ctrlPos are not in the same exon
	       if(areNotInSameExon(\%allSnpGroup, $trgtPos, $snpGroupRef, $ctrlPos)){
	         # a valid target-control SNV pair, now check their allele difference\n";
		 # in the current implementation, one position is assumed to have possibly multiple SNV types separated by ';', so we need to split the tuple first

                 my $targetSnpInfo = undef;
		 if(defined($powSnps{$trgtPos})){
		   $targetSnpInfo = $powSnps{$trgtPos};
		 }elsif(defined($ordSnps{$trgtPos})){
		   $targetSnpInfo = $ordSnps{$trgtPos};
		 }else{
		   die "ERROR: SNP at position $trgtPos not recorded for gene $gene at Chr $i\n";
		 }

                 my @allTargetSnps = split(';', $targetSnpInfo);
		 for(@allTargetSnps){
	           my ($tP, $tPos, $tAlleles, $tSnpId, $tAllel1, $tAllel2) = getSnpInfo($_);
                 
		   # control can only be a powerful SNV
		   my @allControlSnps = split(';', $powSnps{$ctrlPos});
		   for(@allControlSnps){
		     my ($cP, $cPos, $cAlleles, $cSnpId, $cAllel1, $cAllel2) = getSnpInfo($_);
	   	     # Step 2.2 Performing fisher's test on target: $trgtPos VS control: $ctrlPos -- $_

		     #Step 2.3 Use the R object to make a Fisher's exact test
		     #print "testing [$tAllel1, $tAllel2, $cAllel1, $cAllel2]\n";
		     $R->set('x', [$tAllel1, $tAllel2, $cAllel1, $cAllel2]);
		     $R->run('xm = matrix(data = x, nrow = 2)');
		     $R->run('p = fisher.test(xm)$p.value');
		     my $pValue = $R->get('p');
		     #print "fisher test result 1: $pValue\n";
		     
		     #print "testing [$tAllel2, $tAllel1, $cAllel1, $cAllel2]\n";
		     $R->set('x2', [$tAllel2, $tAllel1, $cAllel1, $cAllel2]);
		     $R->run('xm2 = matrix(data = x2, nrow = 2)');
		     $R->run('p2 = fisher.test(xm2)$p.value');
		     my $pValue2 = $R->get('p2');
		     #print "fisher test result 2: $pValue2\n";
		     
		     # give up new p-value scheme (too stringent for exp cases): assume each haplotype is equally probable:
		     #$pValue = ($pValue+$pValue2)/2; # 1/2*$pValue + 1/2*$pValue2
		     my $tRatio = $tAllel1/($tAllel1+$tAllel2);
		     my $cRatio = $cAllel1/($cAllel1+$cAllel2);
		     #Step 2.3 Check if the allelic ratio difference is larger than the threshold
		     my $absRatio = abs($tRatio-$cRatio);
		     my $absRatio2 = abs($tRatio-(1-$cRatio));
		     if($absRatio >= $alleleRatioCutoff || $absRatio2 >= $alleleRatioCutoff){ #allelic difference filter linked with fisher test
		       # just use multiple tests to evaluate the p-value; fisher test and ratio diff should be linked
		       my @ps = ();
		       if($absRatio >= $alleleRatioCutoff){
		         push @ps, $pValue;
			 #print "Add $gene: $trgtPos, $pValue, $absRatio\n";
		       }
		       if($absRatio2 >= $alleleRatioCutoff){
		         push @ps, $pValue2;
			 #print "Add $gene: $trgtPos, $pValue2, $absRatio2\n";
		       }
		       $pValue = 1; #re-init (not a good choice but just not to modify too much following)
		       for(@ps){
		         if($_ < $pValue){  $pValue = $_; } #just use the smaller one, if there are more than one
		       }

		       # record the number (if Bonferroni Correction is used) or all the p-values (if FDR control is used)
		       if(!defined($pValueSnpHash{$trgtPos})){
		         $pValueSnpHash{$trgtPos} = @ps; # counts for BF correction
		         #$pValueSnpHash{$trgtPos} = $pValue; # counts for BF correction
		       }else{
		         $pValueSnpHash{$trgtPos} += @ps; # counts for BF correction
		         #$pValueSnpHash{$trgtPos} .= ",$pValue"; # counts for BF correction
		       }
		       
		       #print "absolute ratio difference found: $tRatio VS $cRatio\n ASARP $gene $trgtPos found!\n";
		       if($pValue <= $asarpPValueCutoff) { #significant
		         #print "significant ($pValue) candidate pair found! $gene: $trgtPos (AI: $altInit AT: $altTerm AS: $altSpInfo) VS $ctrlPos: $powSnps{$ctrlPos}\n";	 
			 my @types = ();
			 if($altInit ne ''){ push @types, "AI:$altInit"; } #alternative 5' initiation
			 if($altTerm ne ''){  push @types, "AT:$altTerm"; } #alternative 3' termination
			 if($altSpInfo){ push @types, "AS:$altSpInfo"; } #alternative splicing
			 my $type = join(',', @types);
			 $asarpGeneControls{$gene} .= "$type;$pValue;$trgtPos;$ctrlPos\t";
			 #$asarpGeneControls{$gene} .= "$type;$trgtPos $tSnpId $tAlleles $tAllel1:$tAllel2;$ctrlPos $cSnpId $cAlleles $cAllel1:$cAllel2;$tRatio $cRatio\t"; 
			 my $snpCheck = $gene.",".$trgtPos;
			 if(!defined($outTabu{$snpCheck})){
			   $asarpGeneHash{$gene} .= "$type;$trgtPos $tSnpId $tAlleles $tAllel1:$tAllel2\t"; 
			   $outTabu{$snpCheck} = 1;
			 }
			 my $snpStub = $gene.",".$tSnpId.",".$tAlleles."\t";
			 if(!defined($asarpSnpHash{$trgtPos}) || !($asarpSnpHash{$trgtPos} =~ /$snpStub/)){
			   $asarpSnpHash{$trgtPos} .= $type.";".$snpStub;
			 }
		       }
		     }
		   }
		 }
	       }
	     }
	   }
	 }
       }
     }
     
     #print "# new p-value scheme candidate post-filtering:\n pValue candidates:\n" if(keys %pValueSnpHash > 0); 
     # first to keep necessary information only: for Bonferroni correction, this step can be skipped 
     for(keys %pValueSnpHash){
       #print "$_\t";
       if(!defined($asarpSnpHash{$_})){ #only keep those in the final target list
         delete $pValueSnpHash{$_};
       }else{
         #my @pList = split(',', $pValueSnpHash{$_});
	 #my $pValueSnpHash{$_} = fdrControl(\@pList, $asarpPValueCutoff, 0); #0/undef--not verbose
       }
     }
     #print "\n# new p-value scheme: correction\n" if(keys %pValueSnpHash > 0);
     for(keys %asarpGeneControls){
       #print "Before: $asarpGeneControls{$_}\n";
       $asarpGeneControls{$_} = pValueCorrection($asarpGeneControls{$_}, \%pValueSnpHash, $asarpPValueCutoff);
       #print "After:  $asarpGeneControls{$_}\n";
     }
     #refine the other results:
     my %newAsarpGeneHash = ();
     my %newAsarpSnpHash = ();
     my %newOutTabu = (); #new tabu list
     for(keys %asarpGeneControls){
       my $gene = $_;
       if($asarpGeneControls{$gene} eq ''){ #all the candidates have been removed
         delete $asarpGeneControls{$gene}; #delete the key as well
	 #print "All control SNVs for gene $gene deleted\n";
	 next; #nothing to be added 
       }

       #print "Getting new $gene results\n";
       #non-empty, get information for newGeneHash and newSnpHash:
       my @snps = split(/\t/, $asarpGeneControls{$_});
       my @targetGeneSnps = split(/\t/, $asarpGeneHash{$_});
       for(@snps){
         my ($type, $pValue, $trgtPos, $ctrlPos) = split(/;/, $_);
	 #print " ($type, $pValue, $trgtPos, $ctrlPos)\n";
	 
	 #now we need to get back extra info of the target SNV position ...
	 my $targetSnpInfo = undef;
	 if(defined($powSnps{$trgtPos})){
	   $targetSnpInfo = $powSnps{$trgtPos};
	 }elsif(defined($ordSnps{$trgtPos})){
	   $targetSnpInfo = $ordSnps{$trgtPos};
	 }else{
	   die "ERROR: SNP at position $trgtPos not recorded for gene $gene at Chr $i\n";
	 }
         my ($tDummyP, $tDummyPos, $tAlleles, $tSnpId, $tAllel1, $tAllel2) = (undef, undef, undef, undef, undef, undef);
	 my @allTargetSnps = split(';', $targetSnpInfo);
	 for(@allTargetSnps){
	   ($tDummyP, $tDummyPos, $tAlleles, $tSnpId, $tAllel1, $tAllel2) = getSnpInfo($_);
	   #just get the first targetSnp
	   last;
	 }

	 #print "Retrieved $trgtPos: ($tDummyP, $tDummyPos, $tAlleles, $tSnpId, $tAllel1, $tAllel2)\n";
	 
	 #re-do similar things in the processASEWithNev sub
	 #but this time the hashes are newAsarp*
	 my $snpCheck = $gene.",".$trgtPos;
	 if(!defined($newOutTabu{$snpCheck})){
	   $newAsarpGeneHash{$gene} .= "$type;$trgtPos $tSnpId $tAlleles $tAllel1:$tAllel2\t"; 
	   $newOutTabu{$snpCheck} = 1;
	 }
	 my $snpStub = $gene.",".$tSnpId.",".$tAlleles."\t";
	 if(!defined($newAsarpSnpHash{$trgtPos}) || !($newAsarpSnpHash{$trgtPos} =~ /$snpStub/)){
	   $newAsarpSnpHash{$trgtPos} .= $type.";".$snpStub;
	 }
       }
     }

     $aseGenes[$i] = \%aseGeneHash; #not affected
     $asarpGenes[$i] = \%newAsarpGeneHash; # updated by Correction
     $asarpControls[$i] = \%asarpGeneControls; # updated by Correction (in place)
     $asarpSnps[$i] = \%newAsarpSnpHash; # updated by Correction

  }
  $R->stop;

  #collect all results
  my %allAsarps = (
   'ASEgene' => \@aseGenes,
   'ASARPgene' => \@asarpGenes,
   'ASARPsnp' => \@asarpSnps,
   'ASARPcontrol' => \@asarpControls,
   'aseSnvCnt' => $aseSnvCnt,
   'powSnvCnt' => $powSnvCnt,
   'non0AseGeneCnt' => $non0AseGeneCnt,
   'powGeneCnt' => $powGeneCnt,
   'aseSnvStr' => $aseSnvStr,
   'powSnvStr' => $powSnvStr,
  );

  return \%allAsarps;
}

# sub-routine to do p-value correction
sub pValueCorrection{
  my ($snpPairs, $pControlRef, $threshold) = @_;

  my %pControl = %$pControlRef;
  #what the format is:
  #$asarpGeneControls{$gene} .= "$type;$pValue;$trgtPos;$ctrlPos\t";
  my $correctedPairs = "";

  my @snps = split(/\t/, $snpPairs);
  for(@snps){
    my ($type, $pValue, $trgtPos, $ctrlPos) = split(/;/, $_);
    #Bonferroni correction:
    #print "$trgtPos: $pValue corrected by $pControl{$trgtPos} tests = ";
    $pValue =$pValue*$pControl{$trgtPos}; # correction done
    #print "$pValue\n";
    if($pValue <= $threshold){
      $correctedPairs .= "$type;$pValue;$trgtPos;$ctrlPos\t"; # filtering
    }
    
    #FDR control correction
    #if($pValue <= $pControl{$trgtPos}){ # FDR: pControl now is the FDR corrected thresholds
    #  $correctedPairs .= "$type;$pValue;$trgtPos;$ctrlPos\t"; # filtering
    #  print "WARNING: FDR threshold $pControl{$trgtPos} < p-value $threshold\n";
    #}

  }
  return $correctedPairs;
}


# sub-routine to merge the +/- ase, asarp results for the strand-specific version

sub mergeAsarpByKey{
  my ($aRef, $aRcRef, $keyword) = @_;
  #print "keyword: $keyword\n";

  my $resRef = $aRef->{$keyword};
  my $resRcRef = $aRcRef->{$keyword};

  #init the results
  my @chrs = ();
  for(my $i=0; $i<=$CHRNUM; $i++){
    push @chrs, {};
  }
  my @resChrs = ();
  my @resRcChrs = ();
  if(defined($resRef)){
    @resChrs = @$resRef; 
  }
  if(defined($resRcRef)){
    @resRcChrs = @$resRcRef; 
  }
  for(my $i=1; $i<=$CHRNUM; $i++){
     #init
     my %asarp = ();
     my %hs = ();
     my %hsRc = ();
     if(defined($resChrs[$i])){
       %hs = %{$resChrs[$i]};
     }
     if(defined($resRcChrs[$i])){
       %hsRc = %{$resRcChrs[$i]};
     }
     # first it's empty
     for(keys %hs){
       if(0){ #$keyword ne 'ASARPsnp'){
         $asarp{$_} = $hs{$_};
       }else{
         # need to add the strand info
	 my @snps = split(/\t/, $hs{$_});
	 for(my $j = 0; $j < @snps; $j++){
	   $snps[$j] .= ';+'; 
	 }
	 $asarp{$_} = join("\t", @snps);
       }
       
     }
     # then it's cautious
     for(keys %hsRc){
       if(0){ #$keyword ne 'ASARPsnp'){ #only this uses the SNV position as the key, and thus $_ may appear in both strands
         if(!defined($asarp{$_})){
           $asarp{$_} = $hsRc{$_}; 
         }else{ # on gene level ($_) the result should be strand specific and therefore they should be exclusive
	   die "ERROR: \n+: $asarp{$_} \nexisting when \n-: $hsRc{$_} \n is to be added\n";
         }
       }else{ #ASARPsnp # need to merge +/-
         # have to split them to add the strand information (which can in fact be determined by the gene)?
	 my @snps = ();
	 if(defined($asarp{$_})){
	   @snps = split(/\t/, $asarp{$_});
	 }
	 my @snpsRc = split(/\t/, $hsRc{$_});
	 for(my $j = 0; $j < @snpsRc; $j++){
	   $snpsRc[$j] .= ';-'; 
	 }
         $asarp{$_} = join("\t", @snps, @snpsRc);
       }
     }

     $chrs[$i] = \%asarp;
  }
  return \@chrs;
}


sub mergeASARP
{
  my ($aRef, $aRcRef) = @_;
  my %as = %$aRef;
  my %rc = %$aRcRef;

  # to store merged results
  my %merged = ();
  # merge all the basic statistics (genes should be exclusive on the two strands)
  my %stats = (
   'aseSnvCnt' => 0,
   'powSnvCnt' => 0,
   'non0AseGeneCnt' => 0,
   'powGeneCnt' => 0,
  );
#   'aseSnvStr' => $aseSnvStr,
#   'powSnvStr' => $pwSnvStr,
  for(keys %stats){
    #print "$_: ";
    if(defined($as{$_})){
      $stats{$_} += $as{$_};
      #print "+: $as{$_}\t";
    }
    if(defined($rc{$_})){
      $stats{$_} += $rc{$_};
      #print "-: $rc{$_}\t";
    }
    #print "total: $stats{$_}\n";
  }
  # no need to merge SNVs at the same coordinate from 2 strands in the current version, each snv on a particular strand is unique!
  for(keys %stats){ #all basic stats have been merged now
    $merged{$_} = $stats{$_};
  }

  # merge the ASE and ASARP results
  for(qw (ASEgene ASARPgene ASARPsnp ASARPcontrol)){
    $merged{$_} =  mergeAsarpByKey($aRef, $aRcRef, $_);
    #print "merged: $merged{$_}\n";
  }

  return \%merged;
}


# sub-routine to get the detailed NEV and Alt Init/Term type information 
# from a potential target SNV position
#	input		Init and Term information (a string) and SNV position
#	output		the strings with detailed info of Init, Term respectively and the target flag
# if the target flag is non-zero (1) here, it means the SNV is a target

# reference for the format
  #$updatedEvents .= "$type,$pos,$nev,$altRegion,$constRegion\t";
  #print "$gene\t$type,$pos,$nev,$altRegion,$constRegion\n";
sub getTargetInitTermInfo{
  my ($geneInfo, $trgtPos) = @_;

  my ($altInit, $altTerm) = ('', '');
  my $targetFlag = 0;
  if($geneInfo =~ /([5|3][+|-]),$trgtPos,([\d|\.]+),/){
       $targetFlag = 1;
       my $endType = $1; #the matched parenthesis
       my $nev = $2;
       #print "nev is: $nev\n";
       if($endType eq '5+' || $endType eq '5-'){
	 $altInit = "$endType($nev)";
       }else{
	 $altTerm = "$endType($nev)";
       }
  }
  return ($altInit, $altTerm, $targetFlag);
}



# sub-routine to get the detailed NEV and (Splicing) type information 
# from a potential target SNV position
#	input		splicing information (a string) and SNV position
#	output		the string with detailed info and the target flag
# if the target flag is non-zero (1) here, it means the SNV is a target
         
# reference for the format:	 
 #print "We want this NEV: $nev, $_\n";
 #$spHash{$gene} .= join(";", $nev, $flankUsed, $snpPos, $eRegion, $lRegion, $rRegion, $strand, $additional, $tag)."\t";
sub getTargetSplicingInfo
{
  my ($geneInfo, $trgtPos) = @_;

  my $altSpInfo = '';
  my %hasType = ();
  my $targetFlag = 0; #false
  #print "sp: $powSps{$gene}\n";
  while($geneInfo =~ /([\d|\.]+);([0|1]);$trgtPos;[\d|\-|\+|;|:]+(ASS|SE|RI|UN)/g){
       $targetFlag = 1;
       my $nev = $1;
       my $isFlanking = $2;
       #print "nev is: $nev ($isFlanking)\n";
       #print "$gene: $trgtPos $1\n";
       if(!defined($hasType{$3}) || $hasType{$3} > $nev){ # just get the smallest one
         $hasType{$3} = $nev; #$isFlanking";
       }
       #print "$trgtPos matched powSps in $gene\n";
       #print "$powSps{$gene}\n";
  }
  if($targetFlag){
    for(keys %hasType){
      $altSpInfo .= "^$_($hasType{$_})";
    }
  }

  return ($altSpInfo, $targetFlag);
}

sub outputRawASARP{
  my ($allAsarpsRef, $key, $outputFile) = (undef, undef, undef);
  ($allAsarpsRef, $key, $outputFile) = @_;
 
  my $file = undef;
  if(defined($outputFile)){
   print "Output results to $outputFile\n";
   open($file, ">", $outputFile) or die "ERROR: Cannot open $outputFile to write the results\n";
   select $file; #redirect all output to the result file
  }

  my $header = '';
  my $isAsarp = 1;
  if($key eq 'ASEgene'){
    $header = "ASE gene level (all powerful SNVs are ASEs)";
    $isAsarp = 0;
  }elsif($key eq 'ASARPgene'){
    $header = "ASARP gene level";
  }elsif($key eq 'ASARPsnp'){
    $header = "ASARP snp level";
  }elsif($key eq 'ASARPcontrol'){
    $header = "ASARP control SNVs";
  }else{
    die "ERROR: Unsupported key for ASE/ASARP\n";
  }

  print $header."\n";
  my @allAsarps = @{$allAsarpsRef->{$key}};
  for(my $i=1; $i<=$CHRNUM; $i++){

    if(defined($allAsarps[$i])){
      my %chrRes = %{$allAsarps[$i]};
      for(keys %chrRes){
         print formatChr($i)."\t$_\n";
	 my @info = split('\t', $chrRes{$_});
	 foreach(@info){
	   print "$_\n";
	 }
	 print "\n";
      }
    }
  }

  if(defined($file)){
    close($file);
    select STDOUT; # back to normal
  }

}


sub outputSnvTracks{
  my ($allAsarpsRef, $outputPrefix, $db, $extraFlag) = (undef, undef, undef, undef);
  ($allAsarpsRef, $outputPrefix, $db, $extraFlag) = @_;

  if(!defined($db)){ die "ERROR: UCSC db name has to be set, e.g. hg19\n"; }
  if(!defined($extraFlag)){ $extraFlag = 0; }
 
  if(defined($outputPrefix)){
   print "Output SNV tracks to files with prefix $outputPrefix\n";
  }else{
    $outputPrefix = "snv.track";
  }
  
  my @allAsarps = @{$allAsarpsRef->{'ASARPgene'}}; #use gene.prediction
  my %catalog = (); # AS, AI, AT
  for(my $i=1; $i<=$CHRNUM; $i++){
    if(defined($allAsarps[$i])){
      my %chrRes = %{$allAsarps[$i]};
      if(keys %chrRes == 0){ next; } #skip this chromosome
      my $chr = formatChr($i);
      # determine the snp track file to write
      my $file ='';
      if(substr($outputPrefix,-1,1) eq '/'){
        $file = $outputPrefix.$chr.".pgSnp"; #no dot in between
      }else{
        $file = "$outputPrefix.$chr.pgSnp";
      }	
      my $content = '';
      my $extraCont = ''; #extra content to use
      my ($startPos, $endPos) = (0, 0); #init
      for(keys %chrRes){
	 my $gene = $_;
	 my @info = split('\t', $chrRes{$gene});
	 foreach(@info){
	   my ($catInfo, $snpInfo) = split(';', $_);
	   # chr	start		end		all	2	reads	qual
	   #chr21	31812007	31812008	T/G	2	21,70	90,70
	   my($pos, $id, $al, $reads) = split(' ', $snpInfo);
	   my ($al1, $al2) = split('>', $al);
	   my ($r1, $r2) = split(':', $reads);
	   my $snpTrack = join("\t", $chr, $pos-1, $pos, "$al1/$al2", 2, "$r1,$r2", "0,0"); #extra #, $gene, $catInfo);
	   $content .= $snpTrack."\n";
	   if($extraFlag){	$extraCont .= join("\t", $chr, $pos, $gene, $catInfo)."\n";  }

	   if($startPos ==0 || $startPos > $pos){	$startPos = $pos;	}
	   if($endPos < $pos){	$endPos = $pos;	}
	 }
      }
      if($endPos ==0 || $startPos > $endPos){
        die "Error getting $chr positions: $startPos-$endPos\n";
      }
      my $header = "track type=pgSnp visibility=3 db=$db name=\"SNV_$chr\" description=\"ASARP SNVs in $chr\"\nbrowser position $chr:$startPos-$endPos\n";
      open(FP, ">", $file) or die "ERROR: cannot open $file to output pgSnp tracks\n";
      print FP $header;
      print FP $content;
      close(FP);
      if($extraFlag){
        open(EP, ">", "$file.aux") or die "ERROR: cannot open $file.aux to output pgSnp auxiluary information\n";
	print EP $extraCont;
	close(EP);
      }
    }
  }

}


sub formatOutputVerNAR{

  my ($allAsarpsRef, $outputFile) = (undef, undef, undef);
  ($allAsarpsRef, $outputFile) = @_;
 
  my $file = undef;
  if(defined($outputFile)){
   print "Output results to $outputFile\n";
  }
  # init the output string structure
  my ($summary, $geneLevel, $snvLevel) = ("", "GENES\n", "SNVs\n");

  #basic statistics
  my ($aseSnvCnt, $powSnvCnt, $non0AseGeneCnt, $powGeneCnt) = 
  ($allAsarpsRef->{'aseSnvCnt'}, $allAsarpsRef->{'powSnvCnt'}, $allAsarpsRef->{'non0AseGeneCnt'}, $allAsarpsRef->{'powGeneCnt'});
  $summary .= "# ASE SNVs: $aseSnvCnt; # Powerful SNVs: $powSnvCnt; Percentage: ";
  if($aseSnvCnt){  $summary .= sprintf("%.1f%%", $aseSnvCnt*100/$powSnvCnt)."\n"; }
  else{ $summary .= "0%\n"; }
  $summary .= "# Genes with ASE SNVs>0: $non0AseGeneCnt; # Genes with Powerful SNVs: $powGeneCnt; Percentage: ";
  if($non0AseGeneCnt){  $summary .= sprintf("%.1f%%", $non0AseGeneCnt*100/$powGeneCnt)."\n"; }
  else{ $summary .= "0%\n"; }
  
  my @allGeneAses = @{$allAsarpsRef->{'ASEgene'}};
  my $aseCount = 0;
  for(my $i=1; $i<=$CHRNUM; $i++){

    if(defined($allGeneAses[$i])){
      my %chrRes = %{$allGeneAses[$i]};
      $aseCount += keys %chrRes;
      for(keys %chrRes){ #gene
	 #my @info = split('\t', $chrRes{$_});
	 $geneLevel .= formatChr($i)."\t".$_."\tAllASE\n";
      }
    }
  }
  if($aseCount>=2){	$summary .= "There are $aseCount genes";
  }else{		$summary .= "There is $aseCount gene"; }
  $summary .= " whose powerful SNVs (>=2) are all ASEs\n";
 
  my ($geneSumRef, $snvSumRef) = formatGeneLevelVerNAR($allAsarpsRef->{'ASARPgene'});
  my %snvHash = %$snvSumRef;
  for(keys %snvHash){
    if($_ eq 'AS'){
      $snvLevel .= "Alternative Splicing\n";
    }elsif($_ eq 'AI'){
      $snvLevel .= "Alternative Initiation\n";
    }elsif($_ eq 'AT'){
      $snvLevel .= "Alternative Termination\n";
    }else{
      die "ERROR: unknown SNV event type: $_\n";
    }
    $snvLevel .= $snvHash{$_}."\n";
  }

  my %geneHash = %$geneSumRef;
  my ($cntI, $cntS, $cntT, $cntComp) = (0, 0, 0, 0);
  my ($cntSIT, $cntSI, $cntST, $cntIT) = (0, 0, 0, 0); #for complex genes
  for my $gene (keys %geneHash){
    my ($chr, $withTypes) = split('\t', $geneHash{$gene});
    my @allEvents = split(';', $withTypes);
    my $text = '';
    if(@allEvents > 1){
      ++$cntComp;
      if($withTypes =~ 'AI;' && $withTypes =~ 'AS;' && $withTypes =~ 'AT;'){
        ++$cntSIT;
      }
      elsif($withTypes =~ 'AI;' && $withTypes =~ 'AS;'){
        ++$cntSI;
      }
      elsif($withTypes =~ 'AS;' && $withTypes =~ 'AT;'){
        ++$cntST;
      }
      elsif($withTypes =~ 'AI;' && $withTypes =~ 'AT;'){
        ++$cntIT;
      }

      $text = "Complex";
    }else{
      if($allEvents[0] eq 'AI'){
        ++$cntI;
	$text = "Initiation";
      }elsif($allEvents[0] eq 'AS'){
        ++$cntS;
	$text = "Splicing";
      }elsif($allEvents[0] eq 'AT'){
        ++$cntT;
	$text = "Termination";
      }
    }
    $geneLevel .= join("\t", $chr, $gene, $text)."\n";
  }

  $summary .= "There are $cntI 5' alternative initiation genes\n";
  $summary .= "There are $cntT 3' alternative termination genes\n";
  $summary .= "There are $cntS alternative splicing genes\n";
  
  $summary .= "There are $cntComp complex evented genes where\n".
  "  There are $cntSIT alternative splicing, initiation and termination evented genes\n".
  "  There are $cntSI alternative splicing and initiation evented genes\n".
  "  There are $cntST alternative splicing and termination evented genes\n".
  "  There are $cntIT alternative initiation and termination evented genes\n";
  my $cntTotalAsarp = $cntI + $cntT + $cntS + $cntComp;

  $summary .= "There are $cntTotalAsarp event related (non-all-ASE) genes\n\n";

  return $summary.$snvLevel.$geneLevel;
}


sub formatGeneLevelVerNAR{
  my ($inputRef) = @_; 
  my @genes = @$inputRef;
  my ($geneCnt, $typeCnt) = 0;
  my %asarps = (); #store all gene level ASEs/ASARPs
  my %geneTypes = (); #store how many different events the genes have

  for(my $i=1; $i<=$CHRNUM; $i++){

    if(defined($genes[$i])){
      my %chrRes = %{$genes[$i]};
      $geneCnt += keys %chrRes;
      for my $gene (keys %chrRes){ #each gene
        my %tabu = (); #just for this gene
	
	my @info = split('\t', $chrRes{$gene});
	foreach(@info){
	 #my ($event, $pAsarp, $target, $control) = split(';', $_);
	 my ($event, $target, $tStrand) = split(';', $_);
	 if(!defined($tStrand)){
	   $tStrand = ""; # no strand information
	 }else{
	   $tStrand = " $tStrand"; # add a space
	 }
	 my @allEvents = split(',', $event);
	 for(@allEvents){
	   my ($alt, $detail) = split(':', $_);
           if(!defined($tabu{$target.$alt})){
	     $asarps{$alt} .= $gene."\t".formatChr($i)." ". $target."$tStrand\n";
	     $tabu{$target.$alt} = $detail;
	   }
	   my $stub = $alt.";";
	   if(!defined($geneTypes{$gene})){
	     $geneTypes{$gene} = formatChr($i)."\t";
	   }
	   if(!($geneTypes{$gene} =~ $stub)){
	     $geneTypes{$gene} .= $stub; #add this new alt type
	   }
	 }
	}
      }

    }
  }
  return (\%geneTypes, \%asarps);
}


sub areNotInSameExon
{
  my ($targetRef, $target, $controlRef, $control, $chrTransRef) = @_;
  my $targetSnpRef = $targetRef->{$target};
  my $controlSnpRef = $targetRef->{$control};
  my %targetSnp = %$targetSnpRef;
  my %controlSnp = %$controlSnpRef;
  
  # Step 1: a quick check with just the intersect exons (pre-screening)
  if(defined($targetSnp{'intron+'}) || defined($targetSnp{'intron-'}) ||defined($controlSnp{'intron+'}) || defined($controlSnp{'intron-'})){
    #print "$target $control one of the 2 snps is in intron\n";
    return 1;
  }
 
  # Step 2: if they are only in exons, have to also check whether there is any transcript in which one is absent
  for my $tag ('exon+', 'exon-'){ # same intron is also considered "NotInSameExon", not need for 'intron+', 'intron-'){

    if(defined($targetSnp{$tag}) && defined($controlSnp{$tag})){
      my ($tS, $tE, $tTransIds) = split(';', $targetSnp{$tag});
      my ($cS, $cE, $cTransIds) = split(';', $controlSnp{$tag});
      my $nonOverlapFlag = 1;
      #print "$target: $targetSnp{$tag}\n";
      #print "$control: $controlSnp{$tag}\n";
      if($tS > $cE || $cS > $tE){ #exon overlaps
	#print "$target not overlaps $control $tag: $tS-$tE $cS-$cE\n";
        return 1;
      }
      
      #check if they are from different transcripts
      my @tIds = split(',', $tTransIds);
      my @cIds = split(',', $cTransIds);
      my $tIdSize = @tIds;
      my $cIdSize = @cIds;
      if($tIdSize != $cIdSize){
        #print "$target $control different transcript IDs: $tIdSize != $cIdSize\n";
	return 1;
      }else{
        @tIds = sort @tIds;
	@cIds = sort @cIds;
	for(my $q = 0; $q < $tIdSize; $q++){
	  if($tIds[$q] ne $cIds[$q]){
	    #print "$target $control different at $tIds[$q] $cIds[$q]\n";
	    return 1;
	  }
	}
	#they are always in the same exon (have to go through + and -)
      }
    }else{ 
      if(defined($targetSnp{$tag}) || defined($controlSnp{$tag})){
        #print "$target and $control one not defined with $tag\n";
        return 1; 
      }
    }
  }
  #print "$target and $control do overlap\n";
  return 0;
}



#################################################
sub filterSnpEventsWithNev
{
  my ($snpRef, $geneSnpRef, $snpEventsRef, $bedF, $spEventsListRef, $nevCutoffLower, $nevCutoffUpper, $STRAND) = @_;
  
  # Preparation of all the splicing and 5'/3' alt init/term events from $snpEventsRef
  my %ss = %$snpEventsRef;
  my ($powAltRef, $snpAltRef, $powSpRef, $snpSpRef) = ($ss{'powSnpAlt'}, $ss{'snpAlt'}, $ss{'powSnpSp'}, $ss{'snpSp'});
  
  my @allPowAlts = @$powAltRef;
  my @allSnpAlts = @$snpAltRef;
  my @allPowSps = @$powSpRef; 
  my @allSnpSps = @$snpSpRef;
 

  # init of the result components
  my @nevPowAlts = ();
  my @nevSnpAlts = ();
  my @nevPowSps = ();
  my @nevSnpSps = ();

  print "Filtering splicing events (alternative splicing, 5'/3' alternations) based on NEV's.\n";
  for(my $i=1; $i<=$CHRNUM; $i++){
     #init
     # get ready for ASARP (splicing and 5'/3' alt init/term events)
     my %powAlts = %{$allPowAlts[$i]};
     my %snpAlts = %{$allSnpAlts[$i]};

     my %powSps = %{$allPowSps[$i]};
     my %snpSps = %{$allSnpSps[$i]};
   
     my $powAltCnt = keys %powAlts; my $snpAltCnt = keys %snpAlts;
     my $powSpCnt = keys %powSps;   my $snpSpCnt = keys %snpSps;
     # read bed file of this chromosome when it is needed
     my $bedRef = undef;
     if($powAltCnt  > 0 || $snpAltCnt > 0 || $powSpCnt > 0 || $snpSpCnt > 0){ #it only makes sense when snps are there
        #readBed
        printChr($i); print "\n";
	$bedRef = readBedByChr($bedF, $i);
     }

     #print "# for 5'/3' alt init/term events\n";
     #update (shortlist) the alt events with NEV values calculated from bed information
     if($powAltCnt  > 0){	$nevPowAlts[$i] = calAltEventNev(\%powAlts, $bedRef, $i, $nevCutoffLower, $nevCutoffUpper, $STRAND);	}
     if($snpAltCnt  > 0){ 	$nevSnpAlts[$i] = calAltEventNev(\%snpAlts, $bedRef, $i, $nevCutoffLower, $nevCutoffUpper, $STRAND);	}

     #print "#for splicing events NEV calculation\n";
     if($powSpCnt  > 0){       $nevPowSps[$i] = calSplicingEventNev(\%powSps, $bedRef, $i, $spEventsListRef, $geneSnpRef, 'gPowSnps', $nevCutoffLower, $nevCutoffUpper, $STRAND);	} 
     if($snpSpCnt  > 0){       $nevSnpSps[$i] = calSplicingEventNev(\%snpSps, $bedRef, $i, $spEventsListRef, $geneSnpRef, 'gSnps', $nevCutoffLower, $nevCutoffUpper, $STRAND); 	}
  }

  my %snpEventsNev = (
    'nevPowSnpAlt' => \@nevPowAlts,
    'nevSnpAlt' => \@nevSnpAlts,
    'nevPowSnpSp' => \@nevPowSps,
    'nevSnpSp' => \@nevSnpSps,
  );
  return \%snpEventsNev;

}


# group the gene's snps
# the input $geneSnps is the value containing all snps under the gene key (a string not yet split)
sub groupGeneSnps
{
  my ($geneSnps) = @_;
  my %groups = ();
  
  #print "Group GeneSnps $geneSnps\n\n";
  # side-track: get all Snps first
  my %snps = ();
  my %snpsNeg = ();
  
  my %snpsIn = ();
  my %snpsNegIn = ();


  # now you can refine the exons that contain this snp (getting the minimal exon overlap)
  my %snpExonStarts = (); #get the intersect exon start for each snp
  my %snpExonEnds = (); #get the interset exon end for each snp
  my %snpIntronStarts = ();
  my %snpIntronEnds = ();

  my %snpNegExonStarts = ();
  my %snpNegExonEnds = ();
  my %snpNegIntronStarts = ();
  my %snpNegIntronEnds = ();
  my @allSnps = split('\t', $geneSnps);
  foreach(@allSnps){
    my ($snpPos, $exonIntronType, $geneName, $txStart, $transId, $exonS, $exonE) =  split(';', $_);
    my @typeInfo = split(':', $exonIntronType); 
    if($typeInfo[-1] eq '+'){ #forward strand
      if($typeInfo[0] eq 'exon'){
        $snps{$snpPos} .= $transId.','; #defined now
        if(!defined($snpExonStarts{$snpPos}) || $snpExonStarts{$snpPos} < $exonS){
          $snpExonStarts{$snpPos} = $exonS;
        }
        if(!defined($snpExonEnds{$snpPos}) || $snpExonEnds{$snpPos} > $exonE){
          $snpExonEnds{$snpPos} = $exonE;
        }
      }
      else{ # intron type
        $snpsIn{$snpPos} .= $transId.','; #defined now
        if(!defined($snpIntronStarts{$snpPos}) || $snpIntronStarts{$snpPos} < $exonS){
          $snpIntronStarts{$snpPos} = $exonS;
        }
        if(!defined($snpIntronEnds{$snpPos}) || $snpIntronEnds{$snpPos} > $exonE){
          $snpIntronEnds{$snpPos} = $exonE;
        }
      }
    }else{ #strand is -
      if($typeInfo[0] eq 'exon'){
	$snpsNeg{$snpPos} .= $transId.','; #defined now
        if(!defined($snpNegExonStarts{$snpPos}) || $snpNegExonStarts{$snpPos} < $exonS){
          $snpNegExonStarts{$snpPos} = $exonS;
        }
        if(!defined($snpNegExonEnds{$snpPos}) || $snpNegExonEnds{$snpPos} > $exonE){
          $snpNegExonEnds{$snpPos} = $exonE;
        }
      }
      else{ # intron type
        $snpsNegIn{$snpPos} .= $transId.','; #defined now
        if(!defined($snpNegIntronStarts{$snpPos}) || $snpNegIntronStarts{$snpPos} < $exonS){
          $snpNegIntronStarts{$snpPos} = $exonS;
        }
        if(!defined($snpNegIntronEnds{$snpPos}) || $snpNegIntronEnds{$snpPos} > $exonE){
          $snpNegIntronEnds{$snpPos} = $exonE;
        }
      }
    }
  }

  for(keys %snps){
    if(!defined($groups{$_})){
      $groups{$_} = {}; #empty hash
    }
    my %hash = %{$groups{$_}};
    if(defined($snpExonStarts{$_})){
      $hash{'exon+'} = join(';', $snpExonStarts{$_}, $snpExonEnds{$_}, $snps{$_});
    }
    $groups{$_} = \%hash; #update the reference
  }

  for(keys %snpsIn){
    if(!defined($groups{$_})){
      $groups{$_} = {}; #empty hash
    }
    my %hash = %{$groups{$_}};
    if(defined($snpIntronStarts{$_})){
      $hash{'intron+'} = join(';', $snpIntronStarts{$_}, $snpIntronEnds{$_}, $snpsIn{$_});
    }
    $groups{$_} = \%hash; #update the reference
  }

  for(keys %snpsNeg){
    if(!defined($groups{$_})){
      $groups{$_} = {}; #empty hash
    }
    my %hash = %{$groups{$_}};
    if(defined($snpNegExonStarts{$_})){
      $hash{'exon-'} = join(';', $snpNegExonStarts{$_}, $snpNegExonEnds{$_}, $snpsNeg{$_});
    }
    $groups{$_} = \%hash;
  }
  
  for(keys %snpsNegIn){
    if(!defined($groups{$_})){
      $groups{$_} = {}; #empty hash
    }
    my %hash = %{$groups{$_}};
    if(defined($snpNegIntronStarts{$_})){
      $hash{'intron-'} = join(';', $snpNegIntronStarts{$_}, $snpNegIntronEnds{$_}, $snpsNegIn{$_});
    }
    $groups{$_} = \%hash;
  }

  return \%groups;
}

sub getGroupedSnpInfoByType{
  my ($groupsRef, $pos, $type) = @_;
  #print  "Getting SNV $pos for type $type\n";
  my %groups = %$groupsRef;
  if(defined($groups{$pos})){
    my %snps = %{$groups{$pos}};
    if(defined($snps{$type})){
       return $snps{$type};
    }
  }
  return ''; #empty info

}

sub calSplicingEventNev
{
  #print "Splicing\n";
  my ($spsRef, $bedRef, $chr, $spEventsListRef, $gSnpRef, $gSnpKey, $nevCutoffLower, $nevCutoffUpper, $setStrand) = @_;
  
  # gene level Snps information for splicing NEV calculation 
  my ($geneSnpChrRef) = getListByKeyChr($gSnpRef, $gSnpKey, $chr);
  my %geneSnps = %$geneSnpChrRef;
  my %spHash = (); #store the updated results  

  my %spEventsList = %$spEventsListRef;
  my %spConstExons = ();

  foreach(keys %spEventsList){
    if(defined($spEventsList{$_}->{'type'})){
      my $tag = $spEventsList{$_}->{'type'};
      if(checkSupportedType($tag)){
        #print "Loading constitutive exons for gene $_ with $tag\n";
        $spConstExons{$tag} = getConstitutiveExonsByChr($spEventsList{$_}, $chr);
      }
    }
    else{ die "ERROR parsing event: no type available for $_\n"; }
  }

  my %allGeneSpSnps = %$spsRef;
  for(keys %allGeneSpSnps){
    my $gene = $_;
    
    my %geneConstRatio = (
     'anno' => undef,
     'rna' => undef,
     'est' => undef,
    );
    #pre-calculate all the geneConstRatios
    for my $t (keys %spConstExons){
      if(defined($spConstExons{$t})){
        my %chrCE = %{$spConstExons{$t}};

        if(defined($chrCE{$gene})){
          my $constExonSet = $chrCE{$gene};   
          #print "$gene ";
	  $geneConstRatio{$t} = calConstRatio($constExonSet, $bedRef, $setStrand);
	  #print "$t: $gene const ratio: $geneConstRatio{$t}\n";

        } #get the constitutive ratio if there is event evidence for this gene
      }
    }
    my $groupsRef = groupGeneSnps($geneSnps{$gene});

    my @allEvents = split('\t', $allGeneSpSnps{$gene});
    foreach(@allEvents){
     my ($snpPos, $eRegion, $lRegion, $rRegion, $strand, $additional, $tag) = split(';', $_);
     my $groupSnpInfo = getGroupedSnpInfoByType($groupsRef, $snpPos,"exon".$strand);
     if($groupSnpInfo eq ''){ #no exon info
       $groupSnpInfo = getGroupedSnpInfoByType($groupsRef, $snpPos,"intron".$strand);
     }

     if($groupSnpInfo ne ''){ #there is information
       # refine (intersect) the event region $eRegion as needed
       my ($s, $e, $transIds) = split(';', $groupSnpInfo);
       my ($eStart, $eEnd) = split(':', $eRegion);
       
       # do intersection if they overlap
       if($eStart <= $e && $eEnd >= $s){
         #print "Intersect! $eStart,$eEnd with $s, $e\n";
	 if($eStart < $s){ $eStart = $s; }
	 if($eEnd > $e){ $eEnd = $e; }
       }
       #print "Effective length: ".($eEnd-$eStart+1)."\n";
       if(!defined($geneConstRatio{$tag})){
         #print "Warning: no const ratio for $tag for $gene\n $allGeneSpSnps{$gene}\n";
         last;
       }
       my ($nev, $flankUsed) = calSpNev($eStart, $eEnd, $lRegion, $rRegion, $bedRef, $geneConstRatio{$tag}, $setStrand); 
       $nev = sprintf("%.4f", $nev);
       if($nev > $nevCutoffLower && $nev < $nevCutoffUpper){
         #print "We want this NEV: $nev, $_\n";
	 $spHash{$gene} .= join(";", $nev, $flankUsed, $snpPos, $eRegion, $lRegion, $rRegion, $strand, $additional, $tag)."\t";
       }else{
         #print "We dont want this NEV: $nev, $_\n";
       }
     }else{
       #print "WARNING: SNP $snpPos $tag event does not match $gene: $_\n";
       #exit;
     }
    }
  }
  return \%spHash;
}

sub calConstRatio
{
  my ($constExonSet, $bedRef, $setStrand) = @_;
  my @allConstExons = split(';', $constExonSet);
  my ($readCount, $effLen) = (0, 0);
  foreach(@allConstExons){
    my ($s, $e) = split('-', $_);
    my ($r, $l) = getEffReadSumLength($bedRef, $s, $e, $setStrand);
    #print "$_: reads: $r, length: $l\n";
    $readCount += $r;
    $effLen += $l;
  }
  if($effLen == 0){ 
    #print "Warning: a gene without constExonSet reads: $constExonSet\n"; 
    return 0;
  }
 
  #print "const read: $readCount const len: $effLen\n";

  return $readCount/$effLen;
}

sub calSpNev
{
  my ($es, $ee, $lRegion, $rRegion, $bedRef, $constExonRatio, $setStrand) = @_;
  my ($ls, $le) = split(':', $lRegion);
  my ($rs, $re) = split(':', $rRegion);
  my ($nev, $nev2, $altReads, $altLen) = (-1, -1, 0, 0);
  my $isFlanking = 0; #check if flanking region is used  
  # possible alt regions
  ($altReads, $altLen) = getEffReadSumLength($bedRef, $es, $ee, $setStrand);

  #use the event only to calculate Nev
  #cal psi based on consituent exon  
  if($constExonRatio > 0 && $altLen >0){
    $nev = ($altReads/$altLen)/$constExonRatio;
  }
  
  #calculate the ratio based on the flanking regions
  if($ls != -1 && $rs !=-1){
    # flanking regions
    my ($c, $l) = getEffReadSumLength($bedRef, $ls, $le, $setStrand); 
    my ($c2, $l2) = getEffReadSumLength($bedRef, $rs, $re, $setStrand); 
    $c += $c2; $l += $l2;
    if($l > 0 && $altLen >0){
      $nev2 = ($altReads/$altLen)/($c/$l);
    }

  }
 
  #get the smaller one (but has to >0)
  if($nev2 > 0 && $nev2 < $nev){
   $isFlanking = 1; # flanking used
   $nev = $nev2;
  }
  return ($nev, $isFlanking);
}

sub calAltEventNev
{
  my ($powAltsRef, $bedRef, $chr, $nevCutoffLower, $nevCutoffUpper, $setStrand) = @_;
  my %allGeneSnps = %$powAltsRef;
  for(keys %allGeneSnps){
     my $gene = $_;
     my $updatedEvents = '';
     my @allEvents = split('\t', $allGeneSnps{$_});
     foreach(@allEvents){
       my ($type,$pos,$altRegion,$constRegion) = split(',', $_);
       my ($altL, $altR) = split('-', $altRegion);
       my ($constL, $constR) = split('-', $constRegion);
       my ($altRead, $altLen) = getEffReadSumLength($bedRef, $altL, $altR, $setStrand);
       my ($constRead, $constLen) = getEffReadSumLength($bedRef, $constL, $constR, $setStrand);
       my ($altRatio, $constRatio) = (0, 0);
       if($altLen >0){
         $altRatio = $altRead/$altLen;
       }
       if($constLen > 0){
         $constRatio = $constRead/$constLen;
       }
       my $nev = 0; #normalized expression value, similar to psi
       if($constRatio>0){
         $nev = $altRatio/$constRatio;
       }
       $nev = sprintf("%.4f", $nev);
       #print "$gene $pos NEV $nev\n";
       if($nev > $nevCutoffLower && $nev < $nevCutoffUpper){
          $updatedEvents .= "$type,$pos,$nev,$altRegion,$constRegion\t";
	  #print "$gene\t$type,$pos,$nev,$altRegion,$constRegion\n";
       }

     }
     $allGeneSnps{$_} = $updatedEvents;
  }
  return \%allGeneSnps;
}

#################################################
sub setSnpEvents{
  my ($geneSnpRef, $altRef, $splicingRef) = @_;
  my ($gPowSnpAltRef, $gSnpAltRef) = setSnpAltEvents($geneSnpRef, $altRef); #match snps with Alt Events
  my ($gPowSnpSpRef, $gSnpSpRef) = setSnpSplicingEvents($geneSnpRef, $splicingRef); #match snps with Splicing Events
  
  my %snpEvents = (
    'powSnpAlt' => $gPowSnpAltRef,
    'snpAlt' => $gSnpAltRef,
    'powSnpSp' => $gPowSnpSpRef,
    'snpSp' => $gSnpSpRef,
  );
  return \%snpEvents;
}
#################################################
## the following sub is to match the Snps with the Alt Events
sub setSnpAltEvents{
  my ($geneSnpRef, $altRef) = @_;

  my @snpEvents = ();
  my @powSnpEvents = ();
  for(my $i=0; $i<=$CHRNUM; $i++){
    push @snpEvents, {};
    push @powSnpEvents, {};
  }
  my $gPowSnpAltEndsRef = \@powSnpEvents;
  my $gSnpAltEndsRef = \@snpEvents;

  ## attributes for each gene in a chromosome of the snpEvents
  # 'AltEnd', 'AltSp', 'GeneLv' = 1, or 0
  $gPowSnpAltEndsRef = snpVsAltEvents($gPowSnpAltEndsRef, $geneSnpRef, 'gPowSnps', $altRef);
  $gSnpAltEndsRef = snpVsAltEvents($gSnpAltEndsRef, $geneSnpRef, 'gSnps', $altRef);

  return ($gPowSnpAltEndsRef, $gSnpAltEndsRef);
}

# matching snps of a particular type (gSnps or gPowSnps) with Alt Ends
sub snpVsAltEvents
{
  my ($snpEventRef, $geneSnpRef, $geneSnpKey, $altEventRef) = @_;
  #print "Gene level SNP ($geneSnpKey) VS Alternative transcript ends (5', 3')\n";
  my @snpEvents = @$snpEventRef;

  for (my $i=1; $i<=$CHRNUM; $i++){
    my ($geneRef) = getListByKeyChr($geneSnpRef, $geneSnpKey, $i); 
    my $geneChrRef = getChrGeneSnpsSorted($geneSnpRef, $geneSnpKey, $i);
    my @chrGenes = @$geneChrRef;

    # Alt 5', 3' transcript ends information
    my %alts = %{getAltEndsListByChr($altEventRef, $i)};
    
    #my %genesPrinted = (); #to store genes printed already in order not to double-print
    if(@chrGenes){
      #printChr($i); print "\n";
      foreach(@chrGenes){
        my @allGenes = split('\t', $_);
        foreach(@allGenes){
	  #get the snp information
	  
          my %snpPosTabu = (); # to avoid duplicate processing of snps for this particular gene
	  my @geneMatches = split('\t', $geneRef->{$_});
	  for(@geneMatches){
	    my ($snpPos, $matchInfo, $geneName, $txStart, $id, $regStart, $regEnd) = split(';', $_);
            if(defined($snpPosTabu{$snpPos})){
	      next; # skip all the following processing
	    }
	    $snpPosTabu{$snpPos}=1; #set the tabu list
	    # just check for printing debug
	    #if(!defined($genesPrinted{$geneName})){
	    #  print $geneName."\n";
	    #  $genesPrinted{$geneName} = 1;
	    #}else{ $genesPrinted{$geneName} += 1; }
            
	    my $snpAltInfoAdd = '';
            if(defined($alts{$geneName})){ #only when there are events defined for this gene
	      $snpAltInfoAdd = matchSnpPoswithAltEvents($snpPos, $alts{$geneName}); 
	    }
	    if($snpAltInfoAdd ne ''){
	       #print "$snpAltInfoAdd\n";
	       my %snpAltHash = ();
               if(defined($snpEvents[$i])){
	         # there is information originally, use it
		 %snpAltHash = %{$snpEvents[$i]};
               }
	       if(!defined($snpAltHash{$geneName})){
	         $snpAltHash{$geneName} = $snpAltInfoAdd;
	       }else{
	         $snpAltHash{$geneName} .= $snpAltInfoAdd;
	       }
	       $snpEvents[$i] = \%snpAltHash; #get back to the reference
	    }
	  } #endof foreach(@geneMatches) 
        } #endof foreach (@allGenes)
      } #endof foreach(@chrGenes)
    }
  } #end of chromosome

  return \@snpEvents;
}

sub matchSnpPoswithAltEvents
{
  my ($pos, $altRef) = @_;
  my %alts = %$altRef;
  my $snpInfoToAdd = '';

  my @allTypes = ('5+', '5-', '3+', '3-');
  for my $type (@allTypes){
    if(defined($alts{$type})){ #exists
      my %hash = %{$alts{$type}};
      for my $altKeyPos (keys %hash){
	
	if($type eq '5+' || $type eq '3-'){ #key is larger than value
          if($altKeyPos >= $pos){ #need to check
            my @altPos = split(';', $hash{$altKeyPos});
	    @altPos = sort { $a <=> $b } @altPos; #sorted in ascending order
	    if(@altPos > 1){ # there are alternative splicing
              my $pi = 0;
	      while($pi<@altPos-1){
	        if($altPos[$pi] <= $pos && $pos <= $altPos[$pi+1]-1){ # a hit
		  #store both alt and const regions
		  $snpInfoToAdd .= "$type,$pos,$altPos[$pi]-".($altPos[$pi+1]-1).",$altPos[-1]-$altKeyPos\t";
		  #print "$snpInfoToAdd as $pos hit key $altKeyPos: @altPos\n";
		  last;
		}
	        $pi++;
	      }

	    } # if @altPos <= 1, no alternative start/end for 5' or 3'
	  }
	
	}else{ #'5-' or '3+': key is smaller than value
          if($altKeyPos <= $pos){ #need to check
            my @altPos = split(';', $hash{$altKeyPos});
	    @altPos = sort { $a <=> $b } @altPos; #sorted in ascending order
	    if(@altPos > 1){ # there are alternative splicing
              my $pi = 1;
	      while($pi<@altPos){
	        if($altPos[$pi-1]+1 <= $pos && $pos <= $altPos[$pi]){ # a hit
		  $snpInfoToAdd .= "$type,$pos,".($altPos[$pi-1]+1)."-$altPos[$pi],$altKeyPos-$altPos[0]\t";
		  #print "$snpInfoToAdd as $pos hit key $altKeyPos: @altPos\n";
		  last;
		}
	        $pi++;
	      }
	    } # if @altPos <= 1, no alternative start/end for 5' or 3'
	  }
	
	}
      }
    }
  }
  return $snpInfoToAdd;
}

#################################################
## the following sub is to match the Snps with the Splicing Events

sub setSnpSplicingEvents{
  my ($geneSnpRef, $splicingRef) = @_;

  my @snpEvents = ();
  my @powSnpEvents = ();
  for(my $i=0; $i<=$CHRNUM; $i++){
    push @snpEvents, {};
    push @powSnpEvents, {};
  }
  my $gPowSnpSplicingRef = \@powSnpEvents;
  my $gSnpSplicingRef = \@snpEvents;

  ## attributes for each gene in a chromosome of the snpEvents
  # 'AltEnd', 'AltSp', 'GeneLv' = 1, or 0
  $gPowSnpSplicingRef = snpVsSplicingEvents($gPowSnpSplicingRef, $geneSnpRef, 'gPowSnps', $splicingRef);
  $gSnpSplicingRef = snpVsSplicingEvents($gSnpSplicingRef, $geneSnpRef, 'gSnps', $splicingRef);

  return ($gPowSnpSplicingRef, $gSnpSplicingRef);
}


# matching snps of a particular type (gSnps or gPowSnps) with internal exon splicing events
sub snpVsSplicingEvents
{
  my ($snpEventRef, $geneSnpRef, $geneSnpKey, $SplicingEventsRef) = @_;
  #print "Gene level SNP ($geneSnpKey) VS Splicing events\n";
  my @snpEvents = @$snpEventRef;

  for (my $i=1; $i<=$CHRNUM; $i++){
    my ($geneRef) = getListByKeyChr($geneSnpRef, $geneSnpKey, $i); 
    my $geneChrRef = getChrGeneSnpsSorted($geneSnpRef, $geneSnpKey, $i);
    my @chrGenes = @$geneChrRef;
   
    my ($spRef, $spIdx) = getListByKeyChr($SplicingEventsRef, 'events', $i);
    my %splicing = %$spRef; 
    #my %genesPrinted = (); #to store genes printed already in order not to double-print

    if(@chrGenes){
      #printChr($i); print "\n";
      foreach(@chrGenes){
        my @allGenes = split('\t', $_);
        foreach(@allGenes){ #each element is a gene name
	  #get the snp information
	    
          my %snpPosTabu = (); # to avoid duplicate processing of snps for this particular gene
	  my @geneMatches = split('\t', $geneRef->{$_});
	  for(@geneMatches){
	    my ($snpPos, $matchInfo, $geneName, $txStart, $id, $regStart, $regEnd) = split(';', $_);
            if(defined($snpPosTabu{$snpPos})){
	      next; # skip all the following processing
	    }
	    $snpPosTabu{$snpPos}=1; #set the tabu list
	    # just check for printing debug
	    #if(!defined($genesPrinted{$geneName})){
	    #  print $geneName."\n";
	    #  $genesPrinted{$geneName} = 1;
	    #}else{ $genesPrinted{$geneName} += 1; }
	    
	    my $snpSpInfoAdd = '';
            if(defined($splicing{$geneName})){ #only when there are events defined for this gene
	      $snpSpInfoAdd = matchSnpPoswithSplicingEvents($snpPos, $splicing{$geneName}); 
	    }
	    if($snpSpInfoAdd ne ''){
	       #print "$snpSpInfoAdd\n";
	       my %snpSpHash = ();
               if(defined($snpEvents[$i])){
	         # there is information originally, use it
		 %snpSpHash = %{$snpEvents[$i]};
               }
	       if(!defined($snpSpHash{$geneName})){
	         $snpSpHash{$geneName} = $snpSpInfoAdd;
	       }else{
	         $snpSpHash{$geneName} .= $snpSpInfoAdd;
	       }
	       $snpEvents[$i] = \%snpSpHash; #get back to the reference
	    }
	  } #endof foreach(@geneMatches) 
        } #endof foreach (@allGenes)
      } #endof foreach(@chrGenes)
    }
  } #end of chromosome

  return \@snpEvents;
}

sub matchSnpPoswithSplicingEvents
{
  my ($pos, $splicing) = @_;
  my $snpInfoToAdd = '';
  my @events = split('\t', $splicing);
  foreach(@events){
    my ($eRegion, $lRegion, $rRegion, $strand, $additional, $tag) = split(';', $_);
    my ($eStart, $eEnd) = split(':', $eRegion);
    if($eStart <= $pos && $pos <= $eEnd){ # a snp match!
      $snpInfoToAdd .= $pos.";".$_."\t";
      #print "Add event: $_ [$eStart, $eEnd] to $pos\n";
    }
    
  }
  return $snpInfoToAdd;
}

# FDR control using the ASE SNVs (i.e. those powerful SNVs with Chi Square Test issued)
# The FDR control is a modified version from BH's method
# The original implementation is in Matlab and now it's done using R and Perl
# the reference: Controlling the proportion of falsely-rejected hypotheses
# when conducting multiple tests with climatological data
#	input:		reference to the p-value list, FDR threshold
#	optional input:	whether verbose print outs are enabled
#	output:		the p-value threshold for the FDR threshold
sub fdrControl{
  my ($pRef, $fdrCutoff, $isVerbose) = @_;
  my $orgFdrCutoff = $fdrCutoff; # fall-back plan when the adjusted FDR fails
  if(!defined($isVerbose)){ $isVerbose = 0;   }
  my @pList = @$pRef;
  my $pListSize = @pList; #size
  @pList = sort{$a<=>$b}@pList; #sorted
  my $pSize = @pList;
  # Create a communication bridge with R and start R
  my $R = Statistics::R->new();

  #$R->set('plist', \@pList);
  my $stepSize = 10000;
  if($pSize > $stepSize){ # huge input
    $R->run("plist <- c()"); #initialize
    for(my $xi = 0; $xi < $pSize; $xi+=$stepSize){ #each time
      my $end = $xi+$stepSize-1;
      if($end >= $pSize){ $end = $pSize-1; }  #the upper bound
      #print "getting plist: $xi to $end\n";
      $R->run("plist <- c(plist, c(".join(",", @pList[$xi .. $end])."))");
    }
  }else{
    $R->run("plist <- c(".join(",",@pList).")");
  }
  #the expected proportion of false discoveries amongst the rejected hypotheses
  #http://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html
  #print "Running R using BY\n";
  $R->run('x <- p.adjust(plist, method="BY")');
  $R->run('rLen <- length(x)');
  my $rSize = $R->get('rLen');
  #print "Getting x from R: size $rSize\n";
  
  my @pAdjust = ();
  if($rSize > $stepSize){
    # need to use the same trick to get parts from Statistics::R
    for(my $xi = 1; $xi <= $pSize; $xi+=$stepSize){ #each time: one-based in R
      #1-based in R!
      my $end = $xi+$stepSize-1;
      if($end > $pSize){ $end = $pSize; }  #the upper bound

      #print "getting x: $xi to $end\n";
      $R->run("xSlice <- x[$xi:$end]");
      my $pSliceRef = $R->get('xSlice');
      push(@pAdjust, @$pSliceRef);
    }
  }else{
    my $pAdjustRef = $R->get('x');
    @pAdjust = @$pAdjustRef;
  }
  $R->stop;
  if($rSize != $pSize || @pAdjust != $pSize){
    print "ERROR: Statistics::R result size: ".(scalar @pAdjust).", or R result size: $pSize different from input size: $pSize\n";
    die "ERROR: the adjusted p-value list from R is inconsistent (see STDOUT for details)! Aborted\n";
  }
  #estimate a new a, default parameters used
  my $aHat = 0;
  my $bigI = 20;
  my $x0 = 0.8;
  my $bigFxi = 0;
  my $start = 0; # because pList is sorted, we can cal a start index for each bin
  for(my $i=0; $i<$bigI; $i++){
    my $xi = $x0+(1-$x0)*$i/$bigI;
    my $end = $start;
    while($end<@pList){
      if($pList[$end] > $xi){ #exceed bin end
        last;
      }
      $end++;
    }
    if($start <= $end-1){
     #print "Bin $i: $xi: from $start+1 ($pList[$start]) to $end ($pList[$end-1])\n";
     $bigFxi += ($end-$start)/@pList; #histogram
     #print "CDF: $bigFxi\n";
    }
    if($bigFxi>$xi){
      $aHat += ($bigFxi-$xi)/(1-$xi);
    }
    if($end >= @pList){ # boundary case
      last;
    }
    $start = $end;
  }
  $aHat /= $bigI;

  print "Estimated alternative percentage (a) out of $pListSize p-values: $aHat\n" if $isVerbose;
  if($aHat >= 1){
    return $pList[-1]; # all cases are estimated to be from the alternative
  }
  $fdrCutoff /= (1-$aHat);
  print "Adjusted FDR (BH): $fdrCutoff\n" if $isVerbose; 
  my $pos = 0;
  while($pos < @pAdjust){
    #if($pList[$pos] > $fdrCutoff/$norm*($pos+1)/@pAdjust){
    if($pAdjust[$pos] > $fdrCutoff){
      last;
    }
    $pos++;
  }
  if($pos == 0){ 
    if($isVerbose){
      print "WARNING: No p-value cutoff can satisfy FDR <= $fdrCutoff out of $pListSize p-values\n";
      print "Set a default: $fdrCutoff (NO FDR control!!) instead. \nWARNING: Recommended to set p-value cutoff instead of FDR in parameter config file\n";
    }
    return $fdrCutoff;
  }
  #print "$pos SNVs out of ".(scalar @pList)." with FDR <= $fdrCutoff (adjusted p: $pAdjust[$pos-1] original p: $pList[$pos-1])\n";
  if($pList[$pos-1] <= $orgFdrCutoff){
    return $pList[$pos-1];
  }
  #fall-back plan
  print "WARNING: The adjusted FDR method does not work (i.e. cutoff > $orgFdrCutoff). Switched to BH method\n" if $isVerbose;
  $pos = 0;
  while($pos < @pAdjust){
    if($pAdjust[$pos] > $orgFdrCutoff){
      last;
    }
    $pos++;
  }
  return $pList[$pos-1];
}

################ minor auxiliary (mainly for print outs and debugs) ####
# print out the gene snps results for certain type (powerful or ordinary snps)
sub printGetGeneSnpsResults
{
  my ($geneSnpRef, $geneSnpKey, $snpRef, $snpKey, $snvPValueCutoff) = @_;
  print "Gene level SNP ($geneSnpKey) VS transcript results\n";

  for (my $i=1; $i<=$CHRNUM; $i++){
    my ($geneRef) = getListByKeyChr($geneSnpRef, $geneSnpKey, $i); 
    my $geneChrRef = getChrGeneSnpsSorted($geneSnpRef, $geneSnpKey, $i);
    my ($snpInfoRef) = getListByKeyChr($snpRef, $snpKey, $i);
    my %snps = %$snpInfoRef; #snps information
    my @genes = @$geneChrRef;
    my %genesPrinted = (); #to store genes printed already in order not to double-print
    if(@genes){
      printChr($i); print "\n";
      #print "@genes\t";
      foreach(@genes){
        my @allGenes = split('\t', $_);
        foreach(@allGenes){
	  #get the snp information
	  my @geneMatches = split('\t', $geneRef->{$_});
	  for(@geneMatches){
	    my ($snpPos, $matchInfo, $geneName, $txStart, $id, $regStart, $regEnd) = split(';', $_);
	    if(!defined($genesPrinted{$geneName})){
	      print $geneName."\n";
	      $genesPrinted{$geneName} = 1;
	    }else{ $genesPrinted{$geneName} += 1; }
	    #print $snps{$snpPos}."\n";
	    my @allSnpInfo = split(';', $snps{$snpPos}); #separate by ;, if there are multiple snps at the same position
            foreach(@allSnpInfo){
	      my $toPrint = $matchInfo."\t".$geneName."\t".$id.":".$regStart."-".$regEnd."\n";
	      if($geneSnpKey eq 'gPowSnps'){
	        #$p."\t".$pos."\t".$alleles."\t".$snpName."\t".$refAl."\t".$altAl.";";
		my ($p, $pos, $alleles, $snpId) = getSnpInfo($_);
	        if($p <= $snvPValueCutoff){
	          print "$snpId,$p,$alleles,$pos\t".$toPrint;
	        }
	      }else{
	        print $_."\t".$toPrint;
	      }
	    }
	  }
        }
      }
      print "\n";
    }
  }

}

sub printSnpEventsResultsByType
{
  my ($snpEventsRef, $key) = @_;
  # print the snpEvents out
  #my $key = 'powSnpAlt';
  my $testEventsRef = $snpEventsRef->{$key};
  my @testEvents = @$testEventsRef;
  for(my $i=1; $i<=$CHRNUM; $i++){
    if(defined($testEvents[$i])){
      printChr($i); print "\n";
      my %geneHash = %{$testEvents[$i]};
      my @genes = keys %geneHash;
      if(@genes>0){
        foreach(@genes){
          print "$_\n";
          print $geneHash{$_}."\n";
        }
      }
    }
  }
}
1;


=head1 NAME

snpParser.pl -- All the sub-routines for SNV (sometimes denoted interchangeably as SNP) handling in the ASARP pipeline.

=head1 SYNOPSIS

	use Statistics::R; #interact with R
	require "fileParser.pl"; #sub's for input annotation files
	require "snpParser.pl"; #sub's for snps
	...

... get all configs, input files (see L<fileParser>)


	# read and parse SNVs
	my ($snpRef, $pRef) = initSnp($snpF, $POWCUTOFF);
        # suggested, get the Chi-Squared Test p-value cutoff from FDR ($FDRCUTOFF)
	$SNVPCUTOFF = fdrControl($pRef, $FDRCUTOFF, 1); #1--verbose
	# match SNVs with gene transcript annotations
	my $geneSnpRef = setGeneSnps($snpRef, $transRef);
	# match gene SNVs with AI/AT and alternative splicing (AS) events
	my ($snpEventsRef) = 
	setSnpEvents($geneSnpRef, $altRef, $splicingRef);

	# calculate NEV and filter the matched gene SNVs with AI/AT/AS events
	my ($snpsNevRef) = 
	filterSnpEventsWithNev($snpRef, $geneSnpRef, $snpEventsRef, $bedF, 
	$allEventsListRef, $NEVCUTOFFLOWER, $NEVCUTOFFUPPER); 

	# process ASE and ASARP
	my ($allAsarpsRef) = 
	processASEWithNev($snpRef, $geneSnpRef, $snpsNevRef, $SNVPCUTOFF, 
	$ASARPPCUTOFF, $ALRATIOCUTOFF);

	# format results to output
	my $outputGene = $outputFile.'.gene.prediction';
	outputRawASARP($allAsarpsRef, 'ASARPgene', $outputGene);
	my $allNarOutput = formatOutputVerNAR($allAsarpsRef);

=head1 REQUIREMENT

C<Statistics::R>: has to be installed. See http://search.cpan.org/~fangly/Statistics-R/lib/Statistics/R.pm 

=head1 DESCRIPTION

This perl file contains all the sub-routines for SNV handling and ASARP processing, as well as result formatting. They are quite procedural and one should first get the input files such as annotations and events using the sub-routines in L<fileParser>.

Basically there are 3 steps:

1. read and parse the individual SNVs

2. match the SNVs to transcripts, and then events, and then filter them based on the PSI like Normalized Expression Value (NEV) calculation

3. process the SNVs with ASE patterns and SNV pairs with other ASARP patterns: AI/AT/AS, and output the formatted results

AI/AT/AS categories are briefly illustrated below (where the red dots represent SNVs with ASARP patterns):

G<img/Types.png>

=head2 SNV List Format

The SNV list input file contains the list of all SNVs covered by RNA-Seq in 
the genes of interest, with the read counts of the reference (Ref) and alternative (Alt) alleles.

This file is space delimited with the following fields for
each SNV:
	
	chromosome
	coordinate
	alleles (reference allele>alternative allele)
	dbsnpID
	RNA-Seq counts 
	(# reads for 
	reference allele:alternative allele:wrong nucleotide)

Example file: F<../data/snp.list.fig3>

=head2 Sub-routines (major)

=over 6

=item C<initSnp>

read and parse SNVs
  
  input: ($snpF, $POWCUTOFF) --SNV file path, powerful SNV cutoff

  output ($snpRef, $pRef) 
  --reference to SNVs, categorized into powerful/non-powerful internally
  --reference to the p-value list, which will be used to control FDR

=item C<setGeneSnps>

match SNVs with gene transcript annotations

  input: ($snpRef, $transRef);
  --reference to SNVs, reference to gene transcripts
  
  output: $geneSnpRef 
  --reference to SNVs matching gene transcripts

=item C<setSnpEvents>

match gene SNVs with AI/AT and alternative splicing (AS) events

  input: ($geneSnpRef, $altRef, $splicingRef)
  --reference to gene SNVs ($geneSnpRef), 
  --reference to AI/AT events ($altRef)
  --reference to AS events ($splicingRef)

  output: ($snpEventsRef) --gene SNVs matching AI/AT/AS events

=item C<filterSnpEventsWithNev>

calculate NEV and filter the matched gene SNVs with AI/AT/AS events
  
  input: ($snpRef, $geneSnpRef, $snpEventsRef, $bedF, 
  $allEventsListRef, $NEVCUTOFFLOWER, $NEVCUTOFFUPPER)
  --reference to SNVs ($snpRef),
  --reference to gene SNVs ($geneSnpRef),
  --reference to SNVs matching AI/AT/AS events ($snpEventsRef),
  --reference to bed folder for mapped reads ($bedF),
  --reference to all parsed events ($allEvetnsListRef),
  --lower and upper cutoffs (excl.) for NEV ($NEVCUTOFFLOWER/UPPER)

  output: ($snpsNevRef) --gene SNVs matching NEV and AI/AT/AS events 

=item C<processASEWithNev>

process ASE and ASARP

  intput: ($snpRef, $geneSnpRef, $snpsNevRef, 
  $SNVPCUTOFF, $ASARPPCUTOFF, $ALRATIOCUTOFF)
  --see above for $snpRef, $geneSnpRef, $snpsNevRef
  --Chi-Squared Test p-value cutoff on individual SNVs for ASE ($SNVPCUTOFF)
  --Fisher's Exact Test p-value cutoff on target-control SNV pairs for
  ASARP ($ASARPPCUTOFF)
  --allelic ratio difference cutoff for target-control SNV pairs for
  ASARP ($ALRATIOCUTOFF)

  output: ($allAsarpsRef) --reference to all ASE and ASARP results 

=item C<outputRawASARP>

format results to output

  input: ($allAsarpsRef, $key, $outputFile)
  --reference to ASARP results ($allAsarpsRef)
  --result type to output ($key) with choices: 
  'ASEgene'--ASE results arranged by genes,
  'ASARPgene'--ASARP results arranged by genes
  'ASARPsnp' --ASARP results arranged by SNVs
  --the output file for the results ($outputFile)

  output: corresponding ASE/ASARP results written to $outputGene

=item C<formatOutputVerNAR>

format results to be like the old version for NAR
  
  input: $allAsarpsRef --see above
  
  output: ($allNarOutput) 
  --text formatted according to the old version

=head1 SEE ALSO

L<asarp>, L<fileParser>, L<MyConstants>

=head1 COPYRIGHT

This pipeline is free software; you can redistribute it and/or modify it given that the related works and authors are cited and acknowledged.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.

=head1 AUTHOR

Cyrus Tak-Ming CHAN

Xiao Lab, Department of Integrative Biology & Physiology, UCLA

=cut

