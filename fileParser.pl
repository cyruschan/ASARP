#!/usr/bin/perl/ -w
use strict;
#use IO::Handle;

use MyConstants qw( $CHRNUM $supportedList $supportedTags );
require "bedHandler.pl"; #sub's for handling bed files

# subroutine to get the arguments to the entry program
sub getArgs{

  my ($outputFile, $inputConfig, $inputParam) = @_;

  #my $configs='test.config';
  my $configs='default.config';
  my $params = 'default.param';

  if(defined($inputConfig)){
    $configs = $inputConfig;
    print "User ";
  }else{	print "Default ";	 }
  print "config file: $configs used\n";

  if(defined($inputParam)){
    $params = $inputParam;
    print "User ";
  }else{	print "Default ";	 }
  print "parameter file: $params used\n";

  if(defined($outputFile)){
    print "User specified output file: $outputFile\n";
  }
  
  return ($outputFile, $configs, $params);

}


# subroutine to read the *.param file and get all input parameters
sub getParameters
{
  my ($config) = @_;
  my %filemap = ();
  open(my $fh, "<", $config) or die "Cannot open parameter file: $config\n";
  print "Reading parameters in file: $config...\n";
  while(<$fh>){ #read each line of the file handle
    chomp;
    if($_ =~/^(\w+)\t\s*([0-9|\.|\+|\-|e|E]+)\s*$/){
      #get the file_var ($1) and file_path ($2)
      my $argName = lc $1;
      $filemap{$argName} = $2;
    }else{
      if(!($_ =~ /^\#/ || $_ =~ /^\s*$/)){
        print "Error: cannot parse the parameter in line: $_\n";
        exit;
      }
    }
  }
  close($fh);

  my ($POWCUTOFF, $SNVPCUTOFF, $ASARPPCUTOFF, $NEVCUTOFFLOWER, $NEVCUTOFFUPPER, $ALRATIOCUTOFF) = (undef, undef, undef, undef, undef, undef);
  for(keys %filemap){
    my ($name, $val) = ($_, $filemap{$_});
    if($name eq 'p_chi_snv'){
      $SNVPCUTOFF = $val; #p-value cutoff for the Chi-squared test
    }
    elsif($name eq 'p_fisher_pair'){
      $ASARPPCUTOFF = $val; #p-value cutoff for the Fisher'exact test
    }
    elsif($name eq 'nev_upper'){
      $NEVCUTOFFUPPER = $val;
    }
    elsif($name eq 'nev_lower'){
      $NEVCUTOFFLOWER = $val;
    }
    elsif($name eq 'ratio_diff'){
      $ALRATIOCUTOFF = $val;
    }
    elsif($name eq 'powerful_snv'){
      $POWCUTOFF = $val;
    }else{
      print "Error parsing the parameter type: $name with value $val\n";
      exit;
    }
  }

  if(!defined($POWCUTOFF)){
    die "Error: cutoff for powerful SNVs (parameter: powerful_snv) not set\n";
  }
  if(!defined($SNVPCUTOFF)){
    die "Error: p-value cutoff for the Chi-Squared Test on alleles of individual SNVs (parameter: p_chi_snv) not set\n";
  }
  if(!defined($ASARPPCUTOFF)){
    die "Error: p-value cutoff for the Fisher Exact Test on target-control SNV pairs (parameter: p_fisher_pair) not set\n";
  }
  if(!defined($NEVCUTOFFUPPER) || !defined($NEVCUTOFFLOWER) || !defined($ALRATIOCUTOFF)){
    die "Error: cutoff(s) for the Normalized Expression Value (NEV) lower/upper bound(s) and/or allele-pair ratio difference (parameters: nev_lower, nev_upper and ratio_diff) not set\n";
  }
  return ($POWCUTOFF, $SNVPCUTOFF, $ASARPPCUTOFF, $NEVCUTOFFLOWER, $NEVCUTOFFUPPER, $ALRATIOCUTOFF);

}

# subroutine to read the annotation_file.config file and get all input file names
sub getRefFileConfig
{
  my ($config) = @_;
  my %filemap = ();
  open(my $fh, "<", $config) or die "Cannot open config file: $config\n";
  print "Reading file names in config: $config...\n";
  while(<$fh>){ #read each line of the file handle
    chomp;
    if($_ =~/^(\w+)\t\s*([\w|\.|\\|\/|\*]+)\s*$/){
      #get the file_var ($1) and file_path ($2)
      my $argName = lc $1;
      $filemap{$argName} = $2;
    }
  }
  close($fh);

  #get all file names, F means file name/path
  my ($snpF, $bedF, $rnaseqF, $xiaoF, $splicingF, $estF) = ('', '', '', '', '', ''); #$genomeF
  
  # match and output all files
  for(keys %filemap){
    my $k = $_;
    #print $k."\n";
    my $curFile = $filemap{$k};
    if($k eq 'snpfile'){
      $snpF=$curFile;
    }
    elsif($k eq 'bedfolder'){
      $bedF=$curFile;
    }
    elsif($k eq 'rnaseqfile'){
      $rnaseqF=$curFile;
    }
    elsif($k eq 'xiaofile'){
      $xiaoF=$curFile;
    }
    elsif($k eq 'splicingfile'){
      $splicingF=$curFile;
    }
    elsif($k eq 'estfile'){
      $estF=$curFile;
    }else{
      die "ERROR: Unknown file_var in config file $config: $k with file/path: $curFile\n";
    }
    #elsif($k eq 'genomepath'){
    #  $genomeF=$curFile;
    #}
  }
  # checking
  if($snpF eq '' || $bedF eq '' || $xiaoF eq '' || $splicingF eq ''){
    die "ERROR: Required user-specific input file(s) not all provided:\nsnpfile=$snpF\nbedfile=$bedF\nxiaofile=$xiaoF\nsplicingfile=$splicingF\n";
  }
  if($rnaseqF eq '' || $estF eq ''){
    print "NOTE: Optional user-specific (RNA-Seq/EST) splicing events file(s) not used:\nrnaseqfile=$rnaseqF\nestfile=$estF\n";
  }
  #die "ERROR: Required annotation file(s) undefined:\nxiaofile=$xiaoF\nsplicingfile=$splicingF\n";

  # print all the file names
  #for(keys %filemap){
  #  print $_, " = \t\t", $filemap{$_}, "\n";
  #}

  #return
  return ($snpF, $bedF, $rnaseqF, $xiaoF, $splicingF, $estF);

}

## a wrapper sub for all splicing related events
sub readAllEvents
{
  my ($splicingF, $rnaseqF, $estF, $transRef, $geneNamesRef) = @_;
  my %spEventsList = ();

  my $annoEventsRef = readDifEvent($splicingF, 'anno', $geneNamesRef, $estF); #debug "");
  $spEventsList{'anno'} = $annoEventsRef;
  #printListByKey($annoEventsRef, 'events');
  
  #optional event annotations
  if($rnaseqF ne ''){ #not empty
    my $rnaEventsRef = readDifEvent($rnaseqF, 'rna', $geneNamesRef);
    #printListByKey($rnaEventsRef, 'events');
    $spEventsList{'rna'} = $rnaEventsRef;
  }
  if($estF ne ''){ #not empty
    my $estEventsRef = readEstEvent($estF, $transRef, $annoEventsRef);
    #printListByKey($estEventsRef, 'events');
    $spEventsList{'est'} = $estEventsRef;
  }

  #for my $testTag (keys %spEventsList){
  #  print "$testTag\n";
  #  my @array = @{$spEventsList{$testTag}->{'constExons'}};
  #  my %coexons = %{$array[1]};
  #  my $c = 0;
  #  for(keys %coexons){
  #    print "$_ $coexons{$_}\n";
  #    $c++;
  #    if($c>=5){ last; }
  #  }
  #  print "REF $spEventsList{$testTag}->{'constExons'} $testTag\n";
  #}
  return \%spEventsList;
}


#subroutine to read *.event files including rnaseq.event and annotation.event
#	input:  file name (including the path) to be read
#	tag:	text to indicate the type of events being parsed
#	output: the event reference, which includes hashes and indices for the events
sub readDifEvent
{
   #array of hashes
   my @chrs=();
   my @constExons = ();
   for(my $i=0; $i<=$CHRNUM; $i++){
     push @chrs, {}; #empty initialization
     push @constExons, {};
   }

   my ($fileName, $tag, $geneNamesRef, $isEstNeeded)=@_;
   if(defined($isEstNeeded) && $isEstNeeded ne ''){
     print "The constitutive exon set will be used by est events: $isEstNeeded\n";
     $isEstNeeded = 1;
   }else{
     $isEstNeeded = 0;
   }
   open(my $fh, "<", $fileName) or die "Cannot open event file: $fileName\n";
   print "Reading alternative splicing events from ", $fileName, "\n";
   my $count=0;
   my $geneName = undef;
   my $constExonList = undef;

   #to handle those genes without annotation events, but may be used in est events
   my $annotationChr = 0; #to get the gene name and confirm the chromosome faster
   my @geneNames = @$geneNamesRef;

   while(<$fh>){
     $count++;
     #if(!($count%10000)){	print $count, "\t";     }
     chomp;
     if($_ =~/^\>/){
       my @geneInfo = split('\t', $_);
       $geneName = substr($geneInfo[0],1); #gene name no >
       $geneName = uc $geneName;

       $constExonList = undef; #for checking
       if(@geneInfo > 1){
         $constExonList = $geneInfo[1];
         #find and set the const exons for the gen even it (possibly) has no events
	 if($isEstNeeded){
	   for(my $x = 0; $x <= $CHRNUM; $x++){
	     if(defined($geneNames[$x]) && defined($geneNames[$x]->{$geneName})){
               #if($x > $currentChr+1){
	       #  print "Warning: gene $geneName is in a different chromosome ($x)\n";
	       #  #print "from the neighboring events ($currentChr); skipped\n";
	       #  last;
	       #}
	       #if($x == $currentChr + 1){
	       #  print "Warning special: $geneName in $x\n";
	       #}
	       $constExons[$x]{$geneName} = $constExonList;
	       $annotationChr = $x;
	       #print "$x\t$geneName\n";
	       last;
	     }
	     #if($x==$CHRNUM){ print  "Warning: cannot find $geneName in any chromosomes\n"; }
	   }
	 }
       }
     }elsif($_ =~/^EVENT\t/){
       #get the line
       my ($dummy, $chrRaw, $geneNameCheck, $strand, $evtRegion, $lRegion, $rRegion)
       = split(/\t/, $_);
       $geneNameCheck = uc $geneNameCheck;
       if($geneNameCheck ne $geneName){
         die "Inconsistent gene name in line $count: $geneNameCheck and $geneName\n";
       }
       my $chrID = getChrID($chrRaw);
       #check numeric, skip all non-cannonical chromosomes
       if(!($chrID=~/^\d+$/)){ 
         #print "$chrID\n"; 
	 next; 
       }
       # store the constitutive exon set list to the gene
       if(!defined($constExonList)){
         die "Error: No const exon info for $geneName EVENT!\n";
       }
       if(!defined($constExons[$chrID]{$geneName})){
         $constExons[$chrID]{$geneName} = $constExonList;
	 #if($isEstNeeded && $annotationChr != $chrID){
	   #print "Warning: Gene $geneName event in ".formatChr($chrID).", annotation in ".formatChr($annotationChr)."\n";
	 #}
       }
       #could make it 0-based by $chrID-=1; here, but not convenient for biologists
       my ($start, $end) = split(/\-/, $evtRegion);

       # try to pack the information
       if($strand == 1){
         $strand = '+';
       }else{
         $strand = '-';
       }
       # determine splicing type
       my ($lStart, $lEnd) = split('-', $lRegion);
       my ($rStart, $rEnd) = split('-', $rRegion);
       my $spliceType = deriveSpliceType($start, $end, $lStart, $lEnd, $rStart, $rEnd);

       $lRegion = $lStart.":".$lEnd;  #join(':', split('-', $lRegion));
       $rRegion = $rStart.":".$rEnd;  #join(':', split('-', $rRegion));
       $chrs[$chrID]{$start}.=$end.";".$geneName.";".$lRegion.";".$rRegion.";".$strand.";".$spliceType."\t"; # derived SP type added to all events
     }
   }
   close($fh);

   #set the indices (such that event starts are sorted) for chrs
   my @chr_idx = ();
   for(my $i=1; $i<=$CHRNUM; $i++){
     #printChr($i); 
     $chr_idx[$i] = [sort {$a<=>$b} keys %{$chrs[$i]}];
   }
  my %eventList = ( 
      'events' => \@chrs,
      'events_idx' => \@chr_idx,
      'type' => $tag,
      'constExons' => \@constExons,
  );  
  return \%eventList;
}

# to compile all the splicing (not alt 5'/3') events according to gene arrangement
sub compileGeneSplicingEvents
{
  my ($geneRef, @eList) = @_;
  my ($rnaRef, $annoRef, $estRef) = (undef, undef, undef);
  my ($rnaIdx, $annoIdx, $estIdx) = (undef, undef, undef);
  
  my @events = ();
  my @events_idx = ();
  for(my $i=0; $i<=$CHRNUM; $i++){
    push @events, {};
    push @events_idx, ();
  }
  
  # get all references and indices
  foreach(@eList){
    my $list = $_;
    my $tag = getListTag($list);
    if($tag eq 'rna' || $tag eq 'est' || $tag eq 'anno'){ #valid lists
      if($tag eq 'rna'){  $rnaRef = $list; }
      elsif($tag eq 'est'){  $estRef = $list;  }
      else{ $annoRef = $list;  } #$tag eq 'anno'
    }else{
      print "Warning: not supported tag $tag for this list (ignored)\n";
    }
  }

  # store information based on gene's index
  for(my $i=1; $i<=$CHRNUM; ++$i){
     if(defined($rnaRef)){   $events[$i] = addEventsByChr($events[$i], $geneRef, $rnaRef, $i, 'rna'); }
     if(defined($estRef)){   $events[$i] = addEventsByChr($events[$i], $geneRef, $estRef, $i, 'est'); }
     if(defined($annoRef)){   $events[$i] = addEventsByChr($events[$i], $geneRef, $annoRef, $i, 'anno'); }
  }

  # build indices
  my @genes = @{$geneRef};
  for(my $i=1; $i<=$CHRNUM; ++$i){
    my $x = 0;
    my @chrGenes = @{$genes[$i]};
    my %chrEvents = %{$events[$i]};
    my @chrEventGenes = ();
    foreach(@chrGenes){
      my ($gene, $s, $e) = split(';', $_);
      if(defined($chrEvents{$gene})){
        $chrEventGenes[$x] = $gene;
	$x++;
      }
    }
    $events_idx[$i] = \@chrEventGenes;
  }
  
  my %geneEvents = (
      'events' => \@events,
      'events_idx' => \@events_idx,
      'type' => 'all combined splicing events arranged by genes',
  );
  return \%geneEvents;
}

# add particular event to the genes on the specified chromosome
sub addEventsByChr
{
  my ($eventChrRef, $geneRef, $rnaRef, $chr, $tag) = @_;
  
  my %eventHash = %$eventChrRef;

  my @geneIdxAll = @{$geneRef};
  my @geneIdx = @{$geneIdxAll[$chr]};

  my ($rnaHashRef, $rnaIdxRef) = getListByKeyChr($rnaRef, 'events', $chr);
  my %rnas = %$rnaHashRef;
  my @rnaIdx = @$rnaIdxRef;
 
  for(my $ei = 0; $ei<@rnaIdx; $ei++){
    my $eStart = $rnaIdx[$ei];
    my @eSet = split('\t', $rnas{$eStart});
    foreach(@eSet){
      my ($eEnd, $gene, $lRegion, $rRegion, $strand, $additional) = split(';', $_);
      #all events have the uniform format now, and can be distinguished by $tag
      $eventHash{$gene} .= join(';', $eStart.":".$eEnd, $lRegion, $rRegion, $strand, $additional, $tag)."\t";
    }
  }
  return \%eventHash;

}

# subroutine to read in the merged transcriptome annotation, i.e. xiaoF
# ID, chr, strand, txStart, txEnd, cdsstart, cdsend, exoncount, exonstarts, 
# exonends, genename, cdsstartstat,cdsendstat
# 	input: 	file name (including the path) to be read
#	output:	the hashes and indices for the events


sub readTranscriptFile
{

  #array of hashes
  my @trans=();
  for(my $i=0; $i<=$CHRNUM; $i++){
    push @trans, {}; #empty initialization
  }

  #keep track of discarded chromosomes
  my %discardedChrs = ();

  my ($fileName)=@_;
  open(my $fh, "<", $fileName) or die "Cannot open merged transcrimptome file: $fileName\n";
  print "Reading the transcriptome file ", $fileName, "\n"; 
  
  my $count = 0;
  while(<$fh>){
    $count++;
    if(!($count%10000)){
      #print $count, "\t";
      #STDOUT->autoflush(1);#need to flush if STDOUT not attached to terminal
    }
    chomp;
    my ($ID, $chrRaw, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd, 
    $exonCount, $exonStarts, $exonEnds, $geneName, $cdsStartStat, $cdsEndStat)
    = split(/\t/, $_);# ID, chr, strand, txStart, txEnd, 
    #cdsstart, cdsend, exoncount, exonstarts, 
    # exonends, genename, cdsstartstat,cdsendstat
    my $chrID = getChrID($chrRaw);
    if(!($chrID=~/^\d+$/)){
      if(!defined($discardedChrs{$chrID})){
        $discardedChrs{$chrID}=1;
      }else{ $discardedChrs{$chrID}+=1;  }
      #print "$chrID\n"; 
      #print "Not valid chr ID: $chrID\n"; 
      next;
    }
    $ID = uc $ID; # to upper
    #change 0-based to 1-based for all starts
    $txStart+=1;
    $cdsStart+=1;
    my @allExonStarts = split(/,/, $exonStarts);
    for(my $i=0; $i< @allExonStarts; $i++){
      $allExonStarts[$i]+=1;
    }
    $exonStarts = join(",", @allExonStarts);
    #but no need to add 1 for ends
    $exonEnds = substr($exonEnds, 0, -1); #because there is an additional "," in the original file

    my $isCoding = 1;
    #if($geneName=~/noncoding/){  $isCoding = 0; }
    if($cdsStart > $cdsEnd){ $isCoding = 0; }
    #if($isCoding==0){   print "Non-coding: $ID, $geneName\n"; }
    # currently skip calculating the 5UTR 3UTR and exon length here
    
    $trans[$chrID]{$txStart}.=$txEnd.";".$cdsStart.";".$cdsEnd.";".$exonStarts.";".$exonEnds.";".$ID.";".uc($geneName).";".$isCoding.";".$strand."\t";
  }
  close($fh);

  #set the indices (such that event starts are sorted) for chrs
  my @trans_idx = ();

  for(my $i=1; $i<=$CHRNUM; $i++){
    #printChr($i); 
    $trans_idx[$i] = [sort {$a<=>$b} keys %{$trans[$i]}];
  }
  #debug print
  #my $chrToPrint = 4;
  #printChr($chrToPrint); print "\n";
  #foreach (@{$trans_idx[$chrToPrint]}){ print $_-1,"\n"; }
  
  #printDiscardedChrs(\%discardedChrs);  
  
  my %transList = (
    'trans' => \@trans,
    'trans_idx' => \@trans_idx,
    'type' => 'transcripts (not events)',
  );
  return \%transList;
}

# at the same time, we can also get the multiple 5' and 3' transcript ends in this sub
sub getGeneAltTransEnds{
  my ($transRef) = @_;
  # this is the init for 5' and 3' transcript ends
  my @allTransEnds = ();
  for(my $i=0; $i<=$CHRNUM; $i++){
    push @allTransEnds, {};
  }

  for(my $i=1; $i<=$CHRNUM; $i++){

    my ($chrTransRef, $chrTransIdxRef) = getListByKeyChr($transRef, 'trans', $i);
    my @trans_idx = @$chrTransIdxRef;
    my %trans = %$chrTransRef; #chromosome specific
    my %transEnds = ();

    foreach(@trans_idx){
      my $txStart = $_;
      my @transList = split('\t',$trans{$txStart});
      foreach(@transList){
        my ($txEnd, $cdsS, $cdsE, $exSs, $exEs, $id, $gene, $isCoding, $strand) = split(';', $_);
	
	######################## for transcript ends
	my @exStSet = split(',', $exSs);
	my @exEnSet = split(',', $exEs);
	my %altHash = (
	    '5+'=>{},
	    '3+'=>{},
	    '5-'=>{},
	    '3-'=>{},
	);
	if(defined($transEnds{$gene})){
	   %altHash = %{$transEnds{$gene}};
	}
	  
	if($strand eq '+'){
	  # for 5+ end, key is the first exon end, and values are all possible txStarts
	  $altHash{'5+'} = addAltKeyValuePair($altHash{'5+'}, $exEnSet[0], $txStart);
	  # for 3+ end, key is the last exon start, and values are all possible txEnds
          $altHash{'3+'} = addAltKeyValuePair($altHash{'3+'}, $exStSet[-1], $txEnd);
	}else{
	  # for 5- end, key is the last exon start in + (i.e. first exon end on -)
	  #, and values are all possible txEnds in + (i.e. all possible txStarts on -)
	  $altHash{'5-'} = addAltKeyValuePair($altHash{'5-'}, $exStSet[-1], $txEnd);
	  # for 3- end, key is the first exon end on +, and values are all possible txStarts
          $altHash{'3-'} = addAltKeyValuePair($altHash{'3-'}, $exEnSet[0], $txStart);	    
	}
	$transEnds{$gene} = \%altHash;

      }
    } #end of foreach trans_idx
    ######################### for transcript ends
    $allTransEnds[$i] = \%transEnds;    
  }

  return \@allTransEnds;
}

# to add an key-value pair to the alternative ends hash
sub addAltKeyValuePair{
  my ($ref, $key, $value) = @_;
  my %hash = %$ref;
  #only add unique value
  my $stub = $value.";";
  if(!defined($hash{$key}) || !($hash{$key} =~/$stub/)){
    $hash{$key} .= $value.";";
  }
  return \%hash;
}

# auxiliary: print out the 5' and 3' alternative starts
sub printAltEnds{
  my ($allTransEndsRef) = @_;
  print "Alternative ends type:\n";
  for(my $i=1; $i<=$CHRNUM; $i++){
    printChr($i); print "\n";
    my %trans = %{getAltEndsListByChr($allTransEndsRef, $i)};
    my @genes = keys %trans;
    if(@genes>0){ 
      printChr($i); print "\n"; 
      foreach(@genes){
        print $_."\n";
        printAltEvent($trans{$_}, '5+');
        printAltEvent($trans{$_}, '5-');
        printAltEvent($trans{$_}, '3+');
        printAltEvent($trans{$_}, '3-');
      }
    }
  }

}

sub getAltEndsListByChr{
  my ($allTransEndsRef, $chr) = @_;
  my @allTransEnds = @$allTransEndsRef;
  return $allTransEnds[$chr];
}

# just to output the alt events under one gene
sub printAltEvent
{
  my ($ref, $altType) = @_;
  my $hsRef = getAltEventByType($ref, $altType);
  my %hash = %$hsRef;
  for(keys %hash){
    my @allAlts = split(';', $hash{$_});
    if(@allAlts>0){
      print " has $altType ".scalar @allAlts." alts: $_ - $hash{$_}\n"; 
    }
  }
}

# get the junction (neighboring) read info for an alt event
sub getAltEventByType
{
  my ($ref, $altType) = @_;
  my %hs = %$ref;
  return $hs{$altType};
  
}

# sub getGeneIndex to get each gene's starting location (min txStart) in each chromosome
sub getGeneIndex{
  my ($transRef) = @_;

  # this is the init for gene indices
  my @geneIndices = ();
  my @geneNames = ();


  for(my $i=0; $i<=$CHRNUM; $i++){
    push @geneIndices, (); #array of arrays of strings
    push @geneNames, ();
  }
  for(my $i=1; $i<=$CHRNUM; $i++){
    my %genes = ();
    my @genesArray = (); # the locations
    my $g = 0; #gene counter


    my ($chrTransRef, $chrTransIdxRef) = getListByKeyChr($transRef, 'trans', $i);
    my @trans_idx = @$chrTransIdxRef;
    my %trans = %$chrTransRef; #chromosome specific

    foreach(@trans_idx){
      my $txStart = $_;
      my @transList = split('\t', $trans{$txStart});
      foreach(@transList){
        my ($txEnd, $cdsS, $cdsE, $exSs, $exEs, $id, $gene, $isCoding, $strand) = split(';', $_);
	######################### for gene indices
	if(!defined($genes{$gene})){
          $genes{$gene} = $txStart; #initial minimal
	  $genesArray[$g] = $gene.';'.$txStart.';'.$txEnd;
	  $g++;
        }elsif($genes{$gene}> $txStart){
          die "This is not possible for $gene at $genes{$gene} as $txStart is sorted\n";
        }

      }
    }
    ######################### for gene indices
    $geneIndices[$i] = \@genesArray;
    $geneNames[$i] = \%genes;
  }

  return (\@geneIndices, \@geneNames);
}

# sub routine to read the est event file (est.event)
sub readEstEventFile{

  my ($fileName)=@_;
  open(my $fh, "<", $fileName) or die "Cannot open est events file: $fileName\n";
  print "Reading est event info from ", $fileName, "\n"; 
  
  #array of hashes
  my @ests=();
  for(my $i=0; $i<=$CHRNUM; $i++){
    push @ests, {}; #empty initialization
  }
  
  my $count = 0;
  while(<$fh>){
    $count++;
    if(!($count%10000)){
      #print $count, "\t";
      #STDOUT->autoflush(1);#need to flush if STDOUT not attached to terminal
    }
    chomp;
    my ($type, $info, $start, $end) = split('\t', $_);
    my ($chrRaw, $pos, $strand) = split(':', $info);
    my $chr = getChrID($chrRaw);

    #temporary indexed by $start
    $ests[$chr]{$start} .= "$end;$type;$pos;$strand\t";
  }
  close($fh);
  
  #create index for genome locations
  my @ests_idx = ();

  for(my $i=1; $i<=$CHRNUM; $i++){
    #printChr($i); 
    $ests_idx[$i] = [sort {$a<=>$b} keys %{$ests[$i]}];
  }

  my %estList = ( 
      'events' => \@ests,
      'events_idx' => \@ests_idx,
      'type' => 'est raw info (not events)',
  );  
  return \%estList;
}
sub readEstEvent{
 
  my $checkStrandConsist = 0; #by defatult: not checking strand consistency between EST and Transcripts
  my $isNonFlankingEventKept = 1; #whether to keep events without flanking regions: no-0; yes-1
  my ($fileName, $transRef, $annoEventsRef, $checkStrandInputFlag)=@_;
  if(defined($checkStrandInputFlag)){
    $checkStrandConsist = $checkStrandInputFlag;
  }
  my $estListRef = readEstEventFile($fileName);
  my @chrs = ();
  for(my $i=0; $i<=$CHRNUM; $i++){
    #push @splice, {}; #empty initialization
    push @chrs, {};
  }

  print "Mapping ests to transcripts to create est events\n";
  #print "Matching splicgraph event to genes\n"; #printout copied from Gang's code
  for(my $i=1; $i<=$CHRNUM; $i++){

    my ($estRef, $estIdxRef) = getListByKeyChr($estListRef, 'events', $i);
    my @chrEsts_idx = @{$estIdxRef};
    my ($chrTransRef, $chrTransIdxRef) = getListByKeyChr($transRef, 'trans', $i);
    my @trans_idx = @$chrTransIdxRef;
    my %trans = %$chrTransRef; #chromosome specific
    my ($ei, $ti) = (0, 0); #est index and transcript index

    #print "SP Events: ".scalar @chrEsts_idx."; Transcripts: ".scalar @trans_idx."\n"; 
    while($ei < @chrEsts_idx && $ti < @trans_idx){
       my $eStart = $chrEsts_idx[$ei];
       my $tPos = $trans_idx[$ti];
       
       # an event must be inside certain transcript (i.e. both the eStart>=$trans_idx[$ti] and eEnd<=txEnd in $trans{$trans_idx[$ti]}
       if($eStart < $trans_idx[$ti]){ $ei++; next; }
       else{ #event left VS transcript left done, i.e., eStart>=$trans_idx[$ti], need further checking
         my $newTi = $ti;
	 while($newTi<@trans_idx && $eStart >= $trans_idx[$newTi]){
	   #print "$ei, $ti ($newTi): $eStart, $trans_idx[$newTi]\n";
	   my $tPos = $trans_idx[$newTi];
	   my @tSet = split('\t', $trans{$tPos});
	   my $maxTxEnd = -1; # to store the largest transcript end of @tSet
	   foreach(@tSet){
	     my ($txEnd, $cdsStart, $cdsEnd, $exonStarts, $exonEnds, $ID, $gene, $isCoding, $txStrand) = split(';', $_);
	     if($txEnd > $maxTxEnd){	$maxTxEnd = $txEnd;	}
	     if($eStart<=$txEnd){ # event left VS this transcript right in @tSet done, need further checking
	       my @exss = split(',', $exonStarts);
	       my @exes = split(',', $exonEnds);
	       my $exNo = @exss;
	       #check all $eEnd's (event right's) in the set with the same $eStart
	       my @evntSet = split('\t', $estRef->{$eStart}); 
	       foreach(@evntSet){
                 my $isInJunction = 0; # to indicate whether the est event is within certain junction
		 my ($eEnd,$type,$pos,$strand) = split(';', $_);
	         # start investigating event right ($eEnd)     
	         #perform a binary search to locate the start against exon ends: boundary cases (i.e. matches) put to the left
                 #corrected: junction should be something in between exons (intronic)
		 my ($loc, $unMatchFlag) = binarySearch(\@exes, $eStart, 0, $exNo-1, 'left');
	       
	         if($loc>0 && $loc<$exNo){ #not the first nor after the last exon in the exon set
	           # no need to do binary search on exes, just check whether the end is < exon start
		   if($eEnd < $exss[$loc]){
		     #event right VS transcript exon end: eStart-$eEnd is some junction
		     
		     if(1){#
		       # a new event for this gene
		       # calculate the junction region
		       $isInJunction = 1; # there is a junction
		       my ($lEnd, $rStart) = ($exes[$loc-1], $exss[$loc]); #left flanking region end and right flanking region start
                       #left flanking region start is the closest exon start
		       my $lStart = $exss[$loc-1]; #initially the previous exon's start
		       # so there may be other exon sets starting with the same start but with different ends
		       for(my $lPos = $loc-1; $lPos>0 && $exss[$lPos-1]>=$exss[$loc-1]; --$lPos){
		         if($lStart < $exss[$lPos-1]){ #find the maximal
		           $lStart = $exss[$lPos-1];
		         }
		       }
		       my $rEnd = $exes[$loc]; #initial
		       for(my $rPos = $loc+1; $rPos<$exNo && $exes[$rPos]<=$exes[$loc]; ++$rPos){
		         if($rEnd > $exes[$rPos]){ #find the minimal
		           $rEnd = $exes[$rPos]; 
		         }
		       }

		       my $spliceType = deriveSpliceType($eStart, $eEnd, $lStart, $lEnd, $rStart, $rEnd);
		       if(1){ #$spliceType ne 'UN') #keeping 'UN' derived type as well
		         # 'D' means 'D'erived from the flanking regions.
		         #$splice[$i]{$gene} .= join(';', ('D', $spliceType, $eStart, $eEnd, $lStart, $lEnd, $rStart, $rEnd))."\t"; 
		         #print $gene." ".$splice[$i]{$gene}."\n";

                         my $lRegion = $lStart.":".$lEnd;
			 my $rRegion = $rStart.":".$rEnd;

			 if(!$checkStrandConsist || $txStrand eq $strand){
                           $chrs[$i]{$eStart}.=$eEnd.";".$gene.";".$lRegion.";".$rRegion.";".$txStrand.";".$spliceType."\t";
                         }
		       }

		     }
		   }
		 }

	         if(!$isInJunction && $isNonFlankingEventKept){
	           #there may be multiple events
		   #add this into the gene: 'F' means the original est.event 'F'ile information is used.
	           #print $gene." ".$splice[$i]{$gene}."\n";
	           if(!$checkStrandConsist || $txStrand eq $strand){
		     $chrs[$i]{$eStart}.=$eEnd.";".$gene.";-1:-1;-1:-1;".$txStrand.";F:".$type."\t";
		   }
	         }
	       } #end of foreach(@evntSet)
	     }
	   } #end of foreach (@tSet)
	   if($eStart > $maxTxEnd && $newTi == $ti){ #no longer to check this again
	     $ti += 1;
	   }
           $newTi++;
	 } #end while $eStrat and $newTi
	 $ei += 1; #$ei has gone through all the $newTi
       
       } #end else
    
    } # end while $ei, $ti

    #all the splicing events should be ready by now for chromosome $i
  }
  
  #set the indices (such that event starts are sorted) for chrs
  my @chr_idx = ();
  for(my $i=1; $i<=$CHRNUM; $i++){
    #printChr($i); 
    $chr_idx[$i] = [sort {$a<=>$b} keys %{$chrs[$i]}];
  }
  my %eventList = ( 
      'events' => \@chrs,
      'events_idx' => \@chr_idx,
      'type' => 'est',
      'constExons' => $annoEventsRef->{'constExons'}, 
  );  
  return \%eventList;
  #return \@splice;
}

# sub routine to determine the splicing type
sub deriveSpliceType{
  my ($eStart, $eEnd, $lStart, $lEnd, $rStart, $rEnd) = @_;
  my $type = '';
  # SE: skipped exon: the event region is in between the left and right flanking regions
  if($lEnd<$eStart-1 && $rStart>$eEnd+1){
    $type = 'SE';
  }elsif($lEnd==$eStart-1 && $rStart==$eEnd+1){
    #RI: the event region just connects the left and right flanking regions
    $type = 'RI';
  }elsif($lEnd==$eStart-1 || $rStart==$eEnd+1){
    #ASS: A5SS (left connected) or A3SS (right connected)
    $type = 'ASS';
  }else{ $type = 'UN'; } #unknown
  return $type;
}

##################################################################################
### the following sub's are about processing 5' alt init and 3' alt term

sub getAlt5Or3Events
{
  my ($transListRef, $geneRef) = @_;
  my @alt5Int = ();
  my @alt5Int_idx = ();
  for(my $i=0; $i<=$CHRNUM; $i++){
    push @alt5Int, {};
    push @alt5Int_idx, ();
  }

  my @allGenes = @$geneRef;
  for(my $i=0; $i<=$CHRNUM; $i++){
    my ($transR, $transIdx) = getListByKeyChr($transListRef, 'trans', $i);
    my %trans = %$transR;
    my @transIdx = @$transIdx;
    my @genes = @{$allGenes[$i]}; #just for this chromosome

    foreach(@genes){
      my ($gene, $txStart, $txEnd) = split(';', $_);

    }
  }
}

##################################################################################
# auxiliary subroutines
# to process the structured data of the annotation files and snps 
# In general, this structured data uses arrays of hashes to store information, 
# and arrays of arrays (suffix _idx) to store the sorted locations in the genome

# the general print procedure for the arrays of hash+idx structures used 
# in the parser files, e.g. fileParser.pl and snpParser.pl
sub printListByKey{
   my ($ref, $key) = @_;
   my ($hsRef, $idxRef) = getListArraysByKey($ref, $key);

   my @hsArry = @{$hsRef}; #hash array
   my @idxArry = @{$idxRef}; #indexing array
   
   print "Key: $key\n";
   for(my $i=1; $i<=$CHRNUM; $i++){
     #see if it is empty
     my @chr_idx = @{$idxArry[$i]};
     my %chr_hash = %{$hsArry[$i]};
     if(@chr_idx){
       printChr($i); print "\n";
       my $x=0;
       foreach(@chr_idx){
         $x++;
         print $_."\t".$chr_hash{$_}."\n";
       }
       print "\n";
     }
   }
}

# get list tag (type of resources)
# input		the list reference of events of a certain type
# output	the tagged type for this list
sub getListTag{
  my ($listRef) = @_;
  return $listRef->{'type'};
}

# get chromosome specific information (both hash and array) from the structured data
sub getListByKeyChr{
  my ($ref, $key, $chr) = @_;
  #get the arrays of all chromosomes first
  my ($hsRef, $idxRef) = getListArraysByKey($ref, $key);
  my @hsArry = @{$hsRef}; #hash array for the whole genome (i.e. all chromosomes)
  my @idxArry = @{$idxRef}; #indexing array

  #get the specific chromosome info
  return ($hsArry[$chr], $idxArry[$chr]); #reference to hash, reference to idx (location array)
}


# get the ref of arrays of hashes and arrays (idx) for the whole genome
sub getListArraysByKey{
  my ($ref, $key) = @_;
    my ($hsArryRef, $idxArryRef) = (undef, undef);
    if(checkSupportedList($key)){
      $hsArryRef = $ref->{$key}; #hash array
      $idxArryRef = $ref->{$key."_idx"}; #indexing array
    } 
    return ($hsArryRef, $idxArryRef);
}

sub getConstitutiveExonsByChr{
  my ($ref, $chr) = @_;
  my $checkTag = $ref->{'type'};
  my @chrsConsExons = ();
  if(checkSupportedType($checkTag)){
     @chrsConsExons = @{$ref->{'constExons'}};
  }else{
    print "cannot get const for $checkTag: $chr\n";
  }
  my %coexons = %{$chrsConsExons[$chr]};
  #my $testRef = \%coexons;
  #print "test REF: $testRef for chr $chr\n";
  #print "chrsConsExons for $chr: $chrsConsExons[$chr]\n";
  #return $testRef;
  return $chrsConsExons[$chr];
}

# checking procedure if any call above provides a valid key
sub checkSupportedList{
  my $key = $_[0];
  # supported list: i.e. what have been implemented in the parsers
  #my $supportedList = " snps; powSnps; trans; events; geneSnps;";
  my $keyStub = ' '.$key.';'; #add the header and footer to prevent excessive matching
  # $supportedList is an external constant in MyConstants.pm
  if($supportedList=~/$keyStub/){
    return 1;
  }else{
    die "Error using structured data: Unknown/unsupported list key: $key\n";
  }
}


# checking procedure if any call above provides a valid key
sub checkSupportedType{
  my $tag = $_[0];
  my $tagStub = ' '.$tag.';'; #add the header and footer to prevent excessive matching
  # $supportedTags is an external constant in MyConstants.pm
  if($supportedTags =~/$tagStub/){
    return 1;
  }else{
    die "Error using structured data: Unknown/unsupported type (used to get events/constitutive exon sets): $tag\n";
  }
}

# assumption: the values in the list are unique
# binary search for insert (or location), including left, right insert (location) cases
sub binarySearch{
  my ($listRef, $x, $imin, $imax, $type) = @_;
  $type = uc $type;
  my @arr = @$listRef;
  while($imin <= $imax){
    my $midPoint = int(($imin+$imax)/2);
    if($arr[$midPoint] < $x){
      $imin = $midPoint+1;
    }elsif($arr[$midPoint] > $x){
      $imax = $midPoint-1;
    }else{ # a match
      if($type eq 'right'){
        return ($midPoint+1, 0); #put it to the right
      }else{
        return ($midPoint, 0); #put it to the left: default
      }
    }
  }
  #this means $imin = $imax+1, and neither $arr[$imin] or $arr[$imax] equals $x
  if($imax<0){ return (0, 1); } # $arr[$imax] when imax == -1 has a different meaning in Perl!!
  if($arr[$imax] < $x){ return ($imin, 1);  }
  return ($imax, 1);

}


################ minor auxiliary functions mainly for debug info ####################################
# parse the chromosome for its numeric ID (not numeric if the chr is special)
sub getChrID
{
  my $chrID = uc $_[0]; 
  if($chrID=~'^CHR'){
    $chrID=substr($chrID, 3); 
  }
  if($chrID eq 'X'){ $chrID = 23; }
  elsif($chrID eq 'Y'){ $chrID = 24; }
  
  return $chrID;
}

# print out the chromosome based on its numeric value 1-$CHRNUM
sub printChr
{
   my $i=$_[0];
   print "Chr ";
   if($i>0 && $i<23){  print "$i";  }
   elsif($i==23){  print "X"; }
   elsif($i==24){  print "Y"; }
   else{ die "Unknown Chr $i\n"; }
}

# print out the chromosome based on its numeric value 1-$CHRNUM
sub formatChr
{
   my $i=$_[0];
   my $str = "chr";
   my $chrX ='';
   if($i>0 && $i<23){  $chrX = $i;  }
   elsif($i==23){  $chrX = "X"; }
   elsif($i==24){  $chrX = "Y"; }
   else{ die "Unknown Chr $i\n"; }

   return $str.$chrX;
}

# print those discarded chromosomes
sub printDiscardedChrs
{ 
  my %hs=%{$_[0]};
  print "Alternative chromosomes discarded with counts:\n";
  foreach(sort keys %hs){
    print $_." (".$hs{$_}.");";
  }
  print "\n";
}

1;

=head1 NAME

filePaser.pl -- All the sub-routines for getting and parsing input files (NOT involving SNV handling) in the ASARP pipeline.

=head1 SYNOPSIS

	require "fileParser.pl";

	# input arguments: $outputFile--output, 
	# $configs--input configuration file, $params--parameter file
	my ($outputFile, $configs, $params) = getArgs(@ARGV); 
	my ($snpF, $bedF, $rnaseqF, $xiaoF, $splicingF, $estF) = getRefFileConfig($configs);
	my ($POWCUTOFF, $SNVPCUTOFF, $ASARPPCUTOFF, $NEVCUTOFFLOWER, $NEVCUTOFFUPPER, $ALRATIOCUTOFF) = getParameters($params);

	# read the transcript annotation file
	my $transRef = readTranscriptFile($xiaoF);
	
	#get alternative initiation/termination (AI/AT) events from transcripts
	my $altRef = getGeneAltTransEnds($transRef); 

	# get indices of gene transcript starts and gene names (prepared also for SNVs)
	my ($genesRef, $geneNamesRef) = getGeneIndex($transRef); 

	# read all annotations, optionally rna-seq and est, events and compile them
	my $allEventsListRef = readAllEvents($splicingF, $rnaseqF, $estF, $transRef, $geneNamesRef);
	my $splicingRef = compileGeneSplicingEvents($genesRef, values %$allEventsListRef);

=head1 DESCRIPTION

This perl file contains all the sub-routines that handle the input arguments, read configuration files, and parse all the annotations and events that are input to be matched with SNVs to discover ASE/ASARP. SNVs are handled separately in L<snpParser>.

=head2 Sub-routines (major)

These are the major (interface) sub-routines that will be used in correlation with SNVs in the whole pipeline. Read them one by one as they are quite procedural.

=over 6

=item C<getRefFileConfig>

get all input annotation/event file/folder paths contained in the configuration file.

 input: $configs --configuration file, check out default.config for the formats.
 
 output: ($snpF, $bedF, $rnaseqF, $xiaoF, $splicingF, $estF) 
 --SNV list file path ($snpF), 
 --the (sam.)bed folder path ($bedF), 
 --rna-seq.event file path ($rnaseqF, optional: '' returned if not provided in $configs), 
 --transcript annotation file path ($xiaoF), 
 --annotation.event path ($splicingF), 
 --est.event file path ($estF, optional: '' returned if not provided in $configs)

The default config file as an example can be found in F<../default.config>.

Lines for RNA-Seq.event and EST.event file paths may be skipped as they are optional.

B<File formats>:

=over 8

=item Annotation file format:
Transcript and gene annotation specified in C<$xiaoF>

The example file by default, F<../data/hg19.merged.to.ensg.all.tx.03.18.2011.txt>, 
was created by merging ensembl Refseq, UCSC knowngene, Gencode
gene, and Vegagene. 

I<Format> (tab delimited):
ID, chr, strand, txStart,
txEnd, cdsstart, cdsend, exoncount, exonstarts, exonends, genename,
cdsstartstat,cdsendstat

I<IMPORTANT>: all coordinates are hg19, 0-based start and 1-based end 
coordinates (UCSC convention) **in this file only???**.

=item Event (splicing event) file formats:

Annotation events specified in C<$splicingF> (example: F<../data/annotation.event>)

The file contains splicing events as annotated in the above file
(C<$xiaoF>).  The format is the same as that for rnaseq_event.
(1-based start and end)

RNA-Seq events specified in C<$rnaseqF> (example: F<../data/rnaseq.event>)

The file contains splicing events as determined by our RNA-seq data.
It lists the events for each gene. The format of the events is, 
EVENT, chromosome, genename, strand, event_region, flanking_region_1, 
flanking_region_2, where *_region are in the format of 
starting_coordinate-ending_coordinate. (1-based start and end)

EST events specified in C<$estF> (example: F<../data/est.event>) 

The file contains splicing events as determined from hg19 EST and cDNA data.
The format is tab-delimited as: event_type, event_name, starting_coordinate,
ending_coordinate. (1-based start and end)

=back

=item C<getParameters>

get all the numeric parameters including p-value cutoffs, NEV lower/upper thresholds, and the allelic ratio difference threshold.

 input: $params --configuration file for the parameters

 output: ($POWCUTOFF, $SNVPCUTOFF, $ASARPPCUTOFF, $NEVCUTOFFLOWER, $NEVCUTOFFUPPER, $ALRATIOCUTOFF)
 --read count cutoff for powerful SNVs ($POWCUTOFF),
 --Chi-Squared Test p-value cutoff for individual SNVs ($SNVPCUTOFF),
 --Fisher's Exact Test p-value cutoff for target-control SNV pairs in ASARP ($ASARPPCUTOFF),
 --NEV lower and upper cutoffs (excl.) ($NEVCUTOFFLOWER, $NEVCUTOFFUPPER),
 --allelic ratio difference cutoff for target-control SNV pairs in ASARP ($ALRATIOCUTOFF)

The default parameter config file as an example can be found in F<../default.param>.
All parameters are required and have to be set.

=item C<readTranscriptFile>

read the transcript annotation file

 input: $xiaoF --file path of the transcript (also gene) annotation file. 

 output: $transRef --reference to an array of the parsed transcripts. 

The results can be printed out using utility sub using key 'trans': C<printListByKey($transRef, 'trans');>


=item C<getGeneIndex>

intermediate sub to get indices of gene transcript starts and gene names

 input: $transRef (see above)

 output: $genesRef --reference to the index for every gene's minimal transcript start

 	 $geneNamesRef --reference to the gene names
	
=item C<getGeneAltTransEnds>

get alternative initiation/termination (AI/AT) events from transcripts

 input: $transRef (see above)

 output: $altRef--reference to the AI/AT events

The AI/AT results can be printed out using utility sub: C<printAltEnds($altRef);>

=item C<readAllEvents>

read all annotations, optionally rna-seq and est, splicing events

 input: ($splicingF, $rnaseqF, $estF, $transRef, $geneNamesRef);
 --event files (see above): anno ($splicingF), rna ($rnaseqF), est ($estF)
 --see above for $transRef, $geneNamesRef
 
 output: $allEventsListRef --the reference to a hash of 
 all 'anno', optionally 'rna', 'est' events parsed
 quoted are keys to access them in the hash
	

=item C<compileGeneSplicingEvents>

compile events from different files and arrange them according to genes

 input: ($genesRef, values %$allEventsListRef)
 --see above for $genesRef
 --the 2nd argument is the event type(s) provided, e.g. 'anno', 'rna', 'est'
 output: $splicingRef --reference to all the splicing events 

=back

=head2 Utility subs

Utility sub-routines are implemented to get or display parsing annotations/events for intermediate usage or checking. 

Parsing results are usually arranged in an internal data structure whose reference is returned. 

=over 6

=item The 1st layer is a hash with 3 keys: the data 'x', its companian index 'x_idx', and 'type' indicating the tag or telling briefly what the data is about. 

=item The 2nd layer is an array with number of elements equal to the chromosome number (see constant $CHRNUM in L<MyConstants>). 

=item The 3rd layer (i.e. each element representing one chromosome) is the reference to a hash containing specific information of that structure. E.g. in the case of transcript annotation, the keys are the starting locations of the transcripts, and the values are the transcript information where multiple transcripts are separated by "\t".

=back

=over 6

=item C<printListByKey>

the general print procedure for the arrays of hash+idx structures used 

  input: ($ref, $key) --the reference to a structure and the corresponding key for the structure.
  output: print out the content of the structure chromosome by chromosome.

e.g. C<printListByKey($transRef, 'trans');> to print out all transcripts; 

or C<printListByKey($snpRef, 'powSnps');> to print out all powerful SNVs (SNPs)

=item C<getListTag>

get list tag (type of resources)

  input: the reference to the structure of events/annotations of a certain type
  output: the tagged type (or simply tag) for this list

=item C<getListByKeyChr>

get the references to hash+idx of a particular chromsome number from an internal structure

  input: ($ref, $key, $chr) 
  --reference to the structure, corresponding key of the structure, the chromosome number specified (1 - $CHRNUM)
  
  output: ($hsArry[$chr], $idxArry[$chr]) 
  --the references to hash+idx of a particular chromsome number

This utility sub is frequently used in the intermediate processing of all kinds of events, annotations and SNVs.

=item C<binarySearch>

serve as a general utility to get the location of an element in an ordered list (a good example is the indices _idx)

binary search for insert (or location), including left, right insert (location) cases, similar to bisectin Python. Assume the elements in the list are sorted in an ascending order and there is no duplicate.

  input: ($listRef, $x, $imin, $imax, $type) 
  --the reference to a list, the element to be searched ($x), 
  --starting index ($imin), e.g. 0, ending index ($imax), e.g. size of the list - 1,
  --$type: 'left' or 'right' mimicing left or right bisect. 'left': if $x matches an element, the return location will be left to the element; 'right': return location right to the element

  output: ($loc, $flag) 
  --location of $x in the list (NOTE: range is 0 to size of list, 0-base)
  --match flag: 1 means match, 0 means no match

=item Other utility subs

Who actually cares? If you do, kindly go to have a look at the source: F<../fileParser.pl>. There will be some comments around. Good luck! : )

=back

=head1 SEE ALSO

L<snpParser>, L<MyConstants>

=head1 COPYRIGHT

This pipeline is free software; you can redistribute it and/or modify it given that the related works and authors are cited and acknowledged.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.

=head1 AUTHOR

Cyrus Tak-Ming CHAN

Xiao Lab, Department of Integrative Biology & Physiology, UCLA

=cut
