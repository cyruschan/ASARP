#!/usr/bin/perl
use warnings;
use strict;
require "anlysUtilities.pl"; #sub's for compSamples

# set autoflush for error and output
select(STDERR);
$| = 1;
select(STDOUT);
$| = 1;

if(@ARGV < 9){

  print <<EOT;

USAGE: perl $0 intron_output1 config1 param1 result1 intron_output2 config2 param2 result2 output [asarp_path]

This pipeline combines several perl scripts to analyze the intronic ASE SNVs
from two samples (datasets), especially in the context where the two samples 
represent different celluar compartments, e.g. cytocol and nucleuse.

In the input, 1, 2 indicate the corresponding in/output for samples (datasets) 
1 and 2, respectively. 

intron_output	the intronic ASE snv intermediate output prefix
		severval files with suffixes: .ase .ase_distri.list, etc.
		will be generated as the intermediate output for one sample
		NOTE: make sure the file names are **DIFFERENT** between 
		intron_output1 and intron_output2
config		the .config file which is required for ASARP pipeline entry
param		the .param file, optional in the ASARP pipeline, is required
		in this pipeline
result		the result (main output summary file) of ASARP
		the corresponding .ase.prediction and .gene.prediction will
		be needed (assumed to be present with the main file)
output		the finaly output file of this analysis pipeline

OPTIONAL:
asarp_path	the -I path for the ASARP pipeline folder if it is run outside
		the original folder. Because this pipeline executes several
		scripts in ASARP: aseSnvs.pl, snp_distri.pl, etc., the ASARP
		path is needed for them.

EOT

exit;
}

# focus on the intron SNVs only and compare them with the 2 sample results
print "Compare intronic ASE SNVs\n";

my ($intron1, $config1, $param1, $result1, $intron2, $config2, $param2, $result2, $output, $iPath) = @ARGV;
if(!defined($iPath)){
  $iPath = ".";
}
if(substr($iPath, -1, 1) ne "/" ){
  $iPath .= "/"; #add the path slash
}

my $specificType = 'INTRON';
print "\n[[I]]. Obtain all intronic ASE SNVs from sample 1\n";
my ($iAseGenesRef1, $iAseSnvsRef1, $iAsarpGenesRef1, $iAsarpSnvsRef1) = specificAsePipeline($intron1, $config1, $param1, $result1, $specificType, $iPath);

print "\n[[II]]. Obtain all intronic ASE SNVs from sample 2\n";
my ($iAseGenesRef2, $iAseSnvsRef2, $iAsarpGenesRef2, $iAsarpSnvsRef2) = specificAsePipeline($intron2, $config2, $param2, $result2, $specificType, $iPath);

# now you can intersect what ever you want
print "\n[[III]]. Intersect ASE Genes with $specificType ASE SNVs\n";
my ($iAse1, $iCom, $iAse2) = intersectHashes($iAseGenesRef1, $iAseGenesRef2);
print "ASE genes containing intronic ASE SNVs:\n";
my %iAseHs1 = %$iAse1; 
my %iAseHs2 = %$iAse2;
my %iComHs = %$iCom;
my $noAG1 = keys %iAseHs1;
my $noAG2 = keys %iAseHs2;
my $noComAG = keys %iComHs;
print "$noAG1\t$noComAG\t$noAG2\n";  

my ($iAseS1, $iComS, $iAseS2) = intersectHashes($iAseSnvsRef1, $iAseSnvsRef2);
print "$specificType ASE SNVs in ASE genes:\n";
my %iAseSnvHs1 = %$iAseS1; 
my %iAseSnvHs2 = %$iAseS2;
my %iComSnvHs = %$iComS;
my $noAS1 = keys %iAseSnvHs1;
my $noAS2 = keys %iAseSnvHs2;
my $noComAS = keys %iComSnvHs;
print "$noAS1\t$noComAS\t$noAS2\n";  

# what to intersect then?
print "\n[[IV]]. Compare Retained Intron (RI) genes with intronic SNVs\n";
# may use part of samplesAnlys.pl to get only common RI results
my $asarpType = 'ASARP'; print STDERR "fake asarp type $asarpType for debugging only\n";
my $riGeneRef1 =$iAsarpGenesRef1->{$asarpType};
my $riGeneRef2 =$iAsarpGenesRef2->{$asarpType};

my ($iRiRef1, $iRiComRef, $iRiRef2) = intersectHashes($riGeneRef1, $riGeneRef2);
print "Intronic ASE SNVs common in genes with RI (Retained Introns):\n";
my %iRiComGenes = %$iRiComRef; 
my $noRiComGenes = keys %iRiComGenes;
print "$noRiComGenes\n";
for(keys %iRiComGenes){
  print "$_\n$iRiComGenes{$_}\n";
}

###################################################################################

sub specificAsePipeline
{
  my ($outputAse, $configs, $params, $result, $specificType, $iPath) = @_;
  #strand-specific setting is more related to file configs (data dependent)
  my ($snpF, $bedF, $rnaseqF, $xiaoF, $splicingF, $estF, $STRANDFLAG) = getRefFileConfig($configs); # input annotation/event files
  my ($POWCUTOFF, $SNVPCUTOFF, $FDRCUTOFF, $ASARPPCUTOFF, $NEVCUTOFFLOWER, $NEVCUTOFFUPPER, $ALRATIOCUTOFF) = getParameters($params); # parameters

  # run aseSnvs to get ASE SNVs only
  my $aseSnvs = "$outputAse.ase";
  print "\n[1]. Get ASE SNVs and output them to $aseSnvs\n\n";
  my $aseSnvCmd = "perl -I $iPath $iPath"."aseSnvs.pl $outputAse $configs";
  if(!defined($params)){
    $aseSnvCmd .=" $params"; #optional parameter
  }
  print "$aseSnvCmd\n";
  if(system($aseSnvCmd)){
    die "FAILED to finish $aseSnvCmd\n";
  }

  # run snp_distri.pl to get the SNV distributions
  my $aseSnvDistri = "$outputAse.ase_distri.list";
  print "\n[2]. Get SNV distributions of ASE SNVs from $aseSnvs and output the detailed list to $aseSnvDistri";
  if($STRANDFLAG){
    print ".plus and .minus";
  }
  print "\n(summary to $outputAse.distri_summary.txt)\n\n";
  # the input will be "$outputAse.ase"
  my $snvDistCmd = "perl -I $iPath $iPath"."snp_distri.pl $outputAse.distri_summary.txt $aseSnvs $xiaoF $STRANDFLAG $POWCUTOFF $aseSnvDistri";
  print "$snvDistCmd\n";
  if(system($snvDistCmd)){
    die "FAILED to finish $snvDistCmd\n";
  }

  # if strand-specific
  #$ordOut = $snvOrdOut.".plus";
  #$ordOutRc = $snvOrdOut.".minus";
  # get only the intronic positions from the list

  # two lists are needed (.plus and .minus) if strand-specific flag is set
  print "\n[3]. Get only $specificType ASE SNV information from $aseSnvDistri";
  if($STRANDFLAG){
    print ".plus and .minus";
  }
  print"\n\n";

  my %introns = ();
  if(!$STRANDFLAG){ #nss
    my @snvs = readSpecificNoStrandInfo($aseSnvDistri, $specificType);
    for(@snvs){ $introns{$_} = 1; }
    my $ttlNo = @snvs;
    print "There are $ttlNo $specificType ASE SNVs\n";
  }else{ #strand-specific
    my @snvs = readSpecificNoStrandInfo("$aseSnvDistri.plus", $specificType); 
    for(@snvs){ $introns{"$_;+"} = 1; }
    my @snvsRc = readSpecificNoStrandInfo("$aseSnvDistri.minus", $specificType); 
    for(@snvsRc){ $introns{"$_;-"} = 1; }
    my $ttlNoPlus = @snvs;
    my $ttlNoMinus = @snvsRc;
    my $ttlNo = $ttlNoPlus + $ttlNoMinus;
    print "There are $ttlNo $specificType ASE SNVs ($ttlNoPlus + and $ttlNoMinus -)\n";
  }
  #return \%introns;
  
  print "\n[4]. Get ASE and ASARP results from $result and keep only $specificType results\n\n";
  my ($aseGeneRef, $aseSnvRef) = getAseAll("$result.ase.prediction");
  my ($iAseGenesRef, $iAseSnvsRef) = getSpecificAse($aseGeneRef, $aseSnvRef, \%introns);

  my ($asarpGeneRef) = getAsarpAll("$result.gene.prediction");
  my ($iAsarpGenesRef, $iAsarpSnvsRef) = getSpecificAsarp($asarpGeneRef, \%introns);

  return ($iAseGenesRef, $iAseSnvsRef, $iAsarpGenesRef, $iAsarpSnvsRef);
}

sub getSpecificAsarp
{
  my ($geneRef, $intronRef) = @_;
  my %asarpGenes = %$geneRef;
  my %intronSnvs = %$intronRef;

  my %intronAsarpGenes = ();
  my %intronAsarpSnvs = ();
  for(keys %asarpGenes){ # for every type in the results
    my %intronAsGenes = ();
    my %intronAsSnvs = ();
    my %asGenes = %{$asarpGenes{$_}};
    for(keys %asGenes){
      my ($chr, $gene) = split(';', $_);
      my $hasIntron = 0;
      my @snvs = split(/\t/, $asGenes{$_});
      for(@snvs){
        #split to get detailed information
	my ($dummyInfo, $snpInfo, $strandInfo) = split(';', $_);
	my($pos, $id, $al, $reads) = split(' ', $snpInfo);
	my $snvKey = "$chr;$pos";
	if(defined($strandInfo)){
	  $snvKey .= ";$strandInfo";
	}
	if(defined($intronSnvs{$snvKey})){
	  $hasIntron = 1;
	  $intronAsSnvs{$snvKey} .= "$_\t"; #which gene(s) the intron SNV belongs to
	}
      }
      if($hasIntron){
        $intronAsGenes{$_} = $asGenes{$_};
      }
    }
    # now only one ASARP type is finished (AS), other types will be finished similarly
    $intronAsarpSnvs{$_} = \%intronAsSnvs;
    $intronAsarpGenes{$_} = \%intronAsGenes;
  }
  return (\%intronAsarpGenes, \%intronAsarpSnvs); 
}

sub getSpecificAse
{
  my ($geneRef, $aseRef, $intronRef) = @_;
  my %aseGenes = %$geneRef;
  my %aseSnvs = %$aseRef;
  my %intronSnvs = %$intronRef;

  my %intronAseGenes = (); #ASE genes that contain the corresponding introns
  my %intronAseSnvs = (); #all ASE SNVs contained by certain genes
  for(keys %aseGenes){
    # check every gene for its Snvs
    my $hasIntron = 0;
    my @allAseSnvs = split(/\t/, $aseGenes{$_});
    my $geneKey = $_; #to get debug info
    for(@allAseSnvs){
      if(defined($intronSnvs{$_})){
        $hasIntron = 1; #print "HIT: $_ is in ASE GENE: $geneKey\n";
	$intronAseSnvs{$_} = $aseSnvs{$_}; # have to save all intronic SNVs;
	#one intron SNV belonging to (even) multiple ASE genes always has the same information
      }
    }
    if($hasIntron){
      $intronAseGenes{$_} = $aseGenes{$_}; # to include this gene
    }
  }

  return (\%intronAseGenes, \%intronAseSnvs);
}

sub readSpecificNoStrandInfo
{
  my @snvs = ();
  my ($aseSnvDistri, $specificType) = @_;
  open(FP, $aseSnvDistri) or die "ERROR: cannot open ASE SNV distribution file from $aseSnvDistri\n";
  my @lines = <FP>;
  close(FP);
  chomp @lines;

  for(@lines){
    if(index($_, "$specificType;") != -1){ 
      print STDERR "now use $specificType, maybe a fake setting for debug!!\n";
      my ($chr, $pos, $type) = split(/\t/, $_);
      push @snvs, "$chr;$pos";
    }
  }

  return @snvs;
}