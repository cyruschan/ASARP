#!/usr/bin/perl

use strict;
use warnings;

# set autoflush for error and output
select(STDERR);
$| = 1;
select(STDOUT);
$| = 1;

if(@ARGV < 4){
  print <<EOT;
USAGE: perl $0 folder_list prefix suffix output_folder

ARGUMENTS:
folder_list 	a file containing a line-separated list of folders 
          	representing mutilple biological replicates
		e.g. data/rep1
		     data/rep2
		     data/rep3
		     ...
prefix		the common chromosome file prefix in these folders
		e.g. proc_chr for proc_chr1.sam, proc_chr2.sam, etc.
suffix		the common file suffix in these folders
		e.g. rmdup.sam for chr1.rmdup.sam, chr2.rmdup.sam, etc.
output_folder	all the merged prefix*.suffix files will be put to that 
		specified folder (assume it exists)

NOTE: 	files matching folder/prefix*.suffix will be
	enumerated for merging, where "folder" is one element
	in the folder_list. Note the added / and *. in between.
	
	The files are simply merged and not sorted.

EXAMPLE: If one wants to merge all "proc_chr*.rmdup.sam" files under
	 data/rep1, data/rep2, data/rep3 (list stored in "rep.lst") 
	 and output the merged files to folder "data.merged", 
	 mergeSam.pl should be run as follows:
	 perl mergeSam.pl rep.lst proc_chr rmdup.sam dadta.merged
EOT
  exit;
}

my ($inFolders, $filePrfx, $fileExtn, $outFolder) = @ARGV;

open(FP, "<", $inFolders) or die "ERROR: cannot open the input folder list: $inFolders\n";
my @inFlds = <FP>;
chomp @inFlds;
close(FP);

my %allChrs = ();
for(@inFlds){
  my $all = $_."/".$filePrfx."*.".$fileExtn;
  my @chrSams = glob("$all");
  if(!@chrSams){
    die "ERROR: cannot get paths/files with $all\n";
  }
  for(@chrSams){
    #if($_ =~ /$_\/$filePrfx(\w+)\.$fileExtn/){ #get the middle
    if($_ =~ /$filePrfx(\w+)/){ #get the middle
      #print $1."\n"; exit;
      if(!defined($allChrs{$1})){
        $allChrs{$1} = $_;
      }else{
        $allChrs{$1} .= "\t".$_;
      }
    }else{
      die "ERROR: cannot get the file name after $filePrfx in path $_\n";
    }
  }
}
for(keys %allChrs){
  my @paths = split('\t', $allChrs{$_});
  my $folderNo = @paths;
  print "Merging chr$_ in $folderNo replicates: lines..."; #$allChrs{$_}\n\n";
  my $outputName = $outFolder."/chr".$_.".sam";
  open(FP, ">", $outputName) or die "ERROR: cannot open $outputName to output merged SAM files\n";
  for(my $i = 0; $i < @paths; $i++){
    open(IP, "<", $paths[$i]) or die "ERROR: cannot open $paths[$i] to merge\n";
    my @lines = <IP>;
    print " ".(scalar @lines)."...";
    print FP @lines;
    close(IP);
  }
  close(FP);
  print " DONE\n";
}

__END__

=head1 NAME

mergeSam.pl -- Merging corresponding SAM files (e.g. all chr1 SAM files) from different independent replicates. Each replicate is indicated by its folder. The SAM files should have been subject to duplicate removal (L<rmDup>) 

=head1 SYNOPSIS

This is part of the full pre-processing:

=over 6

1. rmDup (removing PCR duplicates for SAM files (including Dr. JH Lee's SAM format))

2. B<mergeSam> (merging SAM files if there are independent duplicates)

3. procReads (processing SAM files to get SNV read counts and generate bedgraph files) 

=back

USAGE: 
 perl mergeSam.pl folder_list prefix suffix output_folder

ARGUMENTS:

 folder_list 	a file containing a line-separated list of folders 
          	representing mutilple biological replicates
		e.g. data/rep1
		     data/rep2
		     data/rep3
		     ...
 prefix		the common chromosome file prefix in these folders
		e.g. proc_chr for proc_chr1.sam, proc_chr2.sam, etc.
 suffix		the common file suffix in these folders
		e.g. rmdup.sam for chr1.rmdup.sam, chr2.rmdup.sam, etc.
 output_folder	all the merged prefix*.suffix files will be put to that 
		specified folder (assume it exists)

NOTE:

files matching folder/prefix*.suffix will be
enumerated for merging, where "folder" is one element
in the folder_list. Note the added / and *. in between.

The files are simply merged and not sorted.

EXAMPLE: 

If one wants to merge all "proc_chr*.rmdup.sam" files under
data/rep1, data/rep2, data/rep3 (list stored in "rep.lst") 
and output the merged files to folder "data.merged", 
mergeSam.pl should be run as follows:

 perl mergeSam.pl rep.lst proc_chr rmdup.sam dadta.merged

=head1 SEE ALSO

L<rmDup>, L<procReads>, L<asarp>

=head1 COPYRIGHT

This pipeline is free software; you can redistribute it and/or modify it given that the related works and authors are cited and acknowledged.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.

=head1 AUTHOR

Cyrus Tak-Ming CHAN

Xiao Lab, Department of Integrative Biology & Physiology, UCLA

=cut
