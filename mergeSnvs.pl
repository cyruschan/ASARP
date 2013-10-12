#!/usr/bin/perl
use warnings;
use strict;

require "readUtilities.pl";

if(@ARGV < 4){
  print <<EOT;

USAGE perl $0 prefix suffix mono=0/1 output [isStrandSp]
merge SNVs of different chromosomes into one SNV list
the equivalent can be achieved using cat and grep commands
for strand-specific cases '+\$' and '\\-\$' can get all SNVs from the
+ and - strands respectively
prefix	the folder path and prefix string for all chr* SNV files
	e.g. "/home/N.A+/" for all chr* SNVs in that folder
suffix	the suffix after chr[1..22/X/Y/M]
mono=	options (no space): 
	"mono=0": allow mono-allelic SNVs; 
	"mono=x": discard SNVs with any allele < x reads; 
	useful when genotype is unknown and dbSNPs are used
	where mono-allelic SNVs would be called
output	the output file name for the merged SNVs
	note that output.plus and output.minus are by-product
	outputs for + and - strands for strand-specific cases
	make sure the output and chr* SNV file names do not conflict
OPTIONAL
isStrandSp	set 1 if the SNV files are strand-specific

EOT
  exit;
}

my ($prefix, $suffix, $mono, $output, $isStranded) = @ARGV;
$mono = lc $mono;
my $isMono = 0;
if($mono =~ /mono=(\d+)/){
  $isMono = $1; # matched
}

if($isMono >= 0){
  if($isMono == 0){
    print "Note: Allow mono-allelic SNVs (suppose heterozyous SNVs are known from genomic sequencing)\n";
  }else{
    print "Note: Discard SNVs with any allele (reference or alternative) < $isMono reads (in the case SNVs are called from RNA-Seq)\n";
  }
}else{
  die "ERROR: mono= should be non-negative ($mono input)\n";
}
mergeSnvs($prefix, $suffix, $isMono, $output, $isStranded);

__END__

=head1 NAME

mergeSnvs.pl -- merge SNVs of different chromosomes into one SNV list
the equivalent can be achieved using cat and grep commands
for strand-specific cases '+$' and '\-$' can get all SNVs from the
+ and - strands respectively

After the full pre-processing, one can run ASARP for each chromosome's
SNVs, or merge all chromosomes' SNVs into a single SNV list using this
script (equivalent to using Unix cat and/or grep).

=head1 SYNOPSIS

USAGE:

 perl $0 prefix suffix mono=0/1 output [isStrandSp]

 prefix		the folder path and prefix string for all chr* SNV files
 		e.g. "/home/N.A+/" for all chr* SNVs in that folder
 suffix		the suffix after chr[1..22/X/Y/M]
 mono=		2 options: 
		"mono=1": allow mono-allelic SNVs; 
		"mono=0": discard mono-allelic SNVs; 
		useful when genotype is unknown and dbSNPs are used
 output		the output file name for the merged SNVs
 		note that output.plus and output.minus are by-product
 		outputs for + and - strands for strand-specific cases
 		make sure the output and chr* SNV file names do not conflict

 OPTIONAL
 isStrandSp	set 1 if the SNV files are strand-specific


=head1 SEE ALSO

L<rmDup>, L<mergeSam>, L<procReads>, L<procReadsJ>, L<asarp>

=head1 COPYRIGHT

This pipeline is free software; you can redistribute it and/or modify it given that the related works and authors are cited and acknowledged.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.

=head1 AUTHOR

Cyrus Tak-Ming CHAN

Xiao Lab, Department of Integrative Biology & Physiology, UCLA

=cut
