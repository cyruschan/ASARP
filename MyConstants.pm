#!/usr/bin/perl -w
package MyConstants;
use strict;

use base 'Exporter';
#use Readonly; #this is the better way but requires installation

our @EXPORT = qw();
our @EXPORT_OK = qw( $CHRNUM $supportedList $supportedTags );

# better way to forbid modifying, but requires installation
#Readonly::Scalar our $CHRNUM => 24; #total number of chromosomes handled
#Readonly::Scalar our $TENKB => 10000; #10KB as the bin size
#Readonly::Scalar our $supportedList => ' snps; powSnps; trans; events; geneSnps;';


our $CHRNUM = 24; #total number of chromosomes handled
#each key word should be prefixed by space ' ' and suffixed by ';'
our $supportedList = ' snps; powSnps; trans; events; gSnps; gPowSnps;'; 
our $supportedTags = ' rna; est; anno;';
1;

=head1 NAME

MyConstants -- All the constants used by the ASARP pipeline.

=head1 SYNOPSIS
  
  use MyConstants qw( $CHRNUM $supportedList $supportedTags );
  Then the constants can be used.

=head1 DESCRIPTION

This module does not really contain actual subroutines,
it only provides the constants convenient to use in the 
other related modules/source files of the ASARP pipeline.

=head2 Constants

=over 12

=item C<$CHRNUM>

The total number of chromosomes handled in the pipeline. 24 means 1-22, X and Y.

=item C<$supportedList>

The supported attributes in the internal data structures, which are with hashes and indices for all chromosomes as an array. The space and semicolon surrounding are literally added as the prefix and the suffix respectively.

=item C<$supportedTags>

The supported types (tags) of input splicing event sources for the SNVs. The space and semicolon surrounding are literally added as the prefix and the suffix respectively.


=back

=head1 SEE ALSO

L<fileParser>, L<snpParser>

=head1 COPYRIGHT

Xiao Lab, Department of Integrative Biology & Physiology, UCLA

=head1 AUTHOR

Cyrus Tak-Ming CHAN

=cut
