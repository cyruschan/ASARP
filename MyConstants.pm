#!/usr/bin/perl -w
package MyConstants;
use strict;

use base 'Exporter';
#use Readonly; #this is the better way but requires installation

our @EXPORT = qw();
our @EXPORT_OK = qw( $CHRNUM $supportedList $supportedTags $asarpTags );

# better way to forbid modifying, but requires installation. e.g.
#Readonly::Scalar our $CHRNUM => 24; #total number of chromosomes handled
#Readonly::Scalar our $supportedList => ' snps; powSnps; trans; events; geneSnps;';

our $CHRNUM = 25; #total number of chromosomes handled
#each key word should be prefixed by space ' ' and suffixed by ';'
our $supportedList = ' snps; powSnps; trans; events; gSnps; gPowSnps;'; 
our $supportedTags = ' rna; est; anno;';

# supported ASARP (sub-) types for investigations
our $asarpTags = 'AS AI AT ASS SE RI';
1;

=head1 NAME

MyConstants -- All the constants used by the ASARP pipeline.

=head1 SYNOPSIS
  
  use MyConstants qw( $CHRNUM $supportedList $supportedTags );

Then the constants can be used.

=head1 DESCRIPTION

This module does not really contain any actual subroutines,
it only provides the constants convenient to use in the 
other related modules/source files of the ASARP pipeline.

=head2 Constants

=over 6

=item C<$CHRNUM>

The total number of chromosomes supported and handled in the pipeline. Currently it is set to 24 for chromsomes 1-22, X and Y. Need to modify implementations of L<fileParser> and L<snpParser> to include extra chromosomes.

=item C<$supportedList>

The supported attribute keywords for querying the internal data structures. Each of these structure is formatted inside its own wrapper as an array of hashes and indices, where each array index indicates one chromosome, i.e. 1, 2, ..., 23 (X) or 24 (Y). The space and semicolon surrounding each keyword are literally added as the disambiguating prefix and suffix. See the utility subs section in L<fileParser> for more details.

=item C<$supportedTags>

The supported types (tags) of input splicing event sources for the SNVs, namely 'anno', 'rna' and 'est'. The space and semicolon surrounding each type (tag) are literally added as the disambiguating prefix and suffix. One can extend this by modifying the implementations. See the input file formats in L<fileParser> for more details.

=back

=head1 SEE ALSO

L<fileParser>, L<snpParser>

=head1 COPYRIGHT

This pipeline is free software; you can redistribute it and/or modify it given that the related works and authors are cited and acknowledged.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.

=head1 AUTHOR

Cyrus Tak-Ming CHAN

Xiao Lab, Department of Integrative Biology & Physiology, UCLA

=cut
