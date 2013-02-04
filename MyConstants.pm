#!/usr/bin/perl -w
package MyConstants;
use strict;

use base 'Exporter';
#use Readonly; #this is the better way but requires installation

our @EXPORT = qw();
our @EXPORT_OK = qw( $CHRNUM $TENKB $supportedList $supportedTags );

# better way to forbid modifying, but requires installation
#Readonly::Scalar our $CHRNUM => 24; #total number of chromosomes handled
#Readonly::Scalar our $TENKB => 10000; #10KB as the bin size
#Readonly::Scalar our $supportedList => ' snps; powSnps; trans; events; geneSnps;';


our $CHRNUM = 24; #total number of chromosomes handled
our $TENKB = 10000; #10KB as the bin size
#each key word should be prefixed by space ' ' and suffixed by ';'
our $supportedList = ' snps; powSnps; trans; events; gSnps; gPowSnps;'; 
our $supportedTags = ' rna; est; anno;';
1;
