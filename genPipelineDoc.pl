#!/usr/bin/perl -w
use strict;

# generate the index and the main frame

# all source code
my @source = qw( asarp fileParser snpParser );
for(@source){
  system("perl genHtmlDoc.pl $_.pl doc/$_.html $_");
}

# perl module (only the constants)
my @pms = qw( MyConstants );
for(@pms){
  system("perl genHtmlDoc.pl $_.pm doc/$_.html $_");
}

