#!/usr/bin/perl
use warnings;
use strict;

# generate the index and the main frame

my $menuContent = "";

# title and index page
genTitleFrame();
genIndexFrameset();

# all sections (and in order)
my @sections = qw (INTRO PREPROCESS ASARP ANALYSIS CORE-SUB);
my %layout = ();
for(@sections){
  $layout{$_} = {};
}

# all pre-processing scripts
my @pres = qw( rmDup mergeSam procReads );
for(@pres){
  system("perl genHtmlDoc.pl $_.pl doc/$_.html $_");
}
$layout{"PREPROCESS"} = \@pres; # in order
# main program, all application scripts
my @apps = qw( asarp aseSnvs snp_distri );
for(@apps){
  system("perl genHtmlDoc.pl $_.pl doc/$_.html $_");
}
$layout{"ASARP"} = \@apps; # in order

## core sub-routines (code)
# all source code
my @source = qw( fileParser snpParser );
for(@source){
  system("perl genHtmlDoc.pl $_.pl doc/$_.html $_");
}
# all constants (wrapped in a module)
my @const = qw( MyConstants );
for(@const){
  system("perl genHtmlDoc.pl $_.pm doc/$_.html $_");
}
my @cores = (@const, @source);
$layout{"CORE-SUB"} = \@cores; # in order

# menu page
genMenuPage(\@pres, \@apps, \@source, \@const);

sub genMenuPage{
  my ($preRef, $appRef, $srcRef, $conRef) = @_;

  open(FP, ">", "doc/menu.html") or die "Can't write menu.html";
  print FP <<"EOMENU";
<html>
  <head>
  <base target="content">
    <link rel="stylesheet" href="default.css" type="text/css">
  </head>
  
  <body style="background-color:#eee">
    <!-- Content Goes Here -->
    <div class='pod'>

EOMENU
  
  my %hs = (
   "PRES" => $preRef, 
   "APPS" => $appRef,
   "CODE" => $srcRef,
   "CONST" => $conRef,
  );
  my %hsNames = (
    "PRES" => "Pre-processing",
    "APPS" => "Applications",
    "CODE" => "Core Subs",
    "CONST" => "Constants",
  );

  for(keys %hs){
    print FP "\n\t<h1>$hsNames{$_}</h1>\n\t<ul>\n";

    my @names = @{$hs{$_}};
    for(@names){
      print FP "\t\t<li><a href='$_.html'>$_</a></li>\n";
    }
    print FP "\t</ul>\n";
  }
  print FP "    </div>\n  </body>\n</html>";
  close (FP);
}


sub genTitleFrame{
  open(FP, ">", "doc/title.html") or die "Can't write title.html";
  print FP <<"EOTITLE";
<html>
  <body style="background-color:#006699">
  <h2><font color="white">Documentation of ASARP in Perl</font></h2>
  </body>
</html>

EOTITLE
  close(FP);

}

sub genIndexFrameset{
  open(FP, ">", "doc/index.html") or die "Can't write index.html";
  print FP <<"EOINDEX";
<html>
  <head>
    <link rel="stylesheet" href="default.css" type="text/css">
  </head>
  
  <frameset border="0" frameborder="0" framespacing="0" rows="60px,*">
    <frame src="title.html" noresize scrolling="no">
    <frameset border="0" frameborder="0" framespacing="0" cols="160px,*">
      <frame name="menu" src="menu.html" scrolling="auto" noresize>
      <frame name="content" src="asarp.html" scrolling="auto" noresize>
    </frameset>
  </frameset>
</html>

EOINDEX
  close(FP);
}
