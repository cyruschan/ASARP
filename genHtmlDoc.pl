#!/usr/bin/perl -w
use strict;
use Pod::HtmlEasy;

if(@ARGV < 3){
  die "USAGE: perl $0 input.pod output.html title\n";
}
my ($inputPod, $outputHtml, $title) = @ARGV;
#safety checks
if(!($inputPod =~ /\.[pod|pl|pm|plx]/)){
  print "ERROR: input $inputPod must have a perl-file extension: pod/pl/pm/plx\n";
  exit;
}

if(!($outputHtml =~ /\.html$/)){
  print "ERROR: output must be html: $outputHtml\n";
  exit;
}
if(!defined($title)){
  print "WARNING: no title is given, use $inputPod instead\n";
  $title = $inputPod;
}

# overwrite on_L to support internal html links and image link inclusion
my $podhtml = Pod::HtmlEasy->new(
  on_G => sub {
    my ( $this , $txt ) = @_ ;
    return "<img src='$txt' border=0>" ; #customized link for images
  },

  on_F => sub {
    my ( $this , $txt ) = @_ ;
    return "<i><a href='$txt' target='_blank'>$txt</a></i>"; #customized link for files
  },
  
  on_L => sub {
    my ( $this , $L , $text, $page , $section, $type ) = @_ ;
    if   ( $type eq 'pod' ) {
      if(defined($section)){
        $section = "#$section" if $section ne '' ;
      }
      return "<i><a href='$text.html'>$text</a></i>"; #relative link customized for this internal doc
      #return "<i><a href='http://search.cpan.org/perldoc?$page$section'>$text</a></i>" ;	
    }												 
    elsif( $type eq 'man' ) { return "<i>$text</i>" ;}
    elsif( $type eq 'url' ) { return "<a href='$page' target='_blank'>$text</a>" ;}  
  }
);
my %options = (
css => 'default.css',
output => $outputHtml,
title => $title,
);

$podhtml->pod2html($inputPod, %options);

linkSpaceWorkaround($outputHtml, "$outputHtml");

sub linkSpaceWorkaround{

  #Pod::HtmlEasy will always generates the abbreviated forms of a html/ftp link
  #in order to work around this, one has to add a space between http:// and the actual link,
  #e.g. http://www.abc.com will be http:// www.abc.com to stop this
  #to maintain the link reachable, this work-around simply again
  #replaces all 'html:// ' and 'ftp:// ' with 'html://' and 'ftp://' using sed
  # if outhtml is not input, the final output will be the input itself
  my ($inhtml, $outhtml) = @_;
  if(!defined($outhtml)){
    $outhtml = $inhtml; # in = out: in-place
  }
  #  sed 's/ftp\:\/\/ /ftp\:\/\//g' <temp.txt >temp.new.txt
  # replace ftp:// and then html:// 
  system("sed 's/ftp\\:\\/\\/ /ftp\\:\\/\\//g' <$inhtml >$inhtml.tmp");
  system("sed 's/http\\:\\/\\/ /http\\:\\/\\//g' <$inhtml.tmp >$outhtml");
  system("rm $inhtml.tmp");

}
