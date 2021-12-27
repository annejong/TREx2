#!/usr/bin/env perl

# Anne de Jong
# 
# remove unwanted \" in JSON
# 

use warnings;
use strict;
use Data::Dumper;
use JSON ;
use lib "/usr/molgentools/lib";
use anne_files ;

my $x = '(\[|\]|\{|\}|\:|\,)' ;
my $filename = "123.txt";
my $line = join('', anne_files::read_lines($filename)).'x' ;
my $result =  substr($line, 0, 1) ;
for (my $i=0; $i<length($line)-2; $i++) {
	if (substr($line, $i, 3) !~ m/$x\"$x/g) { $result .= substr($line, $i+1, 1) ; } 
}
print "$result\n\n";

# my $text = Dumper(decode_json($result)) ;
# $text =~ s/\$VAR1 =/     / ;
# print $text ;

undef $/;
while (<>) {
    print $result->pretty->encode($result->decode($_));
}