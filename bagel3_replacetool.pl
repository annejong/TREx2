#!/usr/bin/perl

# when you move the BAGEL3 results to another location, use this tool to change the url's of results files

use strict ;
use Data::Dumper ;
use lib "/usr/molgentools/lib";
use anne_files ;


my $sessiondir = '.';
my $old_url = "http://genome2d.molgenrug.nl/genome2d_results/T-REx/94.212.212.43808lbr9smek3iekcn4iasif4s0916" ;
my $new_url = "http://ngs.molgenrug.nl/zhu/2015-12/Analysis_VC" ;
my $regex	= "html\$";


my $usage = "/usr/molgentools/trex/bagel3_replacetool.pl
				-s Sessiondir [default=current folder]
				-old old url
				-new new url
				";

my $filenames = anne_files::get_files_from_subdirs( $sessiondir, $regex) ;
foreach my $file (@$filenames) {
	print "$file\n";
	my @lines = anne_files::read_lines($file);
	my @result ;
	foreach my $line (@lines) {
		$line =~ s/$old_url/$new_url/g ;
		push @result, $line ;
	}
	anne_files::write_lines($file, @result) ;	
}
			
				
sub parseparam {
    my $var ;
    my @arg = @ARGV ;

    while(@arg) {
        $var = shift(@arg) ;
        $sessiondir = shift(@arg) if($var eq '-s') ;
        $old_url	= shift(@arg) if($var eq '-old') ;
        $new_url	= shift(@arg) if($var eq '-new') ;
    }
}