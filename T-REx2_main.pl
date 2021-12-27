#!/usr/bin/env perl

# Anne de Jong
# 
# T-REx.v2 main perl script
# 2020 March


use warnings;
use strict;
use lib "/data/molgentools/lib";
use anne_files ;
use anne_misc ;


# -----------------------------------------------------  parameters  ---------------------------------------
my $sessiondir      = './';
my $experiment      = 'trex2';
my $counts_file     = 'Counts.txt';
my $factors_file    = 'Factors.txt';
my $contrast_file   = 'Contrasts.txt';
my $class_file      = 'Class.txt';
my $annotation_file = 'Annotation.txt';
my $chkHeatmap      = 'true';
my $chkCluster      = 'true';
my $chkDE           = 'true';
my $chkRep          = 'true';
my $tmp ;
my $trex_programdir = '/data/trex2';
my $url             = 'http://trex.molgenrug.nl/trex2_results';
my $error = 'ok' ;
my $usage = "/data/trex2/T-REx2_main.pl
				-s			Sessiondir [default=current folder]
				-counts		TAB delimted table with fpkm/rpkm/counts values of all experiments [default=$counts_file]
				-factor		Factor file [default=$factors_file]
				-contrasts	Contrasts file [default=$contrast_file]
				-class		class file [default=$class_file]
				-name		Experiment name will be used to prefix the results [default= my_experiment]
				-heatmap    Generate heatmaps                           [default=true]
				-cluster    Perform Clustering                          [default=true]
				-de         Perform Differential Expression analysis    [default=true]
				-rep        Counts table contain replicates             [default=true]
				-url		Sessionfolder for the webserver
		e.g.  /data/trex2/T-REx2_main.pl -rep false 
/data/trex2/T-REx2_main.pl -s /tmp/TREx2/94.212.189.120.6hm59kmil0sacdu390sf8r6n33.887 -heatmap true -cluster true -de true -rep true		
" ;

&parseparam() ;

my @Rcolors = anne_files::read_lines("$trex_programdir/Rcolors.list") ;

$experiment = 'my_experiment' if ($experiment eq '') ;  # if it comes empty from the webserver
$experiment =~ s/\s+/\_/g ; 
$experiment =~ s/[^a-zA-Z0-9 ]//g ; # remove all non a-Z


print "
============================================
counts_file      =$counts_file    
factors_file     =$factors_file   
contrast_file    =$contrast_file  
class_file       =$class_file     
annotation_file  =$annotation_file

sessiondir=$sessiondir
chkHeatmap=$chkHeatmap 
chkCluster=$chkCluster 
chkDE     =$chkDE      
chkRep    =$chkRep    

===========================================\n ";

anne_files::write_string("$sessiondir/00.perl_wrapper.log","
============================================
counts_file      = $counts_file   
factors_file     = $factors_file
contrast_file    = $contrast_file
class_file       = $class_file  
annotation_file  = $annotation_file

sessiondir = $sessiondir
chkHeatmap = $chkHeatmap 
chkCluster = $chkCluster 
chkDE      = $chkDE
chkRep     = $chkRep

===========================================\n ");



# 1. clean the tables from empty lines and DOS chars

	# handle Mac exported data: replace ^M by \n
	print "Clean files\n" ;
	anne_files::mac2linux("$sessiondir/$counts_file") ;
	anne_files::mac2linux("$sessiondir/$factors_file");
	anne_files::mac2linux("$sessiondir/$contrast_file");
	anne_files::mac2linux("$sessiondir/$class_file");
	
	# add proper header
	change_header_counts_file("$sessiondir/$counts_file") ;
	
	# remove illegal chars
	clean_table_4_R("$sessiondir/$counts_file") ;
	clean_table_4_R("$sessiondir/$factors_file");
	clean_table_4_R("$sessiondir/$contrast_file");
	clean_table_4_R("$sessiondir/$class_file");
#	
#	# check Counts input file 
#	print "Check file formats\n" ;
#	report_error() if (!chk_rpkm_duplicates()) ;
#
#	# check factor input file 
#	report_error() if (!chk_factornames()) ;
#	report_error() if (!chk_experiment_order()) ;
#
#	# check the class input file
#	report_error() if (!chk_class_groupsize())  ;
#	report_error() if (!chk_class_colors()) ;
#	report_error() if (!chk_class_double_geneID()) ;
#
#	# check the Contrasts file
#	report_error() if (!chk_contrasts()) ;
#	report_error() if (!chk_contrasts_names()) ;
#

# ----------------------------------------------------------------  functions  -------------------------------------------------------------------------------------------------------
	


# 2. run R, TREx2 
	print "Start R script\n" ;
	my @lines ;
	@lines = anne_files::read_lines("$sessiondir/$class_file") ;
	my $chkClass = 'true' ;
	$chkClass = 'false' if (scalar @lines < 5) ; # disable class if empty or low number of classes
	
	my $command = "R --vanilla --slave --args $sessiondir $chkHeatmap $chkCluster $chkDE $chkRep $chkClass< $trex_programdir/T-REx2.R 2>>$sessiondir/00.R.sessionprogress.log" ;  
	# my $command = "R --vanilla --slave --args $sessiondir $chkHeatmap $chkCluster $chkDE $chkRep $chkClass< $trex_programdir/T-REx2.R" ;  
	print "$command\n";
	anne_files::append_lines("$sessiondir/00.perl_wrapper.log","Start R script\n$command\n");
	system($command) ;



# 3. Check if run was successful
	print "Checking results\n" ;

	@lines = anne_files::read_lines("$sessiondir/00.R.sessionprogress.log") ;
	if (grep {$_ eq 'Execution halted'} @lines) {
		anne_files::write_string("$sessiondir/00.R.ERROR", "ERROR in R") ; # tell the web-server that an error occurred 
		anne_files::append_lines("$sessiondir/00.R.log", "<br>ERROR in R, please check your input files<br><br>ERROR report:<br>") ;
		anne_files::append_lines("$sessiondir/00.R.log", @lines) ;
		
	} else {
		# 4. zip the results
			# - $tmp = "zip -j $sessiondir/results.zip $sessiondir/$experiment*" ;
			# - system($tmp) ;
			# - $tmp = "zip -j $sessiondir/results.zip $sessiondir/00.*" ;
			# - system($tmp) ;

		# 5. Results in html overview tables

			# over view for clusterid links in multiplots
			overview_clusterIDs('ClassHeatmaps') ;
			overview_clusterIDs('ClassHeatmapsSignal') ;
			
			# Cohesion of Contrasts is not written in R but in Perl
			$tmp = "$trex_programdir/T-REx2_cohesion_contrasts.pl -s $sessiondir -i $sessiondir/$experiment.GeneNetwork_EdgeList.txt -o $experiment.contrasts_cohesion" ;
			print "$tmp\n";
			system($tmp) ;
			$tmp = "$trex_programdir/T-REx2_cohesion_contrasts.pl -s $sessiondir -i $sessiondir/$experiment.GeneNetwork_EdgeList.txt -o $experiment.contrasts_cohesion.high -logFC 2 -pvalue 0.01" ;
			print "$tmp\n";
			system($tmp) ;
			
		# 6. write the results overview table
			anne_files::write_string("$sessiondir/00.experiment_name.txt", $experiment) ; # This file is used for downstream scripts 

	}		


# 8. tell the webserver that the session is finished	
	&anne_files::write_string("$sessiondir/sessionstop" , "run done") ;

	


# ---------------------------------------------------------------------------------------- functions ----------------------------------------------------------------------------------------

sub report_error {
	# report an error and terminate the program
	anne_files::write_string("$sessiondir/00.R.ERROR", $error) ; # tell the web-server that an error occurred 
	anne_files::write_string("$sessiondir/00.R.log", "<br>$error<br>") ; # report the error
	anne_files::write_string("$sessiondir/sessionstop" , "run done") ;
	print "$error\n";
	exit ;
}	

sub overview_clusterIDs {
	my $prefix = shift ;
	my $outputfile = "$sessiondir/$experiment.$prefix.clusters.list.html" ;
	my $result = "GeneIDs and their ClusterID can be found via the links below<hr>\n" ;
	my @lines = anne_files::read_lines("$sessiondir/$experiment.$prefix.clusters.filenames.txt");
	foreach my $line (@lines) {
		if ($line =~ m/(.+)\s(.+)/) {
			my $file = $1 ;
			my $name = $2 ;
			$result .= "<a href=$url/$file>Class $name</a><br>\n";
		}	
	}
	anne_files::write_string($outputfile, $result);
	
}

sub json_cleanup {
	# in R some wrong quote marks are inserted, could not solve it in R, so solved it in perl
	my $filename = shift ;
	my $x = '(\[|\]|\{|\}|\:|\,)' ;
	my $line = join('', anne_files::read_lines($filename)).'x' ;
	$line =~ s/\\//g ;
	my $result =  substr($line, 0, 1) ;
	for (my $i=0; $i<length($line)-2; $i++) {
		if (substr($line, $i, 3) !~ m/$x\"$x/g) { $result .= substr($line, $i+1, 1) ; } 
	}
	# some returns for readability
	$result =~ s/\{/\{\n/g ;
	$result =~ s/\}/\n\}/g ;
	$result =~ s/\,/\,\n/g ;	
	anne_files::write_string($filename, $result) ;

}


sub change_header_counts_file {
	# the R script need the word 'key' in the header of the first column of the rpkm file
	my $filename = shift ;
	my @lines = anne_files::read_lines($filename);
	my @items = split (/\t/, $lines[0]);
	$items[0] = 'key' ;
	$lines[0] = join "\t", @items ;
	&anne_files::write_lines($filename, @lines) ;
}


sub clean_table_4_R {
	my $filename = shift ;
	my @lines = anne_files::read_lines($filename);
	my @results ;
	foreach my $line (@lines) {
		#$line =~ s/(\n|\r|\x0d)//g;		 				# remove line breaks etc
		#$line =~ tr/\cM//d ;
		$line =~ s/\s+$// ;								# remove tailing spaces end of line
		$line =~ s/\s+\t/\t/ ;							# remove tailing spaces before tab
		$line =~ s/(\,|\!|\\|\/|\||\:|\"|\?|\<|\>)//g ;	# remove illegal chars for R, also remove the thousand 'comma' in numbers which causes error in R 
		$line =~ s/\ /\_/g ;							# replace space by _
		while ($line =~ m/\t$/ ) { $line =~ s/\t$// ; }	# remove empty column(s) on the end of row
		if ($line ne '') {
			my @items = split ('\t', $line) ;
			#@items = anne_misc::de_space_array(@items) ; # TODO
			chomp (@items) ;	# remove leading and tailing spaces
			push @results, join("\t", @items) ;
		}
	}
	&anne_files::write_lines($filename, @results) ;
	
}	

sub chk_experiment_order {
	# check if the headers are identical except for the first column
	my @rpkm_header = anne_files::get_table_header("$sessiondir/$counts_file") ;
	my @factor_lines = anne_files::read_lines("$sessiondir/$factors_file") ;
	my $len_rpkm = scalar (@rpkm_header) ;
	my $len_factor = scalar (@factor_lines) ;
	my $result = 1 ;
	if ( $len_rpkm != $len_factor ) {
		$result = 0 ;
	} else {	
		for (my $i = 1; $i < $len_rpkm; $i++) {	
			$result = 0 if ($factor_lines[$i] !~ m/^$rpkm_header[$i]/ ) ;  
			print "$rpkm_header[$i] == $factor_lines[$i] $result\n";
		}
	}
	$error = "ERROR: The order of the Experiment names in the Factors file must have the same order as the column headers in the RPKM file. <br>\nPlease correct your Factors and/or the RPKM file<br>\n" ;
	return $result ;

}

sub chk_contrasts_names {
	# get contrasts names
	my @contrasts ;
	my @lines = anne_files::read_lines("$sessiondir/$contrast_file");
	foreach my $line (@lines) {	push @contrasts, split '-', $line ; }
	# get factors names
	my @factors ;
	@lines = anne_files::read_lines("$sessiondir/$factors_file");
	foreach my $line (@lines) {	
		my @items = split '\t', $line ;
		shift @items ;
		push @factors, join '.', @items ; 
	}
	# check if contrasts names match the factors
	my $result = 1 ;
	foreach my $contrast (@contrasts) {
		my $found = grep(/$contrast/, @factors);
		if ($found == 0) { 
			$error = "Error in the Contrasts file: \'$contrast\' not found in Factors file"; 
			$result = 0 ;
			last ;
		}
	}
	return $result ;
}	
	

sub chk_contrasts {
	my @lines = anne_files::read_lines("$sessiondir/$contrast_file");
	# clean the contrast file
	my $result = 1 ;
	my $count ;
	$error = 'none';
	my @cleanlines ;
	foreach my $line (@lines) {	
		# clean the line; remove tabs and spaces
		$line =~ s/(\t+|\s+)//g ;
		push @cleanlines, $line if ($line ne '') ; 
		# check if multiple - are in the line
		if ( ($count = () = $line =~ /\-/g)!=1 ) { 
			$error = "Error in the Contrasts file: there should be one '-' per line: $line"; 
			$result = 0 ;
		}
	}
	anne_files::read_lines(@cleanlines, "$sessiondir/$contrast_file");
	return $result ;

}	
	
sub chk_factornames {
	# factor names starting with a number is not allowed in R
	my @lines = anne_files::read_lines("$sessiondir/$factors_file");
	my $result = 1 ;
	foreach my $line (@lines) {	$result = 0 if ($line =~ m/\t\d+/g) ; }	
	$error = "ERROR: Factor starts with digid number. This is not allowed in R, please correct your Factors and its cognate RPKM table headers<br>\n" ;
	return $result ;
}

	
sub chk_class_groupsize {
	# Check if the groups contain at least 3 items
	my @lines = anne_files::read_lines("$sessiondir/$class_file") ;
	# make a hash of all groups with array of locus tags
	my %groups ;
	my $header = 0;
	foreach my $line (@lines) {	
		my @items = split "\t", $line ;
		if ($header) { 
			push @{ $groups{$items[2]} }, $items[0] ; 
		} else {
			$header = 1; }
	}	
	# check if there is a group with less then 3 locus tags
	my $result = 1 ;
	my @smallgroup ;
	foreach my $key (keys %groups ) {
		if (scalar @{$groups{$key}} < 3) {
			$result = 0 ;
			push @smallgroup, $key ;
		}	
	}	
	if (!$result) {
		$error = "ERROR in CLASS file: Group size must be >=3, please correct the group that contain less than 3 genes: <br>\n" ;
		$error .= join ( "<br>\n", @smallgroup ) ;
	}	
	return $result ;
}	
	

sub chk_class_colors {
	# Check if the color is a valid R-color name
	my @lines = anne_files::read_lines("$sessiondir/$class_file") ;
	my $header = 1;
	my $result = 1 ;
	my @wrongcolors ;
	foreach my $line (@lines) {	
		my @items = split "\t", $line ;
		if ($header) { $header = 0 ; 
		} else {
			if (not grep(/$items[1]/, @Rcolors)) {
				$result = 0   ;
				push @wrongcolors, $items[1] ; 
			}	
		}	
	}
	if (!$result) {
		$error = "ERROR in CLASS file: Wrong color name(s) found, see tutorial for proper color names: <br>\n" ;
		$error .= join ( "<br>\n", @wrongcolors ) ;
	}	
	return $result ;
}

sub chk_class_double_geneID {
	# Check GeneID is not used for more than 1 group
	my $result = 1 ;
	my @lines = anne_files::read_lines("$sessiondir/$class_file") ;
	my $header = 1;
	# get all geneIDs
	my @geneIDs ;
	foreach my $line (@lines) {	
		if ($header) { $header = 0 ; # skip header
		} else {
			my @items = split "\t", $line ;
			if (scalar(@items) > 2 ) { push @geneIDs, $items[0]; }
		}
	}
	my $last = '';
	my @doubles ;
	foreach my $geneID (sort @geneIDs) {
		if ($geneID eq $last) { push @doubles, $geneID ; }
		$last = $geneID ;
	}	
	if ( scalar(@doubles) >0 ) {
		$error  = "ERROR in CLASS file: Some GeneIDs are member of more than 1 group. Please define a combined group or remove the double GeneIDs<br>\n";
		$error .= "The following geneIDs were found more than 1 time: <br>\n" ;
		$error .= join ( "<br>\n", @doubles ) ;
		$result = 0 ;
		
	}	
	return $result ;
}



sub chk_rpkm_duplicates {
	# Check GeneID is not used for more than 1 group
	my $result = 1 ;
	my @lines = anne_files::read_lines("$sessiondir/$counts_file") ;
	my @geneIDs ;
	my @duplicates ;
	foreach my $line (@lines) {	
		my @items = split "\t", $line ;
		if (scalar(@items) > 2 ) { 
			if ( grep { $_ eq $items[0] } @geneIDs ) { push @duplicates, $items[0] ; }
			push @geneIDs, $items[0]; 
		}

	}
	if ( scalar(@duplicates) >0 ) {
		$error  = "ERROR in Counts file: Duplicate geneIDs found. Please correct this.<br>\n";
		$error .= "The following geneIDs were found more than 1 time: <br>\n" ;
		$error .= join ( "<br>\n", @duplicates ) ;
		$result = 0 ;
		
	}	
	return $result ;
}




sub parseparam {
    my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
        die $usage if($var eq '-h' or $var eq '-help')  ;
		$sessiondir 	= shift(@arg) if($var eq '-s') ;
        $counts_file 		= shift(@arg) if($var eq '-counts') ;
		$factors_file	= shift(@arg) if($var eq '-factor') ;
		$class_file		= shift(@arg) if($var eq '-class') ;
		$contrast_file	= shift(@arg) if($var eq '-contrasts') ; 
		$annotation_file= shift(@arg) if($var eq '-annotation') ; 
        $experiment 	= shift(@arg) if($var eq '-name') ;
		
		$chkHeatmap    = shift(@arg) if($var eq '-heatmap') ; 
		$chkCluster    = shift(@arg) if($var eq '-cluster') ; 
		$chkDE         = shift(@arg) if($var eq '-de') ; 
		$chkRep        = shift(@arg) if($var eq '-rep') ; 

		$url			= shift(@arg) if($var eq '-url') ;

		
    }
    
}
