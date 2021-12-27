#!/usr/bin/env perl

# Anne de Jong
# 
# cohesion of contrasts on the basis of the network interaction file
# Downstream processing of the T-REx pipeline
# March 2020

use strict ;
use warnings ;
use lib "/data/molgentools/lib";
use anne_files;

# -----------------------------------------------------  parameters  ---------------------------------------
my $sessiondir = './';
my $inputfile ;
my $output_prefix = 'my_result.contrasts_cohesion' ;
my $url_base = 'http://trex.molgenrug.nl/trex2_results';
my $logFC = 1 ;
my $pvalue = 0.05 ;
my $usage = "./T-REx_cohesion_contrasts.pl
				-s sessiondir [default=current folder]
				-i network table [ tab delimited; Contrast, Gene, logFC, pvalue ]
				-logFC use only values >= logFC [default = $logFC]
				-pvalue use only values <= pvalue [default = $pvalue]
				-o output file [default=$output_prefix] 
		e.g.  ./T-REx_cohesion_contrasts.pl -i GenNetwork.txt

result consists of 3 files
1. my_result.contrasts_cohesion.genes.txt	=> a text file containing the gene names associated with the contrasts combination
2. my_result.contrasts_cohesion.txt			=> contrasts combination in table format

" ;

&parseparam() ;
my $classStyle = "<html>
				<head>
					<style type=text/css>
						.verticalText
						{
						text-align: center;
						vertical-align: middle;
						width: 20px;
						margin: 0px;
						padding: 0px;
						padding-left: 3px;
						padding-right: 3px;
						padding-top: 10px;
						white-space: nowrap;
						-webkit-transform: rotate(-40deg);
						-moz-transform: rotate(-40deg);
						}
						 td#main { border:0px solid black;border-top:0;text-align:right;width:10px;};
					</style>
					<script type=text/javascript src=zepto.js></script>
					<script type=text/javascript src=foundation.min.js></script>
					<script type=text/javascript src=foundation.tooltips.js></script>
				</head>\n
				<body>
				" ;



# -------------------------------------------------------------------------------------- main ------------------------------------------------------------------------------------------


my $url = $sessiondir ;
$url =~ s/\/tmpdrive\/trex2/$url_base/ ;


read_networkfile () ;
cohesion_2_html();



# ---------------------------------------------------------------------------------------- functions ----------------------------------------------------------------------------------------

sub cohesion_2_html {
	# read the Cohesion file
	my @lines = anne_files::read_lines("$sessiondir/$output_prefix.txt") ;
	my $color_fill	= "#4682B4" ;
	my $color_empty	= "#ADD8E6" ;
	my $color_bg = "#4682B4";
	my $header = 1 ;
	my @html = "$classStyle\n<b>Cohesion of Contrasts</b><br>
				<br><br><br><br><br><br>
				<table id=main>";
	foreach my $line (@lines) {
		my @items = split /\t/, $line ;
		my $colcount = scalar @items ;
		if ($header) {  # the header
			$header = 0 ;
			push @html, "<tr>" ;
			for (my $i=1; $i<$colcount-1; $i++) { 
				push @html, "<th><div class=verticalText>$items[$i]</div></th>"; 
			}
			push @html, "<th><div class=verticalText>Download Genes</div></th><th><div class=verticalText>Group_ID</div></th><th align=left>&nbsp;&nbsp;&nbsp;&nbsp;Cohesion_ID</th></tr>\n</div><tbody>\n" ;
		} else {	# the table body
			push @html, "\t<tr bgcolor=$color_bg>" ;
			for (my $i=1; $i<$colcount-2; $i++) {  # all the contrasts
				if ($items[$i] eq '') {
					push @html, "<td bgcolor=$color_empty> </td>";
				} else {	
					push @html, "<td bgcolor=$color_fill><span title=$items[$i]><font color=$color_fill>XXX</font></span></td>";
				}	
			}	
			push @html, "<td bgcolor=$color_empty>$items[$colcount-2]</td>  
						<td bgcolor=$color_empty> <a href=$url/$output_prefix.groupID_$items[$colcount-1].txt download>genes</a> </td>
						<td bgcolor=$color_empty><span title=$items[0]>$items[$colcount-1]</span></td>
						<td bgcolor=$color_empty>$items[0]</td></tr>\n" ;
		}
	}
	push @html, "</TABLE>\n" ;
	anne_files::write_lines("$sessiondir/$output_prefix.html", @html) ;
}


sub read_networkfile {
	my @lines = anne_files::read_lines($inputfile) ;
	my $header = 1 ;
	my @contrasts  ;
	my %network ;
	my %combinations ;

	# 1. Read the tab delimited table: contrast, locus
		foreach my $line (@lines) {
			if ($header) {  # skip header
				$header = 0 ;
			} else {	
				my @items = split "\t", $line ;
				if ( abs($items[2]) >= $logFC and $items[3] <= $pvalue ) {  # use only values between the cutoff thresholds
					my $contrast = $items[0] ;
					my $locus    = $items[1] ;
					push @{$network{$locus}{single_contrast}}, $contrast ; 
					push @contrasts, $contrast ;
				} 	
			}
		}
		@contrasts = sort (anne_files::unique_array(@contrasts)) ;
	
	# 2. merge the contrasts to one node name and count the number of genes in this combination of contrasts
		foreach my $locus ( sort keys %network) {
			$network{$locus}{combined_contrast} =  (join "|", sort @{$network{$locus}{single_contrast}}) ;
			if (defined($combinations{$network{$locus}{combined_contrast}})) {
				$combinations{$network{$locus}{combined_contrast}}{GeneCount}++ ;
			} else {
				$combinations{$network{$locus}{combined_contrast}}{GeneCount} = 1 ;
			}
		}

		
	# 3. Add group_ID (just numbering the combined_contrasts)
		my $group_ID = 0 ;
		my $prev_combined_contrast = '';
		foreach my $Combined_Contrasts (keys %combinations) {
			$group_ID++ if ($prev_combined_contrast ne $Combined_Contrasts) ;
			$combinations{$Combined_Contrasts}{group_ID} = $group_ID ;
			$prev_combined_contrast = $Combined_Contrasts ;
		}

		
# old	# 4. write the locus tags and their combined_contrasts to file 
# old		@lines = "Group_ID\tGenes\tCohesion_ID" ; # add header
# old		foreach my $locus (sort { $network{$a}{combined_contrast} cmp $network{$b}{combined_contrast} || $a cmp $b } keys %network) {  
# old			my $group_ID = $combinations{$network{$locus}{combined_contrast}}{group_ID} ;
# old			push @lines, "$group_ID\t$locus\t$network{$locus}{combined_contrast}" ;
# old		}
# old		anne_files::write_lines("$sessiondir/$output_prefix.genes.txt", @lines) ;

	# 4. write the locus tags and their combined_contrasts to file // April 2020 changed column order
		@lines = "GeneID\tGroupID\tCohesionID" ; # add header
		foreach my $locus (sort { $network{$a}{combined_contrast} cmp $network{$b}{combined_contrast} || $a cmp $b } keys %network) {  
			my $group_ID = $combinations{$network{$locus}{combined_contrast}}{group_ID} ;
			push @lines, "$locus\t$group_ID\t$network{$locus}{combined_contrast}" ;
		}
		anne_files::write_lines("$sessiondir/$output_prefix.genes.txt", @lines) ;

	# 5. write genelists for download
		foreach my $Combined_Contrasts (keys %combinations) {
			my @keys = grep { $network{$_}{combined_contrast} eq $Combined_Contrasts } keys %network;
			#foreach my $key (sort @keys) {	print "$key\t$Combined_Contrasts\t$combinations{$Combined_Contrasts}{group_ID}\n"; 	}
			anne_files::write_lines("$sessiondir/$output_prefix.groupID_$combinations{$Combined_Contrasts}{group_ID}.txt", @keys) ;
		}	
		
		
	# 6. Make a 2D hash of %combinations; Combined_Contrasts x Contrast  
		foreach my $Combined_Contrasts (keys %combinations) {
			foreach my $contrast (@contrasts) {
				if ($contrast =~ m/$Combined_Contrasts/) {
					$combinations{$Combined_Contrasts}{$contrast} = $contrast ;
				} else {
					$combinations{$Combined_Contrasts}{$contrast} = '' ;
				}
			}
		}
		
	
	# 7. Make Overview table including GeneCounts and GroupID
		@lines = () ;
		foreach my $key (sort { $combinations{$b}{GeneCount} <=> $combinations{$a}{GeneCount} } keys %combinations) {  # Reverse sort on GeneCount
			my @row = $key ;
			foreach my $contrast (@contrasts) {
				push @row, $combinations{$key}{$contrast} ; 
			}	
			push @row, $combinations{$key}{GeneCount},$combinations{$key}{group_ID}  ;
			push @lines, (join "\t",@row);
		}
		unshift @lines, "Contrasts_Cohesion\t".(join "\t", @contrasts)."\tGeneCount\tGroup_ID" ;  # add header
		anne_files::write_lines("$sessiondir/$output_prefix.txt", @lines) ;
	
}


sub parseparam {
    my @arg = @ARGV ;
    while(@arg) {
        my $var = shift(@arg) ;
        $sessiondir = shift(@arg) if($var eq '-s') ;
        $inputfile 	= shift(@arg) if($var eq '-i') ;
        $logFC 		= shift(@arg) if($var eq '-logFC') ;
        $pvalue 	= shift(@arg) if($var eq '-pvalue') ;
        $output_prefix = shift(@arg) if($var eq '-o') ;
    }
    die "No inputfile found\n$usage" if (!$inputfile) ;

}



