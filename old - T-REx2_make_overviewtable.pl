#!/usr/bin/env perl

# Anne de Jong
# 
# Trex2 RNAseq analysis 
# Make html overview table
# March 2020

use warnings;
use strict;
use lib "/data/molgentools/lib";
use anne_files ;
use anne_misc ;

my $sessiondir = ".";
my $experiment = "trex2";
my $url = "";
my $output = "00.overviewtable.html";
my $usage= "/data/trex2/T-REx_make_overviewtable.pl -s $sessiondir -url $url -e $experiment -o $output\n
e.g. /data/molgentools/trex/T-REx_make_overviewtable.pl -url 'trex2_results/T-REx' -e my_experiment -o 00.overviewtable.html
" ;



parseparam() ;
print "\nWrite overview table to $sessiondir/$output\n\n";
&anne_files::write_string("$sessiondir/$output" ,overview_table($experiment)) ;


sub overview_table {
	my $experiment = shift ;
	my $examples = "http://trex.molgenrug.nl/examples";
	
	# Style sheet
	my $results = "
	<style type=text/css>
		table.g2d_one {border-collapse:collapse;
			a:link { color: #D8F2F0; } /* unvisited link */
			a:visited {	color: #D8F2F0;	}/* visited link */
			a:hover { color: #D8F2F0; }/* mouse over link */
			a:active { color: #D8F2F0; }/* selected link */
		}
		td.head {
			font-size:18px;
			color: #D8F2F0; 
			background: #0F2D40; /* Dark green background */
			text-align: center; 
			height: 30px;
		}		
		td.a {
			color: #0F2D40;
			background: #113247; /* Dark */
			text-align: left; 
			padding: 4px;
			}
		td.b {
			color: #0F2D40;
			background: #194759; /* Midled */
			text-align: left; 
			padding: 4px;
		}
		td.c {
			color: #D8F2F0;
			background: #296B73; /* Light */
			text-align: left; 
			padding: 4px;
		}
		td.d {
			color: ##0F2D40;
			background: #475D67; /* Light */
			text-align: left; 
			padding: 4px;
		}
		
		a:link { color: #D8F2F0; } /* unvisited link */
		a:visited {	color: #D8F2F0;	}/* visited link */
		a:hover { color: #D8F2F0; }/* mouse over link */
		a:active { color: #D8F2F0; }/* selected link */

		</style> " ;
		
		
		
		# a:link { color: #193441; } /* unvisited link */
		# a:visited {	color: #0F1F26;	}/* visited link */
		# a:hover { color: #0F1F26; }/* mouse over link */
		# a:active { color: #193441; }/* selected link */
	
	# Global analysis
	my $json = "$url/$experiment.Container.json" ; 
	my $weblink = 'http://trex.molgenrug.nl/trex2_results';
	my $AdMiRE = 'http://trex.molgenrug.nl/admire/index.php?settings='."$sessiondir/$experiment.Container.json" ;
	$json =~ s/$weblink/$AdMiRE/ ;
	$AdMiRE =~ s/tmpdrive\/trex2/trex2_results/ ;
	print "Admire==========================================\n $AdMiRE \n\n\n";
	print "url==========================================\n $url \n\n\n";
	
	$results .= "<h1>T-REx2 results of $experiment</h1>
				&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='http://trex.molgenrug.nl/php/files/Tutorial_interpretation_T-REx_results.docx' target=_blank><font color=blue><b>Tutorial for Interpretation of T-REx Results</b></font></a>
				<hr>";
	$results .= " 
		<TABLE class=g2d_one><tr>
				<tr><td class=head><a href=$AdMiRE target=_blank><img src=http://trex.molgenrug.nl/images/Admire.png height=100 width=180> <br> AdMiRE: Interactive mining</a></td></tr>
		</table><br>

		<TABLE class=g2d_one><tr>
				<tr><td class=head colspan=8><b>Global analysis</b></font></td></tr>
				<tr><td class=a>
						<a href=Counts.txt>Counts data</a><br>
						<a href=Factors.txt>Factors</a><br>
						<a href=Contrasts.txt>Contrasts</a><br>
						<a href=Class.txt>Class data</a><br>
						<a href=results.zip>Download all as ZIP</a>
					</td>
					<td class=b><a href=$experiment.Library_size.png> <img src=$experiment.Library_size.png height=100 width=150> <br> Library size</a></td>
					<td class=a><a href=$experiment.colored_boxplot.png><img src=$experiment.colored_boxplot.png height=100 width=180> <br> Box plot of normalized signals</a></td>
					<td class=b><a href=$experiment.PCA_experiments.png> <img src=$experiment.PCA_experiments.png height=100 width=100> <br> PCA of Experiments</a></td>
					<td class=a><a href=$experiment.PCA_factors.png> <img src=$experiment.PCA_factors.png height=100 width=100> <br> PCA of Factors</a></td>
					<td class=b><a href=$experiment.PCA_genes.png> <img src=$experiment.PCA_genes.png height=100 width=100> <br> PCA of Genes</a></td>
				</tr><tr>
					<td class=d></td>
					<td class=d><a href=$experiment.Library_size.txt>Library size data</a></td>
					<td class=d><a href=$experiment.Library_size.txt>Library size data</a><br>
					<td class=d><a href=$experiment.PCA_experiments.txt>Data table</a></td>
					<td class=d><a href=$experiment.PCA_factors.txt>Data table</a><br>
					<td class=d><a href=$experiment.PCA_genes.txt>Data table</a></td>
				</tr><tr><td></td>
				</tr><tr>
					<td class=b><a href=$experiment.Heatmap_Mean_Signals.png><img src=$experiment.Heatmap_Mean_Signals.png height=100 width=120> <br> Heatmap of Mean Signals</a></td>
					<td class=a><a href=$experiment.Heatmap_Rank_Samples.png><img src=$experiment.Heatmap_Rank_Samples.png height=100 width=150> <br> Heatmap Rank-order Samples</a></td>
					<td class=b><a href=$experiment.Heatmap_Rank_Factors_1.png><img src=$experiment.Heatmap_Rank_Factors_1.png height=100 width=120> <br> Heatmap Rank-order Factors (1)</a></td>
					<td class=a><a href=$experiment.Heatmap_Rank_Factors_2.png><img src=$experiment.Heatmap_Rank_Factors_2.png height=100 width=120> <br> Heatmap Rank-order Factors (2)</a></td>
				</tr><tr>
					<td class=d><a href=$experiment.Mean_signals_of_replicates.txt>Mean Signals of Replicates</a><br>
					<td class=d><a href=$experiment.rank_table.txt>Rank-order Samples</a><br>
					<td class=d><a href=$experiment.rank_mean_table.txt>Rank-order Factors Data</a><br>
					<td class=d><a href=$experiment.rank_mean_table.txt>Rank-order Factors Data</a><br>
				</tr>
				</tr><tr>
					<td class=a><a href=$experiment.Heatmap_Mean_Signals.clusters.txt>ClusterIDs</a><br>
					<td class=a><a href=$experiment.Heatmap_Rank_Samples.clusters.txt>ClusterIDs</a><br>
					<td class=a><a href=$experiment.Heatmap_Rank_Factors_1.clusters.txt>ClusterIDs</a><br>
					<td class=a><a href=$experiment.Heatmap_Rank_Factors_2.clusters.txt>ClusterIDs</a><br>
				</tr>
		</TABLE><br> ";
	
	# Contrasts 
	$results .= " 
		<TABLE class=g2d_one><tr>
				<tr><td class=head colspan=8><b>Differential Gene Expression</b></font></td></tr>
				<tr>
					
					<td class=b><a href=$experiment.VolcanoPlots.html><img src=$examples/Example.Volcano_plot.png height=100 width=180> <br> Volcano plots</a></td>
					<td class=a><a href=$experiment.SigChanged_plots.html><img src=$examples/Example.Significant_Changed_Genes.png height=100 width=180> <br> Significant Changed Genes</a></td>
					<td class=b><a href=$experiment.MAplots.html><img src=$examples/Example.MAplot.png height=100 width=100> <br> MA plots</a></td>
				</tr><tr><td></td>	
				</tr><tr>
					<td class=b><a href=$experiment.Contrasts_overview.html><img src=$examples/Example.contrasts_each.png height=70 width=180> <br> Differential Expression<br>Each Contrast</a></td>
					<td class=a><a href=$experiment.TMEV_TOPHITS_logFC.html><img src=$examples/Example.contrasts_all.png height=100 width=180> <br> Differential Expression<br>All Contrasts</a></td>
					<td class=b><a href=$experiment.Heatmap_TopHits.png> <img src=$experiment.Heatmap_TopHits.png height=100 width=100> <br> Heatmap<br>TopHits</a></td>
					<td class=a><a href=$experiment.Heatmap_TopHits_HighFold.png> <img src=$experiment.Heatmap_TopHits_HighFold.png height=100 width=100> <br> Heatmap<br>TopHits, High Fold</a>
				</tr><tr>
					<td class=d><a href=$experiment.Contrasts_overview.txt>Data table</a></td>
					<td class=d><a href=$experiment.TMEV_TOPHITS_logFC.txt>Data table</a></td>
					<td class=d><a href=$experiment.Heatmap_TopHits.txt>-Heatmap TopHits data</a><br>
								<a href=$experiment.Heatmap_TopHits.clusters.txt>-Heatmap Cluster data</a></td>
					<td class=d><a href=$experiment.Heatmap_TopHits_HighFold.txt>-Heatmap High Fold data</a><br>
								<a href=$experiment.Heatmap_TopHits_HighFold.clusters.txt>-Heatmap Cluster data</a></td>
				</tr>				
		</TABLE><br> "; 
	
	# Experiment Analysis 	
		$results .= " 
		<TABLE class=g2d_one><tr>
				<tr><td class=head colspan=8><b>Experiment Analysis</b></font></td></tr>
				<tr>
					<td class=a><a href=$experiment.Correlation_matrix_experiments.png> <img src=$experiment.Correlation_matrix_experiments.png height=100 width=100> <br> Correlation matrix<br>of Experiments</a></td>
					<td class=b><a href=$experiment.Heatmap_EXPERIMENTS.png> <img src=$experiment.Heatmap_EXPERIMENTS.png height=100 width=100> <br>Heatmap<br>of Experiments</a></td>
					<td class=a><a href=$experiment.GeneNetwork_of_Contrasts.png> <img src=$experiment.GeneNetwork_of_Contrasts.png height=100 width=180> <br>Gene Network<br>Significant Changed Genes</a></td>
					<td class=b><a href=$experiment.GeneNetwork_of_Contrasts_high.png> <img src=$experiment.GeneNetwork_of_Contrasts_high.png height=100 width=180> <br>Gene Network<br>High Fold Changed Genes</a></td>
					<td class=a></td>
				</tr><tr>	
					<td class=d><a href=$experiment.Correlation_matrix_experiments.txt> Correlation data</a></td>
					<td class=d><a href=$experiment.Heatmap_EXPERIMENTS.txt>Heatmap data</a></td>
					<td class=d></td>
				</tr><tr><td></td>	
				</tr><tr>	
					<td class=b><a href=$experiment.contrasts_cohesion.html><img src=$examples/Example.contrasts_cohesion_table.png height=100 width=180> <br> Contrasts Cohesion<br>Genes overlap of Contrasts</a></td>
					<td class=a><a href=$experiment.contrasts_cohesion.genes.txt><img src=$examples/Example.contrasts_cohesion_genes.png height=100 width=180> <br> Contrasts Cohesion<br>Gene list</a></td>
					<td class=b><a href=$experiment.contrasts_cohesion.high.html><img src=$examples/Example.contrasts_cohesion_table.png height=100 width=180> <br> Contrasts Cohesion<br>High Fold Change</a></td>
					<td class=a><a href=$experiment.contrasts_cohesion.high.genes.txt><img src=$examples/Example.contrasts_cohesion_genes.png height=100 width=180> <br> Contrasts Cohesion<br>Gene List</a></td>
					<td class=b></td>
				</tr><tr>	
					<td class=d><a href=$experiment.contrasts_cohesion.txt>Matrix data</a></td>
					<td class=d></td>
					<td class=d><a href=$experiment.contrasts_cohesion.high.txt>Matrix data</a></td>
					<td class=d></td>
					<td class=d></td>
				</tr><tr><td></td>	
				</tr><tr>
					<td class=a><a href=$experiment.MedianFold.kmeans.html> <img src=$examples/Example.kmeans.png height=100 width=180> <br> Clustering<br>Median FOLD k-means</a></td>
					<td class=b><a href=$experiment.HighFold.kmeans.html> <img src=$examples/Example.kmeans.png height=100 width=180> <br> Clustering<br>High FOLD k-means</a></td>
					<td class=a><a href=$experiment.MeanSignal.kmeans.html> <img src=$examples/Example.kmeans.png height=100 width=180> <br> Clustering<br>Mean SIGNALS k-means</a></td>
					<td class=b><a href=$experiment.MedianFold.kmeans_Dendrogram.png> <img src=$experiment.MedianFold.kmeans_Dendrogram.png height=100 width=100> <br>Dendrogram<br>Median Fold k-means</a></td>
					<td class=a><a href=$experiment.HighFold.kmeans_Dendrogram.png> <img src=$experiment.HighFold.kmeans_Dendrogram.png height=100 width=100> <br>Dendrogram<br>High Fold k-means</a></td>
				</tr><tr>
					<td class=d><a href=$experiment.MedianFold.ClusterID.txt>ClusterIDs and data</a></td>
					<td class=d><a href=$experiment.HighFold.ClusterID.txt>ClusterIDs and data</a></td>
					<td class=d><a href=$experiment.MeanSignal.ClusterID.txt>ClusterIDs and data</a></td>
					<td class=d></td>
					<td class=d></td>
				</tr>
				
				</TABLE><br> " ;

	# Analysis of the Classes
		$results .= " 
		<TABLE class=g2d_one><tr>
				<tr><td class=head colspan=8><b>Analysis of the Classes</b></font></td></tr>
				<tr>
					<td class=a><a href=$experiment.ClassTophits.kmeans_plot_cluster.html> <img src=$examples/Example.kmeans.png height=100 width=180> <br>k-means RATIO<br> Class TopHits</a></td>
					<td class=b><a href=$experiment.ClassTophitsSignal.kmeans_plot_cluster.html> <img src=$examples/Example.kmeans.png height=100 width=180> <br>k-means SIGNAL<br> Class TopHits</a></td>
					<td class=a><a href=$experiment.Mean_signals_subClass_TOPHITS.html> <img src=$examples/Example.kmeans.png height=100 width=180> <br> Mean SIGNAL Plots<br>Each Class</a></td>
					<td class=b><a href=$experiment.Correlation_matrix_Class_genes.png> <img src=$experiment.Correlation_matrix_Class_genes.png height=100 width=100> <br> Correlation matrix<br>All Classes</a></td>
					<td class=a><a href=$experiment.SubClass_correlation_matrix.html> <img src=$examples/Example.matrices.png height=100 width=180> <br> Correlation matrices<br>Each Class</a></td>
				</tr><tr>               
					<td class=d><a href=$experiment.ClassTophits.ClusterID.txt>ClusterIDs and data</a></td>
					<td class=d><a href=$experiment.ClassTophitsSignal.ClusterID.txt>ClusterIDs and data</a></td>
					<td class=d><a href=$experiment.Mean_signals_log2_scaled_CLASSES.txt>Mean SIGNAL data</a></td>
					<td class=d><a href=$experiment.Correlation_matrix_Class_genes.txt>Correlation data</td>
					<td class=d></td>
				</tr><tr><td></td>	
				</tr><tr>               
					<td class=b><a href=$experiment.Heatmap_Class_TopHits.png> <img src=$experiment.Heatmap_Class_TopHits.png height=100 width=100> <br> Heatmap<br>Class TopHits</a></td>
					<td class=a><a href=$experiment.Heatmap_Class_Signals.png> <img src=$experiment.Heatmap_Class_Signals.png height=100 width=100> <br> Heatmap<br>Class Signals</a></td>
					<td class=b><a href=$experiment.ClassHeatmaps.html> <img src=$examples/Example.ClassHeatmaps.png height=100 width=100> <br> Heatmaps<br>Each Class Ratio</a></td>
					<td class=a><a href=$experiment.ClassHeatmapsSignal.html> <img src=$examples/Example.ClassHeatmapsSignals.png height=100 width=100> <br> Heatmaps<br>Each Class Signals</a></td>
					<td class=b><a href=$experiment.ClassTophits.kmeans_MDS.png> <img src=$experiment.ClassTophits.kmeans_MDS.png height=100 width=100> <br> MDS plot of k-means<br>of ClassTophits</a></td>
				</tr><tr>               
					<td class=d><a href=$experiment.Heatmap_Class_TopHits.txt>- Heatmap data</a><br>
								<a href=$experiment.Heatmap_Class_TopHits.clusters.txt>- ClusterIDs</a></td>
					<td class=d><a href=$experiment.Heatmap_Class_Signals.txt>- Heatmap Signal data</a><br>
								<a href=$experiment.Heatmap_Class_Signals.clusters.txt>- ClusterIDs</a></td>
					<td class=d><a href=$experiment.ClassHeatmaps.clusters.list.html>Link to ClusterIDs</a></td>
					<td class=d><a href=$experiment.ClassHeatmapsSignal.clusters.list.html>Link to ClusterIDs</a></td>
					<td class=d></td>
				</tr>
		</TABLE><br> "; 
				


	return $results ;
}

sub parseparam {
    my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
        $sessiondir = shift(@arg) if($var eq '-s') ;
        $experiment = shift(@arg) if($var eq '-e') ;
        $url		= shift(@arg) if($var eq '-url') ;
		$output		= shift(@arg) if($var eq '-o') ;
    }
    die $usage if (!$experiment) ;
}
