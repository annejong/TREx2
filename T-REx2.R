# Anne de Jong, Aug 2015, T-REx main R script 
# University of Groningen
# Molecular Genetics
# the Netherlands

# 2015 Aug : First release for publication in BMC genomics
# 2015 Aug : Add, MDS plot for average values

# 2020 March 21 update server to version TREx2

  # Non replicate support
  # Class file optional (was obligatory)
  # Heatmap optional for fast analysis





##############################################################################################
##                        Package installer                                                 ##
##############################################################################################

# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# if(!require("statmod"))     { BiocManager::install("statmod",     ask=FALSE, dependencies = TRUE)}
# if(!require("reshape2"))    { BiocManager::install("reshape2",    ask=FALSE, dependencies = TRUE)}
# if(!require("ggcorrplot"))  { BiocManager::install("ggcorrplot",  ask=FALSE, dependencies = TRUE)}
# if(!require("ggplot2"))     { BiocManager::install("ggplot2",     ask=FALSE, dependencies = TRUE)}
# if(!require("tidyr"))       { BiocManager::install("tidyr",       ask=FALSE, dependencies = TRUE)}
# if(!require("cluster"))     { BiocManager::install("cluster",     ask=FALSE, dependencies = TRUE)}
# if(!require("plyr"))        { BiocManager::install("plyr",        ask=FALSE, dependencies = TRUE)}
# if(!require("grid"))        { BiocManager::install("grid",        ask=FALSE, dependencies = TRUE)}
# if(!require("gplots"))      { BiocManager::install("gplots",      ask=FALSE, dependencies = TRUE)}
# if(!require("jsonlite"))    { BiocManager::install("jsonlite",    ask=FALSE, dependencies = TRUE)}
# if(!require("RColorBrewer")){ BiocManager::install("RColorBrewer",ask=FALSE, dependencies = TRUE)}
# if(!require("limma"))       { BiocManager::install("limma",       ask=FALSE, dependencies = TRUE)}
# if(!require("edgeR"))       { BiocManager::install("edgeR",       ask=FALSE, dependencies = TRUE)}
# if(!require("binr"))        { BiocManager::install("binr",        ask=FALSE, dependencies = TRUE)} 



##############################################################################################
##                        Libraries                                                         ##
##############################################################################################

library(statmod)
library(reshape2)
library(ggcorrplot)
library(ggplot2)
library(cluster)
library(plyr)
library(grid)
library(gplots)
library(jsonlite)
library(RColorBrewer)
library(binr)
library(limma)
library(edgeR)
library(tidyr)

# lib_RNAseq <-  "G:\\My Drive\\WERK\\T-REx_local\\T-REx_functions.R"
lib_RNAseq <- '/data/trex2/T-REx2_functions.R'
source(lib_RNAseq)



##############################################################################################
##                        For Local use                                                     ##
##############################################################################################

# sessiondir="G:\\My Drive\\WERK\\PROJECTS\\Anna\\Course"
# experiment="Run01"
# counts_file="Counts.txt"
# factors_file="Factors.txt"
# class_file="Class.txt" 
# contrast_file="Contrasts.txt" 

##############################################################################################
##                        Analysis options                                                  ##
##############################################################################################

option.heatmaps=TRUE
option.classes=TRUE
option.clustering=TRUE
option.contrasts=TRUE
option.replicates=TRUE



##############################################################################################
##                        Web server or Command line                                        ##
##############################################################################################


# Webserver assume default filenames
experiment="trex2"
counts_file="Counts.txt"
factors_file="Factors.txt"
class_file="Class.txt" 
contrast_file="Contrasts.txt" 



# load parameters from command line or call from the webserver
# from perl:    R --vanilla --slave --args $sessiondir $chkHeatmap $chkCluster $chkDE $chkRep $chkClass
# command line: R --vanilla --slave --args true true true true true
sessiondir <- commandArgs()[5]
if (commandArgs()[6] == 'false') { option.heatmaps=FALSE }
if (commandArgs()[7] == 'false') { option.clustering=FALSE }
if (commandArgs()[8] == 'false') { option.contrasts=FALSE }
if (commandArgs()[9] == 'false') { option.replicates=FALSE }
if (commandArgs()[10] == 'false') { option.classes=FALSE }

# make web folder from the local folder
websessiondir <- gsub("/tmpdrive/trex2", "/trex2_results", sessiondir)
#websessiondir <- "trex2_results"




##############################################################################################
##                        graphics settings                                                 ##
##############################################################################################


pointsizes <- 48
# Image resolution
# 720p
x_720p <- 1280
y_720p <- 720
# FullHD
x_1k <- 1920
y_1k <- 1080
# QHD
x_2k <- 2560
y_2k <- 1440
# UHD 4k
x_4k <- 3840
y_4k <- 2160
# FUHD
x_8k <- 7680
y_8k <- 4320


  
##############################################################################################
##                        Init run                                                          ##
##############################################################################################
  
setwd(sessiondir)
#dir()

# comment SINK if you run this script locally
# sink("00.R.remarks.log", append=TRUE, split=FALSE)

# container for arguments and data tables in jSON format
jSON <- NULL
# creating a menus list for JSON
json.settings <- NULL
json.settings$descriptionsFile <- unbox(paste(websessiondir, "Annotation.txt", sep="/"))
json.settings$descriptionSeparator <- unbox("\t")
json.menus <- list()
  
  
errorfile <- "00.R.sessionprogress.log"

logbook <- c()
logbook <- write_log(logbook, "Pipeline started")
logbook <- write_log(logbook, "Functions loaded")





##############################################################################################
##                        Loading the data                                                  ##
##############################################################################################



######################## 1. Load Counts data ################################################

# old - cut2 <- function(x, breaks) {
#   labels <- paste0("(",  breaks[-length(breaks)], ",", breaks[-1L], "]")
#   return(factor(labels[findInterval(x, breaks)], levels=labels))
# }


	
	unlink(errorfile)
	logbook <- write_log(logbook, "# 1. Analyse Count data")
	x <- read.delim(counts_file,sep="\t", row.names="key")
	has.no.na <- apply(x, 1, function(x){any(is.na(x))})  # remove rows with missing values
	x <- x[!has.no.na,]

	filename <- paste(experiment, "00.counts_in R.txt", sep=".")
	write.table(x, file=filename , quote=F, sep="\t", row.names=TRUE, col.names=NA)

	# save quantile table
	filename <- paste(experiment, "quantile_table.txt", sep=".")
	write.table(apply(x,2,quantile), file=filename , quote=F, sep="\t", row.names=TRUE, col.names=NA)

	logbook <- write_log(logbook, "----> Estimate linear range of counts")
	# Make 12 bins with max count value of each bin 
	bins12 <- apply(x,2,function(y){bins(as.integer(y), target.bins = 12, minpts=2)$binhi})
	filename <- paste(experiment, "bins12_table.txt", sep=".")
	write.table(bins12, file=filename , quote=F, sep="\t", row.names=TRUE, col.names=NA)
	filename <- paste(experiment, "bins12_melt_table.txt", sep=".")
	write.table(melt(bins12), file=filename , quote=F, sep="\t", row.names=TRUE, col.names=NA)
	head(melt(bins12))

	# Estimate the linear range of the Count values; second last / second value
	# old - bin12_fold <- c()
	# old - for (i in bins12) {	  bin12_fold = c(bin12_fold, round(i[11]/i[2], digits=1))  }

	bin12_fold <- apply(bins12,2,function(y){ round(y[11]/y[2], digits=1) })  
	filename <- paste(experiment, "bins12_fold.txt", sep=".")
	write.table(bin12_fold, file=filename , quote=F, sep="\t", row.names=TRUE, col.names=NA)

	p <- ggplot(melt(bins12), aes(Var2, value)) + 
	      ggtitle('Counts Distribution over 12 Bins') + 
    	  geom_line(aes(colour=Var2), size=1.2, show.legend = FALSE) + 
    	  geom_point(aes(colour=Var2), size = 10,shape=3, show.legend = FALSE) +
    	  xlab("Samples") +  ylab("Counts") + 
	      ylim(0, NA) +
    	  #scale_x_continuous(limits=c(1, 12), breaks = scales::pretty_breaks(n = 12)) +
    	  theme(
    	    plot.title = element_text(color="#993333", size=24, face="bold.italic"),
    	    axis.title.x = element_text(color="#993333", size=18, vjust = 1.5),
    	    axis.title.y = element_text(color="#993333", size=20),
    	    axis.text.x = element_text(angle = 15, hjust = .9, size=16),
    	    axis.text.y = element_text( hjust = .9, size=16)
    	  )
	# save the graphics
	filename <- paste(experiment, "bins12.png", sep=".")
	png(file = filename, width = 10, height = 10, units = "in", res=300, pointsize = 12, bg = "white")	
	print(p)
	dev.off()


	logbook <- write_log(logbook, "----> Ranking Count data ")
	# Make a ranked dataset of the Count data 
	rank_table <- x
	for (i in 1:ncol(rank_table)) {
	  rank_table <- rank_table[order(rank_table[,i]),] 
	  rank_table[,i] <- seq(1,nrow(x))
	}
	rank_table <- rank_table[order(rank_table[,1]),] 
	head(rank_table)
	filename <- paste(experiment, "rank_table.txt", sep=".")
	write.table(rank_table, file = filename, quote=F, sep="\t", col.names=NA)
	
	if (option.heatmaps) {	# Heatmap of Rank-order Samples
  	logbook <- write_log(logbook, "# 2b. Heatmap of Rank-order Samples" )
  	heatmap_temp <-(t(apply(rank_table,1,as.numeric)))
  	if (nrow(heatmap_temp) > 3) {
  	  colnames(heatmap_temp) <- colnames(rank_table)
  	  Anne_heatmap(log(heatmap_temp+1,2), experiment, "Heatmap Rank-order Samples", "Heatmap_Rank_Samples", "signal")
  	}	
	}
	
	
	# Boxplot of Count values
	logbook <- write_log(logbook, "----> Boxplot of Count values ")
	x.melt <- melt(x)
	colnames(x.melt) =c("Experiment","Value")
	x.melt$log2Value <- log(x.melt$Value,2)
	
	my_barplot_palette <- colorRampPalette(c("#B4AF91",  "#C03000"))(n = length(x))
	filename <- paste(experiment, "colored_boxplot","png", sep=".")
	png(file = filename, width = x_1k, height = y_1k, units = "px", pointsize = 12, bg = "white")	
	p <- ggplot(x.melt,  aes(Experiment,log2Value) ) +	geom_boxplot( aes(fill=Experiment)) + ggtitle("Boxplot of Expression Values of Experiments") + theme_g2d_barplot +	theme(legend.position="none") +	scale_fill_manual(values = my_barplot_palette)
	print(p)
	dev.off()
	
	
	
######################## 2. Load Factors file ##################################################
	
	unlink(errorfile)
	logbook <- write_log(logbook, "# 2. Load Factors file")
	MF <- read.delim(factors_file,sep="\t")
	MF
	
	# generic import of undefined number of factors
	MF.dim <- dim(MF)[2] 
	MF.factors <- MF[,2] # first factor
	# add more factors if there are any (MF.dim >2)
	if (MF.dim >2) {
	  for (i in 3:MF.dim) { MF.factors <- paste(MF.factors,MF[,i],sep=".") }
	}
	MF.group <- factor(MF.factors)
	MF.group.all <- cbind(MF,Group=MF.group)
	MF.group.levels <- as.list(levels(MF.group.all$Group))
	MF.group.all	
	

	
	
######################## 3. Load Class file ##################################################
	
	if (option.classes) {
		unlink(errorfile)
		logbook <- write_log(logbook, "# 3. Load Class file")
		my_class_genes <- read.delim(class_file,sep="\t")
		colnames(my_class_genes) =c("GeneID","color","Group")
	}
	
	
	
##############################################################################################
##                        Global Analysis                                                   ##
##############################################################################################

	
	unlink(errorfile)
	logbook <- write_log(logbook, "# 4. Global Analysis")
	
	
	
######################## Merge Factors and Counts #########################################

	logbook <- write_log(logbook, "	- Merge Factors and Counts and TMM normalization")

	# Merge Factors and Counts
	y <- DGEList(counts=as.matrix(x),group=MF.group)
	y <- calcNormFactors(y)
	#  y$samples
  
   # Create the design matrix
	designMF <- model.matrix(~0+MF.group,  data=y$samples)
	colnames(designMF) <- levels(MF.group)
	designMF


######################## Library Size ##################################################

	logbook <- write_log(logbook, "	- Library Size Boxplot")

	# prepare data for Boxplot
	my_samples <- y$samples
	my_samples$experiment <- row.names(my_samples)
	my_samples$lib.size <- as.numeric(my_samples$lib.size)
	my_samples$experiment <-factor(my_samples$experiment, levels=my_samples[order(my_samples$experiment), "experiment"])
	my_barplot_palette <- colorRampPalette(c("#B4AF91", "#787746", "#40411E", "#32331D", "#C03000"))(n = dim(my_samples)[1])
	bar_chart <- ggplot(my_samples, aes(experiment, lib.size, fill=factor(lib.size))) +
	  ggtitle("Library Sizes") + 
	  theme_g2d_barplot + 
	  geom_bar(stat="identity") +
	  xlab("Experiments") +  
	  ylab("Number of reads") +
	  scale_fill_manual(values = my_barplot_palette) +
	  theme_g2d_barplot +
	  theme(legend.position="none") +
	  scale_y_continuous(expand = c(0,0))
	filename <- paste(experiment, "Library_size","png", sep=".")
	png(file = filename, width = x_720p, height = y_720p, units = "px", pointsize = 12, bg = "white")
	print(bar_chart)
	dev.off()
	filename <- paste(experiment, "Library_size.txt", sep=".")
	write.table(my_samples, file=filename , quote=F, sep="\t", col.names=NA)		
	

	logbook <- write_log(logbook, "	- Library Size Boxplot raw data")
	# Boxplot of Raw Count data
	filename <- paste(experiment, "Box_plot_raw","png", sep=".")
	png(file = filename, width = x_4k, height = y_4k, units = "px", pointsize = 12, bg = "white")
	boxplot(as.data.frame(log(y$counts,2)), main="Raw Count Data", ylab="log2(counts)",  notch=FALSE)
	dev.off()
	

######################## 5. Calculation of Dispersion ###############################################
	

	# 5. Calculation of Dispersion
	unlink(errorfile)
	logbook <- write_log(logbook, "# 5. Calculation of Dispersion")
	if (!option.replicates) {
	  y <- estimateGLMCommonDisp(y,designMF,method="deviance", robust=TRUE, subset=NULL)  # no replicates
	} else {  
	  y <- estimateGLMCommonDisp(y,designMF)
	  y <- estimateGLMTrendedDisp(y,designMF, method="bin.loess")
	  y <- estimateGLMTagwiseDisp(y,designMF)
	} 
	
	


	
######################## 6. Mean values of replicates ######################################

	unlink(errorfile)
	logbook <- write_log(logbook, "# 6. Mean values of replicates")

	# Take mean of replicates if are replicates 
	y_groupheader <- y$counts
	colnames(y_groupheader)<-MF.group 
	mean_table <- c()
	for (i in MF.group.levels) {
		print(i)
		my_cols<-which(colnames(y_groupheader) %in% i)
		print(my_cols)
		if (length(my_cols)>1) {
			my_mean<-rowMeans(y_groupheader[,my_cols])
		} else {
			my_mean<-(y_groupheader[,my_cols])
		}	
		mean_table<-cbind(mean_table,my_mean)
	}

	colnames(mean_table)<-MF.group.levels
	head(mean_table)
	filename <- paste(experiment, "Mean_signals_of_replicates.txt", sep=".")
	write.table(mean_table, file=filename , quote=F, sep="\t", col.names=NA)

	# Heatmaps of Mean Signal data
	logbook <- write_log(logbook, "---> Mean Signal data Heatmap" )
	heatmap_temp <-(t(apply(mean_table,1,as.numeric)))
	if (nrow(heatmap_temp) > 3) {
		colnames(heatmap_temp) <- colnames(mean_table)
		head(heatmap_temp)
		Anne_heatmap(log(heatmap_temp+1,2), experiment, "Heatmap of Mean Signal log2(Signal)", "Heatmap_Mean_Signals", "signal")
	}
	
# Make rank-order table of Mean Signals		
	logbook <- write_log(logbook, "---> Mean Signal Rank Table and Heatmap" )
	rank_mean_table <- mean_table
	for (i in 1:ncol(rank_mean_table)) {
		rank_mean_table <- rank_mean_table[order(rank_mean_table[,i]),] 
		rank_mean_table[,i] <- seq(1,nrow(x))
	}
	rank_mean_table <- rank_mean_table[order(row.names(rank_mean_table)),] 
	head(rank_mean_table)
	filename <- paste(experiment, "rank_mean_table","txt", sep=".")
	write.table(rank_mean_table, file = filename, quote=F, sep="\t", col.names=NA)
	# Heatmap of Rank Mean Signal
	heatmap_temp <-(t(apply(rank_mean_table,1,as.numeric)))
	if (nrow(heatmap_temp) > 3) {
		colnames(heatmap_temp) <- colnames(rank_mean_table)
		head(heatmap_temp)
		Anne_heatmap(heatmap_temp, experiment, "Heatmap Rank-order Factors (1)", "Heatmap_Rank_Factors_1", "ratio")
		Anne_heatmap(heatmap_temp, experiment, "Heatmap Rank-order Factors (2)", "Heatmap_Rank_Factors_2", "signal")
	}

	
# Mean Signals Normalized per gene for heatmaps  (x-min)/max
	logbook <- write_log(logbook, "---> Mean Signals Normalized per gene for heatmaps" )
	mean_table_log <- log(mean_table,2)
	head(mean_table_log)
	filename <- paste(experiment, "Mean_signals_log2.txt", sep=".")
	write.table(mean_table_log, file=filename , quote=F, sep="\t", col.names=NA)	
	
	mean_table_log_scaled <- round((100 * mean_table_log / max(mean_table_log)), digits=3)
	head(mean_table_log_scaled)
	filename <- paste(experiment, "Mean_signals_log2_scaled.txt", sep=".")
	write.table(mean_table_log_scaled, file=filename , quote=F, sep="\t", col.names=NA)
	
	mean_table_log_scaled_per_gene <- round(t(apply(mean_table_log, 1, function(x) (x-min(x))/(max(x)-min(x)))), digits=4)
	filename <- paste(experiment, "Mean_signals_log2_scaled_per_gene.txt", sep=".")
	write.table(mean_table_log_scaled_per_gene, file=filename , quote=F, sep="\t", col.names=NA)
	

	
##############################################################################################
##                        PCA plot and Correlation matrix of Experiments and Factors                               ##
##############################################################################################
	

	unlink(errorfile)
	logbook <- write_log(logbook, "# 7. PCA plot and Correlation matrix of Experiments and Factors")
	
		
	logbook <- write_log(logbook, "----> PCA plot of Experiments and Factors")
	# PCA on experiments
	pca_data <- t(log(y$counts+1)) 
	Anne_PCA_plot(pca_data, "PCA_experiments", "PCA Plot of Experiments", "false")
	
	# PCA on experiments mean signals based on the Factors
	pca_data <- t(log(mean_table+1)) 
	Anne_PCA_plot(pca_data, "PCA_factors", "PCA Plot of Factors", "false")
	
	
	logbook <- write_log(logbook, "----> Make correlation Matrix of experiment")
	logFC <- predFC(y,designMF,prior.count=1,dispersion=0.05)
	corLogFC <- cor(logFC)
	filename <- paste(experiment, "Correlation_matrix_experiments.txt", sep=".")
	write.table(corLogFC, file=filename , quote=F, sep="\t", row.names=F)	
	
	Theme_corrMatrix <- my_theme <- theme(title=element_text(size=18, hjust=0.5),
	        axis.title = element_text(size=14),
	        axis.text  = element_text(size=14))
	p <- ggcorrplot(corLogFC, hc.order = TRUE, outline.col = "white",lab_size = 10 ) +
	  labs(title = 'Correlation Matrix of Experiments') +
	  scale_fill_gradient2(limit = c(min(corLogFC),1),  mid = "white", midpoint = (min(corLogFC))) +
	  Theme_corrMatrix
	filename <- paste(experiment, "Correlation_matrix_experiments.png", sep=".")
	# old -png(file = filename, width = x_1k, height = y_1k, units = "px", pointsize = 12, bg = "white")
	png(file = filename, width = 10, height = 10, units="in", res=300, pointsize = 24, bg = "white")
	print(p)
	dev.off()		
	
	

	
	
##############################################################################################
##                        Analysing Class data (optional)                                   ##
##############################################################################################
	
	if (option.classes) {  # A classification file is optional
    	unlink(errorfile)
    	logbook <- write_log(logbook, "# 8. Analysing Class data")
    
    # Add class (color and groups) to the mean table
    	mean_table_log_scaled_class <- mean_table_log_scaled
    	mean_table_log_scaled_class <- cbind(mean_table_log_scaled_class, row.names(mean_table_log_scaled_class))
    	colnames(mean_table_log_scaled_class) <- c(colnames(mean_table_log_scaled), "GeneID")
    	mean_table_log_scaled_class <- merge(mean_table_log_scaled_class, my_class_genes, by="GeneID",all.x=T)
    	mean_table_log_scaled_class$color <- ifelse(is.na(mean_table_log_scaled_class$color), "grey", as.character(mean_table_log_scaled_class$color))	# add a default color if the color is NA
    	mean_table_log_scaled_class$Group <- ifelse(is.na(mean_table_log_scaled_class$Group), "NA", as.character(mean_table_log_scaled_class$Group))		# add a default Group if the Group is NA
    	head(mean_table_log_scaled_class)
    	filename <- paste(experiment, "Mean_signals_log2_scaled_CLASSES.txt", sep=".")
    	write.table(mean_table_log_scaled_class, file=filename , quote=F, sep="\t", row.names=FALSE)
    
    
    ######################## Line Charts of Classes ######################################
    
    # Line Chart for each subClass on the basis of the group column
    	mean_table_class <- subset(mean_table, row.names(mean_table) %in% my_class_genes$GeneID)
    	mean_table_class_log <- subset(mean_table_log, row.names(mean_table_log) %in% my_class_genes$GeneID)
    	head(mean_table_class)
    	if (nrow(mean_table_class)<5) {
    		logbook <- write_log(logbook, "# ERROR Skipping Class signal plots due to low numbers classes found")
    	} else {
    		#draw a graph for each class color
    		logbook <- write_log(logbook, "-   a) Create Graphics for mean signals of each subClass")
    		Class_signal_plots(mean_table_class, "Mean_Signals_subClass")
    		logbook <- write_log(logbook, "-   b) Create Graphics for mean log2(signal) of each subClass")
    		Class_signal_plots(mean_table_class_log, "Mean_Signals_subClass_log2")
    	}
    
    ######################## Correlation matrix of Classes ######################################
    
    # Correlation matrix of All Class genes, using mean signal of replicates
    	unlink(errorfile)
    	logbook <- write_log(logbook, "----> Correlation matrix of All Class genes, using mean signal of replicates")
    	cor_mean_table_class <- cor(t(mean_table_class))	
    	filename <- paste(experiment, "Correlation_matrix_Class_genes.txt", sep=".")
    	write.table(cor_mean_table_class, file=filename , quote=F, sep="\t", row.names=F)
    
    	filename <- paste(experiment, "Correlation_matrix_Class_genes.png", sep=".")
    	png(file = filename, width = x_2k, height = y_2k, units = "px", pointsize = 12, bg = "white")
    	m1 <- melt(cor_mean_table_class, id = row.names)
    	#qplot(x=X1, y=X2, data=as.data.frame(m1),main = "Correlation matrix of Class genes", xlab="Class genes", ylab="Class genes", fill=value, geom="tile", zlim= c(0,1)) + theme_g2d_plot_cor + scale_fill_gradient2(limits=c(-1, 1))
    	qplot(x=X1, y=X2, data=as.data.frame(m1),main = "Correlation matrix of Class genes", xlab="Class genes", ylab="Class genes", fill=value, geom="tile") + theme_g2d_plot_cor + scale_fill_gradient2(limits=c(-1, 1))
    	dev.off()
    
    # Correlation matrix of each subClass, using mean signal of replicates
    	unlink(errorfile)
    	logbook <- write_log(logbook, "----> Correlation matrix of subClass genes, using mean signal of replicates")
    	SubClass_correlation_plots(cor_mean_table_class, "SubClass_correlation_matrix" )
  }
	
	

	
##############################################################################################
##                        Contrasts analysis                                                ##
##############################################################################################
	
	
  ######################## Load the Contrasts ######################################
	unlink(errorfile)
	logbook <- write_log(logbook, "# 9. Load the Contrasts")
	# contrast from file
	my_contrasts <- readLines(contrast_file)
	my_contrasts
	contrast.matrix <- makeContrasts(contrasts=my_contrasts, levels=designMF)
	contrast.matrix
	contrast.matrix.dim <- dim(contrast.matrix)[2]


  ######################## Fit the data ######################################
	unlink(errorfile)
	logbook <- write_log(logbook, "# 10. Fit the data")
	
	if (option.replicates) { # Data contains replicates
	  logbook <- write_log(logbook, "----> Replicates found using calculated dispersion")
	  fit.design <- glmFit(y, designMF)
	  lrt.design <- glmLRT(fit.design, designMF)
	  filename <- paste(experiment, "Normalised_data_table.txt", sep=".")
	  write.table(lrt.design$fitted.values, file=filename , quote=F, sep="\t", col.names=NA)
	} else {  # NO replicates
	  logbook <- write_log(logbook, "----> No Replicates, estimated dispersion will be used ")
	  bcv <- 0.05      	
	  DGE_TableList <- c()
	  for (i in 1:length(my_contrasts)) {
	    my_list <- as.list(strsplit(my_contrasts[i], "-")[[1]])
	    #OLD: contrast_columns <- strtoi( rownames(MF.group.all[MF.group.all$Group %in% my_list, ]) ) # Get the rownames (in this case the numbers) of the MF.group.all$Group 
	    #OLD: DGE_TableList[[i]] <- exactTest(y[,contrast_columns], dispersion=bcv^2) 
	    DGE_TableList[[i]] <- exactTest(y, pair=rev(my_list), dispersion=bcv^2) # new
		# change PValue to pvalue
	    colnames(DGE_TableList[[i]]$table)[colnames(DGE_TableList[[i]]$table)=="PValue"] <- "pvalue"
	    # add the BH adjusted p-value
	    DGE_TableList[[i]]$table$adj_pvalue <- p.adjust(DGE_TableList[[i]]$table$pvalue, method = "BH", n = length(DGE_TableList[[i]]$table$pvalue))
	    # add the Fold Column
	    DGE_TableList[[i]]$table$Fold <- round( sign(DGE_TableList[[i]]$table$logFC) * abs(DGE_TableList[[i]]$table$logFC)^2 , 3)
	    # add the GeneID column from rownames
	    DGE_TableList[[i]]$table$GeneID <- rownames(DGE_TableList[[i]]$table)
	    # reorder the columns
	    col_order <- c("GeneID","Fold","logFC","logCPM","pvalue","adj_pvalue")
	    DGE_TableList[[i]]$table <- DGE_TableList[[i]]$table[, col_order]
	    
	    filename <- paste(experiment, "Contrast",my_contrasts[i],"txt", sep=".")
	    write.table(DGE_TableList[[i]]$table, file=filename , quote=F, sep="\t", col.names=NA)  
	  }
	  # show first contrast      	
	  head(DGE_TableList[[1]]$table)  
  } 

	
	
	
	
	
	
##############################################################################################
##                        Analyse and plot each Contrast                                    ##
##############################################################################################
	

	unlink(errorfile)
	logbook <- write_log(logbook, "# 11. Analyse and plot each Contrast")
	
		
  # Loop through the Contrasts and draw for each Contrast a: 
  		## - DE (Differential Expression ) table
  		## - Significant Fold Change Plots
  		## - MA plot
  		## - Volcano plot
  		## - Fold change distribution plot
  		## - Count distribution plot
	tophits <- c()
	tophits_highfold <-c()
	tophits_no_background <-c()
	filenames.MAplots <- list()
	filenames.VolcanoPlots <- list()
	filenames.SigChanged_plots <- list()
	filenames.Significant_Changed_Genes_tables <- list()
	my_overview <- c()
	my_GeneNetwork <- c()
	volcano_tables <- NULL
  i <- 1
	for (i in 1:contrast.matrix.dim) {
		logbook <- write_log(logbook, my_contrasts[i] )
		# A. DE (Differential Expression ) table
  		# Replicates
  		if (option.replicates) {
  		  lrt <- glmLRT(fit.design , contrast=contrast.matrix[,i])
  		  my.DGEList.glmFit = topTags(lrt, n=100000, adjust.method="BH")
  		  DGE_Table = cbind(rownames(my.DGEList.glmFit$table), my.DGEList.glmFit$table)
  		  # add fold change
  		  DGE_Table <- cbind(DGE_Table, sign(DGE_Table$logFC)*(2^abs(DGE_Table$logFC)))
  		  colnames(DGE_Table) =c("GeneID","logFC", "logCPM", "LR", "pvalue", "adj_pvalue", "Fold")
  		  
  		} else { # no replicates
  		  DGE_Table = DGE_TableList[[i]]$table    # get it from the ExactTest with estimated dispersion
  		  
  		}  		
  		
  		DGE_Table$minFDR <- -log(DGE_Table$adj_pvalue,2)
  		DGE_Table$minFDR[as.numeric(DGE_Table$minFDR) > 200] <- 200    # use 100 for cutoff
  		filename <- paste(experiment, "Differential_Expression", i, my_contrasts[i] ,"txt", sep=".")
  		write.table(DGE_Table, file=filename , quote=F, sep="\t", row.names=T)
  		
  	# B. For the html table only export the up and down regulated genes
  		tophits_ONLY <- DGE_Table[( abs(DGE_Table$Fold)>=2 & DGE_Table$adj_pvalue<0.05),]
  		tophits_ONLY <- tophits_ONLY[with(tophits_ONLY, order(-Fold)), ]
  		tophits_ONLY$logFC <- sprintf("%.2f", tophits_ONLY$logFC)
  		tophits_ONLY$logCPM <- sprintf("%.2f", tophits_ONLY$logCPM)
  		tophits_ONLY$minFDR <- sprintf("%.2f", tophits_ONLY$minFDR)
  		tophits_ONLY$minFDR[as.numeric(tophits_ONLY$minFDR) > 100] <- 100    # use 100 for cutoff
  		tophits_ONLY$pvalue <- sprintf("%.1e", tophits_ONLY$pvalue) 
  		tophits_ONLY$adj_pvalue <- sprintf("%.1e", tophits_ONLY$adj_pvalue) 
  		tophits_ONLY$Fold <- sprintf("%.1f", tophits_ONLY$Fold)
  		
  		filename <- paste(experiment, "Significant_Changed_Genes", i, my_contrasts[i] ,"txt", sep=".")
  		filenames.Significant_Changed_Genes_tables[i] <- filename	
  		write.table(tophits_ONLY, file=filename , quote=F, sep="\t", row.names=F)	
  		
  	# C. GeneNetwork; Add tophits for each contrast to network compatible table
  		if (nrow(tophits_ONLY) > 0) {
  		  my_GeneNetwork <- rbind(my_GeneNetwork, data.frame(Contrasts=my_contrasts[i], GeneID=tophits_ONLY$GeneID, logFC=tophits_ONLY$logFC, adj_pvalue=tophits_ONLY$adj_pvalue))
  		}
  		
  	# E. Add data to overview table
  		total<- length(DGE_Table$GeneID)
  		up   <- length(DGE_Table[( DGE_Table$Fold >= 2 & DGE_Table$adj_pvalue <= 0.05),]$GeneID)
  		down <- length(DGE_Table[( DGE_Table$Fold <= -2 & DGE_Table$adj_pvalue <= 0.05),]$GeneID)
  		my_overview <- rbind(my_overview, data.frame(Contrasts=my_contrasts[i], Total=total, Up=up, Down=down, TableFilename=filename))	
  		
  		
  	# D. Plot tophits_ONLY for all genes
  		filename <- paste(experiment, "Significant_Changed_Genes", i, my_contrasts[i] ,"png", sep=".")
  		filenames.SigChanged_plots[i] <- filename		
  		png(file = filename, width = x_2k, height = y_2k, units = "px", pointsize = 12, bg = "white")
  		my_title=paste(experiment, "Significant changed genes", my_contrasts[i], sep=" ")
  		tophits_ONLY$logFC <- as.numeric(tophits_ONLY$logFC)
  		tophits_ONLY$GeneID <-factor(tophits_ONLY$GeneID, levels=tophits_ONLY[order(-tophits_ONLY$logFC), "GeneID"])
  		tophits_ONLY$Up_Down <- ifelse(tophits_ONLY$logFC>0, "UP", "DOWN")
  		bar_chart <- ggplot(tophits_ONLY, aes(GeneID, logFC, fill=factor(Up_Down), colour = factor(Up_Down)), main = my_title ) + 
  		  ggtitle(my_title) + 
  		  theme_g2d_barplot + 
  		  geom_bar(stat="identity") +  
  		  xlab("Genes") +  ylab("Fold Change (Log2FC)") + 
  		  scale_colour_manual(breaks = tophits_ONLY$updown, values = c("#FFFAF0","#5494FF")) + 
  		  scale_fill_manual(breaks = tophits_ONLY$updown, values = c("#FFD312","#0608B2"))
  		print(bar_chart)
  		dev.off()
  		
  	# E. Plot All genes of each class
  		# to do
  		
  		
  	# F. add the GeneIDs to the Tophits list for all experiments
  		tophits <- c(tophits, as.character(DGE_Table[( abs(DGE_Table$Fold)>2 & DGE_Table$adj_pvalue<0.05),]$GeneID))
  		tophits_highfold <- c(tophits_highfold, as.character(DGE_Table[( abs(DGE_Table$Fold)>5 & DGE_Table$adj_pvalue<0.01),]$GeneID))
  		tophits_no_background <- c(tophits_no_background, as.character(DGE_Table[( abs(DGE_Table$Fold)>1.4 & DGE_Table$adj_pvalue<0.20),]$GeneID))
  		
  		
  		
  	# G. Add all individual files to one file for TMEV or other downstream programs	
  		my_columns <- cbind(rownames(DGE_Table), DGE_Table$logFC)
  		colnames(my_columns) = c("GeneID","logFC")
  		if (i == 1 ) {
  		  TMEV <- my_columns
  		  TMEV.colnames <- c("GeneID",my_contrasts[i]) ;
  		}
  		if (i > 1 ) {	
  		  TMEV<-merge(TMEV, my_columns, by="GeneID",all.x=T)
  		  TMEV.colnames <- c(TMEV.colnames,my_contrasts[i]) ;
  		}	
  		colnames(TMEV) = TMEV.colnames
  		
  	# H. Select data for colouring classes 
  		if (option.classes) {
    		my_class <- DGE_Table[which(row.names(DGE_Table) %in% my_class_genes$GeneID),]
    		my_class <- merge(my_class, my_class_genes,  by="GeneID",all.x=T)
    		colnames(my_class) =c("GeneID","logFC", "logCPM", "LR", "pvalue", "adj_pvalue", "Fold", "minFDR", "color")
  		} else {
  		  # just take the top 20 and color it red
    		my_class <- DGE_Table[order(abs(DGE_Table$Fold),decreasing=T),][1:20,]
    		my_class$color <- "red"
  		}
      head(my_class)
  		
  		# Select data with 0 expression and mark with X it in the plots
  		#	my_null_class <- DGE_Table[which(row.names(DGE_Table) %in% row.names(lrt$fitted.values)[lrt$fitted.values==0]),]
  		#	my_null_class
  		
  		
  		
  	# J - Volcano Data
  		# Volcano Plot
  		
  		# add class data
  		row.names(my_class_genes) <- my_class_genes$GeneID
  		head(my_class_genes) 
  		#DGE_Table$sfere <- DGE_Table$logCPM - median(DGE_Table$logCPM)
  		DGE_Table$sfere <- ifelse(DGE_Table$logCPM<1, 1, DGE_Table$logCPM)
  		
  		
  		DGE_Table_Class <- merge(DGE_Table, my_class_genes, by="row.names",all.x=FALSE)
  		colnames(DGE_Table_Class)[colnames(DGE_Table_Class)=="GeneID.x"] <- "GeneID"
  		head(DGE_Table_Class)
  		DGE_Table_All <- merge(DGE_Table, my_class_genes, by="row.names",all.x=TRUE)
  		colnames(DGE_Table_All)[colnames(DGE_Table_All)=="GeneID.x"] <- "GeneID"
  		head(DGE_Table_All)
  		
  		filename  <- paste(experiment, "Volcano_plot", i, my_contrasts[i] ,"txt", sep=".")
  		write.table(DGE_Table_All, file=filename , quote=F, sep="\t", row.names=FALSE)
  		volcano_tables[i] <- list(DGE_Table_All) 	# outside the loop we want to filter on tophits_no_background, so we store this table
  		
  		
  		
  		
  		
  	#  Volcano Plot for Admire contain all data
  		xmin <- min(DGE_Table_All$logFC, na.rm=T)
  		xmax <- max(DGE_Table_All$logFC, na.rm=T)
  		ymin <- min(DGE_Table_All$minFDR, na.rm=T)
  		ymax <- max(DGE_Table_All$minFDR, na.rm=T)
  		FDR_threshold <- -log2(0.05)
  		Fold_threshold <- log2(2)
  		
  		
 		filename  <- paste(experiment, "Volcano_plot", i, my_contrasts[i] ,"png", sep=".")
  		filenames.VolcanoPlots[i] <- filename
  		png(file = filename, width = x_1k, height = y_1k,  units = "px", pointsize = 12, bg = "white")
		# Plot only significant values
  		DGE_Table_All_Signifcant <- subset(DGE_Table_All, abs(logFC)>Fold_threshold & minFDR>0.01) 
  		p <- ggplot(DGE_Table_All,  aes(x=logFC,y=minFDR, color=factor(Group))) + 
  		  geom_point(size = DGE_Table_All$sfere, alpha = 1/4) +
  		  annotate("rect", xmin = xmin, xmax = xmax, ymin = 0, ymax = FDR_threshold, alpha = .2) +
  		  annotate("rect", xmin = -Fold_threshold, xmax = Fold_threshold, ymin = FDR_threshold, ymax = ymax, alpha = .2) +
  		  ggtitle(paste(experiment, "\nVolcano Plot ", my_contrasts[i] , sep=" ")) +
  		  theme(plot.title = element_text(hjust=0.5))
  		print(p)
  		dev.off()
  


  	# I. - MA plots
  		xmin <- min(DGE_Table_All$logCPM, na.rm=T)
  		xmax <- max(DGE_Table_All$logCPM, na.rm=T)
  		Fold_threshold <- log2(2)
  		
  		filename <- paste(experiment, "MAplot", i, my_contrasts[i] ,"png", sep=".")
  		filenames.MAplots[i] <- filename
  		plottitle <- paste("MAplot", experiment, "contrast=", my_contrasts[i] , sep=" ")
  		png(file = filename, width = x_1k, height = y_1k,  units = "px", pointsize = 12, bg = "white")
  		p<- ggplot(DGE_Table_All,  aes(x=logCPM,y=logFC, color=factor(Group))) + geom_point() +
  		  geom_text(data=subset(DGE_Table_All, abs(logFC) > 2 & minFDR > 250), aes(label=GeneID),hjust=0, vjust=0) +
  		  annotate("rect", xmin = xmin, xmax = xmax, ymin = -Fold_threshold, ymax = Fold_threshold, alpha = .2) +
  		  ggtitle(paste(experiment, "\nMA Plot of ", my_contrasts[i] , sep=" ")) +
  		  theme(plot.title = element_text(hjust=0.5))          
  		print(p)
  		dev.off()   
  		
  		if (option.classes) {  # MA plot of Class genes
    		filename <- paste(experiment, "MAplot_Class_Genes", i, my_contrasts[i] ,"png", sep=".")
    		plottitle <- paste("MAplot", experiment, "contrast=", my_contrasts[i] , sep=" ")
    		png(file = filename, width = x_1k, height = y_1k,  units = "px", pointsize = 12, bg = "white")
    		p<- ggplot(DGE_Table_Class,  aes(x=logCPM,y=logFC, color=factor(Group))) + geom_point() +
    		  geom_text(data=subset(DGE_Table_Class, abs(logFC) > 2 & minFDR > 250), aes(label=GeneID),hjust=0, vjust=0) +
    		  annotate("rect", xmin = xmin, xmax = xmax, ymin = -Fold_threshold, ymax = Fold_threshold, alpha = .2) +
    		  ggtitle(paste(experiment, "\nMA Plot of Class Genes of ", my_contrasts[i] , sep=" ")) +
    		  theme(plot.title = element_text(hjust=0.5))          
    		print(p)
    		dev.off() 
  		}
  		
	
		
	}
	
  
  
  
  
  
  
##############################################################################################
##                        TMEV: The Multiple Experiment Values                              ##
##############################################################################################
  
  
  unlink(errorfile)
  logbook <- write_log(logbook, "# TMEV: The Multiple Experiment Values; Merge of all contrasts DGE_Tables" )
  

  # 1. All values
  colnames(TMEV) = TMEV.colnames
  TMEV <- as.data.frame(TMEV)
  filename <- paste(experiment, "TMEV.txt", sep=".")
  write.table(TMEV, file=filename , quote=F, sep="\t", row.names=F)
  row.names(TMEV) <- TMEV$GeneID
  head(TMEV)
  
  # 2. Tophits only
  TMEV_TOPHITS <- as.data.frame(TMEV[which(row.names(TMEV) %in% tophits),])
  TMEV_TOPHITS <- as.data.frame(TMEV_TOPHITS[ order(TMEV_TOPHITS[,1],decreasing = TRUE), ])
  TMEV_TOPHITS[1] <- NULL  # remove the GeneID column which is now the row.names
  filenameTopHits <- paste(experiment, "TMEV_TOPHITS_logFC.txt", sep=".")
  write.table(TMEV_TOPHITS, file=filenameTopHits , quote=F, sep="\t", row.names=T, col.names=NA)	
  head(TMEV_TOPHITS)
  
  # 3. High Fold only 
  TMEV_TOPHITS_highfold <- as.data.frame(TMEV[which(row.names(TMEV) %in% tophits_highfold),])
  TMEV_TOPHITS_highfold <- TMEV_TOPHITS_highfold[ order(TMEV_TOPHITS_highfold[,1],decreasing = TRUE), ]
  TMEV_TOPHITS_highfold[1] <- NULL  # remove the GeneID column which is now the row.names
  head(TMEV_TOPHITS_highfold)
  filenameTopHitsHighFold <- paste(experiment, "TMEV_TOPHITS_logFC_highfold.txt", sep=".")
  write.table(TMEV_TOPHITS_highfold, file=filenameTopHitsHighFold , quote=F, sep="\t", row.names=T, col.names=NA)
  
  
  # 4. Class Genes
  if (option.classes) {
    # Tophits only
    CLASS_TOPHITS <- subset(TMEV_TOPHITS, row.names(TMEV_TOPHITS) %in% my_class_genes$GeneID)
    filename <- paste(experiment, "CLASS_TOPHITS_logFC.txt", sep=".")
    write.table(CLASS_TOPHITS, file=filename , quote=F, sep="\t", row.names=T, col.names=NA)	
    
    # High Fold only 
    CLASS_TOPHITS_highfold <- subset(TMEV_TOPHITS_highfold, row.names(TMEV_TOPHITS) %in% my_class_genes$GeneID)
    filename <- paste(experiment, "CLASS_TOPHITS_highfold_logFC.txt", sep=".")
    write.table(CLASS_TOPHITS_highfold, file=filename , quote=F, sep="\t", row.names=T, col.names=NA)	
    
    # Counts data
    CLASS_TOPHITS_signal <- subset(mean_table_class, row.names(mean_table_class) %in% tophits)
    filename <- paste(experiment, "CLASS_Signal.txt", sep=".")
    write.table(CLASS_TOPHITS_highfold, file=filename , quote=F, sep="\t", row.names=T, col.names=NA)	
    
  }
  


##############################################################################################
##                        Store Filenames of Analyse and plot each Contrast                 ##
##############################################################################################
  
  
  
	# save filenames
	filename <- paste(experiment, "MAplots.filenames.txt", sep=".")
	write.table(unlist(filenames.MAplots), file=filename, row.names = FALSE, col.names = FALSE, quote=FALSE)
	
	filename <- paste(experiment, "VolcanoPlots.filenames.txt", sep=".")
	write.table(unlist(filenames.VolcanoPlots), file=filename, row.names = FALSE, col.names = FALSE, quote=FALSE)
	
	filename <- paste(experiment, "SigChanged_plots.filenames.txt", sep=".")
	write.table(unlist(filenames.SigChanged_plots), file=filename, row.names = FALSE, col.names = FALSE, quote=FALSE)
	
	filename <- paste(experiment, "Significant_Changed_Genes_tables.filenames.txt", sep=".")
	write.table(unlist(filenames.Significant_Changed_Genes_tables), file=filename, row.names = FALSE, col.names = FALSE, quote=FALSE)

	# save my_overview table
	filename <- paste(experiment, "Contrasts_overview.txt", sep=".")
	write.table(my_overview, file=filename , quote=F, sep="\t", row.names=FALSE)

	# Store data for volcano plots where at least one gene is above the background
	logbook <- write_log(logbook, "# 16.5 Write filtered Volcano plots data" )
	

	
	
##############################################################################################
##                        Admire data files                                                 ##
##############################################################################################
	

	unlink(errorfile)
	logbook <- write_log(logbook, "Creating data for Admire" )
	

	######################## Admire Volcano Plots ######################################
	json.menu <- NULL
	json.menu$name = unbox("Admire Volcano plots")	
	# Create a items list containing the plots 
	json.items <- list()
	head(as.data.frame(volcano_tables))
	for (i in 1:contrast.matrix.dim) {
		volcano_filtered <- as.data.frame(volcano_tables[i])
		volcano_filtered <- volcano_filtered[which(volcano_filtered$GeneID %in% tophits_no_background),]
		filename  <- paste(experiment, "Volcano_plot_filtered", i, my_contrasts[i] ,"txt", sep=".")
		write.table(volcano_filtered, file=filename , quote=F, sep="\t", col.names=NA)
		
		# JSON: Description of each plot
		json.item <- NULL
		json.item$name <- unbox(my_contrasts[i])  
		json.item$type <- unbox("scatterplot")
		json.item$settingsFile <- unbox('item-settings/molgen-volcanoplot.json')
			json.items.settings <- NULL
			json.items.settings$filename <-  unbox(paste(websessiondir, filename, sep="/"))
		json.item$settings <- json.items.settings
		# add the description to the items list
		json.items[[i]] <- json.item
	}	
	# once we created all the items add these to the menu
	json.menu$items <- json.items
	# and finally we add the menu to the menus list
	json.menus[[1]] <- json.menu 
	

  	
	######################## Admire PCA Plots #########################################
	pca_data <- as.data.frame(apply(t(TMEV_TOPHITS),1,as.numeric))
	rownames(pca_data) <- rownames(TMEV_TOPHITS)
	pca_data <- pca_data[which(rownames(pca_data) %in% tophits_no_background),]				# remove background
	# Add columns for Admire
	if (!is.null(nrow(pca_data))) { 
		# A. Plot a static PCA and data table
		Anne_PCA_plot(pca_data, "PCA_genes", "PCA Plot of Genes", "true")
		# B. Interactive Plot for Admire
		genes.pca <- prcomp(pca_data, center = TRUE, scale. = FALSE) 
		admire.pca <- as.data.frame(genes.pca$x)
		admire.pca$GeneID <- rownames(admire.pca)
		# take the color groups and sfere from a volcanoplot
		admire.pca <- merge(admire.pca, volcano_tables[1], by="GeneID",all.x=T)
		head(admire.pca)
		filename  <- paste(experiment, "PCA_genes" ,"txt", sep=".")
		write.table(admire.pca, file=filename , quote=F, sep="\t", col.names=NA)
		# add info to jSON file
		json.menu <- NULL
		json.menu$name = unbox("Admire PCA plots")	
		# Create a items list containing the plots 
		json.items <- list()
		# JSON: Description of PCA plot PC1 x PC2
			json.item <- NULL
			json.item$name <- unbox("PCA PC1xPC2")  
			json.item$type <- unbox("scatterplot")
			json.item$settingsFile <- unbox('item-settings/molgen-pca.json')
				json.items.settings <- NULL
				json.items.settings$xParameter <- unbox("PC1")
				json.items.settings$yParameter <- unbox("PC2")
				json.items.settings$xLabel <- unbox("PC1")
				json.items.settings$yLabel <- unbox("PC2")
				json.items.settings$filename <-  unbox(paste(websessiondir, filename, sep="/"))
			json.item$settings <- json.items.settings
			# add the description to the items list
			json.items[[1]] <- json.item
		# JSON: Description of PCA plot PC1 x PC3
		if("PC3" %in% colnames(admire.pca)) { 
			json.item <- NULL
			json.item$name <- unbox("PCA PC1xPC3")  
			json.item$type <- unbox("scatterplot")
			json.item$settingsFile <- unbox('item-settings/molgen-pca.json')
				json.items.settings <- NULL
				json.items.settings$xParameter <- unbox("PC1")
				json.items.settings$yParameter <- unbox("PC3")
				json.items.settings$xLabel <- unbox("PC1")
				json.items.settings$yLabel <- unbox("PC3")
				json.items.settings$filename <-  unbox(paste(websessiondir, filename, sep="/"))
			json.item$settings <- json.items.settings
			# add the description to the items list
			json.items[[2]] <- json.item
		}	
		# once we created all the items add these to the menu
		json.menu$items <- json.items
		# and finally we add the menu to the menus list
		json.menus[[2]] <- json.menu 
	}	
	

##############################################################################################
##                        Cohesion of Contrasts                                                 ##
##############################################################################################
#
#

	filename <- paste(experiment, "GeneNetwork_EdgeList.txt", sep=".")
	write.table(my_GeneNetwork, file=filename , quote=F, sep="\t", row.names=FALSE)	

    my_GeneNetwork_high <- my_GeneNetwork[( abs(as.numeric(as.character(my_GeneNetwork$logFC)))>=4 & as.numeric(as.character(my_GeneNetwork$adj_pvalue))<0.01 ),]
    filename <- paste(experiment, "GeneNetwork_EdgeList_high.txt", sep=".")
    write.table(my_GeneNetwork_high, file=filename , quote=F, sep="\t", row.names=FALSE)	
	
##############################################################################################
##                        Network in iGraph                                                 ##
##############################################################################################
#
#	
#	# 17. Draw the network using iGraph
#	unlink(errorfile)
#	logbook <- write_log(logbook, "# 17. Draw Networks" )
#	
#	# The complete network
#	filename <- paste(experiment, "GeneNetwork_EdgeList.txt", sep=".")
#	write.table(my_GeneNetwork, file=filename , quote=F, sep="\t", row.names=FALSE)	
#	
#	my_Error <- tryCatch(
#	  Draw_GeneNetwork(my_GeneNetwork, "GeneNetwork_of_Contrasts", layout.reingold.tilford ),
#	  error=function(e) e
#	)
#	
#	
#	# The high fold change network
#	my_Error <- tryCatch(
#	  {
#	    my_GeneNetwork_high <- my_GeneNetwork[( abs(as.numeric(as.character(my_GeneNetwork$logFC)))>=4 & as.numeric(as.character(my_GeneNetwork$adj_pvalue))<0.01 ),]
#	    # frequency 
#	    freq_table <- table(as.character(my_GeneNetwork_high$GeneID))
#	    
#	    filename <- paste(experiment, "GeneNetwork_EdgeList_high.txt", sep=".")
#	    write.table(my_GeneNetwork_high, file=filename , quote=F, sep="\t", row.names=FALSE)
#	    Draw_GeneNetwork(my_GeneNetwork_high, "GeneNetwork_of_Contrasts_high", layout.reingold.tilford )
#	  },
#	  error=function(e) e
#	)
	
	
	
	
##############################################################################################
##                        Heatmaps                                                          ##
##############################################################################################
	
  if (option.heatmaps) {
	
  	######################## Heatmap settings #########################################
  	my_heatmap_palette_signal <- colorRampPalette(c("#FCFFF5", "#D1DBBD", "#91AA9D", "#3E606F", "#193441"))(n = 1000) 
  	my_heatmap_palette_ratio <- colorRampPalette(c("#4C3200", "#CC8400", "white","#00006B", "#00004C"))(n = 1000)  # orange - blue
  	my_heatmap_margins <- c(14,10)
  	my_heatmap_pointsize <- 24
  
  	unlink(errorfile)
  	logbook <- write_log(logbook, "# Heatmaps" )
	
	
  	######################## Heatmap of replicates #########################################
  
# TODO:  	logbook <- write_log(logbook, "----> Heatmap of replicates" )
# TODO:  	
# TODO:  	n <- length(row.names(lrt.design$fitted.values))
# TODO:  	normalized_signals <- as.data.frame((lrt.design$fitted.values))
# TODO:  	normalized_signals <- normalized_signals[order(normalized_signals[,1]),]
# TODO:  	normalized_signals <- normalized_signals[round(n*.1):round(n-n*.3),]	# remove high and low 20% counts
# TODO:  	# save the heatmap data to file
# TODO:  	filename <- paste(experiment, "Heatmap_EXPERIMENTS", "txt", sep=".")
# TODO:  	write.table(normalized_signals, file=filename , quote=F, sep="\t", row.names=T, col.names=NA)	
# TODO:  	if (nrow(normalized_signals) > 3) {
# TODO:  		filename <- paste(experiment, "Heatmap_EXPERIMENTS", "png", sep=".")
# TODO:  		png(file = filename, width = x_2k, height = y_2k, units = "px", pointsize = my_heatmap_pointsize, bg = "white")	
# TODO:  		par(mar=c(15,6,13,2)) # bottom, left, top, right
# TODO:  		
# TODO:  		heatmap.2(as.matrix(log2(normalized_signals+0.5)),dendrogram="both",Rowv=TRUE, col=my_heatmap_palette_signal, scale="none", key=T, keysize=1, main="Normalized log2(Signals), Experiments vs Genes",margins=my_heatmap_margins, density.info="none", trace="none",cexCol=0.9, labRow=NA)
# TODO:  		dev.off()
# TODO:  	}	
  

	
		######################## Heatmap of TopHits #########################################	
		
		logbook <- write_log(logbook, "----> Heatmap of TopHits" )
		
  	if (ncol(TMEV_TOPHITS) > 1 & nrow(TMEV_TOPHITS) < 2) {
  	  logbook <- write_log(logbook, "# Skipping heatmaps of TopHits: Only one contrasts found" )
  	 } else {
		  logbook <- write_log(logbook, "# Generating Heatmaps of TopHits" )
  		# new json menu for heatmaps
  		json.menu <- NULL
  		json.menu$name = unbox("Heatmaps of TopHits")	
  		json.items <- list()
  			
  		# Heatmap of TopHits
  		heatmap_temp <-as.matrix(apply(t(TMEV_TOPHITS),1,as.numeric))
			if (nrow(heatmap_temp) > 3) {
				rownames(heatmap_temp) <- rownames(TMEV_TOPHITS)
				colnames(heatmap_temp) <- colnames(TMEV_TOPHITS)
				Anne_heatmap(heatmap_temp, experiment, "Heatmap of TopHits Ratio (logFC)", "Heatmap_TopHits", "ratio")
				filename <- paste(experiment, "Heatmap_TopHits.sorted.txt", sep=".")
				# Init: add heatmaps item to json container
					json.item <- NULL
					json.item$name <- unbox("Heatmap of TopHits")  
					json.item$type <- unbox("heatmap")
					json.item$settingsFile <- unbox('item-settings/molgen-heatmap.json')
						json.items.settings <- NULL
						json.items.settings$xLabel <- unbox("Contrasts")
						json.items.settings$yLabel <- unbox("Genes")
						json.items.settings$filename <-  unbox(paste(websessiondir, filename, sep="/"))
					json.item$settings <- json.items.settings
					json.items[[1]] <- json.item
			}		
			
	    ######################## Heatmap of TopHits #########################################	
			
    	logbook <- write_log(logbook, "----> Heatmap of High Fold TopHits" )
    			
  		heatmap_temp <-as.matrix(apply(t(TMEV_TOPHITS_highfold),1,as.numeric))
  		if (nrow(heatmap_temp) > 3) {
  			rownames(heatmap_temp) <- rownames(TMEV_TOPHITS_highfold)
  			colnames(heatmap_temp) <- colnames(TMEV_TOPHITS_highfold)	
  			Anne_heatmap(heatmap_temp, experiment, "Heatmap of High Fold Changed TopHits", "Heatmap_TopHits_HighFold", "ratio")
  			filename <- paste(experiment, "Heatmap_TopHits_HighFold.sorted.txt", sep=".")
  				# add headtmap to json container
  				json.item <- NULL
  				json.item$name <- unbox("Heatmap of TopHits HighFold")  
  				json.item$type <- unbox("heatmap")
  				json.item$settingsFile <- unbox('item-settings/molgen-heatmap.json')
  					json.items.settings <- NULL
  					json.items.settings$xLabel <- unbox("Contrasts")
  					json.items.settings$yLabel <- unbox("Genes")
  					json.items.settings$filename <-  unbox(paste(websessiondir, filename, sep="/"))
  				json.item$settings <- json.items.settings
  				json.items[[2]] <- json.item
  
  			# Final: add heatmap item to json container 
  			json.menu$items <- json.items
  			json.menus[[3]] <- json.menu 
  		}	
		
  		# 25b. Heatmaps of Rank_Table
  		##			logbook <- write_log(logbook, "# 25b. Heatmap of Rank-order Samples" )
  		##			heatmap_temp <-(t(apply(rank_table,1,as.numeric)))
  		##			if (nrow(heatmap_temp) > 3) {
  		##				colnames(heatmap_temp) <- colnames(rank_table)
  		##				Anne_heatmap(log(heatmap_temp+1,2), experiment, "Heatmap Rank-order Samples", "Heatmap_Rank_Samples", "signal")
  		##			}
		}
  } 
	
	
	
	
	
		
##############################################################################################
##                        Heatmaps of Classes                                               ##
##############################################################################################
		
	if (option.classes & option.heatmaps) {	
	  unlink(errorfile)
		logbook <- write_log(logbook, "# Heatmaps of Classes" )

		######################## Heatmap of All Class Gens TopHits ###############################	
		logbook <- write_log(logbook, "----> Heatmap of All Class Genes TopHits" )
		# OLD: CLASS_TOPHITS <- subset(TMEV_TOPHITS, row.names(TMEV_TOPHITS) %in% my_class_genes$GeneID)
		heatmap_temp <- (t(apply(CLASS_TOPHITS,1,as.numeric)))
		if (nrow(heatmap_temp) > 3) {
			colnames(heatmap_temp) <- colnames(CLASS_TOPHITS)	
			Anne_heatmap(heatmap_temp, experiment, "Heatmap of TopHits of all Class genes (logFC)", "Heatmap_Class_TopHits", "ratio")
		}	

	######################## Heatmap of Each Class TopHits ###################################	
		logbook <- write_log(logbook, "----> Heatmap of Each Class TopHits" )
		Class_GroupNames <- levels(my_class_genes$Group)
		filenames.ClassHeatmaps <- list()
		filenames.Clusters <- list()
		TMEV[1] <- NULL
		for (i in 1:length(Class_GroupNames)) {
			Class_Genes <- subset(my_class_genes$GeneID, my_class_genes$Group %in% Class_GroupNames[i])
			Class_Group <- subset(TMEV, row.names(TMEV) %in% Class_Genes)
			if (dim(Class_Group)[1] > 1) {
				filename_prefix   <- paste("Heatmap_Class", Class_GroupNames[i], sep=".")
				filenames.ClassHeatmaps <- c(filenames.ClassHeatmaps, paste(experiment, filename_prefix ,"png", sep="."))
				filenames.Clusters  <- c(filenames.Clusters, paste(experiment, filename_prefix ,"clusters", "txt", sep="."))
				main_title <- paste(experiment, "Heatmap of Contrasts of Class", Class_GroupNames[i] , sep=" ")
				heatmap_temp <-(t(apply(Class_Group,1,as.numeric)))
				if (nrow(heatmap_temp) > 3) {
					colnames(heatmap_temp) <- colnames(Class_Group)	
					Anne_heatmap(heatmap_temp, experiment, main_title, filename_prefix, "ratio")
				}
			}	
		}

		
		######################## Heatmaps of Class Signal data ###################################	
		logbook <- write_log(logbook, "----> Heatmaps of All Class genes Signal data" )
		# Old: CLASS_TOPHITS_signal <- subset(mean_table_class, row.names(mean_table_class) %in% tophits)
		heatmap_temp <-(t(apply(CLASS_TOPHITS_signal,1,as.numeric)))
		if (nrow(heatmap_temp) > 3) {
		  colnames(heatmap_temp) <- colnames(CLASS_TOPHITS_signal)
		  Anne_heatmap(log(heatmap_temp+1,2), experiment, "Heatmap of all Class genes log2(Signal)", "Heatmap_Class_Signals", "signal")
		}
		

		######################## Save Filenames #################################################
		filename <- paste(experiment, "ClassHeatmaps.filenames.txt", sep=".")
		write.table(unlist(filenames.ClassHeatmaps), file=filename, row.names = FALSE, col.names = FALSE, quote=FALSE)
		# OLD removed april 2020 - write.csv(filenames.ClassHeatmaps, file=filename)

		filename <- paste(experiment, "ClassHeatmaps.clusters.filenames.txt", sep=".")
		temp <- NULL
		temp$files <- sapply(filenames.Clusters, paste0, collapse=",") 
		temp$names <- sapply(Class_GroupNames, paste0, collapse=",") 
		write.table(as.data.frame(temp), file=filename, row.names = FALSE, col.names = FALSE, quote=FALSE)
		

		
		
	# 26. Heatmap of Each Class Signal data
		logbook <- write_log(logbook, "# 26. Heatmap of Each Class Signal data" )
		Class_GroupNames <- levels(my_class_genes$Group)
		filenames.ClassHeatmapsSignal <- list()
		filenames.Clusters <- list()
		for (i in 1:length(Class_GroupNames)) {
			Class_Genes <- subset(my_class_genes$GeneID, my_class_genes$Group %in% Class_GroupNames[i])
			Class_Group <- subset(mean_table_class, row.names(mean_table_class) %in% Class_Genes)
			if (dim(Class_Group)[1] > 1) {
				filename_prefix   <- paste("Heatmap_Class_Signal", Class_GroupNames[i], sep=".")
				filenames.ClassHeatmapsSignal <- c(filenames.ClassHeatmapsSignal, paste(experiment, "Heatmap_Class_Signal", Class_GroupNames[i] ,"png", sep="."))
				filenames.Clusters  <- c(filenames.Clusters, paste(experiment, filename_prefix ,"clusters", "txt", sep="."))
				main_title <- paste(experiment, "Heatmap of Class", Class_GroupNames[i], "log2(Signal)" , sep=" ")
				heatmap_temp <-(t(apply(Class_Group,1,as.numeric)))
				colnames(heatmap_temp) <- colnames(Class_Group)	
				Anne_heatmap(log(heatmap_temp+1,2), experiment, main_title, filename_prefix, "signal")
			}	
		}
		# save filenames
		filename <- paste(experiment, "ClassHeatmapsSignal.filenames.txt", sep=".")
		write.table(unlist(filenames.ClassHeatmapsSignal), file=filename, row.names = FALSE, col.names = FALSE, quote=FALSE)
		#write.csv(filenames.ClassHeatmapsSignal, file=filename)
		
		filename <- paste(experiment, "ClassHeatmapsSignal.clusters.filenames.txt", sep=".")
		temp <- NULL
		temp$files <- sapply(filenames.Clusters, paste0, collapse=",") 
		temp$names <- sapply(Class_GroupNames, paste0, collapse=",") 
		write.table(as.data.frame(temp), file=filename, row.names = FALSE, col.names = FALSE, quote=FALSE)
    
	}
	



  	
##############################################################################################
##                        Line Charts of all Classes                                        ##
##############################################################################################

  	
  # Signal data Expression plots of each Class	
  if (option.classes) {	
	  logbook <- write_log(logbook, "# Signal data expression plots of Classes" )
	  # old: CLASS_TOPHITS_signal <- subset(mean_table_class, row.names(mean_table_class) %in% tophits)
	  Class_signal_plots(log(CLASS_TOPHITS_signal,2),"Mean_signals_subClass_TOPHITS") 	
  } 

  	
  	
##############################################################################################
##                        K-Means Clustering                                                ##
##############################################################################################
	
	if (option.clustering) {
    # K-Means Clustering on Ratios
  	logbook <- write_log(logbook, "----> K-Means Clustering of TopHits (ratios)" )
  	if (length(row.names(TMEV_TOPHITS)) > 10) { KMeans_Clustering(TMEV_TOPHITS, "MedianFold", "log2Ratio") }
  	if (length(row.names(TMEV_TOPHITS_highfold)) > 10) { KMeans_Clustering(TMEV_TOPHITS_highfold, "HighFold", "log2Ratio") }
  	
  	if (option.classes) {
  	  logbook <- write_log(logbook, "----> K-Means Clustering of Classes TopHits (ratios)" )
  	  ## Old: CLASS_TOPHITS <- subset(TMEV_TOPHITS, row.names(TMEV_TOPHITS) %in% my_class_genes$GeneID)
    	if (length(row.names(CLASS_TOPHITS)) > 10) { KMeans_Clustering(CLASS_TOPHITS, "ClassTophits", "log2Ratio") }
    	if (length(row.names(CLASS_TOPHITS_signal)) > 10) { KMeans_Clustering(log(CLASS_TOPHITS_signal+0.1,2), "ClassTophitsSignal", "log2Signal") }
  	}  	
  
    # K-means CLustering on Signals
  	logbook <- write_log(logbook, "---->K-Means Clustering of Mean Signals" )
  	# Clustering all genes is overkill, we will use only genes that changed at least weakly in one of the contrasts
  	my_trimmed_mean_table <-  subset(mean_table, row.names(mean_table) %in% tophits_no_background)
  	# take the mean table and remove rows with 1 zero or more
  	my_trimmed_mean_table <- my_trimmed_mean_table[apply(my_trimmed_mean_table, 1, function(x) !any(x==0)),]
  	head(my_trimmed_mean_table)
  	if (length(row.names(my_trimmed_mean_table)) > 10) { KMeans_Clustering(log(my_trimmed_mean_table,2), "MeanSignal", "log2Signal") }
	}


	
##############################################################################################
##                        Export JSON                                                       ##
##############################################################################################
	
	
  # Finally save the file container for downstream processing
	logbook <- write_log(logbook, "# Write JSON containter file" )
	jSON$settings <- json.settings
	jSON$menus <- json.menus
	jSON_output <- prettify(toJSON(jSON)) #create JSON output
	filename <- paste(experiment, "Container.json", sep=".")
	write(jSON_output, file=filename)	

	
	
##############################################################################################
##                        DONE                                                              ##
##############################################################################################
	
		
logbook <- write_log(logbook, "# R-script DONE" )


