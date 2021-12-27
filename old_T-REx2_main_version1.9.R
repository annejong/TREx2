# Anne de Jong, Sept 2018, T-REx2 main R script 
# University of Groningen
# Molecular Genetics
# the Netherlands



# ----------------------------------- load the library--------------------------------------------------
library(statmod)
library(ggplot2)
library(reshape)
library(cluster)
library(plyr)
library(grid)
library(jsonlite)
library(RColorBrewer)
library(limma)
library(edgeR)
library(gplots)

# ---------------------------------------  If you do not have the packages -------------------------------------------------------- 
# install.packages("ggplot2")
# install.packages("gplots")
# install.packages("statmod")
# install.packages("reshape")
# install.packages("cluster")
# install.packages("plyr")
# install.packages("grid")
# install.packages("jsonlite")
# install.packages("igraph")

# source("http://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("limma")
# biocLite("edgeR")

# -------------------------------------------------------------------------------- for local use------------------------------------------------------ 
sessiondir="G:/My Drive/WERK/T-REx2/examples"
experiment="my_experiment"
rpkm_file="RPKM.txt"
factors_file="Factors.txt"
class_file="Class.txt" 
contrast_file="Contrasts.txt" 
REPLICATES="true"

source("G:/My Drive/WERK/T-REx2/T-REx2_functions.R")
setwd(sessiondir)


# -------------------------------------------------------------------------------- for webserver use ------------------------------------------------- 

# load parameters from command line or call from the webserver
# R command line:  R --vanilla --slave --args . robyn -r bacillus_cereus_merged_normalized.rpkm -f Factors.txt -c Contrasts.txt < /usr/molgentools/rnaseq/RNAseq_multifactor.R
#     sessiondir <- commandArgs()[5]
#     experiment <- commandArgs()[6]
#     rpkm_file  <- commandArgs()[7]
#     factors_file  <- commandArgs()[8]
#     contrast_file <- commandArgs()[9]
#     class_file=commandArgs()[10]
#     # make web folder from the local folder
#     websessiondir <- gsub("/tmp/genome2d", "genome2d_results", sessiondir)

# ----------------------------------- custom settings / or parameters from webserver--------------------------------------------------

# comment SINK if you run this script locally
# sink("00.R.remarks.log", append=TRUE, split=FALSE)

# container for arguments and data tables in jSON format
jSON <- NULL
# creating a menus list for JSON
json.settings <- NULL
json.settings$descriptionsFile <- unbox(paste(websessiondir, "Annotation.table", sep="/"))
json.settings$descriptionSeparator <- unbox("\t")
json.menus <- list()



# graphics settings
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

 
errorfile <- "00.R.sessionprogress.log"

logbook <- c()
logbook <- write_log(logbook, "Pipeline started")


#------------------------------------ Functions --------------------------------------------------

logbook <- write_log(logbook, "Functions loaded")

#------------------------------------ Loading data files --------------------------------------------------

# 1. ===================================================================== RPKM / Count data ========================================== 
    unlink(errorfile)
    logbook <- write_log(logbook, "# 1. Reading the RNA-Seq expression values")
    GeneExpressionTable <- read.delim(rpkm_file,sep="\t", row.names="key")
    # remove rows with missing values
    has.no.na <- apply(GeneExpressionTable, 1, function(x){any(is.na(x))})
    GeneExpressionTable <- GeneExpressionTable[!has.no.na,]
    head(GeneExpressionTable)
    
	# save quantile table
	filename <- paste(experiment, "quantile_table.txt", sep=".")
	write.table(apply(GeneExpressionTable,2,quantile), file=filename , quote=F, sep="\t", row.names=TRUE, col.names=NA)
	

	
	
	
    # Make a Gene Expression Ranked Table
        GeneExpression_RankTable <- GeneExpressionTable
        for (i in 1:ncol(GeneExpression_RankTable)) {
          GeneExpression_RankTable <- GeneExpression_RankTable[order(GeneExpression_RankTable[,i]),] 
          GeneExpression_RankTable[,i] <- seq(1,nrow(GeneExpression_RankTable))
        }
        GeneExpression_RankTable <- GeneExpression_RankTable[order(GeneExpression_RankTable[,1]),] 
        head(GeneExpression_RankTable)
        filename <- paste(experiment, "GeneExpression_RankTable","txt", sep=".")
        write.table(GeneExpression_RankTable, file = filename, quote=F, sep="\t", col.names=NA)
    

    # Heatmap of Rank-order Samples
        logbook <- write_log(logbook, "# 2b. Heatmap of Rank-order Samples" )
        heatmap_temp <-(t(apply(GeneExpression_RankTable,1,as.numeric)))
        if (nrow(heatmap_temp) > 3) {
          colnames(heatmap_temp) <- colnames(GeneExpression_RankTable)
          Anne_heatmap(log(heatmap_temp+1,2), experiment, "Heatmap Rank-order Samples", "Heatmap_Rank_Samples", "signal")
        }	
    
    
    
    
    # BoxPlot of all data columns    
        GeneExpressionTable.melt <- melt(GeneExpressionTable)
        colnames(GeneExpressionTable.melt) =c("Experiment","Value")
        GeneExpressionTable.melt$log2Value <- log(GeneExpressionTable.melt$Value,2)
        head(GeneExpressionTable.melt)
        
        my_barplot_palette <- colorRampPalette(c("#B4AF91",  "#C03000"))(n = length(GeneExpressionTable))
        filename <- paste(experiment, "Boxplot_Raw_Expression_values","png", sep=".")
        png(file = filename, width = x_1k, height = y_1k, units = "px", pointsize = 12, bg = "white")	
        my_plot <- ggplot(GeneExpressionTable.melt,  aes(Experiment,log2Value) ) +	
                        geom_boxplot( aes(fill=Experiment)) +
                        ggtitle("Boxplot of the Raw Expression Values") +
                        theme_g2d_barplot +
                        theme(legend.position="none") +
                        scale_fill_manual(values = my_barplot_palette)
        print(my_plot)
        dev.off()







# 2. ===================================================================== Factors ========================================== 


  # 1. Read the (multi-) factorial design
    	unlink(errorfile)
    	logbook <- write_log(logbook, "# 1. Read the factorial design")
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


# 3. ===================================================================== Classes ==========================================

    	unlink(errorfile)
    	logbook <- write_log(logbook, "# 3. Read the class file")
    	my_class_genes <- read.delim(class_file,sep="\t")
    	colnames(my_class_genes) =c("GeneID","color","Group")
	


# 4. ===================================================================== Normalization ==========================================		

	
	# Normalization 
    	unlink(errorfile)
    	logbook <- write_log(logbook, "# 5. Normalization ")
    	y <- DGEList(counts=GeneExpressionTable,group=MF.group)	
    	y <- calcNormFactors(y)
    	
    	
    	# Create the design matrix
    	designMF <- model.matrix(~0+MF.group,  data=y$samples)
    	colnames(designMF) <- levels(MF.group)
    	designMF
    	
    	# 12. Calculation of Dispersion
    	unlink(errorfile)
    	logbook <- write_log(logbook, "# 12. Calculation of Dispersion")
    	y <- estimateGLMCommonDisp(y,designMF,method="deviance", robust=TRUE, subset=NULL)  # no replicates
    	stored_y <- y
    	y <- estimateGLMCommonDisp(y,designMF)
    	y <- estimateGLMTrendedDisp(y,designMF, method="bin.loess")
    	y <- estimateGLMTagwiseDisp(y,designMF)
    	

    	
    	
    	    	
    # Library size
    	unlink(errorfile)
    	logbook <- write_log(logbook, "# 11. Library size")
    	
    	# prepare data for barplot
    	my_samples <- y$samples
    	my_samples$experiment <- row.names(my_samples)
    	my_samples$lib.size <- as.numeric(my_samples$lib.size)
    	my_samples$experiment <-factor(my_samples$experiment, levels=my_samples[order(my_samples$experiment), "experiment"])
    	#my_barplot_palette <- colorRampPalette(c("#0057FF", "#001640"))(n = dim(my_samples)[1])
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
    	

	


    	
  # Make correlation Matrix of experiments
    	unlink(errorfile)
    	logbook <- write_log(logbook, "# 15. Make correlation Matrix of experiment")
    	logFC <- predFC(y,designMF,prior.count=1,dispersion=0.05)
    	corLogFC <- cor(logFC)
    	filename <- paste(experiment, "Correlation_matrix_experiments.txt", sep=".")
    	write.table(corLogFC, file=filename , quote=F, sep="\t", row.names=F)	
    	
    	filename <- paste(experiment, "Correlation_matrix_experiments.png", sep=".")
    	png(file = filename, width = x_1k, height = y_1k, units = "px", pointsize = 12, bg = "white")
    	m1 <- melt(corLogFC, id = row.names)
    	qplot(x=X1, y=X2, data=as.data.frame(m1),main = "Correlation matrix experiments", xlab="Experiments", ylab="Experiments", fill=value, geom="tile") + theme_g2d_plot_cor
    	dev.off()	
    	
    	

# 5. ===================================================================== Mean values of the replicates ==========================================
    	unlink(errorfile)
    	logbook <- write_log(logbook, "# 6. Calculate mean values of the replicates")
    
    	# Take mean of replicates if are are replicates 
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
    	logbook <- write_log(logbook, "# 6b. Mean Signal data Heatmap" )
    	heatmap_temp <-(t(apply(mean_table,1,as.numeric)))
    	if (nrow(heatmap_temp) > 3) {
    		colnames(heatmap_temp) <- colnames(mean_table)
    		head(heatmap_temp)
    		Anne_heatmap(log(heatmap_temp+1,2), experiment, "Heatmap of Mean Signal log2(Signal)", "Heatmap_Mean_Signals", "signal")
    	}
	
	# Make rank-order table of Mean Signals		
    	logbook <- write_log(logbook, "# 6c. Mean Signal Rank Table and Heatmap" )
    	rank_mean_table <- mean_table
    	for (i in 1:ncol(rank_mean_table)) {
    		rank_mean_table <- rank_mean_table[order(rank_mean_table[,i]),] 
    		rank_mean_table[,i] <- seq(1,nrow(GeneExpressionTable))
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

	
	# Mean Signals Scaled to 100
    	mean_table_log <- log(mean_table,2)
    	head(mean_table_log)
    	filename <- paste(experiment, "Mean_signals_log2.txt", sep=".")
    	write.table(mean_table_log, file=filename , quote=F, sep="\t", col.names=NA)	
    	mean_table_log_scaled <- round((100 * mean_table_log / max(mean_table_log)), digits=3)
    	head(mean_table_log_scaled)
    	filename <- paste(experiment, "Mean_signals_log2_scaled.txt", sep=".")
    	write.table(mean_table_log_scaled, file=filename , quote=F, sep="\t", col.names=NA)
    	
    	# BoxPlot of all data columns    
        	mean_table_log.melt <- melt(mean_table_log)
        	colnames(mean_table_log.melt) =c("Gene","Experiment","Value")
        	head(mean_table_log.melt)
        	
        	my_barplot_palette <- colorRampPalette(c("#B4AF91",  "#C03000"))(n = length(GeneExpressionTable))
        	filename <- paste(experiment, "Boxplot_Normalized_Mean_Expression_values","png", sep=".")
        	png(file = filename, width = x_1k, height = y_1k, units = "px", pointsize = 12, bg = "white")	
        	my_plot <- ggplot(mean_table_log.melt,  aes(Experiment,Value) ) +	
        	  geom_boxplot( aes(fill=Experiment)) +
        	  ggtitle("Boxplot Normalized Mean Expression Values") +
        	  theme_g2d_barplot +
        	  theme(legend.position="none") +
        	  scale_fill_manual(values = my_barplot_palette)
        	print(my_plot)
        	dev.off()
    	
    	    	

    	

# 5. ===================================================================== Add class data ==========================================
  
	# add class (color and groups) to the mean table
    	mean_table_log_scaled_class <- mean_table_log_scaled
    	mean_table_log_scaled_class <- cbind(mean_table_log_scaled_class, row.names(mean_table_log_scaled_class))
    	colnames(mean_table_log_scaled_class) <- c(colnames(mean_table_log_scaled), "GeneID")
    	mean_table_log_scaled_class <- merge(mean_table_log_scaled_class, my_class_genes, by="GeneID",all.x=T)
    	mean_table_log_scaled_class$color <- ifelse(is.na(mean_table_log_scaled_class$color), "grey", as.character(mean_table_log_scaled_class$color))	# add a default color if the color is NA
    	mean_table_log_scaled_class$Group <- ifelse(is.na(mean_table_log_scaled_class$Group), "NA", as.character(mean_table_log_scaled_class$Group))		# add a default Group if the Group is NA
    	head(mean_table_log_scaled_class)
    	filename <- paste(experiment, "Mean_signals_log2_scaled_CLASSES.txt", sep=".")
    	write.table(mean_table_log_scaled_class, file=filename , quote=F, sep="\t", col.names=NA)
	

# 7. Graphics for each subClass on the basis of the group column
    	unlink(errorfile)
    	logbook <- write_log(logbook, "# 7. Signal plots of Class genes")
    	mean_table_class <- subset(mean_table, row.names(mean_table) %in% my_class_genes$GeneID)
    	mean_table_class_log <- subset(mean_table_log, row.names(mean_table_log) %in% my_class_genes$GeneID)
    	head(mean_table_class)
    	if (nrow(mean_table_class)<5) {
    		logbook <- write_log(logbook, "# ERROR - 7. Skipping Class signal plots due to low number of genes in classes")
    	} else {
    		#draw a graph for each class color
    		logbook <- write_log(logbook, "-   a) Create Graphics for mean signals of each subClass")
    		Class_signal_plots(mean_table_class, "Mean_Signals_subClass")
    		logbook <- write_log(logbook, "-   b) Create Graphics for mean log2(signal) of each subClass")
    		Class_signal_plots(mean_table_class_log, "Mean_Signals_subClass_log2")
    	}

# 8. Correlation matrix of All Class genes, using mean signal of replicates
    	unlink(errorfile)
    	logbook <- write_log(logbook, "# 8. Correlation matrix of All Class genes, using mean signal of replicates")
    	cor_mean_table_class <- cor(t(mean_table_class))	
    	filename <- paste(experiment, "Correlation_matrix_Class_genes.txt", sep=".")
    	write.table(cor_mean_table_class, file=filename , quote=F, sep="\t", row.names=F)
    
    	filename <- paste(experiment, "Correlation_matrix_Class_genes.png", sep=".")
    	png(file = filename, width = x_2k, height = y_2k, units = "px", pointsize = 12, bg = "white")
    	m1 <- melt(cor_mean_table_class, id = row.names)
    	qplot(x=X1, y=X2, data=as.data.frame(m1),main = "Correlation matrix of Class genes", xlab="Class genes", ylab="Class genes", fill=value, geom="tile", zlim= c(0,1)) + theme_g2d_plot_cor + scale_fill_gradient2(limits=c(-1, 1))
    	dev.off()

# 9. Correlation matrix of subClass genes, using mean signal of replicates
    	unlink(errorfile)
    	logbook <- write_log(logbook, "# 9. Correlation matrix of subClass genes, using mean signal of replicates")
    	SubClass_correlation_plots(cor_mean_table_class, "SubClass_correlation_matrix" )

	
	
	
# 10. PCA plot of Experiments and Factors
    	unlink(errorfile)
    	logbook <- write_log(logbook, "# 10. PCA plot of Experiments and Factors")
    	# PCA on experiments
    	pca_data <- t(log(y$counts+1)) 
    	Anne_PCA_plot(pca_data, "PCA_experiments", "PCA Plot of Experiments", "false")
    
    	# PCA on experiments mean signals based on the Factors
    	pca_data <- t(log(mean_table+1)) 
    	Anne_PCA_plot(pca_data, "PCA_factors", "PCA Plot of Factors", "false")



	
	

# 6. ===================================================================== Differential Gene Expressions ==========================================
	

# 13. Load the Contrasts
      	unlink(errorfile)
      	logbook <- write_log(logbook, "# 13. Load the Contrasts")
      	# contrast from file
      	my_contrasts <- readLines(contrast_file)
      	my_contrasts
      	contrast.matrix <- makeContrasts(contrasts=my_contrasts, levels=designMF)
      	contrast.matrix
      	contrast.matrix.dim <- dim(contrast.matrix)[2]

# 14. Fit the data
      	unlink(errorfile)
      	logbook <- write_log(logbook, "# 14. Fit the data")
      	fit.design <- glmFit(y, designMF)
      	lrt.design <- glmLRT(fit.design, designMF)
        head(lrt.design$fitted.values)
      	filename <- paste(experiment, "Normalised_data_table.txt", sep=".")
      	write.table(lrt.design$fitted.values, file=filename , quote=F, sep="\t", col.names=NA)
  
      	
# 15 # no replicates: Estimates Dispersion will be used bcv = 0.05
      	if (REPLICATES != "true") {
      	    logbook <- write_log(logbook, "# 15. No Replicates. Estimates Dispersion will be used ")
      	    #REPLICATES == "false"  # no replicates
      	    bcv <- 0.05      	
          	i <- 1 # for testing MF.group.all
          	DGE_TableList <- c()
          	for (i in 1:length(my_contrasts)) {
              	my_list <- as.list(strsplit(my_contrasts[i], "-")[[1]])
              	contrast_columns <- strtoi( rownames(MF.group.all[MF.group.all$Group %in% my_list, ]) ) # Get the rownames (in this case the numbers) of the MF.group.all$Group 
              	DGE_TableList[[i]] <- exactTest(y[,contrast_columns], dispersion=bcv^2) 
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
            # show frist contrast      	
          	i <- 1  # for testing
          	DGE_TableList[[i]]$comparison
          	head(DGE_TableList[[i]]$table)
      	}

# ============================================================================================================================================    	
# ==============================================   Drawing graphics for each Contrast  =======================================================    	
# ============================================================================================================================================    	
      	
    # 16. Loop through the Contrasts and draw for each Contrast a: 
    		## - DE (Differential Expression ) table
    		## - Significant Fold Change Plots
    		## - MA plot
    		## - Volcano plot
    		## - Fold change distribution plot
    		## - Count distribution plot
    	unlink(errorfile)
    	logbook <- write_log(logbook, "# 16. Loop through the Contrasts")
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
      i<-1  # for testing
      REPLICATES = 'true'
    	for (i in 1:contrast.matrix.dim) {
          logbook <- write_log(logbook, my_contrasts[i] )
          
          # A. DE (Differential Expression ) table
          # Replicates
              if (REPLICATES == "true") {
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
          		my_class <- DGE_Table[which(row.names(DGE_Table) %in% my_class_genes$GeneID),]
          		my_class <- merge(my_class, my_class_genes,  by="GeneID",all.x=T)
          		colnames(my_class) =c("GeneID","logFC", "logCPM", "LR", "pvalue", "adj_pvalue", "Fold", "minFDR", "color")
          		# Select data with 0 expression and mark with X it in the plots
             		my_null_class <- DGE_Table[which(row.names(DGE_Table) %in% row.names(lrt$fitted.values)[lrt$fitted.values==0]),]
          		my_null_class
          		

          		
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
                  		#ggplot(DGE_Table_All,  aes(x=logFC,y=minFDR, color=factor(abs(Fold)>1 & adj_pvalue<0.05))) + geom_point()
                  		p <- ggplot(DGE_Table_All,  aes(x=logFC,y=minFDR, color=factor(Group))) + geom_point(size = DGE_Table_All$sfere) +
                  		  geom_text(data=subset(DGE_Table_All, abs(logFC) > 2 & minFDR > 250), aes(label=GeneID),hjust=0, vjust=0) +
                  		  annotate("rect", xmin = xmin, xmax = xmax, ymin = 0, ymax = FDR_threshold, alpha = .2) +
                  		  annotate("rect", xmin = -Fold_threshold, xmax = Fold_threshold, ymin = FDR_threshold, ymax = ymax, alpha = .2) +
                  		  ggtitle(paste(experiment, "\nVolcano Plot ", my_contrasts[i] , sep=" ")) +
                  		  theme(plot.title = element_text(hjust=0.5))
                  		print(p)
          	    	dev.off()
          		
          	    	#  Volcano Plot for only the Class data
          	    	filename  <- paste(experiment, "Volcano_Plot_Class_Genes", i, my_contrasts[i] ,"png", sep=".")
          	      png(file = filename, width = x_1k, height = y_1k,  units = "px", pointsize = 12, bg = "white")
              	    	p<- ggplot(DGE_Table_Class,  aes(x=logFC,y=minFDR, color=factor(Group))) + geom_point(size = DGE_Table_Class$sfere) +
              	    	  geom_text(data=subset(DGE_Table_Class, abs(logFC) > 2 & minFDR > 250), aes(label=GeneID),hjust=0, vjust=0) +
              	    	  annotate("rect", xmin = xmin, xmax = xmax, ymin = 0, ymax = FDR_threshold, alpha = .2) +
              	    	  annotate("rect", xmin = -Fold_threshold, xmax = Fold_threshold, ymin = FDR_threshold, ymax = ymax, alpha = .2) +
              	       	ggtitle(paste(experiment, "\nVolcano Plot of Class Genes ", my_contrasts[i] , sep=" ")) +
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
          		
          # K. - Fold change distribution plots
          		filename  <- paste(experiment, "Fold_distribution", i, my_contrasts[i] ,"png", sep=".")
          		plottitle <- paste("Fold change distribution", experiment, "contrast=", my_contrasts[i] , sep=" ")
          		png(file = filename, width = 1280, height = 1024, units = "px", pointsize = 12, bg = "white")	
          		par(mar=c(6,6,3,2)) # bottom, left, top, right
          		plot(rev(DGE_Table$logFC), cex.main=2,cex.axis=1.8,cex.lab=2.5, main=plottitle, pch=20, xlab="genes (sorted on p-value)", ylab="log2(Fold Change)")
          		rect(0,-1, nrow(DGE_Table),1, col = "red", density = 30, border = "transparent")
          		dev.off()
          
          # L - Count distribution plots
          		filename  <- paste(experiment, "Count_distribution", i, my_contrasts[i] ,"png", sep=".")
          		plottitle <- paste("Count distribution", experiment, "contrast=", my_contrasts[i] , sep=" ")
          		png(file = filename, width = 1280, height = 1024, units = "px", pointsize = 12, bg = "white")	
          		par(mar=c(6,6,3,2)) # bottom, left, top, right
          		plot(sort(DGE_Table$logCPM), cex.main=2,cex.axis=1.8,cex.lab=2.5, main=plottitle, pch=20, xlab="genes used for normalization", ylab="log2(Counts)")
          		rect(0,-100, 0.2*nrow(DGE_Table),100, col = "red", density = 30, border = "transparent")
          		rect(0.8*nrow(DGE_Table),-100, nrow(DGE_Table),100, col = "red", density = 30, border = "transparent")
          		dev.off()		
    	
          		
        }

# ============================================================================================================================================    	
# ============================================================================================================================================    	
	
	
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
	

	# creating a menu where items can be added to, we call it volcano plots
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
	
	
	
# 17. Draw the network using iGraph
	unlink(errorfile)
	logbook <- write_log(logbook, "# 17. Draw Networks" )
	
	# The complete network
	filename <- paste(experiment, "GeneNetwork_EdgeList.txt", sep=".")
	write.table(my_GeneNetwork, file=filename , quote=F, sep="\t", row.names=FALSE)	
	
	my_Error <- tryCatch(
      Draw_GeneNetwork(my_GeneNetwork, "GeneNetwork_of_Contrasts", layout.reingold.tilford ),
      error=function(e) e
	)
	

	# The high fold change network
	my_Error <- tryCatch(
	{
		my_GeneNetwork_high <- my_GeneNetwork[( abs(as.numeric(as.character(my_GeneNetwork$logFC)))>=4 & as.numeric(as.character(my_GeneNetwork$adj_pvalue))<0.01 ),]
		# frequency 
		freq_table <- table(as.character(my_GeneNetwork_high$GeneID))
		
		filename <- paste(experiment, "GeneNetwork_EdgeList_high.txt", sep=".")
		write.table(my_GeneNetwork_high, file=filename , quote=F, sep="\t", row.names=FALSE)
		Draw_GeneNetwork(my_GeneNetwork_high, "GeneNetwork_of_Contrasts_high", layout.reingold.tilford )
	},
		error=function(e) e
	)


		
	
# 18. TMEV: Combine the DGE_Table of all the Contrasts to one file that is compatible with TMEV (TIGR Multi Experiment Viewer)
	
	unlink(errorfile)
	logbook <- write_log(logbook, "# 18. Combine the DGE_Table of all the Contrasts to one file that is compatible with TMEV" )
	colnames(TMEV) = TMEV.colnames
	TMEV <- as.data.frame(TMEV)
	filename <- paste(experiment, "TMEV", "txt", sep=".")
	write.table(TMEV, file=filename , quote=F, sep="\t", row.names=F)
	row.names(TMEV) <- TMEV$GeneID
	head(TMEV)

	TMEV_TOPHITS <- as.data.frame(TMEV[which(row.names(TMEV) %in% tophits),])
	TMEV_TOPHITS <- as.data.frame(TMEV_TOPHITS[ order(TMEV_TOPHITS[,1],decreasing = TRUE), ])
	TMEV_TOPHITS[1] <- NULL  # remove the GeneID column which is now the row.names
	head(TMEV_TOPHITS)

	TMEV_TOPHITS_highfold <- as.data.frame(TMEV[which(row.names(TMEV) %in% tophits_highfold),])
	TMEV_TOPHITS_highfold <- TMEV_TOPHITS_highfold[ order(TMEV_TOPHITS_highfold[,1],decreasing = TRUE), ]
	TMEV_TOPHITS_highfold[1] <- NULL  # remove the GeneID column which is now the row.names
	head(TMEV_TOPHITS_highfold)

	# write data TMEV to files
	filenameTopHits <- paste(experiment, "TMEV_TOPHITS_logFC", "txt", sep=".")
	write.table(TMEV_TOPHITS, file=filenameTopHits , quote=F, sep="\t", row.names=T, col.names=NA)	
	filenameTopHitsHighFold <- paste(experiment, "TMEV_TOPHITS_logFC_highfold", "txt", sep=".")
	write.table(TMEV_TOPHITS_highfold, file=filenameTopHitsHighFold , quote=F, sep="\t", row.names=T, col.names=NA)





# 19 # PCA plots for genes will be send to Admire
	unlink(errorfile)
	logbook <- write_log(logbook, "# 19. PCA plots for genes for Admire" )
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
	


	
	
# ------------------------------------------------------------------------------------------------------------ heatmaps -------------------------------------------------------------------------------------------




# Heatmap settings
	my_heatmap_palette_signal <- colorRampPalette(c("#FCFFF5", "#D1DBBD", "#91AA9D", "#3E606F", "#193441"))(n = 1000) 
	my_heatmap_palette_ratio <- colorRampPalette(c("#4C3200", "#CC8400", "white","#00006B", "#00004C"))(n = 1000)  # orange - blue
	my_heatmap_margins <- c(14,10)
	my_heatmap_pointsize <- 24


# 20. Heatmap of experiments to see how replicates perform
	unlink(errorfile)
	logbook <- write_log(logbook, "# 20. Heatmaps" )
	
	n <- length(row.names(lrt.design$fitted.values))
	normalized_signals <- as.data.frame((lrt.design$fitted.values))
	normalized_signals <- normalized_signals[order(normalized_signals[,1]),]
	normalized_signals <- normalized_signals[round(n*.1):round(n-n*.3),]	# remove high and low 20% counts
	# save the heatmap data to file
		filename <- paste(experiment, "Heatmap_EXPERIMENTS", "txt", sep=".")
		write.table(normalized_signals, file=filename , quote=F, sep="\t", row.names=T, col.names=NA)	
	if (nrow(normalized_signals) > 3) {
		filename <- paste(experiment, "Heatmap_EXPERIMENTS", "png", sep=".")
		png(file = filename, width = x_2k, height = y_2k, units = "px", pointsize = my_heatmap_pointsize, bg = "white")	
		par(mar=c(15,6,13,2)) # bottom, left, top, right
		
		heatmap.2(as.matrix(log2(normalized_signals+0.5)),dendrogram="both",Rowv=TRUE, col=my_heatmap_palette_signal, scale="none", key=T, keysize=1, main="Normalized log2(Signals), Experiments vs Genes",margins=my_heatmap_margins, density.info="none", trace="none",cexCol=0.9, labRow=NA)
		dev.off()
	}	
	

	
	
# 21. Heatmaps of TopHits
	unlink(errorfile)
	if (ncol(TMEV_TOPHITS) > 1 & nrow(TMEV_TOPHITS) > 1) {
			
		# 21 Heatmaps of TopHits and HighFold	
			logbook <- write_log(logbook, "# 21. Draw heatmaps of TopHits" )
			# new json menu for heatmaps
			json.menu <- NULL
			json.menu$name = unbox("Heatmaps of TopHits")	
			json.items <- list()
			
			#heatmap_temp <-as.data.frame(t(apply(as.matrix(TMEV_TOPHITS),1,as.numeric)))
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
			
			# Heatmap High fold change
			#heatmap_temp <-(t(apply(TMEV_TOPHITS_highfold,1,as.numeric)))
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
	

	
		

		# 22. Heatmap: Combine Class data with TopHits data and draw plots of sub classes
			unlink(errorfile)
			logbook <- write_log(logbook, "# 22.  Combine Class data with TopHits data and draw plots and heatmaps of sub classes" )
			CLASS_TOPHITS <- subset(TMEV_TOPHITS, row.names(TMEV_TOPHITS) %in% my_class_genes$GeneID)
			heatmap_temp <- (t(apply(CLASS_TOPHITS,1,as.numeric)))
			if (nrow(heatmap_temp) > 3) {
				colnames(heatmap_temp) <- colnames(CLASS_TOPHITS)	
				Anne_heatmap(heatmap_temp, experiment, "Heatmap of TopHits of all Class genes (logFC)", "Heatmap_Class_TopHits", "ratio")
			}	

		# 23. Heatmap of each Class Ratio
			unlink(errorfile)
			logbook <- write_log(logbook, "# 23. Make Ratio Heatmap of each Class" )
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
			# save filenames
			filename <- paste(experiment, "ClassHeatmaps.filenames.txt", sep=".")
			write.table(unlist(filenames.ClassHeatmaps), file=filename, row.names = FALSE, col.names = FALSE, quote=FALSE)
			write.csv(filenames.ClassHeatmaps, file=filename)

			filename <- paste(experiment, "ClassHeatmaps.clusters.filenames.txt", sep=".")
			temp <- NULL
			temp$files <- sapply(filenames.Clusters, paste0, collapse=",") 
			temp$names <- sapply(Class_GroupNames, paste0, collapse=",") 
			write.table(as.data.frame(temp), file=filename, row.names = FALSE, col.names = FALSE, quote=FALSE)
			

			
		# 24. High Fold change
			logbook <- write_log(logbook, "# 24. High Fold changed genes" )
			CLASS_TOPHITS_highfold <- subset(TMEV_TOPHITS_highfold, row.names(TMEV_TOPHITS) %in% my_class_genes$GeneID)
			filename <- paste(experiment, "CLASS_TOPHITS_highfold_logFC", "txt", sep=".")
			write.table(CLASS_TOPHITS_highfold, file=filename , quote=F, sep="\t", row.names=T, col.names=NA)	

		# 25. Heatmaps of Class Signal data
			logbook <- write_log(logbook, "# 25. Signal data Heat Maps" )
			CLASS_TOPHITS_signal <- subset(mean_table_class, row.names(mean_table_class) %in% tophits)
			heatmap_temp <-(t(apply(CLASS_TOPHITS_signal,1,as.numeric)))
			if (nrow(heatmap_temp) > 3) {
				colnames(heatmap_temp) <- colnames(CLASS_TOPHITS_signal)
				Anne_heatmap(log(heatmap_temp+1,2), experiment, "Heatmap of all Class genes log2(Signal)", "Heatmap_Class_Signals", "signal")
			}
		
		# 25b. Heatmaps of GeneExpression_RankTable
##			logbook <- write_log(logbook, "# 25b. Heatmap of Rank-order Samples" )
##			heatmap_temp <-(t(apply(GeneExpression_RankTable,1,as.numeric)))
##			if (nrow(heatmap_temp) > 3) {
##				colnames(heatmap_temp) <- colnames(GeneExpression_RankTable)
##				Anne_heatmap(log(heatmap_temp+1,2), experiment, "Heatmap Rank-order Samples", "Heatmap_Rank_Samples", "signal")
##			}	

			
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


	
	
	} else {
		logbook <- write_log(logbook, "# 26. Skip heatmaps of TopHits: Only one contrasts found" )
	}
	

	
#--------------------------------------------------------------------- end heatmaps -------------------------------------------------------------------------
	
# 27. Signal data Expression plots	
	logbook <- write_log(logbook, "# 27. Signal data expression plots" )
	CLASS_TOPHITS_signal <- subset(mean_table_class, row.names(mean_table_class) %in% tophits)
	Class_signal_plots(log(CLASS_TOPHITS_signal,2),"Mean_signals_subClass_TOPHITS") 	


	
	
# 28-30. K-Means Clustering
	logbook <- write_log(logbook, "# 28. K-Means Clustering MedianFold" )
	if (length(row.names(TMEV_TOPHITS)) > 10) { KMeans_Clustering(TMEV_TOPHITS, "MedianFold", "log2Ratio") }
	logbook <- write_log(logbook, "# 29. K-Means Clustering Class Tophits" )
	CLASS_TOPHITS <- subset(TMEV_TOPHITS, row.names(TMEV_TOPHITS) %in% my_class_genes$GeneID)
	if (length(row.names(CLASS_TOPHITS)) > 10) { KMeans_Clustering(CLASS_TOPHITS, "ClassTophits", "log2Ratio") }
	logbook <- write_log(logbook, "# 30. K-Means Clustering High Fold" )
	if (length(row.names(TMEV_TOPHITS_highfold)) > 10) { KMeans_Clustering(TMEV_TOPHITS_highfold, "HighFold", "log2Ratio") }
	if (length(row.names(CLASS_TOPHITS_signal)) > 10) { KMeans_Clustering(log(CLASS_TOPHITS_signal+0.1,2), "ClassTophitsSignal", "log2Signal") }

# 31 Signal K-means CLustering
	logbook <- write_log(logbook, "# 31. K-Means Clustering of Mean Signals" )
	# Clustering all genes is overkill, we will use only genes that changed at least weakly in one of the contrasts
	my_trimmed_mean_table <-  subset(mean_table, row.names(mean_table) %in% tophits_no_background)
	# take the mean table and remove rows with 1 zero or more
	my_trimmed_mean_table <- my_trimmed_mean_table[apply(my_trimmed_mean_table, 1, function(x) !any(x==0)),]
	head(my_trimmed_mean_table)
	if (length(row.names(my_trimmed_mean_table)) > 10) { KMeans_Clustering(log(my_trimmed_mean_table,2), "MeanSignal", "log2Signal") }


	
# 32. Finally save the file container for downstream processing
	logbook <- write_log(logbook, "# 32. Write Json containter file" )
	jSON$settings <- json.settings
	jSON$menus <- json.menus
	jSON_output <- prettify(toJSON(jSON)) #create JSON output
	filename <- paste(experiment, "Container.json", sep=".")
	write(jSON_output, file=filename)	
	
logbook <- write_log(logbook, "# R-script DONE" )


