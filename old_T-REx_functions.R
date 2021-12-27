
# Anne de Jong, Aug 2015, T-REx R functions 
# University of Groningen
# Molecular Genetics
# the Netherlands

# 2015 Aug : First release for publication in BMC genomics
# 2015 Aug : Add function: Anne_PCA_plot




theme_g2d_plot= theme(axis.text.x =	element_text(size = 20, colour = "blue", angle = 30), 
						axis.text.y = 	element_text(size = 20),
						title = 		element_text(size = 22, face = "bold", hjust = 0),
						axis.title.x = 	element_text(size = 20, hjust = 0.5),
						axis.title.y = 	element_text(size = 20, hjust = 0.5)
					)
theme_g2d_plot_cor = theme(
						axis.text.x =	element_text(size = 20, angle = 90, hjust = 1), 
						axis.text.y = 	element_text(size = 20, hjust = 1),
						title = 		element_text(size = 32, face = "bold", hjust = 0),
						axis.title.x = 	element_text(size = 28, hjust = 0.5),
						axis.title.y = 	element_text(size = 28, hjust = 0.5)
					)
theme_g2d_barplot = theme(
						axis.text.x =	element_text(size = 18, angle = 45, hjust = 1), 
						axis.text.y = 	element_text(size = 16, hjust = 1),
						title = 		element_text(size = 24, face = "bold", hjust = 0),
						axis.title.x = 	element_text(size = 20, hjust = 0.5),
						axis.title.y = 	element_text(size = 20, hjust = 0.5)
					)
theme_g2d_heatmap = theme(
						axis.text.x =	element_text(size = 18, angle = 45, hjust = 1), 
						axis.text.y = 	element_text(size = 16, hjust = 1),
						title = 		element_text(size = 24, face = "bold", hjust = 0),
						axis.title.x = 	element_text(size = 20, hjust = 0.5),
						axis.title.y = 	element_text(size = 20, hjust = 0.5)
					)

theme_PCA_plot = theme(
						axis.text.x =	element_text(size = 20, hjust = 1), 
						axis.text.y = 	element_text(size = 20, hjust = 1),
						title = 		element_text(size = 36, face = "bold", hjust = 0.5),
						axis.title.x = 	element_text(size = 26, hjust = 0.5),
						axis.title.y = 	element_text(size = 26, hjust = 0.5)
					)


Anne_PCA_plot<- function(pca_data, my_prefix, my_title, ratioaxes="true" ) {	
	# general function for 2D PCA plot
	# pca_data is just a data frame with columns 'experiment', 'PC1', 'PC2'
	# my_title <- "Test"
	# ratioaxes <- "true"
	experiments.pca <- prcomp(pca_data, center = TRUE, scale. = FALSE) 
	experiments.pca$x
	plot(experiments.pca, type = "l")
	summary(experiments.pca)			 
	my_pca_data <- NULL
	my_pca_data <- as.data.frame(experiments.pca$x)
	my_pca_data$experiment <- row.names(my_pca_data)
	head(my_pca_data)
	# add the rotation axis
		my_pca_rotation <- NULL
		my_pca_rotation <- as.data.frame(experiments.pca$rotation)
		my_pca_rotation$experiment <- row.names(my_pca_rotation)
		head(my_pca_rotation)
		rotation_axis_scale <- 5
		result <- lapply(seq(1,nrow(my_pca_rotation)), function(x) {	
			c(0, 0, rotation_axis_scale*my_pca_rotation$PC1[x], rotation_axis_scale*my_pca_rotation$PC2[x])
		})
		rotation_axis <- data.frame(rbind(do.call(rbind,result)))
		colnames(rotation_axis) <- c("x","y","xend","yend")
		rotation_axis <- cbind(my_pca_rotation, rotation_axis)
		rotation_axis


	# the basic plot
	experiment_PCAplot <- ggplot(my_pca_data,  aes(PC1,PC2)) + 
				ggtitle(my_title) + 
				theme_PCA_plot +
				xlab("PC1") +  
				ylab("PC2") +
				theme(legend.position="none") +
				scale_y_continuous(expand = c(0.1,0.1))
	#

	# the rotation axis
	if (ratioaxes == "true") {
		experiment_PCAplot <- experiment_PCAplot + geom_point(colour = "darkblue", size=2)
		experiment_PCAplot <- experiment_PCAplot + geom_segment(data=rotation_axis, mapping=aes(x=x, y=y, xend=xend, yend=yend), size=1, colour = "darkgreen")
		experiment_PCAplot <- experiment_PCAplot + geom_text(data=rotation_axis, aes(label=experiment, x=xend, y=yend), size=10, colour = "darkred") 
	} else {
		fontsize <- ifelse(nrow(my_pca_data)<7, 18, 10)
		experiment_PCAplot <- experiment_PCAplot + geom_text(aes(label=experiment), size=fontsize, colour = "darkblue") 
	}
	
	print(experiment_PCAplot)

	# write plot to file
	filename <- paste(experiment, my_prefix, "png", sep=".")
	png(file = filename, width = x_1k, height = y_1k, units = "px")
	par(mar=c(20,20,20,10)) # bottom, left, top, right
	print(experiment_PCAplot)
	dev.off()
	# write the plot data
	filename <- paste(experiment, my_prefix, "txt", sep=".")
	write.table(my_pca_data, file=filename , quote=F, sep="\t", col.names=NA)	
}


write_log <- function(logbook, my_text) {
	# write for every step of the pipeline to follow the progress
	my_text <- paste(my_text, "<br>")
	logbook <- rbind(logbook, Log=my_text)
	write.table(logbook, file="00.R.log" , quote=F, sep="\t", row.names=F)
	return(logbook)
}

Anne_heatmap <- function(my_data, experiment, my_title, my_prefix, my_type) {
	
	# experiment = session name of the experiment
	# my_data = data.matix object
	# my_title = for the heatmap plot
	# my_prefix = file name prefix 
	# my_type = 'signal' or 'ratio'

	## for testing this function use
	# my_data <- heatmap_temp
	# my_title <- "Heatmap of TopHits Ratio (logFC)"
	# my_prefix <- "Heatmap_TopHits"
	# my_type <- 'signal'

	#my_heatmap_palette_signal <- colorRampPalette(c("#FCFFF5", "#D1DBBD", "#91AA9D", "#3E606F", "#193441"))(n = 1000) 
	my_heatmap_palette_signal <- colorRampPalette(c("#FCFFF5", "#D1DBBD", "#193441"))(n = 1000) 
	my_heatmap_palette_ratio <- colorRampPalette(c("#4C3200", "#CC8400", "white","#00006B", "#00004C"))(n = 1000)  # orange - blue
	my_heatmap_margins <- c(14,10)
	my_heatmap_pointsize <- 24

	# QHD
		x_2k <- 2560
		y_2k <- 1440
		
	# set the custom distance and clustering functions, per your example
	hclustfunc <- function(x) hclust(x, method="complete")
	distfunc <- function(x) dist(x, method="euclidean")

	## - Determine number of clusters
	wss <- (nrow(my_data)-1)*sum(apply(my_data,2,var))
	nclusters <- length(row.names(my_data))
	if (nclusters > 30) { nclusters <- 30 }   # first limitation for max kmeans clusters
	for (i in 2:nclusters-1) wss[i] <- sum(kmeans(my_data,  centers=i)$withinss)
	# estimate optimal number of clusters = min + 10% max
	nclusters <- length(wss[wss>min(wss)+0.05*max(wss)])+1  # if nclus is > 5% above minimum relative to the max value
	if (nclusters > 9) { nclusters <- 9 } # max number of clusters

	d <- distfunc(my_data)
	fit <- hclustfunc(d)
	clusters <- cutree(fit, k=nclusters)
	nofclust.height <-  length(unique(as.vector(clusters)))
	selcol <- colorRampPalette(brewer.pal(nclusters,"Set1"))
	clustcol = selcol(nofclust.height)

	filename <- paste(experiment, my_prefix, "png", sep=".")
	png(file = filename, width = x_2k, height = y_2k, units = "px", pointsize = my_heatmap_pointsize, bg = "white")	
	par(mar=c(1,1,1,0)) # bottom, left, top, right
	if (my_type == "signal") {
		heatmap_data <- heatmap.2(my_data, main=my_title, RowSideColors=clustcol[clusters], col=my_heatmap_palette_signal, xlab=list(cex=.05), dendrogram="row",  key=T, keysize=1, trace="none",scale="none", margins=my_heatmap_margins, density.info="none", cexCol=1.5, Colv=FALSE)
	} else {
		heatmap_data <- heatmap.2(my_data, main=my_title, RowSideColors=clustcol[clusters], col=my_heatmap_palette_ratio, xlab=list(cex=.05), dendrogram="both",  key=T, keysize=1, trace="none",scale="none", margins=my_heatmap_margins, density.info="none", cexCol=1.5)
	}	
	legend(x="topright", legend = c(1:nclusters), col = clustcol, lwd = 10 ,y.intersp= 0.8, box.col = "white" , title = "Clusters" )

	dev.off()

	# save the heatmap data to file
		my_data_sorted <- as.data.frame(t(heatmap_data$carpet))
		my_data_sorted <- my_data_sorted[,order(names(my_data_sorted))]
		filename <- paste(experiment, my_prefix, "sorted", "txt", sep=".")
		write.table(my_data_sorted, file=filename , quote=F, sep="\t", row.names=T, col.names=NA)

	# save the raw data to file
		filename <- paste(experiment, my_prefix, "txt", sep=".")
		write.table(my_data, file=filename , quote=F, sep="\t", row.names=T, col.names=NA)
	# save the cluster info
		clusters <- as.data.frame(clusters)
		clusters$GeneID <- rownames(clusters)
		clusters <- clusters[with(clusters, order(clusters)), ]
		clusters <- clusters[c("GeneID","clusters")]
		filename <- paste(experiment, my_prefix, "clusters", "txt", sep=".")
		write.table(clusters, file=filename , quote=F, sep="\t", row.names=FALSE, col.names=TRUE)
}
					
Class_signal_plots <- function(my_data, my_description) {	
	#my_data<-CLASS_TOPHITS_signal
	#my_description<-"Mean_signals_subClass"
	#draw a graph for each class color
	filenames <- list()
	for (i in levels(my_class_genes$Group)) {
		my_subset_genes <- subset(my_class_genes, my_class_genes$Group==i )	# select subset of genes on the basis of color 
		#print(head(my_subset_genes))
		my_subset_table <- subset(my_data, row.names(my_data) %in% my_subset_genes$GeneID) # select subset of data
		#print(head(my_subset_table))
		if (nrow(my_subset_table)>1) {   # do not draw tables less then 2 genes
			m1 <- melt(my_subset_table, id = row.names)
			#print(head(m1))
			min_y = 0
			max_y <- 2 * median(log(m1$value,2), na.rm=TRUE)
			filename <- paste(experiment, my_description ,i,"png", sep=".")
			filenames[i] <- filename
			png(file = filename, width = x_2k, height = y_2k, units = "px", pointsize = 36, bg = "white")
			my_title=paste(my_description," = ",i)
			print(my_title)
			p1 <- qplot(x=m1$X2, y=log(m1$value,2), data=m1, group=m1$X1,colour = factor(m1$X1), main = my_title, xlab="Factors", ylab="Log2(Signal)", geom = "line") + theme_g2d_plot
			# , ylim=c(min_y,max_y)
			print(p1)
			dev.off()
		}
	}
	# save the filenames of the plots
	filename <- paste(experiment, my_description, "filenames.txt", sep=".")
	#write.csv(filenames, file=filename)
	write.table(unlist(filenames), file=filename, row.names = FALSE, col.names = FALSE, quote=FALSE)

}



SubClass_correlation_plots <- function(my_data, my_description) {	


	#my_data<-cor_mean_table_class
	#my_description<-"Correlation_plots_subClasses"
	#draw a correlation matrix for each subclass based on the color
	filenames <- list()
	for (i in levels(my_class_genes$Group)) {
		my_subset_genes <- subset(my_class_genes, my_class_genes$Group==i )	# select subset of genes on the basis of color 
		#print(head(my_subset_genes))
		my_subset_table <- subset(my_data, row.names(my_data) %in% my_subset_genes$GeneID) # select subset of data
		#print(head(my_subset_table))
		if (nrow(my_subset_table)>1) {   # do not draw tables with less then 2 genes
			cor_subset_class <- cor(t(my_subset_table))
			# write a tab delimited table
			filename <- paste(experiment, my_description ,i,"txt", sep=".")
			write.table(cor_subset_class, file=filename , quote=F, sep="\t", row.names=F)
			jSON$SubClass_correlation_plots[i] <- filename
			# plot the correlation matrix
			m1 <- melt(cor_subset_class, id = row.names)
			filename <- paste(experiment, my_description ,i,"png", sep=".")
			print(filename)
			filenames[i] <- filename
			png(file = filename, width = x_1k, height = y_1k, units = "px", pointsize = 12, bg = "white")
			my_title=paste(my_description," = ",i)
			p1 <- qplot(x=X1, y=X2, data=as.data.frame(m1),main = my_title, xlab="Class genes", ylab="Class genes", fill=value, geom="tile") + theme_g2d_plot_cor + scale_fill_gradient2(limits=c(-1, 1))
			print(p1)
			dev.off()
		}
	}
	# save the filenames of the plots
	filename <- paste(experiment, my_description, "filenames.txt", sep=".")
	#write.csv(filenames, file=filename)
	write.table(unlist(filenames), file=filename, row.names = FALSE, col.names = FALSE, quote=FALSE)

}



KMeans_Clustering <- function(my_data, my_description, y_label) {	
# K-Means Clustering
	theme_kmeansplot= theme(axis.text.x =	element_text(size = 20, colour = "blue", angle = 30), 
							axis.text.y = 	element_text(size = 20),
							title = 		element_text(size = 36, face = "bold", hjust = 0),
							axis.title.x = 	element_text(size = 20, hjust = 0.5),
							axis.title.y = 	element_text(size = 20, vjust = 0.5)
						)
	## - Determine number of clusters
	wss <- (nrow(my_data)-1)*sum(apply(my_data,2,var))
	nclusters <- length(row.names(my_data))
	if (nclusters > 30) { nclusters <- 30 } 
	for (i in 2:nclusters-1) wss[i] <- sum(kmeans(my_data,  centers=i)$withinss)
	filename <- paste(experiment, my_description,"kmeans_estimates", "png", sep=".")
	png(file = filename, width = 1280, height = 1280, units = "px", pointsize = 12, bg = "white")	
	par(mar=c(6,6,13,2), col="black") # bottom, left, top, right
	plot(wss, type="b", xlab="Estimate Number of Clusters",  ylab="Within groups sum of squares") 
	# estimate optimal number of clusters = min + 10% max
	nclusters <- length(wss[wss>min(wss)+0.05*max(wss)])+4  # if nclus is > 5% above minimum relative to the max value
	lines(c(nclusters,nclusters),c(-100,100000), col="red")
	dev.off()

	# K-Means Clustering with nclusters clusters
	fit <- kmeans(my_data, nclusters)

	## - Cluster Plot against 1st 2 principal components
	# filename <- paste(experiment, my_description,"kmeans_ClusterPlot", "png", sep=".")
	# png(file = filename, width = 1280, height = 1280, units = "px", pointsize = 12, bg = "white")	
	# par(mar=c(16,16,13,2), col="black") # bottom, left, top, right
	# clusplot(my_data, fit$cluster, color=TRUE, shade=TRUE,   labels=2, lines=0)
	# dev.off()

	## - Plot the Ward Hierarchical Clustering with k-means
	d <- dist(my_data, method = "euclidean") # distance matrix
	fit <- hclust(d, method="ward")
	filename <- paste(experiment, my_description,"kmeans_Dendrogram", "png", sep=".")
	png(file = filename, width = x_4k, height = y_4k, units = "px", pointsize = 12, bg = "white")	
	par(mar=c(6,6,13,2), col="black") # bottom, left, top, right
	plot(fit) # display dendogram
	groups <- cutree(fit, k=nclusters) # cut tree into 5 clusters
	# draw dendogram with red borders around the  clusters
	rect.hclust(fit, k=nclusters, border="red") 
	dev.off() 

	## - Plot MDS of k-means Clusters
	filename <- paste(experiment, my_description,"kmeans_MDS", "png", sep=".")
	png(file = filename, width = 1280, height = 1280, units = "px", pointsize = 12, bg = "white")	
	par(mar=c(6,6,13,2), col="black") # bottom, left, top, right
	mds<- cmdscale(dist(my_data))
	plot(mds, col=groups, pch=19, main="Distance of genes between K-means clusters")
	dev.off()
	
	
	
	## -  Write the Tophits
	f_tophits_Clusters <- as.data.frame(cbind(my_data, clusterID=groups))
	f_tophits_Clusters <- f_tophits_Clusters[ order(f_tophits_Clusters$clusterID), ]
	filename <- paste(experiment, my_description,"ClusterID", "txt", sep=".")
	write.table(f_tophits_Clusters, file=filename , quote=F, sep="\t", row.names=T, col.names=NA)
	jSON$KMeans_Clustering[my_description] <- filename
	
	## plot the k-means clusters
	plots <- list()  # new empty list
	filenames <- list()
	theme_set(theme_grey(base_size = 18)) 

	for (i in 1:nclusters) {
		# Get the cluster i from TOPHITS, unlist and melt
		cluster <- f_tophits_Clusters[f_tophits_Clusters$clusterID == i,]
		cluster$clusterID <- NULL
		cluster2 <- t(apply(cluster, 1, unlist))
		if (nrow(cluster2)>0) {
			
			# print the plot
			filename <- paste(experiment, my_description,"kmeans_plot_cluster",i, "png", sep=".")
			filenames[i] <- filename
			png(file = filename, width = 1280, height = 1280, units = "px", pointsize = 12)	
			p1_title <- paste("k-means cluster ",i," - ",my_description,sep=" ")

			m1 <- melt(cluster2, id = row.names)
			if( mode(my_data) == "numeric") { 
				max_y <- round_any( max(m1$value) , 1, f=ceiling) # round to whole numbers
				min_y <- round_any( min(m1$value) , 1, f=floor) # round to whole numbers
				p1 = qplot(x=X2, y=value, data=m1, group=X1, main = p1_title, colour = factor(X1), xlab="Experiment", ylab=y_label, geom = "line",  ylim=c(min_y,max_y)) + theme_kmeansplot + scale_colour_hue(name = "GeneID")
			} else {
				max_y <- round_any( max(as.numeric(levels(m1$value)[as.numeric(m1$value)])) , 1, f=ceiling) # round to whole numbers
				min_y <- round_any( min(as.numeric(levels(m1$value)[as.numeric(m1$value)])) , 1, f=floor)
				y_values <- as.numeric(levels(m1$value)[as.numeric(m1$value)])
				p1 = qplot(x=X2, y=y_values, data=m1, group=X1, main = p1_title, colour = factor(X1), xlab="Experiment", ylab=y_label, geom = "line",  ylim=c(min_y,max_y)) + theme_kmeansplot + scale_colour_hue(name = "GeneID")
			}	
			
			print(p1)
			dev.off()	

		}
	}
	# save the filenames
	filename <- paste(experiment, my_description,"kmeans_plot_cluster.filenames.txt", sep=".")
	#write.csv(filenames, file=filename)
	write.table(unlist(filenames), file=filename, row.names = FALSE, col.names = FALSE, quote=FALSE)


}


Draw_GeneNetwork <- function(my_data, my_description, my_layout) {
	# layouts: layout.kamada.kawai, layout.reingold.tilford
	if (nrow(my_data) > 3) { 
		filename <- paste(experiment, my_description,"png", sep=".")
		plottitle <- paste(experiment, "GeneNetwork on Contrasts", sep=" ")
		png(file = filename, width = x_2k, height = y_2k, units = "px", pointsize = 12, bg = "white")	
		par(mar=c(6,6,3,2)) # bottom, left, top, right
		my_network_graph=graph.data.frame(my_data,directed=FALSE)
		E(my_network_graph)$weight=as.numeric(as.character(my_data$logFC))
		E(my_network_graph)[ weight < 0 ]$color <- "#672000"
		E(my_network_graph)[ weight > 0 ]$color <- "#16675E"
		V(my_network_graph)$color<-ifelse(V(my_network_graph)$name  %in% levels(my_data$Contrasts), '#DB9E36', '#BD4932')
		V(my_network_graph)$shape<-ifelse(V(my_network_graph)$name  %in% levels(my_data$Contrasts), "square", "circle")
		V(my_network_graph)$size<-ifelse(V(my_network_graph)$name  %in% levels(my_data$Contrasts), 25, 10)
		p <- plot.igraph(my_network_graph, edge.width=E(my_network_graph)$weight, main=plottitle, layout= my_layout )
		print(p)
		dev.off()
	}	
}