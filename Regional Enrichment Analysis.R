#This function is designed to identify unique cell state differences affiliated with condition category
    # It uses nearest neighbor analysis and UMAP visualization to identify regions in high-dimensional space significsantly enriched for a particular condition
    # It then reclassifies by this enrichment status, and can compute signal intensity comparisons across given markers


#Load Library
library(pheatmap)
library(uwot) 
library(FNN)
library(prodlim)
library(dplyr)    
library(ggplot2)  
library(reshape2)
library(scales)  
library(viridis)
library(nnet)   
library(samr)

############## Load Function 1 of 2
regionalEnrichmentAnalysis <- function ( singleCellData, 
                                         saveDirectory, 
                                         saveName, 
                                         markers,
                                         conditions,
                                         cell_subset_number = NA,
                                         by_min_sampleN = FALSE,
                                         groupSize = 15,
                                         return_result = TRUE,

                                         EuclideanDistance = TRUE,
                                         printUMAP_cluster = TRUE,
                                         printUMAP_condition = TRUE,
                                         printUMAP_markers = TRUE,
                                         print_NeighborhoodHeatMap = FALSE,
                                         print_NeighborhoodUMAP = TRUE
                                        ) {
  
  
  
  #~~# Create Save Directory
  
  dir.create(paste0(saveDirectory,Sys.Date(),"_",saveName,"_Output"))
  
  
 #~~# Clean the dataset
  nonMarkerColumns = c("cluster","sampleID", "condition")
          #These elements are non-negotiable. sampleID MUST be the individual identifier
          #NO underscores allowed in the condition element, or ideally any of these
          singleCellData$cluster <- gsub("_",".", singleCellData$cluster)
          singleCellData$sampleID <- gsub("_",".", singleCellData$sampleID)
          singleCellData$condition <- gsub("_",".", singleCellData$condition)
                conditions <- gsub("_",".", conditions)
          singleCellData$sampleID <- paste(singleCellData$condition, singleCellData$sampleID,sep="_")
          
  if(is.null(markers)){
      markers <- colnames(singleCellData)[1:(which(colnames(singleCellData) == "cluster")-1)]
  }#ERROR: depends on cluster being present
  
  clean_SCData <- singleCellData[,c(which(colnames(singleCellData) %in% markers),which(colnames(singleCellData) %in% nonMarkerColumns))]
      
  subset_SCData <- data.frame()
  if(!(is.na(cell_subset_number))){
    for(j in unique(clean_SCData$condition)){
      pullDat <- sample_n(clean_SCData[which(clean_SCData$condition == j),], cell_subset_number)
      subset_SCData <- rbind(subset_SCData, pullDat)
    }
  } else if (by_min_sampleN){
    all <- as.numeric(table(gsub("\\_.*$","",unique(clean_SCData$sampleID)))); min_sampleN = min(all)
    for(p in unique(clean_SCData$condition)){
      sample_N <- length(unique(clean_SCData$sampleID[which(clean_SCData$condition == p)]))
      cell_subset_number <- (min_sampleN/sample_N) * length(which(clean_SCData$condition == p))
      pullDat <- sample_n(clean_SCData[which(clean_SCData$condition == p),], cell_subset_number)
      subset_SCData <- rbind(subset_SCData, pullDat)
    }
  } else{
    subset_SCData <- clean_SCData
  }
  
if(is.null(conditions)){conditions <- unique(subset_SCData$conditions)
    conditions <- conditions[order(conditions, decreasing = TRUE)]
      if(any(grepl("Untreat",conditions,ignore.case=TRUE))){
        conditions <- conditions[c(which(grepl("Untreat",conditions,ignore.case=TRUE)),
                                   which(!(grepl("Untreat",conditions,ignore.case=TRUE))))]
      } else if(any(grepl("Control",conditions,ignore.case=TRUE))){
        conditions <- conditions[c(which(grepl("Control",conditions,ignore.case=TRUE)),
                                   which(!(grepl("Control",conditions,ignore.case=TRUE))))]
      }
} 
  
  
  
  #~~# Verify color inputs
  
  tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", 
                  "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", 
                  "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
  
  color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  color = color[-which(grepl("white", color))]
  
  if(!(length(tol21rainbow) >= length(conditions))){
    condition_colors = sample(color, length(conditions))
  } else{condition_colors = tol21rainbow[1:length(conditions)]}
  condition_colors[1] = "gray27"
  
 
   
#~~# Compute Euclidean Distance between conditions
  
  if(EuclideanDistance) {
    print("Starting euclidean distance analysis")
    dist_results <- subset_SCData[,which(colnames(subset_SCData) %in% nonMarkerColumns)]; dist_results$distance <- NA
        base_Vector <- colMeans(subset_SCData[which(subset_SCData$condition == conditions[1]),-which(colnames(subset_SCData) %in% nonMarkerColumns)])
    
        for(p in 1:nrow(dist_results)){
          dist_results$distance[p] <- as.numeric(dist(rbind(base_Vector, subset_SCData[p,-which(colnames(subset_SCData) %in% nonMarkerColumns)])))
        }
    
    #Hist plot
        dist_results$condition <- factor(dist_results$condition, levels = conditions)
        hist_plot <- ggplot(dist_results, aes(x = distance, color = condition)) + geom_density(lwd = 1) +
          theme_classic() + ylab("Relative Number of Cells") + xlab("Euclidean Distance") +
          scale_color_manual(values = condition_colors) + 
          ggtitle(paste("Histogram of euclidean distances \n from avg.",conditions[1],"marker expression"))
        ggsave(hist_plot, filename = paste0(saveDirectory,Sys.Date(),"_",saveName,"_Output/", 
                                            "Histogram-of-distances-from-",conditions[1],".pdf"), 
                                             width = 6, height = 5, useDingbats = FALSE)
    
    #Boxplot
        dist_results_2 <- aggregate(dist_results[, ncol(dist_results)],list(dist_results$sampleID), mean)
        colnames(dist_results_2) <- c("sampleID", "distance")
        dist_results_2$condition <- factor(as.character(gsub("\\_.*$","",dist_results_2$sampleID)), levels = conditions) 
              stat_Output <- character()
              for(c in 2:length(conditions)){
                test <- wilcox.test(dist_results_2$distance[which(dist_results_2$condition == conditions[1])], 
                                    dist_results_2$distance[which(dist_results_2$condition == conditions[c])]); p <- test$p.value
                stat_Output <- paste(stat_Output, "\n    vs",conditions[c],":",round(p,5))
              }
        box_plot <- ggplot(dist_results_2, aes(x = condition, y = distance, color = condition)) + 
          geom_boxplot(outlier.size = 0, lwd=1) +
          geom_point (size = 2)+
          theme_classic() +
          scale_color_manual(values = condition_colors) +
          ggtitle(paste("Euclidean distances from avg.",conditions[1],stat_Output))
        ggsave(box_plot, filename = paste0(saveDirectory,Sys.Date(),"_",saveName,"_Output/", 
                                           "Euclidean-distances-from-",conditions[1],"-by-mouse.pdf"), 
                                            height = 5, width = 5, useDingbats = FALSE)
        
    print("Completed euclidean distance analysis")
  }
  
  
  
  #~~# UMAP plotting
  print("Starting UMAP analysis")
  Immune_Map <- umap(subset_SCData[,-which(colnames(subset_SCData) %in% nonMarkerColumns)], n_neighbors = groupSize,  learning_rate = 0.5) 
        subset_SCData <- cbind(subset_SCData,Immune_Map)
        colnames(subset_SCData)[c(ncol(subset_SCData)-1,ncol(subset_SCData))] <- c("UMAP1","UMAP2")
        added <- 2
        write.csv(subset_SCData, file = paste0(saveDirectory,Sys.Date(),"_",saveName,"_Output/","UMAP_Results.csv"), row.names = FALSE)
        
        subset_SCData$condition <- factor(subset_SCData$condition, levels = conditions) 
        #By Condition
        if(printUMAP_condition){
          UMAP_condition <- ggplot(subset_SCData, aes(UMAP1, UMAP2, color = condition)) + 
            geom_point(stroke = 0.3, size = 0.3) + 
            theme_classic() + scale_color_manual(values = condition_colors) + 
            ggtitle(paste(saveName, "by Condition"))
          ggsave(UMAP_condition, filename = paste0(saveDirectory,Sys.Date(),"_",saveName,"_Output/", 
                                                   "UMAP_colored-by-condition.pdf"), height = 5, width = 6.5, useDingbats = FALSE) 
          
          p <-  ggplot(subset_SCData, aes(UMAP1, UMAP2, color = condition)) + geom_point(stroke = 0.2, size = 0.2) + 
            theme_classic() + scale_color_manual(values = condition_colors) + facet_wrap(~condition) +
            ggtitle(paste(saveName, "by Condition"))
          ggsave(p,filename = paste0(saveDirectory,Sys.Date(),"_",saveName,"_Output/","UMAP_split-by-condition.pdf"), 
                 width = 6, height = 5, useDingbats=FALSE)
        }
        
        #By Cluster
        if(printUMAP_cluster){
          if(!(length(tol21rainbow) >= length(unique(subset_SCData$cluster)))){
            cluster_colors = sample(color, unique(subset_SCData$cluster))
          } else{cluster_colors = tol21rainbow} # ERROR CHECK 3/24/21 some error here if clusters greater
          
          UMAP_cluster <- ggplot(subset_SCData, aes(UMAP1, UMAP2, color = cluster)) + 
            geom_point(stroke = 0.3, size = 0.3) + 
            theme_classic() + scale_color_manual(values = cluster_colors) + ggtitle(paste(saveName, "by Cluster"))
          ggsave(UMAP_cluster, filename = paste0(saveDirectory,Sys.Date(),"_",saveName,"_Output/", 
                                                 "UMAP-by-cluster.pdf"), height = 5, width = 6, useDingbats = FALSE)
          
          cluster_barplot <-ggplot(subset_SCData, aes(x = cluster, fill = condition)) + geom_bar() + 
            scale_fill_manual(values = condition_colors) + theme_classic() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") +
            ggtitle("Frequency of clusters for each condition")
          ggsave(cluster_barplot, filename = paste0(saveDirectory,Sys.Date(),"_",saveName,"_Output/", 
                                                    "Condition-frequency-within-clusters.pdf"), height = 5, width = 5, useDingbats = FALSE)
        }
        
        #By marker
        if(printUMAP_markers){
        dir.create(paste0(saveDirectory,Sys.Date(),"_",saveName,"_Output/marker-expression-plots/"))
        for(j in markers){
          upper <- as.numeric(quantile(subset_SCData[,which(colnames(subset_SCData) == j)], prob=c(.75)))
          express_plot <- ggplot(subset_SCData, aes(UMAP1, UMAP2, color = get(j))) + geom_point(size = 0.3) + theme_classic() +
            scale_color_viridis(limits = c(0, upper), oob = scales::squish) + ggtitle(j) + labs(color = "Signal Intensity \n (max 4th quart.)")
          ggsave(express_plot, filename = paste0(saveDirectory,Sys.Date(),"_",saveName,"_Output/",
                                                 "marker-expression-plots/","UMAP_",j,".pdf"), height = 5, width = 6, useDingbats = FALSE)
        }
        }
        print("Completed UMAP analysis")
        
        
  #~~# Find Nearest Neighbors
 print("Starting nearest neighbor analysis")
  nn_subset <- subset_SCData[,c(1:(ncol(subset_SCData)-(length(nonMarkerColumns)+added)))]
        nearest_neighbors <- get.knn(nn_subset, k=groupSize, algorithm= "kd_tree")
        nn_results <- nearest_neighbors$nn.index; colnames(nn_results) <- paste0("neighbor",c(1:groupSize))
        
        nn_identities <- data.frame(matrix(nrow = 0, ncol = length(conditions),dimnames=list(c(), conditions)))
        for(p in 1:nrow(nn_results)){
          nn_identities[p,] <- 0
          neighborhood <- as.factor(subset_SCData$condition[as.integer(nn_results[p,])]); sum <- summary(neighborhood)
          for(j in 1:length(names(sum))){
            nn_identities[p,which(colnames(nn_identities) == names(sum)[j])] <- as.integer(sum[j])
          }
        }
        subset_SCData<- cbind(subset_SCData,nn_identities); added <- added+ncol(nn_identities)
        print("Completed neartest neighbor analysis")        
        
        
        
 #~~# Define Enrichment Neighborhood
        print("Calculating neighbor cutoffs")
        enrich_cutoffs <- data.frame(condition = conditions); enrich_cutoffs$cutoff <- NA
        for(c in conditions){
          freq <- length(which(subset_SCData$condition == c)) / nrow(subset_SCData)
          this_s <- numeric()
          for(s in 1:groupSize){this_s <- append(this_s,dbinom(s, size = groupSize, freq))}
          for(L in 1:length(this_s)){
            if(is.na(enrich_cutoffs$cutoff[which(enrich_cutoffs$condition == c)])){
              this_L <- sum(this_s[L:length(this_s)])
              if(this_L < 0.05){
                enrich_cutoffs$cutoff[which(enrich_cutoffs$condition == c)] = L
              }
            }
          }
        }
        nn_identities_binomial <- nn_identities
        for(i in 1:ncol(nn_identities_binomial)){
          e.point <- enrich_cutoffs$cutoff[which(enrich_cutoffs$condition == colnames(nn_identities_binomial)[i]) ]
          nn_identities_binomial[which(nn_identities_binomial[,i] < e.point), i] = 0
          nn_identities_binomial[which(nn_identities_binomial[,i] >= e.point), i] = 1
        }
        neighoborhood_types <- nn_identities_binomial %>% group_by_all() %>% summarise(COUNT = n()) 
        neighoborhood_types$Frequency <- round(neighoborhood_types$COUNT / sum(neighoborhood_types$COUNT)*100,2)
        neighoborhood_types <- neighoborhood_types[order(neighoborhood_types$COUNT, decreasing = TRUE),]
      
        # Ask user to define the group name for each neighborhood
        userDefined <- as.character(readline("Would you like to label enrichment groups? (y or n)    "))
        
        if(userDefined != "y" & userDefined != "n"){
          while(userDefined != "y" & userDefined != "n") {
            print("Please respond: y or n")
            userDefined <- readline("Would you like to label enrichment groups?    ")  
          }
        }
      
        neighoborhood_types$EnrichmentCategory <- NA
        categories <- character()
        if(userDefined == "y"){
        for(n in 1:nrow(neighoborhood_types)){
          if(sum(neighoborhood_types[n,1:length(conditions)]) != 0){
          print(neighoborhood_types[n,])
          inputCategory <- readline("What would you like to name this group?    ")  
          categories <- append(categories, inputCategory)
          neighoborhood_types$EnrichmentCategory[n] <- inputCategory
          }
        }
        print("Group labeling complete")
        } else if(userDefined == "n") {
          categories <- paste0("Group", c(1:nrow(neighoborhood_types)))
          neighoborhood_types$EnrichmentCategory <- categories
        }
        categories <- unique(categories)
        if(any(rowSums(neighoborhood_types[,1:length(conditions)]) == 0)){
          neighoborhood_types$EnrichmentCategory[which(rowSums(neighoborhood_types[,1:length(conditions)]) == 0)] = "Unenriched"
          categories <- append(categories, "Unenriched")
        }
        write.csv(neighoborhood_types, file = paste0(saveDirectory,Sys.Date(),"_",saveName,"_Output/","Enrichment_Categories.csv"), row.names = FALSE)
        
        
        #Append enrichment categories to main dataframe
        subset_SCData$EnrichmentCategory <- NA
        suppressWarnings({ 
          for(d in 1:nrow(neighoborhood_types)){
            subset_SCData$EnrichmentCategory[which(row.match(nn_identities_binomial, neighoborhood_types[d,1:length(conditions)]) == 1)] = neighoborhood_types$EnrichmentCategory[d]
          }
        })
        write.csv(subset_SCData, file = paste0(saveDirectory,Sys.Date(),"_",saveName,"_Output/","Regional-Enrichment-Analysis_Results.csv"), row.names = FALSE)
        
        # Plot summaries of enrichment groups
        neighoborhood_types$EnrichmentCategory <- factor(neighoborhood_types$EnrichmentCategory, levels = categories)
        
        print("Plotting neighborhood summaries")
        REA_colors <- character()
        if(all(categories[-which(categories == "Unenriched")] %in% conditions)){
          for(c in categories){
            REA_colors <- append(REA_colors, condition_colors[which(conditions == c)])
          }
          REA_colors <- append(REA_colors, "#e1f7fc")
        } else{
          if(!(length(tol21rainbow) >= length(categories))){
            REA_colors = sample(color, length(categories))
          } else{REA_colors = tol21rainbow[c((length(tol21rainbow)-length(categories)+1):length(tol21rainbow))]}
          REA_colors[length(REA_colors)] = "#e1f7fc"
        }
        
        identity_plot <- ggplot(neighoborhood_types, aes(x="", y=Frequency, fill = EnrichmentCategory)) +
              geom_bar(width = 2, stat = "identity") + coord_polar("y", start=0) +
              theme_classic() + scale_fill_manual(values = REA_colors) +
              ggtitle("Regional Enrichment Category Breakdown")
         ggsave(identity_plot, filename = paste0(saveDirectory,Sys.Date(),"_",saveName,"_Output/",
                                                "Regional-Enrichment-Category-Frequencies.pdf"), 
                                                height = 5, width = 7, useDingbats = FALSE)

        subset_SCData$EnrichmentCategory <- factor(subset_SCData$EnrichmentCategory, levels = rev(categories))
        identity_barplot <-ggplot(subset_SCData, aes(x = cluster, fill = EnrichmentCategory)) + geom_bar() +
              scale_fill_manual(values = rev(REA_colors)) + theme_classic() +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") +
              ggtitle("Regional Enrichment Analsis by Cluster")
        ggsave(identity_barplot, filename = paste0(saveDirectory,Sys.Date(),"_",saveName,"_Output/",
                                                   "Regional-Enrichment-Category-by-clusters.pdf"),
                                                    height = 5, width = 7, useDingbats = FALSE)

  
        
        
 #~~# Plot Neighborhood Barcodes
        # To do: adjsut the barcode to adjust for where the enrichment cutoff is
        
        if(print_NeighborhoodHeatMap){

          neighbor_grouped <- subset_SCData[,c((ncol(subset_SCData) - (ncol(nn_identities))):ncol(subset_SCData))]
          neighoborhoods_all_ordered <- neighbor_grouped %>% group_by_all()
          neighoborhoods_all_ordered <- neighoborhoods_all_ordered[order(neighoborhoods_all_ordered$EnrichmentCategory),]

          for(p in categories){
            group <- neighoborhoods_all_ordered[which(neighoborhoods_all_ordered$EnrichmentCategory == p),]
              index <- which.is.max(group[1,-c((ncol(group)-1), ncol(group))])
              order <- order( as.numeric( unlist(group[,index])))
              neighoborhoods_all_ordered[which(neighoborhoods_all_ordered$EnrichmentCategory == p),] = group[order,]
          }

         # neighoborhoods_all_ordered$EnrichmentCategory <- factor(neighoborhoods_all_ordered$EnrichmentCategory, levels = categories)

          inputDat <- as.data.frame(t(neighoborhoods_all_ordered[,1:length(conditions)]))
          interpretation <- as.data.frame(neighoborhoods_all_ordered$EnrichmentCategory)
          rownames(interpretation) <- colnames(inputDat); colnames(interpretation) <- "group"
          colfunc <- colorRampPalette(c("white", "#3e48fa")); heatmap_cols <- colfunc(groupSize)
          newCols <- REA_colors; names(newCols) <- categories; newCols <- list(group = newCols)

          pheatmap(inputDat, color = heatmap_cols,
                   annotation = interpretation,
                   annotation_colors = newCols, show_colnames = F,
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   main = "Neighborhoods",
                   filename = paste0(saveDirectory,Sys.Date(),"_",saveName,"_Output/",
                                     "Neighborhood-Heatmap.pdf"), width = 10, heigh = 5, useDingbats = FALSE)
        }
      
      
  
  #~~# Plot Regional Enrichment Map

  if(print_NeighborhoodUMAP){
    classification_UMAP <- ggplot(subset_SCData, aes(UMAP1, UMAP2, color = EnrichmentCategory)) + 
          geom_point(stroke = 0.3, size = 0.3) + 
          theme_classic() + scale_color_manual(values = rev(REA_colors)) + 
          ggtitle("Regional Enrichment Analysis")
    ggsave(classification_UMAP, filename = paste0(saveDirectory,Sys.Date(),"_",saveName,"_Output/",
                                                  "UMAP_Regional-Enrichment-Analysis.pdf"),
                                                   height = 5, width = 7, useDingbats = FALSE)
    if(printUMAP_cluster){
      identity_barplot_2 <- ggplot(subset_SCData, aes(x = cluster, fill = EnrichmentCategory)) + geom_bar() + 
            scale_fill_manual(values = rev(REA_colors)) + theme_classic() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") +
            ggtitle("Regional Enrichment Categories by Cluster")
      ggsave(identity_barplot_2, filename = paste0(saveDirectory,Sys.Date(),"_",saveName,"_Output/",
                                                   "Regional-Enrichment-Categories-by-Cluster.pdf"),
                                                    height = 5, width = 7, useDingbats = FALSE)
    }
  }
 
  print("Regional Enrichment Analysis Complete.")
  if(return_result){return(subset_SCData)}
}








############## Load Function 2 of 2
REA_differentiating_Markers <- function (REA_Results,
                                         saveDirectory, 
                                         saveName, 
                                         
                                         arcsinh_transformed = TRUE, #Yes if the data came out of scaffold clustering
                                         cluster.of.interest = NULL,
                                         markers = NULL,
                                         groups.to.test = NULL, #enrichment groups you are interested in comparing
                                         compare.to = NULL, #Default comparison to unenriched. Can change to a particular group here if so desired.
                                         representation.cutoff = 100
                                         ){
  
  
  
  #~~# Initializing needed function
  # SAMR stat analysis 
    compare.markers.between.neighborhoods <- function(dat,group1,group2, proteins,
                                                      arcsinh_transformed){
      arcsinh_trans.value = 5
      dat1 = dat[which(dat$EnrichmentCategory == group1),]
      dat2 = dat[which(dat$EnrichmentCategory  == group2),]
      
      remove_dat1 = which(!(colnames(dat1) %in% append(proteins, c("EnrichmentCategory","sampleID"))))
      if(length(remove_dat1) > 0) {dat1 = dat1[,-remove_dat1]}
      remove_dat2 = which(!(colnames(dat2) %in% append(proteins, c("EnrichmentCategory","sampleID"))))
      if(length(remove_dat2) > 0) {dat2 = dat2[,-remove_dat2]}
      
      # Order & Untransform data  
        get1 = c(which(colnames(dat1) == "sampleID"), which(colnames(dat1) == "EnrichmentCategory"))
      dat1 = dat1[,c(get1, c(1:ncol(dat1))[-which(1:ncol(dat1) %in% get1)])]
        get2 = c(which(colnames(dat2) == "sampleID"), which(colnames(dat2) == "EnrichmentCategory"))
      dat2 = dat2[,c(get2, c(1:ncol(dat2))[-which(1:ncol(dat2) %in% get2)])]
      
      if(arcsinh_transformed){
        dat1[,3:ncol(dat1)] = sinh(dat1[,3:ncol(dat1)]) * arcsinh_trans.value
        dat2[,3:ncol(dat2)] = sinh(dat2[,3:ncol(dat2)]) * arcsinh_trans.value
      }
      
      Protein = names(dat1[3:ncol(dat1)])
      output = as.data.frame(Protein)
      Median1 = numeric()
      Median2= numeric()
      for(i in 3:length(dat1)) {
        g1 = dat1[[i]]; Median1 = append(Median1, median(g1))
        g2 = dat2[[i]]; Median2 = append(Median2, median(g2))
      }
      output = cbind(output, Group1_MSI = Median1, Group2_MSI = Median2)
      
      
      # # # Read in SAMR sub-functions
      
      samTable = function(dat1, dat2) {
        sumStats1 = data.frame(matrix(nrow=ncol(dat1)-2))
        rownames(sumStats1) = names(dat1[3:ncol(dat1)])
        files = unique(dat1$sampleID)
        for(i in 1:length(files)) {
          x = dat1[which(dat1$sampleID == files[i]), 3:ncol(dat1)]
          sumStats1[i] = apply(x, 2, function(x) median(x))
          names(sumStats1)[i] = paste(files[i], "G1", sep ="_")
        }
        sumStats2 = data.frame(matrix(nrow=ncol(dat2)-2))
        rownames(sumStats2) = names(dat2[3:ncol(dat2)])
        files2 = unique(dat2$sampleID)
        for(i in 1:length(files2)) {
          x = dat2[which(dat2$sampleID == files2[i]), 3:ncol(dat2)]
          sumStats2[i] = apply(x, 2, function(x) median(x))
          names(sumStats2)[i] = paste(files2[i], "G2", sep ="_")
        }
        sumStats.all <- cbind(sumStats1, sumStats2)
        return(sumStats.all)
      }
      
      getSampleIDs = function(names) {
        sampleID = rep.int(0, length(names))
        sampleID[grep("G1", names)] = 1
        sampleID[grep("G2", names)] = 2
        if (!all(sampleID!=0)) {print("There is an unassigned file")}
        return(sampleID)
      }
      
      runStats = function(dat1, dat2, runName) {
        sumStats = samTable(dat1, dat2)
        sampleID = getSampleIDs(colnames(sumStats))
        fullMatrix = as.data.frame(sumStats)
        
        row_without_values = apply(fullMatrix, 1, function(row) all(row==0))
        row_without_values_indecies = which(row_without_values == TRUE)
        if (length(row_without_values_indecies) > 0) {statMatrix = fullMatrix[-row_without_values_indecies,]
        } else {statMatrix = fullMatrix}
        
        if(any(is.na(colSums(statMatrix)))){
          sampleID = sampleID[-which(is.na(colSums(statMatrix)))]
          statMatrix = as.matrix(statMatrix[,-which(is.na(colSums(statMatrix)))])
        } else {statMatrix = as.matrix(statMatrix)
        }
        samResults = SAM(x=statMatrix,y=sampleID,resp.type="Two class unpaired",
                         genenames=rownames(statMatrix), geneid=rownames(statMatrix), nperms=10000)
        qStat = as.data.frame(rep(NA, length=nrow(sumStats)))
        rownames(qStat) = rownames(sumStats)
        colnames(qStat) = "SAMqval"
        qStat$SAMfoldChange = rep(NA, length=nrow(sumStats))
        UP = as.data.frame(samResults$siggenes.table$genes.up, stringsAsFactors = FALSE)
        DOWN = as.data.frame(samResults$siggenes.table$genes.lo, stringsAsFactors = FALSE)
        for(p in 1:nrow(UP)) {
          qStat$SAMqval[which(rownames(qStat) == UP$`Gene ID`[p])] = UP$`q-value(%)`[p]
          qStat$SAMfoldChange[which(rownames(qStat) == UP$`Gene ID`[p])] = UP$`Fold Change`[p]
        }
        for(p in 1:nrow(DOWN)) {
          qStat$SAMqval[which(rownames(qStat) == DOWN$`Gene ID`[p])] = DOWN$`q-value(%)`[p]
          qStat$SAMfoldChange[which(rownames(qStat) == DOWN$`Gene ID`[p])] = DOWN$`Fold Change`[p]
        }
        return(qStat)
      }
      # # # Execution
      SAMresult = runStats(dat1,dat2,runName)
      output <- cbind(output, SAMresult)
      return(output) 
    }
  
  
  
  #~~# Set save directory
  save.space <- paste0(saveDirectory,Sys.Date(),"_Diff.Marker.Analysis_",saveName)
  dir.create(save.space)
  
  
  #~~# Organize dataset
  if(is.null(markers)){
    if("cluster" %in% colnames(REA_Results)){
      markers <- colnames(REA_Results)[1:(which(colnames(REA_Results) == "cluster")-1)]
    } else{
      markers <- colnames(REA_Results)[1:(which(colnames(REA_Results) == "sampleID")-1)]
    }
  }
  
  REA_Results.clean <- REA_Results
  if(!(is.null(cluster.of.interest))){
    if("cluster" %in% colnames(REA_Results.clean)){
      if(all(cluster.of.interest %in% REA_Results.clean$cluster)){
        REA_Results.clean <- REA_Results.clean[which(REA_Results.clean$cluster %in% cluster.of.interest), ]
      } else{
        print("The specified cluster.of.interest was not found. Defaulting to all cells.")
      }
    } else{
      print("No 'cluster' column was found. Defaulting to all cells.")
    }
  }

  REA_Results.clean <- REA_Results.clean[,which(colnames(REA_Results.clean) %in% markers | colnames(REA_Results.clean) == "EnrichmentCategory" |
                                              colnames(REA_Results.clean) == "sampleID")]
  
  if(is.null(compare.to)){
    if(!("Unenriched" %in% REA_Results.clean$EnrichmentCategory)){
      stop("Error. Could not find unenriched cells")
    } else{
      compare.to <- "Unenriched"
    }
  }
  if (!(compare.to %in% REA_Results.clean$EnrichmentCategory)){
    stop("Error. Main comparison group is not properly defined")
  }
  
  if(!(is.null(groups.to.test))){
    if(compare.to %in% groups.to.test){
      groups.to.test <- groups.to.test[-which(groups.to.test == compare.to)]
      if(length(groups.to.test) == 0){stop("Error. Main comparison and groups are not properly defined")}
    }
  } else{
    groups.to.test <- unique(as.character(unlist(REA_Results.clean$EnrichmentCategory)))[-which(unique(as.character(unlist(REA_Results.clean$EnrichmentCategory))) == compare.to)]
    groups.to.test <- groups.to.test[order(groups.to.test, decreasing = TRUE)]
  }

  if(all(groups.to.test %in% REA_Results.clean$EnrichmentCategory)){
    REA_Results.clean <- REA_Results.clean[which(REA_Results.clean$EnrichmentCategory %in% append(groups.to.test, compare.to)), ]
  } else{
    stop("Error. Main comparison and groups are not properly defined")
  } 
  
  REA_Results.clean$EnrichmentCategory <- as.character(REA_Results.clean$EnrichmentCategory)
  
  if(any(as.numeric(table(REA_Results.clean$EnrichmentCategory)) < representation.cutoff)){
    remove_group.s <- names(table(REA_Results.clean$EnrichmentCategory))[which(as.numeric(table(REA_Results.clean$EnrichmentCategory)) < representation.cutoff)]
    if(length(remove_group.s) > 1){
      for(d in remove_group.s){
        if(d == remove_group.s[1]){print("The following groups do not have sufficient cells for statistcial analysis. Removing these groups:")}
        print(d)}
    } else{
      cat(paste("There are not sufficient cells in the following group for statistcial analysis. Removing this group:", "\n",remove_group.s))
      print(remove_group.s)
    }
    
    REA_Results.clean <- REA_Results.clean[-which(REA_Results.clean$EnrichmentCategory %in% remove_group.s),]
  }
  
  #~~# Calculate significantly changed markers via SAMR  https://www.rdocumentation.org/packages/samr/versions/3.0/topics/samr
  
  all_stat.Results <- data.frame()
  for(group2 in groups.to.test){
    res <- compare.markers.between.neighborhoods(dat = REA_Results.clean, group1 = compare.to, group2, proteins = markers, arcsinh_transformed)
    res$comparison <- paste0(compare.to,".vs.",group2); res <- res[,c(ncol(res), c(1:(ncol(res)-1)))]
    all_stat.Results <- rbind(all_stat.Results,res)
  }
      all_stat.Results$SAMqval <- as.numeric(all_stat.Results$SAMqval) / 100
      colnames(all_stat.Results)[which(colnames(all_stat.Results) == "SAMqval")] = "SAMpvalue"
      write.csv(all_stat.Results, file = paste0(save.space,"/Cluster-MicroSubset-Marker-Analysis.csv"), row.names = F)
      all_stat.Results$Sig <- "p>0.05"; all_stat.Results$Sig[which(as.numeric(as.character(all_stat.Results$SAMpvalue)) < 0.05)] <- "Significant"
  
  
  
  #~~# Plotting results
        
  ## histogram plots      
  hist_dat <- melt(REA_Results.clean)
      allLevels <- append(compare.to, groups.to.test)
  
        tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", 
                        "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", 
                        "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
        
        color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
        color = color[-which(grepl("white", color))]
        
        if(!(length(tol21rainbow) >= length(allLevels))){
          enrich_colors = sample(color, length(allLevels))
        } else{enrich_colors = tol21rainbow[1:length(allLevels)]}
        enrich_colors[1] = "lightgray"
  
    hist_dat$EnrichmentCategory <- factor(hist_dat$EnrichmentCategory, levels = allLevels)
    hist_dat$Sig <- "base"
      plot_index <- unique(paste(hist_dat$EnrichmentCategory, hist_dat$variable, sep = "_")); plot_index <- plot_index[-which(grepl(compare.to, plot_index))]
        sig_res_index <- paste(gsub(paste0(compare.to,".vs."),"",all_stat.Results$comparison),all_stat.Results$Protein,sep="_")
      
  for(j in plot_index){hist_dat$Sig[which(paste(hist_dat$EnrichmentCategory, hist_dat$variable, sep = "_") == j)] = all_stat.Results$Sig[which(sig_res_index == j)]}
        
  hist_dat2 <- hist_dat[-which(hist_dat$Sig == "p>0.05"),] 
  hist_dat2$fill_it <- NA; hist_dat2$fill_it[which(hist_dat2$Sig == "base")]<- "baseline"
  
  if(arcsinh_transformed){
    plot <- ggplot(hist_dat2, aes(x = value)) + 
      geom_density(aes(group = EnrichmentCategory, color = EnrichmentCategory, fill = fill_it, y = ..scaled..), size =0.7) + 
      scale_x_continuous(trans = 'asinh',limits = c((min(hist_dat2$value)-1),max(hist_dat2$value))) +
      xlab("arcsinh(Signal Intensity)") + ylab("counts") + facet_wrap(~ variable) + theme_classic() + ylim(c(0,1.1))+
      theme(axis.text.x = element_text(angle = 45)) +
      scale_color_manual(values = enrich_colors) +
      scale_fill_manual(values = enrich_colors) + 
      ggtitle(saveName)
  } else {
    hist_dat2$value <- asinh(hist_dat2$value/5)
    plot <- ggplot(hist_dat2, aes(x = value)) + 
      geom_density(aes(group = EnrichmentCategory, color = EnrichmentCategory, fill = fill_it, y = ..scaled..), size =0.7) + 
      scale_x_continuous(trans = 'asinh',limits = c((min(hist_dat2$value)-1),max(hist_dat2$value))) +
      xlab("arcsinh(Signal Intensity)") + ylab("counts") + facet_wrap(~ variable) + theme_classic() + ylim(c(0,1.1))+
      theme(axis.text.x = element_text(angle = 45)) + 
      scale_color_manual(values = enrich_colors) +
      scale_fill_manual(values = enrich_colors) + 
      ggtitle(saveName)
  }
  
  
  ggsave(plot, filename = paste0(save.space,"/Regional-Enrichment_Marker.Histograms.pdf"), width = 15, height = 10, useDingbats = FALSE)

  
  ## heatmap plot
  heat_dat <- REA_Results.clean
      #generate Z-Score
      for(m in markers){
        dat <- heat_dat[,which(colnames(heat_dat) == m)]
        dat <- (dat - mean(dat)) / sd(dat)
        heat_dat[,which(colnames(heat_dat) == m)] <- dat
      }
  m <- melt(heat_dat)
  m$sst <- paste(m$EnrichmentCategory, m$variable, sep = "_")
      m2 <- aggregate(value ~ sst, data = m, mean)
      m2$EnrichmentCategory <- gsub("\\_.*$","", m2$sst)
          m2$EnrichmentCategory <- factor(m2$EnrichmentCategory, levels = rev(allLevels))
      m2$variable <- gsub(".*\\_","", m2$sst)
          m2$variable <- factor(m2$variable, levels = markers)
      m2$stat <- NA
        for(i in 1:nrow(m2)){
          if(m2$EnrichmentCategory[i] != compare.to){
            index <- which(as.character(all_stat.Results$Protein) == m2$variable[i] & 
                     m2$EnrichmentCategory[i] == gsub(paste0(compare.to,".vs."),"", as.character(all_stat.Results$comparison)))
            m2$stat[i] = all_stat.Results$SAMpvalue[index]
          }
        }
      m2$signif <- "yes"; m2$signif[which(m2$stat > 0.05)] = NA; m2$signif[which(is.na(m2$stat))] = NA
  
  map <- ggplot(m2, aes(x = variable, y = EnrichmentCategory, fill = value)) + geom_tile() + 
    geom_point(aes(shape =signif), color = "black", fill = "white", size = 1.5) +
    scale_fill_viridis(limits = c(-0.5,0.5), oob=squish, option = "C") + theme_minimal() + 
    ggtitle("") + scale_shape_manual(values = c(23)) + xlab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ggtitle(saveName)
  map
  ggsave(map, filename = paste0(save.space,"/Regional-Enrichment_Marker.Zscore.Heatmap.pdf"),
         height = 4, width = 7, useDingbats = FALSE)
}



# Transform data
asinh_trans = function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x/5), 
            inverse = function(x) (sinh(x) * 5))
}







