# Regional_Enrichment_Analysis
Using KNN to define intercluster neighborhoods within clustered single cell datasets

This function is designed to identify unique cell states affilitated with different conditions (applied to different monotherapy or combination immunotherapies in Allen et al 2021, in preparation). It expands upon unsupervised clustering, which imperfectly captures nuanced diffrences in protein expression intensity. REA uses nearets neighbor analysis to identify regions in high-dimensional space that are significantly enriched for a particular condition, allowing for more ganular assessment of unique phenotypes within clusters. Once sub-cluster regions are identified, you can statistcially compare protein expresssion intensities betweeen unenriched and enriched zones. REA is designed to follow unsupervised clustering of user choice (we use Clara clustering) and is optimized for high-dimensional protein expression data. 


## Install required R packages
You will need to install the following packages, using install.packages("..."):  
 * FNN
 * uwot
 * prodlim
 * nnet
 * dyplr
 * pheartmap
 * ggplot2
 * reshape2
 * scales
 * virdis
 * samr


## Install Regional Enrichment Analysis

Once you have successfully installed the packages above, download the code file:  Regional Enrichment Analysis.R  
Open an R session and load the functions into your working directory by running the entire script. You should now have three functions loaded:   
```
regionalEnrichmentAnalysis()  
REA_differentiating_Markers() 
asinh_trans()  
```

## Usage
This is the main function to identify enriched high-dimensional regions.  
```
regionalEnrichmentAnalysis(singleCellData, 
                           saveDirectory, 
                           saveName, 
                           markers,
                           conditions = NULL, #best to specify
                           cell_subset_number = NA,
                           by_min_sampleN = FALSE,
                           groupSize = 15,
                           return_result = TRUE,
                           
                        #Output visuals
                           EuclideanDistance = TRUE,
                           printUMAP_cluster = TRUE,
                           printUMAP_condition = TRUE,
                           printUMAP_markers = TRUE,
                           print_NeighborhoodHeatMap = FALSE,
                           print_NeighborhoodUMAP = TRUE)  
```
**Sample Information**  
 - **singleCellData**: This is tha main data matrix. It must contain individual cells by rows and proteins by column. The protein columns must be followed by three identifier colums (cluster: pre-determined cluster number, sampleID: unique sample idenfiter, condition: treatment or condition status).

Protein 1-N... | cluster | sampleID | condition
------------- | ------------- | ------------- | -------------
Expression Value  | 1 | mouseA | Control
Expression Value  | 2 | mouseB | Untreat
Expression Value  | 2 | mouseC | Treatment1

 - **saveDirectory**: Define the directory for all saved outputs
 - **saveName**: Define the name of a new folder, generated in saveDirectory
 - markers: If NULL, all columns preceeding "cluster" will be used as proteins. Otherwise, input a character vector of protein columns to include. eg markers = c("CD3","CD4","Ki67")
 - **conditions**: If you want conditions to appear graphically in a specific order, input here. It's best to specify and place the control or untreated group first.
 - **cell_subset_number**: Compositioanl. To subset the dataframe by condition before analysis, input number of cells per condition. Subset is random across samples in each conidtion.
 - **by_min_sampleN**: (TRUE/FALSE) Abundance. If you would rather visualize relative number of cells per condition, subset instead to control for the N in each group by selecting TRUE.
 - **groupSize**: K, or the numebr of nearest neighbors to be analyzed. The default is 15.
 - **return_result**: (TRUE/FALSE) To return the resulting dataframe as an object, select TRUE and set regionalEnrichmentAnalysis() to a variable. (eg. output = regionalEnrichmentAnalysis(...) )
 
 **Output visuals**  
 - **EuclideanDistance**: (TRUE/FALSE) Calculates euclidean distance between all cells between conditions and saves plot. Stats caluclated by Wilxocon, but less useful if condition order is not specified above.
 - **printUMAP_cluster**: (TRUE/FALSE) Saves UMAP plot of cells colored by cluster.
 - **printUMAP_condition**: (TRUE/FALSE) Saves UMAP plot of cells colored by conidtion.
 - **printUMAP_markers**: (TRUE/FALSE) Creates sub-directory, and saves UMAP plots colored by intensity of each protein.
 - **print_NeighborhoodHeatMap**: (TRUE/FALSE) Saves heatmap of neighbhoord confirguations called for each enrichment category.
 - **print_NeighborhoodUMAP**: (TRUE/FALSE) Saves UMAP plot of cells colored by Regional Enrichment Analysis.

### Notes: 
 - We perform the analysis on arcsinh transformed values: asinh(X/5)
 - Please **do not** use underscores " _ " in the cluster, sampleID, or condition columns. Replace with "." if needed.
 - If you do not wish to cluster your data beforehand, you can set the cluster column to NA. Make sure to set printUMAP_cluster = FALSE.
 - The function will pause after calculating nearest neighbors and ask if you want to specific the enrichment catagory names. Otherwise it will assign generitc Group1... GroupN. If your specified names all match the conditions, the same color scheme will be used.


## Differentiating Markers, Statistical Analyses
This function identifies the markers distinguishing enriched regions.
Significance Anlaysis of Microarrays () is used to identify significant differences in protein expression intensities, p* < 0.05. The outputs include: 1) sumary csv of median signal intensities for each group, p value, and fold change, 2) histograms for each enrichment group, and 3) a heatmap of protein expression Z-scores across enrichment groups.
```
REA_differentiating_Markers(REA_Results,
                            saveDirectory,
                            saveName, 
                            arcsinh_transformed = TRUE,
                            cluster.of.interest = NULL,
                            markers = NULL,
                            groups.to.test = NULL,
                            compare.to = NULL,
                            representation.cutoff = 100)  
```
- **REA_Results**: The Regional-Enrichment-Analysis_Results output from the above function.
- **saveDirectory**: Directory to save results. Can be same as above.
- **saveName**: Sub-folder name specifying the chosen comparison.
- **arcsinh_transformed**: (TRUE/FALSE) TRUE if you are working with arcsinh transformed data, for histogram graphing.
- **cluster.of.interest**: Specify the cluster for regional comparison, corresponding to your cluster column.
- **markers**: Specify which markers to include in the statistical analysis. Default is NULL, comparing all markers.
- **groups.to.test**: Specify which enrichment groups you wish to compare. Default is NULL, comparing all enrichment groups. 
- **compare.to**: Specify which group to compare every other group. Default is 'Unenriched'. 
- **representation.cutoff**: Note that clusters may not have sufficient cells enriched for every group. Poorly represented groups will be removed by this cutoff, default 100 cells. 
