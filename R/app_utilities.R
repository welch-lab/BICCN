library(dplyr)
library(tidyr)
library(ggplot2)
library(rliger)
library(reshape)

#' @param add_genes provides a space for the user to upload a list of personalized genes
#' Run this function to generate the data which the MarkerGeneBrowser Shiny App can run
#'
#' Example:
#'
#' sus = markergenes(region = "RSP", datasets = c("atac", "sc10Xv2", "sc10Xv3", "sn10Xv3", "smartseq"), analysis_num = 1, filepath = "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region/")
# 'Function to generate the subset normalized/scaled data for each analysis for the viewer
generate_markersobject = function(datasets, region, analysis_num, filepath){
  #Get relevent genes
  genes = readRDS("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Base_Reference_Files/Marker_genes.RDS")
  genes_oi = unique(c(genes$M1, genes$M2, genes$M3, genes$M4, genes$M5, genes$M6))
  #Build a list of your datasets
  dataset_names = c()
  dataset_matrices = c()
  for( i in 1:length(datasets)){
    name = paste0(region, "_", datasets[[i]])
    dataset_names = c(dataset_names, name)
    matrix_filename = paste0(filepath, region, "/Analysis", analysis_num, "_", region, "/", region, "_", datasets[[i]], "_qc.RDS")
    matrix = readRDS(matrix_filename)
    ligs = createLiger(list(newdata = matrix))
    ligs = normalize(ligs)
    ligs = selectGenes(ligs, var.thresh = 0.0)
    genes_needed = subset(ligs@var.genes, ligs@var.genes %in% genes_oi)
    ligs@var.genes = genes_needed
    ligs = scaleNotCenter(ligs)
    dataset_matrices[[i]] = t(ligs@scale.data$newdata)
  }
  names(dataset_matrices) = dataset_names
  output_file = paste0(filepath, region, "/Analysis", analysis_num, "_", region, "/Images/GeneExpression_", region,"_Analysis", analysis_num, ".RDS")
  saveRDS(dataset_matrices, output_file)
}


#Functions used
getMarkerGenes = function(object, verbose = TRUE, chunk = 1000, marker_genes, slot.use = "norm.data"){
  marker_gene_matrix = c()
  if (class(object@raw.data[[1]])[1] == "H5File") {
    # if (verbose) {
    #  message("Start sampling")
    #}
    hdf5_files = names(object@raw.data)
    for (i in 1:length(hdf5_files)) {
      #if (verbose) {
      # message(hdf5_files[i])
      #}
      if (object@h5file.info[[i]][["format.type"]] == "AnnData"){
        genes_i = object@h5file.info[[i]][["genes"]][]$index
      } else {
        genes_i = object@h5file.info[[i]][["genes"]][]}
      genes.use = intersect(genes_i,marker_genes)
      data.subset = Matrix(nrow=length(genes.use),ncol=0,sparse=TRUE)
      chunk_size = chunk
      if (object@h5file.info[[i]][["format.type"]] == "AnnData"){
        barcodes = object@h5file.info[[i]][["barcodes"]][]$cell
        genes = object@h5file.info[[i]][["genes"]][]$index
      } else {
        barcodes = object@h5file.info[[i]][["barcodes"]][]
        genes = object@h5file.info[[i]][["genes"]][]
      }
      num_cells = length(barcodes)
      num_genes = length(genes)
      prev_end_col = 1
      prev_end_data = 1
      prev_end_ind = 0
      num_chunks = ceiling(num_cells/chunk_size)
      if (verbose) {
        pb = txtProgressBar(0, num_chunks, style = 3)
      }
      ind = 0

      while (prev_end_col < num_cells) {
        ind = ind + 1
        if (num_cells - prev_end_col < chunk_size) {
          chunk_size = num_cells - prev_end_col + 1
        }
        if (any(grepl("meth", hdf5_files[i]))){
          slot.use  = "scale.data"
        } else {slot.use  = "norm.data"}


        if (slot.use != "scale.data"){
          start_inds = object@h5file.info[[i]][["indptr"]][prev_end_col:(prev_end_col+chunk_size)]
          row_inds = object@h5file.info[[i]][["indices"]][(prev_end_ind+1):(tail(start_inds, 1))]
          if (slot.use=="norm.data")
          {
            counts = object@norm.data[[i]][(prev_end_ind+1):(tail(start_inds, 1))]
          }
          one_chunk = sparseMatrix(i=row_inds[1:length(counts)]+1,p=start_inds[1:(chunk_size+1)]-prev_end_ind,x=counts,dims=c(num_genes,chunk_size))
          rownames(one_chunk) = genes
          colnames(one_chunk) = barcodes[(prev_end_col):(prev_end_col+chunk_size-1)]
          use_these = colnames(one_chunk)
          one_chunk = one_chunk[genes.use,use_these]
          data.subset = cbind(data.subset,one_chunk)

          num_read = length(counts)
          prev_end_col = prev_end_col + chunk_size
          prev_end_data = prev_end_data + num_read
          prev_end_ind = tail(start_inds, 1)
          setTxtProgressBar(pb, ind)
        } else {
          one_chunk = object@scale.data[[i]][,prev_end_col:(prev_end_col + chunk_size - 1)]
          rownames(one_chunk) = object@var.genes
          colnames(one_chunk) = barcodes[(prev_end_col):(prev_end_col+chunk_size-1)]
          use_these = colnames(one_chunk)
          one_chunk = one_chunk[genes.use,use_these]
          data.subset = cbind(data.subset,one_chunk)

          prev_end_col = prev_end_col + chunk_size
          if (verbose) {
            setTxtProgressBar(pb, ind)
          }
        }
      }
      marker_gene_matrix[[i]] = data.subset

      if (verbose) {
        setTxtProgressBar(pb, num_chunks)
        cat("\n")
      }
    }
    names(marker_gene_matrix) = hdf5_files
  }
  return(marker_gene_matrix)
}


#Function dotplot_generator
#' @param Dataset, what dataset we want expression profiles from
#' @param markerGenes, What dataset we want to draw our marker genes from
#' @param  resolution, can be low or high, indicates that the clusters are drawn from low or high resolution louvain resolution
#' @param CellSubtype, Broad categories, i.e "Exc", "Inh", or "NonN"
#' This returns the general matrix we interested in
dotplot_generator = function(Dataset, markerGenes, resolution, CellSubtype){
  #Load in the list of marker genes, which has the marker genes from atac, methylation, and rna each saved as a named list
  genes = readRDS("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/MarkerGeneBrowser/Marker_genes.RDS")
  #Filter for your specific genes of interest
  genes_of_interest = filter(genes, genes$Dataset == markerGenes)
  genes_of_interest = filter(genes_of_interest, genes_of_interest$MajorType == CellSubtype)
  m_marker_genes = unique(c(genes_of_interest$M1, genes_of_interest$M2, genes_of_interest$M3, genes_of_interest$M4, genes_of_interest$M5, genes_of_interest$M6))
  #Read in and subset the appropriate datasets
  dataset_oi = dataset.vector[[Dataset]]
  dataset_oi = dataset_oi[subset(rownames(dataset_oi), rownames(dataset_oi) %in% m_marker_genes),]
  dataset_oi = as.matrix(dataset_oi)
  dataset_oi = t(dataset_oi)
  dataset_oi = as.data.frame(dataset_oi)
  dataset_oi$Barcode = rownames(dataset_oi)

  results.oi = filter(results.table, results.table$dataset == Dataset)
  clusts = results.oi[,c(resolution, "Barcode")]
  colnames(clusts) = c("resolution", "Barcode")
  totals = cbind(dataset_oi, clusts)
  totals = subset(totals, select = -Barcode)
  totals = subset(totals, select = -Barcode)
  return(totals)
}
PercentAbove <- function(x, threshold){
  100*length(x = x[x > threshold]) / length(x = x)
}
MinMax <- function(data, min, max) {
  data2 <- data
  data2[data2 > max] <- max
  data2[data2 < min] <- min
  return(data2)
}

#Function DotPlot, takes as input
#' @param object meant to be the object returned by the dotplot_generator
#' @param order meant to be the order for the genes, returned by the heatmap of the cell types (i.e marker_table)
#' @param resolution, high or low louvain resolution. Can be specified as lowRcluster or highRcluster
DotPlot = function(object, resolution, cols.use=c("yellow","red"), scale.min = 0, col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 4, scale.by = "radius", scale.max = NA, plot.legend = TRUE, do.return = TRUE, x.lab.rot = TRUE, threshold = 0.0){
  object = object %>% tidyr::pivot_longer(!resolution, names_to = "Gene_Name", values_to = "Expression")
  scale.func <- switch(EXPR = scale.by, size = scale_size, radius = scale_radius)
  colnames(object) = c("Resolution", "Gene_Name", "Expression")
  data.to.plot <- object %>% group_by(Resolution, Gene_Name) %>% summarize(avg.exp = mean(expm1(x = Expression)), pct.exp = PercentAbove(x = Expression, threshold))
  data.to.plot <- data.to.plot %>% ungroup() %>% group_by(Gene_Name) %>% mutate(avg.exp.scale = scale(x = avg.exp))  %>% mutate(avg.exp.scale = MinMax(data = avg.exp.scale,max = col.max, min = col.min))
  data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
  #Reorder the dotplot using a heatmap
  object2 = data.to.plot[, c("Resolution", "Gene_Name", "avg.exp.scale")]
  object2 = data.frame(object2 %>% pivot_wider(names_from = Gene_Name, values_from = avg.exp.scale ))
  rownames(object2) = paste0("Cluster_", object2$Resolution)
  object2 = object2[,-1]
  object2 = as.matrix(object2)
  object2 = t(object2)
  heat = heatmap.2(object2)
  cluster_order = rownames(heat$carpet)
  gene_order = colnames(heat$carpet)
  ###########################################################
  data.to.plot$Gene_Name = factor(data.to.plot$Gene_Name, levels=gene_order) #Reorders the genes to match the cell type markers heatmap
  data.to.plot$Cluster_Name = paste0("Cluster_", data.to.plot$Resolution)

  data.to.plot$Cluster_Name = factor(data.to.plot$Cluster_Name, levels= rev(cluster_order))
  p = ggplot(data = data.to.plot, mapping = aes(x = Cluster_Name, y = Gene_Name )) + geom_point(mapping = aes(size = pct.exp, color = avg.exp.scale)) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text = element_text(size = 20) )

  p <- if (length(x = cols.use) == 1) {
    p <- p + scale_color_distiller(palette = cols.use)
  } else {
    p <- p + scale_color_gradient(low = cols.use[1], high = cols.use[2])
  }
  if (!plot.legend) {
    p <- p + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90,
                                              vjust = 0.5))
  }
  return(p)

}
#marker_table is meant to return a heatmap of genes by cell types, with corresponding entries highlighted to allow easy identification of relationships between marker genes and cell types
#' @param obj meant to be the object returned by dotplot_generator
#' @param markerGenes the dataset by which we want to draw our marker genes from
#' @param CellSubtype the type of broad cell classification we currently wish to inspect

marker_table = function(obj, markerGenes, CellSubtype, resolution, scale.min = 0, col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 4, scale.by = "radius",  threshold = 0.0){
  ##### Get Gene Orderings
  object = obj %>% tidyr::pivot_longer(!resolution, names_to = "Gene_Name", values_to = "Expression")
  scale.func <- switch(EXPR = scale.by, size = scale_size, radius = scale_radius)
  colnames(object) = c("Resolution", "Gene_Name", "Expression")
  data.to.plot <- object %>% group_by(Resolution, Gene_Name) %>% summarize(avg.exp = mean(expm1(x = Expression)), pct.exp = PercentAbove(x = Expression, threshold))
  data.to.plot <- data.to.plot %>% ungroup() %>% group_by(Gene_Name) %>% mutate(avg.exp.scale = scale(x = avg.exp))  %>% mutate(avg.exp.scale = MinMax(data = avg.exp.scale,max = col.max, min = col.min))
  #Reorder the dotplot using a heatmap
  obj2 = data.to.plot[, c("Resolution", "Gene_Name", "avg.exp.scale")]
  obj2 = data.frame(obj2 %>% pivot_wider(names_from = Gene_Name, values_from = avg.exp.scale ))
  rownames(obj2) = paste0("Cluster_", obj2$Resolution)
  obj2 = obj2[,-1]
  obj2 = as.matrix(obj2)
  obj2 = t(obj2)
  heat = heatmap.2(obj2)
  gene_order = colnames(heat$carpet)
  ###########################################################
  cumulative = readRDS("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/MarkerGeneBrowser/Cumulative_marker_genes.RDS")
  inputmatrix = cumulative[[markerGenes]]
  inputmatrix = filter(inputmatrix, inputmatrix$MajorType == CellSubtype)
  reduce_genes = subset(colnames(inputmatrix), colnames(inputmatrix) %in% colnames(obj))
  reduce_genes = c(reduce_genes, "SubType")
  inputmatrix = inputmatrix[,reduce_genes]
  inputmatrix <- inputmatrix %>% pivot_longer(!SubType, values_to = "Relationship")
  inputmatrix$name = factor(inputmatrix$name, levels = gene_order)
  plot_markers =    ggplot(inputmatrix, aes(SubType, name, fill = Relationship)) +
    geom_tile() + theme(axis.text.x = element_text(angle = 90,
                                                   vjust = 0.5)) + xlab("Cell Type") + ylab("Marker Gene") + theme(text = element_text(size = 20))
  return(plot_markers)
}

graph_celltypes = function(resolution, cluster_picked){
  datatypes_present = unique(results.table$dataset)
  if(resolution == "lowRcluster"){
    results.table2 = results.table[,c("dataset", "lowRcluster", "OG_Annotations")]
  } else {
    results.table2 = results.table[,c("dataset", "highRcluster", "OG_Annotations")]
  }
  colnames(results.table2) = c("dataset", "Resolution", "OG_Annotations")

  results.table2 = filter(results.table2, results.table2$Resolution == cluster_picked)
  unos = data.frame(results.table2 %>% group_by(dataset, OG_Annotations) %>% tally())
  final_plot = ggplot(data = unos, mapping = aes(x=OG_Annotations, y=n, fill=OG_Annotations)) +
    geom_bar(stat="identity", position = "dodge") +
    geom_text(aes(label = n), vjust = -0.2, size = 4, position = position_dodge(0.9))+ facet_wrap(facets = vars(dataset))  + labs(x ="Original Annotations", y = "Number of Cells") +  guides(fill=guide_legend(title="Original Annotation")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(text = element_text(size = 20))
  return(final_plot)
}

#`filepath` = "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Sample_Validation/Automated/"
#`analysis_num` = 3
#`region` = "AUD"


#' The generate_umaps function generates labeled umaps given our newly derived LIGER annotations
#' 2nd: It generates a READme with all of the region/analysis umaps in one place
generate_umaps =  function(filepath, analysis_num, region){
  liger.path = paste0(filepath, region, "/Analysis", analysis_num, "_", region, "/onlineINMF_", region, "_object.RDS")
  ligs = readRDS(liger.path)
  results.path = paste0(filepath, region, "/Analysis", analysis_num, "_", region, "/Analysis", analysis_num, "_", region, "_Results_Table.RDS")
  results = readRDS(results.path)
  ligs_clusts = names(ligs@clusters)
  results = results[match(ligs_clusts, results$Barcode), ]
  low_re = results$lowRAnnotations
  names(low_re) = results$Barcode
  high_re = results$highRAnnotations
  names(high_re) = results$Barcode

  ligs_low  = ligs
  ligs_low@clusters = as.factor(low_re)
  ligs_high = ligs
  ligs_high@clusters = as.factor(high_re)

  umaps = data.frame(results$UMAP1, results$UMAP2)
  rownames(umaps) = results$Barcode
  umaps = as.matrix(umaps)

  ligs_low@tsne.coords = umaps
  plot.labeled.low = plotByDatasetAndCluster(ligs_low,return.plots = TRUE, text.size = 6 )
  low_umap2 =paste0(filepath, "/", region, "/Analysis", analysis_num, "_", region, "/Images/Umap2_", region, "_Analysis_", analysis_num, "ligerlabeled.lowresolution.pdf")
  low_umap2png =paste0(filepath, "/", region, "/Analysis", analysis_num, "_", region, "/Images/Umap2_", region, "_Analysis_", analysis_num, "ligerlabeled.lowresolution.png")

  png(low_umap2png, 1000, 800)
  print(plot.labeled.low[[2]])
  dev.off()

  pdf(low_umap2, width = 10, height = 8)
  print(plot.labeled.low[[2]])
  dev.off()

  ligs_high@tsne.coords = umaps
  plot.labeled.high = plotByDatasetAndCluster(ligs_high,return.plots = TRUE, text.size = 6)
  high_umap2 =paste0(filepath, "/", region, "/Analysis", analysis_num, "_", region, "/Images/Umap2_", region, "_Analysis_", analysis_num, "ligerlabeled.highresolution.pdf")
  high_umap2png =paste0(filepath, "/", region, "/Analysis", analysis_num, "_", region, "/Images/Umap2_", region, "_Analysis_", analysis_num, "ligerlabeled.highresolution.png")

  png(high_umap2png, 1000, 800)
  print(plot.labeled.high[[2]])
  dev.off()

  pdf(high_umap2 , width = 10, height = 8)
  print(plot.labeled.high[[2]])
  dev.off()
  ###### Plot the png and pdf
  print("Done! New graphs have been generated and saved!")

}
