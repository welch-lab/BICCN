# Functions required for utilization of the DUMfound deconvolution pipeline on
# fully annotated region atlases

#' @param filepath Path to directory within which the atlas structure was generated
#' @param region The anatomical region with which the analysis is associated
#' @param spatial.data.file A filepath to the spatial expression matrix
#' @param coords An object or filepath to coordinates for the spatial samples
#' @param spatial.data.name A string, the name of the spatial dataset

save_spatial_data = function(filepath,
                        region,
                        spatial.data.file,
                        coords,
                        spatial.data.name
                        ){
  dir_new = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  if(!dir.exists(dir_new)){
    dir.create(dir_new)
    message("Created directory at ", dir_new)
  }
  file.copy(spatial.data.file, paste0(dir_new,"/",spatial.data.name,"_exp.RDS"))
  if(is.character(coords)){
    file.copy(coords, paste0(dir_new,"/",spatial.data.name,"_coords.RDS"))
  } else {
    saveRDS(coords, paste0(dir_new,"/",spatial.data.name,"_coords.RDS"))
  }
}

#' Generate plots of the provided spatial data on the XY, YZ, and XZ planes in
#' a reference grid to allow for indexing.
#'
#' @param filepath Path to directory within which the atlas structure was generated
#' @param region A string corresponding to the name of an anatomical region
#' @param spatial.data.name A string, the name of the spatial dataset
#' @param display.plots A logical, corresponding with if to display requested plots upon
#'    generation
#'
#' @return nothing
#'
#' @import ggplot2, cowplot, viridis
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }
#'

reference_3d_coordinates = function(filepath,
                                    region,
                                    spatial.data.name,
                                    save.plots = FALSE){

  library(ggplot2)
  library(cowplot)
  coords = as.data.frame(readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name,"/",spatial.data.name,"_coords.RDS")))

  minmax = apply(coords, MARGIN = 2, function(x){range(x)})
  ranges = apply(minmax, MARGIN = 2, function(y){return(y[2]-y[1])})

  p1 = ggplot(coords, aes_string(x = colnames(coords)[1], y = colnames(coords)[2]))+
    geom_tile() +
    coord_fixed(ratio = 1)+
    geom_tile(aes(alpha = .1)) +
    coord_fixed(ratio = 1) +
    theme_minimal_grid() +
    #scale_y_continuous(minor_breaks = seq(minmax[1,2] , minmax[2,2], 1), breaks = seq(minmax[1,2] , minmax[2,2], 5)) +
    #scale_x_continuous(minor_breaks = seq(minmax[1,1] , minmax[2,1], 1), breaks = seq(minmax[1,1] , minmax[2,1], 5)) +
    xlab(colnames(coords)[1]) +
    ylab(colnames(coords)[2]) +
    ggtitle(paste0("X-Y silhouette of ",region))+
    theme(legend.position="none",
          text = ggplot2::element_text(size = 5),
          axis.text = ggplot2::element_text(size = 4),
          plot.title = element_text(size=6))
  p2 = ggplot(coords, aes_string(x = colnames(coords)[1], y = colnames(coords)[3])) +
    geom_tile() +
    coord_fixed(ratio = 1)+
    geom_tile(aes(alpha = .1)) +
    coord_fixed(ratio = 1) +
    theme_minimal_grid() +
    #scale_y_continuous(minor_breaks = seq(minmax[1,3] , minmax[2,3], 1), breaks = seq(minmax[1,3] , minmax[2,3], 5)) +
    #scale_x_continuous(minor_breaks = seq(minmax[1,1] , minmax[2,1], 1), breaks = seq(minmax[1,1] , minmax[2,1], 5)) +
    xlab(colnames(coords)[1]) +
    ylab(colnames(coords)[3]) +
    ggtitle(paste0("X-Z silhouette of ",region))+
    theme(legend.position="none",
          text = ggplot2::element_text(size = 5),
          axis.text = ggplot2::element_text(size = 4),
          plot.title = element_text(size=6))
  p3 = ggplot(coords, aes_string(x = colnames(coords)[2], y = colnames(coords)[3])) +
    geom_tile() +
    coord_fixed(ratio = 1)+
    geom_tile(aes(alpha = .1)) +
    coord_fixed(ratio = 1) +
    theme_minimal_grid() +
    #scale_y_continuous(minor_breaks = seq(minmax[1,3] , minmax[2,3], 1), breaks = seq(minmax[1,3] , minmax[2,3], 5)) +
    #scale_x_continuous(minor_breaks = seq(minmax[1,2] , minmax[2,2], 1), breaks = seq(minmax[1,2] , minmax[2,2], 5)) +
    xlab(colnames(coords)[2]) +
    ylab(colnames(coords)[3]) +
    ggtitle(paste0("Y-Z silhouette of ",region))+
    theme(legend.position="none",
          text = ggplot2::element_text(size = 5),
          axis.text = ggplot2::element_text(size = 4),
          plot.title = element_text(size=6))
  print(p1)
  print(p2)
  print(p3)
  message("Plots generated")
  if(save.plots){
    plots_dir = paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots")
    if(!dir.exists(plots_dir)){
      dir.create(plots_dir)
      message("Created directory at ", plots_dir)
    }
    ggplot2::ggsave(paste0(plots_dir, "/x_y_coord_reference.PNG"),
                    p1,
                    width = 500,
                    height = 500,
                    units = "px")
    ggplot2::ggsave(paste0(plots_dir, "/x_z_coord_reference.PNG"),
                    p2,
                    width = 500,
                    height = 500,
                    units = "px")
    ggplot2::ggsave(paste0(plots_dir, "/y_z_coord_reference.PNG"),
                    p3,
                    width = 500,
                    height = 500,
                    units = "px")
    message("Plots saved")
  }
}


#' Subset a spatial dataset by coordinates for analysis
#'
#' @param filepath Path to directory within which the atlas structure was generated
#' @param region The anatomical region with which the analysis is associated
#' @param spatial.data.name A string, the name of the spatial dataset
#' @param subset.specs A list with length equal to the number of axises, with
#'    each entry a vector of length two, with the first element being the
#'    minimum value to include and the second being the maximum, or NaN to
#'    indicate a missing value
#' @param out.filepath Path to directory to save subset data to, if NULL the
#'    expression and coordinates in the atlas directory structure are
#'    overwritten
subset_spatial_data = function(filepath,
                               region,
                               spatial.data.name,
                               subset.specs = list(c(NaN, NaN),
                                                   c(NaN, NaN),
                                                   c(NaN, NaN)),
                               out.filepath = NULL){
  deconv_dir = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/")
  coords = readRDS(paste0(deconv_dir,spatial.data.name,"/",spatial.data.name,"_coords.RDS"))
  spatial.data = readRDS(paste0(deconv_dir,spatial.data.name,"/",spatial.data.name,"_exp.RDS"))
  for(i in 1:ncol(coords)){
    subset.specs[[i]][is.nan(subset.specs[[i]])] = range(coords[,i])[is.nan(subset.specs[[i]])]
    coords = coords[coords[,i] >= subset.specs[[i]][1] & coords[,i] <= subset.specs[[i]][2], ]
  }
  if(nrow(coords) == 0){
    stop("No samples selected -- provide a new range of coordinates.")
  } else {
    message(paste0("Sample subset to ", nrow(coords), " samples."))
  }
  spatial.data = spatial.data[, rownames(coords)]
  new_spatial.data.name = paste0(spatial.data.name,"_subset_",nrow(coords))
  if(!is.null(out.filepath)){
    saveRDS(spatial.data, paste0(out.filepath, "/", new_spatial.data.name,"_exp.RDS"))
    saveRDS(coords, paste0(out.filepath, "/", new_spatial.data.name,"_coords.RDS"))
    message("Saved expression and coordinates to ", out.filepath)
  } else {
    new_dir = paste0(deconv_dir,new_spatial.data.name)
    dir.create(new_dir)
    message("Created directory at ", new_dir)
    saveRDS(spatial.data, paste0(new_dir,"/",new_spatial.data.name,"_exp.RDS"))
    saveRDS(coords, paste0(new_dir,"/",new_spatial.data.name,"_coords.RDS"))
    message("Saved expression and coordinates to ", new_dir)
  }
}

sample_single_cell = function(
    filepath,
    region,
    n.cells = 500,
    known.annotations = NULL,
    naive.clusters = FALSE,
    naive.clusters.remove = NULL,
    rand.seed = 123
){
  
  set.seed(rand.seed)
  
  message("Loading Data")
  
  object_paths = paste0(filepath,"/", region, "/Analysis",c(2,4,5),"_", region, "/onlineINMF_",region, "_object.RDS" )
  objects = lapply(object_paths, function(object_path){readRDS(object_path)})
  
  h5_files = Reduce(c, lapply(objects, function(object){sapply(1:length(object@norm.data), function(i){object@h5file.info[[i]]$file.path})}))
  rna_files = grep(paste0("_(sc10Xv3_|smartseq_|sn10Xv3_|sc10Xv2_)"), h5_files, value = TRUE)
  
  liger_cells = lapply(rna_files, function(i){
    rhdf5::h5read(i, "/matrix/barcodes")#change, extract from H5
  })
  
  liger_genes = lapply(rna_files, function(i){
    rhdf5::h5read(i, "/matrix/features")[[1]] #change, extract from H5
  })
  
  
  
  
  if(is.null(known.annotations)){
    clusters = c()
    for(analysis_num in c(2,4,5)){
      analysis_results = readRDS(paste0(filepath,"/",  region, "/Analysis", analysis_num, "_", region, "/Analysis", analysis_num, "_", region,"_Results_Table.RDS"))
      if("highRcluster" %in% colnames(analysis_results)){
        clust.use = "highRcluster"
      } else {
        clust.use = "lowRcluster"
      }
      if(naive.clusters){
        analysis_clusters = paste0(analysis_num ,"_",analysis_results[,clust.use])
        analysis_clusters[analysis_clusters %in% paste0(analysis_num, "_", naive.clusters.remove[[as.character(analysis_num)]])] = ""
      } else {
        analysis_clusters = as.character(analysis_results$highRAnnotations)
      }
      names(analysis_clusters) = analysis_results$Barcode
      clusters = c(clusters, analysis_clusters)
    }
    clusters = clusters[clusters != ""]
    clusters = as.factor(clusters)
  } else {
    clusters = known.annotations
  }
  

  annotated_cells = intersect(names(clusters),Reduce(union, liger_cells))
  
  clusters = clusters[names(clusters) %in% annotated_cells]
  clusters = droplevels(clusters)
  
  
  
  freq_cells = table(clusters)
  freq_cells = freq_cells[names(freq_cells) != ""]
  
  
  message("Sampling from clusters")
  
  liger_cells_combined = Reduce(c, liger_cells)
  
  sample.cells = Reduce(c, lapply(names(freq_cells), function(cell_type){
    cells = intersect(names(clusters[clusters == cell_type]), liger_cells_combined)
    if(length(cells)> 0){
      return(sample(cells, min(length(cells), n.cells), replace =FALSE))
    } else {
      return(c())
    }
  }))
  
  message("Sample summary -- ")
  print(table(clusters[sample.cells]))
  
  descriptor = as.character(rand.seed)
  
  if(is.null(known.annotations)){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
  }
  
  saveRDS(sample.cells, paste0(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/sampled_cells_", descriptor,".RDS")))
  
  if(is.null(known.annotations)){
     saveRDS(clusters, paste0(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/clusters_",descriptor,".RDS")))
  } else {
     saveRDS(clusters, paste0(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/user_defined_clusters_",descriptor,".RDS")))
  }
          
  
  shared_genes = Reduce(intersect, liger_genes)

  norm.data = lapply(1:length(rna_files), function(i){
    n = rna_files[i]
    out_mat = rliger:::Matrix.column_norm(Matrix::sparseMatrix(
      dims = c(length(liger_genes[[i]]), length(liger_cells[[i]])),
      i = as.numeric(rhdf5::h5read(n, "/matrix/indices")+1),
      p = as.numeric(rhdf5::h5read(n, "/matrix/indptr")),
      x = as.numeric(rhdf5::h5read(n, "/matrix/data"))
    ))
    rownames(out_mat) = liger_genes[[i]]
    colnames(out_mat) = liger_cells[[i]]
    out_mat = out_mat[liger_genes[[i]] %in% shared_genes, liger_cells[[i]] %in% sample.cells]
    gene_means = rhdf5::h5read(rna_files[i], "gene_means")[liger_genes[[i]] %in% shared_genes]
    gene_sum_sq = rhdf5::h5read(rna_files[i], "gene_sum_sq")[liger_genes[[i]] %in% shared_genes]
    if(is.vector(out_mat)){
      out_mat = NULL
    } else {
      root_mean_sum_sq = sqrt(gene_sum_sq/(ncol(out_mat)-1))
      out_mat= sweep(out_mat, 1, root_mean_sum_sq, "/") #liger_cells[[i]] %in% sample.cells
      out_mat[is.na(out_mat)] = 0
      out_mat[out_mat == Inf] = 0
    }
    return(out_mat)
  })
  norm.data = norm.data[!sapply(norm.data, function(x){length(x) == 0})]
  saveRDS(norm.data, paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",descriptor,"_norm_data.RDS"))
}

select_defining_genes = function(
    filepath,
    region,
    clusters.from.atlas = TRUE,
    naive.clusters = FALSE,
    deconv.gene.num = 2000,
    gene.num.tol = 50,
    rand.seed = 123
    ){
  set.seed(rand.seed)
    
  descriptor = as.character(rand.seed)
  
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
  }
      
  norm.data = readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",descriptor,"_norm_data.RDS"))
  
  if(clusters.from.atlas){
     clusters = readRDS(paste0(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/clusters_",descriptor,".RDS")))
  } else {
     clusters = readRDS(paste0(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/user_defined_clusters_",descriptor,".RDS")))
  }
  
  message("Selecting genes with the KW test")
  
  chisq_list = list()
  for(i in 1:length(norm.data)){
    chisq_list[[i]] = matrixTests::row_kruskalwallis(as.matrix(norm.data[[i]]),as.vector(clusters[names(clusters) %in%colnames(norm.data[[i]])]))$statistic
    names(chisq_list[[i]]) = rownames(norm.data[[i]])
  }
  shared_genes = Reduce(c, lapply(norm.data, rownames))
  
  
  var_thresh_start = 0.5
  var_thresh_old = var_thresh_start
  high = 1
  low = 0
  
  gene_vec = shared_genes
  for(i in 1:length(chisq_list)){
    chisq_list[[i]][is.na(chisq_list[[i]])] = 0
    gene_vec = intersect(gene_vec, names(chisq_list[[i]][chisq_list[[i]] > quantile(chisq_list[[i]], var_thresh_start)]))
  }
  
  while(abs(length(gene_vec)-deconv.gene.num) > gene.num.tol){
    if(length(gene_vec)>deconv.gene.num){
      var_thresh_new = (var_thresh_old+high)/2
      low = var_thresh_old
    } else {
      var_thresh_new = (var_thresh_old+low)/2
      high = var_thresh_old
    }
    gene_vec = shared_genes
    for(i in 1:length(chisq_list)){
      gene_vec = intersect(gene_vec, names(chisq_list[[i]][chisq_list[[i]] > quantile(chisq_list[[i]], var_thresh_new, na.rm = TRUE)]))
    }
    var_thresh_old = var_thresh_new
  }
  
  message(paste0(length(gene_vec), " genes found with p = ",var_thresh_old))

  saveRDS(list(chisq_vals = chisq_list, genes_used = gene_vec), paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/gene_selection_",descriptor,".RDS"))
  
}

qc_spatial_data = function(
    filepath,
    region,
    spatial.data.name,
    slide.seq = FALSE,
    clusters.from.atlas = TRUE,
    naive.clusters = FALSE,
    z = 1,
    n.umi.thresh = 150,
    rand.seed = 123
){
  set.seed(rand.seed)
  
  descriptor = as.character(rand.seed)
  
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
    
    dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
    dir_spatial_new = paste0(dir_spatial, "_naive")

    dir.create(dir_spatial_new)
    file.copy(paste0(dir_spatial,"/",spatial.data.name,"_exp.RDS"),dir_spatial_new)
    file.copy(paste0(dir_spatial,"/",spatial.data.name,"_coords.RDS"),dir_spatial_new)
    file.rename(from = paste0(dir_spatial_new,"/",spatial.data.name,"_exp.RDS"), to = paste0(dir_spatial_new,"/",spatial.data.name,"_naive_exp.RDS"))
    file.rename(from = paste0(dir_spatial_new,"/",spatial.data.name,"_coords.RDS"), to = paste0(dir_spatial_new,"/",spatial.data.name,"_naive_coords.RDS"))

    spatial.data.name = paste0(spatial.data.name, "_naive")
  }
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  
  
  
  spatial.data = readRDS(paste0(dir_spatial,"/",spatial.data.name,"_exp.RDS"))
  coords = readRDS(paste0(dir_spatial,"/",spatial.data.name,"_coords.RDS"))

  gene_data = readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/gene_selection_",descriptor,".RDS"))
  gene_vec = gene_data[[2]]
  
  
  if(!slide.seq){
    spatial.data[spatial.data == -1] = NA
    genes_NA = apply(spatial.data, MARGIN = 1, function(x){sum(is.na(x))})
    mean_genes_NA = mean(genes_NA)
    genes_use = rownames(spatial.data)[genes_NA < (mean_genes_NA + z * mean_genes_NA)]
    message(length(genes_use), " genes saved out of ", nrow(spatial.data), " total.")
    spatial.data = spatial.data[rownames(spatial.data) %in% intersect(genes_use, gene_vec),]
    spatial.data[is.na(spatial.data)] = 0
  } else {
    original_dim = dim(spatial.data)
    spatial.data = spatial.data[rownames(spatial.data) %in% gene_vec, ]
    spatial.data = spatial.data[,Matrix::colSums(spatial.data) > n.umi.thresh]
    message(nrow(spatial.data), " genes used out of ", original_dim[1], " and ", ncol(spatial.data), " cells used out of ", original_dim[2])
  }
  
  coords = coords[colnames(spatial.data),]
  
  saveRDS(spatial.data, paste0(dir_spatial,"/",spatial.data.name,"_exp_qc_",descriptor,".RDS"))
  saveRDS(coords, paste0(dir_spatial,"/",spatial.data.name,"_coords_qc_",descriptor,".RDS"))
  saveRDS(rownames(spatial.data), paste0(dir_spatial,"/gene_selection_qc_",descriptor,".RDS"))
}

learn_gene_signatures =function(filepath,
                              region,
                              spatial.data.name,
                              lambda = 1,
                              thresh = 1e-8,
                              max.iters = 100,
                              nrep = 1,
                              rand.seed = 123,
                              print.obj = FALSE,
                              clusters.from.atlas = TRUE,
                              naive.clusters = FALSE,
                              verbose = FALSE){
  set.seed(rand.seed)
  
  descriptor = as.character(rand.seed)
  
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
    spatial.data.name = paste0(spatial.data.name, "_naive")
  }
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  
  norm.data = readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",descriptor,"_norm_data.RDS"))
  gene_vec = readRDS(paste0(dir_spatial,"/gene_selection_qc_",descriptor,".RDS"))
  gene_vec = intersect(gene_vec, rownames(norm.data[[1]]))
  
  sample.cells = readRDS(paste0(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/sampled_cells_", descriptor,".RDS")))
  
  if(clusters.from.atlas){
     clusters = readRDS(paste0(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/clusters_",descriptor,".RDS")))
  } else {
     clusters = readRDS(paste0(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/user_defined_clusters_",descriptor,".RDS")))
  }
  clusters = clusters[sample.cells]
  
  message("Learning gene signatures")
  E = lapply(norm.data, function(x){t(as.matrix(x)[gene_vec,])})
  
  if (!all(sapply(X = E, FUN = is.matrix))) {
    stop("All values in 'object' must be a matrix")
  }
  N <- length(x = E)
  ns <- sapply(X = E, FUN = nrow)
  #if (k >= min(ns)) {
  #  stop('Select k lower than the number of cells in smallest dataset: ', min(ns))
  #}
  tmp <- gc()
  
  clust_levels = levels(clusters)
  numeric_clust = as.numeric(clusters)
  names(numeric_clust) = names(clusters)
  H_cells = lapply(norm.data, colnames)
  clust_vec = rep(0, length(clust_levels))
  N = length(E)
  
  H_indic <- lapply(
    X = H_cells,
    FUN = function(n) {
      numeric_clust_sub = numeric_clust[n]
      mat = do.call(rbind, lapply(X = 1:length(n),
                                  FUN = function(cell){
                                    return(replace(clust_vec,numeric_clust[n[cell]],1))
                                  }))
      rownames(mat) = n
      colnames(mat) = as.character(clust_levels)
      return(mat)
    }
  )
  g <- ncol(x = E[[1]])
  k = length(clust_levels)
  W <- matrix(
    data = abs(x = runif(n = g * k, min = 0, max = 2)),
    nrow = k,
    ncol = g
  )
  
  V <- lapply(
    X = 1:N,
    FUN = function(i) {
      return(matrix(
        data = abs(x = runif(n = g * k, min = 0, max = 2)),
        nrow = k,
        ncol = g
      ))
    }
  )
  tmp <- gc()
  best_obj <- Inf
  V_m = V
  W_m = W
  H_m = H_indic
  
  #E = lapply(E, t)
  for (rep_num in 1:nrep) {
    set.seed(seed = rand.seed + rep_num - 1)
    start_time <- Sys.time()
    delta <- 1
    iters <- 0
    pb <- txtProgressBar(min = 0, max = max.iters, style = 3)
    sqrt_lambda <- sqrt(x = lambda)
    obj0 <- sum(sapply(
      X = 1:N,
      FUN = function(i) {
        return(norm(x = E[[i]] - H_indic[[i]] %*% (W + V[[i]]), type = "F") ^ 2)
      }
    )) +
      sum(sapply(
        X = 1:N,
        FUN = function(i) {
          return(lambda * norm(x = H_indic[[i]] %*% V[[i]], type = "F") ^ 2)
        }
      ))
    tmp <- gc()
    
    while (delta > thresh & iters < max.iters) {
      W <- rliger:::solveNNLS(
        C = do.call(rbind, H_indic),
        B = do.call(rbind, lapply(X = 1:N,
                                  FUN = function(i){
                                    return(E[[i]] - H_indic[[i]] %*% (W + V[[i]]))})
        )
      )
      
      tmp <- gc()
      V <- lapply(
        X = 1:N,
        FUN = function(i) {
          return(rliger:::solveNNLS(
            C = rbind(H_indic[[i]], sqrt_lambda * H_indic[[i]]),
            B = rbind(E[[i]] - H_indic[[i]] %*% W, matrix(data = 0, nrow = ns[[i]], ncol = g))
          ))
        }
      )
      tmp <- gc()
      obj <- sum(sapply(
        X = 1:N,
        FUN = function(i) {
          return(norm(x = E[[i]] - H_indic[[i]] %*% (W + V[[i]]), type = "F") ^ 2)
        }
      )) +
        sum(sapply(
          X = 1:N,
          FUN = function(i) {
            return(lambda * norm(x = H_indic[[i]] %*% V[[i]], type = "F") ^ 2)
          }
        ))
      tmp <- gc()
      delta <- abs(x = obj0 - obj) / (mean(obj0, obj))
      obj0 <- obj
      iters <- iters + 1
      setTxtProgressBar(pb = pb, value = iters)
    }
    setTxtProgressBar(pb = pb, value = max.iters)
    # if (iters == max.iters) {
    #   print("Warning: failed to converge within the allowed number of iterations.
    #         Re-running with a higher max.iters is recommended.")
    # }
    if (obj < best_obj) {
      W_m <- W
      V_m <- V
      best_obj <- obj
      best_seed <- rand.seed + rep_num - 1
    }
    
    
    
    if (verbose) {
      if (print.obj) {
        cat("Objective:", obj, "\n")
      }
      cat("Best results with seed ", best_seed, ".\n", sep = "")
    }
  }
  out <- list()
  out$H <- H_m
  for (i in 1:length(E)) {
    rownames(x = out$H[[i]]) <- rownames(x = E[[i]])
  }
  
  out$V <- V_m
  out$W <- W_m
  
  out$H_refined <- lapply(X = 1:length(E),
                          FUN = function(i){
                            mat = t(x = rliger:::solveNNLS(
                              C = rbind(t(x = W) + t(x = V[[i]]), sqrt_lambda * t(x = V[[i]])),
                              B = rbind(t(x = E[[i]]), matrix(data = 0, nrow = g, ncol = ns[[i]]))
                            )
                            )
                            rownames(mat) = rownames(E[[i]])
                            return(mat)
                          })
  
  names(x = out$V) <- names(x = out$H) <- names(x = out$H_refined) <- names(x = E)
  
  rownames(out$W) = clust_levels
  colnames(out$W) = gene_vec
  
  saveRDS(out, paste0(dir_spatial, "/gene_signature_output_",descriptor,".RDS"))
}

calculate_cell_sizes = function(
    filepath,
    region,
    naive.clusters = FALSE,
    rand.seed = 123,
    plot.hist = F
    ){
      
    object_paths = paste0(filepath,"/", region, "/Analysis",c(2,4,5),"_", region, "/onlineINMF_",region, "_object.RDS" )
    objects = lapply(object_paths, function(object_path){readRDS(object_path)})
    
    h5_files = Reduce(c, lapply(objects, function(object){sapply(1:length(object@norm.data), function(i){object@h5file.info[[i]]$file.path})}))
    rna_files = grep(paste0("_(sc10Xv3_|smartseq_|sn10Xv3_|sc10Xv2_)"), h5_files, value = TRUE)
    
    liger_cells = lapply(rna_files, function(i){
      rhdf5::h5read(i, "/matrix/barcodes")#change, extract from H5
    })
     
    liger_genes = lapply(rna_files, function(i){
      rhdf5::h5read(i, "/matrix/features")[[1]] #change, extract from H5
    })
    
    if(naive.clusters){
      clusters = readRDS(paste0(filepath,"/", region, "/", region, "_Deconvolution_Output/clusters_",rand.seed,"_object_clusters_naive.RDS"))
    } else {
      clusters = readRDS(paste0(filepath,"/", region, "/", region, "_Deconvolution_Output/clusters_",rand.seed,"_object_clusters.RDS"))
    }
    
    size_list = lapply(1:nlevels(clusters),function(x){vector(mode = "integer")})
    names(size_list) = levels(clusters)
    
    for(i in 1:length(rna_files)){
      
      n = rna_files[i]
      raw_data = Matrix::sparseMatrix(
        dims = c(length(liger_genes[[i]]), length(liger_cells[[i]])),
        i = as.numeric(rhdf5::h5read(n, "/matrix/indices")+1),
        p = as.numeric(rhdf5::h5read(n, "/matrix/indptr")),
        x = as.numeric(rhdf5::h5read(n, "/matrix/data"))
      )
      
      colnames(raw_data) = liger_cells[[i]]
      
      for(cell_type in names(size_list)){
       cells_subset = intersect(liger_cells[[i]], names(clusters)[clusters == cell_type])
       if(length(cells_subset)>1){
         size_list[[cell_type]]  = c(size_list[[cell_type]],Matrix::colSums(raw_data[,cells_subset]))
       } else if(length(cells_subset) == 1){
         size_list[[cell_type]] = c(size_list[[cell_type]],sum(raw_data[,cells_subset]))
       }
      }
    }
    
    cell_type_mean = sapply(size_list, mean)
    
    if(plot.hist){
      if(naive.clusters){
        pdf(file = paste0(filepath,"/", region, "/", region, "_Deconvolution_Output/histogram_size_",rand.seed,"_object_clusters_naive.PDF"), width = 6, height = 4)
      } else {
        pdf(file = paste0(filepath,"/", region, "/", region, "_Deconvolution_Output/histogram_size_",rand.seed,"_object_clusters.PDF"), width = 6, height = 4)
      }
      for(i in 1:length(size_list)){
        hist(size_list[[i]],main = paste0("Distribution of cell sizes for ",names(size_list)[i]), xlab = "Counts")
        abline(v = cell_type_mean[i], col = "red")
      }
      dev.off()
    }
    
    if(naive.clusters){
      saveRDS(cell_type_mean, paste0(filepath,"/", region, "/", region, "_Deconvolution_Output/cell_size_",rand.seed,"_object_clusters_naive.RDS"))
    } else {
      saveRDS(cell_type_mean, paste0(filepath,"/",region,"/", region, "_Deconvolution_Output/cell_size_",rand.seed,"_object_clusters.RDS"))
    }
}




deconvolve_spatial = function(filepath,
                              region,
                              spatial.data.name,
                              rand.seed = 123,
                              clusters.from.atlas = TRUE,
                              naive.clusters = FALSE,
                              cell.size = F,
                              W = NULL){
  set.seed(rand.seed)
  
  descriptor = as.character(rand.seed)
  
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
    spatial.data.name = paste0(spatial.data.name, "_naive")
  }
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  
  spatial.data = readRDS(paste0(dir_spatial,"/",spatial.data.name,"_exp_qc_",descriptor,".RDS"))
  
  if(is.null(W)){
    out = readRDS(paste0(dir_spatial, "/gene_signature_output_",descriptor,".RDS"))
    W = out[["W"]]
  } else {
    W = readRDS(W)
    descriptor = paste0(rand.seed, "_custom_w")
  }
  
  gene_vec = intersect(colnames(W), rownames(spatial.data))
  spatial.data = spatial.data[gene_vec,]
  W = W[,gene_vec]
  
  if(cell.size){
    if(naive.clusters){
      cell_size_mean = readRDS(paste0(filepath,"/",region, "/", region, "_Deconvolution_Output/cell_size_",rand.seed,"_object_clusters_naive.RDS"))
    } else {
      cell_size_mean = readRDS(paste0(filepath, "/", region, "/", region, "_Deconvolution_Output/cell_size_",rand.seed,"_object_clusters.RDS"))
    }
    W = t(sapply(1:nrow(W), function(i){return(cell_size_mean[i]*W[i,])}))
    rownames(W) = names(cell_size_mean)
    descriptor = paste0(descriptor, "_size_scaled")
  }
  
  message("Deconvolving spatial data")
  spatial.data = t(scale(t(as.matrix(spatial.data)), center = FALSE))
  spatial.data[spatial.data < 0 ] = 0
  spatial.data[is.nan(spatial.data)] = 0
  deconv_h = t(rliger:::solveNNLS(t(W),spatial.data))
  colnames(deconv_h) = rownames(W)
  deconv_frac = t(apply(deconv_h, MARGIN = 1, function(x){x/sum(x)}))
  rownames(deconv_frac) = rownames(deconv_h) = colnames(spatial.data)
  deconv_frac[is.nan(deconv_frac)] = 0
  saveRDS(list(raw = deconv_h, proportions = deconv_frac), paste0(dir_spatial,"/deconvolution_output_",descriptor,".RDS"))
  message("Deconvolution completed")
  
  dir.create(paste0(dir_spatial,"/",descriptor,"_output"))
  
}

#' Convert proportions generated during the deconvolution into cell-type
#' assignments, for use with single cell spatial modalities (i.e. slideseq)
#'
#' @param filepath Path to directory within which the atlas structure was generated
#' @param region A string corresponding to the name of an anatomical region
#' @param spatial.data.name A string, the name of the spatial dataset
#'
#' @return nothing
#'
#' @import
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }
assign_single_cells = function(
  filepath,
  region,
  spatial.data.name,
  rand.seed = 123,
  clusters.from.atlas = TRUE,
  naive.clusters = FALSE){
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  descriptor = as.character(rand.seed)
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
  }
 
  proportions = readRDS(paste0(dir_spatial,"/deconvolution_output_",descriptor,".RDS"))[[2]]
  cell_types = colnames(proportions)
  max = as.factor(apply(proportions, MARGIN = 1, function(x){cell_types[which.max(x)]}))
  names(max) = rownames(proportions)
  saveRDS(max, paste0(dir_spatial,"/deconvolution_max_factor_assignment_",descriptor,".RDS"))
}


#' Convert single cell spatial modalities into voxels by combining samples at
#' a preselected resolution
#'
#' @param filepath Path to directory within which the atlas structure was generated
#' @param region A string corresponding to the name of an anatomical region
#' @param spatial.data.name A string, the name of the spatial dataset
#' @param voxel.size A numeric, corresponding to the side-length of desired voxels
#' @param out.filepath A string corresponding to the directory to which to save generated data
#' @param verbose A logical, whether to print details on derived data
#' @return nothing
#'
#' @import
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }

voxelize_single_cells = function(
  filepath,
  region,
  spatial.data.name,
  voxel.size,
  out.filepath = NULL,
  verbose = TRUE
){
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  
  spatial.data = readRDS(paste0(dir_spatial,"/",spatial.data.name,"_exp.RDS"))
  coords = readRDS(paste0(dir_spatial,"/",spatial.data.name,"_coords.RDS"))


  minmax = apply(coords, MARGIN = 2, function(x){round(range(x),-1)})
  ranges = apply(minmax, MARGIN = 2, function(y){return(y[2]-y[1])})

  coords_adjusted = sapply(1:3, function(i){
    coords[,i] = voxel.size*round(coords[,i]/voxel.size)
  })
  colnames(coords_adjusted) = colnames(coords)

  coords_out = coords_adjusted[!duplicated(coords_adjusted),]
  voxel_exp = matrix(0L, nrow = nrow(spatial.data), ncol = nrow(coords_out))
  voxel_list = list()

  for(i in 1:nrow(coords_out)){
    voxels_use = coords_out[i,1] == coords_adjusted[,1]&
      coords_out[i,2] == coords_adjusted[,2]&
      coords_out[i,3] == coords_adjusted[,3]
    voxel_data = spatial.data[,voxels_use]
    voxel_list[[i]] = colnames(spatial.data)[voxels_use]
    if(!is.null(dim(voxel_data))){
      voxel_exp[,i] = Matrix::rowSums(voxel_data)
    } else {
      voxel_exp[,i] = voxel_data

    }
  }
  rownames(coords_out) = colnames(voxel_exp) = names(voxel_list) = paste0("voxel_",1:nrow(coords_out))
  rownames(voxel_exp) = rownames(spatial.data)


  if(is.null(out.filepath)){
    new_dir = paste0(dir_spatial, "_", voxel.size)
    dir.create(new_dir)
    message("Saving to new spatial data folder - " , new_dir)
    out.filepath = new_dir
  }
  saveRDS(voxel_exp, paste0(out.filepath,"/",spatial.data.name, "_", voxel.size, "_exp.RDS"))
  saveRDS(coords_out,paste0(out.filepath,"/",spatial.data.name, "_", voxel.size, "_coords.RDS"))
  saveRDS(voxel_list,paste0(out.filepath,"/",spatial.data.name, "_", voxel.size, "_voxels_to_samples.RDS"))

  if(verbose){
    message(paste0("Generated ", nrow(coords_out), " voxels at ", voxel.size, " cubed resolution." ))
    message(paste0("Mean samples per voxel: ", round(mean(sapply(voxel_list, length)))))
    message(paste0("Mean nUMI per voxel: ", round(mean(colSums(voxel_exp)))))
  }
}


#' Generate gifs of cell type distributions derived from deconvolution in space
#'
#' @param filepath Path to directory within which the atlas structure was generated
#' @param region A string corresponding to the name of an anatomical region
#' @param spatial.data.name A string, the name of the spatial dataset
#' @param mat.use A string, either "raw", "proportions", or "assignment",
#'    corresponding to the raw cell type loadings, the normalized loadings, or
#'    cell type assignments for single cell spatial modalities
#' @param cell.types.plot A character vector of cell types to plot
#' @param dims A vector of two integers, corresponding the dimensions of the
#'    gifs generated
#'
#' @return nothing
#'
#' @import rgl
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }

generate_loading_gifs = function(
  filepath,
  region,
  spatial.data.name,
  rand.seed = 123,
  clusters.from.atlas = TRUE,
  naive.clusters = FALSE,
  cell.size = FALSE,
  mat.use = "proportions",#raw, proportions, or assignment
  cell.types.plot = NULL,
  filter = NULL,
  dims = c(500, 500)
){
  set.seed(rand.seed)
  library(rgl)
  
  
  
  descriptor = as.character(rand.seed)
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
    spatial.data.name = paste0(spatial.data.name, "_naive")
  }
  
  coords = readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name,"/",spatial.data.name,"_coords_qc_",descriptor,".RDS"))

  if(cell.size){
     descriptor = paste0(descriptor, "_size_scaled")
  }
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  dir_output = paste0(dir_spatial,"/",descriptor,"_output")

  dir_gifs = paste0(dir_output,"/plots")
  if(!dir.exists(dir_gifs)){
    dir.create(paste0(dir_gifs))
    message("Created directory at ", dir_gifs)
  }
  
  if(mat.use != "assignment"){
    loadings = readRDS(paste0(dir_spatial,"/deconvolution_output_",descriptor,".RDS"))
    cell_types = colnames(loadings[[1]])
  } else {
    assignments = readRDS(paste0(dir_spatial,"/deconvolution_max_factor_assignment_",descriptor,".RDS"))
    cell_types = levels(assignments)
  }
  
  if(!is.null(filter) & mat.use != "assignment"){
    loadings[loadings < filter] = 0
    loadings[loadings >= filter] = 1
    descriptor = paste0(descriptor, "_filter_",filter)
  }
  
  
  grDevices::palette(viridis::viridis(option="A",n=50,direction = -1))
  if(is.null(cell.types.plot)){
    cell.types.plot = cell_types
  } else {
    cell.types.plot = intersect(cell.types.plot, cell_types)
  }
  


  for(cell_type in cell.types.plot){
    cell_type_consistent = sub("/", ".",sub(" ", "_",cell_type))
    if(mat.use != "assignment"){
      colors_view = as.numeric(cut(loadings[[mat.use]][,cell_type],breaks=50))
      try(rgl.close(), silent = TRUE)
      open3d(windowRect = c(0,0, dims[1], dims[2]));
      plot3d(coords[,1],coords[,2],coords[,3],col = colors_view,aspect=c(67,41,58),xlab="Anterior-Posterior",ylab="Inferior-Superior",zlab="Left-Right",size=5, add = TRUE)
      decorate3d(xlab = colnames(coords)[1], ylab = colnames(coords)[2],zlab = colnames(coords)[3], box = FALSE, axes = FALSE)
      axes3d(c("x--","y--","z--"))#axes3d(c("x--","y--","z--"))
      movie3d(spin3d(axis = c(0, 0, 1)), duration = 20, movie = paste0(dir_gifs,"/", region, "_",cell_type_consistent,"_spatial_summary_",descriptor))
    } else {
      colors_view = (assignments == cell_type)*20 + 1
      try(rgl.close(), silent = TRUE)
      open3d(windowRect = c(0,0, dims[1], dims[2]));
      plot3d(coords[,1],coords[,2],coords[,3],col = colors_view,aspect=c(67,41,58),xlab="Anterior-Posterior",ylab="Inferior-Superior",zlab="Left-Right",size=5, type = "p", add = TRUE)
      decorate3d(xlab = colnames(coords)[1], ylab = colnames(coords)[2],zlab = colnames(coords)[3], box = FALSE, axes = FALSE)
      axes3d(c("x--","y--","z--"))#axes3d(c("x--","y--","z--"))
      movie3d(spin3d(axis = c(0, 0, 1)), duration = 20, movie = paste0(dir_gifs, "/", region, "_",cell_type_consistent,"_spatial_summary_",descriptor))
    }

  }
}

generate_label_gifs = function(
  filepath,
  region,
  spatial.data.name,
  labels.plot,
  dims = c(500, 500)
){
  library(rgl)
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  
  dir_gifs = paste0(dir_spatial,"/plots")
  if(!dir.exists(dir_gifs)){
    dir.create(paste0(dir_gifs))
    message("Created directory at ", dir_gifs)
  }
  
  grDevices::palette(viridis::viridis(option="A",n=50,direction = -1))
  coords = readRDS(paste0(dir_spatial,"/",spatial.data.name,"_coords.RDS"))
  
  if(is.vector(labels.plot)){
    coords = coords[names(labels.plot),]
    labels_use = unique(labels.plot)
  } else {
    coords = coords[rownames(labels.plot),]
    labels_use = colnames(labels.plot)
  }
  for(label_unique in labels_use){
    if(is.vector(labels.plot)){
      colors_view = (labels.plot == label_unique)*20 + 1
    } else {
      colors_view = labels.plot[,label_unique]*20 + 1
    }
    try(rgl.close(), silent = TRUE)
    open3d(windowRect = c(0,0, dims[1], dims[2]));
    plot3d(coords[,1],coords[,2],coords[,3],col = colors_view,aspect=c(67,41,58),xlab="Anterior-Posterior",ylab="Inferior-Superior",zlab="Left-Right",size=5, type = "p", add = TRUE)
    decorate3d(xlab = colnames(coords)[1], ylab = colnames(coords)[2],zlab = colnames(coords)[3], box = FALSE, axes = FALSE)
    axes3d(c("x--","y--","z--"))#axes3d(c("x--","y--","z--"))
    label_unique = sub("/", ".",sub(" ", "_",label_unique))
    movie3d(spin3d(axis = c(0, 0, 1)), duration = 20, movie = paste0(dir_gifs, "/", region, "_",label_unique,"_spatial_summary"))
  }
}


#' Convert proportions generated during the deconvolution into cell-type
#' assignments, for use with single cell spatial modalities (i.e. slideseq)
#'
#' @param filepath Path to directory within which the atlas structure was generated
#' @param region A string corresponding to the name of an anatomical region
#' @param layer.list A named list of spatial sample names, corresponding to
#'    groupings to be summarized by the function
#' @param spatial.data.name A string, the name of the spatial dataset
#' @param plot A logical, corresponding with if to make bar plots corresponding
#'    to provided cell type and gene loading
#' @param type A string, either "mean" or "sum" corresponding with how to
#'    summarize across samples
#' @param mat.use A string, either "raw", "proportions", or "assignment",
#'    corresponding to the raw cell type loadings, the normalized loadings, or
#'    cell type assignments for single cell spatial modalities
#' @param use.cell.types A logical, if only the cell types provided to
#'    cell.types.use should be summarized, as opposed to all deconvolved.
#' @param cell.types.use A character vector of cell types to summarize
#' @param genes.use A character vector of genes to summarize
#'
#'
#' @return nothing
#'
#' @import
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }
summarize_by_layer = function(
  filepath,
  region,
  layer.list,
  spatial.data.name,
  plot = FALSE,
  type = "mean",
  mat.use = "proportions",#"assignment
  clusters.from.atlas = T,
  naive.clusters = F,
  use.cell.types = TRUE,
  cell.types.use = NULL,
  cell.size = FALSE,
  genes.use = NULL,
  rand.seed = 123){
  
  set.seed(rand.seed)
  
  descriptor = as.character(rand.seed)
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
    spatial.data.name = paste0(spatial.data.name, "_naive")
  }
  
  if(cell.size){
     descriptor = paste0(descriptor, "_size_scaled")
  }
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  dir_output = paste0(dir_spatial,"/",descriptor,"_output")
  
  if(mat.use == "assignment"){
    assignments = readRDS(paste0(dir_spatial,"/deconvolution_max_factor_assignment_",descriptor,".RDS"))
    sub_vec = rep(0,nlevels(assignments))
    loadings = Reduce(rbind, lapply(assignments, function(x){subbed_vec = sub_vec; subbed_vec[as.numeric(x)] = 1; return(subbed_vec)}))
    colnames(loadings) = levels(assignments)
  } else {
    loadings = readRDS(paste0(dir_spatial,"/deconvolution_output_",descriptor,".RDS"))[[mat.use]]
  }

  if(use.cell.types){
    if(!is.null(cell.types.use)){
      cell.types.use = intersect(cell.types.use, colnames(loadings))
    } else {
      cell.types.use = colnames(loadings)
    }
    cell.types.use = cell.types.use[cell.types.use != ""]
    cell.type.matrix = matrix(0L, nrow = length(layer.list), ncol = length(cell.types.use))
    rownames(cell.type.matrix) = names(layer.list)
    colnames(cell.type.matrix) = cell.types.use
    for(i in 1:length(layer.list)){
      sub_loadings = as.matrix(loadings[rownames(loadings) %in% as.character(layer.list[[i]]), ])
      if(ncol(sub_loadings) == 1){sub_loadings = t(sub_loadings)}
      for(j in 1:length(cell.types.use)){
        if(type == "mean"){
          cell.type.matrix[i,j] = mean(sub_loadings[ ,cell.types.use[j]])
        } else if(type == "sum"){
          cell.type.matrix[i,j] = sum(sub_loadings[ ,cell.types.use[j]])
        }
      }
    }
    saveRDS(cell.type.matrix, paste0(dir_output,"/",spatial.data.name,"_cell_type_layer_summary_",descriptor,".RDS"))
  }
  if(!is.null(genes.use)){
    spatial.data = readRDS(paste0(dir_spatial,"/",spatial.data.name,"_exp.RDS"))
    genes.use = intersect(genes.use, rownames(spatial.data))
    spatial.data[is.na(spatial.data)] = 0
    spatial.data = t(scale(t(spatial.data[genes.use,]), center = FALSE))
    spatial.data[spatial.data < 0 ] = 0
    gene.matrix = matrix(0L, nrow = length(layer.list), ncol = length(genes.use))
    rownames(gene.matrix) = names(layer.list)
    colnames(gene.matrix) = genes.use
    for(i in 1:length(layer.list)){
      sub_loadings = as.matrix(spatial.data[ ,colnames(spatial.data) %in% as.character(layer.list[[i]])])
      if(ncol(sub_loadings) == 1){sub_loadings = t(sub_loadings)}
      for(j in 1:length(genes.use)){
        if(type == "mean"){
          gene.matrix[i,j] = mean(sub_loadings[genes.use[j],])
        } else if(type == "sum"){
          gene.matrix[i,j] = sum(sub_loadings[genes.use[j],])
        }
      }
    }
    saveRDS(gene.matrix, paste0(dir_output,"/",spatial.data.name,"gene_layer_summary_",descriptor,".RDS"))
  }
  if(plot){
    dir_plots = paste0(dir_output,"/plots")
    if(!dir.exists(dir_plots)){
      dir.create(paste0(dir_plots))
      message("Created directory at ", dir_plots)
    }
    ggplot2::theme_set(cowplot::theme_cowplot())
    labels.cell.type = expand.grid(rownames(cell.type.matrix), colnames(cell.type.matrix))
    cell.type.df = data.frame(Layers = as.character(labels.cell.type[,1]),
                              Cell_Types = as.character(labels.cell.type[,2]),
                              Values = as.vector(cell.type.matrix))
 
    if(ncol(cell.type.matrix) > 1){
      overall.cell.type.plot = ggplot2::ggplot(cell.type.df, ggplot2::aes(fill = Cell_Types, y = Values, x = Layers)) +
        ggplot2::theme(text = ggplot2::element_text(size = 10),
                       axis.text = ggplot2::element_text(size = 5),
                       legend.title = ggplot2::element_blank(),
                       legend.text = ggplot2::element_text(size = 3),
                       legend.key.height = ggplot2::unit(3, 'mm'),
                       legend.key.width = ggplot2::unit(1, 'mm')) +
        ggplot2::geom_bar(position = "dodge", stat = "identity") +
        ggplot2::xlab("Layer") +
        ggplot2::ylab("Value") +
        ggplot2::ggtitle(paste0("Distribution of cell types by layer"))

      print(overall.cell.type.plot)

      ggplot2::ggsave(paste0(dir_plots, "/cell_type_layer_distribution_",descriptor,".PNG"),
                      overall.cell.type.plot,
                      width = 1000,
                      height = 800,
                      units = "px")
    }

    for(i in colnames(cell.type.matrix)){
      cell_type_consistent = sub("/", ".",sub(" ", "_",i))
      cell.type.df.sub = cell.type.df[cell.type.df$Cell_Types == i,]
      by.cell.type.plot = ggplot2::ggplot(cell.type.df.sub, ggplot2::aes(fill = Layers, x = Layers, y = Values)) +
        ggplot2::theme(text = ggplot2::element_text(size = 3),
                       axis.text = ggplot2::element_text(size = 2),
                       legend.position="none") +
        ggplot2::geom_bar(position = "dodge", stat = "identity") +
        ggplot2::xlab("Layer") +
        ggplot2::ylab("Value") +
        ggplot2::ggtitle(paste0("Distribution of ",cell_type_consistent, " cells by layer"))
        ggplot2::ggsave(paste0(dir_plots, "/",cell_type_consistent,"_layer_distribution_",descriptor,".PNG"),
                      by.cell.type.plot,
                      width = 500,
                      height = 400,
                      units = "px")
        message("Plotting ", cell_type_consistent)
    }
    if(!is.null(genes.use)){
      labels.genes = expand.grid(rownames(gene.matrix), colnames(gene.matrix))
      gene.df = data.frame(Layers = as.character(labels.genes[,1]),
                           Genes = as.character(labels.genes[,2]),
                           Values = as.vector(gene.matrix))
      if(length(genes.use) > 1){
        overall.gene.plot = ggplot2::ggplot(gene.df, ggplot2::aes(fill = Genes, y = Values, x = Layers)) +
          ggplot2::theme(text = ggplot2::element_text(size = 8),
                         axis.text = ggplot2::element_text(size = 5),
                         legend.title = ggplot2::element_blank(),
                         legend.text = ggplot2::element_text(size = 3),
                         legend.key.height = ggplot2::unit(3, 'mm'),
                         legend.key.width = ggplot2::unit(1, 'mm')) +
          ggplot2::geom_bar(position = "dodge", stat = "identity") +
          ggplot2::xlab("Layer") +
          ggplot2::ylab("Value") +
          ggplot2::ggtitle(paste0("Distribution of gene expression by layer"))

        print(overall.gene.plot)

        ggplot2::ggsave(paste0(dir_plots, "/gene_layer_distribution_",descriptor,".PNG"),
                        overall.gene.plot,
                        width = 1000,
                        height = 800,
                        units = "px")
      }

      for(i in colnames(gene.matrix)){
        gene_consistent = sub("/", ".",sub(" ", "_",i))
        gene.df.sub = gene.df[gene.df$Genes == i,]
        by.gene.plot = ggplot2::ggplot(gene.df.sub, ggplot2::aes(fill = Layers, x = Layers, y = Values)) +
          ggplot2::theme(text = ggplot2::element_text(size = 3),
                         axis.text = ggplot2::element_text(size = 2),
                         legend.position="none") +
          ggplot2::geom_bar(position = "dodge", stat = "identity") +
          ggplot2::xlab("Layer") +
          ggplot2::ylab("Value") +
          ggplot2::ggtitle(paste0("Distribution of ",i, " expression by layer"))
        ggplot2::ggsave(paste0(dir_plots, "/",gene_consistent,"_gene_layer_distribution_",descriptor,".PNG"),
                        by.gene.plot,
                        width = 500,
                        height = 400,
                        units = "px")

      }
    }
  }
}

analyze_gene_signatures = function(filepath,
                                   region,
                                   spatial.data.name,
                                   clusters.from.atlas = TRUE,
                                   naive.clusters = FALSE,
                                   plot = FALSE,
                                   mat.use = "proportions",
                                   cell.types.use = NULL,
                                   cell.size = FALSE,
                                   rand.seed = 123){
  
  library(ggplot2)
  
  set.seed(rand.seed)
  
  descriptor = as.character(rand.seed)
  
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
    spatial.data.name = paste0(spatial.data.name, "_naive")
  }
  
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
    
  gene_sigs = readRDS(paste0(dir_spatial, "/gene_signature_output_",descriptor,".RDS"))
  genes = readRDS(paste0(dir_spatial,"/gene_selection_qc_",descriptor,".RDS"))
  
  if(cell.size){
     descriptor = paste0(descriptor, "_size_scaled")
  }
  
    dir_output = paste0(dir_spatial,"/",descriptor,"_output")

  
  
  signatures = gene_sigs$W
  signatures = signatures[rowSums(signatures) != 0 & rownames(signatures) != "",]
  
  if(!is.null(cell.types.use)){
    signatures = signatures[intersect(rownames(signatures),cell.types.use),]
  }
  
  cos_sim = lsa::cosine(t(signatures))
  
  
  cos_dist= as.dist(1- cos_sim)
  hierarchical_clust <- hclust(cos_dist, method = "ward.D2")
  
  saveRDS(list(cos_sim_signature = cos_sim, dendro_sig = hierarchical_clust),
          paste0(dir_output,"/gene_signature_analysis_summary_",descriptor,".RDS"))
  if(plot){
    
    dir_plots = paste0(dir_output,"/plots")
    if(!dir.exists(dir_plots)){
      dir.create(paste0(dir_plots))
      message("Created directory at ", dir_plots)
    }
    
    heatmap_df = data.frame(expand.grid(Cell_Type_1 = rownames(cos_sim),
                                        Cell_Type_2 = colnames(cos_sim)),
                            cos_sim = as.vector(cos_sim))
    heatmap_plot = ggplot2::ggplot(heatmap_df, ggplot2::aes(x = Cell_Type_1, y = Cell_Type_2, fill = cos_sim)) +
      ggplot2::labs(y = "Cell Types", fill = "", title = "Cosine Similarity for Cell Type Signatures") +
      ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                     axis.title.x = ggplot2::element_blank(),
                     text = ggplot2::element_text(size = 8),
                     axis.text = ggplot2::element_text(size = 5),
                     legend.title = ggplot2::element_blank(),
                     legend.text = ggplot2::element_text(size = 3),
                     legend.key.height = ggplot2::unit(3, 'mm'),
                     legend.key.width = ggplot2::unit(1, 'mm')) +
      ggplot2::geom_tile() +
      viridis::scale_fill_viridis()
    
    dendro_plot = ggdendro::ggdendrogram(hierarchical_clust, rotate = FALSE, size = 2) +
      ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     text = ggplot2::element_text(size = 4),
                     axis.text = ggplot2::element_text(size = 4)) +
      ggplot2::ggtitle("Hierarchical Clustering", subtitle = "By Cell Type Signature")
    
    print(heatmap_plot)
    print(dendro_plot)
   
    ggplot2::ggsave(paste0(dir_plots, "/cell_type_signature_heatmap_",descriptor,".PNG"),
                    heatmap_plot,
                    width = 1500,
                    height = 1200,
                    units = "px")
    ggplot2::ggsave(paste0(dir_plots, "/cell_type_signature_dendrogram_",descriptor,".PNG"),
                    dendro_plot,
                    width = 1000,
                    height = 800,
                    units = "px")
  }
}

cell_type_loading_histogram = function(
  filepath,
  region,
  spatial.data.name,
  rand.seed = 123,
  clusters.from.atlas = TRUE,
  naive.clusters = FALSE,
  cell.size = FALSE,
  mat.use = "proportions",
  cell.types.plot = NULL,
  print.plots = FALSE,
  bin.num = 30){

  library(ggplot2)
  
  set.seed(rand.seed)

  descriptor = as.character(rand.seed)
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
    spatial.data.name = paste0(spatial.data.name, "_naive")
  }
  if(cell.size){
     descriptor = paste0(descriptor, "_size_scaled")
  }

  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  dir_output = paste0(dir_spatial,"/",descriptor,"_output")

  dir_gifs = paste0(dir_output,"/plots")
  if(!dir.exists(dir_gifs)){
    dir.create(paste0(dir_gifs))
    message("Created directory at ", dir_gifs)
  }

  loadings = readRDS(paste0(dir_spatial,"/deconvolution_output_",descriptor,".RDS"))[[mat.use]]
  cell_types = colnames(loadings)

  if(!is.null(cell.types.plot)){
    cell.types.plot = intersect(cell.types.plot, cell_types)
  } 
  cell.types.plot = gsub(" ", ".", gsub("-","_",cell.types.plot))
  colnames(loadings) = gsub(" ", ".", gsub("-","_",colnames(loadings)))

  loadings = as.data.frame(loadings[,cell.types.plot])

  hist_plot = list()
  for(cell_type in colnames(loadings)){
    range = loadings[,cell_type][2]
    hist_plot[[cell_type]] = ggplot(loadings, aes_string(x=cell_type)) + 
      geom_histogram(bins = bin.num, fill = rainbow(bin.num)) +
      labs(y = "Count", x = "",title = paste0("Histogram of ",cell_type, " loading by voxel"))
    if(print.plots){
      print(hist_plot[[cell_type]])
    }
  }


  pdf(file = paste0(dir_gifs, "/",descriptor,"_",mat.use,"_hist.PDF"), width = 7,height = 4)
  for(cell_type in names(hist_plot)){
    print(hist_plot[[cell_type]])
  }
  dev.off()
}

analyze_spatial_correlation = function(filepath,
                                       region,
                                       spatial.data.name,
                                       clusters.from.atlas = TRUE,
                                       naive.clusters = FALSE,
                                       cell.size = FALSE,
                                       plot = FALSE,
                                       mat.use = "proportions",
                                       cell.types.use = NULL,
                                       rand.seed = 123){
  
  library(ggplot2)
  
  set.seed(rand.seed)
  
  descriptor = as.character(rand.seed)
  
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
    spatial.data.name = paste0(spatial.data.name, "_naive")
  }
  
  if(cell.size){
     descriptor = paste0(descriptor, "_size_scaled")
  }
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  dir_output = paste0(dir_spatial,"/",descriptor,"_output")
  
  deconv_out = readRDS(paste0(dir_spatial,"/deconvolution_output_",descriptor,".RDS"))
  loadings = deconv_out[[mat.use]]
  loadings = loadings[, colSums(loadings) != 0 & colnames(loadings)!=""]
  
  if(!is.null(cell.types.use)){
    loadings = loadings[,intersect(colnames(loadings),cell.types.use)]
  }
  
  corr_sim_dist = cor(loadings)
  
  
  
  cor_dist= as.dist(1- corr_sim_dist)
  hierarchical_clust_dist <- hclust(cor_dist, method = "ward.D2")
  
  saveRDS(list(cos_sim_signature = corr_sim_dist, dendro_sig = hierarchical_clust_dist),
          paste0(dir_output,"/spatial_correlation_analysis_summary_",descriptor,".RDS"))
  
  if(plot){
    
    dir_plots = paste0(dir_output,"/plots")
    if(!dir.exists(dir_plots)){
      dir.create(paste0(dir_plots))
      message("Created directory at ", dir_plots)
    }    
    
    heatmap_dist_df = data.frame(expand.grid(Cell_Type_1 = rownames(corr_sim_dist),
                                             Cell_Type_2 = colnames(corr_sim_dist)),
                                 cos_sim = as.vector(corr_sim_dist))
    heatmap_dist_plot = ggplot2::ggplot(heatmap_dist_df, ggplot2::aes(x = Cell_Type_1, y = Cell_Type_2, fill = cos_sim)) +
      ggplot2::labs(y = "Cell Types", fill = "", title = "Correlation for Cell Type Distribution") +
      ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                     axis.title.x = ggplot2::element_blank(),
                     text = ggplot2::element_text(size = 8),
                     axis.text = ggplot2::element_text(size = 5),
                     legend.title = ggplot2::element_blank(),
                     legend.text = ggplot2::element_text(size = 3),
                     legend.key.height = ggplot2::unit(3, 'mm'),
                     legend.key.width = ggplot2::unit(1, 'mm')) +
      ggplot2::geom_tile() +
      viridis::scale_fill_viridis()
    
    dendro_dist_plot = ggdendro::ggdendrogram(hierarchical_clust_dist, rotate = FALSE, size = 2) +
      ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     text = ggplot2::element_text(size = 4),
                     axis.text = ggplot2::element_text(size = 4)) +
      ggplot2::ggtitle("Hierarchical Clustering", subtitle = "By Cell Type Distribution")

    print(heatmap_dist_plot)
    print(dendro_dist_plot)

    
    ggplot2::ggsave(paste0(dir_plots, "/spatial_correlation_heatmap_",descriptor,".PNG"),
                    heatmap_dist_plot,
                    width = 1500,
                    height = 1200,
                    units = "px")
    ggplot2::ggsave(paste0(dir_plots, "/spatial_correlation_dendrogram_",descriptor,".PNG"),
                    dendro_dist_plot,
                    width = 1000,
                    height = 800,
                    units = "px")
  }
}


calculate_wasserstein = function(
    filepath,
    region,
    spatial.data.name,
    clusters.from.atlas = TRUE,
    naive.clusters = FALSE,
    cell.size = FALSE,
    mat.use = "proportions",
    use.cell.types = TRUE,
    cell.types.use = NULL,
    genes.use = NULL,
    p = 2,
    min.samples = 1,
    plot = FALSE,
    rand.seed = 123){
  
  library(ggplot2)
  
  set.seed(rand.seed)
  
  descriptor = as.character(rand.seed)
  
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
    spatial.data.name = paste0(spatial.data.name, "_naive")
  }
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)

  coords = readRDS(paste0(dir_spatial,"/",spatial.data.name,"_coords_qc_",descriptor,".RDS"))
  
  if(cell.size){
     descriptor = paste0(descriptor, "_size_scaled")
  }
  
  dir_output = paste0(dir_spatial,"/",descriptor,"_output")
  
  deconv_out = readRDS(paste0(dir_spatial,"/deconvolution_output_",descriptor,".RDS"))
  loadings = deconv_out[[mat.use]]
  loadings = loadings[, colSums(loadings) != 0 & colnames(loadings)!=""]
  
  
  loadings = loadings[, colSums(loadings) != 0]
  if(use.cell.types){
    if(!is.null(cell.types.use)){
      cell.types.use = intersect(cell.types.use, colnames(loadings))
    } else {
      cell.types.use = colnames(loadings)
    }
    cell.types.use = cell.types.use[cell.types.use != "" & apply(loadings, MARGIN = 2, function(x){sum(x != 0)}) > min.samples]
    loadings = loadings[rownames(loadings) %in% rownames(coords), cell.types.use]
    
    
    if(is.vector(loadings)){
      new_loadings = matrix(loadings)
      rownames(new_loadings) = names(loadings)
      colnames(new_loadings) = cell.types.use
      loadings = new_loadings
      rm(new_loadings)
    }
    
    colnames(loadings) = sub("/",".",sub(" ", "_", cell.types.use))
    
    coords = coords[rownames(loadings),]
  }
  
  if(!is.null(genes.use)){
    exp = t(readRDS(paste0(dir_spatial,"/",spatial.data.name,"_exp_qc_",descriptor,".RDS")))
    exp = exp[, colnames(exp) %in% genes.use]
    exp[exp < 0] = 0
    exp = exp[rownames(loadings),]
    
    loadings = cbind(loadings, exp)
  }
  loadings = scale(loadings, center = F)
  library(transport)
  
  vars_1 = vars_2 = colnames(loadings)
  distance_mat = matrix(0L, nrow = ncol(loadings), ncol = ncol(loadings))
  colnames(distance_mat) = rownames(distance_mat) = vars_1
  loadings = apply(loadings, MARGIN = 2, function(x){x/sum(x)})
  for(var_1 in vars_1){
    mass_1 =  loadings[,var_1]
    distribution_1 = wpp(coords,mass_1)
    for(var_2 in vars_2){
      mass_2 =  loadings[,var_2]
      distribution_2 = wpp(coords,mass_2)
      distance_mat[var_1, var_2] = 
        distance_mat[var_2, var_1] = 
        wasserstein(distribution_1, distribution_2, p)
    }
    vars_2 = vars_2[2:length(vars_2)]
  }
  
  saveRDS(distance_mat,
          paste0(dir_output,"/wasserstein_distance_mat_",descriptor,".RDS"))
  
  if(plot){
    
    dir_plots = paste0(dir_output,"/plots")
    if(!dir.exists(dir_plots)){
      dir.create(paste0(dir_plots))
      message("Created directory at ", dir_plots)
    }    
    
    heatmap_wasserstein_df = data.frame(expand.grid(Cell_Type_1 = rownames(distance_mat),
                                             Cell_Type_2 = colnames(distance_mat)),
                                 wasserstein_dist = as.vector(distance_mat))
    heatmap_wasserstein_plot = ggplot2::ggplot(heatmap_wasserstein_df, ggplot2::aes(x = Cell_Type_1, y = Cell_Type_2, fill = wasserstein_dist)) +
      ggplot2::labs(y = "Cell Types", fill = "", title = "Wasserstein Distance by Cell Type and Gene") +
      ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                     axis.title.x = ggplot2::element_blank(),
                     text = ggplot2::element_text(size = 8),
                     axis.text = ggplot2::element_text(size = 5),
                     legend.title = ggplot2::element_blank(),
                     legend.text = ggplot2::element_text(size = 3),
                     legend.key.height = ggplot2::unit(3, 'mm'),
                     legend.key.width = ggplot2::unit(1, 'mm')) +
      ggplot2::geom_tile() +
      viridis::scale_fill_viridis()
    
    ggplot2::ggsave(paste0(dir_plots, "/wasserstein_heatmap_",descriptor,".PNG"),
                    heatmap_wasserstein_plot,
                    width = 1500,
                    height = 1200,
                    units = "px")
  }
  
}
                                                                                                   
refine_cluster_similarity = function(
    filepath,
    region,
    spatial.data.name,
    clusters.from.atlas = TRUE,
    naive.clusters = FALSE,
    cell.size = FALSE,
    rand.seed){
  library(ggplot2)
  
  set.seed(rand.seed)
  
  descriptor = as.character(rand.seed)
  
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
    spatial.data.name = paste0(spatial.data.name, "_naive")
  }
  
  if(cell.size){
     descriptor = paste0(descriptor, "_size_scaled")
  }
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  dir_output = paste0(dir_spatial,"/",descriptor,"_output")
  
  spatial = readRDS(paste0(dir_output,"/wasserstein_distance_mat_",descriptor,".RDS"))
  
  expression = readRDS(paste0(dir_output,"/gene_signature_analysis_summary_",descriptor,".RDS"))[["cos_sim_signature"]]
  shared_clusters = intersect(colnames(spatial), colnames(expression))
  spatial = spatial[shared_clusters, shared_clusters]
  expression = expression[shared_clusters, shared_clusters]
  placeholder_vec = rep("",length(shared_clusters))
  df_out = data.frame("wasserstein_10" = placeholder_vec,
                      "wasserstein_25" = placeholder_vec,
                      "wasserstein_50" = placeholder_vec,
                      "corr_.9" = placeholder_vec,
                      "corr_.75" = placeholder_vec,
                      "corr_.5" = placeholder_vec)
  rownames(df_out) = shared_clusters
  for(cluster in shared_clusters){
    #come back to this part to figure out how to convert these values into concatenated strings for each output of the clusters denoted.
    clusts_use = setdiff(shared_clusters, cluster)
    wasserstein_vec = spatial[cluster, clusts_use]
    wasserstein_thresh = c(0,as.numeric(quantile(wasserstein_vec, probs=c(.1,.25,.5))))
    wasserstein_strings = sapply(1:3, function(x){paste(clusts_use[wasserstein_vec > wasserstein_thresh[x] & wasserstein_vec <= wasserstein_thresh[x+1]], collapse = ",")})
    corr_thresh = c(.5, .75, .9, 1)
    corr_vec = expression[cluster,clusts_use]
    corr_strings = sapply(3:1, function(x){paste(clusts_use[corr_vec >= corr_thresh[x] & corr_vec < corr_thresh[x+1]], collapse = ",")})
    df_out[cluster,] = c(wasserstein_strings, corr_strings)
  }
  write.csv(df_out,paste0(dir_output,"/",spatial.data.name,"_cluster_similarity_",descriptor,".csv"),row.names = TRUE)    
}


threshold_similar_clusts = function(
    filepath,
    region,
    spatial.data.name,
    clust.compare,
    quantile.wasserstein = .1,
    correlation.val = .9,
    rand.seed = 123){
  library(ggplot2)
  
  set.seed(rand.seed)
  
  descriptor = as.character(rand.seed)
  
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
    spatial.data.name = paste0(spatial.data.name, "_naive")
  }
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  dir_output = paste0(dir_spatial,"/",descriptor,"_output")
  
  spatial = readRDS(paste0(dir_output,"/wasserstein_distance_mat_",descriptor,".RDS"))
  
  expression = readRDS(paste0(dir_output,"/gene_signature_analysis_summary_",descriptor,".RDS"))[["cos_sim_signature"]]
  shared_clusters = intersect(colnames(spatial), colnames(expression))
  clusts_use = setdiff(shared_clusters, clust.compare)
  spatial = spatial[clust.compare, clusts_use]
  expression = expression[clust.compare, clusts_use]
  
  wasserstein_thresh = clusts_use[wasserstein_vec <= quantile(wasserstein_vec, probs=quantile.wasserstein)[1]]
  corr_thresh = clusts_use[expression >= correlation.val]
  thresh_clusts = intersect(wasserstein_thresh, corr_thresh)
  print("Clusters meeting provided thresholds:")
  print(thresh_clusts)
  
  thresh_file = paste0(dir_output, "/thresh_clusters_",descriptor,".RDS")
  if(file.exists(thresh_file)){
    thresh_list = readRDS(thresh_file)
    thresh_file[[paste0("wass_",quantile.wasserstein,"_corr_",correlation.val)]] = thresh_clusts
    saveRDS(thresh_list, thresh_file)
  } else {
    thresh_list = list()
    thresh_list[[paste0("wass_",quantile.wasserstein,"_corr_",correlation.val)]] = thresh_clusts
    saveRDS(thresh_list, thresh_file)
  }
}


describe_voxelized_loading_by_label = function(filepath,
                              region,
                              spatial.data.name,
                              rand.seed = 123,
                              clusters.from.atlas = TRUE,
                              naive.clusters = FALSE,
                              mat.use = "raw",
                              labels.use){
  set.seed(rand.seed)

  descriptor = as.character(rand.seed)
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
    spatial.data.name = paste0(spatial.data.name, "_naive")
  }
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  
  
  voxels_to_samples = readRDS(paste0(dir_spatial, "/",spatial.data.name,"_voxels_to_samples.RDS"))
  raw_loadings = readRDS(paste0(dir_spatial, "/deconvolution_output_",descriptor,".RDS"))[[mat.use]]
  
  unique_labels = unique(labels.use)
  
  voxel_in_subregion = sapply(unique_labels, function(unique_label){
    sapply(names(voxels_to_samples), function(voxel_to_sample){
      any(labels.use[voxels_to_samples[[voxel_to_sample]]] %in% unique_label)
    })
  })
  
  colnames(voxel_in_subregion) = unique_labels
  rownames(voxel_in_subregion) = names(voxels_to_samples)
  
  voxel_in_subregion = voxel_in_subregion[rownames(voxel_in_subregion) %in% rownames(raw_loadings),]
  
  
  proportion_loading_in_subregion = sapply(unique_labels, function(subregion){
    sapply(colnames(raw_loadings), function(cell_type){
      sum(raw_loadings[voxel_in_subregion[,subregion],cell_type])/sum(raw_loadings[,cell_type])
    })
  })
  
  saveRDS(proportion_loading_in_subregion, paste0(dir_spatial, "/",descriptor,"_output/cell_type_loading_by_label_",mat.use,".RDS"))
}
register_voxel_to_label = function(filepath,
                              region,
                              spatial.data.name,
                              labels.use,
                              label.name){

    dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)


    voxels_to_samples = readRDS(paste0(dir_spatial, "/",spatial.data.name,"_voxels_to_samples.RDS"))
  
    unique_labels = unique(labels.use)

    voxel_to_subregion = sapply(unique_labels, function(unique_label){
      label_counts = sapply(names(voxels_to_samples), function(voxel_to_sample){
        sum(labels.use[voxels_to_samples[[voxel_to_sample]]] %in% unique_label)
      })
      return(unique_label[which.max(label_counts)])
    })
  
    names(voxel_to_subregion) = names(voxels_to_samples)
  
    saveRDS(proportion_loading_in_subregion, paste0(dir_spatial, "/voxel_assignment_by_label_",label.name,".RDS"))
  }

mirror_spatial_coords = function(filepath,
                              region,
                              spatial.data.name,
                              axes.flip = c(FALSE,FALSE,FALSE),
                              overwrite = T){
  deconv_dir = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/")
  coords = readRDS(paste0(deconv_dir,spatial.data.name,"/",spatial.data.name,"_coords.RDS"))
  
  descriptor = "mirror"
  axis_des = c("_x","_y", "_z")

  for(i in 1:3){
    if(axes.flip[i]){
      coords_range = range(coords[,i])
      coords[,i] = -1*(coords[,i]-coords_range[1])+coords_range[2]
      descriptor = paste0(descriptor, axis_des[i])
    }
  }
  
  if(overwrite){
    saveRDS(coords, paste0(deconv_dir,spatial.data.name,"/",spatial.data.name,"_coords.RDS"))
  } else {
    save_spatial_data(filepath,
                      region,
                      paste0(deconv_dir,spatial.data.name,"/",spatial.data.name,"_exp.RDS"),
                      coords,
                      paste0(spatial.data.name,"_",descriptor))
  }
  
}

find_common_coords = function(filepath,
                                 region,
                                 spatial.data.name.1,
                                 spatial.data.name.2,
                                 overwrite = T){
  deconv_dir = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/")
  coords_1 = readRDS(paste0(deconv_dir,spatial.data.name.1,"/",spatial.data.name.1,"_coords.RDS"))
  coords_2 = readRDS(paste0(deconv_dir,spatial.data.name.2,"/",spatial.data.name.2,"_coords.RDS"))
  coords_1_range = apply(coords_1, MARGIN = 2, range)
  coords_2_range = apply(coords_2, MARGIN = 2, range)
  coords_2_realigned = sapply(1:3, function(i){
    return((coords_1_range[2,i]-coords_1_range[1,i])*(coords_2[,i]-coords_2_range[1,i])/(coords_2_range[2,i]-coords_2_range[1,i])+coords_1_range[1,i])
  })
  
  rownames(coords_2_realigned) = rownames(coords_2)
  colnames(coords_2_realigned) = colnames(coords_2)
  
  if(overwrite){
    saveRDS(coords_2_realigned, paste0(deconv_dir,spatial.data.name.2,"/",spatial.data.name.2,"_coords.RDS"))
  } else {
    save_spatial_data(filepath,
                      region,
                      paste0(deconv_dir,spatial.data.name.2,"/",spatial.data.name.2,"_exp.RDS"),
                      coords_2_realigned,
                      paste0(spatial.data.name.2,"_aligned"))
  }
}


transform_coords_to_ccf = function(
    filepath,
    region,
    spatial.data.name,
    ish = T,
    overwrite = T){
  scale_factor = if(ish){200}else{25}
  deconv_dir = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/")
  coords = readRDS(paste0(deconv_dir,spatial.data.name,"/",spatial.data.name,"_coords.RDS"))
  coords = coords*scale_factor
  if(ish){
    coords[,3] = -coords[,3]+min(coords[,3])+max(coords[,3])+2400
    coords[,2] = -coords[,2]+min(coords[,2])+max(coords[,2])+5500
  } else {
    coords[,2] = -coords[,2]+min(coords[,2])+max(coords[,2])+3000
  }
  if(overwrite){
    files_dir = list.files(paste0(deconv_dir,spatial.data.name), full.names = T) 
    files_coords = files_dir[grep("coords", files_dir)]
    for(file_coord in files_coords){
      coords_file = readRDS(file_coord)
      saveRDS(coords[rownames(coords) %in% rownames(coords_file),], file_coord)
    }
    message("Rerun downstream plotting functions (generate_loading_gifs, calculate_wasserstein, etc.) to update with transformed coordinates")
  } else {
    saveRDS(coords, paste0("~/",region, "_",spatial.data.name, "_in_ccf.RDS"))
  }
}


summarize_subregions = function(regions, ontology.file = "Downloads/allen_structure_ontology.csv", return = F){
  allen_structure = read.csv(ontology.file,header=T)
  allen_structure_list = lapply(allen_structure$structure_id_path, function(x){strsplit(x, "/")}[[1]])
  
  subregion_vec = vector(mode = "character")
  for(region in regions){
    message(region)
    sub_allen = allen_structure[allen_structure$acronym == region, ]
    if(nrow(sub_allen)>1){
      id = sub_allen$id[which.max(sub_allen$depth)]
    } else {
      id = sub_allen$id[1]
    }
    if(length(id) != 0){
      subregions = allen_structure[sapply(allen_structure_list, function(x){id %in% x}) & grepl(region, allen_structure$acronym),]
      print(subregions[,c("acronym","name")])
      subregion_vec = c(subregion_vec, subregions[,"acronym"])
    } else {
      message("No subregions!")
    }
  }
  if(return){
    return(unique(subregion_vec))
  }
}

summarize_clusters = function(filepath,
                             region,
                             rand.seed,
                             naive = FALSE,
                             return = FALSE){
  deconv_dir = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/")
  if(naive){
    clusts = readRDS(paste0(deconv_dir, "clusters_", rand.seed,"_object_clusters_naive.RDS"))
  } else {
    clusts = readRDS(paste0(deconv_dir, "clusters_", rand.seed,"_object_clusters.RDS"))
  }
  message(paste0(region, " cluster frequency"))
  print(table(clusts))
  if(return){
    return(levels(clusts))
  }
}


overlay_subregion_gifs = function(
    filepath,
    region,
    spatial.data.name,
    rand.seed = 123,
    clusters.from.atlas = TRUE,
    naive.clusters = FALSE,
    cell.size = FALSE,
    mat.use = "proportions",#raw, proportions, or assignment
    cell.types.plot = NULL,
    subregions.plot = NULL,
    filter = NULL,
    dims = c(500, 500)
){
  set.seed(rand.seed)
  library(rgl)
  library(cocoframer)
  
  descriptor = as.character(rand.seed)
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
    spatial.data.name = paste0(spatial.data.name, "_naive")
  }
    
  if(cell.size){
     descriptor = paste0(descriptor, "_size_scaled")
  }
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  
  coords = readRDS(paste0(dir_spatial,"/",spatial.data.name,"_coords_qc_",descriptor,".RDS"))
  
  dir_output = paste0(dir_spatial,"/",descriptor,"_output")
  
  dir_gifs = paste0(dir_output,"/plots")
  if(!dir.exists(dir_gifs)){
    dir.create(paste0(dir_gifs))
    message("Created directory at ", dir_gifs)
  }
  
  if(mat.use != "assignment"){
    loadings = readRDS(paste0(dir_spatial,"/deconvolution_output_",descriptor,".RDS"))
    cell_types = colnames(loadings[[1]])
  } else {
    assignments = readRDS(paste0(dir_spatial,"/deconvolution_max_factor_assignment_",descriptor,".RDS"))
    cell_types = levels(assignments)
  }
  
  if(!is.null(filter) & mat.use != "assignment"){
    loadings[loadings < filter] = 0
    loadings[loadings >= filter] = 1
    descriptor = paste0(descriptor, "_filter_",filter)
  }
  
  
  grDevices::palette(viridis::viridis(option="A",n=50,direction = -1))
  if(is.null(cell.types.plot)){
    cell.types.plot = cell_types
  } else {
    cell.types.plot = intersect(cell.types.plot, cell_types)
  }  
  
  region_mesh = lapply(subregions.plot, function(x){try(ccf_2017_mesh(x), silent = T)})
  names(region_mesh) = subregions.plot
  region_mesh = region_mesh[sapply(region_mesh, function(x){!any(class(x) == "try-error")})]
  
  for(subregion in names(region_mesh)){
    region_mesh[[subregion]]$material$alpha = .1
    region_mesh[[subregion]]$material$color = "gray"
    
    for(cell_type in cell.types.plot){
      cell_type_consistent = sub("/", ".",sub(" ", "_",cell_type))
      if(mat.use != "assignment"){
        colors_view = as.numeric(cut(loadings[[mat.use]][,cell_type],breaks=50))
        try(rgl.close(), silent = TRUE)
        open3d(windowRect = c(0,0, dims[1], dims[2]));
        plot3d(coords[,1],coords[,2],coords[,3],col = colors_view,aspect=c(67,41,58),xlab="Anterior-Posterior",ylab="Inferior-Superior",zlab="Left-Right",size=5, add = TRUE)
        shade3d(region_mesh[[subregion]])
        decorate3d(xlab = colnames(coords)[1], ylab = colnames(coords)[2],zlab = colnames(coords)[3], box = FALSE, axes = FALSE)
        axes3d(c("x--","y--","z--"))#axes3d(c("x--","y--","z--"))
        movie3d(spin3d(axis = c(0, 0, 1)), duration = 20, movie = paste0(dir_gifs,"/", region, "_", subregion,"_",cell_type_consistent,"_spatial_summary_",descriptor))
      } else {
        colors_view = (assignments == cell_type)*20 + 1
        try(rgl.close(), silent = TRUE)
        open3d(windowRect = c(0,0, dims[1], dims[2]));
        plot3d(coords[,1],coords[,2],coords[,3],col = colors_view,aspect=c(67,41,58),xlab="Anterior-Posterior",ylab="Inferior-Superior",zlab="Left-Right",size=1, type = "p", add = TRUE)
        shade3d(region_mesh[[subregion]])
        decorate3d(xlab = colnames(coords)[1], ylab = colnames(coords)[2],zlab = colnames(coords)[3], box = FALSE, axes = FALSE)
        axes3d(c("x--","y--","z--"))#axes3d(c("x--","y--","z--"))
        movie3d(spin3d(axis = c(0, 0, 1)), duration = 20, movie = paste0(dir_gifs,"/", region, "_", subregion,"_",cell_type_consistent,"_spatial_summary_",descriptor))
      }
      
    }
  }
}

wasserstein_test  = function(
    filepath,
    region,
    spatial.data.name,
    n_samples = 200,
    rand.seed = 123,
    clusters.from.atlas = TRUE,
    naive.clusters = FALSE,
    cell.size = FALSE,
    mat.use = "proportions"#raw, proportions, or assignment
){
  set.seed(rand.seed)
  library(transport)
  
  descriptor = as.character(rand.seed)
  if(clusters.from.atlas){
    descriptor = paste0(descriptor, "_object_clusters")
  }
  if(naive.clusters){
    descriptor = paste0(descriptor, "_naive")
    spatial.data.name = paste0(spatial.data.name, "_naive")
  }
    
  if(cell.size){
     descriptor = paste0(descriptor, "_size_scaled")
  }
  
  dir_spatial = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  
  coords = readRDS(paste0(dir_spatial,"/",spatial.data.name,"_coords_qc_",descriptor,".RDS"))
  
  dir_output = paste0(dir_spatial,"/",descriptor,"_output")

  deconv_out = readRDS(paste0(dir_spatial,"/deconvolution_output_",descriptor,".RDS"))
  loadings = deconv_out[[mat.use]]
  loadings = loadings[, colSums(loadings) != 0 & colnames(loadings)!=""]
  raw_deconv_vals = as.vector(loadings)

  bootstrapped_cell_types = sapply(1:n_samples, function(x){vals = sample(raw_deconv_vals, nrow(loadings));return(vals/sum(vals))})
  normalized_deconv = apply(loadings, MARGIN = 2, function(x){x/sum(x)})

  dist_list = list()
  for(var_1 in colnames(loadings)){
    mass_1 =  normalized_deconv[,var_1]
    distribution_1 = wpp(coords,mass_1)
    dist_list[[var_1]] = vector(length = n_samples)
    for(i in 1:n_samples){
      mass_2 =  bootstrapped_cell_types[,i]
      distribution_2 = wpp(coords,mass_2)
      dist_list[[var_1]][i] = wasserstein(distribution_1, distribution_2)
    }
  }
  
  param_list = lapply(dist_list, function(x){return(c(mean(x), sd(x)))})
  wasserstein_distance_mat <- readRDS(paste0(dir_output,"/wasserstein_distance_mat_",descriptor,".RDS"))

  t_score_mat = sapply(1:length(param_list), function(x){
    (wasserstein_distance_mat[,x]-param_list[[x]][1])/param_list[[x]][2]
  })
  colnames(t_score_mat) = rownames(t_score_mat)
  p_val_mat = sapply(1:length(param_list), function(x){
    p = pt(t_score_mat[,x], n_samples-1)
    return(2*(pmin(p, 1-p)))
  })
  colnames(p_val_mat) = rownames(p_val_mat)
  log_p_val_mat = log1p(p_val_mat)
}
  
