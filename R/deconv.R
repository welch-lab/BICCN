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
  coords = readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name,"/",spatial.data.name,"_coords.RDS"))
  
  minmax = apply(coords, MARGIN = 2, function(x){range(x)})
  ranges = apply(minmax, MARGIN = 2, function(y){return(y[2]-y[1])})
  
  p1 = ggplot(coords, aes_string(x = colnames(coords)[1], y = colnames(coords)[2]))+
    geom_tile() +
    coord_fixed(ratio = 1)+
    geom_tile(aes(alpha = .1)) +
    coord_fixed(ratio = 1) +
    theme_minimal_grid() +
    scale_y_continuous(minor_breaks = seq(minmax[1,2] , minmax[2,2], 1), breaks = seq(minmax[1,2] , minmax[2,2], 5)) +
    scale_x_continuous(minor_breaks = seq(minmax[1,1] , minmax[2,1], 1), breaks = seq(minmax[1,1] , minmax[2,1], 5)) +
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
    scale_y_continuous(minor_breaks = seq(minmax[1,3] , minmax[2,3], 1), breaks = seq(minmax[1,3] , minmax[2,3], 5)) +
    scale_x_continuous(minor_breaks = seq(minmax[1,1] , minmax[2,1], 1), breaks = seq(minmax[1,1] , minmax[2,1], 5)) +
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
    scale_y_continuous(minor_breaks = seq(minmax[1,3] , minmax[2,3], 1), breaks = seq(minmax[1,3] , minmax[2,3], 5)) +
    scale_x_continuous(minor_breaks = seq(minmax[1,2] , minmax[2,2], 1), breaks = seq(minmax[1,2] , minmax[2,2], 5)) +
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
  if(save.plots){
    if(!dir.exists(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots"))){
      dir.create(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/plots"))
    }
    ggplot2::ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots/x_y_coord_reference.PNG"),
                    p1,
                    width = 500,
                    height = 500,
                    units = "px")
    ggplot2::ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots/x_z_coord_reference.PNG"),
                    p2,
                    width = 500,
                    height = 500,
                    units = "px")
    ggplot2::ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots/y_z_coord_reference.PNG"),
                    p3,
                    width = 500,
                    height = 500,
                    units = "px")
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
  coords = readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name,"/",spatial.data.name,"_coords.RDS"))
  spatial.data = readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name,"/",spatial.data.name,"_exp.RDS"))
  for(i in 1:ncol(coords)){
    subset.specs[[i]][is.nan(subset.specs[[i]])] = range(coords[i,])[is.nan(subset.specs[[i]])]
    coords = coords[coords[i,] >= subset.specs[[i]][1] & coords[i,] <= subset.specs[[i]][2], ]
  }
  if(nrow(coords) == 0){
    stop("No samples selected -- provide a new range of coordinates.")
  } else {
    message(paste0("Sample subset to ", nrow(coords), " samples."))
  }
  spatial.data = spatial.data[, rownames(coords)]
  if(!is.null(out.filepath)){
    saveRDS(spatial.data, paste0(out.filepath, "/", spatial.data.name,"_exp.RDS"))
    saveRDS(coords, paste0(out.filepath, "/", spatial.data.name,"_coords.RDS"))
  } else {
    saveRDS(spatial.data, paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name,"/",spatial.data.name,"_exp.RDS"))
    saveRDS(coords, paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name,"/",spatial.data.name,"_coords.RDS"))
  }
}

#' Complete the deconvolution, from data preprocessing through feature selection,
#' gene signature calculation, NNLS, and post-processing of loadings
#'
#' @param filepath Path to directory within which the atlas structure was generated
#' @param region The anatomical region with which the analysis is associated
#' @param spatial.data.name A string, the name of the spatial dataset
#' @param slide.seq A logical, if the spatial data was collected with slide-seq
#'  overriding removal of genes from the spatial data
#' @param z The number of standard deviations above the mean number of missing
#'   values at which to remove a gene from the spatial data
#' @param n.umi.thresh An integer, a minimum number of counts on selected
#'   spatial samples for selected genes
#' @param n.cells The maximum number of cells to select from each dataset for
#'   each cluster for use in gene selection and signature calculation
#' @param deconv.gene.num The number of genes to utilize for gene signature
#'   calculation
#' @param gene.num.tol The maximum difference between deconv.gene.num and the
#'   true number of genes used
#' @param lambda Regularization parameter. Larger values penalize dataset-specific effects more
#'   strongly (ie. alignment should increase as lambda increases)
#' @param thresh Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh.
#'   (default 1e-6)
#' @param max.iters Maximum number of block coordinate descent iterations to perform.
#' @param nrep Number of restarts to perform
#' @param rand.seed Random seed to allow reproducible results.
#' @param annotation.level A value, (1,2,3), corresponding to granularity of
#'   known resolution to use in the event that reference annotations are utilized
#' @param known.annotations An object corresponding to cluster labels for single
#'   cells or a filepath to a spreadsheet containing multiple levels of provided
#'   annotations
#' @param print.obj Print objective function values after convergence
#' @param verbose Print progress bar/messages (TRUE by default)
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
deconvolve_spatial = function(filepath,
                              region,
                              spatial.data.name,
                              slide.seq = FALSE,
                              z = 1,
                              n.umi.thresh = 150,
                              n.cells = 500,
                              deconv.gene.num = 2000,
                              gene.num.tol = 50,
                              lambda = 1,
                              thresh = 1e-8,
                              max.iters = 100,
                              nrep = 1,
                              rand.seed = 1,
                              annotation.level = 3,
                              known.annotations = "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Base_Reference_Files/Reference_Annotations.RDS",
                              print.obj = FALSE,
                              verbose = TRUE){
  message("Loading Data")
  object_path = paste0(filepath,"/", region, "/Analysis1_", region, "/onlineINMF_",region, "_object.RDS" )
  object = readRDS(object_path)


  if(!is.null(known.annotations)){
    if(is.character(known.annotations)){
      master = readRDS(known.annotations)
      if (annotation.level == 1){
        clusters = master$Level1
      }
      if (annotation.level == 2){
        clusters = master$Class
      }
      else {
        clusters = master$Type
      }
      clusters = as.factor(clusters)
      names(clusters) = master$Cell_Barcodes
    } else {
      clusters = known.annotations
    }

  } else {
    #prelim, just for testing until april says something
    clusters = c()
    for(analysis_num in c(2,4,5)){
      analysis_results = readRDS(paste0(filepath,"/",  region, "/Analysis", analysis_num, "_", region, "/Analysis", analysis_num, "_", region,"_Results_Table.RDS"))
      analysis_clusters = as.character(analysis_results$lowRannotations)
      names(analysis_clusters) = analysis_results$Barcode
      clusters = c(clusters, analysis_clusters)
    }
    clusters = clusters[clusters != ""]
    clusters = as.factor(clusters)
  }

  h5_files = sapply(1:length(object@norm.data), function(i){object@h5file.info[[i]]$file.path})
  rna_files = grep(paste0("_(sc10Xv3_|smartseq_|sn10Xv3_|sc10Xv2_)"), h5_files, value = TRUE)

  liger_genes = lapply(rna_files, function(i){
    rhdf5::h5read(i, "/matrix/features")[[1]] #change, extract from H5
  })
  liger_cells = lapply(rna_files, function(i){
    rhdf5::h5read(i, "/matrix/barcodes")#change, extract from H5
  })

  message("Preprocessing for gene selection")

  norm.data = lapply(1:length(rna_files), function(i){
    n = rna_files[i]
    rliger:::Matrix.column_norm(Matrix::sparseMatrix(
      dims = c(length(liger_genes[[i]]), length(liger_cells[[i]])),
      i = as.numeric(rhdf5::h5read(n, "/matrix/indices")+1),
      p = as.numeric(rhdf5::h5read(n, "/matrix/indptr")),
      x = as.numeric(rhdf5::h5read(n, "/matrix/data"))
    ))
  })

  spatial.data = readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name,"/",spatial.data.name,"_exp.RDS"))

  if(!slide.seq){
    spatial.data[spatial.data == -1] = NA
    genes_NA = apply(spatial.data, MARGIN = 1, function(x){sum(is.na(x))})
    mean_genes_NA = mean(genes_NA)
    genes_use = rownames(spatial.data)[genes_NA < (mean_genes_NA + z * mean_genes_NA)]
    spatial.data = spatial.data[genes_use,]
    #simplest way to handle this
    spatial.data[is.na(spatial.data)] = 0
  }

  annotated_cells = intersect(names(clusters),Reduce(union, liger_cells))

  shared_genes = intersect(rownames(spatial.data), Reduce(intersect, liger_genes))

  clusters = clusters[names(clusters) %in% annotated_cells]
  clusters = droplevels(clusters)

  freq_cells = table(clusters)
  freq_cells = freq_cells[names(freq_cells) != ""]

  sample.cells = Reduce(c, lapply(names(freq_cells), function(cell_type){
    Reduce(c, lapply(1:length(norm.data), function(i){
      cells = intersect(names(clusters[clusters == cell_type]), liger_cells[[i]])
      if(length(cells)> 0){
        return(sample(cells, min(length(cells), n.cells), replace =TRUE))
      } else {
        return(c())
      }
    }))
  }))
  norm.data = lapply(1:length(norm.data), function(i){
    rownames(norm.data[[i]]) = liger_genes[[i]]
    colnames(norm.data[[i]]) = liger_cells[[i]]
    gene_means = rhdf5::h5read(rna_files[i], "gene_means")[liger_genes[[i]] %in% shared_genes]
    gene_sum_sq = rhdf5::h5read(rna_files[i], "gene_sum_sq")[liger_genes[[i]] %in% shared_genes]
    root_mean_sum_sq = sqrt(gene_sum_sq/(ncol(norm.data[[i]])-1))
    norm.data[[i]]= sweep(norm.data[[i]][liger_genes[[i]] %in% shared_genes, liger_cells[[i]] %in% sample.cells], 1, root_mean_sum_sq, "/") #liger_cells[[i]] %in% sample.cells
    norm.data[[i]][is.na( norm.data[[i]])] = 0
    norm.data[[i]][ norm.data[[i]] == Inf] = 0
    return(norm.data[[i]])
  })
  norm.data = norm.data[!sapply(norm.data, function(x){length(x) == 0})]

  message("Selecting genes with the KW test")
  chisq_list = list()
  for(i in 1:length(norm.data)){
    chisq_list[[i]] = matrixTests::row_kruskalwallis(as.matrix(norm.data[[i]]),as.vector(clusters[names(clusters) %in%colnames(norm.data[[i]])]))$statistic
    names(chisq_list[[i]]) = rownames(norm.data[[i]])
  }


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

  saveRDS(list(chisq_vals = chisq_list, genes_used = gene_vec), paste0(filepath, "/",  region, "/", region,"_Deconvolution_Output/gene_selection_output.RDS"))

  if(slide.seq){
    spatial.data = spatial.data[rownames(spatial.data) %in% gene_vec, ]
    spatial.data = spatial.data[,Matrix::colSums(spatial.data) > n.umi.thresh]
  }

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
    set.seed(seed = rand.seed + i - 1)
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
      best_seed <- rand.seed + i - 1
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

  saveRDS(out, paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/gene_signature_output.RDS"))

  dir_new = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  if(!dir.exists(dir_new)){
    dir.create(dir_new)
  }

  message("Deconvolving spatial data")
  spatial.data = t(scale(t(as.matrix(spatial.data[gene_vec,])), center = FALSE))
  spatial.data[spatial.data < 0 ] = 0
  spatial.data[is.nan(spatial.data)] = 0
  deconv_h = t(rliger:::solveNNLS(t(W),spatial.data))
  colnames(deconv_h) = clust_levels
  deconv_frac = t(apply(deconv_h, MARGIN = 1, function(x){x/sum(x)}))
  rownames(deconv_frac) = rownames(deconv_h) = colnames(spatial.data)
  deconv_frac[is.nan(deconv_frac)] = 0
  saveRDS(list(raw = deconv_h, proportions = deconv_frac), paste0(dir_new, "/deconvolution_output.RDS"))
}


#' Deconvolve new spatial data from generated gene signatures
#'
#' @param filepath Path to directory within which the atlas structure was generated
#' @param region The anatomical region with which the analysis is associated
#' @param spatial.data.name A string, the name of the spatial dataset
#' @param slide.seq A logical, if the spatial data was collected with slide-seq
#'  overriding removal of genes from the spatial data
#' @param z The number of standard deviations above the mean number of missing
#'   values at which to remove a gene from the spatial data
#' @param n.umi.thresh An integer, a minimum number of counts on selected
#'   spatial samples for selected genes
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

deconvolve_new_data = function(filepath,
                    region,
                    spatial.data.name,
                    slide.seq = FALSE,
                    z = 1,
                    n.umi.thresh = 150){
  spatial.data = readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name,"/",spatial.data.name,"_exp.RDS"))

  deconv_out= readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/gene_signature_output.RDS"))
  W = deconv_out$W
  clust_levels = colnames(deconv_out$H[[1]])
  gene_vec = readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/gene_selection_output.RDS"))[[2]]


  if(!slide.seq){
    spatial.data[spatial.data == -1] = NA
    genes_NA = apply(spatial.data, MARGIN = 1, function(x){sum(is.na(x))})
    mean_genes_NA = mean(genes_NA)
    genes_use = rownames(spatial.data)[genes_NA < (mean_genes_NA + z * mean_genes_NA)]
    spatial.data = spatial.data[genes_use,]
    #simplest way to handle this
    spatial.data[is.na(spatial.data)] = 0
  }

  shared_genes = intersect(rownames(spatial.data), gene_vec)

  if(slide.seq){
    spatial.data = spatial.data[rownames(spatial.data) %in% gene_vec, ]
    spatial.data = spatial.data[,Matrix::colSums(spatial.data) > n.umi.thresh]
  }

  W = W[,gene_vec %in% shared_genes]
  dir_new = paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name)
  if(!dir.exists(dir_new)){
    dir.create(dir_new)
  }

  spatial.data = t(scale(t(as.matrix(spatial.data[shared_genes,])), center = FALSE))
  spatial.data[spatial.data < 0 ] = 0
  spatial.data[is.nan(spatial.data)] = 0
  deconv_h = t(rliger:::solveNNLS(t(W),spatial.data))
  colnames(deconv_h) = clust_levels
  deconv_frac = t(apply(deconv_h, MARGIN = 1, function(x){x/sum(x)}))
  rownames(deconv_frac) = rownames(deconv_h) = colnames(spatial.data)
  deconv_frac[is.nan(deconv_frac)] = 0
  saveRDS(list(raw = deconv_h, proportions = deconv_frac), paste0(dir_new, "/deconvolution_output.RDS"))
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
  spatial.data.name){
  proportions = readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name,"/deconvolution_output.RDS"))[[2]]
  cell_types = colnames(proportions)
  max = as.factor(apply(proportions, MARGIN = 1, function(x){cell_types[which.max(x)]}))
  names(max) = rownames(proportions)
  saveRDS(max, paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name,"/deconvolution_max_factor_assignment.RDS"))
}


#' Convert single cell spatial modalities into voxels by combining samples at
#' a preselected resolution
#'
#' @param filepath Path to directory within which the atlas structure was generated
#' @param region A string corresponding to the name of an anatomical region
#' @param spatial.data.name A string, the name of the spatial dataset
#' @param voxel.size A numeric, corresponding to the side-length of desired voxels
#' @param out.filepath A string corresponding to the directory to which to save generated data
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
  out.filepath,
  verbose = TRUE
){
  spatial.data = readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name,"/",spatial.data.name,"_exp.RDS"))
  coords = readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name,"/",spatial.data.name,"_coords.RDS"))


  minmax = apply(coords, MARGIN = 2, function(x){round(range(x),-1)})
  ranges = apply(minmax, MARGIN = 2, function(y){return(y[2]-y[1])})

  coords_adjusted = sapply(1:3, function(i){
    coords[,i] = coords[,i] - minmax[1, i]
    coords[,i] = round(coords[,i]/voxel.size)
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

  saveRDS(voxel_exp, paste0(out.filepath,"/",region,"_generated_voxel_exp.RDS"))
  saveRDS(coords_out, paste0(out.filepath,"/",region,"_generated_voxel_coords.RDS"))
  saveRDS(voxel_list,paste0(out.filepath,"/",region,"_generated_voxels_to_samples.RDS"))
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
  mat.use = "proportions",#raw, proportions, or assignment
  cell.types.plot = NULL,
  dims = c(500, 500)
){
  library(rgl)
  if(mat.use != "assignment"){
    loadings = readRDS(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/deconvolution_output.RDS"))
    cell_types = colnames(loadings[[1]])
  } else {
    assignments = readRDS(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/deconvolution_max_factor_assignment.RDS"))
    cell_types = levels(assignments)
  }
  if(!dir.exists(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/gifs"))){
    dir.create(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/gifs"))
  }
  grDevices::palette(viridis::viridis(option="A",n=50,direction = -1))
  if(is.null(cell.types.plot)){
    cell.types.plot = cell_types
  } else {
    cell.types.plot = intersect(cell.types.plot, cell_types)
  }
  coords = readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name,"/",spatial.data.name,"_coords.RDS"))


  for(cell_type in cell.types.plot){
    if(mat.use != "assignment"){
      colors_view = as.numeric(cut(loadings[[mat.use]][,cell_type],breaks=50))
      try(rgl.close(), silent = TRUE)
      open3d(windowRect = c(0,0, dims[1], dims[2]));
      plot3d(coords[,1],coords[,2],coords[,3],col = colors_view,aspect=c(67,41,58),xlab="Anterior-Posterior",ylab="Inferior-Superior",zlab="Left-Right",size=5, add = TRUE)
      decorate3d(xlab = colnames(coords)[1], ylab = colnames(coords)[2],zlab = colnames(coords)[3], box = FALSE, axes = FALSE)
      axes3d(c("x--","y--","z--"))#axes3d(c("x--","y--","z--"))
      movie3d(spin3d(axis = c(0, 0, 1)), duration = 20, movie = paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/gifs/", region, "_",sub("/", ".",sub(" ", "_",cell_type)),"_spatial_summary"))
    } else {
      colors_view = (assignments == cell_type)*20 + 1
      try(rgl.close(), silent = TRUE)
      open3d(windowRect = c(0,0, dims[1], dims[2]));
      plot3d(coords[,1],coords[,2],coords[,3],col = colors_view,aspect=c(67,41,58),xlab="Anterior-Posterior",ylab="Inferior-Superior",zlab="Left-Right",size=1, type = "p", add = TRUE)
      decorate3d(xlab = colnames(coords)[1], ylab = colnames(coords)[2],zlab = colnames(coords)[3], box = FALSE, axes = FALSE)
      axes3d(c("x--","y--","z--"))#axes3d(c("x--","y--","z--"))
      movie3d(spin3d(axis = c(0, 0, 1)), duration = 20, movie = paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/gifs/", region, "_",sub("/", ".",sub(" ", "_",cell_type)),"_spatial_summary"))
    }

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
  use.cell.types = TRUE,
  cell.types.use = NULL,
  genes.use = NULL){
  if(mat.use == "assignment"){
    paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name,"/deconvolution_max_factor_assignment.RDS")    sub_vec = rep(0,nlevels(assignments))
    loadings = Reduce(rbind, lapply(assignments, function(x){subbed_vec = sub_vec; subbed_vec[as.numeric(x)] = 1; return(subbed_vec)}))
    colnames(loadings) = levels(assignments)
  } else {
    loadings = readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name,"/deconvolution_output.RDS"))[[mat.use]]
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
      sub_loadings = loadings[rownames(loadings) %in% as.character(layer.list[[i]]), ]
      for(j in 1:length(cell.types.use)){
        if(type == "mean"){
          cell.type.matrix[i,j] = mean(sub_loadings[ ,cell.types.use[j]])
        } else if(type == "sum"){
          cell.type.matrix[i,j] = sum(sub_loadings[ ,cell.types.use[j]])
        }
      }
    }
    saveRDS(cell.type.matrix, paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/cell_type_layer_summary.RDS"))
  }
  if(!is.null(genes.use)){
    spatial.data = readRDS(spatial.data.file)
    genes.use = intersect(genes.use, rownames(spatial.data))
    spatial.data[is.na(spatial.data)] = 0
    spatial.data = t(scale(t(spatial.data[genes.use,]), center = FALSE))
    spatial.data[spatial.data < 0 ] = 0
    gene.matrix = matrix(0L, nrow = length(layer.list), ncol = length(genes.use))
    rownames(gene.matrix) = names(layer.list)
    colnames(gene.matrix) = genes.use
    for(i in 1:length(layer.list)){
      sub_loadings = spatial.data[ ,colnames(spatial.data) %in% as.character(layer.list[[i]])]
      for(j in 1:length(genes.use)){
        if(type == "mean"){
          gene.matrix[i,j] = mean(sub_loadings[genes.use[j],])
        } else if(type == "sum"){
          gene.matrix[i,j] = sum(sub_loadings[genes.use[j],])
        }
      }
    }
    saveRDS(cell.type.matrix, paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatiaal.data.name,"gene_layer_summary.RDS"))
  }
  if(plot){
    if(!dir.exists(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots"))){
      dir.create(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots"))
    }
    ggplot2::theme_set(theme_cowplot())
    labels.cell.type = expand.grid(rownames(cell.type.matrix), colnames(cell.type.matrix))
    cell.type.df = data.frame(Layers = as.character(labels.cell.type[,1]),
                              Cell_Types = as.character(labels.cell.type[,2]),
                              Values = as.vector(cell.type.matrix))
    #for(i in 1:nrow(cell.type.df)){
    #  cell.type.df[i, 3] = cell.type.matrix[cell.type.df[i, 1],  cell.type.df[i, 2]]
    #}
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

      ggplot2::ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots/cell_type_layer_distribution.PNG"),
                      overall.cell.type.plot,
                      width = 1000,
                      height = 800,
                      units = "px")
    }

    for(i in colnames(cell.type.matrix)){
      cell.type.df.sub = cell.type.df[cell.type.df$Cell_Types == i,]
      by.cell.type.plot = ggplot2::ggplot(cell.type.df.sub, ggplot2::aes(fill = Layers, x = Layers, y = Values)) +
        ggplot2::theme(text = ggplot2::element_text(size = 3),
                       axis.text = ggplot2::element_text(size = 2),
                       legend.position="none") +
        ggplot2::geom_bar(position = "dodge", stat = "identity") +
        ggplot2::xlab("Layer") +
        ggplot2::ylab("Value") +
        ggplot2::ggtitle(paste0("Distribution of ",i, " cells by layer"))
      ggplot2::ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots/",sub("/", ".",sub(" ", "_",i)),"_layer_distribution.PNG"),
                      by.cell.type.plot,
                      width = 500,
                      height = 400,
                      units = "px")

    }
    if(!is.null(genes.use)){
      labels.genes = expand.grid(rownames(gene.matrix), colnames(gene.matrix))
      gene.df = data.frame(Layers = as.character(labels.genes[,1]),
                           Genes = as.character(labels.genes[,2]),
                           Values = as.vector(gene.matrix))
      #for(i in 1:nrow(gene.df)){
      #  gene.df[i, 3] = gene.matrix[gene.df[i, 1],  gene.df[i, 2]]
      #}
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

        ggplot2::ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots/gene_layer_distribution.PNG"),
                        overall.gene.plot,
                        width = 1000,
                        height = 800,
                        units = "px")
      }

      for(i in colnames(gene.matrix)){
        gene.df.sub = gene.df[gene.df$Genes == i,]
        by.gene.plot = ggplot2::ggplot(gene.df.sub, ggplot2::aes(fill = Layers, x = Layers, y = Values)) +
          ggplot2::theme(text = ggplot2::element_text(size = 3),
                         axis.text = ggplot2::element_text(size = 2),
                         legend.position="none") +
          ggplot2::geom_bar(position = "dodge", stat = "identity") +
          ggplot2::xlab("Layer") +
          ggplot2::ylab("Value") +
          ggplot2::ggtitle(paste0("Distribution of ",i, " expression by layer"))
        ggplot2::ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots/",sub("/", ".",sub(" ", "_",i)),"_gene_layer_distribution.PNG"),
                        by.gene.plot,
                        width = 500,
                        height = 400,
                        units = "px")

      }
    }
  }
}

#' Analyze gene signatures, resulting in analysis of similarity between cell
#' types as well as of covariance in cell type distribution
#'
#' @param filepath Path to directory within which the atlas structure was generated
#' @param region A string corresponding to the name of an anatomical region
#' @param spatial.data.name A string, the name of the spatial dataset
#' @param plot A logical, corresponding with if to make heaatmaps and dendrograms
#'   from the generated data
#' @param mat.use A string, either "raw" or "proportions"
#'   corresponding to the raw cell type loadings or the normalized loadings
#'
#'
#' @return nothing
#'
#' @import ggplot2, lsa, ggdendro
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }

analyze_gene_signatures = function(filepath,
                                   region,
                                   spatial.data.name,
                                   plot = FALSE,
                                   mat.use = "proportions"){

  library(ggplot2)

  gene_sigs = readRDS(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/gene_signature_output.RDS"))
  signatures = gene_sigs$W
  rownames(signatures) = colnames(gene_sigs$H[[1]])
  colnames(signatures) = readRDS(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/gene_selection_output.RDS"))[[2]]
  signatures = signatures[rownames(signatures) != "",]

  cos_sim = lsa::cosine(t(signatures))

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

  cos_dist= as.dist(1- cos_sim)
  hierarchical_clust <- hclust(cos_dist, method = "ward.D2")

  dendro_plot = ggdendro::ggdendrogram(hierarchical_clust, rotate = FALSE, size = 2) +
    ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = 4),
                   axis.text = ggplot2::element_text(size = 4)) +
    ggplot2::ggtitle("Hierarchical Clustering", subtitle = "By Cell Type Signature")

  deconv_out = readRDS(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/deconvolution_output.RDS"))
  loadings = deconv_out[[mat.use]]
  loadings = loadings[, colnames(loadings)!=""]

  corr_sim_dist = cor(loadings)

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

  cor_dist= as.dist(1- corr_sim_dist)
  hierarchical_clust_dist <- hclust(cor_dist, method = "ward.D2")

  dendro_dist_plot = ggdendro::ggdendrogram(hierarchical_clust_dist, rotate = FALSE, size = 2) +
    ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = 4),
                   axis.text = ggplot2::element_text(size = 4)) +
    ggplot2::ggtitle("Hierarchical Clustering", subtitle = "By Cell Type Distribution")
  if(plot){
    print(heatmap_plot)
    print(dendro_plot)
    print(heatmap_dist_plot)
    print(dendro_dist_plot)
    if(!dir.exists(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots"))){
      dir.create(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots"))
    }
    ggplot2::ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots/cell_type_signature_heatmap.PNG"),
                    heatmap_plot,
                    width = 1000,
                    height = 800,
                    units = "px")
    ggplot2::ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots/cell_type_signature_dendrogram.PNG"),
                    dendro_plot,
                    width = 500,
                    height = 400,
                    units = "px")
    ggplot2::ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots/cell_type_distribution_heatmap.PNG"),
                    heatmap_dist_plot,
                    width = 1000,
                    height = 800,
                    units = "px")
    ggplot2::ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots/cell_type_distribution_dendrogram.PNG"),
                    dendro_dist_plot,
                    width = 500,
                    height = 400,
                    units = "px")
  }
  saveRDS(list(cos_sim_signature = cos_sim, cor_sim_distribution = corr_sim_dist, dendro_sig = hierarchical_clust, dendro_dist = hierarchical_clust_dist),
          paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/gene_signature_analysis_summary.RDS"))
  detach("package:ggplot2", unload = TRUE)
}


#' Generate 2D plots of cell type and gene distribution, given a plane.
#'
#' @param filepath Path to directory within which the atlas structure was generated
#' @param region A string corresponding to the name of an anatomical region
#' @param spatial.data.name A string, the name of the spatial dataset
#' @param axis A string corresponding to a column name of the coord object, the
#'    axis along which to plot the slice
#' @param idx An integer or vector of integers, corresponding to either what
#'    slice along the given axis to plot or a range along that axis to flatten
#' @param mat.use A string, either "raw", "proportions", or "assignment",
#'    corresponding to the raw cell type loadings, the normalized loadings, or
#'    cell type assignments for single cell spatial modalities
#' @param use.cell.types A logical, if only the cell types provided to
#'    cell.types.use should be summarized, as opposed to all deconvolved.
#' @param cell.types.use A character vector of cell types to summarize
#' @param genes.use A character vector of genes to summarize
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
plot_layer = function(
  filepath,
  region,
  spatial.data.name,
  axis,
  idx,
  mat.use = "proportions",
  use.cell.types = TRUE,
  cell.types.use = NULL,
  genes.use = NULL,
  display.plots = FALSE){

  library(ggplot2)
  if(!dir.exists(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots"))){
    dir.create(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots"))
  }

  if(mat.use == "assignment"){
    assignments = readRDS(paste0(filepath, "/", region, "/",region, "_Deconvolution_Output/",spatial.data.name,"/deconvolution_max_factor_assignment.RDS"))
    sub_vec = rep(0,nlevels(assignments))
    loadings = Reduce(rbind, lapply(assignments, function(x){subbed_vec = sub_vec; subbed_vec[as.numeric(x)] = 1; return(subbed_vec)}))
    colnames(loadings) = levels(assignments)
    rownames(loadings) = names(assignments)
  } else {
    loadings = readRDS(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/deconvolution_output.RDS"))[[mat.use]]
  }

  theme_set(cowplot::theme_cowplot())

  coords = readRDS(paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/",spatial.data.name,"/",spatial.data.name,"_coords.RDS"))

  if(length(idx) > 1){
    coords = coords[coords[,axis] >= idx[1] & coords[,axis] <= idx[2],
                    setdiff(colnames(coords),axis)]
  } else {
    coords = coords[coords[,axis] == idx, setdiff(colnames(coords),axis)]
  }
  if(use.cell.types){
    if(!is.null(cell.types.use)){
      cell.types.use = intersect(cell.types.use, colnames(loadings))
    } else {
      cell.types.use = colnames(loadings)
    }
    cell.types.use = cell.types.use[cell.types.use != ""]
    loadings = loadings[rownames(loadings) %in% rownames(coords), cell.types.use]

    colnames(loadings) = sub("/",".",sub(" ", "_", cell.types.use))

    coords = coords[rownames(loadings),]
    plotting_df = as.data.frame(cbind(coords, loadings))
    for(i in 1:ncol(loadings)){
      slice_plot = ggplot(plotting_df, aes_string(x = colnames(plotting_df)[1],
                                                  y = colnames(plotting_df)[2],
                                                  fill = colnames(loadings)[i])) +
        geom_tile() +
        coord_fixed(ratio = 1) +
        viridis::scale_fill_viridis() +
        ggtitle(paste0("Distribution of ",cell.types.use[i]),
                subtitle = paste0("In ", axis, " slice ", idx)) +
        theme(legend.title = element_blank(),
              text = element_text(size = 8),
              axis.text = element_text(size = 5))
      if(display.plots){
        print(slice_plot)
      }
      if(length(idx)>1){
        ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots/",
                      colnames(loadings)[i],"_",axis,"_",idx[1],"_to_",idx[2],".PNG"),
               slice_plot,
               width = 1000,
               height = 800,
               units = "px")
      } else {
        ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots/",
                      colnames(loadings)[i],"_",axis,"_",idx,".PNG"),
               slice_plot,
               width = 1000,
               height = 800,
               units = "px")
      }

    }
  }
  if(!is.null(genes.use)){
    spatial.data = readRDS(spatial.data.file)
    genes.use = intersect(genes.use, rownames(spatial.data))
    spatial.data[is.na(spatial.data)] = 0
    spatial.data = scale(t(spatial.data[genes.use,colnames(spatial.data) %in% rownames(coords)]), center = FALSE)
    spatial.data[spatial.data < 0 ] = 0
    plotting_df = as.data.frame(cbind(coords, spatial.data))
    for(i in 1:ncol(spatial.data)){
      slice_plot = ggplot(plotting_df, aes_string(x = colnames(plotting_df)[1],
                                                  y = colnames(plotting_df)[2],
                                                  fill = colnames(spatial.data)[i])) +
        geom_tile() +
        coord_fixed(ratio = 1) +
        viridis::scale_fill_viridis() +
        ggtitle(paste0("Distribution of ",genes.use[i]),
                subtitle = paste0("In ", axis, " slice ", idx)) +
        theme(legend.title = element_blank(),
              text = element_text(size = 8),
              axis.text = element_text(size = 5))
      if(display.plots){
        print(slice_plot)
      }
      if(length(idx)>1){
        ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots/",
                      genes.use[i],"_",axis,"_",idx[1],"_to_",idx[2],".PNG"),
               slice_plot,
               width = 1000,
               height = 800,
               units = "px")
      } else {
        ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/",spatial.data.name,"/plots/",
                      genes.use[i],"_",axis,"_",idx,".PNG"),
               slice_plot,
               width = 1000,
               height = 800,
               units = "px")
      }
    }
  }
  detach("package:ggplot2", unload = TRUE)
}
