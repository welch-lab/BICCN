library(magrittr)
library(stringr)


'%notin%' = Negate('%in%')


create.directories = function(region = "X", desired.filepath = "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/"){
  #Initialize main region directory
  #Add a slash at the end of directory if one is not provided
  if(str_sub(desired.filepath, start = -1L, end = -1L) != "/"){
    desired.filepath = paste0(desired.filepath, "/")
  }
  main_directory = paste0(desired.filepath, region)
  dir.create(main_directory)
  #Create subdirectories
  for (i in 1:5){
    sub_directory_name = paste0(main_directory, "/Analysis", i, "_", region)
    sub_images_name = paste0(main_directory, "/Analysis", i, "_", region, "/Images")
    dir.create(sub_directory_name)
    dir.create(sub_images_name)
  }
  loom_directory = paste0(main_directory, "/", region,"_Loom_Directory")
  dir.create(loom_directory)
  deconv_directory = paste0(main_directory, "/", region,"_Deconvolution_Output")
  dir.create(deconv_directory)
}
################## QC function
#Input: List of filenames and types
#Types should be listed as "sc10Xv3", "smartseq", "atac", "meth", "sn10Xv3", "sc10Xv2"
#If more than one of each type of file exists, they should be listed as "sc10Xv3_1", "sc10Xv3_2", etc.
#output: All files will be out put as .RDS files. QC metrics will be stored in Analysis#_region_Overall_qc.RDS
library(sjmisc)
library(stringr)
library(rliger)
library(dplyr)
library(edgeR)

# Base QC table is available at "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Base_Reference_Files/QC_table.RDS"
# @slot analysis_num Determines directory and output file naming conventions, also ensures that Analysis one will ignore methylation Data

apply_qc = function(filenames, region, analysis_num , qc_table_path, filepath_cytoplasmic = "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Base_Reference_Files/Cytoplasmic_genes.csv",  filepath = "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/", verbose = TRUE){
  if(str_sub(filepath, start = -1L, end = -1L) != "/"){
    filepath = paste0(filepath, "/")
  }
  '%notin%' = Negate('%in%')
  qc_summary = setNames(data.frame(matrix(ncol = 13, nrow = 0)), c("Datatype", "Edition", "OriginalDimensions", "AnalysisDimensions", "FailnUMI", "nUMIThreshold", "FailMito", "MitoThreshold", "FailCyto", "CytoThreshold", "Lost for Gene Counts", "Gene Count Cutoff", "FinalDimensions"))
  for (file.listed in 1:length(filenames)){
    #i.e. for each file 
    # For Analysis != 1, we will need to filter the matrices on the basis on their results
    #i.e. Analysis one should be all cells, not including methylation
    #Analysis 2 should just be non-neuronal cells
    #Analysis 3 should be just neuronal cells + meth
    #Analysis 4 should be just inhib
    #Analysis 5 should be just excit
    #Read in each file
    working_file = readRDS(filenames[[file.listed]])
    
    #Determine type, while accomodating for the fact that some types might have multiple files
    data.type = names(filenames[file.listed])
    if(verbose){
      message("Processing ", data.type)
    }
    
    if(str_contains(data.type,"_")){
      edition = str_split(data.type[[1]], "_")[[1]][2]
      data.type = str_split(data.type[[1]], "_")[[1]][1]
    }else{ edition = 0 }
    #Need to check data for duplicate cellnames
    if (length(unique(colnames(working_file))) != length(colnames(working_file))){
      colnames(working_file) %<>% make.unique()
      warning("Duplicated cell names, making unique")
      
    }
    
    #Determine cell labels 
    if (analysis_num == 1){
      use.cells = colnames(working_file)
    }
    if(analysis_num == 2){ 
      if(data.type != "meth"){
        results_filename = paste0(filepath,region, "/Analysis1_", region, "/Analysis1_", region, "_Results_Table.RDS")
        results = readRDS(results_filename)
        results = filter(results, results$highRAnnotations =="NonN")
        results = subset(results, results$Barcode %in% colnames(working_file))
        use.cells = results$Barcode
      } else {use.cells = colnames(working_file)}
    }
    if(analysis_num == 3){
      if(data.type != "meth"){
        results_filename = paste0(filepath,region, "/Analysis1_", region, "/Analysis1_", region, "_Results_Table.RDS")
        results = readRDS(results_filename)
        results = filter(results, results$highRAnnotations =="Neu")
        results = subset(results, results$Barcode %in% colnames(working_file))
        use.cells = results$Barcode
      } else {use.cells = colnames(working_file)}
    }
    if(analysis_num == 4){
      results_filename = paste0(filepath,region, "/Analysis3_", region, "/Analysis3_", region, "_Results_Table.RDS")
      results = readRDS(results_filename)
      results = filter(results, results$highRAnnotations =="Inh")
      results = subset(results, results$Barcode %in% colnames(working_file))
      use.cells = results$Barcode
    }
    if(analysis_num == 5){
      results_filename = paste0(filepath,region, "/Analysis3_", region, "/Analysis3_", region, "_Results_Table.RDS")
      results = readRDS(results_filename)
      results = filter(results, results$highRAnnotations =="Exc")
      results = subset(results, results$Barcode %in% colnames(working_file))
      use.cells = results$Barcode
    }
    if (length(use.cells) < 30){next}
    ######################### Pre-process 
    #Change the names of the Rik genes. This needs to be done for the sc10Xv3 and smartseq datasets
    if (data.type == "smartseq"){
      #With the new smartseq data, we no longer need the gene converter
      
      #rik_converter = readRDS("/nfs/turbo/umms-welchjd/BRAIN_initiative/Ensembl_conversion/ConvertRik.Tenx.SMART.RDS")
      #Need to raise a warning if there are genes in the smartseq matrix that are not in the conversion table
      #refined_genes = rownames(working_file)
      #potential_outliers = subset(refined_genes, refined_genes %notin% rik_converter$Og_Barcode)
      #if (length(potential_outliers >= 1)) { 
      #  warning("There are genes in the matrix that are not covered in the conversion table. Please update the conversion table to address these genes and then rerun the function.", immediate. = T)
      #} 
      #rowname_new = rik_converter$New_Barcode[match(rownames(working_file), rik_converter$Og_Barcode)]
      #rownames(working_file) = rowname_new 
      qc.matrix = working_file[,use.cells]  #only grab applicable cell populations
      remaining_GeneCounts = remaining_nUMI = remaining_mito = remaining_cyto = colnames(qc.matrix)
      nUMI_cutoff = mito_cutoff = cytoplasmic_cutoff = GeneCount_cutoff = NA
      before_subset = dim(working_file)[[2]] #Gets you original dimensions of matrix
      working_file = working_file[,use.cells] #Gets you a matrix subset for the appropriate cell population
      after_subset = dim(working_file)[[2]]  #Gets you dimensions of matrix after subsetting for the appropriate cell population
      cells_after_subset = colnames(working_file) #Gets you the cell names present after subsetting for cell populations
      
    }
    
    if (data.type == "atac"){
      qc_table = readRDS(qc_table_path)
      #Find the appropriate cutoff to apply
      if (data.type == "atac"){ qc_stats = filter(qc_table, qc_table$DataType == "atac")}
      nUMI_cutoff = qc_stats$nUMI
      mito_cutoff = qc_stats$mito
      ligs = createLiger(list(qc_mat = working_file))
      celldata = ligs@cell.data
      celldata$Mito = getProportionMito(ligs)
      celldata = filter(celldata, celldata$nUMI >= as.numeric(nUMI_cutoff))
      remaining_nUMI = rownames(celldata)
      celldata = filter(celldata, celldata$Mito <= as.numeric(mito_cutoff))
      remaining_mito = rownames(celldata)
      before_subset = dim(working_file)[[2]] #Gets you original dimensions of matrix
      working_file = working_file[,use.cells] #Gets you a matrix subset for the appropriate cell population
      after_subset = dim(working_file)[[2]]  #Gets you dimensions of matrix after subsetting for the appropriate cell population
      cells_after_subset = colnames(working_file) #Gets you the cell names present after subsetting for cell populations
      sub_cells = subset(rownames(celldata), rownames(celldata) %in% colnames(working_file))
      qc.matrix = working_file[, sub_cells]
      remaining_GeneCounts = remaining_cyto = colnames(qc.matrix)
      cytoplasmic_cutoff = GeneCount_cutoff = NA
      
    }
    
    if (data.type == "sc10Xv3" | data.type == "sc10Xv2"){  #i.e. originally called tenx
      qc_table = readRDS(qc_table_path)
      #Find the appropriate cutoff to apply 
      qc_stats = filter(qc_table, qc_table$DataType == "sc10Xv3")
      qc_stats = filter(qc_stats, qc_stats$Analysis == analysis_num)
      if (analysis_num == 1 | analysis_num == 2 | analysis_num == 3){
        qc_stats = filter(qc_stats, qc_stats$Region == "ALL")
      }
      if (analysis_num == 4 | analysis_num == 5){
        if (region != "OLF" & region != "CB"){
          qc_stats = filter(qc_stats, qc_stats$Region == "ALL")
        } else {
          qc_stats = filter(qc_stats, qc_stats$Region == region)
        }
      }
      
      GeneCount_cutoff = qc_stats$GeneCounts
      mito_cutoff = qc_stats$mito
      nUMI_cutoff = NA
      cytoplasmic_cutoff = qc_stats$CytoplasmicScore
      ligs = createLiger(list(qc_mat = working_file))
      celldata = ligs@cell.data
      celldata$Mito = getProportionMito(ligs)
      #Calculate cytoplasmic score 
      if (analysis_num == 3 | analysis_num == 4 | analysis_num == 5){
        celldata$Cytoplasmic_Score = cal_qc_score(working_file, qc.gene.FN = filepath_cytoplasmic)
      }
      remaining_nUMI = rownames(celldata)
      celldata = filter(celldata, celldata$nGene >= as.numeric(GeneCount_cutoff))
      remaining_GeneCounts = rownames(celldata)
     # celldata = filter(celldata, celldata$Mito <= as.numeric(mito_cutoff))
      remaining_mito = rownames(celldata)
      if (analysis_num == 3 | analysis_num == 4 | analysis_num == 5){
        celldata = filter(celldata, celldata$Cytoplasmic_Score > as.numeric(cytoplasmic_cutoff))
      }
      remaining_cyto = rownames(celldata)
      before_subset = dim(working_file)[[2]] #Gets you original dimensions of matrix
      working_file = working_file[,use.cells] #Gets you a matrix subset for the appropriate cell population
      after_subset = dim(working_file)[[2]]  #Gets you dimensions of matrix after subsetting for the appropriate cell population
      cells_after_subset = colnames(working_file) #Gets you the cell names present after subsetting for cell populations
      #If datatype is ATAC sn10Xv3 sc10Xv3, filter for QC cells
      sub_cells = subset(rownames(celldata), rownames(celldata) %in% colnames(working_file))
      qc.matrix = working_file[, sub_cells]
    }
    
    if (data.type == "sn10Xv3" | data.type == "sn10Xv2"){
      #Remember, new macosko data is already pre-filtered. No need to apply any QC.
      qc_table = readRDS(qc_table_path)
      qc_stats = filter(qc_table, qc_table$DataType == "sn10Xv3" )
      nUMI_cutoff = qc_stats$nUMI
      mito_cutoff = qc_stats$mito
      cytoplasmic_cutoff = qc_stats$CytoplasmicScore
      ligs = createLiger(list(qc_mat = working_file))
      celldata = ligs@cell.data
      remaining_nUMI = rownames(celldata)
      remaining_cyto = remaining_mito = rownames(celldata)
      before_subset = dim(working_file)[[2]] #Gets you original dimensions of matrix
      use.cells = as.character(use.cells)
      working_file = working_file[,use.cells] #Gets you a matrix subset for the appropriate cell population
      after_subset = dim(working_file)[[2]]  #Gets you dimensions of matrix after subsetting for the appropriate cell population
      cells_after_subset = colnames(working_file) #Gets you the cell names present after subsetting for cell populations
      #If datatype is ATAC sn10Xv3 sc10Xv3, filter for QC cells
      sub_cells = subset(rownames(celldata), rownames(celldata) %in% colnames(working_file))
      qc.matrix = working_file[, sub_cells]
      
    }
    
    if (data.type != "meth"){
      #Count removed Cells
      #Cell kept after discarding cells that did not pass nUMI threshold
      lost_nUMI = length(subset(cells_after_subset, cells_after_subset %notin% remaining_nUMI))
      lost_mito = length(subset(remaining_nUMI, remaining_nUMI %notin% remaining_mito))
      lost_cyto = length(subset(remaining_mito, remaining_mito %notin% remaining_cyto))
      lost_GeneCounts = length(subset(remaining_cyto, remaining_cyto %notin% remaining_GeneCounts))  
      new.dim = dim(qc.matrix)[[2]]
      if (edition == 0){
        save.filename = paste0(filepath, region, "/Analysis", analysis_num , "_", region, "/", region, "_", data.type, "_qc.RDS")
      } else {save.filename = paste0(filepath, region, "/Analysis", analysis_num , "_", region, "/", region, "_" ,data.type,"_", edition, "_qc.RDS") }
      saveRDS(qc.matrix, save.filename)
      newest_qc = c(data.type, edition, before_subset, after_subset, lost_nUMI, nUMI_cutoff, lost_mito, mito_cutoff, lost_cyto, cytoplasmic_cutoff, lost_GeneCounts,  GeneCount_cutoff, new.dim)
      qc_summary = rbind(qc_summary, newest_qc)
    } 
    
    if( data.type == "meth") {
      if (analysis_num == 1 | analysis_num == 2){
        print("No preprocessing for methylation data necessary, it is not used")
      }
      if (analysis_num != 1 & analysis_num != 2){
        
        before_subset = dim(working_file)[[2]] #Gets you original dimensions of matrix
        use.cells = subset(use.cells, use.cells %in% colnames(working_file))
        qc.matrix = working_file[,use.cells] #Gets you a matrix subset for the appropriate cell population
        after_subset = dim(qc.matrix)[[2]]  #Gets you dimensions of matrix after subsetting for the appropriate cell population
        cells_after_subset = colnames(working_file)
        lost_GeneCounts = lost_nUMI = lost_mito = lost_cyto = 0
        nUMI_cutoff = mito_cutoff = cytoplasmic_cutoff = GeneCount_cutoff = NA
        new.dim = dim(qc.matrix)[[2]]
        
        nUMI_cutoff = mito_cutoff = NA
        if (edition == 0){
          save.filename = paste0(filepath, region, "/Analysis", analysis_num , "_", region, "/", region, "_",data.type, "_qc.RDS")
        } else { save.filename = paste0(filepath, region, "/Analysis", analysis_num , "_", region, "/", region, "_" ,data.type,"_", edition, "_qc.RDS") }
        saveRDS(qc.matrix, save.filename)
        newest_qc = c(data.type, edition, before_subset, after_subset, lost_nUMI, nUMI_cutoff, lost_mito, mito_cutoff, lost_cyto, cytoplasmic_cutoff, lost_GeneCounts,  GeneCount_cutoff, new.dim)
        qc_summary = rbind(qc_summary, newest_qc) }}}
  if(verbose){
    message("Done Processing ", data.type, " ", edition)
  }
  
  if (data.type != "meth" & data.type != "atac" & data.type != "smartseq" & data.type != "sc10Xv3" & data.type != "sn10Xv3" & data.type != "sc10Xv2"){
    warning("Uknown Data Type submitted", immediate. = T)
    
  }
  #Output a .csv to the Analysis directory
  csv_filename = paste0(filepath, region, "/Analysis", analysis_num , "_", region, "/", region, "_Overall_qc.csv")
  colnames(qc_summary) = c("Datatype", "Edition", "OriginalDimensions", "AnalysisDimensions", "FailnUMI", "nUMIThreshold", "FailMito", "MitoThreshold", "FailCyto", "CytoThreshold", "Lost for Gene Counts", "Gene Count Cutoff", "FinalDimensions")
  write.csv(qc_summary, csv_filename)
}
#object is a fully processed object, with clusters from either louvain or max factor assignment, assignment is a factor covering at least k of the cells in the object, k is used for nearest neighbors.
transfer_labels = function(object, annotations, k = 20){
  H.norm = object@H.norm
  previous_clusters = object@clusters
  cells = rownames(object@H.norm)
  
  known_samples = H.norm[rownames(H.norm) %in% names(annotations),]
  unknown_samples = H.norm[setdiff(rownames(H.norm), names(annotations)),]
  if(nrow(known_samples) == 0){
    message("No annotations available")
    return(object)
  }
  annotations = droplevels(annotations[rownames(known_samples)])
  out = RANN::nn2(known_samples, unknown_samples, k = k)
  
  levels_annotations <- levels(annotations)
  
  vote = apply(out$nn.idx, MARGIN = 1, function(unknown){
    (1:nlevels(annotations))[which.max(tabulate(match(annotations[unknown], levels_annotations)))]
  })
  names(vote) = rownames(unknown_samples)
  transfered_annotations = factor(c(annotations, vote))
  levels(transfered_annotations) = levels(annotations)
  cluster_comp = sapply(levels(previous_clusters), function(prev_clust){
    sapply(levels(transfered_annotations), function(trans_clust){
      sum(previous_clusters[cells] == prev_clust & transfered_annotations[cells] == trans_clust)
    })
  })
  cluster_assignments = apply(cluster_comp, MARGIN = 2, function(x){
    levels(transfered_annotations)[which.max(x)]})
  levels(object@clusters) = cluster_assignments
  return(object)
}


annotate_by_modality = function(filepath, 
                                region, 
                                analysis_num, 
                                chunk_size,
                                knownAnnotations = "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Base_Reference_Files/Reference_Annotations.RDS",
                                knn = 10,
                                num_genes = 2500, 
                                gene_num_tolerance = 100, 
                                var_thresh_start = 2, 
                                max_var_thresh = 4, 
                                k = 20,
                                lambda = 5,
                                max.epochs = 5,
                                miniBatch_max_iters = 1,
                                miniBatch_size = 5000,
                                seed = 123,
                                verbose = TRUE){
  
  annot = readRDS(knownAnnotations)
  if(analysis_num == 1){
    annotations = annot$Level1
  } else if(analysis_num %in% c(2,3)){
    annotations = annot$Class 
  } else {
    annotations = annot$Type
  }
  annotations = as.factor(annotations)
  names(annotations) = annot$Cell_Barcodes
  
  liger_name = paste0(filepath, "/", region, "/Analysis", analysis_num, "_", region, "/onlineINMF_",region, "_object.RDS" )
  object = readRDS(liger_name)
  dataset = object@cell.data$dataset
  names(dataset) = rownames(object@cell.data)
  annot_dataset = dataset[names(dataset) %in% names(annotations)]
  freq = table(annot_dataset)
  if(0 == sum(freq[grepl("_(atac_)",names(freq))]) * sum(freq[grepl("_(meth_)",names(freq))]) * sum(freq[grepl("_(sc10Xv3_|smartseq_|sn10Xv3_|sc10Xv2_)",names(freq))])){
    message("Annotations missing for at least one modality --  joint annotation transfer being conducted")
    object = transfer_labels(object, annotations, k)
    return(list("all" = object))
  } else {
    rm(object)
    qc_files = list.files(paste0(filepath, "/", region, "/Analysis", analysis_num , "_", region, "/"))
    qc_files = grep(paste0(region,"_(sc10Xv3_|smartseq_|atac_|meth_|sn10Xv3_|sc10Xv2_).*(qc.RDS)"), qc_files, value = TRUE)
    files = list()
    files[["rna"]] = grep(paste0("_(sc10Xv3_|smartseq_|sn10Xv3_|sc10Xv2_)"), qc_files, value = TRUE)
    files[["atac"]] = grep(paste0("_(atac_)"), qc_files, value = TRUE)
    files[["meth"]] = grep(paste0("_(meth_)"), qc_files, value = TRUE)
    files[length(files) == 0] = NULL
    object_list = list()
    
    for(modality in names(files)){
      if(modality != "meth"){
        rhdf5::h5closeAll()
        
        hdf5_files = paste0(filepath, "/", region, "/Analysis", analysis_num , "_", region, "/",modality,"_subanalysis_",gsub(".RDS", ".H5",files[[modality]]))
        for (i in 1:length(hdf5_files)){
          current_matrix = Matrix::Matrix(readRDS(paste0(filepath, "/", region, "/Analysis", analysis_num , "_", region, "/",files[[modality]][i])), sparse = TRUE)
          rhdf5::h5createFile(hdf5_files[[i]])
          rhdf5::h5createGroup(hdf5_files[[i]], "matrix")
          rhdf5::h5write(current_matrix@Dimnames[[2]], file=hdf5_files[[i]], name="matrix/barcodes") # cell barcodes
          rhdf5::h5createGroup(hdf5_files[[i]], file.path("matrix", "features"))
          rhdf5::h5write(current_matrix@Dimnames[[1]], file=hdf5_files[[i]], name="matrix/features/name")
          rhdf5::h5write(dim(current_matrix), file=hdf5_files[[i]], name="matrix/shape")
          rhdf5::h5createDataset(hdf5_files[[i]],"matrix/data",dims=length(current_matrix@x),storage.mode="double", chunk = min(chunk_size, ncol(current_matrix)))
          rhdf5::h5write(current_matrix@x, file=hdf5_files[[i]], name="matrix/data", index = list(1:length(current_matrix@x)))
          rhdf5::h5createDataset(hdf5_files[[i]],"matrix/indices",dims=length(current_matrix@i),storage.mode="integer", chunk = min(chunk_size, ncol(current_matrix)))
          rhdf5::h5write(current_matrix@i, file=hdf5_files[[i]], name="matrix/indices",index = list(1:length(current_matrix@i))) # already zero-indexed.
          rhdf5::h5createDataset(hdf5_files[[i]],"matrix/indptr",dims=length(current_matrix@p),storage.mode="integer", chunk = min(chunk_size, ncol(current_matrix)))
          rhdf5::h5write(current_matrix@p, file=hdf5_files[[i]], name="matrix/indptr",index = list(1:length(current_matrix@p)))
          rhdf5::h5closeAll()
        }
        data_names = gsub("(_qc.RDS)", "",files[[modality]])
        names(hdf5_files) = data_names
        object = createLiger(raw.data = as.list(hdf5_files))
        object = normalize(object, chunk = chunk_size)
        datasets_use = grep("[sc10Xv3_|smartseq_|sc10Xv2_]",files[[modality]])
        object = selectGenes(object, var.thresh = var_thresh_start, datasets.use = datasets_use)
        var_thresh_old = var_thresh_start
        high = max_var_thresh
        low = 0
        while(abs(length(object@var.genes)-num_genes) > gene_num_tolerance){
          if(length(object@var.genes)>num_genes){
            var_thresh_new = (var_thresh_old+high)/2
            low = var_thresh_old
          } else {
            var_thresh_new = (var_thresh_old+low)/2
            high = var_thresh_old
          }
          object = selectGenes(object, var.thresh = var_thresh_new, datasets.use = datasets_use)
          var_thresh_old = var_thresh_new
        }
        message(paste0(length(object@var.genes), " genes found with var.thresh = ",var_thresh_old))
        object = scaleNotCenter(object, chunk = chunk_size)
        object = online_iNMF(object,h5_chunk_size = chunk_size, k=k, lambda=lambda, max.epochs=max.epochs, miniBatch_max_iters=miniBatch_max_iters, miniBatch_size=miniBatch_size, seed=seed)      
        object = quantile_norm(object, do.center = T)
        object = runUMAP(object, n_neighbors=30, min_dist=0.3, distance ="cosine")
        object = transfer_labels(object, annotations, knn)
        object_list[[modality]] = object
      } else {
        if(length(files[["meth"]])>0){
          hdf5_files = paste0(filepath, "/", region, "/Analysis", analysis_num , "_", region, "/",modality,"_subanalysis_",gsub(".RDS", ".H5",files[["meth"]]))
          data_names = gsub("(_qc.RDS)", "",files[["meth"]])
          
          genes_use = object_list[["rna"]]@var.genes
          for(i in 1:length(hdf5_files)){
            current_matrix = readRDS(paste0(filepath, "/",  region, "/Analysis", analysis_num , "_", region, "/",files[["meth"]][i]))
            current_matrix = current_matrix[rownames(current_matrix) %in% vargenes_rna,]
            genes_use = intersect(genes_use, rownames(current_matrix))
            current_matrix = as.matrix(max(current_matrix) - current_matrix)
            
            met.cell.data = data.frame(dataset = rep(data_names[i], ncol(current_matrix)), nUMI = rep(1,ncol(current_matrix)), nGene = rep(1,ncol(current_matrix)), barcode = colnames(current_matrix))
            met.cell.data$dataset = as.character(met.cell.data$dataset)
            met.cell.data$barcode = as.character(met.cell.data$barcode)
            rhdf5::h5createFile(hdf5_files[i])
            rhdf5::h5createGroup(hdf5_files[i], "matrix")
            rhdf5::h5write(colnames(current_matrix), file=hdf5_files[i], name="matrix/barcodes")
            rhdf5::h5createGroup(hdf5_files[i], file.path("matrix", "features"))
            rhdf5::h5write(rownames(current_matrix), file=hdf5_files[i], name="matrix/features/name")
            rhdf5::h5write(dim(current_matrix), file=hdf5_files[i], name="matrix/shape")
            rhdf5::h5createDataset(hdf5_files[i],"matrix/data",dims=3,storage.mode="double", chunk = 1)
            rhdf5::h5write(1:3, file=hdf5_files[i], name="matrix/data", index = list(1:3))
            rhdf5::h5createDataset(hdf5_files[i],"matrix/indices",dims=3,storage.mode="integer", chunk = 1)
            rhdf5::h5write(1:3, file=hdf5_files[i], name="matrix/indices",index = list(1:3)) # already zero-indexed.
            rhdf5::h5createDataset(hdf5_files[i],"matrix/indptr",dims=3,storage.mode="integer", chunk = 1)
            rhdf5::h5write(1:3, file=hdf5_files[i], name="matrix/indptr",index = list(1:3))
            rhdf5::h5createDataset(hdf5_files[i],"norm.data",dims=3,storage.mode="integer", chunk = 1)
            rhdf5::h5write(1:3, file=hdf5_files[i], name="norm.data",index = list(1:3))
            rhdf5::h5createDataset(hdf5_files[i],"scale.data",dims=c(nrow(current_matrix),ncol(current_matrix)),storage.mode="double", chunk = c(nrow(current_matrix), chunk = c(nrow(meth_1),min(chunk_size, ncol(meth_1)))))
            rhdf5::h5write(current_matrix, file=hdf5_files[i], name="scale.data",index = list(NULL, 1:ncol(current_matrix)))
            rhdf5::h5write(met.cell.data, file=hdf5_files[i], name="cell.data")
            rhdf5::h5closeAll()
          }
          
          names(hdf5_files) = data_names
          object = createLiger(as.list(hdf5_files))
          object@var.genes = genes_use
          object = online_iNMF(object,h5_chunk_size = chunk_size, k=k, lambda=lambda, max.epochs=max.epochs, miniBatch_max_iters=miniBatch_max_iters, miniBatch_size=miniBatch_size, seed=seed)
          object = quantile_norm(object, do.center = T)
          object = runUMAP(object, n_neighbors=30, min_dist=0.3, distance ="cosine")
          object = transfer_labels(object, annotations, knn)
          object_list[[modality]] = object
        }
      }
    }
    return(object_list)
  }
}



####################################### Function to create a master QC file
library(openxlsx)
master_csv = function(region, filepath = "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/"){
  if(str_sub(filepath, start = -1L, end = -1L) != "/"){
    filepath = paste0(filepath, "/")
  }
  output_filepath = paste0(filepath, region,"/Master_QC.xlsx")
  wb <- createWorkbook()
  for (analysis in 1:5){
    input_path = paste0(filepath, region, "/Analysis", analysis , "_", region, "/", region, "_Overall_qc.csv")
    content = read.table(input_path)
    content = data.frame(content)
    current_sheet = paste0("Analysis_", analysis)
    current_sheet
    addWorksheet(wb = wb, sheetName = current_sheet)
    writeData(wb, sheet = current_sheet, x = content)
  }
  saveWorkbook(wb, output_filepath)
}


library(sjmisc)
library(stringr)
library(rliger)
library(dplyr)
library(edgeR)
preprocess_and_run = function(filepath, region, analysis_num, chunk_size, num_genes = 2500, gene_num_tolerance = 100, var_thresh_start = 2, max_var_thresh = 4, customGeneList = NA, return.object = FALSE, qn_ref = NA, knownAnnotations = "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Base_Reference_Files/Reference_Annotations_updated_with_MacoskoLabels.RDS", MaxFactor = FALSE){
  qc_files = list.files(paste0(filepath, "/", region, "/Analysis", analysis_num , "_", region, "/"))
  qc_files = grep(paste0(region,"_(sc10Xv3_|smartseq_|atac_|meth_|sn10Xv3_|sc10Xv2_).*(qc.RDS)"), qc_files, value = TRUE)
  non_meth_files = grep("meth",qc_files, value = TRUE, invert = TRUE)
  rhdf5::h5closeAll()
  print("Writing H5s")
  hdf5_files = paste0(filepath, "/", region, "/Analysis", analysis_num , "_", region, "/",gsub(".RDS", ".H5",non_meth_files))
  for (i in 1:length(hdf5_files)){
    current_matrix = Matrix::Matrix(readRDS(paste0(filepath, "/", region, "/Analysis", analysis_num , "_", region, "/",non_meth_files[i])), sparse = TRUE)
    rhdf5::h5createFile(hdf5_files[[i]])
    rhdf5::h5createGroup(hdf5_files[[i]], "matrix")
    rhdf5::h5write(current_matrix@Dimnames[[2]], file=hdf5_files[[i]], name="matrix/barcodes") # cell barcodes
    rhdf5::h5createGroup(hdf5_files[[i]], file.path("matrix", "features"))
    rhdf5::h5write(current_matrix@Dimnames[[1]], file=hdf5_files[[i]], name="matrix/features/name")
    rhdf5::h5write(dim(current_matrix), file=hdf5_files[[i]], name="matrix/shape")
    rhdf5::h5createDataset(hdf5_files[[i]],"matrix/data",dims=length(current_matrix@x),storage.mode="double", chunk = min(chunk_size, ncol(current_matrix)))
    rhdf5::h5write(current_matrix@x, file=hdf5_files[[i]], name="matrix/data", index = list(1:length(current_matrix@x)))
    rhdf5::h5createDataset(hdf5_files[[i]],"matrix/indices",dims=length(current_matrix@i),storage.mode="integer", chunk = min(chunk_size, ncol(current_matrix)))
    rhdf5::h5write(current_matrix@i, file=hdf5_files[[i]], name="matrix/indices",index = list(1:length(current_matrix@i))) # already zero-indexed.
    rhdf5::h5createDataset(hdf5_files[[i]],"matrix/indptr",dims=length(current_matrix@p),storage.mode="integer", chunk = min(chunk_size, ncol(current_matrix)))
    rhdf5::h5write(current_matrix@p, file=hdf5_files[[i]], name="matrix/indptr",index = list(1:length(current_matrix@p)))
    rhdf5::h5closeAll()
  }
  data_names = gsub("(_qc.RDS)", "", non_meth_files)
  names(hdf5_files) = data_names
  object = createLiger(raw.data = as.list(hdf5_files))
  print("Normalizing data")
  object = normalize(object, chunk = chunk_size)
  datasets_use = grep("[sc10Xv3_|smartseq_|sc10Xv2_]",non_meth_files)
  object = selectGenes(object, var.thresh = var_thresh_start, datasets.use = datasets_use)
  if (is.na(customGeneList)){
    print("Selecting Genes")
    var_thresh_old = var_thresh_start
    high = max_var_thresh
    low = 0
    while(abs(length(object@var.genes)-num_genes) > gene_num_tolerance){
      if(length(object@var.genes)>num_genes){
        var_thresh_new = (var_thresh_old+high)/2
        low = var_thresh_old
      } else {
        var_thresh_new = (var_thresh_old+low)/2
        high = var_thresh_old
      }
      object = selectGenes(object, var.thresh = var_thresh_new, datasets.use = datasets_use)
      var_thresh_old = var_thresh_new
    }
    message(paste0(length(object@var.genes), " genes found with var.thresh = ",var_thresh_old))
  } 
  if (!is.na(customGeneList)){
    print("Using custom gene list")
    customGenes = readRDS(customGeneList)
    object@var.genes = customGenes
  }
  meth_files = setdiff(qc_files, non_meth_files)
  
  if(length(meth_files)>0){
    hdf5_files = paste0(filepath, "/", region, "/Analysis", analysis_num , "_", region, "/",gsub(".RDS", ".H5",meth_files))
    data_names = gsub("(_qc.RDS)", "",meth_files)
    current_matrix = readRDS(paste0(filepath, "/",  region, "/Analysis", analysis_num , "_", region, "/",meth_files[1]))
    meth_genes = rownames(current_matrix)
    new_genes = subset(object@var.genes, object@var.genes %in% meth_genes)
    object@var.genes = new_genes
    message(paste0(length(object@var.genes), " Genes after accounting for those in common with methylation data"))
  }
  
  print("Scaling Object")
  
  object = scaleNotCenter(object, chunk = chunk_size)
  
  if(length(meth_files)>0){
    hdf5_files = paste0(filepath, "/", region, "/Analysis", analysis_num , "_", region, "/",gsub(".RDS", ".H5",meth_files))
    data_names = gsub("(_qc.RDS)", "",meth_files)
    
    for(i in 1:length(meth_files)){
      current_matrix = readRDS(paste0(filepath, "/",  region, "/Analysis", analysis_num , "_", region, "/",meth_files[i]))
      #current_matrix = current_matrix[vargenes,]
      
      current_matrix = current_matrix[object@var.genes,]
      current_matrix = as.matrix(max(current_matrix) - current_matrix)
      met.cell.data = data.frame(dataset = rep(data_names[i], ncol(current_matrix)), nUMI = rep(1,ncol(current_matrix)), nGene = rep(1,ncol(current_matrix)), barcode = colnames(current_matrix))
      met.cell.data$dataset = as.character(met.cell.data$dataset)
      met.cell.data$barcode = as.character(met.cell.data$barcode)
      fname = hdf5_files[[i]]
      rhdf5::h5createFile(fname)
      rhdf5::h5createGroup(fname, "matrix")
      rhdf5::h5write(colnames(current_matrix), file=fname, name="matrix/barcodes")
      rhdf5::h5createGroup(fname, file.path("matrix", "features"))
      rhdf5::h5write(rownames(current_matrix), file=fname, name="matrix/features/name")
      rhdf5::h5write(dim(current_matrix), file=fname, name="matrix/shape")
      rhdf5::h5createDataset(fname,"matrix/data",dims=3,storage.mode="double", chunk = 1)
      rhdf5::h5write(1:3, file=fname, name="matrix/data", index = list(1:3))
      rhdf5::h5createDataset(fname,"matrix/indices",dims=3,storage.mode="integer", chunk = 1)
      rhdf5::h5write(1:3, file=fname, name="matrix/indices",index = list(1:3)) # already zero-indexed.
      rhdf5::h5createDataset(fname,"matrix/indptr",dims=3,storage.mode="integer", chunk = 1)
      rhdf5::h5write(1:3, file=fname, name="matrix/indptr",index = list(1:3))
      rhdf5::h5createDataset(fname,"norm.data",dims=3,storage.mode="integer", chunk = 1)
      rhdf5::h5write(1:3, file=fname, name="norm.data",index = list(1:3))
      rhdf5::h5createDataset(fname,"scale.data",dims=c(nrow(current_matrix),ncol(current_matrix)),storage.mode="double", chunk = c(nrow(current_matrix),min(chunk_size, ncol(current_matrix))))
      rhdf5::h5write(current_matrix, file=fname, name="scale.data",index = list(NULL, 1:ncol(current_matrix)))
      rhdf5::h5write(met.cell.data, file=fname, name="cell.data")
      rhdf5::h5closeAll()
    }
  }
  
  hdf5_files = paste0(filepath, "/", region, "/Analysis", analysis_num , "_", region, "/",gsub(".RDS", ".H5",qc_files))
  data_names = gsub("(_qc.RDS)", "",qc_files)
  names(hdf5_files) = data_names
  rhdf5::h5closeAll()
  
  var.genes = object@var.genes
  rm(object)
  rhdf5::h5closeAll()
  
  object_new = createLiger(as.list(hdf5_files))
  object_new@var.genes = var.genes
  object_new = online_iNMF(object_new, k = 30, lambda = 5, max.epochs = 20, seed = 123)
  if (!is.na(qn_ref)){
    object_new = quantile_norm(object_new, do.center = T, ref_dataset = qn_ref)
  } else{
    object_new = quantile_norm(object_new, do.center = T)
  }
  liger_name = paste0(filepath, "/", region, "/Analysis", analysis_num, "_", region, "/onlineINMF_",region, "_object.RDS" )
  print("Saving LIGER object")
  object_new = runUMAP(object_new,  n_neighbors=30, min_dist=0.3, distance ="cosine")
  saveRDS(object_new, liger_name)
  
  if(MaxFactor == TRUE){
    print("Performing Max Cell Type Assignment")
    liger_low = liger_high = max_factor_assignment(object_new)
  }
  
  if ( MaxFactor== FALSE){
  print("Running low resolution louvain clustering")
  if (analysis_num == 2){
    low_resolution = 0.75
    high_resolution = 1
  }
  if (analysis_num != 2){
    low_resolution = 0.25
    high_resolution = 1
  }
  liger_low = louvainCluster(object_new, resolution = low_resolution, k = 200)
  liger_high = louvainCluster(object_new, resolution = high_resolution, k = 200)
  }
  # #Plot both high and low resolution UMAPs####################
  print("Plotting unlabeled UMAPS....")
  plots_low = plotByDatasetAndCluster(liger_low, return.plots = TRUE, text.size = 6)
  low_umap1 =paste0(filepath, "/",  region, "/Analysis", analysis_num, "_", region, "/Images/Umap1_", region, "_Analysis_", analysis_num, ".pdf")
  low_umap2 =paste0(filepath, "/", region, "/Analysis", analysis_num, "_", region, "/Images/Umap2_", region, "_Analysis_", analysis_num, "unlabeled.lowresolution.pdf")
  low_umap1png =paste0(filepath, "/",  region, "/Analysis", analysis_num, "_", region, "/Images/Umap1_", region, "_Analysis_", analysis_num, ".png")
  low_umap2png =paste0(filepath, "/", region, "/Analysis", analysis_num, "_", region, "/Images/Umap2_", region, "_Analysis_", analysis_num, "unlabeled.lowresolution.png")
  png(low_umap1png, 1000, 800)
  print(plots_low[[1]])
  dev.off()
  png(low_umap2png, 1000, 800)
  print(plots_low[[2]])
  dev.off()
  pdf(low_umap1, width = 10, height = 8)
  print(plots_low[[1]])
  dev.off()
  pdf(low_umap2, width = 10, height = 8)
  print(plots_low[[2]])
  dev.off()
  plots_high = plotByDatasetAndCluster(liger_high, return.plots = TRUE, text.size = 6)
  high_umap2 =paste0(filepath,"/",  region, "/Analysis", analysis_num, "_", region, "/Images/Umap2_", region, "_Analysis_", analysis_num, "unlabeled.highresolution.pdf")
  pdf(high_umap2, width = 10, height = 8)
  print(plots_high[[2]])
  dev.off()
  high_umap2png =paste0(filepath,"/",  region, "/Analysis", analysis_num, "_", region, "/Images/Umap2_", region, "_Analysis_", analysis_num, "unlabeled.highresolution.png")
  png(high_umap2png, 1000, 800)
  print(plots_high[[2]])
  dev.off()
  
  ##########################################################
  # #Generate a saved results table
  result = data.frame(object_new@tsne.coords)
  result$dataset = object_new@cell.data$dataset
  result$lowRcluster = liger_low@clusters
  result$highRcluster = liger_high@clusters
  #
  # #Here we add the old annotations, and also generate a UMAP labeled with these old annotations. Need to discuss with J.S. how to best incorporate nearest neighbors clustering/annotation function
  master = readRDS(knownAnnotations)
  if (analysis_num == 1){
    #Graph neuronal vs. non-neuronal
    annies = master$Level1
    names(annies) = master$Cell_Barcodes
  }
  if (analysis_num == 3){
    #Graph Excitatory vs Inhibitory
    annies = master$Class
    names(annies) = master$Cell_Barcodes
  }
  if (analysis_num == 2 | analysis_num == 4 | analysis_num == 5){
    #   #Graph Type
    annies = master$Type
    names(annies) = master$Cell_Barcodes
    
  }

  clusts = liger_low@clusters
  if (MaxFactor == TRUE){
    clusts2 = data.frame(clusts)
    clusts2$Barcode = rownames(clusts2)
    annies2 = data.frame(annies)
    annies2$Barcode = rownames(annies2)
    full = left_join(clusts2, annies2)
    full2 = full$annies
    names(full2) = full$Barcode
    liger_low@clusters = full2
  }
  if(MaxFactor == FALSE){
  result$ann = annies[names(clusts)]
  liger_low@clusters = as.factor(annies[names(clusts)])
  }
  
  # #Return UMAP with OG labels
  print("Plotting labeled UMAPS....")
  plots_low = plotByDatasetAndCluster(liger_low, return.plots = TRUE, text.size = 6)
  low_umap2 =paste0(filepath,"/", region, "/Analysis", analysis_num, "_", region, "/Images/Umap2_", region, "_Analysis_", analysis_num, "labeled.pdf")
  pdf(low_umap2, width = 10, height = 8)
  print(plots_low[[2]])
  dev.off()
  
  low_umap2.png =paste0(filepath,"/", region, "/Analysis", analysis_num, "_", region, "/Images/Umap2_", region, "_Analysis_", analysis_num, "labeled.png")
  png(low_umap2.png, 1000, 800)
  print(plots_low[[2]])
  dev.off()
  # #Rename Results Table
  names(result) = c("UMAP1","UMAP2","dataset","lowRcluster","highRcluster", "OG_Annotations")
  result$Barcode = rownames(result)
  print("Outputting Annotation CSV")
  cluster_breakdowns_high = result %>% group_by(highRcluster, dataset, OG_Annotations)  %>% tally()
  cluster_breakdowns_low = result %>% group_by(lowRcluster, dataset, OG_Annotations)  %>% tally()
  #
  output_filepath = paste0(filepath,"/",  region, "/Analysis", analysis_num, "_", region, "/Cluster_Breakdowns_",region, "_Analysis_", analysis_num, ".xlsx")
  wb <- createWorkbook()
  addWorksheet(wb = wb, sheetName = "LowRes")
  writeData(wb, sheet = "LowRes", x = cluster_breakdowns_low)
  addWorksheet(wb = wb, sheetName = "HighRes")
  writeData(wb, sheet = "HighRes", x = cluster_breakdowns_high)
  saveWorkbook(wb, output_filepath)
  if (analysis_num == 1 | analysis_num == 3){
    print("Performing max final annotations")
    num_clusts = length(unique(result$highRcluster))
    max_assignments = na.omit(result)
    max_assignments = max_assignments %>% group_by(highRcluster, OG_Annotations)  %>% tally() %>% filter(n == max(n))
    max_assignments =max_assignments[,c("highRcluster", "OG_Annotations")]
    colnames(max_assignments) = c("highRcluster", "highRAnnotations")
    #this is to accomodate for ties in max factor assignment
    if(length(unique(max_assignments$highRcluster)) != length(max_assignments$highRcluster)){
      warning("Clusters with tied OG annotations. Please manually annotate these clusters before proceeding with further analyses")
      duplicated_clusters = max_assignments$highRcluster[duplicated(max_assignments$highRcluster)]
      for (dups in duplicated_clusters){
        max_assignments$highRAnnotations[max_assignments$highRcluster == dups] <- "NA"
      }
    }
    max_assignments = unique(max_assignments)
    max_assign = length(unique(max_assignments$highRcluster))
    result = left_join(result, max_assignments)
    if(max_assign != num_clusts){
      warning("Completely unannotated clusters. Please annotate these clusters before proceeding with further analyses")
    }
  }
  results_filename = paste0(filepath,"/",  region, "/Analysis", analysis_num, "_", region, "/Analysis", analysis_num, "_", region,"_Results_Table.RDS")
  saveRDS(result, results_filename)
  #Also, save cell data information
  cell_data = liger_low@cell.data
  cell_data_filename = paste0(filepath,"/",  region, "/Analysis", analysis_num, "_", region, "/Analysis", analysis_num, "_", region,"_CellData.RDS")
  saveRDS(cell_data, cell_data_filename)
  
}





#' Doublet detection in single-cell RNA sequencing data
#' This function generaetes artificial nearest neighbors from existing single-cell RNA
#' sequencing data. First, real and artificial data are merged. Second, dimension reduction
#' is performed on the merged real-artificial dataset using PCA. Third, the proportion of
#' artificial nearest neighbors is defined for each real cell. Finally, real cells are rank-
#' ordered and predicted doublets are defined via thresholding based on the expected number
#' of doublets.
#'
#' zizheny@alleninstitute.org
#'
#' @param data count matrix
#' @param select.genes : selected genes Default rownames(data) 
#' @param proportion.artificial : The proportion (from 0-1) of the merged real-artificial dataset
#' that is artificial. In other words, this argument defines the total number of artificial doublets.
#' Default is set to 25%, based on optimization on PBMCs (see McGinnis, Murrow and Gartner 2018, BioRxiv).
#' @param k : max number of nearest neighbor RANN:nn2() 
#' @param plot : flag for plots
#' @return doublet.score : doublet score
#' @examples
#' doubletFinder <- function(data, select.genes, proportion.artificial = 0.20,
#'                           k = round(pmin(100, ncol(data) * 0.01)), plot=FALSE)
#'

#library(devtools)
#devtools::install_github("AllenInstitute/scrattch.hicat")
#Note: for successful installation, I first had to install:'impute' and 'GO.db' from Bioconductor
library(scrattch.hicat)
library(RANN)
library(ggplot2)

doubletFinder <- function(data, select.genes, proportion.artificial = 0.20,
                          k = round(pmin(100, ncol(data) * 0.01)), plot=FALSE) {
  
  ## Step 1: Generate artificial doublets from Seurat object input
  print("Creating artificial doublets...")
  real.cells <- colnames(data)
  
  n_real.cells <- length(real.cells)
  n_doublets <- round(n_real.cells/(1-proportion.artificial)-n_real.cells)
  real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  doublets <- (data[ , real.cells1] + data[ , real.cells2])
  colnames(doublets) <- paste("X", 1:n_doublets, sep="")
  
  data_wdoublets <- cbind(data, doublets)  
  norm.dat = logCPM(data_wdoublets)
  
  print("Running PCA")
  if(ncol(data) > 10000){
    sampled.cells = sample(1:ncol(data), pmin(ncol(data),10000))
    rd.dat = rd_PCA(norm.dat, select.genes=select.genes, select.cells= colnames(data_wdoublets), th=0, max.pca=50, sampled.cells= sampled.cells)
  }
  else{
    #sampled.cells = sample(1:ncol(data), ncol(data)*0.6)
    rd.dat = rd_PCA(norm.dat, select.genes=select.genes, select.cells= colnames(data_wdoublets), th=0, max.pca=50)
  }
  
  rd.dat<-rd.dat$rd.dat ##depending on hicat version used rd.dat is embedded or not.
  print("Initialize pANN structure") 
  knn.result = RANN::nn2(rd.dat, k=k)
  knn.idx = knn.result[[1]]
  knn.dist = knn.result[[2]]
  
  num = ncol(data)
  knn.result1 = RANN::nn2(rd.dat[1:num,], rd.dat[(num+1):nrow(rd.dat),], k = 10)
  knn.dist1 = knn.result1[[2]]
  
  dist.th = mean(as.vector(knn.dist1)) + 1.64 * sd(as.vector(knn.dist1))
  
  doublet.freq = knn.idx  > ncol(data) & knn.dist < dist.th
  doublet.score =  doublet.freq[1:ncol(data),]
  row.names(doublet.score) = colnames(data)
  doublet.score = pmax(rowMeans(doublet.freq),rowMeans(doublet.freq[,1:ceiling(k/2)]))
  names(doublet.score) = colnames(data_wdoublets)
  
  if (plot == TRUE) {
    print("plotting")
    ds = pmax(rowMeans(doublet.freq),rowMeans(doublet.freq[,1:ceiling(k/2)])) 
    ds=as.data.frame(ds)
    ds$sample <- colnames(data_wdoublets)
    
    ds$group <- ""
    idx <- startsWith(ds$sample,"X")
    ds[idx, "group"] <- "artifical doublets"
    idx <- !startsWith(ds$sample,"X")
    ds[idx, "group"] <- "samples"
    
    plot.title <- gsub("^.*?-","",ds[1,2])
    
    p=ggplot2::ggplot(ds, aes(x = doublet.score, fill=group, color=group)) +geom_density(alpha=0.4)+scale_color_manual(values=c("#F9627D","#2F3061"))+scale_fill_manual(values=c("#F9627D","#2F3061")) +labs(title=plot.title)
    
    return(list(doublet.score, p))
    
  } else {  
    return(doublet.score) 
  }
  
}

#' QC score based on marker genes
#' 
#' sum of marker genes is used as the qc.score
#'
#' @param data : count matrix
#' @param qc.gene.FN : filename holding QC gene list 
#' @return qc.score : qc score
#' @examples
#' qc.score <- cal_qc_score(data, qc.gene.FN="AIBS_qc_genes_10X.csv") 

cal_qc_score <- function (data, qc.gene.FN="AIBS_qc_genes_10X.csv") {
  
  qc.genes = read.csv(qc.gene.FN)[, "Gene"]
  common.qc.genes = intersect(qc.genes, row.names(data))
  norm.dat = logCPM(data[common.qc.genes,])
  qc.score = apply(norm.dat, 2, sum)
  #qc.score = Matrix::colSums(norm.dat[common.qc.genes,])
  
  return(qc.score)
}

max_factor_assignment = function(object){
  h_factors = object@H.norm
  max_assignments = apply(h_factors, 1, which.max)
  object@clusters = as.factor(max_assignments)
  return(object)
}

# #################### Notes on how we generated the annotation labels
# ann = readRDS("/nfs/turbo/umms-welchjd/BRAIN_initiative/cell_type_annotations_full.RDS")  #4 lists
# annotations = as.factor(c(as.character(ann[[1]]), as.character(ann[[2]]), as.character(ann[[3]]), as.character(ann[[4]])));  # make each one a factorized character and combine into a giant list
# names(annotations) = c(names(ann[[1]]), names(ann[[2]]), names(ann[[3]]), names(ann[[4]]))  #name each element in the giant list
# annot = data.frame(annotations)
# colnames(annot) = c("Type")
# annot$Cell_Barcodes = rownames(annot)
# types_annot = read.csv("/nfs/turbo/umms-welchjd/akriebel/Brain_initiative/cell_type_annot.csv")[1:121,]
#  types_annot$Class = sub("Non-neuron", "NonN", types_annot$Class)
#  types_annot$Class = sub("Excitatory Neuron", "Exc", types_annot$Class)
#  types_annot$Class = sub("Inhibitory Neuron", "Inh", types_annot$Class)
#  master = left_join(annot, types_annot)
#  decoder = data.frame(c("NonN", "Exc", "Inh"), c("NonN", "Neu", "Neu"))
#  colnames(decoder) = c("Class", "Level1")
#  master = left_join(master, decoder)
# saveRDS(master, "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Base_Reference_Files/Reference_Annotations.RDS")
# huh = readRDS("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Base_Reference_Files/Reference_Annotations.RDS")
# unique_cell_types = data.frame(unique(huh$Type))
# write.csv(unique_cell_types, "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Base_Reference_Files/unique_celltypes.csv")
# knownAnnotations = readRDS("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Base_Reference_Files/Reference_Annotations.RDS")
# macosko = read.csv("/scratch/welchjd_root/welchjd0/akriebel/Macosko_QC/fc-secure-1d9abad7-eef0-4161-92fa-4696b6131682/not-uniform/macosko-lab/20220531_draft_clustering_cells.tsv", sep = "\t", header = FALSE)
# library(tidyr)
# cell_annotations = macosko %>% separate("V3", into = c("C1", "C2", "C3", "C4", "C5","C6", "C7", "C8"), remove = FALSE) #3,578,959
# cell_annotations = cell_annotations[,c("V1", "V3", "C1")]
# colnames(cell_annotations) = c("Cell_Barcodes", "Type", "Class")
# cell_annotations$Level1 = cell_annotations$Class
# known = rbind(knownAnnotations, cell_annotations)
# saveRDS(knownAnnotations, "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Base_Reference_Files/Reference_Annotations_updated_with_MacoskoLabels.RDS")

############################ Functions for Shiny App
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
    results_name = paste0(filepath, region, "/Analysis", analysis_num, "_", region, "/Analysis", analysis_num, "_", region, "_Results_Table.RDS")
    results = readRDS(results_name)
    matrix = readRDS(matrix_filename)
    if (datasets[[i]] == "meth_1" | datasets[[i]] == "meth_2" ){
      genes_needed = subset(rownames(matrix), rownames(matrix) %in% genes_oi)
      current_matrix = matrix[genes_needed,]
      current_matrix = as.matrix(max(current_matrix) - current_matrix)
      cells_io = subset(colnames(current_matrix), colnames(current_matrix) %in% results$Barcode)
      dataset_matrices[[i]] = current_matrix[,cells_io]
    } 
    if (datasets[[i]] != "meth_1" & datasets[[i]] != "meth_2" ){
      ligs = createLiger(list(newdata = matrix))
      ligs = normalize(ligs)
      ligs = selectGenes(ligs, var.thresh = 0.0)
      genes_needed = subset(ligs@var.genes, ligs@var.genes %in% genes_oi)
      ligs@var.genes = genes_needed
      ligs = scaleNotCenter(ligs)
      current_matrix = t(ligs@scale.data$newdata)
      cells_io = subset(colnames(current_matrix), colnames(current_matrix) %in% results$Barcode)
      
      
      dataset_matrices[[i]] = current_matrix[,cells_io]
    }
  }
  names(dataset_matrices) = dataset_names
  output_file = paste0(filepath, region, "/Analysis", analysis_num, "_", region, "/Images/GeneExpression_", region,"_Analysis", analysis_num, ".RDS")
  saveRDS(dataset_matrices, output_file)}
#

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

deconvolve_spatial = function(filepath,
  region,
  spatial.data.file,
  annotation.level = 3,
  n.cells = 500,
  known.annotations = "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Base_Reference_Files/Reference_Annotations.RDS",
  deconv.gene.num = 2000,
  gene.num.tol = 50,
  clusters = NULL,
  lambda = 1,
  thresh = 1e-8,
  max.iters = 100,
  nrep = 1,
  rand.seed = 1,
  print.obj = FALSE,
  verbose = TRUE,
  z = 1){
    message("Loading Data")  
    object_path = paste0(filepath, region, "/Analysis1_", region, "/onlineINMF_",region, "_object.RDS" )
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
        analysis_clusters = paste0(analysis_num, "_", analysis_results$lowRannotations)
        names(analysis_clusters) = analysis_results$Barcode
        clusters = c(clusters, analysis_clusters)
      }
      clusters = clusters[clusters != ""]
    
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

    spatial.data = readRDS(spatial.data.file)

    spatial.data[spatial.data == -1] = NA
    genes_NA = apply(spatial.data, MARGIN = 1, function(x){sum(is.na(x))})
    mean_genes_NA = mean(genes_NA)
    genes_use = genes_NA < (mean_genes_NA + z * mean_genes_NA)
    spatial.data = spatial.data[genes_use,]
    #simplest way to handle this
    spatial.data[is.na(spatial.data)] = 0
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
    
    saveRDS(list(chisq_vals = chisq_list, genes_used = gene_vec), paste0(filepath,"/",  region, "/", region,"_Deconvolution_Output/gene_selection_output.RDS"))

    message("Learning gene signatures")
    E = lapply(norm.data, function(x){as.matrix(t(x)[,gene_vec])})
    
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
    
    message("Deconvolving spatial data")
    spatial.data = t(scale(t(spatial.data[gene_vec,]), center = FALSE))
    spatial.data[spatial.data < 0 ] = 0
    deconv_h = t(rliger:::solveNNLS(t(W),spatial.data))
    colnames(deconv_h) = clust_levels
    deconv_frac = t(apply(deconv_h, MARGIN = 1, function(x){x/sum(x)}))
    rownames(deconv_frac) = rownames(deconv_h) = colnames(spatial.data)
    saveRDS(list(raw = deconv_h, proportions = deconv_frac), paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/deconvolution_output.RDS"))
  }

assign_single_cells = function(
  filepath,
  region
){
  proportions = readRDS(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/deconvolution_output.RDS"))[[2]]
  cell_types = colnames(proportions)
  max = as.factor(apply(proportions, MARGIN = 1, function(x){cell_types[which.max(x)]}))
  names(max) = rownames(proportions)
  saveRDS(max, paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/deconvolution_max_factor_assignment.RDS"))
}

generate_loading_gifs = function(
  filepath,
  region,
  coords,
  mat.use = "proportions",
  cell.types.plot = NULL,
  dims = c(500, 500)
){
  library(rgl)
  loadings = readRDS(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/deconvolution_output.RDS"))
  if(!dir.exists(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/gifs"))){
    dir.create(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/gifs"))
  }
  grDevices::palette(viridis::viridis(option="A",n=50,direction = -1))
  if(is.null(cell.types.plot)){
    cell.types.plot = colnames(loadings[[mat.use]])
  } else {
    cell.types.plot = intersect(cell.types.plot, colnames(loadings[[mat.use]]))
  }
  if(is.character(coords)){
    coords = readRDS(coords)
  }
  for(cell_type in cell.types.plot){
    colors_view = as.numeric(cut(loadings[[mat.use]][,cell_type],breaks=50))
    try(rgl.close(), silent = TRUE)
    open3d(windowRect = c(0,0, dims[1], dims[2]));
    plot3d(coords[,1],coords[,2],coords[,3],col = colors_view,aspect=c(67,41,58),xlab="Anterior-Posterior",ylab="Inferior-Superior",zlab="Left-Right",size=5, add = TRUE)
    decorate3d(xlab = colnames(coords)[1], ylab = colnames(coords)[2],zlab = colnames(coords)[3], box = FALSE, axes = FALSE)
    axes3d(c("x--","y--","z--"))#axes3d(c("x--","y--","z--"))
    movie3d(spin3d(axis = c(0, 0, 1)), duration = 20, movie = paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/gifs/", region, "_",sub("/", ".",sub(" ", "_",cell_type)),"_spatial_summary"))
  }
}

#' Convert proportions generated during the deconvolution into cell-type
#' assignments, for use with single cell spatial modalities (i.e. slideseq)
#'
#' @param filepath Path to directory within which the atlas structure was generated
#' @param region A string corresponding to the name of an anatomical region
#' @param layer.list A named list of spatial sample names, corresponding to
#'    groupings to be summarized by the function
#' @param spatial.data.file The location of a file containing spatial expression data
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
  spatial.data.file,
  plot = FALSE,
  type = "mean",
  mat.use = "proportions",#"assignment
  use.cell.types = TRUE,
  cell.types.use = NULL,
  genes.use = NULL){
  if(mat.use == "assignment"){
    assignments = readRDS("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Subset/MOp/MOp_Deconvolution_Output/deconvolution_max_factor_assignment.RDS")
    sub_vec = rep(0,nlevels(assignments))
    loadings = Reduce(rbind, lapply(assignments, function(x){subbed_vec = sub_vec; subbed_vec[as.numeric(x)] = 1; return(subbed_vec)}))
    colnames(loadings) = levels(assignments)
  } else {
    loadings = readRDS(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/deconvolution_output.RDS"))[[mat.use]]
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
    saveRDS(cell.type.matrix, paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/cell_type_layer_summary.RDS"))
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
    saveRDS(cell.type.matrix, paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/gene_layer_summary.RDS"))
  }
  if(plot){
    if(!dir.exists(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/plots"))){
      dir.create(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/plots"))
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

      ggplot2::ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/plots/cell_type_layer_distribution.PNG"),
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
      ggplot2::ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/plots/",sub("/", ".",sub(" ", "_",i)),"_layer_distribution.PNG"),
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

        ggplot2::ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/plots/gene_layer_distribution.PNG"),
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
        ggplot2::ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/plots/",sub("/", ".",sub(" ", "_",i)),"_gene_layer_distribution.PNG"),
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
  plot = FALSE,
  mat.use = "proportions"){
  
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
  
  deconv_out = readRDS(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/deconvolution_output.RDS"))
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
    if(!dir.exists(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/plots"))){
      dir.create(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/plots"))
    }
    ggplot2::ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/plots/cell_type_signature_heatmap.PNG"),
                    heatmap_plot,
                    width = 1000,
                    height = 800,
                    units = "px")
    ggplot2::ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/plots/cell_type_signature_dendrogram.PNG"),
                    dendro_plot,
                    width = 500,
                    height = 400,
                    units = "px")
    ggplot2::ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/plots/cell_type_distribution_heatmap.PNG"),
                    heatmap_dist_plot,
                    width = 1000,
                    height = 800,
                    units = "px")
    ggplot2::ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/plots/cell_type_distribution_dendrogram.PNG"),
                    dendro_dist_plot,
                    width = 500,
                    height = 400,
                    units = "px")
  }
  saveRDS(list(cos_sim_signature = cos_sim, cor_sim_distribution = corr_sim_dist, dendro_sig = hierarchical_clust, dendro_dist = hierarchical_clust_dist),
          paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/gene_signature_analysis_summary.RDS"))
}

plot_layer = function(
  filepath,
  region,
  coords,
  spatial.data.file,
  axis,
  idx,
  mat.use = "proportions",
  use.cell.types = TRUE,
  cell.types.use = NULL,
  genes.use = NULL,
  display.plots = FALSE
){
  
  if(!dir.exists(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/plots"))){
      dir.create(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/plots"))
  }
  theme_set(theme_cowplot())

  if(is.character(coords)){
    coords = readRDS(coords)
  }
  coords = coords[coords[,axis] == idx, setdiff(colnames(coords),axis)]
  if(use.cell.types){
    loadings = readRDS(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/deconvolution_output.RDS"))
    if(!is.null(cell.types.use)){
      cell.types.use = intersect(cell.types.use, colnames(loadings[[mat.use]]))
    } else {
      cell.types.use = colnames(loadings[[mat.use]])
    }
    cell.types.use = cell.types.use[cell.types.use != ""]
    loadings = loadings[[mat.use]][rownames(loadings[[mat.use]]) %in% rownames(coords), cell.types.use]
    
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
          theme(legend.title = ggplot2::element_blank(),
                text = ggplot2::element_text(size = 8), 
                axis.text = ggplot2::element_text(size = 5))
      if(display.plots){
        print(slice_plot)
      }
      ggplot2::ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/plots/",
                    colnames(loadings)[i],"_",axis,"_",idx,".PNG"),
                    slice_plot,
                    width = 1000,
                    height = 800,
                    units = "px")
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
          theme(legend.title = ggplot2::element_blank(),
                text = ggplot2::element_text(size = 8), 
                axis.text = ggplot2::element_text(size = 5))
      if(display.plots){
        print(slice_plot)
      }
      ggplot2::ggsave(paste0(filepath,"/",  region,"/", region,"_Deconvolution_Output/plots/",
                    genes.use[i],"_",axis,"_",idx,".PNG"),
                    slice_plot,
                    width = 1000,
                    height = 800,
                    units = "px")
    }
  }
}




#################################################################################################### Understanding Marker Genes #####################################################################
################################# 
################### Need to generate a reference list object for each datatype
# genes = readRDS("C:/Users/april/OneDrive/Documents/Brain_Initative_v2_2022/Marker_genes.RDS")
# #Change Dataset filtering for atac, rna, and methylation
# genes_of_interest = filter(genes, genes$Dataset == "rna")
# genes_of_interest2 = genes_of_interest[,c("SubType", "M1", "M2", "M3", "M4", "M5", "M6")]
# genes_of_interest2 = genes_of_interest2 %>%  pivot_longer(!SubType, names_to = "Delete", values_to = "Gene")
# genes_of_interest2 = genes_of_interest2[, c("SubType", "Gene")]
# genes_of_interest2 = filter(genes_of_interest2, !is.na(genes_of_interest2$Gene))
# genes_of_interest2$Gene = sub(" ", "", genes_of_interest2$Gene)
# genes_of_interest2$Exists = 1
# genes_of_interest2 = unique(genes_of_interest2)
# gee = pivot_wider(genes_of_interest2, names_from = Gene, values_from  = Exists, values_fill = 0)
# 
# gee = data.frame(gee)
# cat_cells = data.frame(genes_of_interest$MajorType, genes_of_interest$SubType)
# cat_cells = unique(cat_cells)
# colnames(cat_cells) = c("MajorType", "SubType")
# 
# rna = left_join(cat_cells, gee)
# 
# cumulative = list(rna, atac, meth)
# names(cumulative) = c("rna", "atac", "meth")
# saveRDS(cumulative, "C:/Users/april/OneDrive/Documents/Brain_Initative_v2_2022/Cumulative_marker_genes.RDS")
####################################################################################################################################

