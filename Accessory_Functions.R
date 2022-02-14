
library(stringr)

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
}
################## QC function
#Input: List of filenames and types
#Types should be listed as "tenx", "smart", "atac", "meth", "huang"
#If more than one of each type of file exists, they should be listed as "tenx_1", "tenx_2", etc.
#output: All files will be out put as .RDS files. QC metrics will be stored in Analysis#_region_Overall_qc.RDS
library(sjmisc)
library(stringr)
library(rliger)
library(dplyr)
# Base QC table is available at "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Base_Reference_Files/QC_table.RDS"
# @slot analysis_num Determines directory and output file naming conventions, also ensures that Analysis one will ignore methylation Data
apply_qc = function(filenames, region, analysis_num , qc_table_path, filepath = "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/", verbose = TRUE){
  if(str_sub(filepath, start = -1L, end = -1L) != "/"){
    filepath = paste0(filepath, "/")
  }
  '%notin%' = Negate('%in%')
  qc_summary = setNames(data.frame(matrix(ncol = 9, nrow = 0)), c("Datatype", "Edition", "OriginalDimensions", "AnalysisDimensions", "FailnUMI", "nUMIThreshold", "FailMito", "MitoThreshold", "FinalDimensions"))
  
  for (file.listed in 1:length(filenames)){
    #i.e. for each file 
    # For Analysis != 1, we will need to filter the matrices on the basis on their results
    #i.e. Analysis one should be all cells, not including methylation
    #Analysis 2 should just be non-neuronal cells
    #Analysis 3 should be just neuronal cells + methy
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
        results_filename = paste0(filepath,region, "_BICCN/Analysis1_", region, "/Analysis1_", region, "_Results_Table.RDS")
        results = readRDS(results_filename)
        results = filter(results, results$ann =="NN")
        results = subset(results, rownames(results) %in% colnames(working_file))
        use.cells = rownames(results)
      } else {use.cells = colnames(working_file)}
    }
    if(analysis_num == 3){
      if(data.type != "meth"){
        results_filename = paste0(filepath,region, "_BICCN/Analysis1_", region, "/Analysis1_", region, "_Results_Table.RDS")
        results = readRDS(results_filename)
        results = filter(results, results$ann !="NN")
        results = subset(results, rownames(results) %in% colnames(working_file))
        use.cells = rownames(results)
      } else {use.cells = colnames(working_file)}
    }
    if(analysis_num == 4){
      results_filename = paste0(filepath,region, "_BICCN/Analysis3_", region, "/Analysis3_", region, "_Results_Table.RDS")
      results = readRDS(results_filename)
      results = filter(results, results$ann =="Inh")
      results = subset(results, rownames(results) %in% colnames(working_file))
      use.cells = rownames(results)
    }
    if(analysis_num == 5){
      results_filename = paste0(filepath,region, "_BICCN/Analysis3_", region, "/Analysis3_", region, "_Results_Table.RDS")
      results = readRDS(results_filename)
      results = filter(results, results$ann =="Exc")
      results = subset(results, rownames(results) %in% colnames(working_file))
      use.cells = rownames(results)
    }
    
    before_subset = dim(working_file)[[2]]
    working_file = working_file[,use.cells]
    after_subset = dim(working_file)[[2]]
    
    
    ######################### Pre-process   
    if(data.type != "meth"){
      #If datatype is ATAC or Tenx, filter for QC cells
      qc_table = readRDS(qc_table_path)
      if (data.type == "atac" | data.type == "tenx" | data.type == "huang"){
        if (data.type == "atac"){ qc_stats = filter(qc_table, qc_table$DataType == "atac")}
        if(data.type == "tenx" | data.type == "huang"){
          qc_stats = filter(qc_table, qc_table$DataType == "tenx")
          qc_stats = filter(qc_stats, qc_stats$Analysis == analysis_num)
          if (region != "OLF" & region != "CB"){
            qc_stats = filter(qc_stats, qc_stats$Region == "ALL")
          } else {
            qc_stats = filter(qc_stats, qc_stats$Region == region)
          }
        }
        nUMI_cutoff = qc_stats$nUMI
        mito_cutoff = qc_stats$mito
        ligs = createLiger(list(qc_mat = working_file))
        celldata = ligs@cell.data
        celldata$Mito = getProportionMito(ligs)
        celldata = filter(celldata, celldata$nUMI >= nUMI_cutoff)
        lost_nUMI = dim(celldata)[[1]]
        celldata = filter(celldata, celldata$Mito <= mito_cutoff)
        lost_mito = dim(celldata)[[1]]
        qc.matrix = working_file[, rownames(celldata)]
      }
      if (data.type == "smart"){
        qc.matrix = working_file
        lost_nUMI = lost_mito = before_subset - after_subset
        nUMI_cutoff = mito_cutoff = NA
        
      }
      if (data.type == "tenx" | data.type == "smart"){
        rik_converter = readRDS("/nfs/turbo/umms-welchjd/BRAIN_initiative/Ensembl_conversion/ConvertRik.Tenx.SMART.RDS")
        #Need to raise a warning if there are genes in the SMART matrix that are not in the conversion table
        refined_genes = rownames(qc.matrix)
        potential_outliers = subset(refined_genes, refined_genes %notin% rik_converter$Og_Barcode)
        if (length(potential_outliers >= 1)) { 
          warning("There are genes in the matrix that are not covered in the conversion table. Please update the conversion table to address these genes and then rerun the function.", immediate. = T)
        } 
        rowname_new = rik_converter$New_Barcode[match(rownames(qc.matrix), rik_converter$Og_Barcode)]
        rownames(qc.matrix) = rowname_new }
      #Count removed Cells
      lost_nUMI = after_subset - lost_nUMI
      lost_mito = lost_nUMI - lost_mito
      new.dim = dim(qc.matrix)[[2]]
      if (edition == 0){
        save.filename = paste0(filepath, region, "/Analysis", analysis_num , "_", region, "/", region, "_", data.type, "_qc.RDS")
      } else {save.filename = paste0(filepath, region, "/Analysis", analysis_num , "_", region, "/", region, "_" ,data.type,"_", edition, "_qc.RDS") }
      saveRDS(qc.matrix, save.filename)
      newest_qc = c(data.type, edition, before_subset, after_subset, lost_nUMI, nUMI_cutoff, lost_mito, mito_cutoff, new.dim)
      qc_summary = rbind(qc_summary, newest_qc)
    } 
    if( data.type == "meth") {
      #Remember we do not do methylation data for analysis 1 or 2
      if (analysis_num != 1 & analysis_num != 2){
        qc.matrix = as.matrix(max(working_file) - working_file)
        lost_nUMI = lost_mito = 0
        lost_nUMI = after_subset - lost_nUMI
        lost_mito = lost_nUMI - lost_mito
        new.dim = dim(qc.matrix)[[2]]
        
        nUMI_cutoff = mito_cutoff = NA
        if (edition == 0){
          save.filename = paste0(filepath, region, "/Analysis", analysis_num , "_", region, "/", region, "_",data.type, "_qc.RDS")
        } else { save.filename = paste0(filepath, region, "/Analysis", analysis_num , "_", region, "/", region, "_" ,data.type,"_", edition, "_qc.RDS") }
        saveRDS(qc.matrix, save.filename)
        newest_qc = c(data.type, edition, before_subset, after_subset, lost_nUMI, nUMI_cutoff, lost_mito, mito_cutoff, new.dim)
        qc_summary = rbind(qc_summary, newest_qc) }
      
    }
    if(verbose){
      message("Done Processing ", data.type, " ", edition)
    }
  }
  if (data.type != "meth" & data.type != "atac" & data.type != "smart" & data.type != "tenx"){
    warning("Uknown Data Type submitted", immediate. = T)
    
  }
  
  #Output a .csv to the Analysis directory
  csv_filename = paste0(filepath, region, "/Analysis", analysis_num , "_", region, "/", region, "_Overall_qc.csv")
  colnames(qc_summary) = c("Datatype", "Edition", "OriginalDimensions", "AnalysisDimensions", "FailnUMI", "nUMIThreshold", "FailMito", "MitoThreshold", "FinalDimensions")
  
  write.table(qc_summary, csv_filename)
}

preprocess_data = function(filepath, region, analysis_num, chunk_size, num_genes = 2500, gene_num_tolerance = 100, var_thresh_start = 2, max_var_thresh = 4){
  qc_files = list.files(paste0(filepath, region, "/Analysis", analysis_num , "_", region, "/"))
  qc_files = grep(paste0(region,"_(tenx_|smart_|atac_|meth_|huang_).*(qc.RDS)"), qc_files, value = TRUE)
  non_meth_files = grep("[^(meth)]",qc_files, value = TRUE)
  
  rhdf5::h5closeAll()

  hdf5_files = gsub(".RDS", ".H5",non_meth_files)
  for (i in 1:length(hdf5_files)){
    current_matrix = readRDS(paste0(filepath, region, "/Analysis", analysis_num , "_", region, "/",non_meth_files[i]))
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
  data_names = gsub("(_qc.RDS)", "",non_meth_files)
  names(hdf5_files) = data_names
  object = createLiger(raw.data = hdf5_files)
  object = normalize(object, chunk = chunk_size)
  datasets_use = grep("[tenx_|smart_]",non_meth_files)
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
  }
  message(paste0(num_genes, " genes found with var.thresh = ",var_thresh_new))
  object = scaleNotCenter(object, chunk = chunk_size)
  
  meth_files = setdiff(qc_files, non_meth_files)
  hdf5_files = gsub(".RDS", ".H5",meth_files)
  data_names = gsub("(_qc.RDS)", "",meth_files)
  
  
  for(i in 1:length(hdf5_files)){
    current_matrix = readRDS(paste0(filepath, region, "/Analysis", analysis_num , "_", region, "/",meth_file[i]))
    current_matrix = current_matrix[object@var.genes,]
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
    rhdf5::h5createDataset(hdf5_files[i],"scale.data",dims=c(nrow(current_matrix),ncol(current_matrix)),storage.mode="double", chunk = c(nrow(current_matrix),chunk_size))
    rhdf5::h5write(current_matrix, file=hdf5_files[i], name="scale.data",index = list(NULL, 1:ncol(current_matrix)))
    rhdf5::h5write(met.cell.data, file=hdf5_files[i], name="cell.data")
    rhdf5::h5closeAll()
  }
  hdf5_files = gsub(".RDS", ".H5",qc_files)
  data_names = gsub("(_qc.RDS)", "",qc_files)
  names(hdf5_files) = data_names
  rhdf5::h5closeAll()

  object_new = createLiger(hdf5_files)
  object_new@var.genes = object@var.genes
  return(object_new)
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






