
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
apply_qc = function(filenames, region, analysis_num , qc_table_path, filepath = "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/"){
  if(str_sub(filepath, start = -1L, end = -1L) != "/"){
    filepath = paste0(filepath, "/")
  }
  '%notin%' = Negate('%in%')
  qc_summary = data.frame("Datatype", "Edition", "OriginalCells", "RemainingCells", "NumCellsDisgarded")
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
    if(str_contains(data.type,"_")){
      edition = str_split(data.type[[1]], "_")[[1]][2]
      data.type = str_split(data.type[[1]], "_")[[1]][1]
    }else{ edition = 0 }
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
    
    before_subset = dim(working_file)[[1]]
    working_file = working_file[,use.cells]
    
    ######################### Pre-process   
    if(data.type != "meth"){
      #If datatype is ATAC or Tenx, filter for QC cells
      qc_table = readRDS(qc_table_path)
      if (data.type == "atac" | data.type == "tenx" | data.type == "huang"){
        if (data.type == "atac"){ qc_stats = filter(qc_table, qc_table$DataType == "atac")}
        if(data.type == "tenx" | data.type == "hunag"){
          qc_stats = filter(qc_table, qc_table$DataType == "tenx")
          qc_stats = filter(qc_stats, qc_stats$Analysis == analysis_num)
        }
        nUMI_cutoff = qc_stats$nUMI
        mito_cutoff = qc_stats$mito
        ligs = createLiger(list(qc_mat = working_file))
        celldata = ligs@cell.data
        celldata$Mito = getProportionMito(ligs)
        celldata = filter(celldata, celldata$nUMI >= nUMI_cutoff)
        celldata = filter(celldata, celldata$Mito <= mito_cutoff)
        qc.matrix = working_file[, rownames(celldata)]
      }
      if (data.type == "smart"){
        qc.matrix = working_file
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
      og.dim = dim(working_file)[[1]]
      new.dim = dim(qc.matrix)[[1]]
      failing.cells = og.dim - new.dim
      if (edition == 0){
        save.filename = paste0(filepath, region, "/Analysis", analysis_num , "_", region, "/", region, "_", data.type, "_qc.RDS")
      } else {save.filename = paste0(filepath, region, "/Analysis", analysis_num , "_", region, "/", region, "_" ,data.type,"_", edition, "_qc.RDS") }
      saveRDS(qc.matrix, save.filename)
      newest_qc = c(data.type, edition, before_subset, og.dim, new.dim, failing.cells)
      names(newest_qc) = c("Datatype", "Edition", "OriginalCells", "AnalysisSubset" ,"QCPassCells", "NumCellsDisgarded")
      qc_summary = rbind(qc_summary, newest_qc)
    } 
    if( data.type == "meth") {
      #Remember we do not do methylation data for analysis 1 or 2
      if (analysis_num != 1 & analysis_num != 2){
        qc.matrix = as.matrix(max(working_file) - working_file)
        og.dim = dim(working_file)[[1]]
        new.dim = dim(qc.matrix)[[1]]
        failing.cells = og.dim - new.dim
        if (edition == 0){
          save.filename = paste0("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/", region, "/Analysis", analysis_num , "_", region, "/", region, "_",data.type, "_qc.RDS")
        } else { save.filename = paste0("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/", region, "/Analysis", analysis_num , "_", region, "/", region, "_" ,data.type,"_", edition, "_qc.RDS") }
        saveRDS(qc.matrix, save.filename)
        newest_qc = c(data.type, edition, before_subset, og.dim, new.dim, failing.cells)
        names(newest_qc) = c("Datatype", "Edition", "OriginalCells", "AnalysisSubset" ,"QCPassCells", "NumCellsDisgarded")
        qc_summary = rbind(qc_summary, newest_qc) }
      
    }
  }
  if (data.type != "meth" & data.type != "atac" & data.type != "smart" & data.type != "tenx"){
    warning("Uknown Data Type submitted", immediate. = T)
    
  }
  
  #Output a .csv to the Analysis directory
  csv_filename = paste0(filepath, region, "/Analysis", analysis_num , "_", region, "/", region, "_Overall_qc.csv")
  write.table(qc_summary, csv_filename, header = TRUE)
}