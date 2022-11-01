#' Generate atlas directory structure for a region
#'
#' This function generates the required file structure for the 5 analyses,
#' deconvolution output, and associated loom files
#'
#' @param region A string corresponding to the name of an anatomical region
#' @param desired.filepath the name of an existing directory within which to
#'    generate the file structure
#'
#' @return nothing
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }
#'
create.directories = function(region = "X",
                              desired.filepath = "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/"){
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

#' Generate an xlsx file to track QC across analysis levels
#'
#' This function generates an xlsx file with a sheet for each analysis to track
#' differences in QC based on varying thresholds
#'
#' @param region The anatomical region with which the analysis is associated
#' @param filepath Path to directory within which the atlas structure was generated
#'
#' @return Nothing
#'
#' @import openxlsx
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }
#'

master_csv = function(region,
                      filepath = "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/"){
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

#' Apply quality control to provided data
#'
#' This function utilizes preset filters on sequencing counts as well at
#' mitochondrial and cytoplasmic gene expression to subset provided data to only
#' high quality samples
#'
#' @param filenames A list of RDS files containing single cell data
#' @param region The anatomical region with which the samples are associated
#' @param analysis_num What level of the analysis to apply QC to.
#' @param qc_table_path A file path to an RDS object providing QC values
#' @param filepath_cytoplasmic A file path to a CSV supporting QC based on
#'    marker gene expression
#' @param filepath Path to directory within which the atlas structure was generated
#' @param verbose A logical, if the function should display what file is being processed
#'
#' @return Nothing
#'
#' @import
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }
#'
#'
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

#' Generate QCed expression matrices for downstream analysis and save related data
#' @param filenames a named list of file paths to RDS objects containing eexpression matrices
#' @param region A string corresponding to the name of an anatomical region
#' @param analysis_num What level of the analysis to apply QC to.
#' @param qc_table_path A string, the location of an object containing specific QC thresholds
#' @param filepath_cytoplasmic A string, the location of a CSV containing a list of cytoplasmic genes
#'    generation
#' @param filepath Path to directory within which the atlas structure was generated
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
#'
apply_qc = function(filenames,
                    region,
                    analysis_num ,
                    qc_table_path,
                    filepath_cytoplasmic = "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Base_Reference_Files/Cytoplasmic_genes.csv",
                    filepath = "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/",
                    verbose = TRUE){
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
      if(data.type %notin% c("meth","atac")){
        results_filename = paste0(filepath,region, "/Analysis1_", region, "/Analysis1_", region, "_Results_Table.RDS")
        results = readRDS(results_filename)
        results = filter(results, results$highRAnnotations =="NonN")
        results = subset(results, results$Barcode %in% colnames(working_file))
        use.cells = results$Barcode
      } else {use.cells = colnames(working_file)}
    }
    if(analysis_num == 3){
      if(data.type %notin% c("meth","atac")){
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
    
    #if (data.type == "atac"){
    #  qc_table = readRDS(qc_table_path)
    #  #Find the appropriate cutoff to apply
    #  if (data.type == "atac"){ qc_stats = filter(qc_table, qc_table$DataType == "atac")}
    #  nUMI_cutoff = qc_stats$nUMI
    #  mito_cutoff = qc_stats$mito
    #  ligs = createLiger(list(qc_mat = working_file))
    #  celldata = ligs@cell.data
    #  celldata$Mito = getProportionMito(ligs)
    #  celldata = filter(celldata, celldata$nUMI >= as.numeric(nUMI_cutoff))
    #  remaining_nUMI = rownames(celldata)
    #  celldata = filter(celldata, celldata$Mito <= as.numeric(mito_cutoff))
    #  remaining_mito = rownames(celldata)
    #  before_subset = dim(working_file)[[2]] #Gets you original dimensions of matrix
    #  working_file = working_file[,use.cells] #Gets you a matrix subset for the appropriate cell population
    #  after_subset = dim(working_file)[[2]]  #Gets you dimensions of matrix after subsetting for the appropriate cell population
    #  cells_after_subset = colnames(working_file) #Gets you the cell names present after subsetting for cell populations
    #  sub_cells = subset(rownames(celldata), rownames(celldata) %in% colnames(working_file))
    #  qc.matrix = working_file[, sub_cells]
    #  remaining_GeneCounts = remaining_cyto = colnames(qc.matrix)
    #  cytoplasmic_cutoff = GeneCount_cutoff = NA
    #}
    
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
    
    if (data.type %notin% c("meth","atac")){
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
    
    if( data.type == "meth"){
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
        qc_summary = rbind(qc_summary, newest_qc) }}
    if( data.type == "atac") {
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
      qc_summary = rbind(qc_summary, newest_qc) }
    
    if(verbose){
      message("Done Processing ", data.type, " ", edition)
    }
    
    if (data.type != "meth" & data.type != "atac" & data.type != "smartseq" & data.type != "sc10Xv3" & data.type != "sn10Xv3" & data.type != "sc10Xv2"){
      warning("Uknown Data Type submitted", immediate. = T)
    }}
  #Output a .csv to the Analysis directory
  csv_filename = paste0(filepath, region, "/Analysis", analysis_num , "_", region, "/", region, "_Overall_qc.csv")
  colnames(qc_summary) = c("Datatype", "Edition", "OriginalDimensions", "AnalysisDimensions", "FailnUMI", "nUMIThreshold", "FailMito", "MitoThreshold", "FailCyto", "CytoThreshold", "Lost for Gene Counts", "Gene Count Cutoff", "FinalDimensions")
  write.csv(qc_summary, csv_filename)
}



#' Calculates a cytoplasmic score, used in determination of sample quality
#'
#' @param data an expression matrix
#' @param qc.gene.FN a filepath for a csv containing a list of genes to use
#'    in score calculation
#'
#' @return a numeric corresponding to the cytoplasmic score
#'
#' @import scrattch.hicat
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }
#'

cal_qc_score <- function(data, qc.gene.FN="AIBS_qc_genes_10X.csv") {

  qc.genes = read.csv(qc.gene.FN)[, "Gene"]
  common.qc.genes = intersect(qc.genes, row.names(data))
  norm.dat = scrattch.hicat::logCPM(data[common.qc.genes,])
  qc.score = apply(norm.dat, 2, sum)
  #qc.score = Matrix::colSums(norm.dat[common.qc.genes,])

  return(qc.score)
}


#' Assign cells to clusters based on their maximum metagene loading
#'
#' @param object a liger object, having been quantile normalized
#'
#' @return a liger object, with a modified cluster slot
#'
#' @import scrattch.hicat
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }
#'

max_factor_assignment = function(object){
  h_factors = object@H.norm
  max_assignments = apply(h_factors, 1, which.max)
  object@clusters = as.factor(max_assignments)
  return(object)
}

#' Complete a preliminary single cell subanalysis for a region, from
#' preprocessing through preliminary clustering
#'
#' @param filepath Path to directory within which the atlas structure was generated
#' @param region The anatomical region with which the analysis is associated
#' @param analysis_num A number, 1:5, corresponding to which subanalysis to complete
#' @param chunk_size size of chunks in hdf5 files
#' @param num_genes The number of genes to utilize for gene signature
#'   calculation
#' @param gene_num_tolerance The maximum difference between deconv.gene.num and the
#'   true number of genes used
#' @param var_thresh_start
#' @param max_var_thresh
#' @param customGeneList
#' @param return.object a logical, if the processed object should be returned by
#'    the function
#' @param qn_ref a string, corresponding to which dataset to utilize as a
#'    reference during quantile normalization
#' @param knownAnnotations a filepath to a spreadsheet containing multiple levels of provided
#'   annotations
#' @param MaxFactor a logical, whether to use max factor clustering
#'
#' @return Nothing
#'
#' @import rhdf5, Matrix
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }
#'
preprocess_and_run = function(filepath, region, analysis_num, chunk_size, num_genes = 2500, gene_num_tolerance = 100, var_thresh_start = 2, max_var_thresh = 4, customGeneList = NA, return.object = FALSE, qn_ref = NA, knownAnnotations = "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Base_Reference_Files/Reference_Annotations_updated_with_MacoskoLabels.RDS", MaxFactor = FALSE, labels = TRUE, manualH5s = c(), numFactors = 30, qn_k = 20, randomSeed = 45){
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
  if (length(manualH5s) >= 1){
    nonmeth_add = grep("meth", manualH5s, value = TRUE, invert = TRUE)
    data_names = c(data_names, nonmeth_add)
    for (item in 1:length(nonmeth_add)){
      nonmeth_add[[item]] = paste0(filepath, "/", region, "/Analysis", analysis_num , "_", region, "/", nonmeth_add[[item]])
    }
    data_names = sort(data_names)
    hdf5_files = c(hdf5_files, nonmeth_add)
    hdf5_files = sort(hdf5_files)
  }
  
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
  
  if (length(manualH5s) >= 1){
    nonmeth_add = grep("meth", manualH5s, value = TRUE, invert = TRUE)
    data_names = c(data_names, nonmeth_add)
    for (item in 1:length(nonmeth_add)){
      nonmeth_add[[item]] = paste0(filepath, "/", region, "/Analysis", analysis_num , "_", region, "/", nonmeth_add[[item]])
    }
    data_names = sort(data_names)
    hdf5_files = c(hdf5_files, nonmeth_add)
    hdf5_files = sort(hdf5_files)
  }
  
  names(hdf5_files) = data_names
  rhdf5::h5closeAll()
  
  var.genes = object@var.genes
  rm(object)
  rhdf5::h5closeAll()
  
  object_new = createLiger(as.list(hdf5_files))
  object_new@var.genes = var.genes
  object_new = online_iNMF(object_new, k = numFactors , lambda = 5, max.epochs = 20, seed = randomSeed)
  if (!is.na(qn_ref)){
    object_new = quantile_norm(object_new, do.center = T, ref_dataset = qn_ref, knn_k = qn_k, rand.seed = randomSeed)
  } else{
    object_new = quantile_norm(object_new, do.center = T, knn_k = qn_k, rand.seed = randomSeed)
  }
  liger_name = paste0(filepath, "/", region, "/Analysis", analysis_num, "_", region, "/onlineINMF_",region, "_object.RDS" )
  print("Saving LIGER object")
  object_new = runUMAP(object_new,  n_neighbors=30, min_dist=0.3, distance ="cosine", rand.seed = randomSeed)
  #object_new = readSubset(object_new)
  saveRDS(object_new, liger_name)
  max_factor_assignment = function(object){
    h_factors = object@H.norm
    max_assignments = apply(h_factors, 1, which.max)
    object@clusters = as.factor(max_assignments)
    return(object)
  }
  
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
    liger_low = louvainCluster(object_new, resolution = low_resolution, k = 200, random.seed = randomSeed)
    liger_high = louvainCluster(object_new, resolution = high_resolution, k = 200, random.seed = randomSeed)
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
    full2 = as.factor(full2)
    liger_low@clusters = full2
    result$Barcode = rownames(result)
    result = left_join(result, full)
    result = select(result, -c(clusts))
    rownames(result) = result$Barcode
    result = select(result, -c(Barcode))
    
  }
  if(MaxFactor == FALSE){
    result$ann = annies[names(clusts)]
    liger_low@clusters = as.factor(annies[names(clusts)])
  }
  
  if (labels == TRUE){
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
  }
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

#' annotate new clusters using reference labels for a subset of cells
#'
#' @param object a liger object, with a filled clusters slot
#' @param annotations a named vector of annotations for a subset of cells
#' @param k how many nearest neighbors to utilize for cell type assignment
#'
#' @return a liger object with clusters overwritten
#'
#' @import RANN
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }
#'

transfer_labels = function(object,
                           annotations,
                           k = 20){
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

## TODO is this too dumb to use. nowhere to put it in the analysis

#' ana
#'
#' @param data A
#' @param qc.gene.FN
#'
#' @return Nothing
#'
#' @import scrattch.hicat
#'
#' @export
#' @examples
#' \dontrun{
#'
#' }
#'

annotate_by_modality = function(filepath,
                                region,
                                analysis_num,
                                chunk_size,
                                knownAnnotations = "/nfs/turbo/umms-welchjd/BRAIN_initiative/Final_integration_workflow/SupportFiles/Known_Annotationsv3.rds",
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
generate_markergenes = function(region, analysis_num, filepath  ="/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region/", added_genes = c(), baseMarkers = TRUE){
  #Get a list of all relevant files
  direc = paste0(filepath, region, "/", "Analysis", analysis_num, "_", region)
  dirfiles = list.files(direc)
  dirfiles = grep("qc.RDS", dirfiles, value = TRUE)
  #Remove methylation matrices as they are typically uninformative
  #dirfiles = grep("meth", dirfiles, value = TRUE, invert = TRUE)
  methfiles = grep("meth", dirfiles, value = TRUE)
  #Read in files
  matrix_list = c()
  for (rfile in 1:length(dirfiles)){
    fn = paste0(direc, "/", dirfiles[[rfile]])
    mm = readRDS(fn)
    mm = as(mm, "dgCMatrix")
    matrix_list[[rfile]] = mm
    mname = sub("_qc.RDS", "", dirfiles[[rfile]])
    names(matrix_list)[[rfile]] = mname
  }
  #Read in the LIGER object
  obj = readRDS(paste0(direc, "/onlineINMF_", region, "_object.RDS"))
  #Create a useable LIGER object
  ligs = createLiger(matrix_list, remove.missing = FALSE)
  ligs = normalize(ligs, remove.missing = FALSE)
  ligs = selectGenes(ligs)
  ligs@var.genes = obj@var.genes
  if(length(methfiles) > 1){
    for (mf in 1:length(methfiles)){
      ligs@norm.data[[mf]] = matrix_list[[mf]]
    }
  }
  ligs = scaleNotCenter(ligs)
  
  if(dim(ligs@H.norm)[1] != dim(ligs@cell.data)[1]){
    bars_we_have = rownames(ligs@cell.data)
    for (k in 1:length(obj@H)){
      rel_bars = subset(bars_we_have, bars_we_have %in% rownames(obj@H[[k]]))
      obj@H[[k]] = obj@H[[k]][rel_bars,]
    }
    obj@H.norm = obj@H.norm[bars_we_have,]
    obj@tsne.coords = obj@tsne.coords[bars_we_have,]
    obj@clusters = obj@clusters[bars_we_have]
  }
  
  ligs@H = obj@H
  ligs@H.norm = obj@H.norm
  ligs@V = obj@V
  ligs@W = obj@W
  ligs@tsne.coords = obj@tsne.coords
  ligs@clusters = obj@clusters
  
  
  
  
  
  
  #Read in desired genes
  #If base files is TRUE, generate a pdf of the Marker Genes
  if(baseMarkers == TRUE){
    genes_oi = readRDS("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Base_Reference_Files/Marker_genes_vector.RDS")
    pdf_name1 = paste0(direc, "/Images/Analysis", analysis_num, "_", region, "_BaseMarkerPlots.pdf" )
    pdf(pdf_name1)
    for (gene in genes_oi){
      plotGene(ligs, gene)
    }
    dev.off()
  }
  #If added genes is true, plot these additional genes into a seperate pdf
  if (length(added_genes) != 0){
    pdf_name1 = paste0(direc, "/Images/Analysis", analysis_num, "_", region, "_AddedMarkerPlots.pdf" )
    pdf(pdf_name1)
    for (gene in added_genes){
      plotGene(ligs, gene)
    }
    dev.off()
  }
}

