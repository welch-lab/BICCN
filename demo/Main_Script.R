#Main WorkFlow

filepath = "/scratch/welchjd_root/welchjd0/akriebel/BICCN_start"
#First initialize file directory substructure
create.directories(region = "AUD", desired.filepath = filepath)

#Next create a list of filepaths for the files related to this analysis. 
#Please name each file path with it's respective data type, i.e "tenx", "smart", "huang", "atac", or "meth"
#For datasets where there is more than instance of a datatype, please denote each dataset as "tenx_1", "tenx_2", etc.

filenames = c(tenx = "/nfs/turbo/umms-welchjd/BRAIN_initiative/iso_allen/AUD_rna_mat.RDS", smart = "/nfs/turbo/umms-welchjd/BRAIN_initiative/rna/zeng_smart/AUD_smart_mat.RDS", huang = "/nfs/turbo/umms-welchjd/akriebel/Brain_initiative/AUDv2/QC_AUD/QC_data/Huang_filtered_data_AUD.RDS", atac = "/nfs/turbo/umms-welchjd/akriebel/Brain_initiative/AUDv2/QC_AUD/QC_data/ATAC_filtered_data_AUD.RDS", meth_1 = "/nfs/turbo/umms-welchjd/akriebel/Brain_initiative/AUDv2/raw_data/mc_tsv/Methylation_AUD_removeduplicates.RDS", meth_2 = "/nfs/turbo/umms-welchjd/akriebel/Brain_initiative/AUDv2/raw_data/retro_snmc_seq/Methylation_retro_AUD_duplicatesremoved.RDS")
apply_qc(filenames, 
         "AUD", 
         1, 
         "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Base_Reference_Files/QC_table.RDS",
         filepath = filepath)
preprocess_data(filepath, 
                "AUD", 
                1, 
                1000)
run_online(filepath,
           "AUD",
           1)
annotate_by_modality(filepath, 
                "AUD", 
                1, 
                1000)



#Output QC results to main folder. Run this AFTER running all five analyses.
master_csv("AUD")



### Example deconvolution on the VIS

save_spatial_data("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region",
                              "VIS",
                              "/nfs/turbo/umms-welchjd/sodicoff/sub_mat_vis.RDS",
                              "/nfs/turbo/umms-welchjd/sodicoff/sub_coords_vis.RDS",
                              "ISH")

save_spatial_data("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region",
                              "VIS",
                              "/scratch/welchjd_root/welchjd/sodicoff/slide_seq_regions/region_mats/VIS_mat.RDS",
                              "/scratch/welchjd_root/welchjd/sodicoff/slide_seq_regions/region_coords/VIS_coords.RDS",
                              "slideseq")

deconvolve_spatial("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region",
                              "VIS",
                              "ISH",
                              n.cells = 1000,
                              deconv.gene.num = 2000,
                              gene.num.tol = 10,
                              known.annotations = NULL)

deconvolve_new_data("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region",
                              "VIS",
                              "slideseq",
                              slide.seq = TRUE)

analyze_gene_signatures("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region",
                              "VIS",
                              "ISH",
                              plot = TRUE)

generate_loading_gifs ("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region",
                              "VIS",
                              "ISH")

assign_single_cells("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region",
                              "VIS",
                              "slideseq")

generate_loading_gifs ("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region",
                              "VIS",
                              "slideseq",
                              mat.use = "assignment")

voxelize_single_cells("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region",
                              "VIS",
                              "slideseq",
                              10,
                              "~")

save_spatial_data("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region",
                  "VIS","VIS_generated_voxel_exp.RDS",
                  "VIS_generated_voxel_coords.RDS",
                  "slideseq_10")

deconvolve_new_data("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region",
                              "VIS",
                              "slideseq_10",
                              slide.seq = TRUE)

generate_loading_gifs ("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region",
                              "VIS",
                              "slideseq_10")

voxelize_single_cells("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region",
                              "VIS",
                              "slideseq",
                              20,
                              "~")

save_spatial_data("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region",
                  "VIS","VIS_generated_voxel_exp.RDS",
                  "VIS_generated_voxel_coords.RDS",
                  "slideseq_20")

deconvolve_new_data("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region",
                              "VIS",
                              "slideseq_20",
                              slide.seq = TRUE)

generate_loading_gifs ("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region",
                              "VIS",
                              "slideseq_20")

voxelize_single_cells("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region",
                              "VIS",
                              "slideseq",
                              5,
                              "~")

save_spatial_data("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region",
                  "VIS","VIS_generated_voxel_exp.RDS",
                  "VIS_generated_voxel_coords.RDS",
                  "slideseq_5")

deconvolve_new_data("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region",
                              "VIS",
                              "slideseq_5",
                              slide.seq = TRUE)

generate_loading_gifs ("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Analyses_By_Region",
                              "VIS",
                              "slideseq_5")
