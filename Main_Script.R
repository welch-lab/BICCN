#Main WorkFlow

#First initialize file directory substructure
create.directories(region = "ORB")

#Next create a list of filepaths for the files related to this analysis. 
#Please name each file path with it's respective data type, i.e "tenx", "smart", "atac", or "meth"
#For datasets where there is more than instance of a datatype, please denote each dataset as "tenx_1", "tenx_2", etc.

filenames = c(tenx_1 = "/nfs/turbo/umms-welchjd/BRAIN_initiative/iso_allen/AUD_rna_mat.RDS", smart = "/nfs/turbo/umms-welchjd/BRAIN_initiative/rna/zeng_smart/AUD_smart_mat.RDS", tenx_2 = "/nfs/turbo/umms-welchjd/akriebel/Brain_initiative/AUDv2/QC_AUD/QC_data/Huang_filtered_data_AUD.RDS", atac = "/nfs/turbo/umms-welchjd/akriebel/Brain_initiative/AUDv2/QC_AUD/QC_data/ATAC_filtered_data_AUD.RDS")
apply_qc(filenames, region = "ORB", analysis_num = 1, qc_table_path = "/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses/Base_Reference_Files/QC_table.RDS", filepath = "/scratch/welchjd_root/welchjd0/akriebel/Test")
object = preprocess_data("/nfs/turbo/umms-welchjd/BRAIN_initiative/BICCN_integration_Analyses", "ORB", 1, 1000)



#Output QC results to main folder
master_csv("AUD")
