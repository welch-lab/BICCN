#' Generate atlas directory structure for a region
#'
#' This function generates the required file structure for the 5 analyses, 
#' deconvolution output, and associated loom files
#'
#' @param region A string corresponding to the name of an anatomical region
#' @param desired.filepath the name of an existing directory within which to 
#' generate the file structure
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

