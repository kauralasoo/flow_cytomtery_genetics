library("flowCore")
library("flowViz")
library("tidyr")
library("dplyr")
library("devtools")
source("analysis/functions/flow_functions.R")

#Import line metadata
donor_line_map = readRDS("data/compiled_line_metadata.rds") %>% 
  dplyr::select(donor, genotype_id) %>% unique()

#Construct a flow metadata data.frame:
file_names = data.frame(file_name = list.files("data/fcs_files/"), stringsAsFactors = FALSE)
metadata = tidyr::separate(file_names, file_name, c("name", "suffix"), sep = "\\.", remove = FALSE) %>% 
  tidyr::separate(name, into = c("donor", "date", "staining"), sep = "_", remove = FALSE) %>% 
  dplyr::mutate(sample = paste(donor, date, sep = "_")) %>% 
  dplyr::select(name, sample, donor, date, staining, file_name) %>%
  dplyr::left_join(donor_line_map, by = "donor")

#Save sample metadata to disk for uploading
meta_export = dplyr::select(metadata, -name) %>% 
  dplyr::filter(donor != "piun-SN")
write.table(meta_export, "data/flow_sample_metadata.txt", quote = FALSE, row.names = FALSE, sep = "\t")

#Read all FCS files into a list
flow_files_list = as.list(paste("data/fcs_files/", metadata$file_name, sep =""))
names(flow_files_list) = metadata$name
fcs_list = lapply(flow_files_list, read.FCS, alter.names = TRUE)

#Rename columns in FCS files
name_map = data.frame(old_names = c("X670.14..640..A","X450.50..405..A","X586.15..561..A"), 
                      new_names = c("APC.A","Pacific.Blue.A", "PE.A"), stringsAsFactors = FALSE)
fcs_list_renamed = lapply(fcs_list, changeFlowFrameNames, name_map)

#Keep only selected columns and merge into flowSet
selected_columns = c("FSC.A","SSC.A","APC.A","Pacific.Blue.A","PE.A","Time")
fcs_list_selected = lapply(fcs_list_renamed, function(flowframe, sel_cols){return(flowframe[,sel_cols])}, selected_columns)
flow_set = as(fcs_list_selected, "flowSet")
flow_set = addMetadataToFlowset(flow_set, metadata)

#Save the flow set to disk
saveRDS(flow_set, "data/combined_flowSet.rds")
