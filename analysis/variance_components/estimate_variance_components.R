library("lme4")
library("dplyr")
library("ggplot2")

#' Calculate the proportion of variance explaned by different factors in a lme4 model
varianceExplained <- function(lmer_model){
  variance = as.data.frame(lme4::VarCorr(lmer_model))
  var_percent = dplyr::mutate(variance, percent_variance = vcov/sum(vcov)) %>% 
    dplyr::select(grp, percent_variance) %>% 
    dplyr::mutate(type = "gene")
  var_row = tidyr::spread(var_percent, grp, percent_variance)
  return(var_row)  
}

#Import processed data
flow_processed = readRDS("results/processed_flow_cytometry_data.rds")
line_medatada = readRDS("data/compiled_line_metadata.rds")

#Map flow cytometry channels to specifc proteins
channel_marker_map = data_frame(channel = c("APC.A","PE.A","Pacific.Blue.A"), 
                                protein_name = c("CD206","CD16","CD14"))

unique_lines = dplyr::select(line_medatada, line_id, donor, genotype_id) %>% unique()
flow_data = dplyr::left_join(flow_processed, channel_marker_map, by = "channel") %>%
  dplyr::mutate(donor = ifelse(donor == "fpdj", "nibo",donor)) %>% #fpdj and nibo are the same donors
  dplyr::left_join(unique_lines, by = "donor") %>%
  dplyr::mutate(intensity = mean2-mean1) %>%
  dplyr::select(line_id, genotype_id, donor, flow_date, protein_name, purity, intensity)

#Construct a matrix of intensity values
intensity_matrix = dplyr::select(flow_data, line_id, genotype_id, flow_date, protein_name, intensity) %>% 
  tidyr::spread(protein_name, intensity) %>%
  dplyr::mutate(sample_id = paste(line_id, as.character(flow_date), sep = "_"))

#Make a matrix of flow data and perform PCA
flow_matrix = as.matrix(flow_df[,c(4,5,6)])
rownames(flow_matrix) = flow_df$sample_id
pca_res = prcomp(flow_matrix, scale = TRUE, center = TRUE)

#Make a PCA plot
pca_df = dplyr::mutate(as.data.frame(pca_res$x), sample_id = rownames(pca_res$x))
ggplot(pca_df, aes(x = PC1, y = PC2, label = sample_id)) + geom_point() + geom_text()

#Choose outliers based on PCA and remove them
outlier_samples = c("fafq_1_2015-10-16","iill_1_2015-10-20")
flow_df_filtered = dplyr::filter(flow_df, !(sample_id %in% outlier_samples))

#Quick variance component analysis
#Estimate variance explained by different factors
cd14_exp = lmer(CD14 ~ (1|flow_date) + (1|line_id), flow_df_filtered) %>% varianceExplained()
cd16_exp = lmer(CD16 ~ (1|flow_date) + (1|line_id), flow_df_filtered) %>% varianceExplained()
cd206_exp = lmer(CD206 ~ (1|flow_date) + (1|line_id), flow_df_filtered) %>% varianceExplained()

#Identify donors that have more than one measurement
replicated_donors = dplyr::group_by(flow_df_filtered, line_id) %>% 
  dplyr::summarise(n_replicates = length(line_id)) %>% 
  dplyr::filter(n_replicates > 1)
flow_df_replicated = dplyr::filter(flow_df_filtered, line_id %in% replicated_donors$line_id)

#Repeat the variance component analysis
cd14_exp = lmer(CD14 ~ (1|flow_date) + (1|line_id), flow_df_replicated) %>% varianceExplained()
cd16_exp = lmer(CD16 ~ (1|flow_date) + (1|line_id), flow_df_replicated) %>% varianceExplained()
cd206_exp = lmer(CD206 ~ (1|flow_date) + (1|line_id), flow_df_replicated) %>% varianceExplained()



