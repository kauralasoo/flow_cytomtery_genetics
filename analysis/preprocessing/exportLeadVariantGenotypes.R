flow_processed = readRDS("../../results/processed_flow_cytometry_data.rds")
line_medatada = readRDS("data/compiled_line_metadata.rds")

open_access_donors = dplyr::filter(line_medatada, open_access == 1) %>% dplyr::select(genotype_id) %>%
  dplyr::distinct()
write.table(open_access_donors, "data/open_access_donors.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

geno = gdsToMatrix("../macrophage-gxe-study/flow/genotypes/flow_cis_regions.gds")
cd14_lead = geno$genotypes["rs778587",]
cd14_lead_genotype = data_frame(genotype_id = names(cd14_lead), rs778587 = cd14_lead) %>% 
  distinct()
write.table(cd14_lead_genotype, "data/cd14_lead_variant.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
