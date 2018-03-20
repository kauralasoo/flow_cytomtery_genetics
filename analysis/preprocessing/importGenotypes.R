library("SNPRelate")
library("GWASTools")

#' Import genotypes in GDS file to a matrix
#'
#' Open GDS file (created with SNPRelate::snpgdsVCF2GDS) into R and convert it into a matrix of 
#' variant positions and matrix of allele dosage (0,1,2)
#' 
#' @param file Path to the GDS file.
#' @return List containing snpspos and genotypes matrix.
#' @author Kaur Alasoo
#' @export 
gdsToMatrix <- function(gds_file){
  
  #Extract genotypes
  gds <- GWASTools::GdsGenotypeReader(gds_file)
  genotypes = GWASTools::getGenotype(gds)
  sample_ids = GWASTools::getVariable(gds, "sample.id")
  snp_rs_ids = GWASTools::getVariable(gds, "snp.rs.id")
  snp_ids = GWASTools::getVariable(gds, "snp.id")
  
  #Invent id for snps that do not have an rs id
  new_snp_ids = paste("snp",snp_ids[snp_rs_ids == ""], sep = "")
  snp_rs_ids[snp_rs_ids == ""] = new_snp_ids
  colnames(genotypes) = sample_ids
  rownames(genotypes) = snp_rs_ids
  
  #Extract SNP coordinates
  snpspos = dplyr::data_frame(snpid = snp_rs_ids, 
                              chr = GWASTools::getVariable(gds, "snp.chromosome"), 
                              pos = GWASTools::getVariable(gds, "snp.position"))
  GWASTools::close(gds)
  return(list(snpspos = snpspos, genotypes = genotypes))
}

#Convert VCF to GDS
snpgdsVCF2GDS("data/genotypes/open_access_genotypes.vcf.gz", 
              "data/genotypes/open_access_genotypes.gds", method="biallelic.only")
genotypes = gdsToMatrix("data/genotypes/open_access_genotypes.gds")
saveRDS(genotypes, "data/genotypes/open_access_genotypes.rds")

