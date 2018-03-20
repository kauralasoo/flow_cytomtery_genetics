#' A wrapper around MatrixeQTL
#'
#' Prepares SlicedData objects of expression data, genotypes and covariates 
#' and then runs MatrixeQTL on the data.
#' 
#' @param exp_data Matrix of gene expression data (genes in rows, samples in columns).
#' @param geno_data Matrix of genotype data (SNPs in rows, samples in columns).
#' @param snpspos Matrix of SNP coordinates (columns: snpid,chr,pos).
#' @param genepos Matrix of gene coordinates (columns: geneid,chr,left,right).
#' @param covariates Matrix of covariates (samples in columns).
#' @param cisDist cis distance from the gene.
#' @param pvOutputThreshold Maximum p-value to report. Smaller values make the code quicker.
#' @param permute Permute genotype labels before qtl mapping.
#' @param model Specifies which MatrixEQTL model to use. Options: modelLINEAR, modelLINEAR_CROSS and modelANOVA.
#' @return MatrixeQTL result object.
#' @author Kaur Alasoo
#' @export 
runMatrixEQTL <- function(exp_data, geno_data, snpspos, genepos, covariates = NULL, 
                          cisDist = 5e5, pvOutputThreshold = 1e-2, permute = FALSE, model = modelLINEAR){
  #Run matrixeQTL on a prepared data set
  
  #Perform some sanity checks
  if(!all(colnames(exp_data) == colnames(geno_data))){
    stop("Column names of expression and genotype data are not equal.")
  }
  
  #Construct a SlicedData object of the expression data
  expression_sliced = SlicedData$new()
  expression_sliced$CreateFromMatrix(exp_data)
  expression_sliced$ResliceCombined(sliceSize = 2000)
  
  #Create a SlicedData obejct for the genotypes
  if(permute == TRUE){
    #Permute column labels of the genotype data
    genotype_labels = colnames(geno_data)
    geno_data = geno_data[,sample(length(genotype_labels), length(genotype_labels))]
    colnames(geno_data) = genotype_labels
  }
  snps = SlicedData$new()
  snps$CreateFromMatrix(geno_data)
  snps$ResliceCombined()
  
  #Add covariates
  cvrt = SlicedData$new()
  if (!is.null(covariates)){
    if(!all(colnames(exp_data) == colnames(covariates))){
      stop("Column names of expression and covariates data are not equal.")
    }
    cvrt$CreateFromMatrix(covariates)
    cvrt$ResliceCombined()
  }
  
  #RUN
  me = Matrix_eQTL_main(
    snps = snps,
    gene = expression_sliced,
    cvrt = cvrt,
    output_file_name = "",
    pvOutputThreshold = 0,  
    output_file_name.cis = NULL,
    pvOutputThreshold.cis = pvOutputThreshold,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    useModel = model, 
    errorCovariance = numeric(), 
    verbose = TRUE,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
  
  return(me)
}
