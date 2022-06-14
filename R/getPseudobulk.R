#' Construct sample-level pseudobulk gene expression matrix
#'
#' This function takes in raw count gene expression matrix and cell-level metadata as inputs, and generates sample-level normalized pseudobulk gene expression matrix.
#' @param expr              A raw count gene expression matrix
#' @param sample_meta       sample=level metadata
#' @param cell_meta         A data frame for cell-level metadata, where each row is a cell. Must contain these columns: barcode, celltype and sample.
#' @param columns           specify column(s) name to pool pseudobulk (default is sample)
#' @param filter_prop       A number ranges from 0 to 1, to filter low expressed genes across samples (by default: 0.1). Genes with at least this proportion of samples with log2-normalized count greater than 0.01 are retained.
#' @return A sample-level normalized pseudobulk gene expression matrix with genes in the row and samples in the columns
#' @importFrom Matrix Matrix colSums t
#' @export
#' @author Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji
#'

getPseudobulk <- function(expr, cell_meta, sample_meta, columns = 'sample', filter_prop = 0.1){

    barcode <- cell_meta$barcode
    comb <- suppressMessages(inner_join(sample_meta,cell_meta))

    sample_name <- comb[,columns]
    sample_name <- factor(sample_name,levels = unique(sample_name))
    expr_subset <- expr[,barcode]
    expr_subset <- Matrix::Matrix(expr_subset,sparse = T)

    if(length(levels(sample_name))>1){
        mod <- model.matrix(~ 0 + sample_name)
        colnames(mod) <- levels(sample_name)
        ct <- expr_subset %*% mod
    }else{
        ct <- expr_subset %*% matrix(1,nrow = ncol(expr_subset),ncol = 1,dimnames = list(NULL,unique(levels(sample_name))))
    }

    ## filter out genes with 0 count (CPM normalization)
    rc <- Matrix::colSums(ct)
    # modify: scaling factor: median library size
    rc <- rc/median(rc,na.rm = T) # 1e6
    denom0 <- which(rc==0)
    if(length(denom0)>0){
        warning(paste0('There are ', length(denom0),'genes having zero expression across all samples. Proceed with excluding these genes.'))
        ct <- ct[,-denom0,drop=F]
        rc <- rc[-denom0]
    }
    res <- log2(Matrix::t(Matrix::t(ct)/rc + 1)) %>% as.matrix

    ## filter out low express genes (filter_prop)
    res <- res[rowMeans(res > 0.01) >= filter_prop,,drop = F]
}
