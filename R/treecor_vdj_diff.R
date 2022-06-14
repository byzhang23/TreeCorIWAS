#' Differential gene usage analysis
#'
#' @param immunarch_list    A list of immune profiling info in immunarch format. Default is NULL. Users can provide their processed pseudobulk list (e.g. after covariate adjustment) via this parameter. Note that the names of list shall be matched with \code{`id`} extracted from \code{`hierarchy_list`}.
#' @param chain             specify which tcr chain: 'tra' or 'trb'.
#' @param gene  specify one or multiple queried tcr or bcr genes (use '+' to separate).
#' @param norm  an indicator to specify normalization or not. Default is TRUE (normalized to frequency).
#' @param transform specify either 'clr' (centered log ratio) transformation or 'none' (default).
#' @param hierarchy_list a hierarchy list by running \code{'extract_hrchy_string()'} or \code{'extract_hrchy_seurat()'} function.
#' @param sample_meta       Sample-level metadata, where each row is a sample. Must contain 'sample' column and additional variables to be used in the analysis, such as covariates or outcomes of interest.
#' @param response_variable A vector of response variables.
#' @param separate          A TRUE (default) or FALSE indicator, specifying how to evaluate multivariate outcomes.
#'                              \itemize{
#'                              \item TRUE: evaluate multivariate phenotype separately (it is equivalent to run this pipeline for each univariate phenotype).
#'                              \item FALSE: evaluate multivariate phenotype jointly.
#'                              }
#' @param weight            A weight matrix to combine multivariate phenotype. The dimension should be number_phenotype * 1 If none is provided, then PC1 will be used as a joint univariate phenotype.
#' @param formula           An object of class 'formula': a symbolic description of adjustment formula (i.e. only includes covariates other than response variable(s))
#' @param coef              A column number or column name specifying which coefficient to be extracted (by default: 2).
#' @param fdr_cutoff        Cutoff value for FDR. Only genes with lower FDR are listed. Default is 0.05.
#' @param filter_prop       Filter VDJ gene usage (default: 0.1)
#' @param ncores            Number of cores to be used. If ncores > 1, it will be implemented in a parallel mode.
#' @param save_as_csv       An indicator to save identified DEGs in csv files. DEGs for each cell cluster is saved as 'responsevariable_celltype_DEG.csv' and a summary file of combining DEGs from all cell clusters is saved as 'responsevariable_combinedDEG.csv'.
#' @param verbose           Show progress
#' @return A list of three elements:
#' \itemize{
#' \item diff.summary: A summary table of number of differential features for each tree node.
#' \item diff.ls: A comprehensive list of outcome(s)-associated feature set for each tree node. Use \code{`result$dge.ls$response_variable[[celltype]]`} to extract DEGs for a specific cell type
#' \item feature.ls: A list of processed feature set for each cell type
#'}
#' @export
#' @import dplyr parallel limma immunarch
#' @author Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji


treecor_vdj_diff <- function(immunarch_list, chain, gene, norm = T, transform = 'none',
                             hierarchy_list, sample_meta, response_variable, separate = T, differential_method = 'limma',weight = NULL, formula = NULL, coef = 2, fdr_cutoff = 0.05, filter_prop = 0.1,
                             ncores = parallel::detectCores(), save_as_csv = T, verbose = T){

    if (Sys.info()[['sysname']]=='Windows') {
        message('Parallel is disabled for Windows. Running with one core')
        ncores <- 1
    }

    ## Extract info
    label_info <- hierarchy_list$layout
    if(sum(duplicated(label_info$label))) stop('Duplicated cell type name.')
    leaves_info <- hierarchy_list$leaves_info
    unq_id <- unique(leaves_info$id)

    if(verbose){
        message(paste0('Limit to ', chain, ' chain...'))
    }
    dat.ls <- immunarch_list[[chain]]

    ## =============================== ##
    ## 1. Prepare feature list
    ## =============================== ##
    if(verbose){
        message('Extract gene usage...')
    }
    # pb_names <- intersect(as.character(unq_id),names(dat.ls))
    pb_names <- intersect(label_info$label,names(dat.ls))
    if(length(pb_names) != length(unq_id)){
        stop('Invalid feature list, please check if the following nodes exist: \n',
             # paste(hierarchy_list$layout$label[which(hierarchy_list$layout$id %in% as.numeric(setdiff(as.character(unq_id),pb_names)))],collapse = ',')
             paste(setdiff(label_info$label,names(dat.ls)),collapse = ',')
             )
    }
    # get gene usage
    pb.ls <- mclapply(1:length(dat.ls),function(i){
        if(verbose) message(i)
        getGeneUsage(dat.ls[[i]],gene,norm,transform,ncores,filter_prop)
    },mc.cores = ncores)
    names(pb.ls) <- names(dat.ls)

    if(!separate & length(response_variable)>1){
        # multivariate outcome; combine
        if(!is.null(weight)){
            sample_meta$combined_phenotype <- as.matrix(sample_meta[,response_variable]) %*% weight
        }else{
            # default: PC1
            sample_meta$combined_phenotype <- prcomp(sample_meta[,response_variable],scale. = T)$x[,1]
            warning('Multiple response variables are provided. Proceed with 1st PC of response variables and rename as `combined_phenotype`')
        }
        response_variable <- 'combined_phenotype'
    }

    ## =============================== ##
    ## 2. call diff_limma function
    ## =============================== ##
    res <- diff_limma(feature_list = pb.ls,
                      feature_type = gsub('\\+','\\-',gsub(' ','',gene)),
                      hierarchy_list,
                      sample_meta,
                      response_variable,
                      separate,
                      differential_method,
                      formula,
                      coef,
                      fdr_cutoff,
                      'sample',
                      ncores,
                      save_as_csv,
                      verbose)

}
