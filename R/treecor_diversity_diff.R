#' Differential diversity
#'
#' @param immunarch_list    A list of immune profiling info in immunarch format. Default is NULL. Users can provide their processed pseudobulk list (e.g. after covariate adjustment) via this parameter. Note that the names of list shall be matched with \code{`id`} extracted from \code{`hierarchy_list`}.
#' @param chain             specify which tcr chain: 'tra' or 'trb' or 'tra+trb'.
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
#' @param ncores            Number of cores to be used. If ncores > 1, it will be implemented in a parallel mode.
#' @param save_as_csv       An indicator to save identified DEGs in csv files. DEGs for each cell cluster is saved as 'responsevariable_celltype_DEG.csv' and a summary file of combining DEGs from all cell clusters is saved as 'responsevariable_combinedDEG.csv'.
#' @param verbose           Show progress
#' @return A list of three elements:
#' \itemize{}
#' @export
#' @import dplyr parallel limma immunarch
#' @author Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji

treecor_diversity_diff <- function(immunarch_list, chain, hierarchy_list, sample_meta, response_variable, separate = T, differential_method = 'limma', weight = NULL, formula = NULL, coef = 2, fdr_cutoff = 0.05,
                                   method = 'chao1', max.q = 6, min.q = 1,q = 5, step = NA, quantile = c(0.025, 0.975), extrapolation = NA, perc = 50, norm = T, do.norm = NA, laplace = 0,
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
        message('Calculate diversity for each node...')
    }
    pb_names <- intersect(label_info$label,names(dat.ls))
    if(length(pb_names) != length(unq_id)){
        stop('Invalid feature list, please check if the following nodes exist: \n',
             # paste(hierarchy_list$layout$label[which(hierarchy_list$layout$id %in% as.numeric(setdiff(as.character(unq_id),pb_names)))],collapse = ',')
             paste(setdiff(label_info$label,names(dat.ls)),collapse = ',')
        )
    }
    # get diversity
    pb.ls <- mclapply(1:length(dat.ls),function(i){
        if(verbose) message(i)

        # exclude samples with NULL
        dt <- dat.ls[[i]]
        div_data <- repDiversity(.data = dt[which(sapply(dt,function(x) !is.null(x)))],
                                 .method = method,
                                 .col = "aa",
                                 .max.q = max.q,
                                 .min.q = min.q,
                                 .q = q,
                                 .step = step,
                                 .quantile = quantile,
                                 .extrapolation = extrapolation,
                                 .perc = perc,
                                 .norm = norm,
                                 .verbose = verbose,
                                 .do.norm = do.norm,
                                 .laplace = laplace)

        ## different methods - output data frame slightly different (need a function)
        if(method == 'chao1'){
            res <- div_data[,1,drop = F]
            colnames(res) <- 'chao1'
        }else if(method == 'hill'){
            res <- div_data %>% tidyr::spread(Q,Value)
            rownames(res) <- res$Sample
            res <- res[,-1,drop = F]
            colnames(res) <- paste0('hill_Qvalue',colnames(res))
        }else if(method %in% c('div','gini.simp','inv.simp')){
            res <- div_data[,-1,drop = F]
            rownames(res) <- div_data$Sample
            colnames(res) <- method
        }else if(method == 'gini'){
            res <- div_data
            colnames(res) <- 'gini'
        }else if(method %in% c('d50','dxx')){
            res <- div_data[,1,drop = F]
            colnames(res) <- paste0(method,'_',colnames(res))
        }else{
            stop('Other methods are not supported in this version.')
        }
        t(res)
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

    res <- diff_limma(feature_list = pb.ls,
                      feature_type = method,
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
                      verbose) # 1/0 represent diversity different across groups

}
