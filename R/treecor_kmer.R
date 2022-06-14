#' treecor_kmer
#' @param immunarch_list    A list of immune profiling info in immunarch format. Default is NULL. Users can provide their processed pseudobulk list (e.g. after covariate adjustment) via this parameter. Note that the names of list shall be matched with \code{`id`} extracted from \code{`hierarchy_list`}.
#' @param chain             specify which tcr chain: 'tra' or 'trb'.
#' @param k                 kmer
#' @param unq_kmer_per_seq  an indicator to specify whether to count only unique kmer
#' @param remove_residual   specify whether remove residuals and keep only contact region
#' @param start_end_pos     if remove residual is true, then specify starting and ending positions. Default is c(4,-4), which indicates to start at 4th amino acid and ends at 4th from the end.
#' @param norm              normalize to frequency
#' @param hierarchy_list    a hierarchy list by running \code{'extract_hrchy_string()'} or \code{'extract_hrchy_seurat()'} function.
#' @param sample_meta       Sample-level metadata, where each row is a sample. Must contain 'sample' column and additional variables to be used in the analysis, such as covariates or outcomes of interest.
#' @param response_variable A vector of response variables.
#' @param separate          A TRUE (default) or FALSE indicator, specifying how to evaluate multivariate outcomes.
#'                              \itemize{
#'                              \item TRUE: evaluate multivariate phenotype separately (it is equivalent to run this pipeline for each univariate phenotype).
#'                              \item FALSE: evaluate multivariate phenotype jointly.
#'                              }
#' @param differential_method   Specify which type of differential methods: limma (default) or gam
#' @param weight            A weight matrix to combine multivariate phenotype. The dimension should be number_phenotype * 1 If none is provided, then PC1 will be used as a joint univariate phenotype.
#' @param formula           An object of class 'formula': a symbolic description of adjustment formula (i.e. only includes covariates other than response variable(s))
#' @param coef              A column number or column name specifying which coefficient to be extracted (by default: 2).
#' @param fdr_cutoff        Cutoff value for FDR. Only genes with lower FDR are listed. Default is 0.05.
#' @param ncores            Number of cores to be used. If ncores > 1, it will be implemented in a parallel mode.
#' @param save_as_csv       An indicator to save identified DEGs in csv files. DEGs for each cell cluster is saved as 'responsevariable_celltype_DEG.csv' and a summary file of combining DEGs from all cell clusters is saved as 'responsevariable_combinedDEG.csv'.
#' @param verbose           Show progress
#' @return A list of four elements:
#' \itemize{
#' \item diff.summary: A summary table of number of differential features for each tree node.
#' \item diff.ls: A comprehensive list of outcome(s)-associated feature set for each tree node. Use \code{`result$dge.ls$response_variable[[celltype]]`} to extract DEGs for a specific cell type
#' \item feature.ls: A list of processed feature set for each cell type
#'}
#' @export
#' @import dplyr parallel limma immunarch
#' @author Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji
#'

treecor_kmer <- function(immunarch_list, chain, k = 5, unq_kmer_per_seq = T,
                         remove_residual = T, start_end_pos = c(4,-4), norm = T,
                         hierarchy_list, sample_meta, response_variable, separate = T, differential_method = 'limma',weight = NULL, formula = NULL, coef = 2, fdr_cutoff = 0.05,
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

    if(remove_residual){
        if(verbose) message('Remove residuals and keep contact regions...')
        start <- start_end_pos[1]
        end <- start_end_pos[2]
        if(start<0|end>0) stop('Incorrect start or end position. please check again.')

        dat.ls <- parallel::mclapply(dat.ls,function(dt){
            getContactRegion(dt,start,end)
        },mc.cores = ncores)
    }

    if(verbose){
        message('Count kmer per sequence...')
    }

    # if(!unq_kmer_per_seq){
    #     # immunarch
    #     kmer.ls <- mclapply(1:length(dat.ls),function(i){
    #         if(verbose) message(paste0('node ',i))
    #         dt <- dat.ls[[i]]
    #         if(!is.null(dt)){
    #             dt <- dt[sapply(dt,function(x) !is.null(x))]
    #             dk <- getKmers(dt,k,'aa',T)
    #             nk <- as.matrix(dk[,-1,drop = F])
    #             nk[is.na(nk)] <- 0
    #
    #             if(norm){
    #                 nk <- t(t(nk)/colSums(nk))
    #             }
    #             rownames(nk) <- dk[,1] %>% unlist
    #             colnames(nk) <- colnames(dk)[-1]
    #
    #         }else{
    #             nk <- NULL
    #         }
    #         nk
    #     },mc.cores = ncores)
    # }else{

    ##### ++++ SPEED UP: ONLY KMER ON CHILDREN NODE -> COMBINE INTO THE INTERMEDIATE NODE
    kmer.full.ls <- mclapply(1:length(dat.ls),function(i){
        if(verbose) message(paste0('node ',i))
        dt <- dat.ls[[i]]
        if(!is.null(dt)){

            dt <- dt[sapply(dt,function(x) !is.null(x))]
            if(chain == 'tra+trb'){
                full <- getUnqKmers(dt,k,'aa',F,ncores)
            }else{
                full <- getUnqKmers(dt,k,'aa',T,ncores)
            }
            # obtain matrix
            dk <- full$mat
            nk <- as.matrix(dk[,-1,drop = F])
            nk[is.na(nk)] <- 0
            ttl_sum <- sapply(dt,function(x) sum(x$Clones))

            if(norm){
                nk <- t(t(nk)/ttl_sum)
            }
            rownames(nk) <- dk[,1] %>% unlist
            colnames(nk) <- colnames(dk)[-1]

            full_kmer_table <- full$full_kmer

        }else{
            nk <- NULL
            full_kmer_table <- NULL
        }
        list(mat = nk,
             full_kmer_table = full_kmer_table)
    },mc.cores = ncores)
    # }
    kmer.ls <- lapply(kmer.full.ls,function(d) d$mat)
    full.ls <- lapply(kmer.full.ls,function(d) d$full_kmer_table)
    names(kmer.ls) <- names(full.ls) <- names(dat.ls)

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
    if(verbose){
        message('Perform differential analysis...')
    }
    res <- diff_limma(feature_list = kmer.ls,
                      feature_type = 'kmer',
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

    append(res,list(full_kmer_table = full.ls))

}
