#' Sample-level differential analysis
#'
#' Provides feature sets and uses LIMMA to identify differentially features for each cell cluster.
#' @param feature_list    A list of sample-level feature matrices. Default is NULL.
#' @param feature_type      Name of input feature type, will be used as suffix of output file and column names.
#' @param hierarchy_list a hierarchy list by running \code{'extract_hrchy_string()'} or \code{'extract_hrchy_seurat()'} function.
#' @param sample_meta       Sample-level metadata, where each row is a sample. Must contain 'sample' column and additional variables to be used in the analysis, such as covariates or outcomes of interest.
#' @param response_variable A vector of response variables.
#' @param separate          A TRUE (default) or FALSE indicator, specifying how to evaluate multivariate outcomes.
#'                              \itemize{
#'                              \item TRUE: evaluate multivariate phenotype separately (it is equivalent to run this pipeline for each univariate phenotype).
#'                              \item FALSE: evaluate multivariate phenotype jointly.
#'                              }
#' @param differential_method   Specify which type of differential methods: limma (default) or gam
#' @param formula           An object of class 'formula': a symbolic description of adjustment formula (i.e. only includes covariates other than response variable(s))
#' @param coef              A column number or column name specifying which coefficient to be extracted (by default: 2).
#' @param fdr_cutoff        Cutoff value for FDR. Only genes with lower FDR are listed. Default is 0.05.
#' @param column_name       column name of feature matrix. default is 'sample.
#' @param ncores            Number of cores to be used. If ncores > 1, it will be implemented in a parallel mode.
#' @param save_as_csv       An indicator to save identified DEGs in csv files. DEGs for each cell cluster is saved as 'responsevariable_celltype_DEG.csv' and a summary file of combining DEGs from all cell clusters is saved as 'responsevariable_combinedDEG.csv'.
#' @param verbose           Show progress
#' @return A list of three elements:
#' \itemize{
#' \item diff.summary: A summary table of number of differential features for each tree node.
#' \item diff.ls: A comprehensive list of outcome(s)-associated feature set for each tree node. Use \code{`result$dge.ls$response_variable[[celltype]]`} to extract DEGs for a specific cell type
#' \item feature.ls: A list of input feature list
#' \item diff.summary.bygene: A summary table for each differentiated gene.
#'}
#' @export
#' @import dplyr parallel limma immunarch mgcv
#' @importFrom lmtest lrtest
#' @author Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji

diff_limma <- function(feature_list, feature_type,hierarchy_list, sample_meta, response_variable,
                       separate = T, differential_method = 'limma', formula = NULL, coef = 2, fdr_cutoff = 0.05, column_name = 'sample',
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

    if(!is.null(formula)){
        formula = gsub('\\~|\\ ','',formula)
    }

    ## =============================== ##
    ## 1. read in feature list
    ## =============================== ##
    pb.ls <- feature_list
    if(sum(names(pb.ls) == as.character(nrow(label_info))) == 0){
        pb.label <- names(pb.ls)
        pb.id <- suppressMessages(inner_join(data.frame(label = names(pb.ls)),
                                             label_info)) %>% select(id) %>% unlist
        names(pb.ls) <- pb.id
    }

    ## =============================== ##
    ## 2. Differential
    ## =============================== ##
    if(verbose) message(paste0('Conduct differential test using ', differential_method,'...'))
    resp.ls <- lapply(response_variable,function(r){

        if(verbose) message(paste0('## Response variable: ',r))

        dge.ls <- mclapply(unq_id,function(tid){

            if(verbose) message(paste0('node ',tid))
            pb <- pb.ls[[as.character(tid)]]
            if(!is.null(pb)){
                if(nrow(pb) > 0 & ncol(pb)>0){
                    coln <- data.frame(sample = colnames(pb))
                    colnames(coln) <- column_name

                    tmp_meta <- suppressMessages(inner_join(coln,sample_meta))

                    if(nrow(tmp_meta)>1){
                        # drop factor level
                        factor.id <- which(sapply(1:ncol(tmp_meta),function(x) is.factor(tmp_meta[,x])))
                        for(k in factor.id){
                            tmp_meta[,k] <- droplevels(tmp_meta[,k])
                        }

                        ## limma
                        if(differential_method == 'limma'){
                            # design matrix
                            design <- try(model.matrix(as.formula(paste(c(paste0('~',r),formula),collapse = '+')),data = tmp_meta))
                            if(verbose) message(paste0('Checking model: ',paste(c(paste0('~',r),formula),collapse = '+')))

                            # check collinearity/error in building design matrix
                            if(class(design)[1] == 'try-error'){ #|qr(design)$rank<ncol(design)
                                warning(paste0('Error occurs, proceed with dropping the covariates: `formula = ~',r,'`'))
                                design <- model.matrix(as.formula(paste0('~',r)),data = tmp_meta)
                                if(verbose & !(qr(design)$rank<ncol(design))) message(paste0('Final linear model fit: ~',r))
                            }
                            if(verbose){
                                message('Column names of design matrix:\n', paste(colnames(design),collapse  = ' '))
                                if(is.numeric(coef)){
                                    if(coef <= ncol(design)){
                                        message('Extract coefficient by column name: ', colnames(design)[coef])
                                    }else{
                                        stop(paste0('Incorrect column id specifying which coefficient to be extracted. Must be in the range of 1-',ncol(design)))
                                    }
                                }else{
                                    if(sum(colnames(design)==coef)>0){
                                        message('Extract coefficient by column name: ', coef)
                                    }else{
                                        stop(paste0('Incorrect column name specifying which coefficient to be extracted. Must be one of :',paste(colnames(design),collapse  = ' ')))
                                    }
                                }
                            }
                            # limma
                            # if(qr(design)$rank<ncol(design)){
                            #     if(verbose) warning('Check model again. Proceed with num_DEG = 0.')
                            #     diff <- NULL
                            # }else{
                            fit <- lmFit(pb,design)
                            eb <- try(eBayes(fit))
                            if(class(eb)[1] == 'try-error'){
                                diff <- NULL
                            }else{
                                diff <- topTable(eb,coef = coef,n=nrow(pb),sort.by = "P", p.value = fdr_cutoff)
                                if(nrow(diff)>0){
                                    diff <- diff[,c('logFC','t','P.Value','adj.P.Val')]
                                    if(is.numeric(coef)){
                                        colnames(diff) <- paste0(colnames(design)[coef],'.',c('logFC','t','p','fdr'))
                                    }else{
                                        colnames(diff) <- paste0(coef,'.',c('logFC','t','p','fdr'))
                                    }
                                }else{
                                    diff <- NULL
                                }
                            }
                        }else if(differential_method == 'gam'){

                            full.formula <- paste(c(paste0('~s(',r,')'),formula),collapse = '+')
                            reduced.formula <- paste(c(paste0('~1'),formula),collapse = '+')
                            if(verbose) message(paste0('Full model: ',full.formula,'\n',
                                                       'Reduced model: ', reduced.formula, '\n'))

                            diff <- do.call(rbind,mclapply(1:nrow(pb),function(i){
                                y <- pb[i,]
                                full <- mgcv::gam(formula = as.formula(paste0('y',full.formula)),data = tmp_meta)
                                reduced <- mgcv::gam(formula = as.formula(paste0('y',reduced.formula)),data = tmp_meta)
                                lr <- lmtest::lrtest(full,reduced)

                                data.frame(gene = rownames(pb)[i],
                                           chisq = lr$Chisq[2],
                                           p = lr$`Pr(>Chisq)`[2])

                            }, mc.cores = ncores))
                            diff$fdr <- p.adjust(diff$p, method = 'fdr')
                            rownames(diff) <- diff$gene
                            diff <- diff %>% select(-gene) %>% filter(fdr < fdr_cutoff)
                            colnames(diff) <- paste0(r,colnames(diff))
                            if(nrow(diff) == 0) diff <- NULL

                        }

                    }else{
                        diff <- NULL
                    }
                }else{
                    diff <- NULL
                }
            }else{
                diff <- NULL
            }
            diff
        },mc.cores = ncores)
        names(dge.ls) <- label_info$label[unq_id]
        dge.ls

    })
    names(resp.ls) <- response_variable
    names(pb.ls) <- suppressMessages(inner_join(data.frame(id = names(pb.ls) %>% as.numeric),
                                                label_info)) %>% select(label) %>% unlist

    ## =============================== ##
    ## 3. Final result
    ## =============================== ##
    if(save_as_csv){
        for(resp in response_variable){
            rdt <- resp.ls[[resp]]
            for(j in names(rdt)){
                rjdt <- rdt[[j]]
                if(!is.null(rjdt)){
                    write.csv(rjdt,paste0(resp,'_',j,'_diff_',feature_type,'.csv'),row.names = T)
                }
            }
            # all in 1 csv
            rdt.comb <- do.call(rbind,lapply(1:length(rdt),function(j){
                dt <- rdt[[j]]
                if(!is.null(dt)){
                    dt$Gene <- rownames(dt)
                    dt$celltype <- names(rdt)[j]
                }
                dt
            }))
            if(nrow(rdt.comb)>0) write.csv(rdt.comb,paste0(resp,'_combined_diff_',feature_type,'.csv'),row.names = T)
        }
    }

    if(verbose){
        message('Summarize final result...')
    }
    dge.dat <- suppressMessages(Reduce(inner_join,lapply(1:length(resp.ls),function(i){

        dt <- resp.ls[[i]]
        df <- sapply(dt,function(d) ifelse(is.null(d),0,nrow(d))) %>% unlist
        setNames(data.frame(label = names(df),df),
                 c('label',paste0(names(resp.ls)[i],'.num_diff_',feature_type)))


    })))
    res <- suppressMessages(inner_join(dge.dat,label_info))

    # add summarize gene across clusters
    res.wide <- lapply(resp.ls, function(w){
        wide <- do.call(rbind,lapply(names(w),function(i){
            x <- w[[i]]
            if(!is.null(x)){
                a <- data.frame(gene = rownames(x), label = i)
            }else{
                a <- NULL
            }
            a
        }))

        if(is.null(wide)){
            summary <- NULL
        }else{
            summary <- inner_join(wide %>% group_by(gene) %>% mutate(labels = paste0(label, collapse = ",")),
                                  wide %>% group_by(gene) %>% summarise(num_nodes = n_distinct(label))) %>% select(-label) %>% arrange(-num_nodes)
        }
        summary %>% unique
    })


    final.ls <- list(diff.summary = res,
                     diff.ls = resp.ls,
                     feature.ls = pb.ls,
                     diff.summary.bygene = res.wide)
}
