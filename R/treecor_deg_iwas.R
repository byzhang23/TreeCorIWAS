treecor_deg_iwas <- function(expr, hierarchy_list, cell_meta, sample_meta,
                             pool_pseudobulk_columns = 'sample',
                             response_variable, separate = T,differential_method = 'limma',
                             weight = NULL, formula = NULL, # if want to test for interaction, inlucde in formula
                             coef = 2, fdr_cutoff = 0.05,
                             filter_prop = 0.1,
                             pseudobulk_list = NULL,
                             ncores = parallel::detectCores(),
                             save_as_csv = T, verbose = T){
    ## construct pseudobulk
    if(length(pool_pseudobulk_columns) > 1){
        # add new column
        new_col <- paste0(pool_pseudobulk_columns,collapse = '.')
        message(paste0('Multiple columns are specified. Collapse these columns into a new column: ',new_col))
        sample_meta[,new_col] <- apply(sample_meta[,columns],1,function(x) paste0(x,collapse = '.'))
    }else{
        new_col <- pool_pseudobulk_columns
    }

    ## Extract info
    label_info <- hierarchy_list$layout
    if(sum(duplicated(label_info$label))) stop('Duplicated cell type name.')
    leaves_info <- hierarchy_list$leaves_info

    if(!is.null(formula)){
        formula = gsub('\\~|\\ ','',formula)
    }
    if(verbose){
        message('Construct sample-level pseudobulk...')
    }
    unq_id <- unique(leaves_info$id)

    if(is.null(pseudobulk_list)){
        pb.ls <- mclapply(unq_id,function(tid){
            if(verbose) message(paste0('node ',tid))
            node_info <- leaves_info %>% filter(id==tid)
            sub_meta <- cell_meta %>% filter(celltype %in% node_info$children)
            pb <- getPseudobulk(expr = expr,
                                cell_meta = sub_meta,
                                sample_meta = sample_meta,
                                columns = new_col,
                                filter_prop = filter_prop)
        },mc.cores = ncores)
        names(pb.ls) <- unq_id
    }else{
        # input pseudobulk list (user-specify)
        message('Input pseudobulk list:')
        pb_names <- intersect(as.character(unq_id),names(pseudobulk_list))
        if(length(pb_names) != length(unq_id)){
            ## extract label name
            ol <- names(pseudobulk_list)
            nl <- suppressMessages(inner_join(data.frame(label = ol),
                                              label_info)) %>% select(id) %>% unlist
            names(pseudobulk_list) <- nl
            pb_names <- intersect(as.character(unq_id),names(pseudobulk_list))
            if(length(pb_names) != length(unq_id)){
                stop('Invalid pseudobulk list, please check if the following nodes exist: \n', paste(setdiff(as.character(unq_id),pb_names),collapse = ','))
            }
        }
        pb.ls <- pseudobulk_list
    }

    ## DEG
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
                      feature_type = 'deg',
                      hierarchy_list,
                      sample_meta,
                      response_variable,
                      separate,
                      differential_method,
                      formula,
                      coef,
                      fdr_cutoff,
                      new_col,
                      ncores,
                      save_as_csv,
                      verbose)

}
