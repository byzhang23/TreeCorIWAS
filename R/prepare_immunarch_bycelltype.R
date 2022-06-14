#' Prepare immunarch type of data across cell types
#'@param immunarch_data_tra immunarch format for tra
#'@param immunarch_data_trb immunarch format for trb
#'@param input_string string to define tree structure
#'@param cell_meta cell metadata
#'@param filter_chain specify chain to be filtered, can be one of 'tra' or 'trb' or 'tra+trb'.
#'@param min_cells minimum number of cells per cell type
#'@param ncores number of cores
#'@param verbose Default is T
#'@export
#'@import dplyr tidyr TreeCorTreat parallel

prepare_immunarch_by_celltype <- function(immunarch_data_tra = NULL,immunarch_data_trb = NULL,
                                          input_string = NULL,cell_meta,
                                          filter_chain = 'tra', min_cells = 1000, ncores = parallel::detectCores(), verbose = T){

    if(is.null(immunarch_data_tra) & is.null(immunarch_data_trb)) stop('Please specify either TRA or TRB chain.')
    if(is.null(input_string)) stop('Please input a valid input_string; Prepare string in following format: @All(cellltype_1,celltype_2)')

    if (Sys.info()[['sysname']]=='Windows') {
        message('Parallel is disabled for Windows. Running with one core')
        ncores <- 1
    }

    convert_to_long <- function(data){
        long <- do.call(rbind, parallel::mclapply(1:length(data), function(i){

            sample <- names(data)[i]
            if(verbose) message(sample)
            dt <- data[[i]]
            if(sum(colnames(dt) %in% 'sample') == 0) dt$sample <- sample
            new <- dt %>% tidyr::as_tibble() %>% tidyr::separate_rows(Barcode,sep = ';',convert = T)
            new
        },mc.cores = ncores))
    }

    # TRA
    if(!is.null(immunarch_data_tra)){
        if(verbose) message('Process TRA chain...')

        tra.long <- convert_to_long(immunarch_data_tra)
        colnames(tra.long) <- gsub('Barcode','barcode',colnames(tra.long))
        tra <- suppressMessages(inner_join(tra.long,cell_meta))

    }else{
        tra <- summary.tra <- NULL
    }

    # TRB
    if(!is.null(immunarch_data_trb)){
        if(verbose) message('Process TRB chain...')

        trb.long <- convert_to_long(immunarch_data_trb)
        colnames(trb.long) <- gsub('Barcode','barcode',colnames(trb.long))
        trb <- suppressMessages(inner_join(trb.long,cell_meta))

    }else{
        trb <- summary.trb <- NULL
    }

    # TRA and TRB
    if(!is.null(immunarch_data_tra) & !is.null(immunarch_data_trb)){

        if(verbose) message('Process TRA and TRB chains...')

        a1 <- tra; b1 <- trb
        unchanged <- c('barcode','sample','celltype')
        colnames(a1)[-which(colnames(a1) %in% unchanged)] <- paste0('tra.',colnames(a1)[-which(colnames(a1) %in% unchanged)])
        colnames(b1)[-which(colnames(a1) %in% unchanged)] <- paste0('trb.',colnames(b1)[-which(colnames(a1) %in% unchanged)])
        ab1 <- suppressMessages(inner_join(a1,b1))
        common.col <- intersect(gsub('tra.','',colnames(a1)[-which(colnames(a1) %in% unchanged)]),
                                gsub('trb.','',colnames(b1)[-which(colnames(b1) %in% unchanged)]))
        for(col in common.col){
            mod.cols <- c(paste0('tra.',col),paste0('trb.',col))
            ab1[,col] <- do.call(paste, c(ab1[mod.cols], sep = ";"))
        }
        trab <- ab1[,c(unchanged,common.col)]


    }else{
        trab <- summary.trab <- NULL
    }

    if(verbose) message(paste0('Remove cell types < ', min_cells,' cells (', filter_chain,' chain):'))
    if(filter_chain == 'tra' & !is.null(tra)){
        remove.ct <- tra %>% group_by(celltype) %>% summarise(count = n()) %>% filter(count < min_cells) %>% select(celltype) %>% unlist
    }else if(filter.chain == 'trb' & !is.null(trb)){
        remove.ct <- trb %>% group_by(celltype) %>% summarise(count = n()) %>% filter(count < min_cells) %>% select(celltype) %>% unlist
    }else if(filter.chain == 'tra+trb' & !is.null(tra) & !is.null(trb)){
        remove.ct <- trab %>% group_by(celltype) %>% summarise(count = n()) %>% filter(count < min_cells) %>% select(celltype) %>% unlist
    }

    if(length(remove.ct) > 0){
        remove_string <- paste(remove.ct, collapse = '|')
        if(verbose) message(remove_string)

        input_string <- gsub('\\,)',')',
                             gsub('\\(,','(',
                                  gsub(',{2,}',',',
                                       gsub(remove_string,'',input_string))))
        if(verbose) cat(paste0('New input string after filtering:\n', input_string))

    }

    hierarchy_list <- TreeCorTreat::extract_hrchy_string(input_string,special_character = '@', plot = F)
    layout <- hierarchy_list$layout
    leaf_info <- hierarchy_list$leaves_info

    ## filter tra and trb
    if(!is.null(tra)){
        tra <- tra %>% filter(celltype %in% layout$label)
        ## summary
        summary.tra <- suppressMessages(inner_join(cell_meta %>% group_by(celltype) %>% summarise(n_expr_cells = n_distinct(barcode)),
                                                   tra %>% group_by(celltype) %>% summarise(n_tra_cells = n_distinct(barcode))))

    }
    if(!is.null(trb)){
        trb <- trb %>% filter(celltype %in% layout$label)
        ## summary
        summary.trb <- suppressMessages(inner_join(cell_meta %>% group_by(celltype) %>% summarise(n_expr_cells = n_distinct(barcode)),
                                                   trb %>% group_by(celltype) %>% summarise(n_trb_cells = n_distinct(barcode))))

    }
    if(!is.null(trab)){
        trab <- trab %>% filter(celltype %in% layout$label)
        ## summary
        summary.trab <- suppressMessages(inner_join(cell_meta %>% group_by(celltype) %>% summarise(n_expr_cells = n_distinct(barcode)),
                                                    trab %>% group_by(celltype) %>% summarise(n_trab_cells = n_distinct(barcode))))
    }
    summary <- list(summary.tra,summary.trb,summary.trab)
    summary <- Reduce(inner_join,summary[!sapply(summary,is.null)])

    ## split by cell type
    # leaf
    if(!is.null(tra)) tra.leaf <- split(tra,factor(tra$celltype))
    if(!is.null(trb)) trb.leaf <- split(trb,factor(trb$celltype))
    if(!is.null(trab)) trab.leaf <- split(trab,factor(trab$celltype))

    # nonleaf
    tra.nonleaf <- trb.nonleaf <- trab.nonleaf <- list()
    nonleaf.label <- layout$label[!layout$leaf]
    for(i in 1:length(nonleaf.label)){

        children <- leaf_info %>% filter(label == nonleaf.label[i]) %>% select(children) %>% unlist
        if(!is.null(tra)){
            tra.nonleaf[[i]] <- do.call(rbind,lapply(tra.leaf[children],function(x){
                x$celltype <- nonleaf.label[i]
                x
            }))
            rownames(tra.nonleaf[[i]]) <- NULL
        }

        if(!is.null(trb)){
            trb.nonleaf[[i]] <- do.call(rbind,lapply(trb.leaf[children],function(x){
                x$celltype <- nonleaf.label[i]
                x
            }))
            rownames(trb.nonleaf[[i]]) <- NULL
        }

        if(!is.null(trab)){
            trab.nonleaf[[i]] <- do.call(rbind,lapply(trab.leaf[children],function(x){
                x$celltype <- nonleaf.label[i]
                x
            }))
            rownames(trab.nonleaf[[i]]) <- NULL
        }
    }
    if(!is.null(tra)) names(tra.nonleaf) <- nonleaf.label
    if(!is.null(trb)) names(trb.nonleaf) <- nonleaf.label
    if(!is.null(trab)) names(trab.nonleaf) <- nonleaf.label

    if(!is.null(tra)) tra.ls <- append(tra.nonleaf, tra.leaf)
    if(!is.null(trb)) trb.ls <- append(trb.nonleaf, trb.leaf)
    if(!is.null(trab)) trab.ls <- append(trab.nonleaf, trab.leaf)

    ## use immunarch format
    listBySample <- function(dat){
        pattern <- unique(dat$chain)
        dat.ls <- parallel::mclapply(unique(dat$sample),function(s){
            if(verbose) message(s)
            sub <- dat %>% filter(sample==s)
            ttl <- length(unique(sub$barcode))
            suppressMessages(sub %>%
                                 group_by(chain,sample,CDR3.nt,CDR3.aa,V.name,D.name,J.name) %>%
                                 summarise(Clones = n_distinct(barcode),
                                           Proportion = n_distinct(barcode)/ttl,
                                           Barcode = paste0(barcode,collapse = ';'))) %>% data.frame
        },mc.cores = ncores)
        names(dat.ls) <- unique(dat$sample)
        dat.ls
    }

    tra.proc <- trb.proc <- trab.proc <- NULL
    if(!is.null(tra)) tra.proc <- lapply(tra.ls,function(dt) listBySample(dt))
    if(!is.null(trb)) trb.proc <- lapply(trb.ls,function(dt) listBySample(dt))
    if(!is.null(trab)) trab.proc <- lapply(trab.ls,function(dt) listBySample(dt))

    comb <- list(tra = tra.proc,
                 trb = trb.proc,
                 'tra+trb' = trab.proc)

    res <- list(summary = summary,
                immunarch_by_celltype = comb)

}
