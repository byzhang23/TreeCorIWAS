#' Get gene usage
#' @param data  a list of immune profiling info for each sample. in `immunarch` format. It corresponds to extract [['data']] slot in `immunarch`.
#' @param gene  specify one or multiple queried tcr or bcr genes (use '+' to separate).
#' @param norm  an indicator to specify normalization or not. Default is TRUE (normalized to frequency).
#' @param transform specify either 'clr' (centered log ratio) transformation or 'none' (default).
#' @param ncores number of cores
#' @export
#' @import immunarch dplyr parallel
#' @importFrom compositions clr
#' @author Boyang Zhang <bzhang34@jhu.edu>


# unexplored parameters: ambig, type, quant (get_genes does not contain all valid VDJ genes)
# library(dplyr)
# library(immunarch)
# library(compositions)

getGeneUsage <- function(data,gene,norm = T,transform = 'none',ncores = parallel::detectCores(), filter_prop = 0.1){

    gene <- gsub(' ','',tolower(gene[1]))

    # VDJ gene usage
    multi_genes <- strsplit(gene,'\\+') %>% unlist
    if(sum(multi_genes %in% c('v','d','j')) != length(multi_genes)){
        stop('Incorrect gene usage, please use either `v` or `d` or `j`.')
    }

    # 1 sample
    if(class(data) != 'list'){
        data <- list(data)
        if(sum(colnames(data)=='sample')>0){
            names(data) <- data[[1]]$sample %>% unique
        }else{
            warning('Arbitrary assign a sample name: s1')
            names(data) <- 's1'
        }

    }

    # # check input and data
    # imrc_type_data <- sub('\\V.*|\\v.*','',data[[1]]$V.name[1]) %>% tolower
    # imrc_type_input <- substr(multi_genes[1],1,3)
    # if(imrc_type_data != imrc_type_input) stop(paste0('Inconsistent immune receptor type: data is ',imrc_type_data,', but input is ',imrc_type_input))

    # gene_db <- expand.grid(lapply(paste0('hs.',multi_genes),function(g) immunarch::get_genes(g)))
    # gene_db$Names <- do.call(paste, c(gene_db, sep = ';'))
    # column_nm <- paste0(substr(multi_genes,4,4) %>% toupper(),'.name')

    column_nm <- paste0(toupper(multi_genes),'.name')
    gene_db <- do.call(rbind,lapply(data,function(x) x[,column_nm,drop=F])) %>% unique
    if(length(column_nm)>1){
        gene_db$Names <- do.call(paste, c(gene_db, sep = '+'))
    }else{
        gene_db$Names <- gene_db[,column_nm]
    }

    if(sum(column_nm %in% 'D.name')>0 & length(setdiff(unique(data[[1]]$D.name),c('None','none',NA,'NA')))==0) warning('`D.name` contains no valid D genes. Suggest to exclude D gene usage from this analysis.')

    summary.ls <- mclapply(1:length(data),function(i){
        message(paste0('sample ', i))
        dt <- data[[i]]
        if(!is.null(dt)){
            dt$Names <- do.call(paste, c(dt[column_nm], sep = '+'))
            dt$Names[dt$Names == paste(rep('NA',length(column_nm)),collapse = '+')] <- NA
            bysamp <- suppressMessages(dt %>% group_by(Names) %>% summarise(count = sum(Clones)))
            res <- suppressMessages(left_join(gene_db %>% select(Names),bysamp))
            colnames(res)[2] <- names(data)[i]
        }else{
            res <- NULL
        }
        res
    },mc.cores = ncores)
    summary.ls <- summary.ls[which(sapply(summary.ls,function(x) !is.null(x)))]
    summary <- suppressMessages(Reduce(inner_join,summary.ls)) %>% data.frame(check.names = F) %>% filter(!is.na(Names))


    summary[is.na(summary)] <- 0
    dat <- summary[,-1] %>% as.matrix
    rownames(dat) <- summary$Names

    dat <- dat[rowSums(dat)>0,]
    if(norm){
        dat <- t(t(dat)/colSums(dat))
        if(transform == 'clr'){
            t.dat <- t(dat) # row: sample, col: genes
            dat <- t(sapply(colnames(t.dat),function(c) unclass(compositions::clr(t.dat[,c])),USE.NAMES = TRUE))
        }
    }
    # filter prop
    dat <- dat[rowMeans(abs(dat) > 0) >= filter_prop,,drop = F]
    if(nrow(dat) == 0) dat <- NULL

    dat
}
