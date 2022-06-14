#' Identify exapnded clones
#' @param data  a list of immune profiling info for each sample. in `immunarch` format. It corresponds to extract [['data']] slot in `immunarch`.
#' @param match   match by 'clonotype' or 'clonotype+vdj'
#' @param type  use 'count' or 'frequency' to define expanded clones
#' @param threshold  specify a minimum threshold for identifying expanded clones for each sample. Default threshold is 3. Clones shared among at least `threshold` number of cells would be identified as expanded clones.
#' @param aggregate  specify either aggregate at sample-level ('sample') or combine all samples together ('all')
#' @param ncores    number of cores
#' @export
#' @import immunarch dplyr parallel
#' @author Boyang Zhang <bzhang34@jhu.edu>


getExpandedClone <- function(data, match = 'clonotype', type = 'count', threshold = 3, aggregate = 'sample',ncores = parallel::detectCores()){

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
    # raw data
    raw <- sapply(data,function(x){
        if(!is.null(x)){
            r <- nrow(x)
        }else{
            r <- NA
        }
        r
    }) %>% unlist

    # aggregate at which level
    if(aggregate == 'all'){

        agg <- do.call(rbind,data)
        if(match == 'clonotype+vdj'){
            agg <- agg %>% mutate(key = paste0(CDR3.aa,'_',V.name,'+',D.name,'+', J.name))
        }else if (match == 'clonotype'){
            agg <- agg %>% mutate(key = CDR3.aa)
        }
        comb <- agg %>% group_by(key) %>% summarise(Clones = sum(Clones))
        comb$Freq <- comb$Clones/sum(comb$Clones)
        if(type == 'count'){
            expanded.clone <- comb %>% filter(Clones >= threshold)
        }else if(type == 'frequency'){
            expanded.clone <- comb %>% filter(Freq >= threshold)
        }
        expanded.CDR3.key <- expanded.clone$key

        # by each sample
        expand.ls <- mclapply(1:length(data),function(i){
            dt <- data[[i]]
            if(!is.null(dt)){
                if(match == 'clonotype'){
                    expanded.clone <- dt %>% mutate(key = CDR3.aa) %>%
                        group_by(key) %>% summarise(Clones = sum(Clones)) %>%
                        filter(key %in% expanded.CDR3.key)
                }else if(match == 'clonotype+vdj'){
                    expanded.clone <- dt %>%
                        mutate(key = paste0(CDR3.aa,'_',V.name,'+',D.name,'+', J.name)) %>%
                        group_by(key) %>% summarise(Clones = sum(Clones)) %>%
                        filter(key %in% expanded.CDR3.key)
                }

                res <- expanded.clone %>% dplyr::select(key, Clones) %>% dplyr::rename(Key = key)
            }else{
                res <- NULL
            }
            res
        },mc.cores = ncores)
        names(expand.ls) <- names(data)


    }else{
        expand.ls <- mclapply(1:length(data),function(i){

            dt <- data[[i]]
            if(!is.null(dt)){
                if(match == 'clonotype'){
                    comb <- dt %>% mutate(key = CDR3.aa) %>% group_by(key) %>% summarise(Clones = sum(Clones))
                }else if(match == 'clonotype+vdj'){
                    comb <- dt %>% mutate(key = paste0(CDR3.aa,'_',V.name,'+',D.name,'+', J.name)) %>% group_by(key) %>% summarise(Clones = sum(Clones))
                }
                comb$Freq <- comb$Clones/sum(comb$Clones)

                if(type == 'count'){
                    expanded.clone <- comb %>% filter(Clones >= threshold)
                }else if(type == 'frequency'){
                    expanded.clone <- comb %>% filter(Freq >= threshold)
                }
                res <- expanded.clone %>% dplyr::select(key, Clones) %>% dplyr::rename(Key = key)
            }else{
                res <- NULL
            }
            res
        },mc.cores = ncores)
        names(expand.ls) <- names(data)
    }

    # proc
    proc <- sapply(expand.ls,function(x){
        if(!is.null(x)){
            r <- nrow(x)
        }else{
            r <- NA
        }
        r
    }) %>% unlist

    summary <- suppressMessages(inner_join(data.frame(sample = names(raw),total_clones = raw),
                                           data.frame(sample = names(proc),expanded_clones = proc))) %>%
        mutate(expanded_prop = expanded_clones/total_clones)

    summary[is.na(summary)] <- 0

    # return
    res <- list(summary = summary,
                expanded_clone_list = expand.ls)
}
