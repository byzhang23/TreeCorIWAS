#' get unique kmer per sequence
#' @param .data immunarch data
#' @param .k specify kmer length
#' @param .col specify amino acid or nucleotide
#' @param .coding coding region
#' @param ncores number of cores
#' @importFrom data.table data.table
#' @importFrom dplyr as_tibble
#' @import immunarch parallel
#' @export

getUnqKmers <- function(.data, .k, .col = c("aa", "nt"), .coding = TRUE, ncores = parallel::detectCores()) {
    seq_col <- paste0('CDR3.',.col[1])# immunarch::switch_type(.col[1])

    if (.coding) {
        .data <- coding(.data)
    }

    if (class(.data) == 'list') {

        res <- do.call(rbind,
                       mclapply(.data,function(s){
                           tbl <- data.table(split_to_unq_kmers(collect(select(s, seq_col), n = Inf)[[1]], .k = .k))
                           if(nrow(tbl) > 0){
                               res <- suppressMessages(merge(tbl,
                                                             s %>% select(Clones,sample,Barcode) %>% mutate(id = 1:nrow(s)) %>% data.table, by = 'id', all = T) %>% as_tibble %>%
                                                           group_by(Kmer,sample) %>% summarise(Count = sum(Clones),
                                                                                               Barcode = paste(Barcode,collapse = ';')))
                           }else{
                               res <- NULL
                           }
                           res
                       },mc.cores = ncores))


        mat <- res %>% select(-Barcode) %>% dplyr::filter(!is.na(Kmer)) %>% spread(sample,Count)

    } else {
        tbl <- data.table(split_to_unq_kmers(collect(select(s, seq_col), n = Inf)[[1]], .k = .k))
        if(nrow(tbl) > 0){
            res <- merge(data.table(split_to_unq_kmers(collect(select(.data, seq_col), n = Inf)[[1]], .k = .k)),
                         s %>% select(Clones,sample,Barcode) %>% mutate(id = 1:nrow(s)) %>% data.table, by = 'id', all = T) %>% as_tibble %>%
                group_by(Kmer,sample) %>% summarise(Count = sum(Clones),
                                                    Barcode = paste(Barcode,collapse = ';'))
            mat <- res %>% select(-Barcode) %>% dplyr::filter(!is.na(Kmer)) %>% spread(sample,Count)
        }else{
            mat <- res <- NULL
        }

    }

    list(mat = as_tibble(as.data.frame(mat)),
         full_kmer = as_tibble(as.data.frame(res)))
}


#' Count unique kmers
#' @param .data list of sequences
#' @param .k Integer. Size of kmers.
#' @param ncores number of cores
#' @importFrom dplyr as_tibble
#' @importFrom stringr str_sub
#' @import immunarch parallel
#' @export

split_to_unq_kmers <- function(.data, .k, ncores = parallel::detectCores()) {
    max_len <- max(nchar(.data))
    str.ls <- strsplit(.data, ";", fixed = TRUE)

    # deal with multiple chains
    tab <- do.call(rbind,mclapply(1:length(str.ls),
                                  function(id){
                                      x <- str.ls[[id]]
                                      sub <- do.call(rbind,lapply(1:(max_len - .k + 1), function(i) substr(x, i, i + .k - 1)))
                                      x1 <- unique(sub[nchar(sub[,1]) >= .k,1])
                                      if(ncol(sub)>1){
                                          x2 <- unique(sub[nchar(sub[,2]) >= .k,2])
                                          if(length(x1)>0 & length(x2)>0){
                                              m <- expand.grid(x1 = x1, x2 = x2)
                                          }else{
                                              m <- expand.grid(x1,x2)
                                          }
                                      }else{
                                          # check one chain
                                          m <- unique(sub[nchar(sub[,1]) >= .k,1,drop = F])
                                      }
                                      if(nrow(m) >0){
                                          return(as_tibble(as.data.frame(data.frame(Kmer = sapply(1:nrow(m),function(i) paste(unlist(m[i,]),collapse = ';')), id = id), stringsAsFactors = FALSE)))
                                      }else{
                                          return(NULL)
                                      }

                                  },mc.cores = ncores))


    # tab <- do.call(rbind,lapply(1:(max_len - .k + 1), function(i) substr(.data, i, i + .k - 1)))
    # tab <- table(unlist(sapply(1:ncol(tab),function(x) unique(tab[,x]))))
    # tab <- tab[nchar(names(tab)) == .k]
    # if(length(tab) > 1){
    #     tab <- as_tibble(as.data.frame(tab, stringsAsFactors = FALSE))
    # }else{
    #     tab <- as_tibble(as.data.frame(data.frame(Kmer = names(tab),Count = tab), stringsAsFactors = FALSE))
    # }
    # names(tab) <- c("Kmer", "Count")
    tab
}
