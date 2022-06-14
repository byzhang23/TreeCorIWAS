#' Get Contact Region (remove residues)
#'@param .data data frame
#'@param .start starting position
#'@param .end ending position
#'@export
#'@import dplyr
#'@importFrom stringr str_sub

getContactRegion <- function (.data,.start = 4,.end = -4){
    if (sum(class(.data) %in% "list")>0) {
        lapply(.data, getContactRegion)
    }
    else {
        dt_flag <- FALSE
        if (!is.null(.data)){
            if (sum(class(.data) %in% "data.table")>0) {
                dt_flag <- TRUE
                .data <- .data %>% lazy_dt()
            }
            d <- collect(.data, n = Inf)
            d <- d %>% filter(!is.na(CDR3.aa))

            # deal with multiple chains
            d$CDR3.aa <- vapply(strsplit(d$CDR3.aa, ";", fixed = TRUE),
                                function(x) paste(stringr::str_sub(x,start = .start,end = .end), collapse = ";"),
                                character(1L))
            # d$CDR3.aa <- stringr::str_sub(d$CDR3.aa,start = .start,end = .end)
            if (dt_flag) {
                d <- data.table(d)
            }
            d <- d %>% filter(CDR3.aa != '')
        }else{
            d <- NULL
        }
        d
    }
}
