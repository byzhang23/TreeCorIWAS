#' Perform GO analysis
#' @param hierarchy_list    hierarchy list
#' @param deg_list          differential expressed genes, can be obtained by running `treecor_deg_iwas` or `treecor_deg`
#' @param background_gene   specify a set of background genes
#' @param ontology          specify ontology to be used in GO analysis. Default is 'BP' (biological process).
#' @param annotation_function specify an annotation function. Default is 'annFUN.org'.
#' @param mapping           specify mapping database, default is 'org.Hs.eg.db'.
#' @param algorithm         algorithm to be used in runTest
#' @param statistic         statistic to be used in runTest
#' @param fdr_cutoff        FDR cutoff (default = 0.05)
#' @param fold_change_cutoff  Fold change cutoff (default = 1.5)
#' @param min_DEGs          proceed with GO analysis with at least `min_DEGs` (default = 50)
#' @param min_annotated_genes  retain GO terms with at least `min_annotated_genes` (default = 10)
#' @param ncores            multiple cores
#' @param verbose           show progress
#' @export
#' @import topGO dplyr org.Hs.eg.db
#' @importFrom tidyr spread

treecor_GO <- function(hierarchy_list, deg_list, background_gene, ontology_name = 'BP', annotation_function = 'annFUN.org', mapping = 'org.Hs.eg.db',  algorithm = 'classic', statistic = 'fisher',
                       min_DEGs = 50, min_annotated_genes = 10, fdr_cutoff = 0.05, fold_change_cutoff = 1.5,
                       ncores = parallel::detectCores(), verbose = T){

    label_info <- hierarchy_list$layout
    if(verbose){
        message('Perform GO enrichment analysis...')
        message(paste0('Using ontology = ',ontology_name,'; annotation_function = ',annotation_function))
        message(paste0('Using algorithm = ',algorithm,'; statistic = ',statistic))
    }
    res.ls <- mclapply(1:length(deg_list),function(i){

        if(verbose) message(paste0('node ',i))
        dt <- deg_list[[i]]
        if(is.null(dt)){
            sigres <- NULL
        }else{

            gl <- rownames(dt)
            bg <- background_gene

            if(length(gl) >= min_DEGs){
                if(length(setdiff(gl,bg)) > 0) stop(paste0('Error. There are ',length(setdiff(gl,bg)), ' genes not in the background gene set. Please check it again!'))

                geneList <- factor(as.integer(bg %in% gl))
                names(geneList) <- bg

                suppressMessages({
                    GOdata <- new("topGOdata",
                                  ontology = ontology_name,
                                  allGenes = geneList,
                                  geneSel=function(a) {a},
                                  annot = get(annotation_function),
                                  mapping = mapping,
                                  ID = "Symbol")
                    resultFisher <- runTest(GOdata, algorithm = algorithm, statistic = statistic)

                    sigres <- GenTable(GOdata, classic_fisher = resultFisher, topNodes = length(resultFisher@score),orderBy="classic_fisher",numChar=1000)
                })

                sigres$classic_fisher[sigres$classic_fisher=="< 1e-30"] <- 1e-30
                sigres$classic_fisher <- as.numeric(sigres$classic_fisher)
                sigres <- sigres[sigres$Annotated >= min_annotated_genes,]
                sigres$FDR <- p.adjust(sigres$classic_fisher,method="fdr")
                colnames(sigres)[which(colnames(sigres) == 'classic_fisher')] <- paste0(algorithm,'_',statistic)
                fc <- (sigres[,"Significant"]/sum(GOdata@allScores[GOdata@feasible]==1))/(sigres[,"Annotated"]/sum(GOdata@feasible))
                sigres <- data.frame(sigres,fold_change = fc) %>% mutate(FDR = as.numeric(FDR)) %>% arrange(FDR,-fold_change)

                if(!is.null(sigres)){
                    if(nrow(sigres) > 0){
                        dt <- sigres %>% dplyr::filter(FDR <= fdr_cutoff & fold_change >= fold_change_cutoff)
                        if(nrow(dt) > 0){

                            fisher.go <- dt$GO.ID
                            fisher.ann.genes <- lapply(genesInTerm(GOdata, whichGO=fisher.go), function(x) paste(intersect(x,gl),collapse = ';'))
                            fisher.ann.dt <- data.frame(`GO.ID` = names(fisher.ann.genes),
                                                        Significant_Annotated_Genes = unlist(fisher.ann.genes))
                            sigres <- suppressMessages(left_join(sigres,fisher.ann.dt))
                        }
                    }
                }

            }else{
                sigres <- NULL
            }
        }
        sigres
    },mc.cores = ncores)

    sign.ls <- lapply(res.ls,function(dt){
        if(!is.null(dt)){
            dt <- dt %>% dplyr::filter(FDR <= fdr_cutoff & fold_change >= fold_change_cutoff)
            if(nrow(dt) == 0){
                dt <- NULL
            }
        }else{
            dt <- NULL
        }
        dt
    })

    names(res.ls) <- names(sign.ls) <- names(deg_list)

    ## summary as a wide table
    if(verbose) message('Reshape the significant GO terms into a table...')
    sign.union <- sapply(sign.ls[which(sapply(sign.ls,function(x) !is.null(x)))],function(x) paste0(x$GO.ID,'\n',x$Term)) %>% unlist %>% unique
    sub.long <- do.call(rbind,lapply(1:length(res.ls),function(i){
        if(verbose) message(paste0('node ',i))
        dt <- res.ls[[i]]
        if(!is.null(dt)){
            dt <- dt %>% select(GO.ID, Term, fold_change, FDR) %>% mutate(GO.term = paste0(GO.ID,'\n',Term)) %>%
                dplyr::filter(GO.term %in% sign.union) %>% mutate(label = names(res.ls)[i]) %>%
                select(c(GO.term, fold_change, FDR, label))
            if(nrow(dt) == 0){
                dt <- data.frame(GO.term = sign.union, fold_change = NA, FDR = NA, label = names(res.ls)[i])
            }
        }else{
            dt <- data.frame(GO.term = sign.union, fold_change = NA, FDR = NA, label = names(res.ls)[i])
        }
        dt
    }))

    sub.wide.fdr <- sub.long %>% select(label, GO.term, FDR) %>% tidyr::spread(GO.term,FDR)
    colnames(sub.wide.fdr)[-1] <- paste0(colnames(sub.wide.fdr)[-1],'.fdr')

    sub.wide.fc <- sub.long %>% select(label, GO.term, fold_change) %>% tidyr::spread(GO.term,fold_change)
    colnames(sub.wide.fc)[-1] <- paste0(colnames(sub.wide.fc)[-1],'.fold_change')

    sub.wide <- suppressMessages(inner_join(label_info,inner_join(sub.wide.fc,sub.wide.fdr)))


    list(summary = sub.wide,
         signGO = sign.ls,
         allGO = res.ls)

}
