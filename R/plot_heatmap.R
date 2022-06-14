#' plot_heatmap
#'
#' @export
#' @import pheatmap dplyr RColorBrewer

plot_heatmap <- function(mat, sample_meta, scale = 'row',sort_bycol = NA,
                         cluster_rows = F,cluster_cols = F,show_colnames = F, show_rownames = T,
                         annotation_col = NA, annotation_colors = NA, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),breaks = NA,fontsize = 12){

    # annotation
    if(!is.na(annotation_col[1])){
        annot <- suppressMessages(inner_join(data.frame(sample = colnames(mat)),sample_meta[,c('sample',annotation_col)]) %>% arrange(across(annotation_col)))
        rownames(annot) <- annot$sample
        annot <- annot %>% select(-sample)

        if(!is.na(sort_bycol[1])){
            annot <- annot %>% arrange(across(sort_bycol))
            mat <- mat[,rownames(annot)]
        }

    }else{
        annot <- NA
    }

    if(!is.na(breaks[1])){
        if(scale=='none'){
            if(min(breaks) > min(mat)) breaks <- c(min(mat),breaks)
            if(max(breaks) < max(mat)) breaks <- c(breaks,max(mat))
        }else if(scale=='col'){
            mat_scale <- scale(mat)
            if(min(breaks) > min(mat_scale)) breaks <- c(min(mat_scale),breaks)
            if(max(breaks) < max(mat_scale)) breaks <- c(breaks,max(mat_scale))

        }else if(scale=='row'){
            mat_scale <- t(scale(t(mat)))
            if(min(breaks) > min(mat_scale)) breaks <- c(min(mat_scale),breaks)
            if(max(breaks) < max(mat_scale)) breaks <- c(breaks,max(mat_scale))

        }
    }

    nbreaks <- length(breaks)
    ncolor <- length(color)
    if(nbreaks != (ncolor + 1)){
        warning(paste(paste0('Length of `breaks` vector: ', nbreaks),
                      paste0('Length of `color` vector: ', ncolor),
                      paste0('length(breaks) ', ifelse(nbreaks<ncolor+1,'<','>'),' (length(color) + 1)'),
                      sep = '\n'))
    }
    # heatmap
    pheatmap(mat,
             color = color,
             breaks = breaks,
             scale = scale,
             show_colnames = show_colnames,
             show_rownames = show_rownames,
             cluster_rows = cluster_rows,
             cluster_cols = cluster_cols,
             annotation_col = annot,
             annotation_colors = annotation_colors,
             fontsize = fontsize)
}
