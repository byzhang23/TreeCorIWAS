#' Visualize TreeCorTreat via Concantenation (Version 2)
#'
#' Generate TreeCorTreat plot to summarize and visualize the result.
#' @param hierarchy_list        A hierarchy list by running \code{'extract_hrchy_text()'} or \code{'extract_hrchy_seurat()'} function. Contains 'edges', 'layout', 'immediate_children' and 'leaves_info' as 4 elements.
#' @param annotated_df          A data frame contains summary statistic (canonical correlation or number of differential expressed genes, etc) and/or adjusted p-value for each tree node. Can be obtained by running \code{'treecor_ctprop()$canonical_corr'} or \code{'treecor_expr()$canonical_corr'} or \code{'treecor_deg()$dge.summary'}.
#' @param response_variable     A vector of response variables.
#' @param separate              A TRUE (default) or FALSE indicator, specifying whether response variables are evaluated separately (TRUE) or jointly (FALSE). This should be the same 'separate' parameter specified in \code{'treecor_ctprop()$canonical_corr'} or \code{'treecor_expr()$canonical_corr'}.
#' @param num_cancor_components Number of canonical components to be extracted (by default: 1). Only works when separate = F.
#' @param color_variable        Set color aesthetics to a variable. This variable should be one of column names in \code{annotated_df}.
#' @param size_variable         Set size aesthetics to a variable. This variable should be one of column names in \code{annotated_df}.
#' @param alpha_variable        Set alpha (transparency) aesthetics to a variable. This variable should be one of column names in \code{annotated_df}.
#' @param advanced_list         A list for more advanced figure configuration. Default is NULL. Users can specify label color or manually specify your own set of mappings from levels in the data to aesthetic values. Names of this list must be one of 'label_info','size_info','color_info', 'alpha_info','palette','plot_layout'.
#'                              \itemize{
#'                              \item label_info: a data frame that contains two columns: 'label' and 'label.color'. The label name should be same as the ones provided in \code{annotated_df}, which will be used to merged with \code{annotated_df} data frame.
#'                              \item color_info: a data frame that contains two columns: 'breaks' and 'values'. This will be used to override default setting using \code{scale_color_manual()} function (currently only support discrete variables)
#'                              \item size_info: a data frame that contains two columns: 'breaks' and 'values'. This will be used to override default setting using \code{scale_size_manual()} function (currently only support discrete variables)
#'                              \item alpha_info: a data frame that contains two columns: 'breaks' and 'values'. This will be used to override default setting using \code{scale_alpha_manual()} function (currently only support discrete variables)
#'                              \item palette: a string or an index specify palette when plotting continuous variable. Default is 'Spectral'.
#'                              \item layout_widths: a vector of two numbers to specify plot layout (left:right). The relative widths of non-leaf (left) and leaf part (right) in the grid. Default is c(10,5).
#'                              }
#' @param annotate_number       A TRUE or FALSE (default) indicator, specifying whether to annotate numbers in the plot.
#' @param annotate_number_column Specify column names to annotate text
#' @param annotate_number_fontsize Specify font size of annotated numbers
#' @param annotate_number_color Specify color of annotated text if \code{annotate_number = T}. Default is 'black'.
#' @param x_axis_text_angle     Specify angle of x axis. Default is 0.
#' @param font_size             Font size
#' @param nonleaf_label_pos     A number from 0 to 1, controlling the relative positions of non-leaf labels (default: 0.6). Default tree height between adjacent hierarchy levels is 1.
#' @param edge_path_type        Specify an edge type, must be one of:
#'                              \itemize{
#'                              \item link (default): Use straight lines to connect start and end nodes by applying \code{geom_edge_link()} in \code{ggraph} package. Recommended only for pre-set approach (i.e. generate hierarchy_list using \code{extract_hrchy_string()}).
#'                              \item diagonal: Use diagonal bezier curves to connect start and end nodes by applying \code{geom_edge_diagonal()} in \code{ggraph} package. Recommended only for pre-set approach (i.e. generate hierarchy_list using \code{extract_hrchy_string()}).
#'                              \item elbow: Draw classical dendrograms by applying \code{geom_edge_elbow()} in \code{ggraph} package. Recommended only for data-driven approach (i.e. generate hierarchy_list using \code{extract_hrchy_seurat()}).
#'                              }
#' @param plot_type             Specify a type to display the leaf part, must be one of:
#'                              \itemize{
#'                              \item balloon (default): balloon plot, where size of a dot reflects magnitude of result.
#'                              \item bar: bar plot, where height of each bar reflects magnitude of result.
#'                              \item heatmap: heatmap, where color of each grid reflects magnitude of result.
#'                              }
#' @param line_color            Specify line color for hierarchical tree structure (default: black).
#' @param line_type             Specify line type for hierarchical tree structure, can be one of 'solid' (default), 'dashed', 'dotted', 'dotdash', 'longdash' and 'twodash'.
#' @param plot                  Print the plot
#' @return A ggplot object
#' @export
#' @import dplyr ggplot2 ggraph RColorBrewer patchwork
#' @importFrom tidyr gather
#' @importFrom igraph graph_from_data_frame
#' @importFrom ggrepel geom_label_repel
#' @author Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji
#' @examples
#' treecortreatplot(hierarchy_list, annotated_df, response_variable = 'severity', color_variable = 'adjp.sign', size_variable = 'cancor', alpha_variable = 'adjp.sign', plot = T)


treecortreatplot_v2 <- function(hierarchy_list, annotated_df, response_variable, separate = T, num_cancor_components = 1,color_variable = NULL, size_variable = NULL, alpha_variable = NULL, advanced_list = NULL, annotate_number = F, annotate_number_column = NULL, annotate_number_fontsize = 5,annotate_number_color = 'black', font_size = 10,
                                x_axis_text_angle = 0, nonleaf_label_pos = 0.1, edge_path_type = 'link', plot_type = 'balloon',line_color = 'black', line_type = 'solid',plot = T){

    #############################
    ## 1. Define outcome variable
    #############################

    if(!separate){
        outcome_variable <- paste(response_variable,collapse = '_')
        if(num_cancor_components>1){
            outcome_variable <- paste0(outcome_variable,'_Component',1:num_cancor_components)
        }

    }else{
        outcome_variable <- response_variable
    }
    #############################
    ## 2. Combine annotated info with original skeleton
    #############################
    edges_dt <- hierarchy_list$edges %>% dplyr::rename(from = parent, to = children)
    mygraph <- graph_from_data_frame(edges_dt)
    layout_all <- create_layout(mygraph,layout = 'dendrogram')

    #############################
    ## 3. Basic info for ggplot2: size/color/alpha range/categories
    #############################
    # Fix colors for categorical variables
    gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
    }

    ## FORM A FUNCTION
    collect_info <- function(var,color = F){
        info <- annotated_df %>% dplyr::select(any_of(paste0(outcome_variable,'.',var))) %>% unlist %>% unique
        if(is.character(info)|is.factor(info)){
            ## categorical (discrete)
            type <- 'discrete'
            limit <- sort(info)
            ls <- list(type = type,limit = limit)
            if(color){
                fix_colors <- gg_color_hue(n = length(limit))
                ls <- list(type = type,limit = limit,fix_colors = fix_colors)
            }
        }else{
            ## numeric (continuous)
            type <- 'continuous'
            limit <- c(min(info,na.rm = T),max(info,na.rm = T))
            ls <- list(type = type,limit = limit)
        }
        ls
    }

    # (1) Size
    if(!is.null(size_variable)){
        size_info <- collect_info(size_variable)
    }
    # (2) Color
    if(!is.null(color_variable)){
        color_info <- collect_info(color_variable,color = T)
    }
    # (3) Alpha
    if(!is.null(alpha_variable)){
        alpha_info <- collect_info(alpha_variable)
    }

    if(sum(names(advanced_list) == 'label_info')>0){
        annotated_df <- suppressMessages(inner_join(annotated_df,advanced_list$label_info))
    }

    if(sum(colnames(annotated_df)=='label.color')==0){
        annotated_df$`label.color` <- 'black'
    }
    text_variable <- annotate_number_column

    #############################
    ## 4. Organize data for ggplot
    #############################
    organize_leaf <- function(data,variable,name){

        data %>%
            dplyr::select(c(id,x,y),label,any_of(paste0(outcome_variable,'.',variable))) %>%
            tidyr::gather(group,value,-c(id:label)) %>%
            dplyr::mutate(group = gsub(paste0('\\.',variable),'',group))  %>%
            dplyr::rename(!!name:=value)
    }


    distinct_classes <- unique(layout_all$y) %>% sort %>% as.numeric
    layout_ls <- lapply(distinct_classes,function(l){
        nm <- layout_all %>% dplyr::filter(y>=l) %>% dplyr::select(name) %>% unlist
        sub <- edges_dt %>% dplyr::filter(to %in% nm)
        if(nrow(sub) == 0){
            # root node
            sub <- data.frame(from = nm, to = nm)
        }
        layout <- create_layout(graph_from_data_frame(sub),layout = 'dendrogram')
    })
    names(layout_ls) <- distinct_classes

    combine_leaf_ls <- lapply(layout_ls,function(layout){
        leaf_node <- suppressMessages(inner_join(annotated_df,layout %>% dplyr::filter(leaf==T) %>% dplyr::rename(label=name) %>% dplyr::select(label)))
        leaf_node$label <- factor(leaf_node$label,levels = leaf_node$label[order(leaf_node$x)]) # must; order label

        combine_leaf <- lapply(c('size','color','alpha','text'),function(s){
            if(exists(paste0(s,'_variable'))){
                nm <- eval(as.name(paste0(s,'_variable')))
                if(!is.null(nm)){
                    dt <- organize_leaf(leaf_node,nm,s)
                }else{
                    dt <- NULL
                }
            }else{
                dt <- NULL
            }
            dt
        })
        combine_leaf <- suppressMessages(Reduce(inner_join,combine_leaf[unlist(sapply(combine_leaf, function(x) !is.null(x)))])) %>%
            mutate(group = factor(group,levels = outcome_variable))
    })


    #############################
    ## 5. Plot (ggplot2)
    #############################
    # (1) Base plot for non-leaf
    custom_nonleaf_plot <- function(layout,combine_leaf){

        # (1) Base plot for non-leaf
        nonleaf_base <- ggraph(layout, layout = 'dendrogram', circular = FALSE) +
            geom_node_point() +
            theme_minimal(base_size = font_size) +
            labs(y='') +
            coord_flip() +
            scale_x_discrete(breaks = seq(0,max(combine_leaf$x),1)) +
            theme(axis.ticks = element_blank(),
                  axis.text.x = element_blank(),
                  axis.title.y = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank()) +
            scale_y_reverse(expand = expansion(mult = c(0.1, 0.05))) +
            expand_limits(y=c(0,max(layout_all$y)))

        ## modify label position
        x_add <- 0 # nonleaf_label_pos
        y_add <- nonleaf_label_pos # nonleaf_point_gap
        hjust <- 1
        vjust <- 0.5

        if(nrow(layout)>1){
            nonleaf_node <- layout %>% filter(leaf == F) %>% dplyr::rename(label = name)
        }else{
            nonleaf_node <- layout %>% dplyr::rename(label = name)
        }
        nonleaf_data <- suppressMessages(inner_join(annotated_df %>% dplyr::select(-c(x,y)),nonleaf_node))

        if(edge_path_type=='link'){
            # label on the left
            nonleaf_base <- nonleaf_base +
                geom_edge_link(color = line_color,linetype = line_type) +
                geom_text(data = nonleaf_data,
                          aes(x = x + x_add,y = y + y_add,label = label),
                          color = nonleaf_data$`label.color`,
                          size = font_size*0.35,
                          hjust = hjust,vjust = vjust,show.legend  = FALSE)
        }else if(edge_path_type=='diagonal'){
            # label on the left
            nonleaf_base <- nonleaf_base +
                geom_edge_diagonal(color = line_color,linetype = line_type) +
                geom_text(data = nonleaf_data,
                          aes(x = x + x_add,y = y + y_add,label = label,color = `label.color`),
                          size = font_size*0.35,
                          color = nonleaf_data$`label.color`,
                          hjust = hjust,vjust = vjust,show.legend  = FALSE)
        }else if(edge_path_type=='elbow'){
            # label on the top
            x_add <- nonleaf_label_pos
            y_add <- nonleaf_point_gap
            hjust <- 0.5
            vjust <- 0
            nonleaf_base <- nonleaf_base +
                geom_edge_elbow(color = line_color,linetype = line_type) +
                geom_text(data = nonleaf_data,
                          aes(x = x + x_add,y = y + y_add,label = label), # y + nonleaf_label_pos
                          color = nonleaf_data$`label.color`,
                          size = font_size*0.35,
                          hjust = hjust,vjust = vjust,show.legend  = FALSE)

        }
        return(nonleaf_base)
    }

    # (2) Base plot for leaf
    custom_leaf_plot <- function(layout,combine_leaf){

        leaf_data <- annotated_df %>% filter(label %in% unique(combine_leaf$label))
        color <- ifelse(!is.null(color_variable),ifelse(sum(!is.na(combine_leaf$color))==0,'NULL','color'),'NULL')
        alpha <- ifelse(!is.null(alpha_variable),ifelse(sum(!is.na(combine_leaf$alpha))==0,'NULL','alpha'),'NULL')
        size <- ifelse(!is.null(size_variable),ifelse(sum(!is.na(combine_leaf$size))==0,'NULL','size'),'NULL')

        if(size == 'NULL' & alpha == 'NULL' & color == 'NULL'){
            leaf_base <- ggplot(data = combine_leaf,
                                aes_string(x = 'group', y = 'label',
                                           fill = ifelse(!is.null(color_variable),ifelse(sum(!is.na(combine_leaf$color))==0,'NULL','color'),'NULL'),
                                           alpha = ifelse(!is.null(alpha_variable),ifelse(sum(!is.na(combine_leaf$alpha))==0,'NULL','alpha'),'NULL'))) +
                geom_blank() +
                theme_minimal(base_size = font_size) +
                labs(x=NULL,y=NULL) +
                theme(axis.text = element_text(size = font_size,color = 'black'),
                      axis.text.y = element_text(hjust = 0),
                      axis.text.x = element_text(angle = x_axis_text_angle))
        }else{
            if(plot_type == 'heatmap'){
                # (2.1) heatmap
                leaf_base <- ggplot() +
                    geom_tile(data = combine_leaf,
                              aes_string(x = 'group', y = 'label',
                                         fill = ifelse(!is.null(color_variable),ifelse(sum(!is.na(combine_leaf$color))==0,'NULL','color'),'NULL'),
                                         alpha = ifelse(!is.null(alpha_variable),ifelse(sum(!is.na(combine_leaf$alpha))==0,'NULL','alpha'),'NULL')),
                              color = 'grey') +
                    theme_minimal(base_size = font_size) +
                    labs(x=NULL,y=NULL) +
                    theme(axis.text = element_text(size = font_size,color = 'black'),
                          axis.text.y = element_text(hjust = 0),
                          axis.text.x = element_text(angle = x_axis_text_angle))
            }else if(plot_type == 'bar'){
                # (2.2) barplot
                leaf_base <- ggplot(data = combine_leaf,
                                    aes_string(x = 'label', y = 'size',
                                               fill = ifelse(!is.null(color_variable),ifelse(sum(!is.na(combine_leaf$color))==0,'NULL','color'),'NULL'),
                                               alpha = ifelse(!is.null(alpha_variable),ifelse(sum(!is.na(combine_leaf$alpha))==0,'NULL','alpha'),'NULL'))) +
                    geom_bar(stat = 'identity', position = 'dodge') +
                    theme_minimal(base_size = font_size) +
                    labs(x=NULL,y=NULL) +
                    coord_flip() +
                    facet_wrap(~group,strip.position = 'bottom',nrow = 1) +
                    theme(axis.text = element_text(size = font_size,color = 'black'),
                          axis.text.y = element_text(hjust = 0),
                          axis.text.x = element_blank(),
                          strip.text = element_text(size = font_size, angle = x_axis_text_angle)) +

                    expand_limits(y=size_info$limit)

            }else{
                # (2.3) balloon point
                leaf_base <- ggplot() +
                    geom_point(data = combine_leaf,
                               aes_string(x = 'group', y = 'label',
                                          size = ifelse(!is.null(size_variable),ifelse(sum(!is.na(combine_leaf$size))==0,'NA','size'),'NULL'),
                                          color = ifelse(!is.null(color_variable),ifelse(sum(!is.na(combine_leaf$color))==0,'NULL','color'),'NULL'),
                                          alpha = ifelse(!is.null(alpha_variable),ifelse(sum(!is.na(combine_leaf$alpha))==0,'NULL','alpha'),'NULL'))) +
                    theme_minimal(base_size = font_size) +
                    labs(x=NULL,y=NULL) +
                    theme(axis.text = element_text(size = font_size,color = 'black'),
                          axis.text.y = element_text(hjust = 0),
                          axis.text.x = element_text(angle = x_axis_text_angle))
            }
        }

        leaf_base <- suppressWarnings(leaf_base +
                                          theme(axis.text.y = element_text(color = leaf_data %>% arrange(x) %>% dplyr::select(`label.color`) %>% unlist)))

        # (optional) add text
        text <- 'text'
        if(length(text)>0 & annotate_number){

            is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
            if(sum(is.wholenumber(combine_leaf[,text]),na.rm = T) == 0){

                combine_leaf$text <- round(combine_leaf[,text],2)
                leaf_scentific_idx <- which(abs(combine_leaf[,text])<1e-3|combine_leaf[,text]>1e3)
                if(length(leaf_scentific_idx)>0) combine_leaf$text[leaf_scentific_idx] <- formatC(combine_leaf[leaf_scentific_idx,text], format = "e", digits = 2)


            }else{
                combine_leaf$text <- combine_leaf[,text]
            }

            if(plot_type != 'bar'){
                leaf_base <- leaf_base +
                    geom_text(data = combine_leaf,aes(x = group, y = label,label = text),size = annotate_number_fontsize, color = annotate_number_color, alpha = 1, show.legend = F)
            }else{
                leaf_base <- leaf_base +
                    geom_text(data = combine_leaf,aes(x = label, y = size,label = text),size = annotate_number_fontsize,hjust = 0, color = annotate_number_color, alpha = 1, show.legend = F)
            }
        }

        # (3) Add scale options
        # a. size
        if(!is.null(size_variable)){
            if(size_info$type == 'continuous'){
                leaf_plot_size <- leaf_base +
                    scale_size_continuous(size_variable,
                                          range = c(1,12),
                                          limits = size_info$limit)

            }else{
                if(sum(names(advanced_list) == 'size_info')){

                    leaf_plot_size <- leaf_base +
                        scale_size_manual(size_variable,
                                          breaks = advanced_list$size_info$breaks,
                                          values = advanced_list$size_info$values,
                                          limits = size_info$limit)
                }else{

                    leaf_plot_size <- leaf_base +
                        scale_size_discrete(size_variable,
                                            range = c(1,12),
                                            breaks = size_info$limit,
                                            limits = size_info$limit)
                }
            }
        }else{

            leaf_plot_size <- leaf_base
        }


        # b. color
        if(!is.null(color_variable)){
            if(color_info$type == 'continuous'){
                if(plot_type == 'balloon'){
                    leaf_plot_color <- leaf_plot_size +
                        scale_color_distiller(color_variable,
                                              palette = ifelse(sum(names(advanced_list)=='palette')>0,advanced_list$palette[1],'Spectral'),
                                              na.value = "white",
                                              limits = color_info$limit)
                }else{
                    leaf_plot_color <- leaf_plot_size +
                        scale_fill_distiller(color_variable,
                                             palette = ifelse(sum(names(advanced_list)=='palette')>0,advanced_list$palette[1],'Spectral'),
                                             na.value = "white",
                                             limits = color_info$limit)
                }


            }else{
                if(sum(names(advanced_list) == 'color_info')>0){
                    color_breaks <- advanced_list$color_info$breaks
                    color_values <- advanced_list$color_info$values
                }else{
                    color_breaks <- color_info$limit
                    color_values <- color_info$fix_colors
                }

                if(plot_type == 'balloon'){
                    leaf_plot_color <- leaf_plot_size +
                        scale_color_manual(color_variable,
                                           breaks = color_breaks,
                                           limits = color_info$limit,
                                           values = color_values,
                                           na.value = "white")
                }else{
                    leaf_plot_color <- leaf_plot_size +
                        scale_fill_manual(color_variable,
                                          breaks = color_breaks,
                                          limits = color_info$limit,
                                          values = color_values,
                                          na.value = "white")
                }


            }
        }else{
            leaf_plot_color <- leaf_plot_size
        }


        # c. alpha
        if(!is.null(alpha_variable)){
            if(alpha_info$type == 'continuous'){

                leaf_plot <- leaf_plot_color +
                    scale_alpha_continuous(alpha_variable,
                                           limits = alpha_info$limit,
                                           range = c(0.1,1))


            }else{
                if(sum(names(advanced_list) == 'alpha_info')){

                    leaf_plot <- leaf_plot_color +
                        scale_alpha_manual(alpha_variable,
                                           limits = alpha_info$limit,
                                           breaks = advanced_list$alpha_info$breaks,
                                           values = advanced_list$alpha_info$values)
                }else{

                    leaf_plot <- leaf_plot_color +
                        scale_alpha_discrete(alpha_variable,
                                             limits = alpha_info$limit,
                                             breaks = alpha_info$limit,
                                             range = c(0.1,1))
                }

            }
        }else{

            leaf_plot <- leaf_plot_color
        }
        leaf_plot + theme(legend.position = "bottom",
                          legend.margin = margin(t=0,unit = 'lines'),
                          legend.box.just = 'left',
                          legend.box = 'horizontal')
    }


    # (4) Merge two panels as TreeCorTreat plot
    generate_plot <- function(nonleaf_plot,leaf_plot){
        comb_plot <- nonleaf_plot + leaf_plot
        if(sum(names(advanced_list)=='layout_widths')>0){
            new_plot <- comb_plot + patchwork::plot_layout(guides = "collect",widths = c(advanced_list$layout_widths[1],advanced_list$layout_widths[2]))
        }else{
            new_plot <- comb_plot + patchwork::plot_layout(guides = "collect",widths = c(1,6))
        }
        return(new_plot)
    }

    # left and right
    left.ls <- lapply(1:length(layout_ls),function(i){
        nonleaf_plot <- custom_nonleaf_plot(layout_ls[[i]],combine_leaf_ls[[i]])
    })
    right.ls <- lapply(1:length(layout_ls),function(i){
        leaf_plot <- custom_leaf_plot(layout_ls[[i]],combine_leaf_ls[[i]])
    })

    # plot at each level
    eachlevel.ls <- lapply(1:length(layout_ls),function(i){
        generate_plot(left.ls[[i]],right.ls[[i]]) &
            theme(legend.position = "bottom",
                  legend.margin = margin(t=0,unit = 'lines'),
                  legend.box.just = 'left',
                  legend.box = 'horizontal',
                  axis.text.y = element_text(hjust = 1))
    })
    names(eachlevel.ls) <- names(layout_ls)


    # combine plot
    ratio <- table(layout_all$y) %>% unlist
    final_left <- wrap_plots(left.ls,ncol = 1, heights = ratio) &
        theme(legend.position = 'right', # "bottom",
              legend.margin = margin(t=0,unit = 'lines'),
              legend.box.just = 'left',
              legend.box = 'vertical', # 'horizontal',
              axis.text.y = element_text(hjust = 1),
              plot.margin = margin(t=1,r=3,b=-1,l=5.5,unit = 'pt'))
    # remove x axis
    right.ls <- lapply(1:length(right.ls),function(i){
        g <- right.ls[[i]]
        if(i != length(right.ls)){
            gedit <- g + theme(axis.text.x = element_blank(),
                               strip.text = element_blank())
        }else{
            gedit <- g
        }
        gedit
    })
    final_right <- wrap_plots(right.ls,ncol = 1,guides = "collect", heights = ratio) &
        theme(legend.position = 'right', # "bottom",
              legend.margin = margin(t=0,unit = 'lines'),
              legend.box.just = 'left',
              legend.box = 'vertical', # 'horizontal',
              axis.text.y = element_text(hjust = 1))

    if(sum(names(advanced_list)=='layout_widths')>0){
        final_plot <- (final_left|final_right) + plot_layout(widths = c(advanced_list$layout_widths[1],advanced_list$layout_widths[2]))
    }else{
        final_plot <- (final_left|final_right)
    }


    if(plot){
        print(final_plot)
    }
    plots <- list(by_each_level = eachlevel.ls,
                  combine = final_plot)
}
