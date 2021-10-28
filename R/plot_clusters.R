#' @export
#' @import RColorBrewer
#' @import tidyverse

# plot on the t-SNE axes and color by ideal FlowSOM number of clusters
plot_clusters <- function (x, y,clusters,xlab,ylab,legendtitle,title){

  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                             rownames(qual_col_pals)))
  max_x = max(c(max(x),abs(min(x))))
  max_y = max(c(max(y),abs(min(y))))
  plot <- data.frame(x = x, y = y, col = as.factor(clusters))
  FlowSOM_clusters_on_viSNE <- ggplot(plot) + geom_point(aes(x=x, y=y, color=col)) + labs(color = legendtitle, x = xlab, y = ylab, title = title) +
    scale_color_manual(values = col_vector)+ guides(colour = guide_legend(override.aes = list(size=5)))+theme_bw() + theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())+ylim(-max_y,max_y) +xlim(-max_x,max_x)
  return(FlowSOM_clusters_on_viSNE)}
