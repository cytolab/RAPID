#' @export
#' @import RColorBrewer
#' @import tidyverse

# plot on the t-SNE axes and color by ideal FlowSOM number of clusters
plot_clusters <- function (x, y,clusters,xlab,ylab,legendtitle,title){

  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                             rownames(qual_col_pals)))
  col_vector = col_vector[-c(4,17,19,27,29:45)]
  max_x = max(c(max(x),abs(min(x))))
  max_y = max(c(max(y),abs(min(y))))
  plot <- data.frame(x = x, y = y, col = clusters)
  values = col_vector
  cols <- c("1" =values[1] ,"2" =values[2],"3" = values[3],"4" = values[4],"5" = values[5],"6" = values[6],"7" = values[7],"8" = values[8],"9" = values[9],
            "10" =values[10] ,"11" =values[11], "12" =values[12],"13" = values[13],"14" = values[14],"15" = values[15],"16" = values[16],"17" = values[17],
            "18" = values[18],"19" = values[19],"20" =values[20] ,"21" =values[21], "22" =values[22],"23" = values[23],"24" = values[24],"25" = values[25],
            "26" = values[26],"27" = values[27],"28" = values[28],"29" = values[29],"30" =values[30] ,"31" =values[31] ,"32" =values[32],"33" = values[33],
            "34" = values[34],"35" = values[35],"36" = values[36],"37" = values[37],"38" = values[38],"39" = values[39],"40" =values[40],"41" =values[41] ,"42" =values[42],"43" = values[43],
            "44" = values[44],"45" = values[45],"46" = values[46],"47" = values[47],"48" = values[48],"39" = values[49],"50" =values[50])

  plot <- data.frame(x = x, y = y, col = as.factor(clusters))
  FlowSOM_clusters_on_viSNE <- ggplot(plot) + geom_point(aes(x=x, y=y, color=col)) + labs(color = legendtitle, x = xlab, y = ylab, title = title) +
    scale_color_manual(breaks = c(1:50),values = cols)+ guides(colour = guide_legend(override.aes = list(size=5)))+theme_bw() + theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())+ylim(-max_y,max_y) +xlim(-max_x,max_x)
  return(FlowSOM_clusters_on_viSNE)}
