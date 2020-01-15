#' @export
#' @import tidyverse

plot_sig_clusters <- function (x, y, clusters, survival_stats, xlab, ylab, legendtitle, title){

  data = as.data.frame(cbind(x,y,as.factor(clusters)))
  colnames(data)[ncol(data)] <- "cluster"

  GNP05 <- list()
  GNP01 <- list()
  GNP1 <- list()
  GPP05 <- list()
  GPP01 <- list()
  GPP1 <- list()

  countgnp01 = 0
  countgpp01 = 0
  countgnp05 = 0
  countgpp05 = 0
  countgnp1 = 0
  countgpp1 = 0

  for (z in 1:length(survival_stats)){
    survival_stats[[z]] = survival_stats[[z]]$coefficients
    if (survival_stats[[z]][,c(5)] <= 0.01 & survival_stats[[z]][,c(2)] <= 1){
      countgpp01 = countgpp01 + 1
      GPP01[countgpp01] = z}

    else if (survival_stats[[z]][,c(5)] <= 0.05 & survival_stats[[z]][,c(2)] <= 1){
      countgpp05 = countgpp05 + 1
      GPP05[countgpp05] = z}

    else if (survival_stats[[z]][,c(5)] <= 0.1 & survival_stats[[z]][,c(2)] <= 1){
      countgpp1 = countgpp1 + 1
      GPP1[countgpp1] = z}

    if (survival_stats[[z]][,c(5)] <= 0.01 & survival_stats[[z]][,c(2)] >= 1){
      countgnp01 = countgnp01 + 1
      GNP01[countgnp01] = z}

    else if (survival_stats[[z]][,c(5)] <= 0.05 & survival_stats[[z]][,c(2)] >= 1){
      countgnp05 = countgnp05 + 1
      GNP05[countgnp05] = z}

    else if (survival_stats[[z]][,c(5)] <= 0.1 & survival_stats[[z]][,c(2)] >= 1){
      countgnp1 = countgnp1 + 1
      GNP1[countgnp1] = z}}

  data$status <- "p>0.1"
  if (length(GNP01) >= 1)
  {for (X in (1:length(GNP01))){
    data$status[data$cluster == GNP01[X]] <- "p<0.01 HR>1"}}

  if (length(GPP01) >= 1)
  {for (Y in (1:length(GPP01))){
    data$status[data$cluster == GPP01[Y]] <- "p<0.01 HR<1"}}

  if (length(GPP05) >= 1)
  {for (Z in (1:length(GPP05))){
    data$status[data$cluster == GPP05[Z]] <- "p<0.05 HR<1"}}

  if (length(GNP05) >= 1)
  {for (W in (1:length(GNP05))){
    data$status[data$cluster == GNP05[W]] <- "p<0.05 HR>1"}}

  if (length(GPP1) >= 1)
  {for (V in (1:length(GPP1))){
    data$status[data$cluster == GPP1[V]] <- "p<0.1 HR<1"}}

  if (length(GNP1) >= 1)
  {for (H in (1:length(GNP1))){
    data$status[data$cluster == GNP1[H]] <- "p<0.1 HR>1"}}

  cols <- c("p>0.1" ="lightgray" ,"p<0.01 HR>1" = "darkred","p<0.05 HR>1" = "red","p<0.1 HR>1" = "darkgray","p<0.01 HR<1" = "darkblue","p<0.05 HR<1" = "blue","p<0.1 HR<1"= "darkgray")

  xvalues = c(max(x),abs(min(x)))
  max_x = max(xvalues)
  yvalues = c(max(y),abs(min(y)))
  max_y = max(yvalues)

  plot <- data.frame(x = data[,1], y = data[,2], col = as.factor(data$status))
  clinical_clusters <- ggplot(plot)+ geom_point(aes(x=data[,1], y=data[,2], col = as.factor(data$status))) + scale_color_manual(breaks = c("p<0.01 HR<1","p<0.05 HR<1","p<0.1 HR<1", "p>0.1","p<0.1 HR>1","p<0.05 HR>1","p<0.01 HR>1"),values = cols)+ guides(colour = guide_legend(override.aes = list(size=5)))+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(color = legendtitle, x = xlab, y = ylab, title = title) +ylim(-max_y,max_y) +xlim(-max_x,max_x)

  return(clinical_clusters)}
