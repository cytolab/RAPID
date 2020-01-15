#' @export


# plot GNP and GPP cell subsets back onto the t-SNE axes with color represent degree to which the clusters stratify patient outcome as well as the degree of risk for patients with higher abundances of those specific subsets
significant_clusters <- function (survival_stats,cluster_data){
  GNP <- vector()
  GPP <- vector()
  countgnp = 0
  countgpp = 0

  dir.create(file.path(getwd(), "output files"), showWarnings = FALSE)

  sink(paste("./output files/",strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),"_prognostic_subset_values.txt",sep=""))

  for (z in 1:length(survival_stats)){
    survival_stats[[z]] = survival_stats[[z]]$coefficients
    if (survival_stats[[z]][,5] <= 0.05 & survival_stats[[z]][,2] <= 1){
      print(paste("Positive Prognostic Subset #", z))
      countgpp = countgpp + 1
      GPP[countgpp] = z}

    if (survival_stats[[z]][,5] <= 0.05 & survival_stats[[z]][,2] >= 1){
      print(paste("Negative Prognostic Subset #",z))
      countgnp = countgnp + 1
      GNP[countgnp] = z}
  }

  sink()

  cluster_data = as.data.frame(cluster_data)
  colnames(cluster_data)[1]<- "Cluster"
  cluster_data$Prog_Status <-1
  if (length(GNP) >= 1)
  {for (x in (1:length(GNP))){
    cluster_data$Prog_Status[cluster_data[,1] == GNP[x]] <- 3}}

  if (length(GPP) >= 1)
  {for (y in (1:length(GPP))){
    cluster_data$Prog_Status[cluster_data[,1] == GPP[y]] <- 2}}

  return(cluster_data)}
