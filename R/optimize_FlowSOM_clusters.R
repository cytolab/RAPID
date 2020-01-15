#' @export
#' @import FlowSOM
#' @import matrixStats
#' @import Biobase

optimize_FlowSOM_clusters <- function (FlowSOM_input, marker_data, N, seed) {

  # create flowFrame from input data
  mat.input <- as.matrix(FlowSOM_input)
  metadata <- data.frame(name = dimnames(mat.input)[[2]], desc = dimnames(mat.input)[[2]])
  metadata$range <- apply(apply(mat.input, 2, range), 2, diff)
  metadata$minRange <- apply(mat.input, 2, min)
  metadata$maxRange <- apply(mat.input, 2, max)
  input.flowframe <- new("flowFrame", exprs=mat.input,parameters = AnnotatedDataFrame(metadata))
  # create empty data frames for the for loop to store values
  flowSOM.clusters <- list()
  data.with.clusters <- list()
  subset.data <- list()
  subset.IQR.markers <- list()
  mean.subset.IQR.markers <- list()

  TEST.CLUSTER.NUM = N-4
  # for loop to implement FlowSOM on each tSNE data set
  for  (k in 1:TEST.CLUSTER.NUM) {
    CLUSTER.NUMBER = k+4

    # implement the FlowSOM function on the data
    fSOM <- FlowSOM(input.flowframe, scale = TRUE, colsToUse = c(1:length(FlowSOM_input)), nClus    = CLUSTER.NUMBER, seed = seed)

    # create a cluster matrix
    flowSOM.clusters[[k]] <- as.matrix(fSOM[[2]][fSOM[[1]]$map$mapping[,1]])
    data.with.clusters[[k]] = cbind(marker_data, flowSOM.clusters[[k]])
    colnames(data.with.clusters[[k]])[colnames(data.with.clusters[[k]]) ==  'flowSOM.clusters[[k]]'] <- 'cluster'
    subset.data[[k]] <- split(data.with.clusters[[k]], data.with.clusters[[k]]$cluster)
    FlowSOM.subset.data = subset.data[[k]]

    # finding average IQR for each marker in each subset
    for (l in 1:length(FlowSOM.subset.data)){
      subset.matrix = as.matrix(FlowSOM.subset.data[[l]][-c(length(FlowSOM.subset.data[[l]]))])
      subset.IQR.markers[[l]] = colIQRs(subset.matrix)}
    mean.subset.IQR.markers[[k]] = rowMeans(as.data.frame(subset.IQR.markers))
  }

  # minimization of sums of average IQR for each marker
  df = as.data.frame(t(data.frame(mean.subset.IQR.markers)))
  colnames(df)<-colnames(marker_data)
  chosen.index <- list()

  # find the optimal or ideal number of subsets, or FlowSOM clusters
  for (m in 1:ncol(df)){
    chosen.index[[m]] = which.min((df[,c(m)]))}

  chosen.indexes = do.call(rbind, chosen.index)
  IDEAL.INDEX = round(median(chosen.indexes))
  IDEAL.NUM.SUBSETS = IDEAL.INDEX + 4
  print(paste("Ideal number of FlowSOM clusters was calculated to be ", IDEAL.NUM.SUBSETS, sep = ""))
  data.with.clusters[[IDEAL.INDEX]]$cluster = as.numeric(as.vector(data.with.clusters[[IDEAL.INDEX]]$cluster))
  return(data.with.clusters[[IDEAL.INDEX]]$cluster)}
