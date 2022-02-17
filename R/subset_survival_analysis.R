#' @export
#' @import survival
#' @import survminer
#' @import plyr

subset_survival_analysis <- function (cluster_and_patient_df, survival_data, clinical_col, status_col){

  cluster_and_patient_df[,1] <- as.factor(cluster_and_patient_df[,1])
  cluster_and_patient_df[,2] <- as.factor(cluster_and_patient_df[,2])
  patient.subsets <- split(cluster_and_patient_df,cluster_and_patient_df[,2])

  subset.abundance <- list()
  # find subset abundances for each patient
  for (j in 1:length(patient.subsets)){
    subset.abundance[[j]] = (summary(patient.subsets[[j]][,1]))*100/nrow(patient.subsets[[j]])}

  # this variable is the percentages of GBM cells in each subset for every patient (each column represents one of the subsets)
  all.subset.abundances = ldply(subset.abundance,rbind)
  all.subset.abundances[is.na(all.subset.abundances)]<-0

  # split data frame by column, that is to say, separate out all abundances for all patients per subset
  single.subset.abundance.data <- list()
  for(d in 1:ncol(all.subset.abundances)){
    single.subset.abundance.data[[d]] = all.subset.abundances[,c(d)]}

  # find abundance IQR for each subset
  abundance.IQR <- list()
  for (b in 1:ncol(all.subset.abundances)){
    abundance.IQR[[b]] = IQR(all.subset.abundances[,c(b)])}

  Abundance_IQR = as.vector(unlist(abundance.IQR))
  rownames(all.subset.abundances)<- c(1:nrow(all.subset.abundances))

  # divide patients into high and low groups based on subset abundance
  abundance.groups <- list()
  for(e in 1:ncol(all.subset.abundances)){
    abundance.groups[[e]] <- (single.subset.abundance.data[[e]]>abundance.IQR[[e]])
    abundance.groups[[e]][abundance.groups[[e]]==TRUE] = 1}             # 1 = high, 0 = low

  cox.summary<- list()
  count = 1
  high.low.groups <- list()
  low.median = vector()
  high.median = vector()
  for (s in 1:ncol(all.subset.abundances)){
    Group <- factor(abundance.groups[[s]], levels = c(0,1), labels = c("Low", "High"))
    survival_data1 = cbind(survival_data,Group)
    high.low.groups[[s]] = split(survival_data1,survival_data1$Group)
    low.median[[s]] = median(high.low.groups[[s]][["Low"]][[clinical_col]])
    high.median[[s]] = median(high.low.groups[[s]][["High"]][[clinical_col]])
    model2 <- coxph(Surv(survival_data[,clinical_col], survival_data[,status_col]) ~ Group, data=survival_data1)
    cox.summary[[colnames(all.subset.abundances)[s]]] = summary(model2)}

  Median_High_Group = high.median
  Median_Low_Group = low.median
  all.abundance.data = rbind(all.subset.abundances,Abundance_IQR, Median_High_Group,Median_Low_Group)
  colnames(all.abundance.data) <- paste("Subset", colnames(all.abundance.data), sep = "_")
  rownames(all.abundance.data)[(nrow(all.abundance.data)-2):(nrow(all.abundance.data))]<-c("Abundance_IQR", "Median_High_Group","Median_Low_Group")


  dir.create(file.path(getwd(), "output files"), showWarnings = FALSE)
  write.csv(all.abundance.data, file = paste("./output files/",strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),"_subset_abundances_withIQR.csv",sep=""))

  return(cox.summary)}
