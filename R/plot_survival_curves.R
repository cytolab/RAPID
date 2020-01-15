#' @export
#' @import survival
#' @import survminer
#' @import plyr

plot_survival_curves <- function(survival_stats, cluster_and_patient_df, survival_data, clinical_col, status_col){

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

  # divide patients into high and low groups based on subset abundance
  abundance.groups <- list()
  for(e in 1:ncol(all.subset.abundances)){
    abundance.groups[[e]] <- (single.subset.abundance.data[[e]]>abundance.IQR[[e]])
    abundance.groups[[e]][abundance.groups[[e]]==TRUE] = 1}             # 1 = high, 0 = low

  survival_plot <- list()
  count = 1
  colnames(survival_data)[clinical_col] <- "clinicalstat"
  colnames(survival_data)[status_col] <- "status"

  dir.create(file.path(getwd(), "output files"), showWarnings = FALSE)

  sink(paste("./output files/",strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),"_subset_survival_stats.txt",sep=""))

  for (s in 1:ncol(all.subset.abundances)){
    survival_stat = survival_stats[[s]]$coefficients
    CI = survival_stats[[s]]$conf.int
    Group = factor(abundance.groups[[s]], levels = c(0,1), labels = c("Low", "High"))
    survival_data1 = cbind(survival_data,Group)
    model1 <- survfit(Surv(clinicalstat, status) ~ Group, data=survival_data1)
    print(paste("Subset #", colnames(all.subset.abundances)[s], " (p = ", round(survival_stat[,5],3), ", HR = ",round(survival_stat[,2],3),", CI[",round(CI[,3],3),",",round(CI[,4],3),"])", sep=""))
    if (survival_stat[,5] <= 0.05){
      survival_plot[[count]] <- ggsurvplot(model1, data=survival_data1, title = paste("Subset #", colnames(all.subset.abundances)[s], " (p = ", round(survival_stat[,5],3), ", HR = ",round(survival_stat[,2],3),", CI[",round(CI[,3],3),",",round(CI[,4],3),"])", sep=""),conf.int=F, pval=F, risk.table=T, tables.y.text = FALSE, legend.labs = c("Low", "High"), legend.title = "Group",censor.shape=124)
      count = count+1}}

  sink()
  return(survival_plot)}
