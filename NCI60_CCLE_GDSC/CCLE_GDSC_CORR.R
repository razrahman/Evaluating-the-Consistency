cat("\014") # Clear history
rm(list = ls() ) # Clear Environment variables
# Check if libraries are installed
paket <- function(pak){
  new_pak <- pak[!(pak %in% rownames(installed.packages()))]
  if (length(new_pak)) 
    install.packages(new_pak, dependencies = TRUE)
  sapply(pak, library, character.only = TRUE)
}
#"data.table", "ggplot2", "CORElearn", "h2o" , "lubridate", "Rcpp", "VennDiagram" , "extrafont"
listOfPackages <- c("data.table", "VennDiagram","ggplot2")

paket(listOfPackages)

CCLE_ALL <- fread("CCLE_ccle-gdsc2.csv", na.strings = c(""," "), 
                           nrows=-1L, verbose=getOption("datatable.verbose"),
                           header=T, stringsAsFactors=FALSE)

GDSC_ALL <- fread("GDSC_ccle-gdsc2.csv", na.strings = c(""," "), 
                  nrows=-1L, verbose=getOption("datatable.verbose"),
                  header=T, stringsAsFactors=FALSE)

GDSC_drugs <- fread("GDSC_ALL_drugs.csv", na.strings = c(""," "), 
                  nrows=-1L, verbose=getOption("datatable.verbose"),
                  header=T, stringsAsFactors=FALSE)

ccle_gdsc_cells <- fread("ccle_gdsc_cells.csv", na.strings = c(""," "), 
                    nrows=-1L, verbose=getOption("datatable.verbose"),
                    header=T, stringsAsFactors=FALSE)

ccle_gdsc_drugs <- fread("ccle_gdsc_drugs.csv", na.strings = c(""," "), 
                         nrows=-1L, verbose=getOption("datatable.verbose"),
                         header=T, stringsAsFactors=FALSE)

Direct_CORR=matrix(data = NA, ncol=2,nrow=length(ccle_gdsc_drugs[[2]]))
Direct_CORR2=matrix(data = NA, ncol=2,nrow=length(ccle_gdsc_drugs[[2]]))
Direct_CORR3=matrix(data = NA, ncol=2,nrow=length(ccle_gdsc_drugs[[2]]))
setkey(GDSC_ALL, Drug_Name)
setkey(CCLE_ALL,Compound)
for(i in seq_len(length(ccle_gdsc_drugs[[2]]))){
  CC1 = which(ccle_gdsc_drugs[[2]][i]==GDSC_drugs, arr.in=TRUE)
  gdsc_drugs_i = GDSC_drugs[[2]][CC1[1]]
  GDSC_AA = GDSC_ALL[.(gdsc_drugs_i), LN_IC50:Cell_Line]
  GDSC_M = matrix(c(GDSC_AA[,Cell_Line],GDSC_AA[,MAX_CONC_MICROMOLAR*10^-6],GDSC_AA[,exp(LN_IC50)*10^-6]),ncol=3)
  GDSC_M2 = GDSC_M[complete.cases(GDSC_M), ]
  if (mean(as.numeric(GDSC_M2[,2]))*length(GDSC_M2[,2])!=sum(as.numeric(GDSC_M2[,2]))) stop("Maximum concentration is not same")
  
  CCLE_AA = CCLE_ALL[.(ccle_gdsc_drugs[i][[2]]), Cell_Line:IC50]
  CCLE_M = matrix(c(CCLE_AA[,Cell_Line], CCLE_AA[,IC50*10^-6]),ncol=2)
  CCLE_M2 = CCLE_M[complete.cases(CCLE_M), ]
  
  if (length(CCLE_M2[,1])>length(ccle_gdsc_cells[[1]]) || length(CCLE_M2[,1])>length(ccle_gdsc_cells[[1]]))
      stop("Number of available cell line is higher")
  
  Common_Cell = intersect(CCLE_M2[,1], GDSC_M2[,1])
  CCLE_index = match(Common_Cell,CCLE_M2[,1])
  GDSC_index = match(Common_Cell,GDSC_M2[,1])
  IC50 = matrix(as.numeric(c(CCLE_M2[CCLE_index,2],GDSC_M2[GDSC_index,3])),ncol=2)
  Direct_CORR[i,] =  c(cor(IC50[,1],IC50[,2],method = "pearson"),cor(IC50[,1],IC50[,2],method = "spearman"))
  #names(IC50) <- c("CCLE_IC50","GDSC_IC50")
  topConc = min(mean(as.numeric(GDSC_M2[,2])),8*10^(-6))
  lowConc = max(min(IC50[,2]),0.0025*10^(-6))
  IC50_2 = pmax(pmin(IC50,topConc),lowConc)
  Direct_CORR2[i,] =  c(cor(IC50_2[,1],IC50_2[,2],method = "pearson"),cor(IC50_2[,1],IC50_2[,2],method = "spearman"))
  
  IC50_3 = (-(log10(IC50_2)-log10(topConc)))/max(-(log10(IC50_2)-log10(topConc)))
  Direct_CORR3[i,] =  c(cor(IC50_3[,1],IC50_3[,2],method = "pearson"),cor(IC50_3[,1],IC50_3[,2],method = "spearman"))
}
Direct_CORR <- matrix(c(matrix(ccle_gdsc_drugs[[2]],ncol=1),Direct_CORR),ncol=3)
colnames(Direct_CORR)<- c("Drug Name","Pearson Correlation","Spearman Correlation")
Direct_CORR2 <- matrix(c(matrix(ccle_gdsc_drugs[[2]],ncol=1),Direct_CORR2),ncol=3)
colnames(Direct_CORR2)<- c("Drug Name","Pearson Correlation","Spearman Correlation")
Direct_CORR3 <- matrix(c(matrix(ccle_gdsc_drugs[[2]],ncol=1),Direct_CORR3),ncol=3)
colnames(Direct_CORR3)<- c("Drug Name","Pearson Correlation","Spearman Correlation")
fwrite(data.table(Direct_CORR),file="Direct_CORR.csv", 
       showProgress = getOption("datatable.showProgress"),
       verbose = getOption("datatable.verbose"),
       sep = ",", col.names = T, row.names = T)