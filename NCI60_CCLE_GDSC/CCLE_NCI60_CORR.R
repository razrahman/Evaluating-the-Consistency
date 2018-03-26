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

CCLE_ALL <- fread("CCLE_nci-ccle2.csv", na.strings = c(""," "), 
                  nrows=-1L, verbose=getOption("datatable.verbose"),
                  header=T, stringsAsFactors=FALSE)

NCI60_ALL <- fread("NCI60_nci-ccle2.csv", na.strings = c(""," "), 
                  nrows=-1L, verbose=getOption("datatable.verbose"),
                  header=T, stringsAsFactors=FALSE)

NCI60_drugs <- fread("nsc_names_nci60.csv", na.strings = c(""," "), 
                    nrows=-1L, verbose=getOption("datatable.verbose"),
                    header=T, stringsAsFactors=FALSE)

ccle_NCI60_cells <- fread("nci_ccle_cells.csv", na.strings = c(""," "), 
                         nrows=-1L, verbose=getOption("datatable.verbose"),
                         header=T, stringsAsFactors=FALSE)

ccle_NCI60_drugs <- fread("nci_ccle_drugs.csv", na.strings = c(""," "), 
                         nrows=-1L, verbose=getOption("datatable.verbose"),
                         header=T, stringsAsFactors=FALSE)

Direct_CORR=matrix(data = NA, ncol=2,nrow=length(ccle_NCI60_drugs[[2]]))
Direct_CORR2=matrix(data = NA, ncol=2,nrow=length(ccle_NCI60_drugs[[2]]))
Direct_CORR3=matrix(data = NA, ncol=2,nrow=length(ccle_NCI60_drugs[[2]]))
aa1=NULL
aa2=NULL
aa3=NULL
setkey(NCI60_ALL, NSC)
setkey(CCLE_ALL,Compound)
setkey(NCI60_drugs,Drug_Name)
for(i in seq_len(length(ccle_NCI60_drugs[[2]]))){
  CC1 = which(ccle_NCI60_drugs[[2]][i]==NCI60_drugs[,Drug_Name], arr.in=TRUE)
  NCI60_drugs_i = NCI60_drugs[[2]][CC1[1]]
  NCI60_AA = NCI60_ALL[.(NCI60_drugs_i), NSC:NLOGGI50]
  NCI60_M = NCI60_AA[complete.cases(NCI60_AA), ]
  if (unique(NCI60_M[,2])!="M") 
    stop("Maximum concentration unit is not same")
  CC2=NULL
  for (ii in 1:length(unique(NCI60_M[,5][[1]]))){
    CC2[ii]= match(unique(NCI60_M[,5][[1]])[ii],NCI60_M[,5][[1]])
  }
  NCI60_M1 = NCI60_M[CC2,]
  NCI60_M2 = matrix(c(NCI60_M1[,CELL],NCI60_M1[,10^LCONC],NCI60_M1[,10^-NLOGGI50]),ncol=3)
  #if (mean(as.numeric(NCI60_M2[,2]))*length(NCI60_M2[,2])!=sum(as.numeric(NCI60_M2[,2]))) 
   # stop("Maximum concentration is not same")
  
  CCLE_AA = CCLE_ALL[.(ccle_NCI60_drugs[i][[2]]), Cell_Line:IC50]
  CCLE_M = matrix(c(CCLE_AA[,Cell_Line], CCLE_AA[,IC50*10^-6]),ncol=2)
  CCLE_M2 = CCLE_M[complete.cases(CCLE_M), ]
  
  if (length(CCLE_M2[,1])>length(ccle_NCI60_cells[[1]]) || length(CCLE_M2[,1])>length(ccle_NCI60_cells[[1]]))
    stop("Number of available cell line is higher")
  
  Common_Cell = intersect(CCLE_M2[,1], NCI60_M2[,1])
  CCLE_index = match(Common_Cell,CCLE_M2[,1])
  NCI60_index = match(Common_Cell,NCI60_M2[,1])
  IC50 = matrix(as.numeric(c(CCLE_M2[CCLE_index,2],NCI60_M2[NCI60_index,3])),ncol=2)
  Direct_CORR[i,] =  c(cor(IC50[,1],IC50[,2],method = "pearson"),cor(IC50[,1],IC50[,2],method = "spearman"))
  #names(IC50) <- c("CCLE_IC50","NCI60_IC50")
  topConc = min(min(as.numeric(NCI60_M2[,2])),8*10^(-6))
  lowConc = max(min(IC50[,2]),0.0025*10^(-6))
  IC50_2 = pmax(pmin(IC50,topConc),lowConc)#pmin(IC50,topConc)#pmax(pmin(IC50,topConc),lowConc)
  Direct_CORR2[i,] =  c(cor(IC50_2[,1],IC50_2[,2],method = "pearson"),cor(IC50_2[,1],IC50_2[,2],method = "spearman"))
  
  IC50_3 = (-(log10(IC50_2)-log10(topConc)))/max(-(log10(IC50_2)-log10(topConc)))
  Direct_CORR3[i,] =  c(cor(IC50_3[,1],IC50_3[,2],method = "pearson"),cor(IC50_3[,1],IC50_3[,2],method = "spearman"))
  
  aa1=rbind(aa1,IC50)
  aa2=rbind(aa2,IC50_2)
  aa3=rbind(aa3,IC50_3)
}
c(cor(aa1[,1],aa1[,2],method = "pearson"),cor(aa1[,1],aa1[,2],method = "spearman"))
c(cor(aa2[,1],aa2[,2],method = "pearson"),cor(aa2[,1],aa2[,2],method = "spearman"))
c(cor(aa3[,1],aa3[,2],method = "pearson"),cor(aa3[,1],aa3[,2],method = "spearman"))
Direct_CORR <- matrix(c(matrix(ccle_NCI60_drugs[[2]],ncol=1),Direct_CORR),ncol=3)
colnames(Direct_CORR)<- c("Drug Name","Pearson Correlation","Spearman Correlation")
Direct_CORR2 <- matrix(c(matrix(ccle_NCI60_drugs[[2]],ncol=1),Direct_CORR2),ncol=3)
colnames(Direct_CORR2)<- c("Drug Name","Pearson Correlation","Spearman Correlation")
Direct_CORR3 <- matrix(c(matrix(ccle_NCI60_drugs[[2]],ncol=1),Direct_CORR3),ncol=3)
colnames(Direct_CORR3)<- c("Drug Name","Pearson Correlation","Spearman Correlation")
fwrite(data.table(Direct_CORR3),file="Direct_CORR.csv", 
       showProgress = getOption("datatable.showProgress"),
       verbose = getOption("datatable.verbose"),
       sep = ",", col.names = T, row.names = T)