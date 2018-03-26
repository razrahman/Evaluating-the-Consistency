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

GDSC_ALL <- fread("GDSC_nci-gdsc2.csv", na.strings = c(""," "), 
                  nrows=-1L, verbose=getOption("datatable.verbose"),
                  header=T, stringsAsFactors=FALSE)

GDSC_drugs <- fread("GDSC_ALL_drugs.csv", na.strings = c(""," "), 
                    nrows=-1L, verbose=getOption("datatable.verbose"),
                    header=T, stringsAsFactors=FALSE)

NCI60_ALL <- fread("NCI60_nci-gdsc2.csv", na.strings = c(""," "), 
                   nrows=-1L, verbose=getOption("datatable.verbose"),
                   header=T, stringsAsFactors=FALSE)

NCI60_drugs <- fread("nsc_names_nci60.csv", na.strings = c(""," "), 
                     nrows=-1L, verbose=getOption("datatable.verbose"),
                     header=T, stringsAsFactors=FALSE)

nci60_gdsc_cells <- fread("nci_gdsc_cells.csv", na.strings = c(""," "), 
                         nrows=-1L, verbose=getOption("datatable.verbose"),
                         header=T, stringsAsFactors=FALSE)

nci60_gdsc_drugs <- fread("nci_gdsc_drugs.csv", na.strings = c(""," "), 
                         nrows=-1L, verbose=getOption("datatable.verbose"),
                         header=T, stringsAsFactors=FALSE)

Direct_CORR=matrix(data = NA, ncol=2,nrow=length(nci60_gdsc_drugs[[2]]))
Direct_CORR2=matrix(data = NA, ncol=2,nrow=length(nci60_gdsc_drugs[[2]]))
Direct_CORR3=matrix(data = NA, ncol=2,nrow=length(nci60_gdsc_drugs[[2]]))
aa1=NULL
aa2=NULL
aa3=NULL
setkey(GDSC_ALL, Drug_Name)
setkey(NCI60_ALL, NSC)
setkey(NCI60_drugs,Drug_Name)
for(i in seq_len(length(nci60_gdsc_drugs[[2]]))){
  ## GDSC data ##
  CC1 = which(nci60_gdsc_drugs[[2]][i]==GDSC_drugs, arr.in=TRUE)
  gdsc_drugs_i = GDSC_drugs[[2]][CC1[1]]
  GDSC_AA = GDSC_ALL[.(gdsc_drugs_i), LN_IC50:Cell_Line]
  GDSC_M = matrix(c(GDSC_AA[,Cell_Line],GDSC_AA[,MAX_CONC_MICROMOLAR*10^-6],GDSC_AA[,exp(LN_IC50)*10^-6]),ncol=3)
  GDSC_M2 = GDSC_M[complete.cases(GDSC_M), ]
  #if (mean(as.numeric(GDSC_M2[,2]))*length(GDSC_M2[,2])!=sum(as.numeric(GDSC_M2[,2]))) stop("Maximum concentration is not same")
  
  ## NCI60 data ##
  CC2 = which(nci60_gdsc_drugs[[2]][i]==NCI60_drugs[,Drug_Name], arr.in=TRUE)
  NCI60_drugs_i = NCI60_drugs[[2]][CC2[1]]
  NCI60_AA = NCI60_ALL[.(NCI60_drugs_i), NSC:NLOGGI50]
  NCI60_M = NCI60_AA[complete.cases(NCI60_AA), ]
  N_molar=1
  if (unique(NCI60_M[,2])!="M") {
    if (unique(NCI60_M[,2])=="u" || unique(NCI60_M[,2])=="V"){
      N_molar=10^(-6)
    }else stop("Maximum concentration unit is not same")
  }
    
  CC3=NULL
  for (ii in 1:length(unique(NCI60_M[,5][[1]]))){
    CC3[ii]= match(unique(NCI60_M[,5][[1]])[ii],NCI60_M[,5][[1]])
  }
  NCI60_M1 = NCI60_M[CC3,]
  NCI60_M2 = matrix(c(NCI60_M1[,CELL],NCI60_M1[,(10^LCONC)*N_molar],NCI60_M1[,(10^-NLOGGI50)*N_molar]),ncol=3)
  #if (mean(as.numeric(NCI60_M2[,2]))*length(NCI60_M2[,2])!=sum(as.numeric(NCI60_M2[,2]))) 
  # stop("Maximum concentration is not same")
  
  if (length(NCI60_M2[,1])>length(nci60_gdsc_cells[[1]]) || length(NCI60_M2[,1])>length(nci60_gdsc_cells[[1]]))
    stop("Number of available cell line is higher")
  
  Common_Cell = intersect(NCI60_M2[,1], GDSC_M2[,1])
  NCI60_index = match(Common_Cell,NCI60_M2[,1])
  GDSC_index = match(Common_Cell,GDSC_M2[,1])
  IC50 = matrix(as.numeric(c(NCI60_M2[NCI60_index,3],GDSC_M2[GDSC_index,3])),ncol=2)
  aa1=rbind(aa1,IC50)
  Direct_CORR[i,] =  c(cor(IC50[,1],IC50[,2],method = "pearson"),cor(IC50[,1],IC50[,2],method = "spearman"))
  #names(IC50) <- c("NCI60_IC50","GDSC_IC50")
  topConc = min(min(as.numeric(GDSC_M2[,2])),min(as.numeric(NCI60_M2[,2])))
  lowConc = max(min(IC50[,2]),min(IC50[,1]))
  IC50_2 = pmax(pmin(IC50,topConc),lowConc)
  Direct_CORR2[i,] =  c(cor(IC50_2[,1],IC50_2[,2],method = "pearson"),cor(IC50_2[,1],IC50_2[,2],method = "spearman"))
  aa2=rbind(aa2,IC50_2)
  
  IC50_3 = (-(log10(IC50_2)-log10(topConc)))/max(-(log10(IC50_2)-log10(topConc)))
  Direct_CORR3[i,] =  c(cor(IC50_3[,1],IC50_3[,2],method = "pearson"),cor(IC50_3[,1],IC50_3[,2],method = "spearman"))
  aa3=rbind(aa3,IC50_3)
}
c(cor(aa1[,1],aa1[,2],method = "pearson"),cor(aa1[,1],aa1[,2],method = "spearman"))
c(cor(aa2[,1],aa2[,2],method = "pearson"),cor(aa2[,1],aa2[,2],method = "spearman"))
c(cor(aa3[,1],aa3[,2],method = "pearson"),cor(aa3[,1],aa3[,2],method = "spearman"))
Direct_CORR <- matrix(c(matrix(nci60_gdsc_drugs[[2]],ncol=1),Direct_CORR),ncol=3)
colnames(Direct_CORR)<- c("Drug Name","Pearson Correlation","Spearman Correlation")
Direct_CORR2 <- matrix(c(matrix(nci60_gdsc_drugs[[2]],ncol=1),Direct_CORR2),ncol=3)
colnames(Direct_CORR2)<- c("Drug Name","Pearson Correlation","Spearman Correlation")
Direct_CORR3 <- matrix(c(matrix(nci60_gdsc_drugs[[2]],ncol=1),Direct_CORR3),ncol=3)
colnames(Direct_CORR3)<- c("Drug Name","Pearson Correlation","Spearman Correlation")
fwrite(data.table(Direct_CORR3),file="Direct_CORR.csv", 
       showProgress = getOption("datatable.showProgress"),
       verbose = getOption("datatable.verbose"),
       sep = ",", col.names = T, row.names = T)