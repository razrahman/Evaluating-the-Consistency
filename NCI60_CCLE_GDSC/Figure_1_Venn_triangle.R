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

######## Cells Lines #####
  g1 = draw.triple.venn(
    area1 = 60, # 10 O'clock - NCI60
    area2 = 501,  # 2 O'clock - CCLE
    area3 = 1027, # 6 O'clock - GDSC
    n12 = 30, # 12 O'clock
    n23 = 378, # 4 O'clock
    n13 = 59,# 8 O'clock
    n123 = 30, # Center
    category = c("NCI60", "CCLE", "GDSC"),
    fill = c("cyan", "magenta", "yellow"),
    lty = "blank",
    cex = rep(7,7), # Size of the areas # A vector (length 7) giving the size of the areas' labels
    cat.cex = rep(4,3), # Size of the names
    # Fonts: https://cran.r-project.org/web/packages/extrafont/README.html
    cat.fontfamily = rep("Cambria-Bold", 3),
    fontfamily = rep("Lucida Console",7),
    ext.text = TRUE,
    cat.col = c("cornflower blue", "purple4", "salmon4"),
    scaled=TRUE, ind = FALSE,
    alpha = rep(0.5,3) )
  #setwd(mainPath)
  png(file = "cells_Venn_all_datasets_color2.png", width = 1024, height = 1024, units = "px")
  
  # vennTitle <- "Gene intersection diagram of the Gene Expression data in the datasets: NCI60,CCLE and GDSC"
  grid.draw(g1)
  # grid.arrange(gTree(children=g), main= textGrob(paste(vennTitle,sep=""),
  #                                               gp=gpar(fontsize=16,font=2)) )
  dev.off()
  grid.newpage()
 
  #### drugs #########
  # including the gdsc synonyms
  g2 = draw.triple.venn(
    area1 = 34496, # NCI60 -> 34,496 (From GI50) - 10 O'clock 52,204
    area2 = 24,  # CCLE ->  - 24 - 2 O'clock
    area3 = 250, # GDSC -> 250 - 6 O'clock 
    n12 = 10, # NCI/CCLE ->  10 - 12 O'clock
    n23 = 16, # CCLE/GDSC -> 9 - 4 O'clock 
    n13 = 132, # NCI/GDSC -> 47 - 8 O'clock
    n123 = 9, # 6 - Center
    category = c("NCI60", "CCLE", "GDSC"),
    fill = c("cyan", "magenta", "yellow"),
    lty = "blank",
    cex = rep(7,7), # Size of the areas
    cat.cex = rep(4,3), # Size of the names
    # Fonts: https://cran.r-project.org/web/packages/extrafont/README.html
    cat.fontfamily = rep("Cambria-Bold", 3),
    fontfamily = rep("Lucida Console",7),
    ext.text = TRUE,
    cat.col = c("cornflower blue", "purple4", "salmon4"),
    scaled=TRUE, ind = FALSE,
    alpha = rep(0.5,3) )
  #setwd(mainPath)
  png(file = "drugs_Venn_all_datasets_color2.png", width = 1024, height = 1024, units = "px")
  # To display Commas
  idx <- sapply(g2, function(i) grepl("text", i$name))
  for(i in 1:7){
    g2[idx][[i]]$label <- 
      format(as.numeric(g2[idx][[i]]$label), big.mark=",", scientific=FALSE)}
  # vennTitle <- "Gene intersection diagram of the Gene Expression data in the datasets: NCI60,CCLE and GDSC"
  grid.draw(g2)
  # grid.arrange(gTree(children=g), main= textGrob(paste(vennTitle,sep=""),
  #                                              gp=gpar(fontsize=16,font=2)) )
  dev.off()
  grid.newpage()
  
  ##### Gene #####
  
    g3 = draw.triple.venn(
      area1 = 23686, # NCI60 -> 23,686 10 O'clock
      area2 = 54336,  # CCLE -> 54,336 genes - 2 O'clock
      area3 = 56697, # GDSC -> 56,697 genes - 6 O'clock 
      n12 = 21197, # 21,197 - 12 O'clock
      n23 = 16996, # 16,996 - 4 O'clock
      n13 = 15966,# 15,996 - 8 O'clock
      n123 = 15579, # 15,579 Center
      category = c("NCI60", "CCLE", "GDSC"),
      fill = c("cyan", "magenta", "yellow"),
      lty = "blank",
      cex = rep(7,7), # Size of the areas
      cat.cex = rep(4,3), # Size of the names
      # Fonts: https://cran.r-project.org/web/packages/extrafont/README.html
      cat.fontfamily = rep("Cambria-Bold", 3),
      fontfamily = rep("Lucida Console",7),
      ext.text = TRUE,
      cat.col = c("cornflower blue", "purple4", "salmon4"),
      scaled=TRUE, ind = FALSE,
      alpha = rep(0.5,3) )
    #setwd(mainPath)
    png(file = "genes_Venn_all_datasets_color2.png", width = 1024, height = 1024, units = "px")
    # To display Commas
    idx <- sapply(g3, function(i) grepl("text", i$name))
    for(i in 1:7){
      g3[idx][[i]]$label <- 
        format(as.numeric(g3[idx][[i]]$label), big.mark=",", scientific=FALSE)}
    vennTitle <- "Gene intersection diagram of the Gene Expression data in the datasets: NCI60,CCLE and GDSC"
    grid.draw(g3)
    # grid.arrange(gTree(children=g), main= textGrob(paste(vennTitle,sep=""),
    #                                              gp=gpar(fontsize=16,font=2)) )
    dev.off()
    grid.newpage()
  
  