#############################################
#        _                 _       _        #
#     __| |_   _  ___ _ __| | __ _| |__     #
#    / _` | | | |/ _ \ '__| |/ _` | '_ \    #
#   | (_| | |_| |  __/ |  | | (_| | |_) |   #
#    \__,_|\__, |\___|_|  |_|\__,_|_.__/    #
#          |___/                            #
#                                           #
#############################################

#' Routines used to take raw ABI msAFLP data from Chitra's project and 
#'  process them for input into gstudio.
#'  
rm(list=ls())
require(binner) 
require(gWidgets)
require(gWidgetsRGtk2)


if (0) {
  files1 <- list.files(path="./msp/",pattern="*.fsa",full.names = TRUE)
  fsa1 <- readFSA( files1, dye = c("FAM") )
  save(fsa1,file="msp.fsa.rda")
  
  files2 <- list.files(path="./hpa/",pattern="*.fsa",full.names = TRUE)
  fsa2 <- readFSA( files2, dye = c("FAM") )
  save(fsa2,file="hpa.fsa.rda")
} else {
  
  load("fsa1.rda")
  load("fsa2.rda")
}



fsa1 <- fsaNormalize( fsa1 )
fsa2 <- fsaNormalize( fsa2 )

pt <- fsa2PeakTab( fsa1 , dye="FAM")
bins <- fsaRGbin(pt, mnbin=1, mxbin=1.5, verbose=TRUE, start=50, end=495 )

#' use this to visualize the bins
# scanGel(pt,bins) ##Some black/green bands...samples were tagged with FAM dye and have the LIZ sizing control.

#' Lets call alleles now.
pt1 <- fsa2PeakTab(fsa1, dye="FAM")
pt2 <- fsa2PeakTab(fsa2, dye="FAM")
msp <- binSet(pt1,bins,pref="msp")
hpa <- binSet(pt2,bins,pref="hpa")


# make the msp data frame and clean it up
data <- data.frame(msp[,,"alleles"])
samples <- rownames(data)
pop_id <- matrix( unlist(strsplit(substring(samples, first=2), split=".", fixed=T)), ncol=2,byrow=T)
data$Population <- pop_id[,1]
data$ID <- pop_id[,2]
#genos.msp <- data[,c(708,709,1:707)]
#Error in `[.data.frame`(data, , c(708, 709, 1:707)) :undefined columns selected
genos.msp <- data[,c(657,658,1:656)]
#genos.msp$msp.NaN <- genos.msp$msp.NaN.1 <- NULL
genos.msp <- genos.msp[ order(genos.msp$Population, genos.msp$ID),]
save(genos.msp,file="genos.msp1.rda")


# make the hpa data frame and clean it up
data <- data.frame(hpa[,,"alleles"])
samples <- rownames(data)
pop_id <- matrix( unlist(strsplit(substring(samples, first=2), split=".", fixed=T)), ncol=2,byrow=T)
data$Population <- pop_id[,1]
data$ID <- pop_id[,2]
#genos.hpa <- data[,c(708,709,1:707)]
#Error in `[.data.frame`(data, , c(708, 709, 1:707)) : undefined columns selected
genos.hpa <- data[,c(657,658,1:656)]
genos.hpa$hpa.NaN <- genos.hpa$hpa.NaN.1 <- NULL
genos.hpa <- genos.hpa[ order(genos.hpa$Population, genos.hpa$ID),] 
save(genos.hpa, file="genos.hpa1.rda")

