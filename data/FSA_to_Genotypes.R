#' This function takes the FSA trace files and create a data.frame 
#'   of genotypes.
library( tidyverse )

# Check to see if there is a binner package and if not, install it.
if( !require(binner) ){
  if( !require(devtools)){
    install.packages("devtools")
  }
  devtools::install_github("plantarum/binner")
}
library(binner)

filz <-  list.files(path="./data", pattern=".rda")

###########   MSP DATA Generation #############################################
#                                                                             #
# Check to see if there is no msp data present, and if so, make it.           #
#                                                                             #
###############################################################################

if( !("msp_fsa.rda" %in% filz)) {
  # Read in raw files, load them in and save as a data.frame 
  msp_files <- list.files( path="data/traces/MSP-raw", pattern="*.fsa", full.names=TRUE)
  msp_fsa <- readFSA( msp_files, dye=c("FAM"), bin.width = 1, min.peak.height = 50, CORES=4)
  
  # Three sizing errors were reported
  drop <- c("M161.3","M168.3","M89.4")
  idx <- which( names(msp_fsa$ep) %in% drop ) 
  msp_fsa <- fsaDrop( msp_fsa, idx )
  
  # Save for output
  save(msp_fsa, file="data/msp_fsa.rda")
} else {
  load("data/msp_fsa.rda")
}


# Normalize the traces and then set up a peak table and calculate default AFLP bin sizes.
msp_fsa <- fsaNormalize(msp_fsa)
pts.msp <- fsa2PeakTab( msp_fsa , dye="FAM")
bins <- fsaRGbin(pt, mnbin=1, mxbin=1.5, verbose=FALSE )



###########   HPA DATA Generation #############################################
#                                                                             #
# Check to see if there is no hpa data present, and if so, make it.           #
#                                                                             #
###############################################################################

if( !("hpa_fsa.rda" %in% filz) ){
  hpa_files <- list.files( path="data/traces/HPA-raw/", pattern="*.fsa", full.names=TRUE )
  hpa_fsa <- readFSA( hpa_files, dye=c("FAM"), bin.width=1, min.peak.height = 50, CORES=4)
  save(hpa_fsa, file="data/hpa_fsa.rda")
} else {
  load("data/hpa_fsa.rda")
}
hpa_fsa <- fsaNormalize(hpa_fsa)
pts.hpa <- fsa2PeakTab( hpa_fsa, dye="FAM" )



## Now turn the pts data.frame into set of gentypes based upon pts and bins
#   from binner
ID = sort(unique( pts$sample.name) ) 
m <- matrix(0,nrow=length(ID),ncol=nrow(bins))
colnames(m) <- paste("loc",seq(1,nrow(bins)),sep="-")
rownames(m) <- ID
df <- data.frame( ID )
for( i in 1:nrow(bins)){
  bands <-  pts$sample.name[ pts$bp >= bins[i,1] & 
                               pts$bp < bins[i,2] ]
  col <- rep(0,length(ID))
  col[ ID %in% bands ] <- 1
  df[[ paste("loc",i,sep="-")]] <- locus( col, type="aflp")
}





