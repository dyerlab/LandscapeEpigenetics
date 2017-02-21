#msap analysis

rm(list=ls())
load('genos.hpa.rda')
load('genos.msp.rda')

data <- paste(c("","","",names(genos.hpa)[4:450]),collapse=",")


for( i in 1:nrow(genos.hpa)){
  pop <- as.character(genos.hpa$Population[i])
  ind <- as.character(genos.hpa$ID[i])
  
  hpa <- as.numeric( genos.hpa[i,4:450])
  msp <- as.numeric( genos.msp[i,4:450])
  

  hpa.rep <- paste(c(pop,ind,"HPA",as.character(hpa)),collapse=",")
  msp.rep <- paste(c(pop,ind,"MSP",as.character(msp)),collapse=",") 
  
  data <- c( data, hpa.rep, msp.rep )
  
  
}

write.csv(data,file="data.csv", quote=FALSE, row.names = FALSE)


msap("data.csv",name="Araptus")