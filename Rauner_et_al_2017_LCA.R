# ------------------------------------------------------------------------------
# Program Name: Rauner_et_al_2017_LCA.R
# Author(s): Sebastian Rauner 
# Date Last Updated: Mar 20, 2017 
# script developed to perform the LCA in the publication Sebastian Rauner et al 2017 Environ. Res. Lett. https://doi.org/10.1088/1748-9326/aa914d
# Program Purpose: This script reads the spares ecoinvent matricies and perfomes the LCA for the technologies in the lists technology_list_variable.csv and technology_list_fixed.csv
# additonally it perfomes the spatial contribution analysis
# Input Files: A (process), B (emission) and C (impact) matrices & the technologies lists (demand vectors)
# Output Files: s, g, h matrices
# Notes: 
#setwd("Y:/Home/rauner/LCA/ecoinvent_32_cutoff_universal_18_05_16/cutoff")

library(Matrix)
library(miscTools)

# chose the ecoinvent model (cutoff, APOS, consequential)
set_ecoinvent_model <- 'cutoff'

setwd(paste0("Y:/Home/rauner/LCA/ecoinvent_32_cutoff_universal_18_05_16/",set_ecoinvent_model)

#read in row, columns and A,B,C matrices
#watch out for commas and semicolons
A.technology.row <- as.data.frame(read.csv("ie_index_ordered.csv",header=T,sep=",", dec = "."))
A.technology.column <- as.data.frame(read.csv("ie_index_ordered.csv",header=T,sep=",", dec = "."))
A.technology.readin <- as.data.frame(read.csv("A_confidential_hidden.csv",header=T,sep=";", dec = "."))
B.intervention.row <- as.data.frame(read.csv("ee_index_ordered.csv",header=T,sep=",", dec = "."))
B.intervention.column <- as.data.frame(read.csv("ie_index_ordered.csv",header=T,sep=",", dec = "."))
B.intervention.readin <- as.data.frame(read.csv("B_confidential_hidden.csv",header=T,sep=";", dec = "."))
C.characterization.row <- as.data.frame(read.csv("LCIA_index_ordered.csv",header=T,sep=",", dec = "."))
C.characterization.column <- as.data.frame(read.csv("ee_index_ordered.csv",header=T,sep=",", dec = "."))
C.characterization.readin <- as.data.frame(read.csv("C.csv",header=T,sep=";", dec = "."))
technology.list.variable <- as.data.frame(read.csv("technology_list_variable.csv",header=T,sep=",", dec = "."))
technology.list.fixed <- as.data.frame(read.csv("technology_list_fixed.csv",header=T,sep=",", dec = "."))

#assign rownames
row.names(technology.list.fixed)    <- technology.list.fixed[,"fixed"]
row.names(technology.list.variable) <- technology.list.variable[,"variable"]

#generate sparse matrices, no units
A.technology.final.sparse       <- sparseMatrix(i = A.technology.readin$row, j = A.technology.readin$column,
                                                x = A.technology.readin$coefficient, dim=c(nrow(A.technology.row),nrow(A.technology.column)),dimnames = list(A.technology.row$label,A.technology.column$label),index1=FALSE)

B.intervention.final.sparse     <- sparseMatrix(i = B.intervention.readin$row, j = B.intervention.readin$column,
                                                x = B.intervention.readin$coefficient, dim=c(nrow(B.intervention.row),nrow(B.intervention.column)),dimnames = list(B.intervention.row$label,B.intervention.column$label),index1=FALSE)

C.characterization.final.sparse <- sparseMatrix(i = C.characterization.readin$row, j = C.characterization.readin$column,
                                                x = C.characterization.readin$coefficient, dim=c(nrow(C.characterization.row),nrow(C.characterization.column)),dimnames = list(C.characterization.row$label,C.characterization.column$label),index1=FALSE)



#A matrix of technology list -> which and how many of one process goes into another process
A.technology.final.technology.list.variable <- A.technology.final.sparse[,colnames(A.technology.final.sparse)%in%technology.list.variable[,1]] 
A.technology.final.technology.list.fixed    <- A.technology.final.sparse[,colnames(A.technology.final.sparse)%in%technology.list.fixed[,1]] 

# see book "The Computational Structure of Life Cycle Assessment" page 219
# calculation of the inverse of A, then the scaling vector s and finally the inventory vector g.

# the demand vector is constructed for all variable and fixed parts of the relevant processes
# finally the impact vector h can be calculated with h=C*g (page 168)

# s=A^-1*f, g=B*s, h=C*g



#inverse of A
A.technology.final.inverse.sparse <- solve(A.technology.final.sparse)

dimnames(A.technology.final.inverse.sparse)  <- dimnames(A.technology.final.sparse)


#construct demand vector of the single process
#demand vector is the columns of A.technology.final.technology.list
#------------------------------demand vector f---------------------------------------------------------------
f_process_only <- Matrix(0, nrow(A.technology.final.technology.list.variable), ncol(A.technology.final.technology.list.variable), sparse=TRUE)
dimnames(f_process_only) <- dimnames(A.technology.final.technology.list.variable)
f_variable <- f_process_only
dimnames(f_variable) <- dimnames(A.technology.final.technology.list.variable)
f_fixed_unit   <- Matrix(0, nrow(A.technology.final.technology.list.fixed), ncol(A.technology.final.technology.list.fixed), sparse=TRUE)
dimnames(f_fixed_unit) <- dimnames(A.technology.final.technology.list.fixed)
f_fixed_share  <- f_fixed_unit 
f_fixed_kW     <- f_fixed_unit
f_elec_mix     <- Matrix(0, nrow(A.technology.final.technology.list.variable), 1, sparse=TRUE)
rownames(f_elec_mix)     <- rownames(A.technology.final.technology.list.variable)



x<-1
while(x <= ncol(A.technology.final.technology.list.variable))
{ 
  y<-1
  while(y <= nrow(A.technology.final.technology.list.variable))
  {
  if(rownames(A.technology.final.technology.list.variable)[y] == colnames(A.technology.final.technology.list.variable)[x])
  {f_process_only[y,x] <- A.technology.final.technology.list.variable[y,x]}
  
  y<-y+1
  
  }
  x<-x+1
  
}

x<-1
while(x <= ncol(A.technology.final.technology.list.fixed))
{ 
  y<-1
  while(y <= nrow(A.technology.final.technology.list.fixed))
  {
    if(rownames(A.technology.final.technology.list.fixed)[y] == colnames(A.technology.final.technology.list.fixed)[x])
    { f_fixed_unit[y,x]  <- 1
      f_fixed_kW[y,x]    <- A.technology.final.technology.list.fixed[y,x]/technology.list.fixed[colnames(A.technology.final.technology.list.fixed)[x],"amount_fixed"]
      #variable process to fixed process mapping with typ
      process_variable<- technology.list.fixed[rownames(A.technology.final.technology.list.fixed)[y],"Typ"]
      z<-1
      while(z <= nrow(technology.list.variable))
      {
      if(toString(technology.list.variable[z,"Typ"])== toString(process_variable)) {
        variable_process_col <-  toString(technology.list.variable[z,"variable"])
        
      }
      z<-z+1
    }
      
     f_fixed_share[y,x] <- A.technology.final.sparse[rownames(A.technology.final.technology.list.fixed)[y],variable_process_col]

    }
    
    y<-y+1
    
  }
  x<-x+1
  
}

#-------------------------------------------Mapping with the GESOP technology names


#insert one row at position 2
#go through the columns and fill with mapping process
# if name = name_mapping then mapping

x<-1
while(x <= ncol(A.technology.final.technology.list.variable))
{
  y<-1
  while(y <= nrow(A.technology.final.technology.list.fixed))
  {
    if(colnames(f_process_only)[x] == toString(technology.list.variable[y,"variable"])) {colnames(f_process_only)[x]<-toString(technology.list.variable[y,"Typ"])}
    if(colnames(f_variable)[x] == toString(technology.list.variable[y,"variable"])) {colnames(f_process_only)[x]<-toString(technology.list.variable[y,"Typ"])}
    
    y<-y+1
  }
  x<-x+1  
}

x<-1
while(x <= ncol(A.technology.final.technology.list.fixed))
{
  y<-1
  while(y <= nrow(technology.list.fixed))
  {
    
    
    if(colnames(f_fixed_unit)[x] == toString(technology.list.fixed[y,"fixed"])) {colnames(f_fixed_unit)[x]<-toString(technology.list.fixed[y,"Typ"])}       
    if(colnames(f_fixed_share)[x] == toString(technology.list.fixed[y,"fixed"])) {colnames(f_fixed_share)[x]<-toString(technology.list.fixed[y,"Typ"])}  
    if(colnames(f_fixed_kW)[x] == toString(technology.list.fixed[y,"fixed"])) {colnames(f_fixed_kW)[x]<-toString(technology.list.fixed[y,"Typ"])}
       
    
    y<-y+1
  }
  
  x<-x+1
  
}

#condense the fixed demand vectors

f_fixed_unit<-f_fixed_unit %*% sapply(unique(colnames(f_fixed_unit)),"==", colnames(f_fixed_unit))
f_fixed_share<-f_fixed_share %*% sapply(unique(colnames(f_fixed_share)),"==", colnames(f_fixed_share))
f_fixed_kW<-f_fixed_kW %*% sapply(unique(colnames(f_fixed_kW)),"==", colnames(f_fixed_kW))

colnames(f_variable) <-  colnames(f_process_only)

x<-1
while(x <= ncol(f_variable))
{
  y<-1
  while(y <= ncol(f_variable))
  {
    
  
    if(colnames(f_variable)[x] == colnames(f_fixed_share)[y]) {f_variable[,x] <- f_process_only[,x] - abs(f_fixed_share[,y])}       
    
    y<-y+1
  }
    
    
  x<-x+1
  
}

# write the different demand vectors

#f_process_only#only the process with a 1
#f_fixed_unit  #only the unit with a 1
#f_fixed_share #share of a unit in 1kWh
#f_fixed_kW    # one kW installed capacity
#f_variable    # f_process_only - abs(f_fixed_share)


write.csv(as.matrix(f_process_only), file="f_process_only.csv")
write.csv(as.matrix(f_fixed_unit), file="f_fixed_unit.csv")
write.csv(as.matrix(f_fixed_share), file="f_fixed_share.csv")
write.csv(as.matrix(f_fixed_kW), file="f_fixed_kW.csv")
write.csv(as.matrix(f_variable), file="f_variable.csv")

#------------------------------INITIALIZE S,G and H----------------------------------------------------------------------

s_process_only<- f_process_only
g_process_only<-Matrix(0, nrow(B.intervention.final.sparse), ncol(f_process_only), sparse=TRUE)
colnames(g_process_only)  <- colnames(f_process_only)
rownames(g_process_only) <- rownames(B.intervention.final.sparse)
h_process_only<-Matrix(0, nrow(C.characterization.final.sparse), ncol(f_process_only), sparse=TRUE)
colnames(h_process_only)  <- colnames(f_process_only)
rownames(h_process_only) <- rownames(C.characterization.final.sparse)

s_variable<-f_variable
g_variable<-Matrix(0, nrow(B.intervention.final.sparse), ncol(f_process_only), sparse=TRUE)
colnames(g_variable)  <- colnames(f_process_only)
rownames(g_variable) <- rownames(B.intervention.final.sparse)
h_variable<-Matrix(0, nrow(C.characterization.final.sparse), ncol(f_process_only), sparse=TRUE)
colnames(h_variable)  <- colnames(f_process_only)
rownames(h_variable) <- rownames(C.characterization.final.sparse)

s_fixed_unit<-f_fixed_unit
g_fixed_unit<-Matrix(0, nrow(B.intervention.final.sparse), ncol(f_fixed_unit), sparse=TRUE)
colnames(g_fixed_unit)  <- colnames(f_fixed_unit)
rownames(g_fixed_unit) <- rownames(B.intervention.final.sparse)
h_fixed_unit<-Matrix(0, nrow(C.characterization.final.sparse), ncol(f_fixed_unit), sparse=TRUE)
colnames(h_fixed_unit)  <- colnames(f_fixed_unit)
rownames(h_fixed_unit) <- rownames(C.characterization.final.sparse)

s_fixed_share<-f_fixed_share
g_fixed_share<-Matrix(0, nrow(B.intervention.final.sparse), ncol(f_fixed_share), sparse=TRUE)
colnames(g_fixed_share)  <- colnames(f_fixed_share)
rownames(g_fixed_share) <- rownames(B.intervention.final.sparse)
h_fixed_share<-Matrix(0, nrow(C.characterization.final.sparse), ncol(f_fixed_share), sparse=TRUE)
colnames(h_fixed_share)  <- colnames(f_fixed_share)
rownames(h_fixed_share) <- rownames(C.characterization.final.sparse)

s_fixed_kW<-f_fixed_kW
g_fixed_kW<-Matrix(0, nrow(B.intervention.final.sparse), ncol(f_fixed_kW), sparse=TRUE)
colnames(g_fixed_kW)  <- colnames(f_fixed_kW)
rownames(g_fixed_kW) <- rownames(B.intervention.final.sparse)
h_fixed_kW<-Matrix(0, nrow(C.characterization.final.sparse), ncol(f_fixed_kW), sparse=TRUE)
colnames(h_fixed_kW)  <- colnames(f_fixed_kW)
rownames(h_fixed_kW) <- rownames(C.characterization.final.sparse)

f_elec_mix["electricity, high voltage//[DE] market for electricity, high voltage",1] <-  1
h_elec_mix<-Matrix(0, nrow(C.characterization.final.sparse), ncol(f_elec_mix), sparse=TRUE)

#calculate s,g and h of the above demand vectors, here for process only, below for the rest
x<-1
while(x <= ncol(f_process_only))
{ 
  
  s_process_only[,x]<- A.technology.final.inverse.sparse %*% f_process_only[,x]
  g_process_only[,x]<- B.intervention.final.sparse %*% s_process_only[,x]
  h_process_only[,x]<- C.characterization.final.sparse %*% g_process_only [,x]
  
  s_variable[,x]<- A.technology.final.inverse.sparse %*% f_variable[,x]
  g_variable[,x]<- B.intervention.final.sparse %*% s_variable[,x]
  h_variable[,x]<- C.characterization.final.sparse %*% g_variable[,x]
  
  x<-x+1
  
}

x<-1
while(x <= ncol(f_fixed_unit))
{
  
  s_fixed_unit[,x]<- A.technology.final.inverse.sparse %*% f_fixed_unit[,x]
  g_fixed_unit[,x]<- B.intervention.final.sparse %*% s_fixed_unit[,x]
  h_fixed_unit[,x]<- C.characterization.final.sparse %*% g_fixed_unit [,x]
  
  s_fixed_share[,x]<- A.technology.final.inverse.sparse %*% f_fixed_share[,x]
  g_fixed_share[,x]<- B.intervention.final.sparse %*% s_fixed_share[,x]
  h_fixed_share[,x]<- C.characterization.final.sparse %*% g_fixed_share [,x]
  
  s_fixed_kW[,x]<- A.technology.final.inverse.sparse %*% f_fixed_kW[,x]
  g_fixed_kW[,x]<- B.intervention.final.sparse %*% s_fixed_kW[,x]
  h_fixed_kW[,x]<- C.characterization.final.sparse %*% g_fixed_kW [,x]
  
  

x<-x+1

}

#####################contribution analysis -> see what processes contribute how much to the overall impact score
  s_elec_mix    <- A.technology.final.inverse.sparse %*% f_elec_mix
  ss_elec_mix    <- as.matrix(s_elec_mix)
  sss_elec_mix  <- matrix( rep( t(ss_elec_mix) , nrow(C.characterization.final.sparse) ) , ncol =  ncol(t(ss_elec_mix)) , byrow = TRUE)
  h_elec_mix   <- (C.characterization.final.sparse %*%  B.intervention.final.sparse) * sss_elec_mix

  write.csv(as.matrix(h_elec_mix), file="h_elec_mix.csv")
  
#####################contribution analysis -> see what processes contribute how much to the overall impact score

#assigne units, technology units missing
#s <- [s[,1:1] ]
s_variable <- cbind(as.matrix(A.technology.row$unit),as.matrix(s_variable))
s_fixed_kW <- cbind(as.matrix(A.technology.row$unit),as.matrix(s_fixed_kW))

h_variable <- cbind(as.matrix(C.characterization.row[,3]), as.matrix(h_variable))
h_fixed_kW <- cbind(as.matrix(C.characterization.row[,3]), as.matrix(h_fixed_kW))


#------------------------------write_out of variable and fixed_kW as an input for the model-------------------------------------------------------


write.csv(as.matrix(s_variable), file="s_variable.csv")
write.csv(as.matrix(g_variable), file="g_variable.csv")
write.csv(as.matrix(h_variable), file="h_variable.csv")

write.csv(as.matrix(s_fixed_kW), file="s_fixed_kW.csv")
write.csv(as.matrix(g_fixed_kW), file="g_fixed_kW.csv")
write.csv(as.matrix(h_fixed_kW), file="h_fixed_kW.csv")

