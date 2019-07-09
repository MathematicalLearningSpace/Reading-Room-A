#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#

#--------------------Work In Progress-------------------------------------#

#------------------------------R API----------------------------------#
library(gdata);library(bio3d);library(igraph);library(sna);library(ips);
library(phangorn);library(proteomics)
library(dcGOR);library(MDplot);library(UniProt.ws);
library(circlize);library(BioPhysConnectoR);library(protr)
library(seqinr);library(Biostrings);library(Peptides);
library(PearsonDS);library(xtable)
library(rcdk);library(BioMedR);library(ChemmineR);
library(Matrix):library(fingerprint)
library(igraph);library(readr);library(leaps);library(caret); library(GA) 
library(ggplot2);library(kohonen);library(pROC)
#------Scientific  Visualization-----#
library(corrplot);library(plot3D);library(scatterplot3d);library(rgl)

#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#
W<-data.frame();X<-data.frame();Y<-data.frame();Z<-data.frame();
#------------Initial Value-------------------#
#------------Parameter Model-----------------#

#---------------------------------------------------------------------#
#------------------------------Functions------------------------------#
#---------------------------------------------------------------------#
Review.Notes<-function(X)
 {
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
 
 
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 
 }
test.Review.Notes.1<-Review.Notes.1("1")
test.Review.Notes.1
#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#
Ribosome.Model.1<-function(X,Visualization=FALSE)
 {
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
 Table.4.df<-as.data.frame(Table.4)
 Table.5.df<-as.data.frame(Table.5);Table.6.df<-as.data.frame(Table.6);Table.7.df<-as.data.frame(Table.7);Table.8.df<-as.data.frame(Table.8)
 Table.9.df<-as.data.frame(Table.9);Table.10.df<-as.data.frame(Table.10)
 
  setwd("Ribosome Model")
  protein.model.files.pdb<-list.files();protein.model.list<-list();K<-length(protein.model.files.pdb)
  #----------------------------------------------------------------------------------------------------------------------------#
  #---------------------------------------PDB Files----------------------------------------------------------------------------#
  #----------------------------------------------------------------------------------------------------------------------------#
  for(i in 1:k){protein.model.list[[i]]<-read.pdb(protein.model.files.pdb[i])}
  bio.Unit.1 <- biounit(protein.model.list[[3]]);names(bio.Unit.1)	
  #---------------------------------------#
  #-------Chain Sequence Modeling---------#
  #---------------------------------------#
  
  #-----------------------------------------------------#
  #-----------------------Residue Selection-------------#
  #-----------------------------------------------------#
 
  #---------------------------------------#
  #-------Network Dependency Modeling-----#
  #---------------------------------------#
  Protein.X.modes<-list()
  for(i in 1:length(protein.model.list))
  {
    Protein.X.modes[[i]] <- nma(protein.model.list[[i]])
  }
  #---------------------------------------#
  #-------Deformation energies------------# 
  #---------------------------------------#	

  #---------------------------------------#
  #----------Fluctuations-----------------# 
  #---------------------------------------#
	
  #---------------------------------------#
  #-------Torsion Analysis By Chain-------#
  #---------------------------------------#

  #---------------------------------------#
  #---cross-correlations------------------#
  #---------------------------------------#
  
  #---------------------------------------------#
  #---------------Graph Theory Analysis---------#
#---------------------------------------------#
 if(Visualization)
	 {
	 png(file = stringr::str_c('Figures//Example_',1,'_Figure_',1,'.png'))
	 Figure.1<-plot(Protein.X.modes[[1]], spread=TRUE)
	 dev.off()
 	}
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 
 }
test.Ribosome.Model.1<-Ribosome.Model.1("1",TRUE)
test.Ribosome.Model.1

Ribosome.Model.2<-function(X)
 {
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
 
 
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 
 }
test.Ribosome.Model.2<-Ribosome.Model.2("1")
test.Ribosome.Model.2

#---------Signal Model-------------#
Signal.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Signal.Model.1<-Signal.Model.1("1")
test.Signal.Model.1

#---------Spline Model-------------#
Spline.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Spline.Model.1<-Spline.Model.1("1")
test.Spline.Model.1

#---Nonlinear Tridiagonal mRNA translation Ribosome Model-------------#

#--------------------Solutions--------------------#
Solution.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
   return(output)
 }
test.Solution.Model.1<-Solution.Model.1("1")
test.Solution.Model.1
#---------------------------------------------------------------------#
#------------------------------Analysis-------------------------------#
#---------------------------------------------------------------------#
Analysis.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Analysis.Model.1<-Analysis.Model.1("1")
test.Analysis.Model.1


#----------------------------------Optimization-----------------------#
Optimization.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Optimization.Model.1<-Optimization.Model.1("1")
test.Optimization.Model.1
#---------------------------------------------------------------------#
#------------------------------Tables---------------------------------#
#---------------------------------------------------------------------#

#----------Table 1-------------#
#----------Table 2-------------#
#----------Table 3-------------#
#----------Table 4-------------#
#----------Table 1------#
Table.1.TeX<-xtable::xtable(Table.1.df)
#----------Table 2------#
Table.2.TeX<-xtable::xtable(Table.2.df)
#----------Table 3------#
Table.3.TeX<-xtable::xtable(Table.3.df)
#----------Table 4------#
Table.4.TeX<-xtable::xtable(Table.4.df)
#---------------------------------------------------------------------#
#------------------------------Figures--------------------------------#
#---------------------------------------------------------------------#

#---Figure 1 for Classroom Presentation------#
#---Figure 2 for Classroom Presentation------#
#---Figure 3 for Classroom Presentation------#
#---Figure 4 for Classroom Presentation------#
png(file = stringr::str_c('Figures//Example_',1,'_Figure_',1,'.png'))
op <- par(mfrow = c(2,2),mar=c(3,3,3,3))
hist(W, main="Title 1",xlab="X Value")
text(4, 9, expression(hat(theta) == (W^t))
legend("topright", legend = paste(seq(1:7),LETTERS[1:7]),lty = 1, cex = .8, y.intersp = 1)
hist(X, main="Title 2",xlab="Note Value")
legend("topright", legend = paste(seq(1:7),LETTERS[1:7]),lty = 1, cex = .8, y.intersp = 1)
hist(Y, main="Title 3",xlab="Note Value")
legend("topright", legend = paste(seq(1:7),LETTERS[1:7]),lty = 1, cex = .8, y.intersp = 1)
hist(Z, main="Title 4",xlab="Note Value")
legend("topright", legend = paste(seq(1:7),LETTERS[1:7]),lty = 1, cex = .8, y.intersp = 1)
par(op)
dev.off()


#-----------------------------References------------------------------#






