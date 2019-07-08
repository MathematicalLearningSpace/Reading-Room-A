#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#

#--------------------Work In Progress-------------------------------------#


#------------------------------R API----------------------------------#
library(gdata);library(bio3d);library(igraph);library(sna);library(ips);
library(phangorn);library(proteomics);
library(dcGOR);library(MDplot);library(UniProt.ws);
library(circlize);library(BioPhysConnectoR);library(protr)
library(seqinr);library(Biostrings);library(Peptides);
library(PearsonDS);library(xtable)
library(rcdk);library(BioMedR);library(ChemmineR);
library(Matrix):library(fingerprint)
library(readr);library(leaps);library(caret); library(GA) 
library(ggplot2);library(kohonen);library(pROC)
#------Scientific  Visualization-----#
library(corrplot);library(plot3D);library(scatterplot3d);library(rgl)

#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#
Table.1.df<-as.data.frame(Table.1);Table.2.df<-as.data.frame(Table.2);Table.3.df<-as.data.frame(Table.3);
Table.4.df<-as.data.frame(Table.4)Table.5.df<-as.data.frame(Table.5);Table.6.df<-as.data.frame(Table.6);
Table.7.df<-as.data.frame(Table.7);Table.8.df<-as.data.frame(Table.8)Table.9.df<-as.data.frame(Table.9);
Table.10.df<-as.data.frame(Table.10)

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

#---------Translational and Transcriptional Error Models I------------#
Error.Translation.Transcriptional.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Error.Translation.Transcriptional.Model.1<-Error.Translation.Transcriptional.Model.1("1")
test.Error.Translation.Transcriptional.Model.1
#---------Translational and Transcriptional Error Models II-----------#
Error.Translation.Transcriptional.Model.2<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Error.Translation.Transcriptional.Model.2<-Error.Translation.Transcriptional.Model.2("1")
test.Error.Translation.Transcriptional.Model.2
#---------Translational and Transcriptional Error Models III----------#
Error.Translation.Transcriptional.Model.3<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Error.Translation.Transcriptional.Model.3<-Error.Translation.Transcriptional.Model.3("1")
test.Error.Translation.Transcriptional.Model.3
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
#---------------------------------------------------------------------#
#------------------------------Tables---------------------------------#
#---------------------------------------------------------------------#

#-----------Table 1-------------------#
#-----------Table 2-------------------#
#-----------Table 3-------------------#
#-----------Table 4-------------------#

#----------Table 1------#
Table.1.TeX<-xtable::xtable(Table.1.df)
#----------Table 2------#
Table.2.TeX<-xtable::xtable(Table.2.df)
#----------Table 3------#
Table.3.TeX<-xtable::xtable(Table.3.df)
#----------Table 4------#
Table.4.TeX<-xtable::xtable(Table.4.df)
#----------Table 5------#
Table.1.TeX<-xtable::xtable(Table.5.df)
#----------Table 6------#
Table.2.TeX<-xtable::xtable(Table.6.df)
#----------Table 7------#
Table.3.TeX<-xtable::xtable(Table.7.df)
#----------Table 8------#
Table.4.TeX<-xtable::xtable(Table.8.df)

#---------------------------------------------------------------------#
#------------------------------Figures--------------------------------#
#---------------------------------------------------------------------#

#------------Figure 1---------------------#
#------------Figure 2---------------------#
#------------Figure 3---------------------#

png(file = stringr::str_c('Figures//Example_',1,'_Figure_',3,'.png'))
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
