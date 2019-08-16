#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#

#-----Work In Progress----------August 2019---------------------------#

#------------------------------R API----------------------------------#
library(gdata);library(bio3d);library(igraph);library(sna);library(ips);
library(phangorn);library(proteomics)
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
W<-data.frame();X<-data.frame();Y<-data.frame();Z<-data.frame();

Gastric.Cancer.Genes<-c("CDX2", "MUC2", "REG4", "CDH17", "MDR1", "SHH","p53", "p21", 
                        "BAX", "p48", "GADD45", "BAK", "POLK","Retinoic.Acid", "RAR.Beta", "RXR",
                        "DV1", "GSK.3Beta", "Beta.Catenin", "Axin", "APC", "CK1.alpha","GBP","EFG", 
                        "ERBB2", "SHC", "GRB2", "SOS", "RAS", "RAF", "MEK", "ERK.1","PI3K", "PIP3", 
                        "AKT", "mTOR", "p53", "S6K", "BCL2","TGF.Beta","TGF.BetaRI","TGF.BetaRII", "SMAD.2", "SMAD.4", "p15", "p21",
                        "HGF", "c.MET", "GRB2", "SOS", "RAS", "RAF", "MEK", "ERK.1",
                        "FGF", "FGFR2", "GAB1", "PI3K", "PIP3", "AKT", "mTOR", "GSK.3Beta")

Cancer.Models.<-c("Nutrition Model","Cancer Model","Stomach Model","Intestine Model","Digestion Model","Protein Model","Carbohydrates Model","Fats Model",
                  "Metabolism Model","Signal Transduction Model","Ribosome Model","Chaperonin Model","Proteasome Model")

A.Drug.Resistance<-c(CDX2, MUC2, REG4, CDH17, MDR1, SHH)
B.Genomic.Instability<-c(p53, p21, BAX, p48, GADD45, BAK, POLK)
C.Tumor.Progression<-c(Retinoic.Acid, RAR.Beta, RXR)
D.Intestinal.Metaplasia<-c(DV1, GSK.3Beta, Beta.Catenin, Axin,APC, CK1.alpha,GBP)
E.Dysplasia.Path.1<-c(EFG, ERBB2, SHC, GRB2, SOS, RAS, RAF, MEK, ERK.1)
F.Dysplasia.Path.2<-c(PI3K, PIP3, AKT, mTOR, p53, S6K, BCL2)
G.Dysplasia.Path.3<-c(TGF.Beta,TGF.BetaRI,TGF.BetaRII, SMAD.2, SMAD.4, p15, p21)
H.Normal.Gastic.Muscosa.1<-c(HGF, c.MET, GRB2, SOS, RAS, RAF, MEK, ERK.1)
I.Normal.Gastic.Muscosa.Survival.Path.1<-c(FGF, FGFR2, GAB1, PI3K, PIP3, AKT, mTOR, GSK.3Beta)



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

#--------Normal Mode Analysis-----------------#
NMA.Model.1<-function(X,visualization=FALSE)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  if(visualization)
   {
   #------------------------Correlation Matrix ------------------------------------
    #png(file = stringr::str_join("Figures//1//Example_",1,"_Figure_","A",".png"))
    #car::scatterplotMatrix(Feature.Matrix.1[,k1:k2])
    #dev.off()
   }
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.NMA.Model.1<-NMA.Model.1("1")
test.NMA.Model.1
#--------Compound Discovery Design------------#
CCD.Model.1<-function(X,visualization=FALSE)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  if(visualization)
   {
   #------------------------Correlation Matrix ------------------------------------
    #png(file = stringr::str_join("Figures//1//Example_",1,"_Figure_","A",".png"))
    #car::scatterplotMatrix(Feature.Matrix.1[,k1:k2])
    #dev.off()
   }
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.CCD.Model.1<-CCD.Model.1("1")
test.CCD.Model.1
#--------QSAR Models--------------------------#
QSAR.Model.1<-function(X,visualization=FALSE)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  Feature.Matrix.1<-X 
 
  if(visualization)
   {
   #------------------------Correlation Matrix ------------------------------------
    #png(file = stringr::str_join("Figures//1//Example_",1,"_Figure_","A",".png"))
    #car::scatterplotMatrix(Feature.Matrix.1[,k1:k2])
    #dev.off()
   }
 
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.QSAR.Model.1<-QSAR.Model.1("1")
test.QSAR.Model.1
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

#----------------------------------Optimization--------------------------------------#
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

#---------------------------------------------------------------------#
#------------------------------Figures--------------------------------#
#---------------------------------------------------------------------#

#------------Figure 1---------------------#
#------------Figure 2---------------------#
#------------Figure 3---------------------#
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
