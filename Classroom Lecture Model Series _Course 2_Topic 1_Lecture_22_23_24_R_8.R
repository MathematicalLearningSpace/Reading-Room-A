#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#

#--------------------Work In Progress-------------------------------------#

#------------------------------R API----------------------------------#
library(deSolve);library(ReacTran);library(rootSolve);
library(rcellminer);library(rcellminerData);library(bioCancer);
library(trena);library(ctSGE);library(anamiR) 
library(srnadiff);library(RmiR);library(CancerSubtypes);
library(GRENITS);library(compEpiTools);library(trigger);
library(CVE);library(qgraph);library(igraph)
library(sgnesR);library(rcdk);
library(fda);library(phaseR)
library(pracma);library(GA)
#----------------------------Scientific  Visualization----------------#
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
#-----------------------Review Notes----------------------------------#
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
#------------------------------Functions------------------------------#
#---------------------------------------------------------------------#
Transformation.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Transformation.1<-Transformation.1("1")
test.Transformation.1
#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#

#---------------Digestion I-----------------#
Digestion.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  setwd("Digestion Model")
 
  DGM.1.Name="Salivary secretion";DGM.1<-xml2::read_xml("hsa04970.xml")
  DGM.2.Name="Gastric acid secretion";DGM.2<-xml2::read_xml("hsa04971.xml")
  DGM.3.Name="Pancreatic secretion";DGM.3<-xml2::read_xml("hsa04972.xml")
  DGM.4.Name="Bile secretion";DGM.4<-xml2::read_xml("hsa04976.xml")
  DGM.5.Name="Carbohydrate digestion and absorption";DGM.5<-xml2::read_xml("hsa04973.xml")
  DGM.6.Name="Protein digestion and absorption";DGM.6<-xml2::read_xml("hsa04974.xml")
  DGM.7.Name="Fat digestion and absorption";DGM.7<-xml2::read_xml("hsa04975.xml")
  DGM.8.Name="Cholesterol metabolism";DGM.8<-xml2::read_xml("hsa04979.xml")
  DGM.9.Name="Vitamin digestion and absorption";DGM.9<-xml2::read_xml("hsa04977.xml")
  DGM.10.Name="Mineral absorption";DGM.10<-xml2::read_xml("hsa04978.xml")
 
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Digestion.Model.1<-Digestion.Model.1("1")
test.Digestion.Model.1
#---------------Digestion II----------------#
Digestion.Model.2<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Digestion.Model.2<-Digestion.Model.2("1")
test.Digestion.Model.2
#---------------Sugar-----------------------#
Sugar.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Sugar.Model.1<-Sugar.Model.1("1")
test.Sugar.Model.1
#---------------Insulin---------------------#
Insulin.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Insulin.Model.1<-Insulin.Model.1("1")
test.Insulin.Model.1
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
#----------------------------------Optimization-------------------------#
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

#----------Table 1 Group A---------#

#----------Table 2 Group B---------#

#----------Table 3 Group C---------#

#----------Table 4 Group D---------#
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

#----------Figure 1 Group A---------#

#----------Figure 1 Group B---------#

#----------Figure 1 Group C---------#

#----------Figure 1 Group D---------#
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
