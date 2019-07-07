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

library(easyPubMed);library(readr);library(CHNOSZ);
library(stringr);library(seqinr);library(seqLogo);library(msa);library(ape);
library(odseq);library(rphast)
library(plyr)

#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#

#--------------------Review Notes------------------------#
#-----Signal Transduction Networks-----------------------#
A.xml<-fetch_pubmed_data(get_pubmed_ids("Ras signaling pathway")) 
B.xml<-fetch_pubmed_data(get_pubmed_ids("PI3K-Akt signaling pathway"))  
C.xml<-fetch_pubmed_data(get_pubmed_ids("Hippo signaling pathway")) 
D.xml<-fetch_pubmed_data(get_pubmed_ids("Wnt signaling pathway")) 
E.xml<-fetch_pubmed_data(get_pubmed_ids("p53 signaling pathway"))  
F.xml<-fetch_pubmed_data(get_pubmed_ids("Thyroid hormone signaling pathway")) 
G.xml<-fetch_pubmed_data(get_pubmed_ids("FoxO signaling pathway"))
H.xml<-fetch_pubmed_data(get_pubmed_ids("mTOR signaling pathway")) 
I.xml<-fetch_pubmed_data(get_pubmed_ids("ErbB signaling pathway"))  
J.xml<-fetch_pubmed_data(get_pubmed_ids("MAPK signaling pathway"))  
K.xml<-fetch_pubmed_data(get_pubmed_ids("Hippo signaling pathway"))  
L.xml<-fetch_pubmed_data(get_pubmed_ids("Apoptosis"))  
M.xml<-fetch_pubmed_data(get_pubmed_ids("Cell cycle")) 
N.xml<-fetch_pubmed_data(get_pubmed_ids("Rap1 signaling pathway"))  
O.xml<-fetch_pubmed_data(get_pubmed_ids("Chemokine signaling pathway"))
p.xml<-fetch_pubmed_data(get_pubmed_ids("Phosphatidylinositol signaling system"))
Q.xml<-fetch_pubmed_data(get_pubmed_ids("T cell receptor signaling pathway"))
R.xml<-fetch_pubmed_data(get_pubmed_ids("Estrogen signaling pathway"))
S.xml<-fetch_pubmed_data(get_pubmed_ids("VEGF signaling pathway"))

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

#---------------Signal Transduction------------#
Signal.Transduction.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Signal.Transduction.Model.1<-Signal.Transduction.Model.1("1")
test.Signal.Transduction.Model.1
#---------------QSAR Model I-------------------#
QSAR.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.QSAR.Model.1<-QSAR.Model.1("1")
test.QSAR.Model.1
#---------------QSAR Model II------------------#
QSAR.Model.2<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.QSAR.Model.2<-QSAR.Model.1("1")
test.QSAR.Model.2
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

#------------Table 1-----------------#
#------------Table 2-----------------#
#------------Table 3-----------------#
#------------Table 4-----------------#
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

#------------Figure 1-----------------#
#------------Figure 2-----------------#
#------------Figure 3-----------------#
#------------Figure 4-----------------#
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
