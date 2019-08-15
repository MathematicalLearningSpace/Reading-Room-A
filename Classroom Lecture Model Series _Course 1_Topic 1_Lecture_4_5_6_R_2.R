#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#

#--------------------Work In Progress-------------------------------------#
#------------------------------R API----------------------------------#
library(deSolve);library(ReacTran);library(rootSolve);
library(fda);library(phaseR)
library(pracma);library(GA);library(igraph);
library(Sim.DiffProc);library(bvpSolve);library(odeintr);library(scaRabee)
library(yuima);library(Sim.DiffProc);library(fptdApprox);library(rpgm);
#-------------------Combinatorics-------------------------------------#
library(combinat);library(GeomComb)
#-------------------Numerical Derivatives-----------------------------#
library(Deriv);library(numDeriv)
#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#

A.Drug.Resistance<-c(CDX2, MUC2, REG4, CDH17, MDR1, SHH)
B.Genomic.Instability<-c(p53, p21, BAX, p48, GADD45, BAK, POLK)
C.Tumor.Progression<-c(Retinoic.Acid, RAR.Beta, RXR)
D.Intestinal.Metaplasia<-c(DV1, GSK.3Beta, Beta.Catenin, Axin,APC, CK1.alpha,GBP)
E.Dysplasia.Path.1<-c(EFG, ERBB2, SHC, GRB2, SOS, RAS, RAF, MEK, ERK.1)
F.Dysplasia.Path.2<-c(PI3K, PIP3, AKT, mTOR, p53, S6K, BCL2)
G.Dysplasia.Path.3<-c(TGF.Beta,TGF.BetaRI,TGF.BetaRII, SMAD.2, SMAD.4, p15, p21)
H.Normal.Gastic.Muscosa.1<-c(HGF, c.MET, GRB2, SOS, RAS, RAF, MEK, ERK.1)
I.Normal.Gastic.Muscosa.Survival.Path.1<-c(FGF, FGFR2, GAB1, PI3K, PIP3, AKT, mTOR, GSK.3Beta)

#-----------------------------Parameter Model-------------------------#
params.1<-c(a11=0.1,a12=0.1,a13=0.1,a14=0.1,a15=0.1,a16=0.1,
            a21=0.1,a22=0.1,a23=0.1,a24=0.1,a25=0.1,a26=0.1,
            a31=0.1,a32=0.1,a33=0.1,a34=0.1,a35=0.1,a36=0.1,
            a41=0.1,a42=0.1,a43=0.1,a44=0.1,a45=0.1,a46=0.1,
            a51=0.1,a52=0.1,a53=0.1,a54=0.1,a55=0.1,a56=0.1,
            a61=0.1,a62=0.1,a63=0.1,a64=0.1,a65=0.1,a66=0.1)

params.epsilon.1<-c(epsilon1,epsilon2,epsilon3)


#---------------------------------------------------------------------#
#-----------------------Review Notes----------------------------------#
#---------------------------------------------------------------------#

Review.Notes<-function(X)
 {
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  #--------------------Review Notes------------------------#
  #-----Signal Transduction Networks-----------------------#
#A.xml<-fetch_pubmed_data(get_pubmed_ids("Ras signaling pathway")) 
#B.xml<-fetch_pubmed_data(get_pubmed_ids("PI3K-Akt signaling pathway"))  
#C.xml<-fetch_pubmed_data(get_pubmed_ids("Hippo signaling pathway")) 
#D.xml<-fetch_pubmed_data(get_pubmed_ids("Wnt signaling pathway")) 
#E.xml<-fetch_pubmed_data(get_pubmed_ids("p53 signaling pathway"))  
#F.xml<-fetch_pubmed_data(get_pubmed_ids("Thyroid hormone signaling pathway")) 
#G.xml<-fetch_pubmed_data(get_pubmed_ids("FoxO signaling pathway"))
#H.xml<-fetch_pubmed_data(get_pubmed_ids("mTOR signaling pathway")) 
#I.xml<-fetch_pubmed_data(get_pubmed_ids("ErbB signaling pathway"))  
#J.xml<-fetch_pubmed_data(get_pubmed_ids("MAPK signaling pathway"))  
#K.xml<-fetch_pubmed_data(get_pubmed_ids("Hippo signaling pathway"))  
#L.xml<-fetch_pubmed_data(get_pubmed_ids("Apoptosis"))  
#M.xml<-fetch_pubmed_data(get_pubmed_ids("Cell cycle")) 
#N.xml<-fetch_pubmed_data(get_pubmed_ids("Rap1 signaling pathway"))  
#O.xml<-fetch_pubmed_data(get_pubmed_ids("Chemokine signaling pathway"))
#p.xml<-fetch_pubmed_data(get_pubmed_ids("Phosphatidylinositol signaling system"))
#Q.xml<-fetch_pubmed_data(get_pubmed_ids("T cell receptor signaling pathway"))
#R.xml<-fetch_pubmed_data(get_pubmed_ids("Estrogen signaling pathway"))
#S.xml<-fetch_pubmed_data(get_pubmed_ids("VEGF signaling pathway"))
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

#----------------------Parameter Values-------------------------------#


#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#

#----------------------Drift------------------------------------------#
Drift.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Drift.Model.1<-Drift.Model.1("1")
test.Drift.Model.1
#----------------------Diffusion--------------------------------------#
Diffusion.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Diffusion.Model.1<-Diffusion.Model.1("1")
test.Diffusion.Model.1
#--------------------- Boundary Model---------------------------------#
Boundary.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Boundary.Model.1<-Boundary.Model.1("1")
test.Boundary.Model.1
#----------------------Estimation-------------------------------------#
Estimation.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Estimation.1<-Estimation.1("1")
test.Estimation.1
#----------------------Joint density----------------------------------#
Density.Joint.1<-function(X)
{
 Table.1.df<-data.frame();
 
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  return(output)
}
test.Density.Joint.1<-Density.Joint.1(X)
test.Density.Joint.1
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

