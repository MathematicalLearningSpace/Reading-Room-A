#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#
#--------------------Work In Progress-------------------------------------#
#------------------------------R API----------------------------------#
library(deSolve);library(ReacTran);library(rootSolve);
library(fda);library(phaseR)
library(pracma);library(GA);library(igraph)
library(sde);library(yuima);library(MsdeParEst);library(rugarch)

#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#

Gastric.Cancer.Genes<-c("CDX2", "MUC2", "REG4", "CDH17", "MDR1", "SHH","p53", "p21", 
                        "BAX", "p48", "GADD45", "BAK", "POLK","Retinoic.Acid", "RAR.Beta", "RXR",
                        "DV1", "GSK.3Beta", "Beta.Catenin", "Axin", "APC", "CK1.alpha","GBP","EFG", 
                        "ERBB2", "SHC", "GRB2", "SOS", "RAS", "RAF", "MEK", "ERK.1","PI3K", "PIP3", 
                        "AKT", "mTOR", "p53", "S6K", "BCL2","TGF.Beta","TGF.BetaRI","TGF.BetaRII", "SMAD.2", "SMAD.4", "p15", "p21",
                        "HGF", "c.MET", "GRB2", "SOS", "RAS", "RAF", "MEK", "ERK.1",
                        "FGF", "FGFR2", "GAB1", "PI3K", "PIP3", "AKT", "mTOR", "GSK.3Beta")


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

#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#

Model.SDE.1<-function(X,s = 1,point.vector=c(x=2,y=-2,z=-2))
{
drift.Model.1 <- expression(x*y , 
                            y , 
                            z*y) 
diffusion.Model.1 <- rep(expression(10^(-1),3)
SDE.1 <- Sim.DiffProc::snssde3d(x0=point.vector,drift=drift.Model.1,diffusion=diffusion.Model.1,M=10^4)
Model.Summary.1<-summary.1(SDE.1, at = s)
output<-list()
output$Model.Summary<-Model.Summary.1  
  return(output)
}  

Model.SDE.2<-function(X)
{
  output<-list()
  
  return(output)
}  

Model.SDE.3<-function(X)
{
  output<-list()
  
  return(output)
}  


#--------MLE, Lasso  and Bayesian Models------------------------------#
MLE.Model.1<-function(X)
 {
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
           
            
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 
 }
MLE.Model.1<-MLE.Model.1("1")
test.MLE.Model.1

LASSO.Model.1<-function(X,p,k)
 {
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
 p<-nrow(X)
 k<-ncol(X)           
 Y<-matrix(1,nrow=10^p,ncol=1) 
 I<-diag(10^k)           
 beta.1<-solve(t(X)%*%X+lambda.1*I)%*%t(X)%*%Y
 Y.Pred.1<-t(X)%*%beta.1     
  output<-list()
  output$X<-X
  output$Y<-Y
  output$beta.1<-beta.1
  output$Y.Pred.1<-Y.Pred.1
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 
 }
LASSO.Model.1<-LASSO.Model.1("1")
test.LASSO.Model.1

Bayesian.Model.1<-function(X)
 {
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 
 }
Bayesian.Model.1<-Bayesian.Model.1("1")
test.Bayesian.Model.1
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
test.Estimation.Model.1<-Estimation.1("1")
test.Estimation.Model.1
#----------------------Joint density----------------------------------#
Density.Joint.1<-function(X)
 {
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 
 }
test.Density.Joint.1<-Density.Joint.1("1")
test.Density.Joint.1
#-------Maximum Likelihood Estimation---------------------------------#
Estimation.MLE.1<-function(X)
 {
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 
 }
test.Estimation.MLE.1<-Estimation.MLE.1("1")
test.Estimation.MLE.1
#-------Bayesian Model with MCMC Estimation --------------------------#
Estimation.Bayesian.MCMC.1<-function(X)
 {
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 
 }
test.Estimation.Bayesian.MCMC.1<-Estimation.Bayesian.MCMC.1("1")
test.Estimation.Bayesian.MCMC.1
#--------------Prior Specification: Probability Density Function------#
Estimation.Bayesian.Prior.1<-function(X)
 {
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 
 }
test.Estimation.Bayesian.Prior.1<-Estimation.Bayesian.Prior.1("1")
test.Estimation.Bayesian.Prior.1
#---------------------Model Testing-----------------------------------#
Testing.Model.1<-function(X)
 {
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Testing.Model.1<-Testing.Model.1("1")
test.Testing.Model.1

#--------Lead - Lag Estimation----------------------------------------#
Estimation.Lead.Lag.1<-function(X)
 {
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Estimation.Lead.Lag.1<-Estimation.Lead.Lag.1("1")
test.Estimation.Lead.Lag.1

#---------------------------------------------------------------------#
#------------------------------Analysis-------------------------------#
#---------------------------------------------------------------------#
Analysis.1<-function(X)
 {
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Analysis.1<-Analysis.1("1")
test.Analysis.1
#---------------------------------------------------------------------#
#------------------------------Tables---------------------------------#
#---------------------------------------------------------------------#

#----------------Table 1---------------------#
#----------------Table 2---------------------#
#----------------Table 3---------------------#
#----------------Table 4---------------------#
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

#----------------Figure 1---------------------#
#----------------Figure 2---------------------#
#----------------Figure 3---------------------#
#----------------Figure 4---------------------#

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


