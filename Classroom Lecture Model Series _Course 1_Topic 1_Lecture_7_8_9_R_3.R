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

#---Parameter Systems-------------#


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

LASSO.Model.1<-function(X)
 {
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
 
  output<-list()
  output$X<-X
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

#--------Lead - Lag Estimation----------------------------------------#

#---------------------------------------------------------------------#
#------------------------------Analysis-------------------------------#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#------------------------------Tables---------------------------------#
#---------------------------------------------------------------------#

#----------------Table 1---------------------#
#----------------Table 2---------------------#
#----------------Table 3---------------------#
#----------------Table 4---------------------#

#---------------------------------------------------------------------#
#------------------------------Figures--------------------------------#
#---------------------------------------------------------------------#

#----------------Figure 1---------------------#
#----------------Figure 2---------------------#
#----------------Figure 3---------------------#
#----------------Figure 4---------------------#




