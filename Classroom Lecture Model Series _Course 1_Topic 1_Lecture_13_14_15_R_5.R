#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#

#--------------------Work In Progress-------------------------------------#

#------------------------------R API----------------------------------#
library(deSolve);library(ReacTran);library(rootSolve);
library(fda);library(phaseR)
library(pracma);library(GA);library(igraph);
library(markovchain);library(HMM);library(adaptMCMC);library(mcmc);library(MCMCpack)
library(combinat);library(coda);library(rbenchmark);library(stringi);library(stringr)
library(Matrix);library(corrplot);library(xtable); library(coda)

library(qgraph);library(qlcMatrix);library(rcrossref);library(PearsonDS);library(catnet)

#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#-----------------------Review Notes----------------------------------#
#---------------------------------------------------------------------#


#---------------------------------------------------------------------#
#------------------------------Functions------------------------------#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#

#-----------Markov Model Parameters-----------------------------------#
Model.Markov.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Model.Markov.1<-Model.Markov.1("1")
test.Model.Markov.1
#-----------Generate Matrices-----------------------------------------#

#-----------Dirichlert Prior------------------------------------------#

#-----------Markov Chains---------------------------------------------#
Model.Markov.Chain.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Model.Markov.Chain.1<-Model.Markov.Chain.1("1")
test.Model.Markov.Chain.1

#---------------------------------------------------------------------#
#------------------------------Analysis-------------------------------#
#---------------------------------------------------------------------#

#---------Multinomial Confidence Intervals----------------------------#

#---------Path Probabilities------------------------------------------#

#------------Graph Theory:Distance------------------------------------#
Graph.Theory.Distance.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Graph.Theory.Distance.1<-Graph.Theory.Distance.1("1")
test.Graph.Theory.Distance.1
#------------Graph Theory:Connectivity--------------------------------#
Graph.Theory.Connectivity.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Graph.Theory.Connectivity.1<-Graph.Theory.Connectivity.1("1")
test.Graph.Theory.Connectivity.1
#------------Graph Theory:Spectra-------------------------------------#
Graph.Theory.Spectra.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Graph.Theory.Spectra.1<-Graph.Theory.Spectra.1("1")
test.Graph.Theory.Spectra.1
#----Correlation Analysis---------------------------------------------#


#---------------------------------------------------------------------#
#------------------------------Tables---------------------------------#
#---------------------------------------------------------------------#

#------------Table 1---------------#
#------------Table 2---------------#
#------------Table 3---------------#
#------------Table 4---------------#


#---------------------------------------------------------------------#
#----------Figures for Classroom Presentation-------------------------#
#---------------------------------------------------------------------#

#-----------Figure 1---------------#
#-----------Figure 2---------------#
#-----------Figure 3---------------#
#-----------Figure 4---------------#
#-----------Figure 5---------------#



