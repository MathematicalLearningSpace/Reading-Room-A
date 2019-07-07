#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#

#--------------------Work In Progress-------------------------------------#

#------------------------------R API----------------------------------#
library(deSolve);library(ReacTran);library(rootSolve);
library(fda);library(phaseR)
library(pracma);library(GA);library(igraph)
library(NMOF);library(xtable)

#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#-----------------------Review Notes----------------------------------#
#---------------------------------------------------------------------#


#---------------------------------------------------------------------#
#------------------------------Functions------------------------------#
#---------------------------------------------------------------------#

#--------------------------Objective Functions------------------------#
OF.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.OF.Model.1<-OF.Model.1("1")
test.OF.Model.1
#--------------------------Fitness Metrics----------------------------#
FM.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.FM.Model.1<-FM.Model.1("1")
test.FM.Model.1
#--------------------------Loss Functions-----------------------------#
Loss.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Loss.Model.1<-Loss.Model.1("1")
test.Loss.Model.1
#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#

#-----------Genetic Algorithm Optimization----------------------------#
Optimization.GA.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Optimization.GA.Model.1<-Optimization.GA.Model.1("1")
test.Optimization.GA.Model.1
#-----------Differential Evolution Optimization------------------------#
Optimization.DE.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Optimization.DE.Model.1<-Optimization.DE.Model.1("1")
test.Optimization.DE.Model.1
#-----------Particle Swarm Optimization-------------------------------#
Optimization.PS.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Optimization.PS.Model.1<-Optimization.PS.Model.1("1")
test.Optimization.PS.Model.1
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

#---------------References--------------------------------------------#

