#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#
#--------------------Work In Progress-------------------------------------#

#------------------------------R API----------------------------------#
library(deSolve);library(ReacTran);library(rootSolve);
library(fda);library(phaseR)
library(pracma);library(GA);library(igraph)
library(NMOF);library(xtable)
library(tm);library(topicmodels);library(wordcloud);

#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#-----------------------Review Notes----------------------------------#
#---------------------------------------------------------------------#

Literature.Search.Example<-function(Topics,Language.Model)
  {
    Review.Notes.Formatted<-list()
    Review.Notes.Unformatted<-list()
    #-------------------Designed in the Classroom---------------------#
    output<-list()
    output$Topics<-Topics
    output$Review.Notes.Formatted<-Review.Notes.Formatted
    output$Review.Notes.Unformatted<-Review.Notes.Unformatted
    return(output)
  }
test.Literature.Search.Example<-Literature.Search.Example("Equations","English")
test.Literature.Search.Example
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

#-------Fitness, Distance Metrics and Loss Functions I, II,III--------#

Metric.Fitness.1<-function(X){return(X)}
Metric.Fitness.2<-function(X){return(X)}
Metric.Fitness.3<-function(X){return(X)}
Metric.Fitness.4<-function(X){return(X)}


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
  operator.population<-function(x){ 
    #-----------------------------Dynamics-----------#
  rate.survival<-0;rate.birth<-0;rate.migration<-0
    #--------------Birth rates
    #--------------Survival rates
    #--------------Migration rates
    
    return(x)}
  operator.mutation<-function(x){return(x)}
  operator.crossover<-function(x){return(x)}
  operator.selection<-function(x){return(x)}
  operator.fitness<-function(x){return(x)}
  operator.sensitivity<-function(x){return(x)}
  operator.stability<-function(x){return(x)}
  
  if(visualization){
     png(file = stringr::str_c('Figures/1/Example_',1,'_Figure_',1,'.png'))
    plot(solution.1[,-1], type = "l")
    dev.off()
    png(file = stringr::str_c('Figures/1/Example_',2,'_Figure_',2,'.png'))
    plot(solution.1[,-1], type = "l")
    dev.off()
    png(file = stringr::str_c('Figures/1/Example_',3,'_Figure_',3,'.png'))
    plot(solution.1[,-1], type = "l")
    dev.off()
    png(file = stringr::str_c('Figures/1/Example_',4,'_Figure_',4,'.png'))
    plot(solution.1[,-1], type = "l")
    dev.off()
    }
  
  output<-list()
  output$X<-X
  output$Solution.1<-solution.1
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  output$Selection<-operator.selection
  output$Crossover<-operator.crossover
  output$Mutation<-operator.mutation
  output$Fitness<-operator.fitness
  output$Sensitivity<-operator.sensitivity
  output$Stability<-operator.stability
  output$Population<-population
  output$Rates.Survival<-rate.survival
  output$Rates.Birth<-rate.birth
  output$Rates.Migration<-rate.migration
  return(output)
 }
test.Optimization.GA.Model.1<-Optimization.GA.Model.1("1")
test.Optimization.GA.Model.1
#-----------Differential Evolution Optimization------------------------#
Optimization.DE.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  
   if(visualization){
    png(file = stringr::str_c('Figures/1/Example_',1,'_Figure_',1,'.png'))
    plot(solution.1[,-1], type = "l")
    dev.off()
    png(file = stringr::str_c('Figures/1/Example_',2,'_Figure_',2,'.png'))
    plot(solution.1[,-1], type = "l")
    dev.off()
    png(file = stringr::str_c('Figures/1/Example_',3,'_Figure_',3,'.png'))
    plot(solution.1[,-1], type = "l")
    dev.off()
    png(file = stringr::str_c('Figures/1/Example_',4,'_Figure_',4,'.png'))
    plot(solution.1[,-1], type = "l")
    dev.off()
    }
  output<-list()
  output$X<-X
  output$Solution.1<-solution.1
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
#----------References-------------------------------------------------#

