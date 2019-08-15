#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#
#----------------------------------R Source Files------------------------------------#
#-------------------------R API----------------------#
library(easyPubMed);library(bio3d);library(readr);library(CHNOSZ);
library(stringr);library(Peptides);library(Biostrings)
library(seqinr);library(seqLogo);library(msa);library(ape);
library(dtw);library(dtwclust);library(odseq);library(rphast)
library(plyr)
#----------------------------------Data----------------------------------------------#
W<-data.frame();X<-data.frame();Y<-data.frame();Z<-data.frame();
#----------------------------------Transformations-----------------------------------#
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
#----------------------------------User Defined Modules and Functions----------------#
F.1<-function(x)
 {
   return(x)
 }
#----------------Ribosomes-----------------#
Ribosome.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
 
 system.equation.model.1<-function(times, variables.intitial.test, parameters.test)
{
  with(as.list(c(parameters.test, variables.intitial.test)), 
       {
#--------------Group 1---------------------
d.1.X1.d.t.1<--a11*X1 + a12*X2
#--------------Group 2---------------------
d.1.X2.d.t.1<-a72*X2 - a81*X1
res <- c(d.1.X1.d.t.1,d.1.X2.d.t.1)
list(res)
})
}
 
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Ribosome.Model.1<-Ribosome.Model.1("1")
test.Ribosome.Model.1
#----------------ER Rough and Smooth-------#
ER.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
 system.equation.model.1<-function(times, variables.intitial.test, parameters.test)
{
  with(as.list(c(parameters.test, variables.intitial.test)), 
       {
#--------------Group 1---------------------
d.1.X1.d.t.1<--a11*X1 + a12*X2
#--------------Group 2---------------------
d.1.X2.d.t.1<-a72*X2 - a81*X1
res <- c(d.1.X1.d.t.1,d.1.X2.d.t.1)
list(res)
})
}
 
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.ER.Model.1<-ER.Model.1("1")
test.ER.Model.1
#----------------Tridiagonal Systems-------#
Tridiagonal.Systems.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
 
 system.equation.model.1<-function(times, variables.intitial.test, parameters.test)
{
  with(as.list(c(parameters.test, variables.intitial.test)), 
       {
#--------------Group 1---------------------
d.1.X1.d.t.1<--a11*X1 + a12*X2
#--------------Group 2---------------------
d.1.X2.d.t.1<-a72*X2 - a81*X1
res <- c(d.1.X1.d.t.1,d.1.X2.d.t.1)
list(res)
})
}
 
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Tridiagonal.Systems.Model.1<-Tridiagonal.Systems.Model.1("1")
test.Tridiagonal.Systems.Model.1
#----------------------------------Network Designs-----------------------------------#

#----------------------------------Equation Systems----------------------------------#
Tridiagonal.Systems.Model.2<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
 system.equation.model.1<-function(times, variables.intitial.test, parameters.test)
{
  with(as.list(c(parameters.test, variables.intitial.test)), 
       {
#--------------Group 1---------------------
d.1.X1.d.t.1<--a11*X1 + a12*X2
#--------------Group 2---------------------
d.1.X2.d.t.1<-a72*X2 - a81*X1
res <- c(d.1.X1.d.t.1,d.1.X2.d.t.1)
list(res)
})
}
 
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Tridiagonal.Systems.Model.2<-Tridiagonal.Systems.Model.2("1")
test.Tridiagonal.Systems.Model.2
#----------------------------------Parameter Tables----------------------------------#

#----------------------------------Network Analysis----------------------------------#
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
Model.Optimization.GA.2 <- ga(type = "real-valued", fitness =  function(x) - F.1(x[1]),min = 0, max =1, popSize = 10, maxiter = 3)
Model.Optimization.GA.Summary<-summary(Model.Optimization.GA.2)
Model.Optimization.GA.Summary$popSize
Model.Optimization.GA.Summary$maxiter
Model.Optimization.GA.Summary$elitism
Model.Optimization.GA.Summary$pcrossover
Model.Optimization.GA.Summary$pmutation
Model.Optimization.GA.Summary$fitness
Model.Optimization.GA.Summary$solution
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Optimization.Model.1<-Optimization.Model.1("1")
test.Optimization.Model.1
#----------------------------------Natural Language Description----------------------#

#----------------------------------Tables--------------------------------------------#

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
#----------------------------------Figures-------------------------------------------#

#----------Figure 1 Group A---------#
#png(file = stringr::str_join("Figures//1//Example_",1,"_Figure_","0A",".png"))
plot(system.equation.model.test.solution.1[,2], type="l", main="Time and Phase")
lines(system.equation.model.test.solution.1[,3],lty = 1)
lines(system.equation.model.test.solution.1[,4],lty = 2)
lines(system.equation.model.test.solution.1[,5],lty = 3)
lines(system.equation.model.test.solution.1[,6],lty = 4)
lines(system.equation.model.test.solution.1[,7],lty = 5)
lines(system.equation.model.test.solution.1[,8],lty = 6)
legend("bottomleft", legend = paste(seq(1:8),Model.Names[2:9]),lty = 1:8, cex = .5, y.intersp = 1)
#dev.off()
#png(file = stringr::str_join("Figures//1//Example_",2,"_Figure_","0A",".png"))
plot(Model.Optimization.GA.Summary$solution, type="l", main="Optimal Path")
legend("bottomleft", legend = paste("A-","Phase"),lty = 1:8, cex = .5, y.intersp = 1)
#dev.off()
#png(file = stringr::str_join("Figures//1//Example_",3,"_Figure_","0A",".png"))
plot(Model.1, type="l", main="Time Variation")
legend("bottomleft", legend = paste("A-","Time"),lty = 1:1, cex = .5, y.intersp = 1)
#dev.off()
#------------------------Correlation Matrix ------------------------------------
#png(file = stringr::str_join("Figures//1//Example_",4,"_Figure_","0A",".png"))
#car::scatterplotMatrix(system.equation.model.test.solution.1[,4:8])
#dev.off()
#----------Figure 2 Group B---------#

#----------Figure 3 Group C---------#

#----------Figure 4 Group D---------#
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
        
#----------------------------------Discussion----------------------------------------#
