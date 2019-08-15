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
F.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.F.1<-F.1("1")
test.F.1
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

#-------------------Semi-Markov Models-------------------------#
Model.Markov.Semi.1<-function(X,p,k,visualization=FALSE)
{
  X1<-matrix(0,nrow=10^p,ncol=10^k)
  X2<-matrix(0,nrow=10^p,ncol=10^k)
  X3<-matrix(0,nrow=10^p,ncol=10^k)
  I<-diag(10^k)
  Table.1.df<-data.frame(X1)
  for(i in 1:nrow(X1))
  {
    X<-t(X) 
  }
  
  if(visualization){
    png(file = stringr::str_c('Figures/1/Example_',1,'_Figure_',1,'.png'))
    plot(yout[,-1], type = "l")
    dev.off()
    png(file = stringr::str_c('Figures/1/Example_',2,'_Figure_',2,'.png'))
    plot(yout[,-1], type = "l")
    dev.off()
    png(file = stringr::str_c('Figures/1/Example_',3,'_Figure_',3,'.png'))
    plot(yout[,-1], type = "l")
    dev.off()
    png(file = stringr::str_c('Figures/1/Example_',4,'_Figure_',4,'.png'))
    plot(yout[,-1], type = "l")
    dev.off()
  }
  Table.1.TEX<-xtable::xtable(Table.1.df)
  output<-list()
  output$Table.1<-Table.1.df
  output$Table.1.TEX<-Table.1.TEX
  return(output)
}
test.Model.Markov.Semi.1<-Model.Markov.Semi.1("1",10,1,FALSE)
test.Model.Markov.Semi.1

#-------------------Model Graph Exponential-------------------------#
Model.Graph.Exponential.1<-function(X,p,k,visualization=FALSE)
{
  A<-matrix(0,nrow=10^p,ncol=10^k)
  B<-matrix(0,nrow=10^p,ncol=10^k)
  C<-matrix(0,nrow=10^p,ncol=10^k)
  I<-diag(10^k)
  Table.1.df<-data.frame(X)
  if(visualization){Figure.1<-plot(X)}
  Table.1.TEX<-xtable::xtable(Table.1.df)
  output<-list()
  output$Table.1<-Table.1.df
  output$Table.1.TEX<-Table.1.TEX
  return(output)
}
test.Model.Graph.Exponential.1<-Model.Graph.Exponential.1("1",10,1,FALSE)
test.Model.Graph.Exponential.1


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

#---------Multinomial Confidence Intervals----------------------------#
Confidence.Intervals.Multinomial.1<-function(X)
 {
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 
 }
test.Confindence.Intervals.Multinomial.1<-Confindence.Intervals.Multinomial.1("1")
test.Confindence.Intervals.Multinomial.1

#---------Path Probabilities------------------------------------------#
Path.Probabilities.1<-function(X)
 {
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 
 }
test.Path.Probabilties.1<-Path.Probabilities.1("1")
test.Path.Probabilities.1

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
#----------Table 1------#
Table.1.TeX<-xtable::xtable(Table.1.df)
#----------Table 2------#
Table.2.TeX<-xtable::xtable(Table.2.df)
#----------Table 3------#
Table.3.TeX<-xtable::xtable(Table.3.df)
#----------Table 4------#
Table.4.TeX<-xtable::xtable(Table.4.df)

#---------------------------------------------------------------------#
#----------Figures for Classroom Presentation-------------------------#
#---------------------------------------------------------------------#

#-----------Figure 1---------------#
#-----------Figure 2---------------#
#-----------Figure 3---------------#
#-----------Figure 4---------------#
#-----------Figure 5---------------#
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



