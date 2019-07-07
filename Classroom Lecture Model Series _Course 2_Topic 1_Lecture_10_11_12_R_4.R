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

#------------MAPK-----------------#
#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
MAPK.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
MAPK.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, func = system.equation.model.1, parms = Params.Equation) 

#------------Calcium--------------#
#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
Calcium.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
Calcium.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, func = system.equation.model.1, parms = Params.Equation) 

#------------cAMP-----------------#
#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
cAMP.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
cAMP.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, func = system.equation.model.1, parms = Params.Equation) 

#------------Chemokine------------#
#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
Chemokine.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
Chemokine.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, func = system.equation.model.1, parms = Params.Equation) 

#------------nf-Kappa-------------#
#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
nf.kappa.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
nf.kappa.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, func = system.equation.model.1, parms = Params.Equation) 

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
