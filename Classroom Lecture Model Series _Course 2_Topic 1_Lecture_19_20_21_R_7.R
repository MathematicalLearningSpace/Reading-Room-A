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
setwd("Signal Transduction Models")
  #Signal.Transduction.Name.1="Two-component system";Signal.Transduction.1<-xml2::read_xml("ko02020.xml")
  Signal.Transduction.Name.2="Ras signaling pathway";Signal.Transduction.2<-xml2::read_xml("hsa04014.xml")
  Signal.Transduction.Name.3="Rap1 signaling pathway";Signal.Transduction.3<-xml2::read_xml("hsa04015.xml")
  Signal.Transduction.Name.4="MAPK signaling pathway";Signal.Transduction.4<-xml2::read_xml("hsa04010.xml")
  #Signal.Transduction.Name.5="MAPK signaling pathway - fly";Signal.Transduction.5<-xml2::read_xml("dme04013.xml")
  #Signal.Transduction.Name.6="MAPK signaling pathway - plant";Signal.Transduction.6<-xml2::read_xml("ath04016.xml")
  #Signal.Transduction.Name.7="MAPK signaling pathway - yeast";Signal.Transduction.7<-xml2::read_xml("sce04011.xml")
  Signal.Transduction.Name.8="ErbB signaling pathway";Signal.Transduction.8<-xml2::read_xml("hsa04012.xml")
  Signal.Transduction.Name.9="Wnt signaling pathway";Signal.Transduction.9<-xml2::read_xml("hsa04310.xml")
  Signal.Transduction.Name.10="Notch signaling pathway";Signal.Transduction.10<-xml2::read_xml("hsa04330.xml")
  Signal.Transduction.Name.11="Hedgehog signaling pathway";Signal.Transduction.11<-xml2::read_xml("hsa04340.xml")
  #Signal.Transduction.Name.12="Hedgehog signaling pathway - fly";Signal.Transduction.12<-xml2::read_xml("dme04341.xml")
  Signal.Transduction.Name.13="TGF-beta signaling pathway";Signal.Transduction.13<-xml2::read_xml("hsa04350.xml")
  Signal.Transduction.Name.14="Hippo signaling pathway";Signal.Transduction.14<-xml2::read_xml("hsa04390.xml")
  #Signal.Transduction.Name.15="Hippo signaling pathway - fly";Signal.Transduction.15<-xml2::read_xml("dme04391.xml");
  #Signal.Transduction.Name.16="Hippo signaling pathway - multiple species";Signal.Transduction.16<-xml2::read_xml("hsa04392.xml")
  Signal.Transduction.Name.17="VEGF signaling pathway";Signal.Transduction.17<-xml2::read_xml("hsa04370.xml")
  Signal.Transduction.Name.18="Apelin signaling pathway";Signal.Transduction.18<-xml2::read_xml("hsa04371.xml");
  Signal.Transduction.Name.19="Jak-STAT signaling pathway";Signal.Transduction.19<-xml2::read_xml("hsa04630.xml");
  Signal.Transduction.Name.20="NF-kappa B signaling pathway";Signal.Transduction.20<-xml2::read_xml("hsa04064.xml")
  Signal.Transduction.Name.21="TNF signaling pathway";Signal.Transduction.21<-xml2::read_xml("hsa04668.xml")
  Signal.Transduction.Name.22="HIF-1 signaling pathway";Signal.Transduction.22<-xml2::read_xml("hsa04066.xml")
  Signal.Transduction.Name.23="FoxO signaling pathway";Signal.Transduction.23<-xml2::read_xml("hsa04068.xml")
  Signal.Transduction.Name.24="Calcium signaling pathway";Signal.Transduction.24<-xml2::read_xml("hsa04020.xml")
  Signal.Transduction.Name.25="Phosphatidylinositol signaling system";Signal.Transduction.25<-xml2::read_xml("hsa04070.xml")
  Signal.Transduction.Name.26="Phospholipase D signaling pathway";Signal.Transduction.26<-xml2::read_xml("hsa04072.xml")
  Signal.Transduction.Name.27="Sphingolipid signaling pathway";Signal.Transduction.27<-xml2::read_xml("hsa04071.xml")
  Signal.Transduction.Name28="cAMP signaling pathway";Signal.Transduction.28<-xml2::read_xml("hsa04024.xml")
  Signal.Transduction.Name.29="cGMP-PKG signaling pathway";Signal.Transduction.29<-xml2::read_xml("hsa04022.xml")
  Signal.Transduction.Name.30="PI3K-Akt signaling pathway";Signal.Transduction.30<-xml2::read_xml("hsa04151.xml")
  Signal.Transduction.Name.31="AMPK signaling pathway";Signal.Transduction.31<-xml2::read_xml("hsa04152.xml")
  Signal.Transduction.Name.32="mTOR signaling pathway";Signal.Transduction.32<-xml2::read_xml("hsa04150.xml")
  Signal.Transduction.Name.33="Plant hormone signal transduction";Signal.Transduction.33<-xml2::read_xml("ath04075.xml")

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

#------------Apelin-----------------#
#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
Apelin.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
Apelin.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, func = Apelin.system.equation.model.1, parms = Params.Equation) 

#------------Hippo------------------#
#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
Hippo.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
Hippo.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, func = Hippo.system.equation.model.1, parms = Params.Equation) 

#------------Jak-STAT---------------#
#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
Jak.STAT.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
Jak.STAT.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, func = Jak.STAT.system.equation.model.1, parms = Params.Equation) 

#------------Insulin----------------#
#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
Insulin.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
Insulin.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, func = Insulin.system.equation.model.1, parms = Params.Equation) 

#------------WDR--------------------#
#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
WDR.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
WDR.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, func = WDR.system.equation.model.1, parms = Params.Equation) 

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











