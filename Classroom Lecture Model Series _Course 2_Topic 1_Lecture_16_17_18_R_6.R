#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#

#--------------------Work In Progress-------------------------------------#

#------------------------------R API-----------------------------------#
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
Signal.Transduction.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Signal.Transduction.Model.1<-Signal.Transduction.Model.1("1")
Signal.Transduction.Model.1
#------------Longevity-----------------#
#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
Longevity.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
Longevity.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, 
						 func = Longevity.system.equation.model.1, parms = Params.Equation) 

#------------Immune System-------------#

Immune.System.Model<-function(X)
{
 Immune.System.df<-as.data.frame(rbind(c("04640","Hematopoietic cell lineage"),
                       c("04610","Complement and coagulation cascades"),
                       c("04611","Platelet activation"),
                       c("04620","Toll-like receptor signaling pathway"),
                       c("04624","Toll and Imd signaling pathway"),
                       c("04621","NOD-like receptor signaling pathway"),
                       c("04622","RIG-I-like receptor signaling pathway"),
                       c("04623","Cytosolic DNA-sensing pathway"),
                       c("04625","C-type lectin receptor signaling pathway"),
                       c("04612","Antigen processing and presentation"),
                       c("04660","T cell receptor signaling pathway"),
                       c("04658","Th1 and Th2 cell differentiation"),
                       c("04659","Th17 cell differentiation"),
                       c("04657","IL-17 signaling pathway"),
                       c("04662","B cell receptor signaling pathway"),
                       c("04664","Fc epsilon RI signaling pathway"),
                       c("04666","Fc gamma R-mediated phagocytosis"),
                       c("04670","Leukocyte transendothelial migration"),
                       c("04672","Intestinal immune network for IgA production"),
                       c("04062","Chemokine signaling pathway")))

 setwd("Immune System Model")
  Immune.System.1.Name=""; Immune.System.1<-xml2::read_xml("hsa04640.xml")
  Immune.System.2.Name=""; Immune.System.2<-xml2::read_xml("hsa04610.xml")
  Immune.System.3.Name=""; Immune.System.3<-xml2::read_xml("hsa04620.xml")
  #Immune.System.4.Name=""; Immune.System.4<-xml2::read_xml("hsa04624.xml")
  Immune.System.5.Name=""; Immune.System.5<-xml2::read_xml("hsa04621.xml")
  Immune.System.6.Name=""; Immune.System.6<-xml2::read_xml("hsa04622.xml")
  Immune.System.7.Name=""; Immune.System.7<-xml2::read_xml("hsa04623.xml")
  Immune.System.8.Name=""; Immune.System.8<-xml2::read_xml("hsa04625.xml")
  Immune.System.9.Name=""; Immune.System.9<-xml2::read_xml("hsa04612.xml")
  Immune.System.10.Name=""; Immune.System.10<-xml2::read_xml("hsa04660.xml")
  Immune.System.11.Name=""; Immune.System.11<-xml2::read_xml("hsa04658.xml")
  Immune.System.12.Name=""; Immune.System.12<-xml2::read_xml("hsa04659.xml")
  Immune.System.13.Name=""; Immune.System.13<-xml2::read_xml("hsa04657.xml")
  Immune.System.14.Name=""; Immune.System.14<-xml2::read_xml("hsa04662.xml")
  Immune.System.15.Name=""; Immune.System.15<-xml2::read_xml("hsa04664.xml")
  Immune.System.16.Name=""; Immune.System.16<-xml2::read_xml("hsa04666.xml")
  Immune.System.17.Name=""; Immune.System.17<-xml2::read_xml("hsa04670.xml")
  Immune.System.18.Name=""; Immune.System.18<-xml2::read_xml("hsa04672.xml")
  Immune.System.19.Name=""; Immune.System.19<-xml2::read_xml("hsa04062.xml")
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
	
}	


#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
Immune.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
Immune.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, 
					      func = Immune.system.equation.model.1, parms = Params.Equation) 

#------------wNT-----------------------#
#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
Wnt.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
Wnt.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, 
					   func = Wnt.system.equation.model.1, parms = Params.Equation) 

#------------Notch---------------------#
#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
Notch.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
Notch.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, 
					     func = Notch.system.equation.model.1, parms = Params.Equation) 

#------------HedgeHog------------------#
#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
Hedgehog.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
Hedgehog.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, 
						func = Hedgehog.system.equation.model.1, parms = Params.Equation) 

#------------TGF-Beta------------------#
#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
TGF.Beta.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
TGF.Beta.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, 
						func = TGF.Beta.system.equation.model.1, parms = Params.Equation) 

#------------VEGF----------------------#
#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
VEGF.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
VEGF.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, 
					    func = VEGF.system.equation.model.1, parms = Params.Equation) 

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

