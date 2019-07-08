#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#

#--------------------Work In Progress-------------------------------------#

#------------------------------R API----------------------------------#
library(deSolve);library(ReacTran);library(rootSolve);
library(fda);library(phaseR)
library(pracma);library(GA);library(igraph)
library(NMOF);library(xtable);

library(OptimalCutpoints);library(Rsolnp);library(NMOF);
library(PearsonDS)
library(rcellminer)
library(cvxclustr)
library(HistogramTools)

#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#

selectedOncogenes <- c("ABL1", "ALK", "BRAF", "CCND1", "CCND3", "CCNE1", "CCNE2", 
                       "CDC25A", "EGFR", "ERBB2", "EZH2", "FOS", "FOXL2", "HRAS", 
                       "IDH1", "IDH2", "JAK2", "KIT", "KRAS", "MET", "MOS", "MYC", 
                       "NRAS", "PDGFB", "PDGFRA", "PIK3CA", "PIK3CB", "PIK3CD", 
                       "PIK3CG", "PIM1", "PML", "RAF1", "REL", "RET", "SRC", "STK11", 
                       "TP63", "WNT10B", "WNT4", "WNT2B", "WNT9A", "WNT3", "WNT5A", 
                       "WNT5B", "WNT10A", "WNT11", "WNT2", "WNT1", "WNT7B", "WISP1", 
                       "WNT8B", "WNT7A", "WNT16", "WISP2", "WISP3", "FZD5", "FZD1")


drugAct <- exprs(getAct(rcellminerData::drugData))
molData <- getMolDataMatrices()
plots <- c("mut", "drug", "cop", "xai", "pro")

Data.design.1 <- function(n, K, constant = TRUE, sigma = 2, oFrac = 0.1) 
{
pIVpars <- list(m=5.1, nu=1, location=0.5, scale=2);
X <- array(rpearsonIV(n * k,params=pIVpars), dim = c(n, k));
b <- rpearsonIV(k,params=pIVpars)
Y<- X %*% b + rpearsonIV(n,params=pIVpars)*0.25 
 output<-list()
 output$X<-X
 output$Y<-Y
 return(output)
} 
test.Data.design.1<-Data.design.1(100,3)
test.Design.design.1
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
Optimization.DE.Model.1<-function(X,Y)
 {
  n <- 100;k <- 5;
  popsize <- 100;generations <- 500; 
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
   de.parms <- list(min = rep(-1,k), max = rep( 1,k), nP = popsize, nG = generations, F = 0.7, 
           CR = 0.9, loopOF = FALSE, printBar = FALSE, printDetail = FALSE) 
 Solution.DE.1 <- DEopt(OF = Y$Specification.1, algo = de.parms, Data =X) 

  output<-list()
  output$X<-X
  output$Solution.DE<-Solution.DE.1
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Optimization.DE.Model.1<-Optimization.DE.Model.1("1")
test.Optimization.DE.Model.1
#-----------Particle Swarm Optimization-------------------------------#
Optimization.PS.Model.1<-function(X,Y)
 {
 popsize <- 100;generations <- 500; 
 n <- 100;k <- 5;
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  ps.parms <- list(min = rep(-1,k), max = rep( 1,k), c1 = 0.9, c2 = 0.9, iner = 0.9, initV = 1, 
           nP = popsize, nG = generations, maxV = 5, loopOF = FALSE, printBar = FALSE, printDetail = FALSE) 
  Solution.PS.1 <- PSopt(OF = Y$Specification.1, algo = de.parms, Data =X) 
 
 output<-list()
  output$X<-X
  output$Solution.PS<-Solution.PS.1
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

png(file = stringr::str_join("Figures//Example_",1,"_Figure_","1",".png"))
Figure.1<-plotCellMiner(drugAct, molData, plots=plots, nsc="94600", gene="", verbose=TRUE)
dev.off()

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
#---------------References--------------------------------------------#

