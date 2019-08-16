#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#

#--------------------Work In Progress-------------------------------------#
#------------------------------R API----------------------------------#
library(gdata);library(bio3d);library(igraph);library(sna);library(ips);
library(phangorn);library(proteomics)
library(dcGOR);library(MDplot);library(UniProt.ws);library(circlize);
library(BioPhysConnectoR);library(protr)
library(seqinr);library(Biostrings);library(Peptides);
library(PearsonDS);library(xtable);
library(quanteda);library(data.table);library(readr)
library(fuzzyjoin);library(tidytext);library(dplyr);library(tm)
#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#

#\section{Apiaceae}
#\section{Araceae}
#\section{Arecaceae}
#\section{Asteraceae}
#\section{Brassicaceae}
#\section{Commelinaceae}
#\section{Ericaceae}
#\section{Liliaceae}
#\section{Amaryllidoideae}
#\section{Clusiaceae}
#\section{Euphorbiaceae}
#\section{Phyllanthaceae}
#\section{Picrodendraceae}
#\section{Apocynaceae}
#\section{Orchidaceae}
#\section{Poaceae}
#\section{Rosaceae}
#\section{Rubiaceae}
#\section{Rutaceae}
#\section{Sapindaceae}
#-------------------------------------------------------------------
Japan.df<-as.data.frame(BIEN_list_country("Japan"))
China.df<-as.data.frame(BIEN_list_country("China"))
Nepal.df<-as.data.frame(BIEN_list_country("Nepal"))
Italy.df<-as.data.frame(BIEN_list_country("Italy")) 
France.df<-as.data.frame(BIEN_list_country("France"))
England.df<-as.data.frame(BIEN_list_country("England"))
Arizona.df<-as.data.frame(BIEN_list_country("Arizona"))
Mexico.df<-as.data.frame(BIEN_list_country("Mexico"))
Argentina.df<-as.data.frame(BIEN_list_country("Argentina"))
Brazil.df<-as.data.frame(BIEN_list_country("Brazil"))
Chile.df<-as.data.frame(BIEN_list_country("Chile"))
Peru.df<-as.data.frame(BIEN_list_country("Peru"))
Venezuela.df<-as.data.frame(BIEN_list_country("Venezuela"))
#------------------------------------------------------------------
W<-data.frame();X<-data.frame();Y<-data.frame();Z<-data.frame();
#-----------------Phytochemicals----------------------
Phytochemicals.Quantity.80.g<-0 
Phytochemicals.Calories.(kcal)<-0
Phytochemicals.Protein.(g)<-0 
Phytochemicals.Fat.(g)<-0 
Phytochemicals.Fibre.(g)<-0
Phytochemicals.Beta.carotene.(mcg)<-0 
Phytochemicals.Vitamin.A.equivalent.(mcg)<-0 
Phytochemicals.Vitamin.B1.(mg)<-0 
Phytochemicals.Vitamin.B6.(mg)<-0 
Phytochemicals.Vitamin.C.(mg)<-0 
Phytochemicals.Vitamin.E.(mg)<-0 
Phytochemicals.Folate.(mcg)<-0
Phytochemicals.Vitamin.K.(mcg)<-0 
#----------------Mineral Compositions------------------
Mineral.composition.Quantity.80g<-0 
Mineral.composition.Calcium.(mg)<-0 
Mineral.composition.Iodine.(mcg)<-0 
Mineral.composition.Iron.(mg)<-0
Mineral.composition.Magnesium.(mg)<-0 
Mineral.composition.Manganese.(mg)<-0 
Mineral.composition.Phosphorus.(mg)<-0
Mineral.composition.Potassium.(mg)<-0 
Mineral.composition.Zinc.(mg)<-0
Mineral.composition.Selenium.(mcg)<-0
Mineral.composition.Sodium.(mg.100g)<-0 
Mineral.composition.Copper.(mg.100g)<-0

#---------------------------------------------------------------------#
#-----------------------Review Notes----------------------------------#
#---------------------------------------------------------------------#
Review.Notes.1<-function(X)
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
Model.Botany.1<-function(X,p,k,visualization=FALSE)
{
 
 #Logit Model
 #Multinomial Probit Model
 #Multinomial Logit Model
 #Dirchlet-Multinomial Mixture
 #Markov Chain Monte Carlo Algorithm 1
 #Markov Chain Monte Carlo Algorithm 2
 
  A<-matrix(0,nrow=10^p,ncol=10^k)
  B<-matrix(0,nrow=10^p,ncol=10^k)
  C<-matrix(0,nrow=10^p,ncol=10^k)
  I<-diag(10^k)
if(visualization){
    
    Figure.1<-plot(X.Model.Solution.1$Solution[,2], type="l", lty=1,col="black", ylab="Simulated Value", xlab="Temporal Position")
    lines(X.Model.Solution.1$Solution[,3], lty=2,col="green")
    lines(X.Model.Solution.1$Solution[,4], lty=3,col="blue")
    rug(X.Model.Solution.1$Solution[,2], side=4, col="black")
    rug(X.Model.Solution.1$Solution[,3], side=4, col="green")
    rug(X.Model.Solution.1$Solution[,4], side=4, col="blue")
    legend("bottomleft",c("x1","x2","x3"),inset = .01,col=c("black","green","blue"),lwd=2,cex=0.5)
    
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
test.Model.Botany.1<-Model.Botany.1("1",10,1,FALSE)
test.Model.Botany.1
  
Model.Botany.2<-function(X,p,k,visualization=FALSE)
{
  A<-matrix(0,nrow=10^p,ncol=10^k)
  B<-matrix(0,nrow=10^p,ncol=10^k)
  C<-matrix(0,nrow=10^p,ncol=10^k)
  I<-diag(10^k)
if(visualization){
    
    Figure.1<-plot(X.Model.Solution.1$Solution[,2], type="l", lty=1,col="black", ylab="Simulated Value", xlab="Temporal Position")
    lines(X.Model.Solution.1$Solution[,3], lty=2,col="green")
    lines(X.Model.Solution.1$Solution[,4], lty=3,col="blue")
    rug(X.Model.Solution.1$Solution[,2], side=4, col="black")
    rug(X.Model.Solution.1$Solution[,3], side=4, col="green")
    rug(X.Model.Solution.1$Solution[,4], side=4, col="blue")
    legend("bottomleft",c("x1","x2","x3"),inset = .01,col=c("black","green","blue"),lwd=2,cex=0.5)
    
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
test.Model.Botany.2<-Model.Botany.2("1",10,1,FALSE)
test.Model.Botany.2
  
Model.Shikimate.1<-function(X,p,k,visualization=FALSE)
{
  A<-matrix(0,nrow=10^p,ncol=10^k)
  B<-matrix(0,nrow=10^p,ncol=10^k)
  C<-matrix(0,nrow=10^p,ncol=10^k)
  I<-diag(10^k)
if(visualization){
    
    Figure.1<-plot(X.Model.Solution.1$Solution[,2], type="l", lty=1,col="black", ylab="Simulated Value", xlab="Temporal Position")
    lines(X.Model.Solution.1$Solution[,3], lty=2,col="green")
    lines(X.Model.Solution.1$Solution[,4], lty=3,col="blue")
    rug(X.Model.Solution.1$Solution[,2], side=4, col="black")
    rug(X.Model.Solution.1$Solution[,3], side=4, col="green")
    rug(X.Model.Solution.1$Solution[,4], side=4, col="blue")
    legend("bottomleft",c("x1","x2","x3"),inset = .01,col=c("black","green","blue"),lwd=2,cex=0.5)
    
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
test.Model.Shikimate.1<-Model.Shikimate.1("1",10,1,FALSE)
test.Model.Shikimate.1
  

Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Model.1<-Model.1("1")
test.Model.1
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
#------Evolutionary Learning of Globally Optimal Classification and Regression Trees 
X1 <- rep(seq(0.25, 1.75, 0.5), each = 8)
X2 <- rep(seq(0.25, 1.75, 0.5), 8)
Y.1 <- rep(1, 32)
#----------------Rule 1
Y.1[(X1 < 1 & X2 < 1) | (X1 > 1 & X2 > 1)] <- 2
#----------------Labels of Factors
Y.1 <- factor(Y.1, labels = c("A", "not A"))
Z.df <- data.frame(Y.1, X1, X2)
#----------------Learning Space
Tree.Model.1<-evtree::evtree(Y.1 ~ ., data = Z.df, minbucket = 1, minsplit = 2)
#----------------Presentation
output<-list()
output$X<-X
output$Tree.Model.1<-Tree.Model.1
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
