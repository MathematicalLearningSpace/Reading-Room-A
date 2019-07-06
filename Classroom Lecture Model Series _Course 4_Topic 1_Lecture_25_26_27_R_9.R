#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#

#--------------------Work In Progress-------------------------------------#

#--------------------------R API-------------------------------------#
#--------------------------Deep Learning-----------------------------#
library(deepnet);library(darch);library(FCNN4R);
library(ForecastCombinations);
library(EbayesThresh);library(HMM);library(markovchain)
#------------------------Data Processing-----------------------------#
library(XML);library(xml2);library(XML2R);library(seqinr);
library(methods);library(dplyr);
library(Matrix);library(pracma);library(expm)

library(recommenderlab);library(readr);library(GA)


#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#
W<-data.frame();X<-data.frame();Y<-data.frame();Z<-data.frame();
#---------------------------------------------------------------------#
#------------------------------Functions------------------------------#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#

#------------Item-based Collaborative Filtering Recommender Model-----#
Recommender.Model.Jazz.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Recommender.Model.Jazz.1<-Recommender.Model.Jazz.1("1")
test.Recommender.Model.Jazz.1

#------------Composition Models I------------------#
Composition.Model.Jazz.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Composition.Model.Jazz.1<-Composition.Model.Jazz.1("1")
test.Composition.Model.Jazz.1
#------------Composition Models II-----------------#
Composition.Model.Jazz.2<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Composition.Model.Jazz.2<-Composition.Model.Jazz.1("1")
test.Composition.Model.Jazz.2
#------------Composition Models III----------------#

Composition.Model.Jazz.3<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Composition.Model.Jazz.3<-Composition.Model.Jazz.1("1")
test.Composition.Model.Jazz.3

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
Table.4.TeX<-xtable::xtable(T

#---------------------------------------------------------------------#
#------------------------------Figures--------------------------------#
#---------------------------------------------------------------------#

#----------Figure 1------#
#----------Figure 2------#
#----------Figure 3------#
#----------Figure 4------#
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
