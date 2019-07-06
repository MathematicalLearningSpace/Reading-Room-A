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


#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#------------------------------Functions------------------------------#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#

#-----------Composition.Jazz.1------------------#

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

#-----------Composition.Jazz.2------------------#
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
test.Composition.Model.Jazz.2<-Composition.Model.Jazz.2("1")
test.Composition.Model.Jazz.2
#-----------Composition.Romantic.1--------------#
Composition.Model.Romantic.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Composition.Model.Romantic.1<-Composition.Model.Romantic.1("1")
test.Composition.Model.Romantic.1
#---------------------------------------------------------------------#
#------------------------------Analysis-------------------------------#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#------------------------------Tables---------------------------------#
#---------------------------------------------------------------------#

#----------Table 1------#
#----------Table 2------#
#----------Table 3------#
#----------Table 4------#

#---------------------------------------------------------------------#
#------------------------------Figures--------------------------------#
#---------------------------------------------------------------------#

#----------Figure 1------#
#----------Figure 2------#
#----------Figure 3------#
#----------Figure 4------#

