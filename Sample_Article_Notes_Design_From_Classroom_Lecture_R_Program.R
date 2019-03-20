#-----------R Program for InClass Lecture For Sample Journal Article Notes---------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------------#
library(pracma);library(xtable);library(deSolve);library(ReacTran);library(rootSolve);library(fda);library(igraph);library(boot)
library(sampling);library(PearsonDS);library(KEGGgraph);library(keggorthology);library(KEGGREST);library(Rgraphviz)
library(igraph);library(jsonlite);library(tsDyn);library(tseriesChaos);library(wavelets);library(waveslim);library(readr);library(randtests)
library(changepoint);library(car);library(tuneR);library(diffEq);library(phaseR)
#----------------------------------------------------------------------------------------------------------------------------------#

Table.1.df<-data.frame()
Table.2.df<-data.frame()
Table.3.df<-data.frame()
Table.4.df<-data.frame()
Table.5.df<-data.frame()
Table.6.df<-data.frame()

Gene.Study<-c("LD","ATRp","p53","p53p","Mdm2c","Mdm2cp","Mdm2n","Aktp","PIP3","PTEN","DDB2","p21tot",
              "p21CE","CycEtot","E2F1","Rbp","Bax","CytoC","Apaf1","Apops","Casp9","Casp3")
Gene.Study.Literature<-c("ATR","p53","p53","Mdm2","Akt","PIP3","PTEN","DDB2","p21",
                         "CycE","E2F1","Rbp","Bax","CytoC","Apaf1","Apops","Casp9","Casp3")
#-------------------Noise Processes------------------------
noise.white<-noise(kind = c("white"))
noise.pink<-noise(kind = c("pink"))
noise.red<-noise(kind = c("red"))
#--------------------Power law-----------------------------
noise.4.4<-noise(kind=c("power"),alpha=(4/4))
noise.5.4<-noise(kind=c("power"),alpha=(5/4))
noise.6.4<-noise(kind=c("power"),alpha=(6/4))
noise.7.4<-noise(kind=c("power"),alpha=(7/4))
noise.8.4<-noise(kind=c("power"),alpha=(8/4))
noise.test<-c(white=noise.white,pink=noise.pink,noise.1=noise.5.4,red.2=noise.7.4,red=noise.red,red.3=noise.8.4)
noise.test.parameters <-c(k.noise.1=0,k.noise.2=0,k.noise.3=0,k.noise.4=0,k.noise.5=0,k.noise.6=0,k.noise.7=0,k.noise.8=0,
                          k.noise.9=0,k.noise.10=0,k.noise.11=0,k.noise.12=0,k.noise.13=0,k.noise.14=0,k.noise.15=0,k.noise.16=0,
                          k.noise.17=0,k.noise.18=0,k.noise.19=0,k.noise.20=0,k.noise.21=0,k.noise.22=0)

#----------------------------------------------------Algorithm Development-------------------------------------------
Differential.Equation.Algorithms = c("lsoda", "lsode", "lsodes", "lsodar", "vode", "daspk", "euler", "rk4", "ode23", "ode45", "radau", "bdf", "bdf_d", "adams", "impAdams", "impAdams_d", 
           "iteration")

Algorithm.1<-function(X)
{
  output<-list()
  return(output)
}
  
Algorithm.2<-function(X)
{
   
  output<-list()
  return(output)
}
  
Algorithm.3<-function(X)
{
  
  output<-list()
  return(output)
}
