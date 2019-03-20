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
#-------------------Noise Processes------------------------#
noise.white<-noise(kind = c("white"))
noise.pink<-noise(kind = c("pink"))
noise.red<-noise(kind = c("red"))
#--------------------Power law-----------------------------#
noise.4.4<-noise(kind=c("power"),alpha=(4/4))
noise.5.4<-noise(kind=c("power"),alpha=(5/4))
noise.6.4<-noise(kind=c("power"),alpha=(6/4))
noise.7.4<-noise(kind=c("power"),alpha=(7/4))
noise.8.4<-noise(kind=c("power"),alpha=(8/4))
noise.test<-c(white=noise.white,pink=noise.pink,noise.1=noise.5.4,red.2=noise.7.4,red=noise.red,red.3=noise.8.4)
noise.test.parameters <-c(k.noise.1=0,k.noise.2=0,k.noise.3=0,k.noise.4=0,k.noise.5=0,k.noise.6=0,k.noise.7=0,k.noise.8=0,
                          k.noise.9=0,k.noise.10=0,k.noise.11=0,k.noise.12=0,k.noise.13=0,k.noise.14=0,k.noise.15=0,k.noise.16=0,
                          k.noise.17=0,k.noise.18=0,k.noise.19=0,k.noise.20=0,k.noise.21=0,k.noise.22=0)

#----------------------------------------------------Algorithm Development-------------------------------------------#
Differential.Equation.Algorithms = c("lsoda", "lsode", "lsodes", "lsodar", "vode", "daspk", "euler", "rk4", 
                                     "ode23", "ode45", "radau", "bdf", "bdf_d", "adams", "impAdams", "impAdams_d", "iteration")

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
#------------------------------------------------Parameter Matrix-----------------------------------------------------#

#------------------------------------------------Equation Systems------------------------------------------------------#


#------------------------------------------------DE Solvers, Optimizers and Algorithms--------------------------------#


#------------------------------------------------Network--------------------------------------------------------------#
Network.A<-c("UV","ATR","ATRp","RPA")
Network.B<-c("p21","p21CE","DDb2","Rb","Rbp","Cyclin E","EF21")
Network.C<-c("p53","p53p","Pten","Mdm2n","Mdm2c","Mdm2cp","Akt","Aktp","PIP2","PIP3")
Network.D<-c("BAX","CytoC","Apops","Apaf-1","CASP9","Procasp9","CASP3","Procasp3","PARP1","PARP3")
Networks<-c(Network.A,Network.B,Network.C,Network.D)
#------------------------------------------------Tables--------------------------------------------------------------#


#------------------------------------------------Figures-------------------------------------------------------------#

png(file = stringr::str_join("Figures//Example_",1,"_Figure_","1",".png"))
Figure.1<-plot(system.equation.model.test.solution.1[,2],xlab="T",ylab="Value")
lines(system.equation.model.test.solution.1[,3],col=2)
lines(system.equation.model.test.solution.1[,4],col=3)
legend("topleft", legend = Network.A, col = seq(Network.A),lty = 1, cex = .8)
dev.off()
png(file = stringr::str_join("Figures//Example_",1,"_Figure_","2",".png"))
Figure.2<-plot(system.equation.model.test.solution.1[,5],xlab="T",ylab="Value")
lines(system.equation.model.test.solution.1[,6],col=2)
lines(system.equation.model.test.solution.1[,7],col=3)
lines(system.equation.model.test.solution.1[,8],col=4)
lines(system.equation.model.test.solution.1[,9],col=5)
lines(system.equation.model.test.solution.1[,10],col=6)
lines(system.equation.model.test.solution.1[,11],col=7)
legend("topleft", legend = Network.B, col = seq(Network.B),lty = 1, cex = .8)
dev.off()

png(file = stringr::str_join("Figures//Example_",1,"_Figure_","3",".png"))
Figure.3<-plot(system.equation.model.test.solution.1[,14],xlab="T",ylab="Value")
lines(system.equation.model.test.solution.1[,15],col=2)
lines(system.equation.model.test.solution.1[,16],col=3)
lines(system.equation.model.test.solution.1[,17],col=4)
lines(system.equation.model.test.solution.1[,18],col=5)
lines(system.equation.model.test.solution.1[,19],col=6)
lines(system.equation.model.test.solution.1[,20],col=7)
lines(system.equation.model.test.solution.1[,21],col=8)
lines(system.equation.model.test.solution.1[,22],col=9)
lines(system.equation.model.test.solution.1[,23],col=10)
legend("topleft", legend = Network.C, col = seq(Network.C),lty = 1, cex = .8)
dev.off()
png(file = stringr::str_join("Figures//Example_",1,"_Figure_","4",".png"))
Figure.4<-plot(system.equation.model.test.solution.1[,26],xlab="T",ylab="Value")
lines(system.equation.model.test.solution.1[,27],col=2)
lines(system.equation.model.test.solution.1[,28],col=3)
lines(system.equation.model.test.solution.1[,29],col=4)
lines(system.equation.model.test.solution.1[,30],col=5)
lines(system.equation.model.test.solution.1[,31],col=6)
legend("topleft", legend = Network.D, col = seq(Network.D),lty = 1, cex = .8)
dev.off()

#-------------------------------------Network A------------------------
#Network.A<-c("UV","ATR","ATRp","RPA")
png(file = stringr::str_join("Figures//Example_",1,"_Figure_","5",".png"))
Figure.1.A<-scatterplot.matrix(system.equation.model.test.solution.1[,2:4])
dev.off()
#-------------------------------------Network B--------------------------
#Network.B<-c("p21","p21CE","DDb2","Rb","Rbp","Cyclin E","EF21")
png(file = stringr::str_join("Figures//Example_",1,"_Figure_","6",".png"))
Figure.1.B<-scatterplot.matrix(system.equation.model.test.solution.1[,5:13])
dev.off()
#-------------------------------------Network C--------------------------
#Network.C<-c("p53","p53p","Pten","Mdm2n","Mdm2c","Mdm2cp","Akt","Aktp","PIP2","PIP3")
png(file = stringr::str_join("Figures//Example_",1,"_Figure_","7",".png"))
Figure.1.C<-scatterplot.matrix(system.equation.model.test.solution.1[,14:25])
dev.off()
#-------------------------------------Network D---------------------------
#Network.D<-c("BAX","CytoC","Apops","Apaf-1","CASP9","Procasp9","CASP3","Procasp3","PARP1","PARP3")
png(file = stringr::str_join("Figures//Example_",1,"_Figure_","8",".png"))
Figure.1.D<-scatterplot.matrix(system.equation.model.test.solution.1[,26:31])
dev.off()
#------------------------------------Histograms and QQ Plots-----------------------------------------------------------#
for(i in 2:31)
{
op <- par(mfrow = c(2, 2))
png(file = stringr::str_join("Figures//Example_",1,"_Figure_",i+10,".png"))
Figure<-hist(system.equation.model.test.solution.1[,i],main=stringr::str_c(Networks[i]),xlab=stringr::str_c(Networks[i]),
                 breaks = 12, col = "green", border = "blue")
dev.off()
}
for(i in 2:31)
{
png(file = stringr::str_join("Figures//Example_",1,"_Figure_",i+40,".png"))
Figure<-qqplot(system.equation.model.test.solution.1[,i], qchisq(ppoints(system.equation.model.test.solution.1[,i]), df = 4),
                   main=stringr::str_c(Networks[i]),xlab=stringr::str_c(Networks[i])); 
abline(0, 1, col = 2, lty = 2)
dev.off()
}


#------------------------------------------------References----------------------------------------------------------#
