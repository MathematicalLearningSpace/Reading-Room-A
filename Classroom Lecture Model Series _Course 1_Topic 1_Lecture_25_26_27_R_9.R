#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#
#--------------------Work In Progress-------------------------------------#
#------------------------------R API----------------------------------#
library(deSolve);library(ReacTran);library(rootSolve);
library(fda);library(phaseR)
library(pracma);library(GA);library(igraph)

library(tseries);library(costat);library(locits);library(wbsts);
library(forecast);library(tsoutliers);library(jmotif)
library(TSclust);library(TSMining);library(ggplot2);
library(tsDyn);library(tseriesChaos);library(yuima);library(DescTools)
library(xtable);library(PearsonDS);library(fitdistrplus);library(psych)

#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#
A.Drug.Resistance<-c(CDX2, MUC2, REG4, CDH17, MDR1, SHH)
B.Genomic.Instability<-c(p53, p21, BAX, p48, GADD45, BAK, POLK)
C.Tumor.Progression<-c(Retinoic.Acid, RAR.Beta, RXR)
D.Intestinal.Metaplasia<-c(DV1, GSK.3Beta, Beta.Catenin, Axin,APC, CK1.alpha,GBP)
E.Dysplasia.Path.1<-c(EFG, ERBB2, SHC, GRB2, SOS, RAS, RAF, MEK, ERK.1)
F.Dysplasia.Path.2<-c(PI3K, PIP3, AKT, mTOR, p53, S6K, BCL2)
G.Dysplasia.Path.3<-c(TGF.Beta,TGF.BetaRI,TGF.BetaRII, SMAD.2, SMAD.4, p15, p21)
H.Normal.Gastic.Muscosa.1<-c(HGF, c.MET, GRB2, SOS, RAS, RAF, MEK, ERK.1)
I.Normal.Gastic.Muscosa.Survival.Path.1<-c(FGF, FGFR2, GAB1, PI3K, PIP3, AKT, mTOR, GSK.3Beta)

#------------------------------Noise----------------------------------#
noise.white<- noise(kind = c("white"))
noise.pink<-noise(kind = c("pink"))
noise.5.4<-noise(kind = c("power"), alpha=(5/4))
noise.7.4<-noise(kind=c("power"),alpha=(7/4))
noise.red<-noise(kind = c("red"))
noise.8.4<-noise(kind=c("power"),alpha=(8/4))

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
Model.1<-function (Time.Sequence, variables.intitial, Params.1)
  {
    with(as.list(c(Params.1, variables.intitial)), 
         {
           #------------Equations Designed in the Classroom-----#
           d.1.x.dt.1<-a0
           d.1.y.dt.1<-a1
           d.1.z.dt.1<-a2
           d.1.w.dt.1<-a3
   res <- c(d.1.x.dt.1,
           d.1.y.dt.1,
           d.1.z.dt.1,
           d.1.w.dt.1)
  list(res)
         })
}                                                
#----------------Parameters-----------#   
initial.vector<-c(x=1,y=0,z=0,w=0)
Params.1<-c(a0<-1,a1<-1,a2<-1,a3<-1,a4<-1,a5<-1)                                                    
Time.Sequence<- seq(0, 48, by = 48/120) 

#---------------------------------------------------------------------#
#------------------------------Solutions------------------------------#
#---------------------------------------------------------------------#
model.system.solution.1 <- ode(y = initial.vector, times = Time.Sequence, 
func = Model.1, parms = Parms.1)
#---------------------------------------------------------------------#
#------------------------------Analysis-------------------------------#
#---------------------------------------------------------------------#
options(digits = 3)
model.system.1.summary<-summary(model.system.solution.1)
#-------------------Topology------------------------------------------#
Topology.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Topology.1<-Topology.1("1")
test.Topology.1

#-------------------Equilibrium---------------------------------------#
Equilibrium.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Equilibrium.1<-Equilibrium.1("1")
test.Equilibrium.1

#-------------------Stability-----------------------------------------#
Stability.1<-function(X,Point.Vector.Collection)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();Z<-NULL
  Path.Matrix.df<-as.data.frame(X)
  Jacobian.List<-list()
  Hessian.List<-list()
  for(j in 1:nrow(Point.Vector.Collection))
  {
    Z[j]<-Point.Vector.Collection[j,]
    J[[j]]<-jacobian(Z[j])
    H[[j]]<-hessian(Z[j])
  }
  output<-list()
  output$X<-X
  output$Z<-Z
  output$J<-J
  output$H<-H
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Stability.1<-Stability.1("1")
test.Stability.1

#---------------------------------------------------------------------#
#------------------------------Tables---------------------------------#
#---------------------------------------------------------------------#

#---------Table 1-----------------------#
#---------Table 2-----------------------#
#---------Table 3-----------------------#
#---------Table 4-----------------------#
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

#---------Figure 1---------------------#
png(file = stringr::str_c("Figures//Example_",1,"_Figure_",0,".png"))
Figure.1<-matplot(model.system.solution.1[,2:5], type = "l", lty = 1, lwd = c(2, 1, 1,1),
                  col = c("darkred", "darkblue", "darkgreen","red"),
                  xlab = "Time [min]", ylab= "Y",
                  main = "Solution Example")
grid()
dev.off() 
#---------Figure 2---------------------#
png(file = stringr::str_c("Figures//Example_",1,"_Figure_",1,".png"))
Figure.1<-scatterplotMatrix(model.system.solution.1[,2:5])
dev.off() 
#---------Figure 3---------------------#

#---------Figure 4---------------------#

#---------Figure 5---------------------#

#---------Figure 6---------------------#
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

