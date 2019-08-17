#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#

#--------------------Work In Progress-------------------------------------#
library(TSA);library(tseries);library(costat);library(locits);library(wbsts);
library(forecast);library(tsoutliers);library(jmotif)
library(TSclust);library(TSMining);library(ggplot2);library(tsDyn);
library(tseriesChaos);library(yuima);library(DescTools)
library(xtable);library(PearsonDS);library(fitdistrplus);library(psych);
library(BNPTSclust)
library(boot);library(sampling);library(RandomFields);
library(dtwclust);library(dtw);library(TSMining)
library(HDMD);library(xtable);library(pdc);library(lattice);library(fracdiff);library(tseriesEntropy)
#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#
W<-data.frame();X<-data.frame();Y<-data.frame();Z<-data.frame();
TS.grp.1 <- replicate(10, arima.sim(n = 48, list(ar = c(0.8, -0.4), ma = c(-0.22, 0.2)),sd = sqrt(0.1)) )
TS.grp.2 <- replicate(10, arima.sim(n = 48, list(ar = c(-0.7, 0.1), ma = c(0.9, -0.1)),sd = sqrt(0.09)) )
#--------------------long tailed distribution
TS.grp.3 <- replicate(10, arima.sim(n = 48, list(ar = c(-0.5, 0.25), ma = c(0.5, 0.25)),rand.gen = function(n, ...) sqrt(0.1) * rt(n, df = 10 ))
#--------------------fractional 
TS.grp.4<-replicate(10,fracdiff.sim(48, ar = .2, ma = .4, d = .3))
#-----------------------------Parameter Model-------------------------#
params.1<-c(a11=0.1,a12=0.1,a13=0.1,a14=0.1,a15=0.1,a16=0.1,
            a21=0.1,a22=0.1,a23=0.1,a24=0.1,a25=0.1,a26=0.1,
            a31=0.1,a32=0.1,a33=0.1,a34=0.1,a35=0.1,a36=0.1,
            a41=0.1,a42=0.1,a43=0.1,a44=0.1,a45=0.1,a46=0.1,
            a51=0.1,a52=0.1,a53=0.1,a54=0.1,a55=0.1,a56=0.1,
            a61=0.1,a62=0.1,a63=0.1,a64=0.1,a65=0.1,a66=0.1)

params.epsilon.1<-c(epsilon1,epsilon2,epsilon3)

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

#----------------------Transformations--------------------------------#
Transformation.1<-function(X,grp1,grp2,grp3,grp4)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  #--------------------------------Minimum Entropic Heuristic-----------------------------
            heuristic.1 <-  entropyHeuristic(grp1[,1] )
            heuristic.delay.1 <-  entropyHeuristic(grp1[,1], t.min=1, t.max=6 )
            heuristic.2 <-  entropyHeuristic(grp2[,1] )
            heuristic.delay.2 <-  entropyHeuristic(grp2[,1], t.min=1, t.max=6 )
            heuristic.3 <-  entropyHeuristic(grp3[,1] )
            heuristic.delay.3 <-  entropyHeuristic(grp3[,1], t.min=1, t.max=6 )
            heuristic.4 <-  entropyHeuristic(grp4[,1]$series )
            heuristic.delay.4 <-  entropyHeuristic(grp4[,1]$series, t.min=1, t.max=6 )
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Transformation.1<-Transformation.1("1",TS.grp.1,TS.grp.2,TS.grp.3,TS.grp.4)
test.Transformation.1
#-----SAX Representations for Motifs----------------------------------#
Motif.SAX.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Motif.SAX.1<-Motif.SAX.1("1")
test.Motif.SAX.1
#-----Hierarchical clustering of time series for Classification-------#
Classifier.Cluster.Hierarchical.1<-function(X,grp1,grp2,grp3,grp4,visualization=TRUE)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
 #----------Permutation Distribution Clustering for time series---------------------------
Group.X <- cbind(grp1,grp2,grp3)
D <- pdcDist(Group.X,3)
clustering <- pdclust(Group.X)
x <- codebook(grp2[,1],m=4)
y <- codebook(grp4[,1]$series,m=4)
D.hellinger<-hellingerDistance(x,y)
D.hellinger.squared<-squaredHellingerDistance(x,y)
Divergence<-symmetricAlphaDivergence(x,y)
  #----------Hierarchical cluster analysis for time series---------------------------
hc.1<-hclust(dist(Group.X))
memb.1 <- cutree(hc.1, k = 20)
memb.2 <- cutree(hc.1, k = 15)
memb.3 <- cutree(hc.1, k = 10)
cent <- NULL
for(k in 1:20){cent <- rbind(cent, colMeans(Group.X[memb.1 == k, , drop = FALSE]))}
hc.2 <- hclust(dist(cent)^2, method = "cen", members = table(memb.1))
cent <- NULL
for(k in 1:15){cent <- rbind(cent, colMeans(Group.X[memb.2 == k, , drop = FALSE]))}
hc.3 <- hclust(dist(cent)^2, method = "cen", members = table(memb.2))
cent <- NULL
for(k in 1:10){cent <- rbind(cent, colMeans(Group.X[memb.3 == k, , drop = FALSE]))}
hc.4 <- hclust(dist(cent)^2, method = "cen", members = table(memb.3))
#----------KMeans clustering------------------------------------------------------
clustering.Kmeans<-kmeans(Group.X,3)
metric.1<- function(x) sum(scale(x, scale = FALSE)^2)
clustering.Kmeans.fitted <- fitted(clustering.Kmeans)
head(clustering.Kmeans.fitted)
clustering.Kmeans.resid<- Group.X - fitted(clustering.Kmeans)
diagnostics.df<-cbind(clustering.Kmeans[c("betweenss", "tot.withinss", "totss")],
      c(metric.1(clustering.Kmeans.fitted), metric.ss(clustering.Kmeans.resid), metric.1(Group.X)))
  Table.1<-xtable(diagnostics.df)
 if(visualization)
    {  
#-----Additional Table Designs provided by students--------------------------------------                                 
#---------------Figures to be presented in the Classroom----------------------------------------------
par(mfrow = c(2,1))
Figure.1<-plot(grp1[,1], type="l", lty=1,col="black",ylab="Simulated Value", xlab="Temporal Position")
lines(grp1[,2], lty=2,col="red")
lines(grp1[,3], lty=3,col="blue")
lines(grp2[,1], lty=4,col="green")
lines(grp2[,2], lty=5,col="yellow")
lines(grp2[,3], lty=6,col="orange")

Figure.2<-plot(grp3[,1], type="l", lty=1,col="black",ylab="Simulated Value", xlab="Temporal Position")
lines(grp3[,2], lty=2,col="red")
lines(grp3[,3], lty=3,col="blue")
lines(grp3[,4], lty=4,col="green")
lines(grp3[,5], lty=5,col="yellow")
lines(grp3[,6], lty=6,col="orange")
#----------------------------------Figure Groups 1------------------------------------
par(mfrow = c(2,2))
Figure.3A<-plot(heuristic.1)
Figure.3B<-plot(heuristic.2)
Figure.3C<-plot(heuristic.3)
Figure.3D<-plot(heuristic.4)
#----------------------------------Figure Groups 2------------------------------------			
par(mfrow = c(2,2))
Figure.4A<-plot(heuristic.delay.1)
Figure.4B<-plot(heuristic.delay.2)
Figure.4C<-plot(heuristic.delay.3)
Figure.4D<-plot(heuristic.delay.4)
Figure.5<-levelplot(as.matrix(D), col.regions=grey.colors(20,start=0.1, end=0.95))
#----------------------------------Figure Groups 4------------------------------------	
par(mfrow = c(1,2))
Figure.6<-plot(clustering, labels=1:30, cols=c(rep("red",10),rep("blue",10),rep("green",10)))
Figure.6A<-plot(clustering, cols=c(rep("red",10),rep("blue",10),("green",10)))
#----------------------------------Figure Groups 4------------------------------------	
par(mfrow = c(2,2))
Figure.7A<-plot(hc.1)
Figure.7B<-plot(hc.2)
Figure.7C<-plot(hc.3)
Figure.7D<-plot(hc.4)
Figure.8<-plot(Group.X, col = clustering.Kmeans$cluster)
points(clustering.Kmeans$centers, col = 1:3, pch = 8, cex = 2)         
            } 
            
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Classifier.Cluster.Hierarchical.1<-Classifier.Cluster.Hierarchical.1("1",TS.grp.1,TS.grp.2,TS.grp.3,TS.grp.4)
test.Classifier.Cluster.Hierarchical..1
#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#

#----------------------Nonlinear Tests--------------------------------#
Tests.Nonlinear.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Tests.Nonlinear.1<-Tests.Nonlinear.1("1")
test.Tests.Nonlinear.1
#---------------------Stationarity-------------------------------------#
Tests.Stationarity.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Tests.Stationarity.1<-Tests.Stationarity.1("1")
test.Tests.Stationarity.1
#---------------------Prediction--------------------------------------#
Prediction.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Prediction.Model.1<-Prediction.Model.1("1")
test.Prediction.Model.1
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
#----------Table 1------#
Table.1.TeX<-xtable::xtable(Table.1.df)
#----------Table 2------#
Table.2.TeX<-xtable::xtable(Table.2.df)
#----------Table 3------#
Table.3.TeX<-xtable::xtable(Table.3.df)
#----------Table 4------#
Table.4.TeX<-xtable::xtable(Table.4.df)

#----------Table 5------#
Table.5.TeX<-xtable::xtable(Table.5.df)
#----------Table 6------#
Table.6.TeX<-xtable::xtable(Table.6.df)
#----------Table 7------#
Table.7.TeX<-xtable::xtable(Table.7.df)
#----------Table 8------#
Table.8.TeX<-xtable::xtable(Table.8.df)

#---------------------------------------------------------------------#
#------------------------------Figures--------------------------------#
#---------------------------------------------------------------------#

#---------Figure Group--A--------------------------------#
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

#----------Figure 1-------------#
#----------Figure 2-------------#
#----------Figure 3-------------#

#---------Figure Group--B--------------------------------#
png(file = stringr::str_c('Figures//Example_',1,'_Figure_',2,'.png'))
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

#----------Figure 1-------------#
#----------Figure 2-------------#
#----------Figure 3-------------#

#---------Figure Group--C--------------------------------#
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

#----------Figure 1-------------#
#----------Figure 2-------------#
#----------Figure 3-------------#




