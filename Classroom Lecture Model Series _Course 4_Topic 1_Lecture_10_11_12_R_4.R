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
#-----------------------Signal Processing----------------------------#
library(audio);library(tuneR);library(dtw);library(wavelets);
library(dtwclust);library(TSclust);library(TSdist)
library(wavethresh);library(waveslim);library(wavemulcor)
library(PearsonDS);library(xtable);library(psych);library(adwave);library(biwavelet)
library(wmtsa);library(wavelets);library(wavemulcor)
library(PearsonDS);library(xtable);library(psych);library(adwave);library(BootWPTOS)
library(liftLRD);library(mwaved);library(mvLSW);library(unbalhaar)
#----------------------Parallel Processing---------------------------#
library(parallel);library(microbenchmark)

#---------------------Data Files in Music XML------------------------#
album.1<-c("A.xml","B.xml","C.xml","D.xml","E.xml","F.xml","G.xml","H.xml","I.xml","J.xml")
W<-data.frame();X<-data.frame();Y<-data.frame();Z<-data.frame();
pearson.N<-512
pIIIpars <- list(shape=1, location=1, scale=1)
error.pearson.3<-rpearsonIII(pearson.N,params=pIIIpars)
experimental.data.1<-window(error.pearson.3, end=256)

cusps.1 <- function(x) -0.5*abs(x-1)^0.5 -0.5* abs(x-2)^0.5 + 0.5*x + 3.5
cusps.2 <- function(x) -1*abs(x-1)^0.75 -0.5* abs(x-8)^0.75 + 1*x + 3     

x <- seq(0, 10, length=512)
y.1 <- splus2R::signalSeries(cusps.1(x), x)
y.2 <- splus2R::signalSeries(cusps.2(x), x)
y.3<- makeHeaviSine(pearson.N)
#-------------------------------------------------#
#---Individual Tracks of Music Compositions---#
#-------------------------------------------------#
music.track.1<-read_xml(album.1[1])
music.track.2<-read_xml(album.1[2])
music.track.3<-read_xml(album.1[3])
music.track.4<-read_xml(album.1[4])
music.track.5<-read_xml(album.1[5])
music.track.6<-read_xml(album.1[6])
music.track.7<-read_xml(album.1[7])
music.track.8<-read_xml(album.1[8])
music.track.9<-read_xml(album.1[9])
music.track.10<-read_xml(album.1[10])

#---------------------------------------------------------------------#
#------------------------------Functions------------------------------#
#---------------------------------------------------------------------#
#---------------------------------------Transforms---------------------------------------
x.cwt.1 <- wavCWT( as.ts(experimental.data.1), wavelet="gaussian1")
x.cwt.2 <- wavCWT( as.ts(experimental.data.1), wavelet="gaussian2")
x.cwt.3 <- wavCWT( as.ts(experimental.data.1), wavelet="Haar")
x.cwt.4 <- wavCWT( as.ts(experimental.data.1), wavelet="morlet")
#-------------------------------------Trees----------------------------------------------
W.tree.1 <- wavCWTTree(x.cwt.1)
W.tree.2 <- wavCWTTree(x.cwt.2)
W.tree.3 <- wavCWTTree(x.cwt.3)
W.tree.4 <- wavCWTTree(x.cwt.4)
#-----------------------------------Holder Exponent---------------------------------------
holder.1 <- holderSpectrum(W.tree.1)
holder.2 <- holderSpectrum(W.tree.2)
holder.3 <- holderSpectrum(W.tree.3)
holder.4 <- holderSpectrum(W.tree.4)
#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#

Model.Wavelet.Analysis.1<-function(X,visualization=TRUE)
 {
Z<-as.matrix(rnorm(10000),nrow=10000,ncol=1)
scales <- seq(1, 64, 3)
W.1.Coefs <- MassSpecWavelet::cwt(X[1:10000], scales=scales, wavelet='mexh')
W.1.Coefs.localMax <- MassSpecWavelet::getLocalMaximumCWT(W.1.Coefs)
if(visualization){Figure.1<-MassSpecWavelet::plotLocalMax(W.1.Coefs.localMax)}

 output<-list()
 output$W.1.Coefs<-W.1.Coefs
 output$W.1.Coefs.localMax<-W.1.Coefs.LocalMax
 return(output)
}
X.1<-as.matrix(rnorm(10000),nrow=10000,ncol=1)
test.Model.Wavelet.Analysis.1<-Model.Wavelet.Analysis.1(X.1)
test.Model.Wavelet.Analysis.1

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
png(file = stringr::str_c('Figures//Example_',1,'_Figure_',1,'.png'))
Figure.1<-plot(W.tree.1$`3`$extrema, type="l", lty=1, col="black",xlab="Observations",ylab="Extrema")
lines(W.tree.2$`3`$extrema,lty=2, col="green")
lines(W.tree.3$`3`$extrema,lty=3, col="blue")
lines(W.tree.4$`3`$extrema,lty=4, col="red")
rug(W.tree.1$`3`$extrema, side=4, col="black")
rug(W.tree.2$`3`$extrema, side=4, col="green")
rug(W.tree.3$`3`$extrema, side=4, col="blue")
rug(W.tree.4$`3`$extrema, side=4, col="red")
legend("topright",c("Branch 5_1=Gaussian 1","Branch 5_2=Gaussian 2","Branch 5_3=Haar","Branch 5_4=Morlet"),
       inset = .01,col=c("black","green","blue","red"),lwd=2,cex=0.8)
dev.off()

png(file = stringr::str_c('Figures//Example_',1,'_Figure_',2,'.png'))
op <- par(mfrow = c(3, 3),mar=c(1,1,1,1))
Figure.2<-plot(W.tree.1, xlab="(a)")
plot(W.tree.2, xlab="(b)")
plot(W.tree.3, xlab="(c)")
plot(W.tree.4, xlab="(d)")
par(op)
dev.off()

png(file = stringr::str_c('Figures//Example_',1,'_Figure_',3,'.png'))
Figure.2<-densityBy(as.data.frame(as.data.frame(y.1)),xlab="CUSP Expression 1", ylab="Simulated",col="green")
dev.off()

png(file = stringr::str_c('Figures//Example_',1,'_Figure_',4,'.png'))
Figure.3<-plot(holder.1$exponent, type="l",ylim=c(-4,1),lty=1,xlab="Observations", ylab="Exponent",col="black")
lines(holder.2$exponent,lty=2, col="green")
lines(holder.3$exponent,lty=3, col="blue")
lines(holder.4$exponent,lty=4, col="red")
rug( holder.1$exponent, side=4, col="black")
rug( holder.2$exponent, side=4, col="green")
rug( holder.3$exponent, side=4, col="blue")
rug( holder.4$exponent, side=4, col="red")
legend("bottomright",c("Holder_1=Gaussian 1","Holder_2=Gaussian 2","Holder_3=Haar","Holder_4=Morlet"),inset = .01,col=c("black","green","blue","red"),lwd=2,cex=0.8)
dev.off()

png(file = stringr::str_c('Figures//Example_',1,'_Figure_',5,'.png'))
Figure.4<-plot(experimental.data.1, type="l",lty=1,ylim<-c(-1,5),xlab="Observations", ylab="value",col="black")
lines(y.1,lty=2, col="green")
lines(y.2,lty=3, col="blue")
lines(y.3,lty=4,col="red")
rug(experimental.data.1, side=4, col="black")
rug(y.1, side=4, col="green")
rug(y.2, side=4, col="blue")
rug(y.3, side=4, col="red")
legend("topleft",c("Pearson III","CUSP 1","CUSP 2","HeaviSine"),inset = .01,col=c("black","green","blue","red"),lwd=2,cex=0.8)
dev.off()
#----------Figure 1-----------------#
#----------Figure 2-----------------#
#----------Figure 3-----------------#
#----------Figure 4-----------------#
png(file = stringr::str_c('Figures//Example_',1,'_Figure_',6,'.png'))
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
