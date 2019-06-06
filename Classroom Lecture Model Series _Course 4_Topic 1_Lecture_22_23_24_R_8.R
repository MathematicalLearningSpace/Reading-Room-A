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
library(Matrix);library(pracma);library(expm);

library(seewave);library(tuneR);library(soundecology);library(glogis);
library(nsRFA)

#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#
Music_Notes <- read_csv("Music_Notes.txt")
Music.Collection.1.1<-readWave("Composition1.wav")
Music.Collection.1.2<-readWave("Composition2.wav")
Music.Collection.1.3<-readWave("Composition3.wav")
Music.Collection.1.4<-readWave("Composition4.wav")
Music.Collection.1.5<-readWave("Composition5.wav")
Music.Collection.1.6<-readWave("Composition6.wav")
Music.Collection.1.7<-readWave("Composition7.wav")
Music.Collection.1.8<-readWave("Composition8.wav")
Music.Collection.1.9<-readWave("Composition9.wav")
Music.Collection.1.10<-readWave("Composition10.wav")

#---------------------------------------------------------------------#
#------------------------------Functions------------------------------#
#---------------------------------------------------------------------#

Acoustic.Index.Prototype<-function(X)
{
  #-----------Work in Progress for the Classroom----------#
  Table.1<-data.frame()
  X.pspectrum<- powspec(X@left, X@samp.rate)
  X.duration<- length(X@left/X@samp.rate)
  X.ACI<-ACI(X, nbwindows=(X/5))
  X.ACI.B<-acoustic_complexity(X,j=5) # j=5
  X.AR<-AR(X)
  X.ADI<-acoustic_diversity(X)
  X.AEI<-acoustic_evenness(X)
  X.BAI<-bioacoustic_index(X)
  X.Entropy<-H(X@left,X@samp.rate)
  X.spec1<-spec(X,f=22050,at=0.2,plot=TRUE)
  ## Shannon spectral entropy
  X.Entropy.Shannon.Spectral<-sh(X.spec1)
  ## Renyi spectral entropy
  X.Entropy.Renyi.Spectral.2<-sh(X.spec1, alpha=2)
  X.Entropy.Renyi.Spectral.3<-sh(X.spec1, alpha=3)
  X.roughness<-roughness(X.spec1)
  X.rugosity<-rugo(X@left/max(X@left))
  Music.Index<-cbind(X.duration,X.ACI.A,X.AR,X.Entropy)
  output<-list()
  output$Music.Index<-Music.Index
  output$Table.1<-Table.1
  return(output)
}


#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#

#-----------Markov Models I--------------#

#-----------Semi-Markov Models I---------#

#-----------Music Machine Learning-------#

#---------------------------------------------------------------------#
#------------------------------Analysis-------------------------------#
#---------------------------------------------------------------------#

#------------------Music Collection Statistics---------#


#---------------------------------------------------------------------#
#------------------------------Tables---------------------------------#
#---------------------------------------------------------------------#

#--------------Tables 1-----------#
#--------------Tables 2-----------#
#--------------Tables 3-----------#
#--------------Tables 4-----------#


#---------------------------------------------------------------------#
#------------------------------Figures--------------------------------#
#---------------------------------------------------------------------#

#-------------Figures 1------------#
#-------------Figures 2------------#
#-------------Figures 3------------#
#-------------Figures 4------------#
#-------------Figures 5------------#
#-------------Figures 6------------#
#-------------Figures 7------------#
