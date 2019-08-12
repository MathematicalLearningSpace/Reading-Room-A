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
W<-data.frame();X<-data.frame();Y<-data.frame();Z<-data.frame();
#---------Introduction
Track.1.midi<-readMidi("Composition_1.mid")
#---------Part I
Track.2.midi<-readMidi("Composition_2.mid")
Track.3.midi<-readMidi("Composition_3.mid")
Track.4.midi<-readMidi("Composition_4.mid")
#---------Part II
Track.5.midi<-readMidi("Composition_5.mid")
Track.6.midi<-readMidi("Composition_6.mid")
Track.7.midi<-readMidi("Composition_7.mid")
#---------Part III
Track.8.midi<-readMidi("Composition_8.mid")
Track.9.midi<-readMidi("Composition_9.mid")
#----------Conclusion
Track.10.midi<-readMidi("Composition_10.mid")
#----------------------------------Midi Processing
Track.1.midi.Notes<-getMidiNotes(Track.1.midi)
Track.2.midi.Notes<-getMidiNotes(Track.2.midi)
Track.3.midi.Notes<-getMidiNotes(Track.3.midi)
Track.4.midi.Notes<-getMidiNotes(Track.4.midi)
Track.5.midi.Notes<-getMidiNotes(Track.5.midi)
Track.6.midi.Notes<-getMidiNotes(Track.6.midi)
Track.7.midi.Notes<-getMidiNotes(Track.7.midi)
Track.8.midi.Notes<-getMidiNotes(Track.8.midi)
Track.9.midi.Notes<-getMidiNotes(Track.9.midi)
Track.10.midi.Notes<-getMidiNotes(Track.10.midi)
#---------------------------------------------------------------------#
#------------------------------Functions------------------------------#
#---------------------------------------------------------------------#
Transformations.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Transformations.1<-Transformations.1("1")
test.Transformations.1
#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#
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
