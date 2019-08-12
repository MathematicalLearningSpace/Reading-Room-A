#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#
#--------------------Work In Progress-------------------------------------#
#-------------------------------------------------------------------------#
library(XML);library(xml2);library(XML2R);library(seqinr);
library(methods);library(dplyr);
library(Matrix);library(pracma);library(expm);
library(dplyr);library(audio);library(seewave)
library(sound);library(soundecology);library(soundgen)
#-----------------------------------Design Matrix---------------------
require(HMM);require(msm);require(ggraph);require(RJaCGH);
require(SemiMarkov);require(surveillance);require(depmixS4)
require(markovchain);require(XML);require(xml2)
#----------------------------------R Source Files------------------------------------#
#----------------------------------Data----------------------------------------------#
W<-data.frame();X<-data.frame();Y<-data.frame();Z<-data.frame();
 #---------------------Data Files in Music XML------------------------#
 album.1<-c("Composition 1.xml","Composition 2.xml","Composition 3.xml","Composition 4.xml","Composition 5.xml",
            "Composition 6.xml","Composition 7.xml","Composition 8.xml","Composition 9.xml","Composition 10.xml")
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
#----------------------------------Transformations-----------------------------------#
#----------------------------------User Defined Modules and Functions----------------#
Module.Music.1<-function(music.track)
{
  tag.remove <- function(tag) {
    return(gsub("<.*?>", "", tag))
  }
  notes.1<-tag.remove(xml_find_all(music.track, ".//note//step"))
  output<-list()
  #output$Model.1<-Model.1
  #output$Model.1A<-Model.2
  #output$Chain.1<-chain.1
  #output$Chain.1A<-chain.1A
  #output$Model.1.LogLikelihood<-Model.1$logLikelihood
  #output$Transition.Matrix<-Transition.Matrix
  #output$Table.1<-xtable::xtable(as.data.frame(Table.1))
  #output$Transition.Matrix.Table.1<-xtable::xtable(as.data.frame(Model.1$estimate@transitionMatrix))
  #output$Transition.Matrix.Table.1A<-xtable::xtable(as.data.frame(Model.1A$estimate@transitionMatrix))
  #output$Transition.Matrix.Table.2<-xtable::xtable(as.data.frame(Transition.Matrix))
  return(output)
}
test.Module.Music.1<-Module.Music.1(music.track.1)
test.Module.Music.1

Motif.Design.1<-function(X)
{

  #---------------------------------------------------------#
  #-----------------------Tables----------------------------#
  #---------------------------------------------------------#
  #Table.1<-table(notes)
  #-----------------------Presentation and Organization--------------
  output<-list()
  #output$Table.1.TeX<-xtable::xtable(as.data.frame(Table.1))
  #output$Transition.Matrix.Table.1<-xtable::xtable(as.data.frame(Model.Chain.1.Overall))
  #output$Table.1.df<-Table.1.df
  #output$Nbr.Measures<-Nbr.Measures
  #output$Nbr.Notes<-Nbr.Notes
  #output$Notes.By.Measure<-Notes.By.Measure
  #output$Notes.Duration<-Notes.Duration
  #output$Note.Duration.Time.Series<-unlist(Notes.Duration)
  #output$Notes.By.Step<-Notes.By.Step
  #output$Notes.By.Octave<-Notes.By.Octave
  #output$Notes.By.Octave<-Notes.By.Pitch
  #output$Notes.Time.Series<-as.numeric(unlist(Notes.By.Step))
  #output$Notes.Octave.Time.Series<-as.numeric(unlist(Notes.By.Octave))
  #output$Notes.Pitch.Time.Series<-unlist(Notes.By.Pitch)
  #output$Model.Chain.1<-Model.Chain.1
  #output$Piano.Keys.Time.Series<-piano.keys
  #output$Piano.Frequencies<-piano.frequencies
  return(output)
}
test.Motif.Design.1<-Motif.Design.1(music.track.1)
test.Motif.Design.1
#
#----------------------------------Network Designs-----------------------------------#
Network.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Network.Model.1<-Network.Model.1("1")
test.Network.Model.1
#----------------------------------Equation Systems----------------------------------#

#----------------------------------Parameter Tables----------------------------------#

#----------------------------------Network Analysis----------------------------------#
Network.Analysis.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Network.Analysis.Model.1<-Network.Analysis.Model.1("1")
test.Network.Analysis.Model.1
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
#----------------------------------Natural Language Description----------------------#
Language.Natural.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Language.Natural.Model.1<-Language.Natural.Model.1("1")
test.Language.Natural.Model.1
#----------------------------------Tables--------------------------------------------#
#----------Table 1------#
Table.1.df<-as.data.frame(matrix(0,nrow=3,ncol=3)
colnames(Table.1.df)<-Letters[1:3]

Table.1.TeX<-xtable::xtable(Table.1.df)
#----------Table 2------#
Table.2.TeX<-xtable::xtable(Table.2.df)
#----------Table 3------#
Table.3.TeX<-xtable::xtable(Table.3.df)
#----------Table 4------#
Table.4.TeX<-xtable::xtable(Table.4.df)
#----------------------------------Figures-------------------------------------------#
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
#----------------------------------Discussion----------------------------------------#
