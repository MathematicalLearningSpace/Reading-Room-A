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
library(dplyr);library(seewave)
library(sound);library(soundecology);library(soundgen)
library(dtwclust);library(TSclust);library(TSdist)
library(wavethresh);library(waveslim);library(wavemulcor)
library(PearsonDS);library(xtable);library(psych);library(adwave);
library(biwavelet)
#----------------------Parallel Processing---------------------------#
library(parallel);library(microbenchmark)

#---------------------Data Files in Music XML------------------------#
album.1<-c("a.xml","b.xml","c.xml","d.xml","e.xml","f.xml","g.xml","h.xml","i.xml","j.xml")
W<-data.frame();X<-data.frame();Y<-data.frame();Z<-data.frame();
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
  Z.dwt.1 <- wavelets::dwt(as.numeric(Z), filter = "la8")
  Z.dwt.2 <- wavelets::dwt(as.numeric(Z), filter = "d4")
  Z.dwt.3 <- wavelets::dwt(as.numeric(Z), filter = "haar")
  Z.dwt.4 <- wavelets::dwt(as.numeric(Z), filter = "c6")
  Z.dwt.list <- list( Z.dwt.1,  Z.dwt.2,  Z.dwt.3,  Z.dwt.4)
  mra.out <- wavelets::mra(as.numeric(Z), n.levels=3, boundary="reflection")
#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#
Note.Alphabet.Value<-c(A=0,B=2); Sequence.Note<-"A B A B A B"; Sequence.Duration<-c(1,1,1,1,1,1)
compose.wave.sine<-function(frequency,duration,tempo)
{
  amplitude<-1;sample.rate <- 44100;  fade <- seq(0, 1, 50 / sample.rate)
  wave <- amplitude*sin(seq(0, duration / tempo * 60, 1 / sample.rate) *frequency * 2 * pi)
  wave.sine.1<-wave*c(fade,rep(1, length(wave) - 2 * length(fade)), rev(fade))
  return(wave.sine.1)
}
Track.Synthetic.1<-Track.Synthetic.1 %>%
mutate(Track.octave = substring(Sequence.Note, nchar(Sequence.Note)) %>% ifelse(is.na(.), 4, .),
         Track.note = Note.Alphabet.Value[substr(Sequence.Note, 1, 1)],
         Track.note = Track.note + grepl("#", Sequence.Note) - grepl("b", Sequence.Note) + Track.octave * 12 + 12 * (Track.note < 3),
         freq = 2 ^ ((Track.note - 60) / 12) * 440)

Track.1.wav<-mapply(compose.wave.sine, Track.1$freq, Track.1$Seq.Duration,120) 

Composition.Model.Classical.1<-function(X,Design.Pattern,Nbr.Measures,Note.Model,Duration.Model,Instruments)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  k<-4
  Composition.Matrix.1.df<-as.data.frame(matrix(0,nrow=10^k,col=length(Model.Instruments)))
  
  for(i in 1:Nbr.Measures){}
  
  output<-list()
  output$Composition.1<-Composition.Matrix.1.df
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Composition.Model.Classical.1<-Composition.Model.Classical.1("1")
test.Composition.Model.Classical.1

Composition.Model.Jazz.1<-function(X,Design.Pattern,Nbr.Measures,Note.Model,Duration.Model,Instruments)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  for(i in 1:Nbr.Measures){}
  
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

#--------------Figure 1----------------------------#
png(file = stringr::str_c('Figures//Example_',1,'_Figure_',1,'.png'))
plot(Track.1.wav[1000:3000,1])
dev.off()
#--------------Figure 2----------------------------#
#--------------Figure 3----------------------------#
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


