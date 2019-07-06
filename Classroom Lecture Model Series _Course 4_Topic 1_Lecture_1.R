#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#

#----------------------------------R Source Files------------------------------------#

#-----------------------R API------
library(tm);library(stringi);library(stringr);library(utils);
library(jsonlite);library(rjson);library(xtable)
#----------------------------------Data----------------------------------------------#

#----------------------------------Transformations-----------------------------------#

#----------------------------------User Defined Modules and Functions----------------#
Acoustic.Complexity.Index<-function(X,Window.Width)
  {
  Acoustic.Complexity.Local<-0
  Acoustic.Complexity.Global<-0
 
  output<-list()
 	output$Acoustic.Complexity.Local<-Acoustic.Complexity.Local
 	output$Acoustic.Complexity.Local<-Acoustic.Complexity.Global
  return(output)
  }
test.Acoustic.Complexity.Index<-Acoustic.Complexity.Index(X,3)
test.Acoustic.Complexity.Index

#-----------------------------------Music Model--------------------------------------#
Music.Model.1<-function(W,X,Y,Z,Visualization=TRUE)
{
setwd("Music Model/Music Sequences/XML Data")
Music.Composition.1<-"";Music.Composition.2<-"";Music.Composition.3<-"";Music.Composition.4<-"";Music.Composition.5<-"";
Music.Composition.6<-"";Music.Composition.7<-"";Music.Composition.8<-"";Music.Composition.9<-"";Music.Composition.10<-"";
Notes.By.Measure<-list();Notes.Duration<-list();Notes.By.Step<-list();
Notes.By.Octave<-list();Notes.By.Pitch<-list();Model.1<-list(); Model.Chain.1<-list();notes.series<-list();
 Music.files<-list.files()
 Music.Composition.1<-read_xml(Music.files[1])
 Music.Composition.2<-read_xml(Music.files[2])
 Music.Composition.3<-read_xml(Music.files[3])
 Music.Composition.4<-read_xml(Music.files[4])
 Music.Composition.5<-read_xml(Music.files[5])
 Music.Composition.6<-read_xml(Music.files[6])
 Music.Composition.7<-read_xml(Music.files[7])
 Music.Composition.8<-read_xml(Music.files[8])
 Music.Composition.9<-read_xml(Music.files[9])
 Music.Composition.10<-read_xml(Music.files[10])
 
 music.track.summary<-function(X)
{
   #--------------Musical Objects-------------------------------------#
  music.categories<-c("chord","pitch","step","alter","octave","duration","voice","type","stem","staff")
  Notes.By.Measure<-list();Notes.Duration<-list();Notes.By.Step<-list();Notes.By.Octave<-list();Notes.By.Pitch<-list()
  Model.1<-list();Model.Chain.1<-list();notes.series<-list()
  
   #-----------------------Presentation and Organization--------------#
  output<-list()
  output$X<-X
  return(output)
}
 
 
 
 #-----------------------Presentation and Organization--------------#
output<-list()
output$Music.Composition.1<-Music.Composition.1
output$Music.Composition.2<-Music.Composition.2
output$Music.Composition.3<-Music.Composition.3
output$Music.Composition.4<-Music.Composition.4
output$Music.Composition.5<-Music.Composition.5
output$Music.Composition.6<-Music.Composition.6
output$Music.Composition.7<-Music.Composition.7
output$Music.Composition.8<-Music.Composition.8
output$Music.Composition.9<-Music.Composition.9
output$Music.Composition.10<-Music.Composition.9
  return(output)
}
test.Music.Model.1<-Music.Model.1("1","1","1","1",FALSE)
test.Music.Model.1
#----------------------------------Network Designs-----------------------------------#

#----------------------------------Equation Systems----------------------------------#

#----------------------------------Parameter Tables----------------------------------#

#----------------------------------Network Analysis----------------------------------#

#----------------------------------Optimization--------------------------------------#

#----------------------------------Natural Language Description----------------------#

#----------------------------------Tables--------------------------------------------#


#-----------Table 1-------------------#
Table.1.TeX<-xtable::xtable(Table.1.df)
#-----------Table 2-------------------#
Table.2.TeX<-xtable::xtable(Table.2.df)
#-----------Table 3-------------------#
Table.3.TeX<-xtable::xtable(Table.3.df)
#-----------Table 4-------------------#
Table.4.TeX<-xtable::xtable(Table.4.df)

#----------------------------------Figures-------------------------------------------#

#-----------Figure 1-------------------#
#-----------Figure 2-------------------#
#-----------Figure 3-------------------#
#-----------Figure 4-------------------#
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
