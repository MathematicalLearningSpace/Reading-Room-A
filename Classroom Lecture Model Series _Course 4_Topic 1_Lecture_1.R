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
#-----------Table 2-------------------#
#-----------Table 3-------------------#
#-----------Table 4-------------------#


#----------------------------------Figures-------------------------------------------#

#-----------Figure 1-------------------#
#-----------Figure 2-------------------#
#-----------Figure 3-------------------#
#-----------Figure 4-------------------#
        
        
#----------------------------------Discussion----------------------------------------#
