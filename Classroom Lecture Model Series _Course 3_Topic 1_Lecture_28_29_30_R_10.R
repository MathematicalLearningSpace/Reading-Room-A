#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#

#--------------------Work In Progress-------------------------------------#
#------------------------------R API----------------------------------#
library(gdata);library(bio3d);library(igraph);library(sna);library(ips);
library(phangorn);library(proteomics)
library(dcGOR);library(MDplot);library(UniProt.ws);
library(circlize);library(BioPhysConnectoR);library(protr)
library(seqinr);library(Biostrings);library(Peptides);

library(easyPubMed);library(readr);library(CHNOSZ);
library(stringr);library(seqinr);library(seqLogo);library(msa);library(ape);
library(dtw);library(dtwclust);library(odseq);library(rphast)
library(plyr)

#------Scientific  Visualization-----#
library(corrplot);library(plot3D);library(scatterplot3d);library(rgl)

#-------------------------Application I-------------------------------#

#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#
setwd("A Model")
  Table.1.df<-as.data.frame(Table.1);Table.2.df<-as.data.frame(Table.2);Table.3.df<-as.data.frame(Table.3);Table.4.df<-as.data.frame(Table.4)
  Table.5.df<-as.data.frame(Table.5);Table.6.df<-as.data.frame(Table.6);Table.7.df<-as.data.frame(Table.7);Table.8.df<-as.data.frame(Table.8)
  Table.9.df<-as.data.frame(Table.9);Table.10.df<-as.data.frame(Table.10)

#--------------Review Notes-------#
#-------------------------------------------------------#
#------------Article Keywords for Search----------------#
#-------------------------------------------------------#

Article.keywords<-c("A","B","C","D");Search.Results<-list();
require(easyPubMed);
for(i in 1:length(Article.keywords))
  {
    search.string.1<-stringr::str_c(Article.keywords[i],"  AND (2018[PDAT]:2019[PDAT])")
    Search.Results[[i]]<-fetch_pubmed_data(get_pubmed_ids(search.string.1))
    saveXML(Search.Results[[i]],file=stringr::str_c("Search_Results_X_",i,".xml"))
    for(j in 1:10000){}
  }
setwd("Model/Abstracts")
require(XML);require(xml2);require(seqinr)
X.Models.Korpus.Files<-list.files()
X.Models.Korpus.XML<-list();X.Models.Korpus.Titles<-list();X.Models.Korpus.Abstracts<-list();
X.Models.Dictionary<-data.frame()
  
for(i in 1:length(X.Models.Korpus.Files))
{
    X.Models.Korpus.XML[[i]]<-read_xml(X.Models.Korpus.Files[i])
    X.Models.Korpus.Titles[[i]]<-xml_text(xml_find_all(X.Models.Korpus.XML[[i]], "//ArticleTitle"))
    X.Models.Korpus.Abstracts[[i]]<-xml_text(xml_find_all(X.Models.Korpus.XML[[i]], ".//AbstractText"),trim=TRUE) 
 }

#---------------------------------------------------------------------#
#------------------------------Functions------------------------------#
#---------------------------------------------------------------------#
Topic.Model.1<-function(X)
{
require(tm);require(topicmodels);require(slam);require(lattice);
}
test.Topic.Model.1<-Topic.Model.1(X.Models.Korpus.Abstracts)
test.Topic.Model.1
#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#

Model.1<-function(X,Visualization=FALSE)
{
  Table.1.df<-data.frame();Table.2.df<-data.frame();Table.3.df<-data.frame();
  
  if(Visualization){plot(X)};
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
} 
test.Model.1<-Model.1("1",TRUE)
test.Model.1
#---------------------------------------------------------------------#
#------------------------------Analysis-------------------------------#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#------------------------------Tables---------------------------------#
#---------------------------------------------------------------------#

#-----------Table 1-------------------#
#-----------Table 2-------------------#
#-----------Table 3-------------------#
#-----------Table 4-------------------#

#---------------------------------------------------------------------#
#------------------------------Figures--------------------------------#
#---------------------------------------------------------------------#

#-----------Figure 1-------------------#
require(circlize)
for(i in 1:3)
{
png(file = stringr::str_c('Figures//Example_',i,'_Figure_',2,'.png'))
circos.initializeWithIdeogram(plotType = c("labels", "axis"))
X = generateRandomBed(nr = 10^2, nc = 2^2)
color.choices = colorRamp2(c(-1, 0, 1), c("blue", "red", "green"))
circos.genomicHeatmap(X, color.choices, side = "inside", border = "black")
circos.genomicHeatmap(X, color.choices, side = "outside", line_col = as.numeric(factor(X[[1]])))
dev.off()
}
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

#-------------------------Application II------------------------------#

#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#

#--------------Review Notes-------#

#---------------------------------------------------------------------#
#------------------------------Functions------------------------------#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#
Model.2<-function(X,Visualization=FALSE)
{
  Table.1.df<-data.frame();Table.2.df<-data.frame();Table.3.df<-data.frame();
  
  if(Visualization){plot(X)};
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
} 
test.Model.2<-Model.2("1",TRUE)
test.Model.2
#---------------------------------------------------------------------#
#------------------------------Analysis-------------------------------#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#------------------------------Tables---------------------------------#
#---------------------------------------------------------------------#

#-----------Table 1-------------------#
#-----------Table 2-------------------#
#-----------Table 3-------------------#
#-----------Table 4-------------------#

#---------------------------------------------------------------------#
#------------------------------Figures--------------------------------#
#---------------------------------------------------------------------#

#-----------Figure 1-------------------#
#-----------Figure 2-------------------#
#-----------Figure 3-------------------#
#-----------Figure 4-------------------#
png(file = stringr::str_c('Figures//Example_',2,'_Figure_',1,'.png'))
op <- par(mfrow = c(2,2),mar=c(3,3,3,3))
hist(W, main="Title 1",xlab="X Value")
legend("topright", legend = paste(seq(1:7),LETTERS[1:7]),lty = 1, cex = .8, y.intersp = 1)
hist(X, main="Title 2",xlab="Note Value")
legend("topright", legend = paste(seq(1:7),LETTERS[1:7]),lty = 1, cex = .8, y.intersp = 1)
hist(Y, main="Title 3",xlab="Note Value")
legend("topright", legend = paste(seq(1:7),LETTERS[1:7]),lty = 1, cex = .8, y.intersp = 1)
hist(Z, main="Title 4",xlab="Note Value")
legend("topright", legend = paste(seq(1:7),LETTERS[1:7]),lty = 1, cex = .8, y.intersp = 1)
par(op)
dev.off()
#-------------------------Application III-----------------------------#

#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#

#--------------Review Notes-------#

#---------------------------------------------------------------------#
#------------------------------Functions------------------------------#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#
Model.3<-function(X,Visualization=FALSE)
{
  Table.1.df<-data.frame();Table.2.df<-data.frame();Table.3.df<-data.frame();
  
  if(Visualization){plot(X)};
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
} 
test.Model.3<-Model.3("1",TRUE)
test.Model.3
#---------------------------------------------------------------------#
#------------------------------Analysis-------------------------------#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#------------------------------Tables---------------------------------#
#---------------------------------------------------------------------#

#-----------Table 1-------------------#
#-----------Table 2-------------------#
#-----------Table 3-------------------#
#-----------Table 4-------------------#


#---------------------------------------------------------------------#
#------------------------------Figures--------------------------------#
#---------------------------------------------------------------------#

#-----------Figure 1-------------------#
#-----------Figure 2-------------------#
#-----------Figure 3-------------------#
#-----------Figure 4-------------------#
png(file = stringr::str_c('Figures//Example_',3,'_Figure_',1,'.png'))
op <- par(mfrow = c(2,2),mar=c(3,3,3,3))
hist(W, main="Title 1",xlab="X Value")
legend("topright", legend = paste(seq(1:7),LETTERS[1:7]),lty = 1, cex = .8, y.intersp = 1)
hist(X, main="Title 2",xlab="Note Value")
legend("topright", legend = paste(seq(1:7),LETTERS[1:7]),lty = 1, cex = .8, y.intersp = 1)
hist(Y, main="Title 3",xlab="Note Value")
legend("topright", legend = paste(seq(1:7),LETTERS[1:7]),lty = 1, cex = .8, y.intersp = 1)
hist(Z, main="Title 4",xlab="Note Value")
legend("topright", legend = paste(seq(1:7),LETTERS[1:7]),lty = 1, cex = .8, y.intersp = 1)
par(op)
dev.off()

