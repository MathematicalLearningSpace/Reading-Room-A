#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#

#--------------------Work In Progress-----------August 2019---------------#
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
Cancer.Application.1<-function(X)
 {
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
 system.equation.model.test<-function(times, variables.intitial.test, parameters.test)
{
  with(as.list(c(parameters.test, variables.intitial.test)), 
       {

#--------------Group 1---------------------
d.1.X1.d.t.1<--a11*X1 + a12*X2
#--------------Group 2---------------------
d.1.X2.d.t.1<-a21*X2
#--------------Group 3---------------------
d.1.X3.d.t.1<-a72*X2 - a81*X1
res <- c(d.1.X1.d.t.1,
         d.1.X2.d.t.1,
         d.1.X3.d.t.1)
list(res)
})
}

 
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 
 }
test.Cancer.Application.1<-Cancer.Application.1("1")
test.Cancer.Application.1

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
  require(bio3d);
  Table.1.df<-data.frame();Table.2.df<-data.frame();Table.3.df<-data.frame();
  setwd("Model/Protein Models")
  protein.model.files.pdb<-list.files();protein.model.list<-list();N.1<-length(protein.model.files.pdb);
  for(i in 1:3){protein.model.list[[i]]<-read.pdb(protein.model.files.pdb[i])}
  bio <- biounit(protein.model.list[[3]]);names(bio)
  #--------------------------------------------------------------------------------#
  #-----------------Sequence Modeling Global, Local, Motif Optimization------------#
  #--------------------------------------------------------------------------------#
  sequence.Pattern.Recognition<-function(Pattern.1,Y,Visualization=FALSE)
  {
    aa.seq<-pdbseq(Y)
    Search.Result<-motif.find(Pattern.1, aa.seq)
    #---------------------------------------------------------------------------#
    #----------------Alignment--------------------------------------------------#
    #---------------------------------------------------------------------------#
    Protein.blast <- blast.pdb()
    Table.1<-Protein.blast$hit.tbl
    if(Visualization){top.hits <- plot(Protein.blast)}
    Table.2<-top.hits$hits
    output<-list()
    output$Pattern.1<-Pattern.1
    output$Model.1<-Y
    output$Table.1<-Table.1
    output$Table.2<-Table.2
    return(output)
  }
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
Analysis.Model.1<-function(X,Visualization=FALSE)
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
test.Analysis.Model.1<-Analysis.Model.1("1",TRUE)
test.Analysis.Model.1
#---------------------------------------------------------------------#
#------------------------------Tables---------------------------------#
#---------------------------------------------------------------------#

#-----------Table 1-------------------#
Table.1.TeX<-xtable::xtable(Table.1.df)
#-----------Table 2-------------------#
Table.2.TeX<-xtable::xtable(Table.2.df)
#-----------Table 3-------------------#
Table.3.TeX<-xtable::xtable(Table.3.df)
#-----------Table 4-------------------#
Table.4.TeX<-xtable::xtable(Table.4.df)
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
Cancer.Application.2<-function(X)
 {
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
 system.equation.model.test<-function(times, variables.intitial.test, parameters.test)
{
  with(as.list(c(parameters.test, variables.intitial.test)), 
       {

#--------------Group 1---------------------
d.1.X1.d.t.1<--a11*X1 + a12*X2
#--------------Group 2---------------------
d.1.X2.d.t.1<-a21*X2
#--------------Group 3---------------------
d.1.X3.d.t.1<-a72*X2 - a81*X1
res <- c(d.1.X1.d.t.1,
         d.1.X2.d.t.1,
         d.1.X3.d.t.1)
list(res)
})
}
 
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 
 }
test.Cancer.Application.2<-Cancer.Application.2("1")
test.Cancer.Application.2
#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#

#--------------Review Notes-------#
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
Analysis.Model.2<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Analysis.Model.2<-Analysis.Model.2("1")
test.Analysis.Model.2
#---------------------------------------------------------------------#
#------------------------------Tables---------------------------------#
#---------------------------------------------------------------------#

#-----------Table 1-------------------#
#-----------Table 2-------------------#
#-----------Table 3-------------------#
#-----------Table 4-------------------#
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
Cancer.Application.3<-function(X)
 {
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
 system.equation.model.test<-function(times, variables.intitial.test, parameters.test)
{
  with(as.list(c(parameters.test, variables.intitial.test)), 
       {

#--------------Group 1---------------------
d.1.X1.d.t.1<--a11*X1 + a12*X2
#--------------Group 2---------------------
d.1.X2.d.t.1<-a21*X2
#--------------Group 3---------------------
d.1.X3.d.t.1<-a72*X2 - a81*X1
res <- c(d.1.X1.d.t.1,
         d.1.X2.d.t.1,
         d.1.X3.d.t.1)
list(res)
})
}
 
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 
 }
test.Cancer.Application.3<-Cancer.Application.3("1")
test.Cancer.Application.3
#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#

#--------------Review Notes-------#
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
Analysis.Model.3<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Analysis.Model.3<-Analysis.Model.3("1")
test.Analysis.Model.3
#---------------------------------------------------------------------#
#------------------------------Tables---------------------------------#
#---------------------------------------------------------------------#

#-----------Table 1-------------------#
#-----------Table 2-------------------#
#-----------Table 3-------------------#
#-----------Table 4-------------------#

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

