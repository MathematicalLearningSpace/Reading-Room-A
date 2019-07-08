#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#

#------------------------------------------------------------------------------------#
#----------------------------------R Source Files------------------------------------#
#------------------------------------------------------------------------------------#

#-------------------------R API-----------------#
library(KEGG.db);library(KEGGgraph);library(KEGGprofile);
library(KEGGREST);library(rentrez)
library(xtable);library(stringi);library(readr);
library(readxl);library(Matrix);library(igraph);library(visNetwork);
library(rcellminer);library(rcellminerData);library(bioCancer);
library(trena);library(ctSGE);library(anamiR) 
library(srnadiff);library(RmiR);library(CancerSubtypes);library(GRENITS);
library(compEpiTools);
library(trigger);library(CVE);library(qgraph)
library(sgnesR);library(rcdk);library(sqldf)
library(tm);library(stringi);library(stringr);library(utils);
library(CHNOSZ);library(PearsonDS);library(xtable)
library(bio3d);library(Peptides);library(jsonlite);
library(rjson);library(seqinr)
library(easyPubMed);library(bio3d);library(readr);library(CHNOSZ);
library(stringr);library(Peptides);library(Biostrings)
library(seqinr);library(seqLogo);library(msa);library(ape);
library(dtw);library(dtwclust);library(odseq);library(rphast)
library(plyr)
#----------------------Parallel Processing and Benchmarking--------------------------#
library(parallel);library(microbenchmark);

#------------------------------------------------------------------------------------#
#----------------------------------Data----------------------------------------------#
#------------------------------------------------------------------------------------#
Table.1.df<-as.data.frame(Table.1);Table.2.df<-as.data.frame(Table.2);Table.3.df<-as.data.frame(Table.3);
Table.4.df<-as.data.frame(Table.4)Table.5.df<-as.data.frame(Table.5);Table.6.df<-as.data.frame(Table.6);
Table.7.df<-as.data.frame(Table.7);Table.8.df<-as.data.frame(Table.8);Table.9.df<-as.data.frame(Table.9);
Table.10.df<-as.data.frame(Table.10)

Oncogenes.Sample.1 <- c("ABL1", "ALK", "BRAF", "CCND1", "CCND3", "CCNE1", "CCNE2", 
                       "CDC25A", "EGFR", "ERBB2", "EZH2", "FOS", "FOXL2", "HRAS", 
                       "IDH1", "IDH2", "JAK2", "KIT", "KRAS", "MET", "MOS", "MYC", 
                       "NRAS", "PDGFB", "PDGFRA", "PIK3CA", "PIK3CB", "PIK3CD", 
                       "PIK3CG", "PIM1", "PML", "RAF1", "REL", "RET", "SRC", "STK11", 
                       "TP63")
Oncogenes.Sample.2 <- c("WNT10B", "WNT4", "WNT2B", "WNT9A", "WNT3", "WNT5A", 
                       "WNT5B", "WNT10A", "WNT11", "WNT2", "WNT1", "WNT7B", "WISP1", 
                       "WNT8B", "WNT7A", "WNT16", "WISP2", "WISP3", "FZD5", "FZD1")


KEGG.organisms<-keggList("organism")
KEGG.Brite<-keggList("brite")
KEGG.organisms.df<-as.data.frame(KEGG.organisms)

#--------------Cell Lines--------------------------#
drug.Act <- exprs(getAct(rcellminerData::drugData))
mol.Data <- getMolDataMatrices()

PlotTypes <- c("drug", "cop", "exp", "xai", "exo", "mut", "mir", "pro", "mda")
cop.1 <- "cop";exp.1 <- "exp";xai.1 <- "xai";exo.1<- "exo";mut.1 <- "mut";mir.1 <- "mir";pro.1 <- "pro"
#------Filters------
cop.data <- mol.Data[[cop.1]];exp.data <- mol.Data[[exp.1]];xai.data <- mol.Data[[xai.1]];exo.data <- mol.Data[[exo.1]]
mut.data <- mol.Data[[mut.1]]; mir.data <- mol.Data[[mir.1]];pro.data <- mol.Data[[pro.1]]

#--------------Gene Ontology-----------------------#
Gene.Ontology.Enrichment.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Gene.Ontology.Enrichment.Model.1<-Gene.Ontology.Enrichment..Model.1("1")
test.Gene.Ontology.Enrichment.Model.1
#------------------------------------------------------------------------------------#
#----------------------------------Transformations-----------------------------------#
#------------------------------------------------------------------------------------#
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
#------------------------------------------------------------------------------------#
#----------------------------------User Defined Modules and Functions----------------#
#------------------------------------------------------------------------------------#
 Cancer.types<-c("Colorectal cancer","Pancreatic cancer","Hepatocellular carcinoma","Gastric cancer","Glioma","Thyroid cancer",
  "Acute myeloid leukemia","Chronic myeloid leukemia","Basal cell carcinoma","Melanoma","Renal cell carcinoma",
  "Bladder cancer","Prostate cancer","Endometrial cancer","Breast cancer","Small cell lung cancer","Non-small cell lung cancer")

#------------------------------------------------------------------------------------------------#
#-------------------------------Cancer Models----------------------------------------------------#
#------------------------------------------------------------------------------------------------#
  setwd("Cancer Models")
  Cancer.Model.1.Name<-c("Colorectal cancer");Cancer.Model.1<-xml2::read_xml("hsa05210.xml")}
  Cancer.Model.2.Name<-c("Pancreatic cancer");Cancer.Model.2<-xml2::read_xml("hsa05212.xml")
  Cancer.Model.3.Name<-c("Hepatocellular carcinoma");Cancer.Model.3<-xml2::read_xml("hsa05225.xml")
  Cancer.Model.4.Name<-c("Gastric cancer");Cancer.Model.4<-xml2::read_xml("hsa05226.xml")
  Cancer.Model.5.Name<-c("Glioma");Cancer.Model.5<-xml2::read_xml("hsa05214.xml")
  Cancer.Model.6.Name<-c("Thyroid cancer");Cancer.Model.6<-xml2::read_xml("hsa05216.xml")
  Cancer.Model.7.Name<-c("Acute myeloid leukemia");Cancer.Model.7<-xml2::read_xml("hsa05221.xml")
  Cancer.Model.8.Name<-c("Chronic myeloid leukemia");Cancer.Model.8<-xml2::read_xml("hsa05220.xml")
  Cancer.Model.9.Name<-c("Basal cell carcinoma");Cancer.Model.9<-xml2::read_xml("hsa05217.xml")
  Cancer.Model.10.Name<-c("Melanoma");Cancer.Model.10<-xml2::read_xml("hsa05218.xml")
  Cancer.Model.11.Name<-c("Renal cell carcinoma");Cancer.Model.11<-xml2::read_xml("hsa05211.xml")
  Cancer.Model.12.Name<-c("Bladder cancer");Cancer.Model.12<-xml2::read_xml("hsa05219.xml")
  Cancer.Model.13.Name<-c("Prostate cancer");Cancer.Model.13<-xml2::read_xml("hsa05215.xml")
  Cancer.Model.14.Name<-c("Endometrial cancer");Cancer.Model.14<-xml2::read_xml("hsa05213.xml")
  Cancer.Model.15.Name<-c("Breast cancer");Cancer.Model.15<-xml2::read_xml("hsa05224.xml")
  Cancer.Model.16.Name<-c("Small cell lung cancer");Cancer.Model.16<-xml2::read_xml("hsa05222.xml")
  Cancer.Model.17.Name<-c("Non-small cell lung cancer");Cancer.Model.17<-xml2::read_xml("hsa05223.xml")
#------------------------------------------------------------------------------------#
#----------------------------------Network Designs-----------------------------------#
#------------------------------------------------------------------------------------#
Network.Design.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Network.Design.Model.1<-Network.Design.Model.1("1")
test.Network.Design.Model.1
#------------------------------------------------------------------------------------#
#----------------------------------Equation Systems----------------------------------#
#------------------------------------------------------------------------------------#
Equation.System.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Equation.System.1<-Equation.System.1("1")
test.Equation.System.1
#------------------------------------------------------------------------------------#
#----------------------------------Parameter Tables----------------------------------#
#------------------------------------------------------------------------------------#
Matrix.Parameter.1<-matrix(1,nrow=10,ncol=10)
Matrix.Parameter.2<-matrix(1,nrow=10,ncol=10)
Matrix.Parameter.3<-matrix(1,nrow=10,ncol=10)
#------------------------------------------------------------------------------------#
#----------------------------------Network Analysis----------------------------------#
#------------------------------------------------------------------------------------#
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
#------------------------------------------------------------------------------------#
#----------------------------------Optimization--------------------------------------#
#------------------------------------------------------------------------------------#
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
test.Optimization.Model.1<-Analysis.Model.1("1")
test.Optimization.Model.1
#------------------------------------------------------------------------------------#
#----------------------------------Natural Language Description----------------------#
#------------------------------------------------------------------------------------#
NLP.Description.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.NLP.Description.Model.1<-NLP.Description.Model.1("1")
test.NLP.Description.Model.1
#------------------------------------------------------------------------------------#
#----------------------------------Tables--------------------------------------------#
#------------------------------------------------------------------------------------#

#----------Table 1 Group A---------#

#----------Table 2 Group B---------#

#----------Table 3 Group C---------#

#----------Table 4 Group D---------#
#----------Table 1------#
Table.1.TeX<-xtable::xtable(Table.1.df)
#----------Table 2------#
Table.2.TeX<-xtable::xtable(Table.2.df)
#----------Table 3------#
Table.3.TeX<-xtable::xtable(Table.3.df)
#----------Table 4------#
Table.4.TeX<-xtable::xtable(Table.4.df)
#----------------------------------Figures-------------------------------------------#

#----------Figure 1 Group A---------#

#----------Figure 1 Group B---------#

#----------Figure 1 Group C---------#

#----------Figure 1 Group D---------#
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

#------------------------------------------------------------------------------------#  
#----------------------------------Discussion----------------------------------------#
#------------------------------------------------------------------------------------#

#----------------References----------------------#
