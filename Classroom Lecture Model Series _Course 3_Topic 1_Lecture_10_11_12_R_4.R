#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#

#--------------------Work In Progress-------------------------------------#
#----------------------R API----------------------------------------------#
library(bio3d);library(xtable);library(Peptides);library(stringi);
library(gdata);library(igraph);library(sna);library(ips);
library(phangorn);library(proteomics)
library(dcGOR);library(MDplot);library(UniProt.ws);
library(circlize);library(BioPhysConnectoR);library(protr)
library(seqinr);library(Biostrings);library(Peptides);

library(easyPubMed);library(bio3d);library(readr);library(CHNOSZ);
library(stringr);library(Biostrings)
library(seqinr);library(seqLogo);library(msa);library(ape);
library(dtw);library(dtwclust);library(odseq);library(rphast)
library(plyr)

WDR<-c(WDR1, WDR5,WDR10, WDR12, WDR13, WDR16, 
       WDR17,WDR18, WDR19, WDR20, WDR21A, 
       WDR21C, WDR22, WDR23, WDR24, WDR25, 
       WDR26, WDR27, WDR3, WDR31, WDR32, WDR33, 
       WDR34, WDR35, WDR36, WDR37, WDR38, WDR4, 
       WDR40A, WDR40B, WDR40C, WDR41, WDR42A, WDR42B, WDR43, 
       WDR44, WDR46, WDR47, WDR48, WDR49, WDR5, WDR51A, WDR51B, 
       WDR52, WDR53, WDR54, WDR55, WDR57, WDR59, WDR5B, 
       WDR6, WDR60, WDR61, WDR62, WDR63, WDR64, 
       WDR65, WDR66, WDR67, WDR68, WDR69, WDR7, 
       WDR70, WDR72, WDR73, WDR74, WDR75, WDR76, 
       WDR77, WDR78, WDR79, WDR8, WDR81, WDR82, 
       WDR85, WDR86, WDR88, WDR89, WDR90, WDR91,WDR92)
#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#

#------------Example 1-----------------------#
WDR5.pdb<-read.pdb("WDR5_HUMAN_1.pdb")

#---------------------------------------------------------------------#
#------------------------------Functions------------------------------#
#---------------------------------------------------------------------#
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
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#

#------------------WDR Models-------------------------------#
WDR.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(outputWDR
 }
test.WDR.Model.1<-WDR.Model.1("1")
test.WDR.Model.1

#------------------G Quadruplexes---------------------------#
#-----------------------------G Quadruplexes-------------------------------------------
Intramolecular.G.Quadruplexes.Basket_xml<-fetch_pubmed_data(get_pubmed_ids("Intramolecular G Quadruplexes Basket"))
Intramolecular.G.Quadruplexes.Chair_xml<-fetch_pubmed_data(get_pubmed_ids("Intramolecular G Quadruplexes Chair"))
Intramolecular.G.Quadruplexes.Propeller_xml<-fetch_pubmed_data(get_pubmed_ids("Intramolecular G Quadruplexes Propeller"))
Intermolecular.G.Quadruplexes.Hairpin.dimmer_xml<-fetch_pubmed_data(get_pubmed_ids("Intermolecular G-Quadruplexes Hairpin dimmer"))
Intermolecular.G.Quadruplexes.Tetrameric_xml<-fetch_pubmed_data(get_pubmed_ids("Intermolecular G-Quadruplexes Tetrameric")
                                                                
G.QUAD.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  setwd("Immune System Model/Protein Models")
  protein.model.files.pdb<-list.files();protein.model.list<-list();K<-length(protein.model.files.pdb)
  #----------------------------------------------------------------------------------------------------------------------------#
  #---------------------------------------PDB Files----------------------------------------------------------------------------#
  #----------------------------------------------------------------------------------------------------------------------------#
  for(i in 1:k){protein.model.list[[i]]<-read.pdb(protein.model.files.pdb[i])}
  bio.Unit.1 <- biounit(protein.model.list[[3]]);names(bio.Unit.1)	
  #---------------------------------------#
  #-------Chain Sequence Modeling---------#
  #---------------------------------------#
  
  #-----------------------------------------------------#
  #-----------------------Residue Selection-------------#
  #-----------------------------------------------------#
 
  #---------------------------------------#
  #-------Network Dependency Modeling-----#
  #---------------------------------------#	

  #---------------------------------------#
  #-------Deformation energies------------# 
  #---------------------------------------#	

  #---------------------------------------#
  #----------Fluctuations-----------------# 
  #---------------------------------------#
	
  #---------------------------------------#
  #-------Torsion Analysis By Chain-------#
  #---------------------------------------#

  #---------------------------------------#
  #---cross-correlations------------------#
  #---------------------------------------#
  
  #---------------------------------------------#
  #---------------Graph Theory Analysis---------#
#---------------------------------------------#
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.G.QUAD.Model.1<-G.QUAD.Model.1("1")
test.G.QUAD.Model.1
#-----------------Heat Shock Proteins-----------------------#
HSP.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.HSP.Model.1<-HSP.Model.1("1")
test.HSP.Model.1
#-----------------Protein To Protein Interaction Models-----#
PIP.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.PIP.Model.1<-PIP.Model.1("1")
test.PIP.Model.1
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

#-----------Table 1-------------------#
#-----------Table 2-------------------#
#-----------Table 3-------------------#
#-----------Table 4-------------------#
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
#-----------Figure 1------------------#
#-----------Figure 2------------------#
#-----------Figure 3------------------#
#-----------Figure 4------------------#

