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
library(PearsonDS);library(xtable)
library(rcdk);library(BioMedR);library(ChemmineR);
library(Matrix):library(fingerprint)
library(readr);library(leaps);library(caret); library(GA) 
library(ggplot2);library(kohonen);library(pROC)
#------Scientific  Visualization-----#
library(corrplot);library(plot3D);library(scatterplot3d);library(rgl)

library(easyPubMed);library(readr);library(CHNOSZ);
library(stringr);library(seqinr);library(seqLogo);library(msa);library(ape);
library(odseq);library(rphast)
library(plyr)

#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#

KEGG.Proteasome<-keggGet("br:ko03051")
KEGG.Proteasome.20s<-keggGet("path:pxb03050")
KEGG.Ubiquitin.System<-keggGet("br:ko04121")
KEGG.Ubiquitin.System.1<-keggGet("path:pxb04120")
KEGG.Ubiquitin.System.2<-keggGet("path:pxb03420")

#--------------Review Notes--------------#
Proteasome_xml <- fetch_pubmed_data(get_pubmed_ids("Proteasome"))

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

#----------------Proteasome Model - I----------#

#---------------------------------------------B: Standard 20s Proteasome AT----------------------------------

KEGG.Proteasome.20s<-keggGet("path:pxb03050")

KEGG.Proteasome.20s.Description<-KEGG.Proteasome.20s[[1]]$DESCRIPTION
KEGG.Proteasome.20s.Gene<-KEGG.Proteasome.20s[[1]]$GENE
KEGG.Proteasome.20s.Articles<-KEGG.Proteasome.20s[[1]]$REFERENCE

number.even.1<-seq(0,length(KEGG.Proteasome.20s[[1]]$GENE),2)
KEGG.Proteasome.20s.Gene<-KEGG.Proteasome.20s[[1]]$GENE[number.even.1]

KEGG.Proteasome.20s.Articles.References<-NULL
for(i in 1:length(KEGG.Proteasome.20s.Articles))
{
  KEGG.Proteasome.20s.Articles.References[i]<-stri_join(KEGG.Proteasome.20s.Articles[[i]]$REFERENCE,";",
                                                        KEGG.Proteasome.20s.Articles[[i]]$AUTHORS,";",
                                                        KEGG.Proteasome.20s.Articles[[i]]$TITLE,";",
                                                        KEGG.Proteasome.20s.Articles[[i]]$JOURNAL)
}
KEGG.Proteasome.20s.Articles.References.df<-as.data.frame(KEGG.Proteasome.20s.Articles.References)
colnames(KEGG.Proteasome.20s.Articles.References.df)<-c("Reference")

Proteasome.20s.df<-data.frame()
Protesome.core.particle.alpha<-c("alpha.1","alpha.2","alpha.3","alpha.4","alpha.5","alpha.6","alpha.7")
Protesome.core.particle.beta<-c("beta.1","beta.2","beta.3","beta.4","beta.5","beta.6","beta.7")
#Reverse the sequence
Protesome.core.particle.alpha.reverse<-c("alpha.7","alpha.6","alpha.5","alpha.4","alpha.3","alpha.6","alpha.7")
Protesome.core.particle.beta.reverse<-c("beta.7","beta.6","beta.5","beta.4","beta.3","beta.2","beta.1")
Proteasome.20s.Middle<-c(Protesome.core.particle.alpha,
                         Protesome.core.particle.beta,
                         Protesome.core.particle.beta.reverse,
                         Protesome.core.particle.alpha.reverse)
Proteasome.20s.Assembled.df<-as.data.frame(c(Proteasome.20s.Middle))



Proteasome.Model.1<-function(X)
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
  
  Sequence.Chain.Model.1<-function(X){output<-list();return(output)}
  
  #-----------------------------------------------------#
  #-----------------------Residue Selection-------------#
  #-----------------------------------------------------#
  
  Operator.Selector.Residue<-function(X){ output<-list(); return(output)}
 
  #---------------------------------------#
  #-------Network Dependency Modeling-----#
  #---------------------------------------#	

  Graph.Dependent.Model<-function(X){ output<-list(); return(output)}
  
  #---------------------------------------#
  #-------Deformation energies------------# 
  #---------------------------------------#	
  Energy.Defomration.Model<-function(X){ output<-list(); return(output)}
  
  #---------------------------------------#
  #----------Fluctuations-----------------# 
  #---------------------------------------#
  
  Energy.Fluctuations.Model<-function(X){ output<-list(); return(output)}

  #---------------------------------------#
  #-------Torsion Analysis By Chain-------#
  #---------------------------------------#
  
  Chain.Torsion.Analysis.Model<-function(X){ output<-list(); return(output)}

  #---------------------------------------#
  #---cross-correlations------------------#
  #---------------------------------------#
  
  Moment.1<-function(X){output<-list(); return(output)}
  Moment.2<-function(X){output<-list(); return(output)}
  Moment.3<-function(X){output<-list(); return(output)}
  Moment.4<-function(X){output<-list(); return(output)}

  #---------------------------------------------#
  #---------------Graph Theory Analysis---------#
  #---------------------------------------------#
  
Graph.Analysis.1<-function(X){output<-list(); return(output)}
  
	
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Proteasome.Model.1<-Proteasome.Model.1("1")
test.Proteasome.Model.1
#----------------Proteasome Model - II----------#
#--------26s Proteasome Assembly-----------------

KEGG.Proteasome.26s.RPN.13.publications<-entrez_search(db="pubmed", term="26S proteasome regulatory subunit RPN13-like", retmax=40)
Proteasome.26s.df<-data.frame()
Proteasome.26s.lid<-c("RPN.3","RPN.5","RPN.6","RPN.7","RPN.8","RPN.9","RPN.11","RPN.12","RPN.15")
Proteasome.26s.base<-c("RPN.1","RPN.2","RPN.13","RPT.1","RPT.2","RPT.6","RPT.4","RPT.5","RPT.3")
Proteasome.26s.Top<-c(Proteasome.26s.lid,Proteasome.26s.base)
Proteasome.26s.Bottom<-c(Proteasome.26s.base,Proteasome.26s.lid)
Proteasome.26s.Assembled.df<-as.data.frame(c(Proteasome.26s.Top,Proteasome.20s.Middle,Proteasome.26s.Bottom))


Proteasome.Model.2<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Proteasome.Model.2<-Proteasome.Model.2("1")
test.Proteasome.Model.2
#----------------Proteasome Model - III---------#
#--------Immunoproteasome-----

Immunoproteasome.20s.df<-data.frame()
Immunoproteasome.20s.core.Regulatory.Particle<-c("PA28.alpha","PA28.beta") # hetro hexamer or heptamer
Immunoproteasome.20s.core.Regulatory.Particle.1<-c("PA28.gamma") #, homo hexamer
Immunoproteasome.20s.core.particle.beta<-c("beta.5i","beta.2i","beta.1i")
Immunoproteasome.20s.core.Regulatory.Particle.reverse<-c("PA28.alpha","PA28.beta")
Immunoproteasome.20s.core.particle.beta.reverse<-c("beta.1i","beta.2i","beta.5i")

Immunoproteasome.20s.Assembled.df<-as.data.frame(c(Immunoproteasome.20s.core.Regulatory.Particle,
                                                   Immunoproteasome.20s.core.particle.beta,
                                                   Proteasome.20s.Middle,
                                                   Immunoproteasome.20s.core.particle.beta.reverse,
                                                   Immunoproteasome.20s.core.Regulatory.Particle.reverse))

Thymoproteasome.df<-data.frame()
Thymoproteasome.core.particle.beta<-c("beta.5t")
Thymoproteasome.Assembled.df<-as.data.frame(c(Proteasome.26s.Top,Proteasome.20s.Middle,Proteasome.26s.Bottom))

Proteasome.Model.3<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Proteasome.Model.3<-Proteasome.Model.3("1")
test.Proteasome.Model.3
#----------------Proteasome Model - IV----------#

#-----------------------------------------------Ubiquitin mediated proteolysis------------------------------------------
Ubiquitin.enzyme.activation.E1<-c("UB.E1","UB.L.E1.A","UB.L.E1.B","UB.E1.C")
Ubiquitin.enzyme.conjugation.E2<-c("UB.E2.A","UB.E2.B","UB.E2.C","UB.E2.D",
                                   "UB.E2.E","UB.E2.F","UB.E2.G.1","UB.E2.G.2",
                                   "UB.E2.H","UB.E2.I","UB.E2.J.1","UB.E2.J.2",
                                   "UB.E2.L.3","UB.E2.L.6","UB.E2.M","UB.E2.N",
                                   "UB.E2.O","UB.E2.Q","UB.E2.R","UB.E2.S",
                                   "UB.E2.U","UB.E2.W","UB.E2.Z","HIP2","APOLLON")

Ubiquitin.ligase.E3.HECT<-c("UB.3.B","UB.E.3.C")
Ubiquitin.ligase.E3.U-Box<-c("UB.E4.A","UB.E4.B","CHIP","CYC4","PRP19","UIP5")
Ubiquitin.ligase.E3.RING.Finger.Single<-c("Mdm2","BRCA1","COP1")
Ubiquitin.ligase.E3.RING.Finger.Subunit.Multi<-c("")
Ubiquitin.APC.subunits.other<-c("Apc.1","Apc.3","Apc.4","Apc.5","Apc.6","Apc.7","Apc.8","Apc.9","Apc.10","Apc.12","Apc.13")
#
#-----------------------------------------------Algorithm---------------------------------------------------------------------
#
#Step 1. Ub connect E1. Result:ATP->AMP
#Step 2. Ub enters domain:E2,E3. Result: Ub connect E2
#Step 3. Ub enters domain:target recognition. Result: UB connect Target, Ub replication
#Step 4. Ub Polyubiquitination Result: Ub replication
#Step 5. Ub enters Proteasome Result: Ub separated Ub.s
#Step 6. Ub.s recycled

Polyubiquitination.publications<-entrez_search(db="pubmed", term="Polyubiquitination", retmax=40)
Ubiquitin.mediated.proteolysis.publications<-entrez_search(db="pubmed", term="Ubiquitin mediated proteolysis", retmax=40)

Proteasome.Model.4<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Proteasome.Model.4<-Proteasome.Model.4("1")
test.Proteasome.Model.4
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
#---------------------------------------------------------------------#
#------------------------------Tables---------------------------------#
#---------------------------------------------------------------------#

Table.1<-xtable(Genetic.Information.Processing.df)
Table.2<-xtable(Proteasome.26s.Assembled.df)
Table.3<-xtable(Proteasome.20s.Assembled.df)
Table.4<-xtable(Immunoproteasome.20s.Assembled.df)
Table.5<-xtable(Thymoproteasome.Assembled.df)

#----------Table 1 Group A---------#

#----------Table 2 Group B---------#

#----------Table 3 Group C---------#

#----------Table 4 Group D---------#
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


