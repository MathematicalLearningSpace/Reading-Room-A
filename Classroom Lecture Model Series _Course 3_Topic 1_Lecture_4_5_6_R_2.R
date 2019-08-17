#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#
#--------------------Work In Progress--------August---2019----------------#
#------------------------------R API----------------------------------#
library(gdata);library(bio3d);library(igraph);library(sna);library(ips);
library(phangorn);library(proteomics);
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
library(triplex);library(BSgenome.Hsapiens.UCSC.hg19);library(heatmaps);library(BSgenome);library(GenomeGraphs);library(biomaRt)
#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#
Table.1.df<-as.data.frame(Table.1);Table.2.df<-as.data.frame(Table.2);Table.3.df<-as.data.frame(Table.3);
Table.4.df<-as.data.frame(Table.4)Table.5.df<-as.data.frame(Table.5);Table.6.df<-as.data.frame(Table.6);
Table.7.df<-as.data.frame(Table.7);Table.8.df<-as.data.frame(Table.8)Table.9.df<-as.data.frame(Table.9);
Table.10.df<-as.data.frame(Table.10)

pearson.N<-10^2
#-------------------------------------Moment Parameters for the Distributions-------------------------------------------------
p0pars <- list(mean=1, sd=1)
pIpars <- list(a=1, b=1, location=1, scale=1) 
pIIpars <- list(a=1, location=1, scale=1) 
pIIIpars <- list(shape=1, location=1, scale=1)
pIVpars <- list(m=3, nu=1, location=1, scale=1)
pVpars <- list(shape=1, location=1, scale=1) 
pVIpars <- list(a=1, b=1, location=1, scale=1) 
pVIIpars <- list(df=10, location=1, scale=1)
#-------------------------------------Generate Random variables from the Distributions----------------------------------------
error.pearson.0<-rpearson0(pearson.N,params=p0pars)
error.pearson.1<-rpearsonI(pearson.N,params=pIpars)
error.pearson.2<-rpearsonII(pearson.N,params=pIIpars)
error.pearson.3<-rpearsonIII(pearson.N,params=pIIIpars)
error.pearson.4<-rpearsonIV(pearson.N,params=pIVpars)
error.pearson.5<-rpearsonV(pearson.N,params=pVpars)
error.pearson.6<-rpearsonVI(pearson.N,params=pVIpars)
error.pearson.7<-rpearsonVII(pearson.N,params=pVIIpars)

strand.types<-c("Parallel-First-0","Parallel-First-1",
                "Parallel-Second-2","Parallel-Second-3",
                "AntiParallel-Second-4", "AntiParallel-Second-5",
                "AntiParallel-First-6", "AntiParallel-First-7")
genome <- BSgenome.Hsapiens.UCSC.hg19
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
#---------Translational and Transcriptional Error Models I------------#
Error.Translation.Transcriptional.Model.0<-function(X)
 {
require(sysBio)
Model.1 <- newModel("Translation Example") 
addMAreaction(Model.1, "DNA -> DNA + mRNA", r1="v1", name="Transcription")
addMAreaction(Model.1, "mRNA -> mRNA + protein", r1="v2", name="Translation")
addMAreaction(Model.1, "DNA + protein -> DNA_protein", r1="v3", name="Binding")
addMAreactRate(Model.1, "v1", "assigned", "a11*DNA")
addMAreactRate(Model.1, "v2", "assigned", "a21*mRNA")
addMAreactRate(Model.1, "v3", "assigned", "a31*DNA*protein")
addMAreactRate(Model.1, "v4", "assigned", "a41*DNA_protein")
addMAreactRate(Model.1, "v5", "assigned", "a51*mRNA")
addSpecies(Model.1, "DNA", 50)
addSpecies(Model.1, "mRNA", 0)
addSpecies(Model.1, "protein", 0)
addSpecies(Model.1, "DNA_protein", 0)
addParameters(Model.1, "a11", 0.2)
addParameters(Model.1, "a21", 20)
addParameters(Model.1, "a32", 0.2)
addParameters(Model.1, "a41", 1)
addParameters(Model.1, "a51", 1.5)
addParameters(Model.1, "a61", 1)
#------------Biology Rules
#------------Botany Rules
#------------Chemistry Rules
addRule(Model.1, "rule 1","ODEs","DNA_protein=0.5*DNA")
makeModel(Model.1)
#----------------------Simulate the Model in the Classroom-------------------------------------------
simResults.1 <-simulateModel(Model.1, times = seq(0, 100, by = 0.1)) 
simResults.1.stoch <- solveStoch(Model.1,10,method = "D", simName = "", tau = 0.3,f = 10, epsilon = 0.03, nc = 10) 
distribution.moments.1.df<-data.frame()
distribution.moments.1.df<-rbind(c("ES_1","Protein",empMoments(simResults.1$protein)),
                                 c("ES_1 Stochastic","Protein",empMoments(simResults.1.stoch$protein)))
colnames(distribution.moments.1.df)<-c("Equation","Variable","Moment 1", "Moment 2","Moment 3","Moment 4")
Table.1<-xtable(distribution.moments.1.df)
if(visualization)
  {
Figure.1<-plotResults(simResults.1, title="Model 1") 
Figure.2<-plotResults(simResults.1.stoch, title="Model 1 Stochastic")
  }
 output<-list()
  output$X<-X
  output$Table.1<-Table.1
  return(output)
 }
test.Error.Translation.Transcriptional.Model.0<-Error.Translation.Transcriptional.Model.0("1")
test.Error.Translation.Transcriptional.Model.0


Error.Translation.Transcriptional.Model.1<-function(X)
 {
 require(seqinr);require(ape);require(protr);require(Peptides)
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
#---------------Generate Theoretical Data Set Based on Alphabet--------------------------------------------------------
nucleotides <- c("A", "C", "G", "T")
probabilities.1 <-c(0.2, 0.3, 0.3, 0.2) 
probabilities.2<- c(0.2, 0.3, 0.3, 0.2)
probabilities.3<- c(0.2, 0.3, 0.3, 0.2)
probabilities.4<- c(0.2, 0.3, 0.3, 0.2)
seqlength <- 100
DNA.Sequence.1<-sample(nucleotides, seqlength, rep=TRUE, prob=probabilities.1)
DNA.Sequence.2<-sample(nucleotides, seqlength, rep=TRUE, prob=probabilities.2)
DNA.Sequence.3<-sample(nucleotides, seqlength, rep=TRUE, prob=probabilities.3)
DNA.Sequence.4<-sample(nucleotides, seqlength, rep=TRUE, prob=probabilities.4)
#----------------------Translate Data to Amino Acids Data Type----------------------------------------------------
AA.DNA.Sequence.1<-c2s(translate(as.character(DNA.Sequence.1),frame=0,sens='F',numcode=1))
AA.DNA.Sequence.2<-c2s(translate(as.character(DNA.Sequence.2),frame=0,sens='F',numcode=1))
AA.DNA.Sequence.3<-c2s(translate(as.character(DNA.Sequence.3),frame=0,sens='F',numcode=1))
AA.DNA.Sequence.4<-c2s(translate(as.character(DNA.Sequence.4),frame=0,sens='F',numcode=1))

  system.equation.model.test<-function(times, variables.intitial.test, parameters.test)
{
  with(as.list(c(parameters.test, variables.intitial.test)), 
       {
#--------------Group 1---------------------
d.1.X1.d.t.1<--a11*X1 + a12*X2
#--------------Group 2---------------------
d.1.X2.d.t.1<-a72*X2 - a81*X1
res <- c(d.1.X1.d.t.1,d.1.X2.d.t.1)
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
test.Error.Translation.Transcriptional.Model.1<-Error.Translation.Transcriptional.Model.1("1")
test.Error.Translation.Transcriptional.Model.1
#---------Translational and Transcriptional Error Models II-----------#
Error.Translation.Transcriptional.Model.2<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
 system.equation.model.test<-function(times, variables.intitial.test, parameters.test)
{
  with(as.list(c(parameters.test, variables.intitial.test)), 
       {

#--------------Group 1---------------------
d.1.X1.d.t.1<--a11*X1 + a12*X2
#--------------Group 2---------------------
d.1.X2.d.t.1<-a72*X2 - a81*X1
res <- c(d.1.X1.d.t.1,d.1.X2.d.t.1)
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
test.Error.Translation.Transcriptional.Model.2<-Error.Translation.Transcriptional.Model.2("1")
test.Error.Translation.Transcriptional.Model.2
#---------Translational and Transcriptional Error Models III----------#
Error.Translation.Transcriptional.Model.3<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
 system.equation.model.test<-function(times, variables.intitial.test, parameters.test)
{
  with(as.list(c(parameters.test, variables.intitial.test)), 
       {

#--------------Group 1---------------------
d.1.X1.d.t.1<--a11*X1 + a12*X2
#--------------Group 2---------------------
d.1.X2.d.t.1<-a72*X2 - a81*X1
res <- c(d.1.X1.d.t.1,d.1.X2.d.t.1)
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
test.Error.Translation.Transcriptional.Model.3<-Error.Translation.Transcriptional.Model.3("1")
test.Error.Translation.Transcriptional.Model.3
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

#----------Table 1------#
Table.1.TeX<-xtable::xtable(Table.1.df)
#----------Table 2------#
Table.2.TeX<-xtable::xtable(Table.2.df)
#----------Table 3------#
Table.3.TeX<-xtable::xtable(Table.3.df)
#----------Table 4------#
Table.4.TeX<-xtable::xtable(Table.4.df)
#----------Table 5------#
Table.1.TeX<-xtable::xtable(Table.5.df)
#----------Table 6------#
Table.2.TeX<-xtable::xtable(Table.6.df)
#----------Table 7------#
Table.3.TeX<-xtable::xtable(Table.7.df)
#----------Table 8------#
Table.4.TeX<-xtable::xtable(Table.8.df)

#---------------------------------------------------------------------#
#------------------------------Figures--------------------------------#
#---------------------------------------------------------------------#

#------------Figure 1---------------------#
#------------Figure 2---------------------#
#------------Figure 3---------------------#

png(file = stringr::str_c('Figures//Example_',1,'_Figure_',3,'.png'))
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
