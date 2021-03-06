#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#
#--------------------Work In Progress----August 2019----------------------#

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

#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#
W<-data.frame();X<-data.frame();Y<-data.frame();Z<-data.frame();
Table.1.df<-as.data.frame(Table.1);Table.2.df<-as.data.frame(Table.2);Table.3.df<-as.data.frame(Table.3);
Table.4.df<-as.data.frame(Table.4)Table.5.df<-as.data.frame(Table.5);Table.6.df<-as.data.frame(Table.6);
Table.7.df<-as.data.frame(Table.7);Table.8.df<-as.data.frame(Table.8)Table.9.df<-as.data.frame(Table.9);
Table.10.df<-as.data.frame(Table.10)

Cell.Cycle.df<- as.data.frame(read_delim("Participating Molecules [R-HSA-1640170].tsv", "\t", escape_double = FALSE, trim_ws = TRUE))
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

#------------------------Feature Matrix Design-----------------------#
Feature.1<-unlist(lapply(Z.2.df,function(x){Peptides::instaIndex(x)}))
Feature.2<-unlist(lapply(Z.2.df,function(x){Peptides::lengthpep(x)}))
Feature.3<-unlist(lapply(Z.2.df,function(x){Peptides::boman(x)}))
Feature.4<-unlist(lapply(Z.2.df,function(x){Peptides::charge(x,pH=5)}))
Feature.5<-unlist(lapply(Z.2.df,function(x){Peptides::charge(x,pH=7)}))
Feature.6<-unlist(lapply(Z.2.df,function(x){Peptides::charge(x,pH=9)}))
Feature.7<-unlist(lapply(Z.2.df,function(x){Peptides::aIndex(x)}))
Feature.8<-unlist(lapply(Z.2.df,function(x){Peptides::pI(x)}))
Feature.9<-unlist(lapply(Z.2.df,function(x){Peptides::mw(x)}))
Feature.10<-unlist(lapply(Z.2.df,function(x){Peptides::hmoment(x,angle=100,window=11)}))
Feature.11<-unlist(lapply(Z.2.df,function(x){Peptides::hmoment(x,angle=160,window=11)}))
#----------------------Feature Design Matrix
Response.Vector<-ifelse(Feature.1>= 40, 1, 0)
Matrix.Feature<-cbind(Matrix.Label,
                      Feature.1,
                      Feature.2,
                      Feature.3,
                      Feature.4,
                      Feature.5,
                      Feature.6,
                      Feature.7,
                      Feature.8,
                      Feature.9,
                      Feature.10,
                      Feature.11)
Matrix.Feature
#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#

#----------------Machine Learning Models I-----------#
Machine.Learning.Model.1<-function(X,Visualization=FALSE)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  require(ggRandomForests)
  Model.1 <- rfsrc(X$Dependent.1 ~ ., data = X)
  Model.ROC.1 <- gg_roc(Model.1, which.outcome=1)
  Model.ROC.2 <- gg_roc(Model.1, which.outcome=2)
  Model.ROC.3 <- gg_roc(Model.1, which.outcome=3)
  Model.ROC.4 <- gg_roc(Model.1, which.outcome=4)

Table.1.df<-cbind(Model.ROC.1$sens,Model.ROC.1$spec,Model.ROC.1$pct)
Table.2.df<-cbind(Model.ROC.2$sens,Model.ROC.2$spec,Model.ROC.2$pct)
Table.3.df<-cbind(Model.ROC.3$sens,Model.ROC.3$spec,Model.ROC.3$pct)
Table.4.df<-cbind(Model.ROC.4$sens,Model.ROC.4$spec,Model.ROC.4$pct)

 if(Visualization)
  {
   png(file = stringr::str_c('Figures//Example_',"MLM_1",'_Figure_',1,'.png'))
   op <- par(mfrow = c(2,2),mar=c(1,1,1,1))
   plot(Model.ROC.1);plot(Model.ROC.2);plot(Model.ROC.3);plot(Model.ROC.4);
   par(op)
   dev.off()
 } 
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  output$Table.4<-Table.4.df
  return(output)
 }
test.Machine.Learning.Model.1<-Machine.Learning.Model.1("1",TRUE)
test.Machine.Learning.Model.1
#----------------Machine Learning Models II----------#
Machine.Learning.Model.2<-function(X)
 {
require(rcdk);require(BioMedR);require(ChemmineR);require(Matrix);require(fingerprint)
require(igraph);require(readr);require(leaps);require(caret);require(GA) 
Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
 #--------------------------Molecular Descriptors-------------------#
 Molecular.Categories <- get.desc.categories()
 dcn.1<-get.desc.names(Molecular.Categories[1])
 dcn.2<-get.desc.names(Molecular.Categories[2])
 dcn.3<-get.desc.names(Molecular.Categories[3])
 dcn.4<-get.desc.names(Molecular.Categories[4])
 dcn.5<-get.desc.names(Molecular.Categories[5])
 Molecular.Categories.Names<- unique(unlist(sapply(get.desc.categories(), get.desc.names)))
 
 Z.mols.ids<-c("")
 Z<-BMgetDrugSmiKEGG(Z.mols.ids[1])
 Z.mols <- parse.smiles(Z)
 Z.mols.1 <- eval.desc(Z.mols, dcn.1)
 #-----------------Topological Metrics------------------------------#
 Z.Complexity<-extrDrugFragmentComplexity(Z.mols, silent = TRUE)
 
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Machine.Learning.Model.2<-Machine.Learning.Model.2("1")
test.Machine.Learning.Model.2
#----------------Machine Learning Models III---------#
Machine.Learning.Model.3<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  require("protr");require("randomForest");require("pROC")
 # calculate APseAAC descriptors
 #x1 <- t(sapply(Y.df$Y.FASTA[1:100], extractAPAAC))
 #x2 <- t(sapply(Y.df$Y.FASTA[100:200], extractAPAAC))
X.Feature.1<-X[1:50,]
X.Feature.2<-X[51:100,]
x <- rbind(X.Feature.1, X.Feature.2)
labels.1 <- as.factor(c(rep(0, nrow(X.Feature.1)), rep(1, nrow(X.Feature.2))))
train.idx <- c(
  sample(1:nrow(X.Feature.1), round(nrow(X.Feature.1) * 0.75)),
  sample(nrow(X.Feature.1) + 1:nrow(X.Feature.2), round(nrow(X.Feature.2) * 0.75))
)
test.idx <- setdiff(1:nrow(x), train.idx)
x.train <- x[train.idx, ]
x.test <- x[test.idx, ]
y.train <- labels.1[train.idx]
y.test <- labels.1[test.idx]
model.rf.fit <- randomForest(x.train, y.train, cv.fold = 5)
model.rf.pred <- predict(model.rf.fit, newdata = x.test, type = "prob")[, 1]
plot.roc(y.te, rf.pred, grid = TRUE, print.auc = TRUE)
 output<-list()
 output$X<-X
 output$Model.Predictions<-model.rf.pred
 output$Table.1<-Table.1.df
 output$Table.2<-Table.2.df
 output$Table.3<-Table.3.df
  return(output)
 }
test.Machine.Learning.Model.3<-Machine.Learning.Model.3("1")
test.Machine.Learning.Model.3
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

#------------Table 1---------------------#
#------------Table 2---------------------#
#------------Table 3---------------------#
Table.1.TeX<-xtable::xtable(Table.1.df)
#----------Table 2------#
Table.2.TeX<-xtable::xtable(Table.2.df)
#----------Table 3------#
Table.3.TeX<-xtable::xtable(Table.3.df)
#----------Table 4------#
Table.4.TeX<-xtable::xtable(Table.4.df)
#---------------------------------------------------------------------#
#------------------------------Figures--------------------------------#
#---------------------------Work in Progress--------------------------#

#------------Figure 1---------------------#
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
#------------Figure 2---------------------#
for(i in 1:3)
{
png(file = stringr::str_c('Figures//Example_',i,'_Figure_',3,'.png'))
heatmap.2(exprs(esetSel), col=topo.colors(75), scale="none", ColSideColors=patientcolors,
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
dev.off()
}
#------------Figure 3---------------------#
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
