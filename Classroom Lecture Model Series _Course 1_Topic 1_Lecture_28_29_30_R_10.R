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

#------Scientific  Visualization-----#
library(corrplot);library(plot3D);library(scatterplot3d);library(rgl)

#-------------------------Application I-------------------------------#
library(koRpus);library(topicmodels);library(tm);library(XML);library(slam)
library(lattice);library(dplyr);library(tidyr);require(ggplot2);require(igraph)

Application.1<-function(X)
 {
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
 
 
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 
 }
test.Application.1<-Application.1("1")
test.Application.1
#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#

article.1.files <- list.files(patt='publications_Topic_*.*csv$')
article.1.files
articles.1.df<-data.frame()
articles.1.df<-read_csv(article.1.files[1])

#-------------Review Notes-------------------#

Search.terms<-c("low dose radiation induced bystander effect",
"genomic instability",
"radiation hypersensitivity", 
"hormesis",
"radioadaptive transgenerational responses", 
"intra signaling",
"intercellular signaling", 
"reactive oxygen species transient persistent signaling", 
"cytokines release bystander effect",
"epigenetic changes", 
"translesional responses", 
"DNA repair capacity",
"TK6 cells") 

X<-fetch_pubmed_data(get_pubmed_ids(Search.terms[3]))
#---------------------------------------------------------------------#
#------------------------------Functions------------------------------#
#---------------------------------------------------------------------#
search.model.LDA.Topics<-function(X.df,search.term,threshold=0.1, k =3)
{
  remove.HTM<-function(s)
  # Search function Completed in the Classroom
  collection <- Corpus(VectorSource(sapply(X.df[, "A"],remove.HTML)))
  X.dtm <- DocumentTermMatrix(collection, control = list(stemming = TRUE, stopwords = TRUE, 
                                                          minWordLength = 3,removeNumbers = TRUE, 
                                                          removePunctuation = TRUE))
  X.tdm <- TermDocumentMatrix(collection, control = list(stemming = TRUE, stopwords = TRUE, 
                                                         minWordLength = 3,removeNumbers = TRUE, 
                                                         removePunctuation = TRUE))
  v <- as.vector(X.dtm[,search.term]>1)
  X.dtm.Filter <- X.dtm[v, ]
  X.dtm.Filter.Terms<-inspect(X.dtm.Filter[, search.term])
  X.dtm.summary<-summary(col_sums(X.dtm))
  X.term_tfidf <- tapply(X.dtm$v/row_sums(X.dtm)[X.dtm$i], X.dtm$j, mean) *log2(nDocs(X.dtm)/col_sums(X.dtm > 0))
  X.dtm.term.summary<-summary(X.term_tfidf)
  X.dtm <- X.dtm[,X.term_tfidf >= threshold]
  X.dtm <- X.dtm[row_sums(X.dtm) > 0,]
  #--------------------------LDA Model Example-------------------------#
  X.TM <- list(VEM = LDA(X.dtm, k = k, control = list(seed =10)),
               VEM_fixed = LDA(X.dtm, k = k, control = list(estimate.alpha = FALSE, seed = 10)),
        Gibbs = LDA(X.dtm, k = k, 
                    method = "Gibbs",control = list(seed = 10, burnin = 1000, thin = 100, iter = 1000)),
               CTM = CTM(X.dtm, k = k, control = list(seed =10, var = list(tol = 10^-4), em = list(tol = 10^-3))))
  sapply(X.TM[1:2], slot, "alpha")
  methods <- c("VEM", "VEM_fixed", "Gibbs", "CTM")
  DF <- data.frame(posterior = unlist(lapply(X.TM, 
                                             function(x) apply(posterior(x)$topics, 1, max))),
                                             method = factor(rep(methods,each = nrow(posterior(X.TM$VEM)$topics)), methods))
  X.posterior.topics<-sapply(X.TM, function(x) mean(apply(posterior(x)$topics, 1, function(z) - sum(z * log(z)))))
  X.Topic <- topics(X.TM[["VEM"]], 1)
  X.Terms.1 <- terms(X.TM[["VEM"]], 5)
  search.output<-list()
  search.output$corpus<-corpus
  search.output$X.dtm.summary<-X.dtm.summary
  search.output$X.dtm.term.summary<-X.dtm.term.summary
  search.output$X.dtm<-X.dtm
  search.output$X.tdm<-X.tdm
  search.output$DF<-DF
  search.output$X.Topic<-X.Topic
  search.output$X.Terms.1<-X.Terms.1 
  search.output$X.posterior.topics<-X.posterior.topics
  search.output$X.dtm.Filter.Terms<- X.dtm.Filter.Terms
  return(search.output)
}
  
#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#
Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Model.1<-Model.1("1")
test.Model.1
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

#------------------Table 1--------------#
#------------------Table 2--------------#
#------------------Table 3--------------#
#------------------Table 4--------------#
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

#------------------Figure 1--------------#
#------------------Figure 2--------------#
#------------------Figure 3--------------#
#------------------Figure 4--------------#
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
#---------------------------------------------------------------------#
#------------References Design Pattern--------------------------------#
#---------------------------------------------------------------------#


#-------------------------Application II------------------------------#

require(PearsonDS);require(stats);require(graphics);require(phaseR)
require(deSolve);require(car);require(xtable);require(tuneR);require(yuima)                                   
#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#
article.2.files <- list.files(patt='publications_Topic_*.*csv$')
article.2.files
articles.2.df<-data.frame()
articles.2.df<-read_csv(article.2.files[1])
#------------------------------Noise----------------------------------#
noise.white<- noise(kind = c("white"))
noise.pink<-noise(kind = c("pink"))
noise.5.4<-noise(kind = c("power"), alpha=(5/4))
noise.7.4<-noise(kind=c("power"),alpha=(7/4))
noise.red<-noise(kind = c("red"))
noise.8.4<-noise(kind=c("power"),alpha=(8/4))
                                                                                                          
#---------------------------------------------------------------------#
#------------------------------Functions------------------------------#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#
Model.1<-function (Time.Sequence, variables.intitial, Params.1)
  {
    with(as.list(c(Params.1, variables.intitial)), 
         {
           #------------Equations Designed in the Classroom-----#
           d.1.x.dt.1<-a0
           d.1.y.dt.1<-a1
           d.1.z.dt.1<-a2
           d.1.w.dt.1<-a3
   res <- c(d.1.x.dt.1,
           d.1.y.dt.1,
           d.1.z.dt.1,
           d.1.w.dt.1)
  list(res)
         })
}                                                
#----------------Parameters-----------#   
initial.vector<-c(x=1,y=0,z=0,w=0)
Params.1<-c(a0<-1,a1<-1,a2<-1,a3<-1,a4<-1,a5<-1)                                                    
Time.Sequence<- seq(0, 48, by = 48/120) 

#---------------------------------------------------------------------#
#------------------------------Solutions------------------------------#
#---------------------------------------------------------------------#
model.system.solution.1 <- ode(y = initial.vector, times = Time.Sequence, 
                                             func = Model.1, parms = Parms.1) 
                                                          
#---------------------------------------------------------------------#
#------------------------------Analysis-------------------------------#
#---------------------------------------------------------------------#
options(digits = 3)
model.system.1.summary<-summary(model.system.solution.1)
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

#------------------Table 1--------------#
Table.1<-xtable(t(model.system.1.summary))             
#------------------Table 2--------------#
#------------------Table 3--------------#
#------------------Table 4--------------#
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

#------------------Figure 1--------------#
png(file = stringr::str_c("Figures//Example_",1,"_Figure_",0,".png"))
Figure.1<-matplot(model.system.solution.1[,2:5], type = "l", lty = 1, lwd = c(2, 1, 1,1),
                  col = c("darkred", "darkblue", "darkgreen","red"),
                  xlab = "Time [min]", ylab= "Y",
                  main = "Solution Example")
grid()
dev.off()                               
#------------------Figure 2--------------#
png(file = stringr::str_c("Figures//Example_",1,"_Figure_",1,".png"))
Figure.1<-scatterplotMatrix(model.system.solution.1[,2:5])
dev.off() 
#------------------Figure 3--------------#
#------------------Figure 4--------------#
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
#---------------------------------------------------------------------#
#------------References Design Pattern--------------------------------#
#---------------------------------------------------------------------#



#-------------------------Application III-----------------------------#
require(DiffusionRimp);require(Langevin);
require(HMM);require(markovchain)
#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#
article.3.files <- list.files(patt='publications_Topic_*.*csv$')
article.3.files
articles.3.df<-data.frame()
articles.3.df<-read_csv(article.3.files[1])
#---------------------------------------------------------------------#
#------------------------------Functions------------------------------#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#
Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Model.1<-Model.1("1")
test.Model.1
#---Model Estimation of Drift and Diffusion Vectors---#
#---Estimate Model Coefficients-----------------------#
#---Generate Temporal Sequence Estimated Coefficients---#

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
#-------------------Condition Testing----------------#

#---------------------------------------------------------------------#
#------------------------------Tables---------------------------------#
#---------------------------------------------------------------------#

#------------------Table 1--------------#
#------------------Table 2--------------#
#------------------Table 3--------------#
#------------------Table 4--------------#
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


#------------------Figure 1--------------#
#------------------Figure 2--------------#
#------------------Figure 3--------------#
#------------------Figure 4--------------#
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

#---------------------------------------------------------------------#
#------------References Design Pattern--------------------------------#
#---------------------------------------------------------------------#

