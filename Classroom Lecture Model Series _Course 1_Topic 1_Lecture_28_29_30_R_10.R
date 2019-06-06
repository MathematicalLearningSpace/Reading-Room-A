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

#---------------------------------------------------------------------#
#------------------------------Analysis-------------------------------#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#------------------------------Tables---------------------------------#
#---------------------------------------------------------------------#

#------------------Table 1--------------#
#------------------Table 2--------------#
#------------------Table 3--------------#
#------------------Table 4--------------#

#---------------------------------------------------------------------#
#------------------------------Figures--------------------------------#
#---------------------------------------------------------------------#

#------------------Figure 1--------------#
#------------------Figure 2--------------#
#------------------Figure 3--------------#
#------------------Figure 4--------------#

#---------------------------------------------------------------------#
#------------References Design Pattern--------------------------------#
#---------------------------------------------------------------------#


#-------------------------Application II------------------------------#

#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#

article.2.files <- list.files(patt='publications_Topic_*.*csv$')
article.2.files
articles.2.df<-data.frame()
articles.2.df<-read_csv(article.2.files[1])
#---------------------------------------------------------------------#
#------------------------------Functions------------------------------#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#------------------------------Analysis-------------------------------#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#------------------------------Tables---------------------------------#
#---------------------------------------------------------------------#

#------------------Table 1--------------#
#------------------Table 2--------------#
#------------------Table 3--------------#
#------------------Table 4--------------#


#---------------------------------------------------------------------#
#------------------------------Figures--------------------------------#
#---------------------------------------------------------------------#

#------------------Figure 1--------------#
#------------------Figure 2--------------#
#------------------Figure 3--------------#
#------------------Figure 4--------------#

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

#---Model Estimation of Drift and Diffusion Vectors---#
#---Estimate Model Coefficients-----------------------#
#---Generate Temporal Sequence Estimated Coefficients---#

#---------------------------------------------------------------------#
#------------------------------Analysis-------------------------------#
#---------------------------------------------------------------------#

#-------------------Condition Testing----------------#

#---------------------------------------------------------------------#
#------------------------------Tables---------------------------------#
#---------------------------------------------------------------------#

#------------------Table 1--------------#
#------------------Table 2--------------#
#------------------Table 3--------------#
#------------------Table 4--------------#

#---------------------------------------------------------------------#
#------------------------------Figures--------------------------------#
#---------------------------------------------------------------------#


#------------------Figure 1--------------#
#------------------Figure 2--------------#
#------------------Figure 3--------------#
#------------------Figure 4--------------#


#---------------------------------------------------------------------#
#------------References Design Pattern--------------------------------#
#---------------------------------------------------------------------#

