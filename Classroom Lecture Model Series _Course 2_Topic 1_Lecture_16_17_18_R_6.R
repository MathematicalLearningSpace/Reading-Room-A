#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#

#--------------------Work In Progress-------------------------------------#

#------------------------------R API-----------------------------------#
library(deSolve);library(ReacTran);library(rootSolve);
library(rcellminer);library(rcellminerData);library(bioCancer);
library(trena);library(ctSGE);library(anamiR) 
library(srnadiff);library(RmiR);library(CancerSubtypes);
library(GRENITS);library(compEpiTools);library(trigger);
library(CVE);library(qgraph);library(igraph)
library(sgnesR);library(rcdk);
library(fda);library(phaseR)
library(pracma);library(GA)
#----------------------------Scientific  Visualization----------------#
library(corrplot);library(plot3D);library(scatterplot3d);library(rgl)
#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#
W<-data.frame();X<-data.frame();Y<-data.frame();Z<-data.frame();
#---------------------------------------------------------------------#
#-----------------------Review Notes----------------------------------#
#-----------------Work In Progress------------------------------------#
Review.Notes<-function(X)
 {
 Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  setwd("Immune System Model/Abstracts")
  require(XML);require(xml2);require(seqinr)
  Immune.Models.Korpus.Files<-list.files()
  Immune.Models.Korpus.XML<-list();Immune.Models.Korpus.Titles<-list();Immune.Models.Korpus.Abstracts<-list();
  Immune.Models.Dictionary<-data.frame()
  for(i in 1:length(Immune.Models.Korpus.Files))
  {
    Immune.Models.Korpus.XML[[i]]<-read_xml(Immune.Models.Korpus.Files[i])
    Immune.Models.Korpus.Titles[[i]]<-xml_text(xml_find_all(Immune.Models.Korpus.XML[[i]], "//ArticleTitle"))
    Immune.Models.Korpus.Abstracts[[i]]<-xml_text(xml_find_all(Immune.Models.Korpus.XML[[i]], ".//AbstractText"),trim=TRUE)
  }	
#--------------------------------------------------------------------------------#
#-----------------------Literature NLP Topic Modeling----------------------------#
#--------------------------------------------------------------------------------#
  require(tm);require(topicmodels);require(slam);require(lattice)
  Immune.System.Model.Dictionary<-c("cell","cancer");
  Immune.System.Model.Dictionary.Correlation.Limit<-c(0.9,0.9)
  strsplit_space_tokenizer <- function(x) unlist(strsplit(as.character(x), "[[:space:]]+"))
  #-----------------------------Titles---------------------------------------#
  Immune.System.Models.Korpus<-Corpus(VectorSource(Immune.Models.Korpus.Titles[1:3]))
  #-----------------------------Abstracts------------------------------------#
  Immune.System.Models.Abstracts.Korpus<-Corpus(VectorSource(Immune.Models.Korpus.Abstracts[1:3]))
  Immune.System.Models.Korpus.DTM <- DocumentTermMatrix(Immune.System.Models.Korpus, control = list(tokenize = strsplit_space_tokenizer,
                                                                                                    stemming = TRUE, 
                                                                                                    stopwords = union(stopwords("en"),c("and")),
                                                                                                    removeWords=c("and"),
                                                                                                    minWordLength = 3,
                                                                                                    removeNumbers = TRUE,
                                                                                                    bounds=NULL,
                                                                                                    dictionary=NULL,
                                                                                                    weighting =function(x) weightTfIdf(x, normalize =FALSE),
                                                                                                    removePunctuation = TRUE))
  Immune.System.Models.Korpus.DTM
  Immune.System.Models.Korpus.TDM<-TermDocumentMatrix(Immune.System.Models.Korpus)
  Feature.1<-Docs(Immune.System.Models.Korpus.TDM)
  Feature.2<-nDocs(Immune.System.Models.Korpus.TDM)
  Feature.3<-nTerms(Immune.System.Models.Korpus.TDM)
  Feature.4<-Terms(Immune.System.Models.Korpus.TDM)
  word.Letter.1<-NULL;word.dictionary<-list();j<-1;words<-list()
  for(i in 1:length(Feature.4))
  {word.Letter.1<-s2c(Feature.4[i])
    for(k in 1:26){
    if(word.Letter.1[1]==letters[k] || word.Letter.1[1]==LETTERS[k] ){word.dictionary[[k]]<-Feature.4[i];j<-length(word.dictionary[[k]])+1;}
    }
  }
 Features<-list(Feature.1,Feature.2,Feature.3,Feature.4)
 Immune.System.Models.Korpus.TDM.tf <- termFreq(PlainTextDocument(Immune.Models.Korpus.Titles[1:3]))
 Immune.System.Models.Korpus.TDM.Most.Frequent<-findMostFreqTerms(Immune.System.Models.Korpus.TDM.tf)
 Immune.System.Models.Korpus.TDM.Most.Frequent
 Immune.System.Models.Korpus.TDM.Associations<-findAssocs(Immune.System.Models.Korpus.TDM,Immune.System.Model.Dictionary,
                                                           Immune.System.Model.Dictionary.Correlation.Limit)
Immune.System.Models.Korpus.TDM.Associations
Immune.System.Models.Korpus.DTM.Summary<-summary(col_sums(Immune.System.Models.Korpus.DTM))
  term_tfidf <- tapply(Immune.System.Models.Korpus.DTM$v/row_sums(Immune.System.Models.Korpus.DTM)[Immune.System.Models.Korpus.DTM$i], 
                       Immune.System.Models.Korpus.DTM$j, mean) *log2(nDocs(Immune.System.Models.Korpus.DTM)/col_sums(Immune.System.Models.Korpus.DTM > 0))
  Immune.System.Models.Korpus.DTM.Summary.1<-summary(term_tfidf)
  Immune.System.Models.Korpus.DTM <- Immune.System.Models.Korpus.DTM[,term_tfidf >= 0.1]
  Immune.System.Models.Korpus.DTM <- Immune.System.Models.Korpus.DTM[row_sums(Immune.System.Models.Korpus.DTM) > 0,]
  Immune.System.Models.Korpus.DTM.Summary.3<-summary(col_sums(Immune.System.Models.Korpus.DTM))
  k <- 30;SEED <- 2010
  Immune.System.Models.TM<- list(VEM = LDA(Immune.System.Models.Korpus.DTM, k = k, control = list(seed = SEED)),
                                 VEM_fixed = LDA(Immune.System.Models.Korpus.DTM, k = k, control = list(estimate.alpha = FALSE, seed = SEED)),
                                 Gibbs = LDA(Immune.System.Models.Korpus.DTM, k = k, method = "Gibbs",control = list(seed = SEED, burnin = 1000, thin = 100, iter = 1000)),
                                 CTM = CTM(Immune.System.Models.Korpus.DTM, k = k, control = list(seed = SEED, var = list(tol = 10^-4), 
                                                                                                  em = list(tol = 10^-3))))
  
  methods <- c("VEM", "VEM_fixed", "Gibbs", "CTM")
  Immune.System.Models.TM.DF <- data.frame(posterior = unlist(lapply(Immune.System.Models.TM, function(x) apply(posterior(x)$topics, 1, max))),
                                           method = factor(rep(methods,each = nrow(posterior(Immune.System.Models.TM$VEM)$topics)), methods))
  sapply(Immune.System.Models.TM, function(x) mean(apply(posterior(x)$topics, 1, function(z) - sum(z * log(z)))))
  Immune.System.Models.TM.Topic <- topics(Immune.System.Models.TM[["VEM"]], 1)
  Immune.System.Models.TM.Terms <- terms(Immune.System.Models.TM[["VEM"]], 10)
  Immune.System.Models.Most.Frequent<- which.max(tabulate(Immune.System.Models.TM.Topic))
  Immune.System.Models.Most.Frequent.Terms<-terms(Immune.System.Models.TM[["VEM"]], 10)[,Immune.System.Models.Most.Frequent]
  Immune.System.Models.Most.Frequent.Terms	
  output<-list()
  output$X<-X
  output$Article.keywords<-Article.keywords
  output$Immune.Models.Korpus.Titles<-Immune.Models.Korpus.Titles
  output$Immune.System.Model.Dictionary<-Immune.System.Model.Dictionary
  output$Immune.System.Models.Korpus.DTM<-Immune.System.Models.Korpus.DTM
  output$Immune.System.Models.Korpus.TDM.Associations<-Immune.System.Models.Korpus.TDM.Associations
  output$Immune.System.Models.TM.Terms<-Immune.System.Models.TM.Terms
  output$Immune.System.Models.TM.Topic<-Immune.System.Models.TM.Topic
  output$Immune.System.Models.Most.Frequent.Terms<-Immune.System.Models.Most.Frequent.Terms
  output$Immune.System<-Immune.System.df							 
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
Transformation.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Transformation.1<-Transformation.1("1")
test.Transformation.1
#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#
Signal.Transduction.Model.1<-function(X)
 {
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  return(output)
 }
test.Signal.Transduction.Model.1<-Signal.Transduction.Model.1("1")
Signal.Transduction.Model.1
#------------Longevity-----------------#
#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
Longevity.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
Longevity.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, 
						 func = Longevity.system.equation.model.1, parms = Params.Equation) 

#------------Immune System-------------#

Immune.System.Model<-function(X)
{
require(WGCNA)
data(ImmunePathwayLists)
	
 Immune.System.df<-as.data.frame(rbind(c("04640","Hematopoietic cell lineage"),
                       c("04610","Complement and coagulation cascades"),
                       c("04611","Platelet activation"),
                       c("04620","Toll-like receptor signaling pathway"),
                       c("04624","Toll and Imd signaling pathway"),
                       c("04621","NOD-like receptor signaling pathway"),
                       c("04622","RIG-I-like receptor signaling pathway"),
                       c("04623","Cytosolic DNA-sensing pathway"),
                       c("04625","C-type lectin receptor signaling pathway"),
                       c("04612","Antigen processing and presentation"),
                       c("04660","T cell receptor signaling pathway"),
                       c("04658","Th1 and Th2 cell differentiation"),
                       c("04659","Th17 cell differentiation"),
                       c("04657","IL-17 signaling pathway"),
                       c("04662","B cell receptor signaling pathway"),
                       c("04664","Fc epsilon RI signaling pathway"),
                       c("04666","Fc gamma R-mediated phagocytosis"),
                       c("04670","Leukocyte transendothelial migration"),
                       c("04672","Intestinal immune network for IgA production"),
                       c("04062","Chemokine signaling pathway")))

 setwd("Immune System Model")
  Immune.System.1.Name=""; Immune.System.1<-xml2::read_xml("hsa04640.xml")
  Immune.System.2.Name=""; Immune.System.2<-xml2::read_xml("hsa04610.xml")
  Immune.System.3.Name=""; Immune.System.3<-xml2::read_xml("hsa04620.xml")
  #Immune.System.4.Name=""; Immune.System.4<-xml2::read_xml("hsa04624.xml")
  Immune.System.5.Name=""; Immune.System.5<-xml2::read_xml("hsa04621.xml")
  Immune.System.6.Name=""; Immune.System.6<-xml2::read_xml("hsa04622.xml")
  Immune.System.7.Name=""; Immune.System.7<-xml2::read_xml("hsa04623.xml")
  Immune.System.8.Name=""; Immune.System.8<-xml2::read_xml("hsa04625.xml")
  Immune.System.9.Name=""; Immune.System.9<-xml2::read_xml("hsa04612.xml")
  Immune.System.10.Name=""; Immune.System.10<-xml2::read_xml("hsa04660.xml")
  Immune.System.11.Name=""; Immune.System.11<-xml2::read_xml("hsa04658.xml")
  Immune.System.12.Name=""; Immune.System.12<-xml2::read_xml("hsa04659.xml")
  Immune.System.13.Name=""; Immune.System.13<-xml2::read_xml("hsa04657.xml")
  Immune.System.14.Name=""; Immune.System.14<-xml2::read_xml("hsa04662.xml")
  Immune.System.15.Name=""; Immune.System.15<-xml2::read_xml("hsa04664.xml")
  Immune.System.16.Name=""; Immune.System.16<-xml2::read_xml("hsa04666.xml")
  Immune.System.17.Name=""; Immune.System.17<-xml2::read_xml("hsa04670.xml")
  Immune.System.18.Name=""; Immune.System.18<-xml2::read_xml("hsa04672.xml")
  Immune.System.19.Name=""; Immune.System.19<-xml2::read_xml("hsa04062.xml")
  Table.1.df<-data.frame(); Table.2.df<-data.frame(); Table.3.df<-data.frame();
  
  output<-list()
  output$X<-X
  output$Table.1<-Table.1.df
  output$Table.2<-Table.2.df
  output$Table.3<-Table.3.df
  output$Table.4<-Table.4.df
  output$Table.5<-Table.5.df
  output$Table.6<-Table.6.df
  return(output)
	
}	


#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
Immune.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
Immune.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, 
					      func = Immune.system.equation.model.1, parms = Params.Equation) 

#------------wNT-----------------------#
#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
Wnt.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
Wnt.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, 
					   func = Wnt.system.equation.model.1, parms = Params.Equation) 

#------------Notch---------------------#
#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
Notch.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
Notch.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, 
					     func = Notch.system.equation.model.1, parms = Params.Equation) 

#------------HedgeHog------------------#
#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
Hedgehog.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
Hedgehog.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, 
						func = Hedgehog.system.equation.model.1, parms = Params.Equation) 

#------------TGF-Beta------------------#
#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
TGF.Beta.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
TGF.Beta.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, 
						func = TGF.Beta.system.equation.model.1, parms = Params.Equation) 

#------------VEGF----------------------#
#------------Parameter Models----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)
#------------Equation Systems----------#
VEGF.system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}
#-------Solutions----------------------#
VEGF.system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, 
					    func = VEGF.system.equation.model.1, parms = Params.Equation) 

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

