#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#
#----------------------Work In Progress--------------------#

library(deSolve);library(ReacTran);library(rootSolve);
library(fda);library(phaseR);library(pracma);
library(xtable);library(GA);library(igraph);library(pracma};library(readr);require(BioMedR)
#----------------------------------R Source Files------------------------------------#
#----------------------------------Data----------------------------------------------#
#---Example from PubChem for the Classroom-------------------------------------------#
Gene_Targets_Patents.df <- as.data.frame(read_csv("Gene_Targets_Patents.csv"))
Protein_Targets_Patents.df <- as.data.frame(read_csv("Protein_Targets_Patents.csv"))
Enzyme_Targets_Patents.df <- as.data.frame(read_csv("Enzyme_Targets_Patents.csv"))

Review.Article.Module.1<-function(X,search.terms,Z)
{
ids<-unique(grep(search.terms,X$patentabstract))
Article.Content<-""
Article.Bibliography<-""
for(i in 1:length(ids))
{
  Article.Content<-stringr::str_c(Article.Content,
                                    X$patentabstract[ids[i]],
                                    "\\cite{",
                                     X$pid[ids[i]],"}"," \\\n "
  )
  Article.Bibliography<-stringr::str_c(Article.Bibliography,
                                   X$pid[ids[i]], " & ",
                                   X$patenttitle[ids[i]]," \\ "
                                   )
}
Article.Review<-stringr::str_c(Article.Content,Article.Bibliography)
title<-stringr::str_c("Working_Review_Article_",search.terms,"_",Z,"_Target_Patents",".txt")
base::write(Article.Review,file=title,append=FALSE)
}
test.Review.Article.Module.1<-Review.Article.Module.1(Gene_Targets_Patents.df,"cancer","Gene")
test.Review.Article.Module.2<-Review.Article.Module.1(Protein_Targets_Patents.df,"cancer","Protein")
test.Review.Article.Module.3<-Review.Article.Module.1(Enzyme_Targets_Patents.df,"cancer","Enzyme")
                                                    
#-------------------Parameter Models-----------#
sequence.time <- seq(0, 10^2, by = 0.1)
Params.Initial<-c(X<-1)
Params.Equation<-c(a11<-1)

#------------------------------------------------------------------------------------#
#----------------------------------Transformations-----------------------------------#
#------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------#						    
#----------------------------------User Defined Modules and Functions----------------
#------------------------------------------------------------------------------------#

#----------------------------------Network Designs-----------------------------------#
#------------------------------------------------------------------------------------#
					    
#----------------------------------Equation Systems----------------------------------#

system.equation.model.1<-function(sequence.time,Params.Initial, Params.Equation)
{
  with(as.list(c(Params.Equation, Params.Initial)), 
       {
	d.1.X.dt.1 = a11*X
res <- c(d.1.X.dt.1)
 list(res)
       })
}

#-------Solutions----------------------------#
system.equation.model.1.solution<- ode(y = Params.Initial, times = sequence.time, func = system.equation.model.1, parms = Params.Equation)                                                     
                                             
#----------------------------------Network Analysis----------------------------------#

#----------------------------------Optimization--------------------------------------#

#----------------------------------Natural Language Description----------------------#

#----------------------------------Tables--------------------------------------------#
#--------------Table 1----------#
#--------------Table 2----------#
#--------------Table 3----------#	

#----------------------------------Figures-------------------------------------------#
#-----------------------Figure 1-------------------#
#-----------------------Figure 2-------------------#
      
#----------------------------------Discussion----------------------------------------#
