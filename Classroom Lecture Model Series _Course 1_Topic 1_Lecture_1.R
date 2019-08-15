#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#
#----------------------Work In Progress--------------------#
library(deSolve);library(ReacTran);library(rootSolve);
library(fda);library(phaseR);library(pracma);
library(xtable);library(GA);library(igraph);library(pracma};library(readr);require(BioMedR)
#----------------------------------R Source Files------------------------------------#
#----------------------------------Data-------------------------------------------#
#----------------------------------Cancer Example: Gastric Cancer----------------#
						    
A.Drug.Resistance<-c(CDX2, MUC2, REG4, CDH17, MDR1, SHH)
B.Genomic.Instability<-c(p53, p21, BAX, p48, GADD45, BAK, POLK)
C.Tumor.Progression<-c(Retinoic.Acid, RAR.Beta, RXR)
D.Intestinal.Metaplasia<-c(DV1, GSK.3Beta, Beta.Catenin, Axin,APC, CK1.alpha,GBP)
E.Dysplasia.Path.1<-c(EFG, ERBB2, SHC, GRB2, SOS, RAS, RAF, MEK, ERK.1)
F.Dysplasia.Path.2<-c(PI3K, PIP3, AKT, mTOR, p53, S6K, BCL2)
G.Dysplasia.Path.3<-c(TGF.Beta,TGF.BetaRI,TGF.BetaRII, SMAD.2, SMAD.4, p15, p21)
H.Normal.Gastic.Muscosa.1<-c(HGF, c.MET, GRB2, SOS, RAS, RAF, MEK, ERK.1)
I.Normal.Gastic.Muscosa.Survival.Path.1<-c(FGF, FGFR2, GAB1, PI3K, PIP3, AKT, mTOR, GSK.3Beta)

Cancer.Models.<-c("Nutrition Model","Cancer Model","Stomach Model","Intestine Model","Digestion Model","Protein Model","Carbohydrates Model","Fats Model",
"Metabolism Model","Signal Transduction Model","Ribosome Model","Chaperonin Model","Proteasome Model")
						    
#-------------------Parameter Models						    
params.1<-c(a11=0.1,a12=0.1,a13=0.1,a14=0.1,a15=0.1,a16=0.1,
            a21=0.1,a22=0.1,a23=0.1,a24=0.1,a25=0.1,a26=0.1,
            a31=0.1,a32=0.1,a33=0.1,a34=0.1,a35=0.1,a36=0.1,
            a41=0.1,a42=0.1,a43=0.1,a44=0.1,a45=0.1,a46=0.1,
            a51=0.1,a52=0.1,a53=0.1,a54=0.1,a55=0.1,a56=0.1,
            a61=0.1,a62=0.1,a63=0.1,a64=0.1,a65=0.1,a66=0.1)

params.epsilon.1<-c(epsilon1,epsilon2,epsilon3)
					    					    
#---------------------------Review Notes--------------------------#
						    
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

#-------------------Gene/Protein/Enzyme Study--#					   
Gastric.Cancer.Genes<-c("CDX2", "MUC2", "REG4", "CDH17", "MDR1", "SHH","p53", "p21", 
                        "BAX", "p48", "GADD45", "BAK", "POLK","Retinoic.Acid", "RAR.Beta", "RXR",
                        "DV1", "GSK.3Beta", "Beta.Catenin", "Axin", "APC", "CK1.alpha","GBP","EFG", 
                        "ERBB2", "SHC", "GRB2", "SOS", "RAS", "RAF", "MEK", "ERK.1","PI3K", "PIP3", 
                        "AKT", "mTOR", "p53", "S6K", "BCL2","TGF.Beta","TGF.BetaRI","TGF.BetaRII", "SMAD.2", "SMAD.4", "p15", "p21",
                        "HGF", "c.MET", "GRB2", "SOS", "RAS", "RAF", "MEK", "ERK.1",
                        "FGF", "FGFR2", "GAB1", "PI3K", "PIP3", "AKT", "mTOR", "GSK.3Beta")
					    
#------------------------------------------------------------------------------------#
#----------------------------------Transformations-----------------------------------#
#------------------------------------------------------------------------------------#
Transformation.1<-function(X,T.Category)
{
if(T.Category==LETTERS[1]){Y<-X}
if(T.Category==LETTERS[2]){Y<-1}
if(T.Category==LETTERS[3]){Y<-2}
if(T.Category==LETTERS[4]){Y<-3}
	
output<-list()
output$X<-X
output$Y<-Y
output$T.Category<-T.Category
return(output)
}
test.Transformation.1<-Transformation.1(c(1,1),"A")
test.Transformation.1			    
						    
#------------------------------------------------------------------------------------#						    
#----------------------------------User Defined Modules and Functions----------------
#------------------------------------------------------------------------------------#
Module.1<-function(X)
{
	F.1<-function(X){output<-list();output$X<-X; return(output)}
	F.2<-function(X){output<-list();output$X<-X; return(output)}
	F.3<-function(X){output<-list();output$X<-X; return(output)}
	F.4<-function(X){output<-list();output$X<-X; return(output)}
	F.5<-function(X){output<-list();output$X<-X; return(output)}
	
	Z1<-F.1(X)$X
	Z2<-F.2(Z1)$X
	Z3<-F.3(Z2)$X
	Z4<-F.4(Z3)$X
	Z5<-F.5(Z4)$X
	
output<-list()
output$X<-X
output$Z5<-Z5
return(output)
}
test.Module.1<-Module.1(c(1,1))
test.Module.1
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
Analysis.Graph.1<-function(X)
{
output<-list()
return(output)
}	
#----------------------------------Optimization--------------------------------------#
Optimization.Model.1<-function(X)
{
output<-list()
return(output)
}	
#----------------------------------Natural Language Description----------------------#
Language.Natural.1<-function(X)
{

output<-list()
return(output)
}
#----------------------------------Tables--------------------------------------------#
#--------------Table 1----------#
						    
#--------------Table 2----------#
						    
#--------------Table 3----------#	

#----------------------------------Figures-------------------------------------------#
#-----------------------Figure 1-------------------#
png(file = stringr::str_c('Figures/1/Example_',1,'_Figure_',1,'.png'))
#plot(yout[,-1], type = "l")
dev.off()
png(file = stringr::str_c('Figures/1/Example_',1,'_Figure_',2,'.png'))
#plot(yout[,-1], type = "l")
dev.off()
png(file = stringr::str_c('Figures/1/Example_',1,'_Figure_',3,'.png'))
#plot(yout[,-1], type = "l")
dev.off()
png(file = stringr::str_c('Figures/1/Example_',1,'_Figure_',4,'.png'))
#plot(yout[,-1], type = "l")
 dev.off()
#-----------------------Figure 2-------------------#
      
#----------------------------------Discussion----------------------------------------#
