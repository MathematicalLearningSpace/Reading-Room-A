
library(WikipediR)
library(textclean)
library(xtable)
library(readr)
#-------------Lecture Management System-----------------#
#-------------Data--------------------------------------#

Folders<-c("Album_Designs",
           "Journal_Designs",
           "Lecture_Designs",
           "Music_Composition_Designs",
           "Reading_List_Designs")



Classes<-c("Fractional Stochastic Differential Equation Theory",
           "Biological Signal Transduction Network Theory",
           "Combinatorial Optimization and Combinatorial Topology",
           "Music Theory and Composition and Mathematics")
#---------------------Course 1-----------------------

Course.1.Title<-c("Differential Equations in Mathematical Biology, Botany and Oncology")
Course.1.Outline<-c("section{Group 1}",
                    "Ordinary Differential Equations ODE",
                    "Delayed Differential Equations DDE",
                    "Stochastic Differential Equations SDE",
                    "Fractional Differential Equations FDE",
                    "DDE SDE",
                    "FDE SDE",
                    "SDE Bridges and Diffusion",
                    "Stochastic Noise Models",
                    "Markov Models",
                    "section{Group 2}",
                    "Newton and Euler",
                    "Runge Kutta ODE23 ODE45",
                    "Milstein",
                    "Stochastic Search and Integration",
                    "Baum-Welch and Viterbi",
                    "Combinatorial Optimization",
                    "Heuristic-Metaheuristic Optimization",
                    "Eigen Values, Modes, Harmonics",
                    "Distance, Connectivity and Toplogical Metrics for Adjacency Lists and Matrices",
                    "section{Group 3}",
                    "Equilibrium Analysis I",
                    "Equilibrium Analysis II",
                    "Global and Local Stability",
                    "Visualization I -Phase Portraits, Null Clines",
                    "Visualization II",
                    "Toplogical Dynamics I",
                    "Topological Dynamics II",
                    "Graph Models",
                    "section{Group 4}",
                    "Lecture Test Day 1",
                    "Lecture Test Day 2",
                    "Lecture Test Day 3")

Course.1.Reading.List<-c("")

#---------------------Course 2------------------------
Course.2.Title<-c("Title: Design of Signal Transduction Networks and Gene Expression in Oncology Cell Lines")
Course.2.Outline<-c("section{Group 1}",
                    "NCI-60 Cell Line Model, Gene Expression, Compound Alignments", 
                    "Biotheorems on Gene Expression, Boundaries and Similarties",
                    "Motif Clustering Algorithms and Sequence Alignments",
                    "Signal Transduction, Projection and Compression",
                    "Amino Acid Based Chemical Modeling I and II",
                    "Gene Inhibitors and Silencers",
                    "Botany Models I - Species Modeling",
                    "Botany Models II - Shikimate Pathway I",
                    "Botany Models III - ShikiMate Pathway II",
                    "section{Group 2}",
                    "MAPK_MAPK14, ErbB_EGFR, Rap1_RRAS2",
                    "Calcium_MAPK14,cGMP_PRKACA, cAMP_PRKG1",
                    "Chemokine_FGR4, NFKappa_BIRC2, HIF1_EGFR",
                    "FOX0_Ubiquitin",
                    "Phosphatadylinositol_PLCD3,Sphingolipid_C06124",
                    "Cell.Cycle_CDC25A, Apoptsis_TP53, p53_MDM2",
                    "mTOR_MTOR,PI3K_AKT_AKT3,AMPK_PRKAG2",
                    "Longevity_SIRT1, Immune Systems",
                    "Wnt_APC2, Notch_RBPJL, Hedgehog_GLI1",
                    "TGF_Beta_SMAD6, VEGF_KDR",
                    "Apelin_GNB5, Hippo_YAP1",
                    "JakStat_JAK1,IL17.TAB3",
                    "Insulin_AKT3,Glucagon_PRKACA",
                    "WDR Proteins Model and Signal Network Reviews",
                    "section{Group 3}",
                    "Digestion Models and Dynamics I",
                    "Digestion Models and Dynamics II",
                    "Sugar and Insulin Models I",
                    "Sugar and Insulin Models II",
                    "Metabolic Rate Modeling",
                    "Obesity Models",
                    "Pancreatic Signaling I", 
                    "Pancreatic Signaling II",
                    "Pancreatic Signaling III",
                    "section{Group 4}",
                    "Lecture Test Day 1",
                    "Lecture Test Day 2",
                    "Lecture Test Day 3")

Course.2.Reading.List<-c("")

#-----------------Course 3 ---------------------------
Course.3.Title<-c("Classroom Lecture Model Series 3 Oncological Molecular Machine Learning with Combinatorial Topological Dynamics for Compound Discovery")
Course.3.Outline<-c("section{Group 1}",
                    "Ribosome Design, Engineering and Simulation",
                    "EndoPlasmic Rectilium Smooth and Rough Design",
                    "Translation and Transcription Error Models I",
                    "Translation and Transcription Error Models II",
                    "Translation and Transcription Error Models III",
                    "Tridiagonal Models", 
                    "Peptide Models I",
                    "Peptide Models II",
                    "Mitochrondria Models",
                    "section{Group 2}",
                    "Modal Analysis and Hamiltonian Dynamics", 
                    "Protein Interaction Models",
                    "Stochastic Search and Optimization in Compound Discovery",
                    "Molecular Machine Learning Models",
                    "Structure Activity Relationship Matrix Designs",
                    "EigenFrequency Classifiers",
                    "WDR Models",
                    "Molecular Machines I",
                    "Molecular Machines II",
                    "section{Group 3}",
                    "Chaperone Design",
                    "Chaperone Engineering",
                    "Chaperone Simulation",
                    "Proteasome Design",
                    "Proteasome Engineering",
                    "Proteasome Simulation",
                    "System  Models: Ribosome, Chaperone, Proteasome I",
                    "System  Models: Ribosome, Chaperone, Proteasome II",
                    "Course 1 2 3 Review",
                    "section{Group 4}",
                    "Lecture Test Day 1",
                    "Lecture Test Day 2",
                    "Lecture Test Day 3")

Course.3.Reading.List<-c("")

#----------------------------------Reading List Sources------------------------------------------------------#


#Reading_List_AIMS <- read_delim("Reading_List_AIMS.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
#Reading_List_Open_Access_Articles_Lectures <- read_delim("Reading_List_Open_Access_Articles_Lectures.txt", "\t", escape_double = FALSE, trim_ws = TRUE)


#-------------------------------------------------------------------------------------------------------------#

Data.1
Data.2
Data.3


f.1<-function(search.term, save.notes=FALSE)
{
  notes<-page_content("en","wikipedia", page_name =search.term)
  notes.filtered<-gsub("<[^>]+>", "", notes$parse$text$`*`)
  notes.content<-unlist(strsplit(notes.filtered," . "))
  temp<-unlist(notes.content)
  notes<-""
  for(j in 1:length(temp))
  {
    if(TRUE){notes<-stringr::str_c(notes,temp[j])}
    
  }
  notes.title<-stringr::str_c("Notes_For_Students_",search.term,".txt")
  
  if(save.notes)
  {
    write(notes, file=notes.title,append=FALSE)
  }
  output<-list()
  output$notes<-notes
  
  return(output)
  
}
test.f.1<-f.1("Ordinary Differential Equations",TRUE)
test.f.1


#-------------------------------------Function Library--------------------------------------------------------#

f.1<-function(X){ output<-list();output$X<-X;return(output)}
f.2<-function(X){ output<-list();output$X<-X;return(output)}
f.3<-function(X){ output<-list();output$X<-X;return(output)}


test.f.1<-f.1(X);test.f.1
test.f.2<-f.2(X);test.f.2
test.f.3<-f.3(X);test.f.3


#-------------------------------------Models------------------------------------------------------------------#

Classes.Titles<-list(Course.1.Title,Course.2.Title,Course.3.Title)

Lecture.Models.content<-list(Course.1.Outline,Course.2.Outline,Course.3.Outline)


#------------------------------------Write Models-------------------------------------------------------------#
setwd("Lecture_Designs/Course 1")

for(j in 1:length(Course.1.Outline))
{
  Lecture.Model.title<-stringr::str_c("Classroom_Lecture_Model_", Course.1.Outline[j])
  write(Lecture.Models.content[[1]],file=Lecture.Model.title)
}

setwd("Lecture_Designs/Course 2")

for(j in 1:length(Course.2.Outline))
{
  Lecture.Model.title<-stringr::str_c("Classroom_Lecture_Model_", Course.2.Outline[j])
  write(Lecture.Models.content[[2]],file=Lecture.Model.title)
}

setwd("Lecture_Designs/Course 3")

for(j in 1:length(Course.3.Outline))
{
  Lecture.Model.title<-stringr::str_c("Classroom_Lecture_Model_", Course.3.Outline[j])
  write(Lecture.Models.content[[3]],file=Lecture.Model.title)
}


#------------------------------------Lecture Management-------------------------------------------------------

Classes.Titles<-list(Course.1.Title,Course.2.Title,Course.3.Title)
for(i in 1:length(Classes.Titles))
{
  write(Classes.Titles[[i]],file=Classes.Titles)
}



#-------------------------------------Tables------------------------------------------------------------------#

Table.1
Table.2
Table.3
Table.4

#-------------------------------------Figures-----------------------------------------------------------------#

Figure.1
Figure.2
Figure.3
Figure.4




