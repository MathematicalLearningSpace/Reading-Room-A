#-------------------------------------------------------------------------#
#--------------------Classroom Lecture Model Series-----------------------#
#-------------------------------------------------------------------------#

#--------------------Work In Progress-------------------------------------#

#------------------------------R API----------------------------------#
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

library(easyPubMed);library(readr);library(CHNOSZ);
library(stringr);library(seqinr);library(seqLogo);library(msa);library(ape);
library(odseq);library(rphast)
library(plyr)

#---------------------------------------------------------------------#
#------------------------------Data-----------------------------------#
#---------------------------------------------------------------------#

#--------------------Review Notes------------------------#
#-----Signal Transduction Networks-----------------------#
A.xml<-fetch_pubmed_data(get_pubmed_ids("Ras signaling pathway")) 
B.xml<-fetch_pubmed_data(get_pubmed_ids("PI3K-Akt signaling pathway"))  
C.xml<-fetch_pubmed_data(get_pubmed_ids("Hippo signaling pathway")) 
D.xml<-fetch_pubmed_data(get_pubmed_ids("Wnt signaling pathway")) 
E.xml<-fetch_pubmed_data(get_pubmed_ids("p53 signaling pathway"))  
F.xml<-fetch_pubmed_data(get_pubmed_ids("Thyroid hormone signaling pathway")) 
G.xml<-fetch_pubmed_data(get_pubmed_ids("FoxO signaling pathway"))
H.xml<-fetch_pubmed_data(get_pubmed_ids("mTOR signaling pathway")) 
I.xml<-fetch_pubmed_data(get_pubmed_ids("ErbB signaling pathway"))  
J.xml<-fetch_pubmed_data(get_pubmed_ids("MAPK signaling pathway"))  
K.xml<-fetch_pubmed_data(get_pubmed_ids("Hippo signaling pathway"))  
L.xml<-fetch_pubmed_data(get_pubmed_ids("Apoptosis"))  
M.xml<-fetch_pubmed_data(get_pubmed_ids("Cell cycle")) 
N.xml<-fetch_pubmed_data(get_pubmed_ids("Rap1 signaling pathway"))  
O.xml<-fetch_pubmed_data(get_pubmed_ids("Chemokine signaling pathway"))
p.xml<-fetch_pubmed_data(get_pubmed_ids("Phosphatidylinositol signaling system"))
Q.xml<-fetch_pubmed_data(get_pubmed_ids("T cell receptor signaling pathway"))
R.xml<-fetch_pubmed_data(get_pubmed_ids("Estrogen signaling pathway"))
S.xml<-fetch_pubmed_data(get_pubmed_ids("VEGF signaling pathway"))

#---------------------------------------------------------------------#
#------------------------------Functions------------------------------#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#------------------------------Models---------------------------------#
#---------------------------------------------------------------------#

#---------------Signal Transduction------------#

#---------------QSAR Model I-------------------#

#---------------QSAR Model II------------------#


#---------------------------------------------------------------------#
#------------------------------Analysis-------------------------------#
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#------------------------------Tables---------------------------------#
#---------------------------------------------------------------------#

#------------Table 1-----------------#
#------------Table 2-----------------#
#------------Table 3-----------------#
#------------Table 4-----------------#



#---------------------------------------------------------------------#
#------------------------------Figures--------------------------------#
#---------------------------------------------------------------------#

#------------Figure 1-----------------#
#------------Figure 2-----------------#
#------------Figure 3-----------------#
#------------Figure 4-----------------#

