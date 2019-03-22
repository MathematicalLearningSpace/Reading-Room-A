#---------------Reading List Design from aRXiv for Classroom Student Presentations-----------------#
library(aRxiv)

#---------------Examples from Lectures------------------------------------#
ode<-arxiv_search(query = 'ti:ordinary differential equations', limit=100)
dde<-arxiv_search(query = 'ti:delayed differential equations', limit=100)
sde<-arxiv_search(query = 'ti:stochastic differential equations', limit=100)
fde<-arxiv_search(query = 'ti:fractional differential equations', limit=100)
