#---------------Reading List Design from aRXiv for Classroom Student Presentations-----------------#
library(aRxiv)

#---------------Examples from Lectures-------------------------------------------#
ode<-arxiv_search(query = 'ti:ordinary differential equations', limit=100)
dde<-arxiv_search(query = 'ti:delayed differential equations', limit=100)
sde<-arxiv_search(query = 'ti:stochastic differential equations', limit=100)
fde<-arxiv_search(query = 'ti:fractional differential equations', limit=100)
#-------------------aRxiv Categories---------------------------------------------#
data(arxiv_cats)
#------------------Query categories----------------------------------------------#
#prefix 	explanation
#ti 	Title
#au 	Author
#abs 	Abstract
#co 	Comment
#jr 	Journal Reference
#cat 	Subject Category
#rn 	Report Number
#id 	Id (use id_list instead)
#all 	All of the above
#-------------------Pearson Distribution-----------------------------------------#
Pearson.Type.0<-arxiv_search(query = 'abs:Normal Distribution', limit=100)
Pearson.Type.I<-arxiv_search(query = 'abs:Beta Distribution',limit=100)
Pearson.Type.II<-arxiv_search(query = 'abs: Symmetric Beta Distribution',limit=100)
Pearson.Type.III<-arxiv_search(query = 'abs:Gamma Distribution',limit=100)
Pearson.Type.IV<-arxiv_search(query = 'abs: Pearson Type IV Distribution',limit=100) 
Pearson.Type.V<-arxiv_search(query = 'abs: Inverse Gamma Distribution',limit=100)
Pearson.Type.VI<-arxiv_search(query = 'abs: Beta Prime Distribution',limit=100)
Pearson.Type.VII<-arxiv_search(query = 'abs: Student t Distribution',limit=100)

Tikhonov.regularization<-arxiv_search(query = 'abs:Generalized Tikhonov regularization AND cat:math.ST',limit=100)
Gene.Expression<-arxiv_search(query = 'abs: Gene Expression AND cat:q-bio.BM',limit=100)
