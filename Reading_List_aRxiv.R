#--------------------------------------------------------------------------------------------------#
#---------------Reading List Design from aRXiv for Classroom Student Presentations-----------------#
#--------------------------------------------------------------------------------------------------#
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
#---Mathematical Categories for Table Designs and Student Classroom Projects and Presentations----#
\hline 
math.AG & Mathematics - Algebraic Geometry & \\
math.AT & Mathematics - Algebraic Topology & \\
math.AP & Mathematics - Analysis of PDEs & \\
math.CT & Mathematics - Category Theory & \\
math.CA & Mathematics - Classical Analysis and ODEs & \\
math.CO & Mathematics - Combinatorics & \\
math.AC & Mathematics - Commutative Algebra & \\
math.CV & Mathematics - Complex Variables & \\
math.DG & Mathematics - Differential Geometry & \\
math.DS & Mathematics - Dynamical Systems & \\
math.FA & Mathematics - Functional Analysis & \\
math.GM & Mathematics - General Mathematics & \\
math.GN & Mathematics - General Topology & \\
math.GT & Mathematics - Geometric Topology & \\
math.GR & Mathematics - Group Theory & \\
math.HO & Mathematics - History and Overview & \\
math.IT & Mathematics - Information Theory & \\
math.KT & Mathematics - K-Theory and Homology & \\
math.LO & Mathematics - Logic & \\
math.MP & Mathematics - Mathematical Physics & \\
math.MG & Mathematics - Metric Geometry & \\
math.NT & Mathematics - Number Theory & \\
math.NA & Mathematics - Numerical Analysis & \\
math.OA & Mathematics - Operator Algebras & \\
math.OC & Mathematics - Optimization and Control & \\
math.PR & Mathematics - Probability & \\
math.QA & Mathematics - Quantum Algebra & \\
math.RT & Mathematics - Representation Theory & \\
math.RA & Mathematics - Rings and Algebras & \\
math.SP & Mathematics - Spectral Theory & \\
math.ST & Mathematics - Statistics & \\
\hline 

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


#----------------Student Assignments-----------------------------------------------#
#----------------To be Designed by Students----------------------------------------#
Module.1.Search<-function(X)
{
  
  output<-list()
  output$X<-X
  return(output)
}
test.Module.1.Search("1")
test.Module.1.Search
