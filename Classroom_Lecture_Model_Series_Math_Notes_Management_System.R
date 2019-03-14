#-------------------------------------Classroom Lecture Model Math Notes Management System-----------------------


#-----------------------------------------R Package Management for Students in the Classroom--------------------------------------

#-----------------------------------------Data------------------------------------------------------------------------------------
db <- hsearch_db()
R.package.Description<-lapply(unique(db$Concepts$Package), function(x) packageDescription(x,fields = c("Title", "Description")))
R.package.Description


#----------------------------------------Tables---------------------------------------------------------------------------------
Table.1.df<-cbind(db$Concepts$Concept,db$Concepts$Package)
Table.1.df



#----------------------------------------Figures--------------------------------------------------------------------------------