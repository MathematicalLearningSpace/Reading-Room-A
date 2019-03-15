#-------------------------------------Classroom Lecture Model Math Notes Management System---------------------------------------
#-------------------------------------R Package Management for Students in the Classroom-----------------------------------------

#-----------------------------------------R Data------------------------------------------------------------------------------------
db <- hsearch_db()
R.package.Description<-lapply(unique(db$Concepts$Package), function(x) packageDescription(x,fields = c("Title", "Description")))
R.package.Description
#----------------------------Function To Be Developed in the Classroom------------------------------------------------------------
search.1<-function(keywords,R.package.Description)
{
R.Titles<-NULL
R.Description<-NULL
R.Package<-NULL
R.References<-""
R.Citations<-""
for(i in 1:length(R.package.Description))
{
  R.Titles[i]<-R.package.Description[[i]]$Title
  R.Description[i]<-R.package.Description[[i]]$Description
  R.Package[i]<-R.package.Description[[i]]$Package
}
R.Titles<-R.Titles[grep(keywords,R.Description)]
R.Description<-R.Description[grep(keywords,R.Description)]
R.Package<-R.Package[grep(keywords,R.Description)]
for(j in 1:length(R.Package))
{
  R.References<-stringr::str_c(R.References,
                               "\\Library(",
                              R.Package[j],
                              ")")
  R.Citations<-stringr::str_c(R.Citations,R.Package[j])
}
output<-list()
output$Title<-R.Titles
output$Description<-R.Description
output$Package<-R.Package
output$References<-R.References
output$Citations<-R.Citations
return(output)
}
test.search.1<-search.1("differential equations",R.package.Description)
test.search.1 
#----------------------------------------Tables---------------------------------------------------------------------------------
Table.1.df<-cbind(db$Concepts$Concept,db$Concepts$Package)
Table.1.df
#----------------------------------------Figures--------------------------------------------------------------------------------
