
Arrhenius.equation.models.prototype<-function(X,PM,T.1)
{
  
  #--------------------------------------A------------------------------------------
  A.E.1<-NULL;A.E.2<-NULL;A.E.3<-NULL;A.E.4<-NULL;A.E.5<-NULL;A.E.6<-NULL;
  A.E.7<-NULL;A.E.8<-NULL;A.E.9<-NULL;A.E.10<-NULL;A.E.11<-NULL;A.E.12<-NULL;
  A.E.13<-NULL;A.E.14<-NULL;A.E.15<-NULL;A.E.16<-NULL;A.E.17<-NULL;A.E.18<-NULL;
  A.E.19<-NULL;A.E.20<-NULL;
  
  Table.1<-data.frame()
  Table.2<-data.frame()
  Table.3<-data.frame()
   
  h.1<-function(X)
  {
    Table.TeX<-""
    column.1<-""
    row.1<-""
    
    for(i in 1:nrow(X))
    {
      for(j in 1:ncol(X))
      {
        column.1<-stringr::str_c(column.1,X[i,j]," & ")
      }
      row.1<-stringr::str_c(row.1,column.1," \\")
    }
    Table.TeX<-row.1
    return(Table.TeX)
  }
  #----------------------------------------B----------------------------------
  for(t in 1:T.1)
  {
    A.E.1[t]<-PM[1,1]+PM[1,2]*exp(-PM[1,3]*t)^(PM[1,4])
    A.E.2[t]<-PM[2,1]+PM[2,2]*exp(-PM[2,3]*t)^(PM[2,4])
    A.E.3[t]<-PM[3,1]+PM[3,2]*exp(-PM[3,3]*t)^(PM[3,4])
    A.E.4[t]<-PM[4,1]+PM[4,2]*exp(-PM[4,3]*t)^(PM[4,4])
    A.E.5[t]<-PM[5,1]+PM[5,2]*exp(-PM[5,3]*t)^(PM[5,4])
    A.E.6[t]<-PM[6,1]+PM[6,2]*exp(-PM[6,3]*t)^(PM[6,4])
    A.E.7[t]<-PM[7,1]+PM[7,2]*exp(-PM[7,3]*t)^(PM[7,4])
    A.E.8[t]<-PM[8,1]+PM[8,2]*exp(-PM[8,3]*t)^(PM[8,4])
    A.E.9[t]<-PM[9,1]+PM[9,2]*exp(-PM[9,3]*t)^(PM[9,4])
    A.E.10[t]<-PM[10,1]+PM[10,2]*exp(-PM[10,3]*t)^(PM[10,4])
    A.E.11[t]<-PM[11,1]+PM[11,2]*exp(-PM[11,3]*t)^(PM[11,4])
    A.E.12[t]<-PM[12,1]+PM[12,2]*exp(-PM[12,3]*t)^(PM[12,4])
    A.E.13[t]<-PM[13,1]+PM[13,2]*exp(-PM[13,3]*t)^(PM[13,4])
    A.E.14[t]<-PM[14,1]+PM[14,2]*exp(-PM[14,3]*t)^(PM[14,4])
    A.E.15[t]<-PM[15,1]+PM[15,2]*exp(-PM[15,3]*t)^(PM[15,4])
    A.E.16[t]<-PM[16,1]+PM[16,2]*exp(-PM[16,3]*t)^(PM[16,4])
    A.E.17[t]<-PM[17,1]+PM[17,2]*exp(-PM[17,3]*t)^(PM[17,4])
    A.E.18[t]<-PM[18,1]+PM[18,2]*exp(-PM[18,3]*t)^(PM[18,4])
    A.E.19[t]<-PM[19,1]+PM[19,2]*exp(-PM[19,3]*t)^(PM[19,4])
    A.E.20[t]<-PM[20,1]+PM[20,2]*exp(-PM[20,3]*t)^(PM[20,4])
  }
  #-------------------------------------C-------------------------------------------
  output<-list()
  output$X
  output$A.E.1<-A.E.1
  output$A.E.2<-A.E.2
  output$A.E.3<-A.E.3
  output$A.E.4<-A.E.4
  output$A.E.5<-A.E.5
  output$A.E.6<-A.E.6
  output$A.E.7<-A.E.7
  output$A.E.8<-A.E.8
  output$A.E.9<-A.E.9
  output$A.E.10<-A.E.10
  output$A.E.11<-A.E.11
  output$A.E.12<-A.E.12
  output$A.E.13<-A.E.13
  output$A.E.14<-A.E.14
  output$A.E.15<-A.E.15
  output$A.E.16<-A.E.16
  output$A.E.17<-A.E.17
  output$A.E.18<-A.E.18
  output$A.E.19<-A.E.19
  output$A.E.20<-A.E.20
  
  Table.1<-rbind(format(output$A.E.1,digits=2),
                 format(output$A.E.2,digits=2),
                 format(output$A.E.3,digits=2),
                 format(output$A.E.4,digits=2),
                 format(output$A.E.5,digits=2),
                 format(output$A.E.6,digits=2),
                 format(output$A.E.7,digits=2),
                 format(output$A.E.8,digits=2),
                 format(output$A.E.9,digits=2),
                 format(output$A.E.10,digits=2),
                 format(output$A.E.11,digits=2),
                 format(output$A.E.12,digits=2),
                 format(output$A.E.13,digits=2),
                 format(output$A.E.14,digits=2),
                 format(output$A.E.15,digits=2),
                 format(output$A.E.16,digits=2),
                 format(output$A.E.17,digits=2),
                 format(output$A.E.18,digits=2),
                 format(output$A.E.19,digits=2),
                 format(output$A.E.20,digits=2))
  
  Table.1.TeX<-h.1(Table.1)
  Table.2.TeX<-h.1(PM)
  output$Table.1<-Table.1
  output$Table.2<-Table.2
  output$Table.3<-Table.3
  output$Table.1.TeX<-Table.1.TeX
  output$Table.2.TeX<-Table.2.TeX
  output$Reference<-stringr::str_c("\\bibitem[1]{key100} Ilknur ALIBAS (2014)",
                                   "\\newblock Mathematical modeling of microwave dried celery leaves and determination of the effective moisture diffusivities and activation energy", 
                                   "\\newblock Food Sci. Technol, Campinas, 34(2): 394-401.")
  
  return(output)
  
}
#---------------------------------------------------Test Functions in the Mathematical Learning Space------------------------------------------------
X<-runif(10^2,0,10^0)
parameter.matrix<-matrix(1,20,6)
Arrhenius.equation.models.prototype(X,parameter.matrix,6)
png(file = stringr::str_join("Figures//Figure","1",".png"))
plot(Arrhenius.equation.models(X,parameter.matrix,6)$A.E.3)
dev.off()
#-----------------------------------------------------A------------------------------------------------------------------------------------
logistic.1<-function(T.1,A=1,K=0,C=1,Q=1,B=1,M=1,v=1)
{
  Table.1<-""
  logistic<-NULL
  #\begin{equation}f(x) = \frac 1 {1+e^{-x}}= \frac {e^{x}} {1+e^{x}}\\end{equation} \cite{key1}
  #\begin{equation}\frac{d}{dx}f(x) = \frac{e^{x} \cdot (1+e^{x})-e^{x} \cdot e^{x}}{(1 + e^{x})^2} =\end{equation} \cite{key2}
  #\begin{equation}\frac{e^{x}}{(1 + e^{x})^2} = f(x)(1-f(x))\end{equation}\cite{key3}
  #\begin{equation}\frac{d}{dx}f(x) = f(x)(1-f(x)) \end{equation} \cite{key4}
  
  for(t in 1:T.1)
  {
    logistic[t]<-A + (K-A)/(C + Q*exp(-B*(t-M))^(1/v))
  }
  
  output<-list()
  output$A<-A
  output$K<-K
  output$C<-C
  output$Q<-Q
  output$B<-B
  output$M<-M
  output$v<-v
  output$logistic<-logistic
  output$Table.1<-Table.1
  output$Reference<-stringr::str_c("\\bibitem[1]{key100}Wikipedia contributors.", 
                                   "\\newblock Generalised logistic function.",
                                   "\\newblock Wikipedia, The Free Encyclopedia. Wikipedia, The Free Encyclopedia, 23 Jun. 2018. Web. 1 Nov. 2018.") 
  return(output)
}

#-------------------------------------------------------------Test---------------------------------------------------
logistic.1(10)

plot_colors <- c("blue","black", "green", "orange", "pink")
png(file = stringr::str_join("Figures//Figure","1",".png"))
plot(seq(0,100,1),logistic.1(101)$logistic,type="l")
legend(x = "top",inset = 0,
       legend = c("Logistic"), 
       col=plot_colors, lwd=5, cex=.5, horiz = TRUE)
dev.off()
