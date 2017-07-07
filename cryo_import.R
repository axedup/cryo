library(survival)
library(gdata)
library(Hmisc)
library(compareGroups)
library(stargazer)
library(xtable)
library(dplyr)
library(ggplot2)
library(cmprsk)


format.pv <- function(p, text=F)
{
  if(p<0.0001) return("<0.0001")
  if(p>=0.0001&p<0.00095) ifelse(text==F,return(sprintf("%.4f", p)),return(paste("=",sprintf("%.4f", p),sep="")))
  if(p>=0.00095&p<=0.0095) ifelse(text==F,return(as.character(signif(p,1))),return(paste("=",as.character(signif(p,1)),sep="")))
  if(p>0.0095&p<0.0995) ifelse(text==F,return(sprintf("%.3f", signif(p,2))),return(paste("=",sprintf("%.3f", signif(p,2)),sep="")))
  if(p>=0.0995) ifelse(text==F,return(sprintf("%.2f", signif(p,2))),return(paste("=",sprintf("%.2f", signif(p,2)),sep="")))
}
format.hr <- function(z)
{
  if(z<0.05) return(sprintf("%.3f",z))
  if(z<=9.95&z>=0.05) return(sprintf("%.2f",z))
  if(z>9.95) return(sprintf("%.1f",z))
}
or.calc <- function(obj, alpha=0.05)
{
  se <- sqrt(diag(obj$var))
  res <- data.frame(OR=exp(obj$coef), lower=exp(obj$coef-qnorm(1-alpha/2)*se), upper=exp(obj$coef+qnorm(1-alpha/2)*se), p.value=(1-pchisq((obj$coef/se)^2,1)))
  return(res)
}


or.calc2 <- function(obj, alpha=0.05)
{
  se <- obj$coefficients[3,2]
  res <- data.frame(OR=round(exp(obj$coefficients[3,1]),2), lower=round(exp(obj$coefficients[3,1]-qnorm(1-alpha/2)*se),2), upper=round(exp(obj$coefficients[3,1]+qnorm(1-alpha/2)*se),2), p.value=as.numeric(as.character(format.pv(obj$coefficients[3,4]))))
  return(res)
}




# nomsvar <- vecnom[i]
# quali <- quali[i]
# noms <- valeurs[i]

# la version originale de Raphaël
# avec un paramètre qui RES
disp.logist.u <- function(x, y, dat, res, quali=0, nomvar=NULL, noms=NULL, tex=T)
{
  require(rms)
  if(is.null(nomvar))
    nomvar <- x
  if(tex==TRUE){
    if(quali==0)
    {
      y <- lrm(formula(paste(y,"~",x,sep="")),data=dat)
      or <- or.calc(y)
      res2 <- rbind(c(nomvar, paste(format.hr(or[2,1])," (",format.hr(or[2,2]),"--", format.hr(or[2,3]), ")",sep=""), format.pv(or[2,4])))
   pg<-c("-")
       res2<-cbind(res2,pg)
      }
    if(quali==1)
    {
      y <- lrm(formula(paste(y,"~",factor(x),sep="")),data=dat)
      # tat<-model.matrix(~dat[,x])
      # colnames(tat)<-c("","tu","su")
      # dat2<-cbind(dat,tat)
      # d<-glm(formula(paste(yy,"~",dat2$tu,"+",dat2$su,sep="")),data=dat2,family="binomial")
      # ano<-anova(d,test="Chisq")
      
      if(is.null(noms))
      {
        aa<-table(dat[,match(x,names(dat))])
        noms<- rownames(aa)
      }
      res2 <- rbind(c(nomvar,"",""),c(paste("\\quad ", noms[1], sep=""),"1",""))
      or <- or.calc(y)
      for(i in 2:length(y$coef))
      {
        resi <- c(paste("\\quad ", noms[i], sep=""),paste(format.hr(or[i,1])," (",format.hr(or[i,2]),"--", format.hr(or[i,3]), ")",sep=""), format.pv(or[i,4]))
        #resi <- c(paste("\\quad ", noms[i], sep=""),paste(format.hr(or[i,1])," (",format.hr(or[i,2]),"--", format.hr(or[i,3]), ")",sep=""), round(ano[i,5],2))
        res2 <- rbind(res2,resi)
      }
      
      pg<-c(round(y$stats[5],2),rep("",nrow(res2)-1))
      res2<-cbind(res2,pg)
    }
    
    
    colnames(res2) <- c("Variable","OR (95\\%CI)","\\emph{P}","P")
    return(res2)
  }
  if(tex==FALSE){
    print(x)
    if(quali==0)
    {
      y <- lrm(formula(paste(y,"~",x,sep="")),data=dat)
      or <- or.calc(y)
      res2 <- rbind(c(nomvar, paste(format.hr(or[2,1])," (",format.hr(or[2,2]),"-", format.hr(or[2,3]), ")",sep=""), format.pv(or[2,4])))
    }
    if(quali==1)
    {
      y <- lrm(formula(paste(y,"~",factor(x),sep="")),data=dat)
      print(y)
      if(is.null(noms))
      {
        aa<-table(dat[,match(x,names(dat))])
        noms<- rownames(aa)
      }
      res2 <- rbind(c(nomvar,"",""),c(paste(" ", noms[1], sep=""),"1",""))
      or <- or.calc(y)
      for(i in 2:length(y$coef))
      {
        resi <- c(paste(" ", noms[i], sep=""),paste(format.hr(or[i,1])," (",format.hr(or[i,2]),"-", format.hr(or[i,3]), ")",sep=""), format.pv(or[i,4]))
        res2 <- rbind(res2,resi)
      }
    }
    
    
    #colnames(res2) <- c("Variable","OR (95%CI)","P")
    return(res2)
  }
}




cryo<-read.csv2("C:/Users/adupont/Documents/projetstlouis/cryotherapie/cryo.csv",na.strings = c("unknown","?"))

cryo$clavien2<-as.factor(ifelse(cryo$clavien %in% c(0,1),"0-1","2-3"))
cryo$nbr_af<-as.factor(cryo$nbr_a)
cryo$nbr_afc<-as.factor(ifelse(cryo$nbr_a<median(cryo$nbr_a),paste0("< ",round(median(cryo$nbr_a),0),""),
                                    paste0(">= ",round(median(cryo$nbr_a),0),"")))


cryo$histof<-as.factor(ifelse(cryo$histo %in% c("CCR"),"CCR","Non-CCR"))

cryo$gradef<-as.factor(cryo$grade)
cryo$renal_scoref<-as.factor(cryo$renal_score)
cryo$renal_scorec<-ifelse(cryo$renal_score<7, "simple","interm�diaire")
cryo$renal_scorec<-ifelse(cryo$renal_score>9, "complexe",cryo$renal_scorec)
cryo$renal_scorec<-as.factor(cryo$renal_scorec)


cryo$mrenal_scoref<-as.factor(cryo$mrenal_score)

cryo$mrenal_scorec<-ifelse(cryo$mrenal_score<7, "simple","interm�diaire")
cryo$mrenal_scorec<-ifelse(cryo$mrenal_score>9, "complexe",cryo$mrenal_scorec)
cryo$mrenal_scorec<-as.factor(cryo$mrenal_scorec)





cryo$volume<-as.numeric(as.character(cryo$volume))
cryo$volumec<-as.factor(ifelse(cryo$volume<median(cryo$volume),paste0("< ",round(median(cryo$volume),0),""),
                               paste0(">= ",round(median(cryo$volume),0),"")))

cryo$taille_maxc<-as.factor(ifelse(cryo$taille_max<median(cryo$taille_max),paste0("< ",round(median(cryo$taille_max),0),""),
                               paste0(">= ",round(median(cryo$taille_max),0),"")))
cryo$dissectionc<-as.factor(ifelse(cryo$dissection %in% c("0"),"non","oui"))
table(cryo$dissectionc,exclude=NULL)

cryo$tt_incom<-as.factor(cryo$tt_incom)
cryo$recidiven<-ifelse(cryo$recidive %in% c("OUI", "OUI "),1,0)
cryo$agec<-as.factor(ifelse(cryo$age<median(cryo$age),paste0("< ",round(median(cryo$age),0),"ans"),paste0(">= ",round(median(cryo$age),0),"ans")))

cryo$margesc<-as.factor(ifelse(cryo$marges >0,"oui","non"))
table(cryo$margesc)
levels(cryo$compli)<-c("non", "non", "non", "oui")


levels(cryo$renal_sinus)<-c("No", "No", "Yes", "Yes")
levels(cryo$vasc)<-c("Hyper", "Hyper", "Hypo", "Hypo")

e<-cryo$nip[cryo$recidive=="OUI"]
cryo[cryo$nip %in% e,]
# patient id : 3214055737 2 cryo dont une tumeur qui récidive on supprime celle qui ne récidive pas #

patients_cryo<-cryo[order(cryo$nip,cryo$age),]
length(unique(cryo$nip)) # 169 patients
patients_cryo<-patients_cryo %>% select(everything())%>%  distinct(nip,.keep_all=TRUE)

table(cryo$nip)
table(patients_cryo$nip)
patient_cryo_descri<-TABKRIS(baz=patients_cryo,vect.var = c("age","sex"),
                                vect.quali = c(0,1),
                                varint=NULL,valvarint =NULL,
                                nomvarint = "Histopathologic subtypes",
                                test=NULL,
                                vecnoms=c("Age","Sexe"),valeurs=NULL,
                                vecrefs=NULL,varassoc=NULL,
                                codassoc=NULL,pres=c("mean_mm",""),langue="en",digits=0)

tumeurs_cryo<-TABKRIS(baz=cryo,vect.var = c("cote","histo","histof","gradef","pole","renal_score",
                                            "renal_scoref","renal_scorec","mrenal_score","mrenal_scoref","mrenal_scorec",
                                            "location","renal_sinus","volume"
                                            ,"taille_max","vasc","T"),
                      vect.quali = c(rep(1,5),0,1,1,0,1,1,1,1,0,0,1,1),
                      varint=NULL,valvarint =NULL,
                      nomvarint = "Histopathologic subtypes",
                      test=NULL,
                      vecnoms=c("Cote","Histo","Histo","Grade","Pole",
                                "renal score","renal score","renal score","mrenal score","mrenal score","mrenal score",
                                "location","renal sinus"
                                ,"volume","taille max","Hyper ou hypovasculaire
","T"),valeurs=NULL,
                      vecrefs=NULL,varassoc=NULL,
                      codassoc=NULL,pres=c(rep("",5),"mean_mm","","","mean_mm","","","","","mean_mm",
                                           "mean_mm","","")
                        ,langue="en",digits=0)


tumeurs_cryo_tt<-TABKRIS(baz=cryo,vect.var = c("nbr_a","nbr_af","nbr_afc","dissection","marges","margesc","clavien2"),
                      vect.quali = c(0,1,1,1,0,1,1),
                      varint=NULL,valvarint =NULL,
                      nomvarint = "Histopathologic subtypes",
                      test=NULL,
                      vecnoms=c("Nombre aiguilles","Nombre aiguilles","Nombre aiguilles","Dissection","Marges","Marges","Clavien"),valeurs=NULL,
                      vecrefs=NULL,varassoc=NULL,
             
                     codassoc=NULL,pres=c("mean_mm",1,1,1,"mean_mm",1,1),langue="en",digits=0)


tumeurs_cryo_out<-TABKRIS(baz=cryo,vect.var = c("compli","recidiven"),
                         vect.quali = c(1,1),
                         varint=NULL,valvarint =NULL,
                         nomvarint = "Histopathologic subtypes",
                         test=NULL,
                         vecnoms=c("Complications","Recidive"),valeurs=NULL,
                         vecrefs=NULL,varassoc=NULL,
                         
                         codassoc=NULL,pres=NULL,langue="en",digits=0)




biv_patients<-TABKRIS(baz=cryo,vect.var = c("age","agec","sex","cote","histo","histof","gradef","pole","renal_score",
                                                     "renal_scoref","renal_scorec","mrenal_score","mrenal_scoref",
                                            "mrenal_scorec",
                                                     "location","renal_sinus","volume"
                                                     ,"volumec","taille_max","taille_maxc","vasc","T","nbr_a",
                                            "nbr_af","nbr_afc","dissection","marges","margesc","clavien2","compli"
                                            ),
                             vect.quali = c(0,1,1,rep(1,5),0,1,1,0,1,1,1,1,0,1,0,1,1,1,0,1,1,1,0,1,1,1),
                             varint=c("recidiven"),valvarint =c("Pas d'echec","Echec"),
                             nomvarint = "",
                             test=NULL,
                             vecnoms=c("Age","Age","Sexe","Cote","Histo","Histo","Grade","Pole",
                                       "renal score","renal score","renal score","mrenal score","mrenal score","mrenal score",
                                       "location","renal sinus"
                                       ,"volume","volume","taille max","taille max","Hyper ou hypovasculaire
                                ","T","Nombre aiguilles","Nombre aiguilles","Nombre aiguilles","Dissection",
                                       "Marges","Marges","Clavien","Complication"),valeurs=NULL,
                             vecrefs=NULL,varassoc=NULL,
                             codassoc=NULL,pres=NULL,langue="en",digits=2)


varlist<-c("agec","sex","nbr_af","cote","histo","gradef","pole",
           "renal_scoref","mrenal_scoref","location","renal_sinus"
           ,"volume","taille_max","vasc","T","dissection","margesc","clavien2","nbr_af","dissection","margesc","clavien2","compli")
models <- lapply(varlist, function(x) {
  summary(glm(substitute(recidiven~ i, list(i = as.name(x))), data=cryo,family=binomial ))
})

capture.output(models, file=paste("C:/Users/adupont/Documents/projetstlouis/cryotherapie/","modeleul.txt",sep=""))



test<-glm(recidiven~ age, data=cryo,family=binomial )
1 - pchisq(112.76, 179)
test2<-glm(recidiven~ factor(age), data=cryo,family=binomial )
anova(test,test2,test="Chisq")

test<-glm(recidiven~ renal_score, data=cryo,family=binomial )
1 - pchisq(112.77, 180)
test2<-glm(recidiven~ factor(age), data=cryo,family=binomial )
anova(test,test2,test="Chisq")




# b<-c(0.01,0.48,0.71,0.95,1.19,0.01,0.48,1.44,0.71,1.96)
#  a<-c(127.6,124.0,110.8,103.9,101.5,130.1,122.0,92.3,113.1,83.7)
#  corrosion<-data.frame(a,b)
# t<-lm(a~ b, data=corrosion )
# p<-lm(a~ factor(b), data=corrosion )
# y<-glm(a~ b, data=corrosion )
# u<-glm(a~ factor(b), data=corrosion )
# 
# anova(y,u,test="F")
# anova(t,p)

# G<-NULL
# for (i in 1:24){
#   var<-c("age","agec","sex","cote","histof","pole","renal_score"
#          ,"renal_scorec","mrenal_score","mrenal_scorec",
#          "location","renal_sinus","volume"
#          ,"volumec","taille_max","taille_maxc","vasc","T","nbr_a","nbr_afc","dissectionc","marges","clavien2","compli")[i]
#   t<-c(0,1,1,rep(1,3),0,1,0,1,1,1,0,1,0,1,1,1,0,1,1,1,1,1)[i]
#   nom<-c("Age","Age","Sexe","Cote","Histo","Pole",
#          "renal score","renal score","mrenal score","mrenal score",
#          "location","renal sinus"
#          ,"volume","volume","taille max","taille max","Hyper ou hypovasculaire
#          ","T","Nombre aiguilles","Nombre aiguilles","Dissection","Marges","Clavien","Complications")[i]
#   g<-disp.logist.u(y="recidiven",x=var, dat=cryo,quali=t,tex=TRUE,nomvar=nom, noms=NULL)
# G<-rbind(G,g)
#   
#   }
#   
# m<-glm(recidiven~ margesc, data=cryo,family=binomial )
# u<-summary(glm(recidiven~ margesc, data=cryo,family=binomial ))
# marges<-data.frame(Variable=c("marges","\\quad oui","\\quad non"),OR=
#                      c("","1", paste(round(exp(u$coefficients)[2,1],2),"(",
#                                                   round(exp(confint(m))[2,1],3) ,"--", 
#                                                   
#                                                   round(exp(confint(m))[2,2],3),")"
#                                                   )),p=c("","","<0.0001"))     
#                      
# colnames(marges)<-colnames(G)
# G<-rbind(G,marges)


# lmm <- lmer(recidiven ~ sex  + (1 | nip), data = cryo,REML = TRUE)
# 
# cryo$nip2<-c(1,rep(1:90,2))
# cryo$recideven5<-c(1,rep(0,95),rep(1,85))
# glmer(recidiven ~ agec  + sex+(1 | nip2), data = cryo, family = binomial, control = glmerControl(optimizer = "bobyqa"))
# 
# glm(recidiven~ sex+agec, data=cryo,family=binomial )
# 

# "nbr_af","cote","histo","pole",
# "renal_scoref","mrenal_scoref","location","renal_sinus"
# ,"volume","taille_max","vasc","T","dissection","margesc","nbr_af","dissection","margesc","clavien2","compli"
# 
# 
# 
# ,"Nb af","Cote","Histo","Pole",
# "Renal score","mrenal score","location","Renal sinus"
# ,"Volume","Taille max","Vascularisation","T","Dissection","Marges","Nbr af","Dissection","Marges","Clavien","Complication"



G2<-NULL
for (i in 1:18){
  var<-c("age","sex","cote","histof","pole","renal_scorec"
         ,"mrenal_scorec",
         "location","renal_sinus","volume",
         "taille_max","vasc","T","nbr_a","nbr_afc","dissectionc","clavien2","compli")[i]
  t<-c(0,1,1,rep(1,2),1,1,1,1,0,0,1,1,0,1,1,1,1,1)[i]
  nom<-c("Age","Sexe","Cote","Histo","Pole",
         "renal score","mrenal score",
         "location","renal sinus"
         ,"volume","taille max","Hyper ou hypovasculaire
         ","T","Nombre aiguilles","Nombre aiguilles","Dissection","Marges","Clavien","Complications")[i]
  g<-disp.logist.u(y="recidiven",x=var, dat=cryo,quali=t,tex=TRUE,nomvar=nom, noms=NULL)
  G2<-rbind(G2,g)
  
}

m<-glm(recidiven~ margesc, data=cryo,family=binomial )
u<-summary(glm(recidiven~ margesc, data=cryo,family=binomial ))
marges<-data.frame(Variable=c("marges","\\quad oui","\\quad non"),OR=
                     c("","1", paste(round(exp(u$coefficients)[2,1],2),"(",
                                     round(exp(confint(m))[2,1],3) ,"--", 
                                     
                                     round(exp(confint(m))[2,2],3),")"
                     )),p=c("","","<0.0001"),pp=c("<0.0001","",""))     

colnames(marges)<-colnames(G2)
G2<-rbind(G2,marges)


colnames(G2)<-c("Variable","OR (95\\%CI)", "\\emph{P}","P maximum vrais" )

sexe<-summary(glm(recidiven~ margesc+sex, data=cryo,family=binomial ))
cote<-summary(glm(recidiven~ margesc+cote, data=cryo,family=binomial ))
histo<-summary(glm(recidiven~ margesc+histo, data=cryo,family=binomial ))
renal<-summary(glm(recidiven~ margesc+renal_sinus, data=cryo,family=binomial ))
vasc<-summary(glm(recidiven~ margesc+vasc, data=cryo,family=binomial ))
T<-summary(glm(recidiven~ margesc+T, data=cryo,family=binomial ))
aig<-summary(glm(recidiven~ margesc+nbr_afc, data=cryo,family=binomial ))
diss<-summary(glm(recidiven~ margesc+dissection, data=cryo,family=binomial ))
clavien<-summary(glm(recidiven~ margesc+clavien2, data=cryo,family=binomial ))
pole<-summary(glm(recidiven~ margesc+pole, data=cryo,family=binomial ))

age<-summary(glm(recidiven~ margesc+age, data=cryo,family=binomial ))
1-pchisq(112.770 - 19.682, df=2)


renal_score<-summary(glm(recidiven~ margesc+renal_score, data=cryo,family=binomial ))
1-pchisq(112.770 - 19.77, df=2)


mrenal_score<-summary(glm(recidiven~ margesc+mrenal_score, data=cryo,family=binomial ))
1-pchisq(112.770 - 19.77, df=2)

volume<-summary(glm(recidiven~ margesc+volume, data=cryo,family=binomial ))
1-pchisq(112.770 - 19.77, df=2)

taille<-summary(glm(recidiven~ margesc+taille_max, data=cryo,family=binomial ))
1-pchisq(112.770 - 19.77, df=2)


### 
cryo$Recurrence<-as.factor(cryo$recidiven)
levels(cryo$Recurrence)<-c("No recurrence","Recurrence")
histog<- ggplot(cryo, aes(as.factor(marges),fill=Recurrence))+geom_bar()+
  scale_fill_manual("legend", values = c("No recurrence" = "#0066CC", "Recurrence" = "red"))+
xlab("Margins (mm)")+ylab("Number")
  theme(legend.title=element_blank())

histog
