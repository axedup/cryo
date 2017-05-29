
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
    }
    if(quali==1)
    {
      y <- lrm(formula(paste(y,"~",factor(x),sep="")),data=dat)
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
        res2 <- rbind(res2,resi)
      }
    }
    
    
    colnames(res2) <- c("Variable","OR (95\\%CI)","\\emph{P}")
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
cryo$renal_scorec<-as.factor(ifelse(cryo$renal_score<median(cryo$renal_score),paste0("< ",round(median(cryo$renal_score),0),""),
                                   paste0(">= ",round(median(cryo$renal_score),0),"")))


cryo$mrenal_scoref<-as.factor(cryo$mrenal_score)

cryo$mrenal_scorefc<-as.factor(ifelse(cryo$mrenal_score<median(cryo$mrenal_score),paste0("< ",round(median(cryo$mrenal_score),0),""),
                                    paste0(">= ",round(median(cryo$mrenal_score),0),"")))
cryo$volume<-as.numeric(as.character(cryo$volume))
cryo$volumec<-as.factor(ifelse(cryo$volume<median(cryo$volume),paste0("< ",round(median(cryo$volume),0),""),
                               paste0(">= ",round(median(cryo$volume),0),"")))

cryo$taille_maxc<-as.factor(ifelse(cryo$taille_max<median(cryo$taille_max),paste0("< ",round(median(cryo$taille_max),0),""),
                               paste0(">= ",round(median(cryo$taille_max),0),"")))


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
patients_cryo<-patients_cryo %>% distinct(nip, .keep_all = TRUE)

table(cryo$nip)
patient_cryo_descri<-TABKRIS(baz=patients_cryo,vect.var = c("age","agec","sex"),
                                vect.quali = c(0,1,1),
                                varint=NULL,valvarint =NULL,
                                nomvarint = "Histopathologic subtypes",
                                test=NULL,
                                vecnoms=c("Age","Age","Sexe"),valeurs=NULL,
                                vecrefs=NULL,varassoc=NULL,
                                codassoc=NULL,pres=NULL,langue="en",digits=2)

tumeurs_cryo<-TABKRIS(baz=cryo,vect.var = c("cote","histo","histof","gradef","pole","renal_score",
                                            "renal_scoref","renal_scorec","mrenal_score","mrenal_scoref","mrenal_scorefc",
                                            "location","renal_sinus","volume"
                                            ,"volumec","taille_max","taille_maxc","vasc","T"),
                      vect.quali = c(rep(1,5),0,1,1,0,1,1,1,1,0,1,0,1,1,1),
                      varint=NULL,valvarint =NULL,
                      nomvarint = "Histopathologic subtypes",
                      test=NULL,
                      vecnoms=c("Cote","Histo","Grade","Pole",
                                "renal score","renal score","renal score","mrenal score","mrenal score","mrenal score",
                                "location","renal sinus"
                                ,"volume","volume","taille max","taille max","Hyper ou hypovasculaire
","T"),valeurs=NULL,
                      vecrefs=NULL,varassoc=NULL,
                      codassoc=NULL,pres=NULL,langue="en",digits=2)


tumeurs_cryo_tt<-TABKRIS(baz=cryo,vect.var = c("nbr_a","nbr_af","nbr_afc","dissection","marges","margesc","clavien2","compli"),
                      vect.quali = c(0,1,1,1,0,1,1,1),
                      varint=NULL,valvarint =NULL,
                      nomvarint = "Histopathologic subtypes",
                      test=NULL,
                      vecnoms=c("Nombre aiguilles","Nombre aiguilles","Nombre aiguilles","Dissection","Marges","Marges","Clavien","Complications"),valeurs=NULL,
                      vecrefs=NULL,varassoc=NULL,
             
                     codassoc=NULL,pres=NULL,langue="en",digits=2)


varlist<-c("agec","sex","nbr_af","cote","histo","gradef","pole",
           "renal_scoref","mrenal_scoref","location","renal_sinus"
           ,"volume","taille_max","vasc","T","dissection","margesc","clavien2","nbr_af","dissection","margesc","clavien2","compli")
models <- lapply(varlist, function(x) {
  summary(glm(substitute(recidiven~ i, list(i = as.name(x))), data=cryo,family=binomial ))
})

capture.output(models, file=paste("C:/Users/adupont/Documents/projetstlouis/cryotherapie/","modeleul.txt",sep=""))




G<-NULL
for (i in 1:27){
  var<-c("age","agec","sex","cote","histof","pole","renal_score",
         "renal_scoref","renal_scorec","mrenal_score","mrenal_scoref","mrenal_scorefc",
         "location","renal_sinus","volume"
         ,"volumec","taille_max","taille_maxc","vasc","T","nbr_a","nbr_af","nbr_afc","dissection","marges","clavien2","compli")[i]
  t<-c(0,1,1,rep(1,3),0,1,1,0,1,1,1,1,0,1,0,1,1,1,0,1,1,1,0,1,1)[i]
  nom<-c("Age","Age","Sexe","Cote","Histo","Pole",
         "renal score","renal score","renal score","mrenal score","mrenal score","mrenal score",
         "location","renal sinus"
         ,"volume","volume","taille max","taille max","Hyper ou hypovasculaire
         ","T","Nombre aiguilles","Nombre aiguilles","Nombre aiguilles","Dissection","Marges","Clavien","Complications")[i]
  g<-disp.logist.u(y="recidiven",x=var, dat=cryo,quali=t,tex=TRUE,nomvar=nom, noms=NULL)
G<-rbind(G,g)
  
  }
  


"nbr_af","cote","histo","pole",
"renal_scoref","mrenal_scoref","location","renal_sinus"
,"volume","taille_max","vasc","T","dissection","margesc","nbr_af","dissection","margesc","clavien2","compli"



,"Nb af","Cote","Histo","Pole",
"Renal score","mrenal score","location","Renal sinus"
,"Volume","Taille max","Vascularisation","T","Dissection","Marges","Nbr af","Dissection","Marges","Clavien","Complication"