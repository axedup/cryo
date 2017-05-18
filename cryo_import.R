cryo<-read.csv2("C:/Users/adupont/Documents/projetstlouis/cryotherapie/cryo.csv",na.strings = c("unknown","?"))

cryo$clavien2<-ifelse(cryo$clavien %in% c(0,1),"0-1","2-3")
cryo$nbr_af<-as.factor(cryo$nbr_a)
cryo$gradef<-as.factor(cryo$grade)
cryo$renal_scoref<-as.factor(cryo$renal_score)
cryo$mrenal_scoref<-as.factor(cryo$mrenal_score)
cryo$volume<-as.numeric(as.character(cryo$volume))
cryo$tt_incom<-as.factor(cryo$tt_incom)
cryo$recidiven<-ifelse(cryo$recidive %in% c("OUI", "OUI "),1,0)


cryo$margesc<-ifelse(cryo$marges >0,"oui","non")
table(cryo$margesc)
levels(cryo$compli)<-c("non", "non", "non", "oui")

e<-cryo$nip[cryo$recidive=="OUI"]
cryo[cryo$nip %in% e,]
# patient id : 3214055737 2 cryo dont une tumeur qui récidive on supprime celle qui ne récidive pas #

patients_cryo<-cryo[order(cryo$nip,cryo$age),]
length(unique(cryo$nip)) # 169 patients
patients_cryo<-patients_cryo %>% distinct(nip, .keep_all = TRUE)

table(cryo$nip)
patient_cryo_descri<-TABKRIS(baz=patients_cryo,vect.var = c("age","sex"),
                                vect.quali = c(0,1),
                                varint=NULL,valvarint =NULL,
                                nomvarint = "Histopathologic subtypes",
                                test=NULL,
                                vecnoms=c("Age","Sexe"),valeurs=NULL,
                                vecrefs=NULL,varassoc=NULL,
                                codassoc=NULL,pres=NULL,langue="en",digits=2)

tumeurs_cryo<-TABKRIS(baz=cryo,vect.var = c("nbr_af","cote","histo","gradef","pole",
                                            "renal_scoref","mrenal_scoref","location","renal_sinus"
                                            ,"volume","taille_max","vasc","T","dissection","margesc","clavien2"),
                      vect.quali = rep(1,16),
                      varint=NULL,valvarint =NULL,
                      nomvarint = "Histopathologic subtypes",
                      test=NULL,
                      vecnoms=c("Nombre aiguilles","Cote","Histo","Grade","Pole",
                                "renal score ","mrenal score","location","renal sinus"
                                ,"volume","taille max","Hyper ou hypovasculaire
","T","dissection","marges","clavien"),valeurs=NULL,
                      vecrefs=NULL,varassoc=NULL,
                      codassoc=NULL,pres=NULL,langue="en",digits=2)


tumeurs_cryo_tt<-TABKRIS(baz=cryo,vect.var = c("nbr_af","dissection","margesc","clavien2","compli"),
                      vect.quali = c(1,1,1,1,1),
                      varint=NULL,valvarint =NULL,
                      nomvarint = "Histopathologic subtypes",
                      test=NULL,
                      vecnoms=c("nombre aiguilles","dissection","marges","clavien","complications"),valeurs=NULL,
                      vecrefs=NULL,varassoc=NULL,
                      codassoc=NULL,pres=NULL,langue="en",digits=2)

