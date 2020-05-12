createCLcdis<-function(CLc,dbedis,taxon0=""){
		CLcdis1<-CLcdis2<-CLc
		CLcdis1@cl<-dbedis@totalW$estim%>%transmute(landCtry="FRA", vslFlgCtry="FRA",
							time,space,technical, taxon=taxon0,
							landCat="HUC",commCatScl="EU",
							commCat="NA",unallocCatchWt=0,
							misRepCatchWt=0,landWt=value,landMult=1,
							landValue=0)
		return(CLcdis1)
}

createCScsex<-function(CSc){
	CSc@sl$spp<-paste(CSc@sl$spp,CSc@sl$sex)
	CSc@hl$spp<-paste(CSc@hl$spp,CSc@hl$sex)
	if(nrow(CSc@ca)>1){
	CSc@ca$spp<-paste(CSc@ca$spp,CSc@ca$sex)
	}
	return(CSc)
}

createdbe<-function(info,myStr){
	#define dbe object
        dbelan<-dbeObject(species=info$species,taxon=info$species,catchCat='LAN',
			  strataDesc=myStr)
	wdbelan <-dbeObject(species=info$species,taxon=info$species,catchCat='LAN',
			    strataDesc=myStr,methodDesc='analytical',
			    param='weight',adjust=F)
	ldbelan <-dbeObject(species=info$species,taxon=info$species,catchCat='LAN',
			    strataDesc=myStr,methodDesc='analytical',
			    param='length',adjust=F)
        dbedis<-dbeObject(species=info$species,taxon=info$species,catchCat='DIS',
			  strataDesc=myStr,methodDesc="analytical")
	wdbedis<-dbeObject(species=info$species,taxon=info$species,catchCat='DIS',
			    strataDesc=myStr,methodDesc='analytical',
			    param='weight',adjust=F)
	ldbedis<-dbeObject(species=info$species,taxon=info$species,catchCat='DIS',
			    strataDesc=myStr,methodDesc='analytical',
			    param='length',adjust=F)
	return(list(dbelan,wdbelan,ldbelan,dbedis,wdbedis,ldbedis))
}


dofinaltab<-function(SIfinal,SDfinal,SDagefinal,allverif,param){

#,results="asis"}
recheck1<-SIfinal%>%select(Year,Season,FishingArea,Fleet,CatchCategory)%>%
	distinct()%>%mutate(ices1="catch")
recheck2<-SDfinal%>%select(Year,Season,FishingArea,Fleet,CatchCategory)%>%
	distinct()%>%mutate(ices2="len")
recheck3<-SDagefinal%>%select(Year,Season,FishingArea,Fleet,CatchCategory)%>%
	distinct()%>%mutate(ices3="age")
recheck<-recheck1%>%full_join(recheck2)%>%full_join(recheck3)
recheck[is.na(recheck)]<-"."
recheck<-recheck%>%unite("ices",ices1,ices2,ices3,sep="/",remove=F)
disrate<-full_join(
	allverif%>%filter(CatchCategory=="L")%>%transmute(Year,Season,FishingArea,Fleet,L=CATON),
	allverif%>%filter(CatchCategory=="D")%>%transmute(Year,Season,FishingArea,Fleet,D=CATON))%>%mutate(disrate=D/(L+D))	

pipo<-allverif%>%
	left_join(recheck)%>%
	arrange(desc(Year),Season,desc(Fleet),FishingArea,CatchCategory,ices)%>%
	left_join(disrate%>%select(-L,-D))%>%
	mutate(Time=paste(Year,Season),
	       Time=cell_spec(Time,
			      bold=ifelse(Year==param$currentyear,T,F),
			      color=ifelse(is.na(ices),"black","blue")),
	       A=gsub("27.","",FishingArea),
	       A=cell_spec(A,
			       bold=ifelse(!is.na(ices),T,F),
			       color=ifelse(is.na(ices),"black","white"),
			       background=ifelse(is.na(ices),"white","blue")),
	       Fleet=cell_spec(Fleet,
			       bold=ifelse(!is.na(ices),T,F),
			       color=ifelse(is.na(ices),"black","white"),
			       background=ifelse(is.na(ices),"white","blue")),
	       C=CatchCategory,
	       C=cell_spec(C,
			       bold=ifelse(!is.na(ices),T,F),
			       color=ifelse(is.na(ices),"black",ifelse(C=="L","white","red")),
			       background=ifelse(is.na(ices),"white","blue")),
	       sn=paste0(s,"/",n),
	       checknbsamp=ifelse(is.na(check2),F,check2),
	       checkdist=ifelse(is.na(checknbsampnbfish),F,checknbsampnbfish),
	       checksoplennona=ifelse(is.na(checksoplen),F,checksoplen),
	       checksopagenona=ifelse(is.na(checksopage),F,checksopage),
	       soplshort=ifelse(is.finite(sopl),round(sopl,2),"."),
	       sopashort=ifelse(is.finite(sopa),round(sopa,2),"."),
	       sopla=paste0(soplshort,"/",sopashort),
	       icesnona=ifelse(is.na(ices),".",ices),
	       w_t=round(CATON/1000,2),
	       w_t=cell_spec(w_t,"latex",color="black",
			     background=ifelse(checkLanDis,"green","white")),
	       OffLan_t=ifelse(OffLandings>=0,round(OffLandings/1000,2),"."),
	       OffLan_t=cell_spec(OffLan_t,"latex",color="black",
			     background=ifelse(checkLanDis,"green","white")),
	       sopl=cell_spec(soplshort,"latex",color="black",
			     background=ifelse(checksoplennona,"green","white")),
	       sopa=cell_spec(sopashort,"latex",color="black",
			     background=ifelse(checksopagenona,"green","white")),
	       s=cell_spec(s,"latex",color="black",
			     background=ifelse(checknbsamp,"green","white")),
	       sn=cell_spec(sn,"latex",color="black",
			     background=ifelse(checkdist,"green","white")),
	       ices=cell_spec(icesnona,"latex",bold=T,
			      color=ifelse(!is.na(ices),"white","black"),
			     background=ifelse(!is.na(ices),"blue","yellow")),
	       disrate=round(disrate,2)
	       ) %>%
		transmute(Time,Fleet,A,C,w=w_t,OffLan=OffLan_t,disrate,s,sn,sopl,sopa,ices)

		#w_t=ifelse(checkLanDis,cell_spec(round(CATON/1000,2),"latex",color="green"),".",),
#		check=paste(checkLanDis,checknbsampnbfish,checksoplen,checksopage))


	pipo[is.na(pipo)]<-"."
	finaltab<-pipo
	return(finaltab)
}


#summarise ices data according to allverifop
getfinaldata<-function(icesdat,SDtmp,SDagetmp,allverif,allverifop,param){
	#filter data
	SIfinal<-left_join(icesdat$SI,allverifop)%>%
		filter(checkLanDis)#%>%
	SIfinal<-SIfinal[,names(icesdat$SI)]%>%select(-n,-s)
	HIfinal<-semi_join(icesdat$HI,SIfinal%>%select(-RecordType))
	SDfinal<-left_join(SDtmp[,names(icesdat$SD)],allverifop)%>%
		filter(checkLanDis & checknbsampnbfish&checksoplen)%>%
		select(names(icesdat$SD))
	SDagefinal<-left_join(SDagetmp[,names(icesdat$SD)],allverifop)%>%
		filter(checkLanDis & checknbsampnbfish&checksopage)%>%
		select(names(icesdat$SDage))

	#summarise ices transmitted info
	nSDfinal<-SDfinal%>%select(Year,Season,Fleet,FishingArea,CatchCategory)%>%distinct()%>%nrow()
	nSDagefinal<-SDagefinal%>%select(Year,Season,Fleet,FishingArea,CatchCategory)%>%distinct()%>%nrow()
	vices<-c(catchCat="ices",nbtot=NA,nb=NA,nbok=nrow(SIfinal),
			  nblen=NA,nblenok=NA,nblensopok=nSDfinal,
			  nbagesopok=nSDagefinal)
	nSDfinalnow<-SDfinal%>%select(Year,Season,Fleet,FishingArea,CatchCategory)%>%
		filter(Year==param$currentyear)%>%distinct()%>%nrow()
	nSDagefinalnow<-SDagefinal%>%select(Year,Season,Fleet,FishingArea,CatchCategory)%>%
		filter(Year==param$currentyear)%>%distinct()%>%nrow()
	vicesnow<-c(catchCat="ices",nbtot=NA,nb=NA,
		    nbok=nrow(SIfinal%>%filter(Year==param$currentyear)),
			  nblen=NA,nblenok=NA,nblensopok=nSDfinalnow,
			  nbagesopok=nSDagefinalnow)
	#overall table on check
	v1<-fctchecklan(allverif,"L")
	v1bms<-fctchecklan(allverif,"B")
	v1regdis<-fctchecklan(allverif,"R")
	v2<-fctcheckdis(allverif)
	v12<-rbind(v1,v2)
	vtot<-c(catchcat="C",apply(v12[,2:8],2,sum))
	vtot<-rbind(vtot,v1bms,v1regdis)#%>%tapply(2,sum)
	v1now<-fctchecklan(allverif%>%filter(Year==param$currentyear))
	v1nowbms<-fctchecklan(allverif%>%filter(Year==param$currentyear),"B")
	v1nowregdis<-fctchecklan(allverif%>%filter(Year==param$currentyear),"R")
	v2now<-fctcheckdis(allverif%>%filter(Year==param$currentyear))
	v12now<-rbind(v1now,v2now)#%>%tapply(2,sum)
	vtotnow<-c(catchcat="C",apply(v12now[,2:8],2,sum))
	vtotnow<-rbind(vtotnow,v1nowbms,v1nowregdis)#%>%tapply(2,sum)
	tabverif<-rbind(v1now,v2now,vtotnow,vicesnow,v1,v2,vtot,vices)
	pipo<-data.frame(time=c(rep(param$currentyear,6),rep(paste(param$currentyear,param$currentyear-param$nbyear,sep="-"),6)))
	tabverif<-cbind(pipo,tabverif)

	return(list(HIfinal,SIfinal,SDfinal,SDagefinal,tabverif))

}


	#summarise check
	fctchecklan<-function(tmp,catchtype="L"){
		#nb de strate
		nbtot<-nrow(tmp)
		#subset by catchcat
		vtmp<-tmp[tmp$CatchCategory==catchtype,]
		#nb strata for this catchcat
		nb<-nrow(vtmp)
		#nb strata for this catchcat ok
		nbok<-length(vtmp$check1[vtmp$check1])
		#nb strata with len distrib
		nblen<-length(vtmp$check1[is.finite(vtmp$CATONSD)])
		#nb strata with len distrib with enough samp and fish nb
		nblenok<-length(vtmp$check1[is.finite(vtmp$CATONSD)&vtmp$check4])
		#nb strata with len distrib with enough samp and fish nb and sop len ok
		nblensopok<-length(vtmp$check1[is.finite(vtmp$CATONSD)&vtmp$check4&vtmp$check5])
		#nb strata with len distrib with enough samp and fish nb and sop age ok
		nbagesopok<-length(vtmp$check1[is.finite(vtmp$CATONSDage)&vtmp$check4&vtmp$check6])
		pipo<-data.frame(catchcat=catchtype,nbtot,nb,nbok,nblen,nblenok,nblensopok,nbagesopok)
		return(pipo)
	}
	#summarise check
	fctcheckdis<-function(allverif){
		#nb de strate
		nbtot<-nrow(allverif)
		#subset by catchcat
		vtmp<-allverif[allverif$CatchCategory=="D",]
		#nb strata for this catchcat
		nb<-nrow(vtmp)
		#nb strata for this catchcat ok
		nbok<-length(vtmp$check2[vtmp$check2])
		#nb strata with len distrib
		nblen<-length(vtmp$check2[vtmp$check2&is.finite(vtmp$CATONSD)])
		#nb strata with len distrib with enough samp and fish nb
		nblenok<-length(vtmp$check2[vtmp$check2&is.finite(vtmp$CATONSD)&vtmp$check3])
		#nb strata with len distrib with enough samp and fish nb and sop len ok
		nblensopok<-length(vtmp$check2[vtmp$check2&is.finite(vtmp$CATONSD)&vtmp$check3&vtmp$check5])
		#nb strata with len distrib with enough samp and fish nb and sop age ok
		nbagesopok<-length(vtmp$check2[vtmp$check2&is.finite(vtmp$CATONSDage)&vtmp$check3&vtmp$check6])
		pipo<-data.frame(catchcat="D",nbtot,nb,nbok,nblen,nblenok,nblensopok,nbagesopok)
		return(pipo)
	}


#check all verif results
checkallverif<-function(allverif,param){
	#1. si diff CATON and landWt CLc > 1e-5 -> stop
	allverif<-allverif%>%
		mutate(check1=ifelse(diffCATONcl<1e-5,TRUE,FALSE))
	if(any(!na.omit(allverif$check1))){
		print(allverif%>%filter(!check1))
		stop("Difference in landings estimated and registered !")
	}
	#2. remove strata with low sampling number for discards
	allverif<-allverif%>%
		mutate(check2=ifelse(s>=param$seuilsampledis,TRUE,FALSE))%>%
		mutate(check2=ifelse(CatchCategory=="D",check2,NA))
	#3. select strata with s and n ok for discards length distrib
	allverif<-allverif%>%
		mutate(check3=ifelse(n>=param$seuilfish & s>=param$seuilsampledis,TRUE,FALSE))%>%
		mutate(check3=ifelse(CatchCategory=="D",check3,NA))
	#4. as 3 but for landings
	allverif<-allverif%>%
		mutate(check4=ifelse(n>=param$seuilfish & s>=param$seuilsample,TRUE,FALSE))%>%
		mutate(check4=ifelse(CatchCategory=="L",check4,NA))
	#5. sop for length distrib
	allverif<-allverif%>%
		mutate(check5=ifelse(0.8<=sopl & sopl<=1.2,TRUE,FALSE))%>%
		mutate(check5=ifelse(is.finite(CATONSD),check5,NA))
	#6. sop for age distrib
	allverif<-allverif%>%
		mutate(check6=ifelse(0.8<=sopa & sopa<=1.2,TRUE,FALSE))%>%
		mutate(check6=ifelse(is.finite(CATONSDage),check6,NA))
	#7. harmonize check
	allverif<-allverif%>%
		mutate(checkLanDis=ifelse(CatchCategory=="L",check1,check2))%>%
		mutate(checkLanDis=ifelse(CatchCategory%in%c("B","R"),T,checkLanDis))%>%
		mutate(checknbsampnbfish=ifelse(CatchCategory=="L",check4,check3))%>%
		mutate(checksoplen=check5,checksopage=check6)#%>%
		#select(-check1,-check2,-check3,-check4,-check5,-check6)
	allverifop<-allverif%>%select(Year,Season,Fleet,FishingArea,CatchCategory,
				      checkLanDis:checksopage)
	return(list(allverif,allverifop))
}



#aggregate rtp if needed
rtpagg<-function(rtp,myStr,param){
	#agregate rtp in space if myStr agregates space
	spacefromto<-data.frame(area=myStr@spRec$from,areanew=myStr@spRec$to)
	rtpagg<-left_join(rtp,spacefromto,by="area")%>%
	mutate(area=areanew)%>%select(-areanew)
	rtpagg<-rtpagg%>%group_by(spp,area,quarter,sex)%>%
	summarise(a=median(a),
	b=median(b))%>%ungroup()
	rtp<-rtpagg
	#if year remove season in rtp
	if(param$timestratif=="Year"&exists("rtp")){
		rtp<-rtp%>%mutate(quarter="")%>%
		group_by(spp,area,quarter,sex)%>%
		summarise(a=median(a),
		b=median(b))%>%ungroup()
	}
	return(rtp)
}

#shoren a species name
shortspp<-function(spp){
	spp1<-unlist(strsplit(spp," "))
	spp1<-paste(substr(spp1,1,3),collapse=" ")
	return(spp1)
}


#source("data2icesxls.R")
dislantot<-function(dbedis,CSr,info,param,myStr){
	if(F){
		load("test.rdata")
		library(COSTdbe)
		library(ggplot2)
		library(dplyr)
	}
	#load comp data
	load("./data/CLrcompall.rdata")
	 y3<-c(param$currentyear:(param$currentyear-param$nbyear))
	 #met lan with samp by area              
	 CLrcomp<-subset(CLrcomp,year%in%y3,table="cl")

	load("./data/CSrcompall.rdata")
	#subset CSrcomp to match the quality check of CSr
	CSrcomp<-subCS(CSrcomp,CSr)
	#consolidate and validate
	CSvcomp <- csDataVal(CSrcomp); CSccomp <- csDataCons(CSvcomp,myStr)
	#idem for CLcomp
	CLvcomp <- clDataVal(CLrcomp); CLccomp <- clDataCons(CLvcomp,myStr)
	#put species in taxon to avoid COST bug in CLc
	CLccomptmp<-CLccomp
	CLccomptmp@cl$taxon<-as.character(CLccomptmp@cl$taxon)
	CLccomptmp@cl$taxon[CLccomptmp@cl$taxon==info$taxon]<-info$species

	#remove MIS_MIS in CSc
	CSccomp<-subset(CSccomp,!grepl("MIS_MIS",technical),table="hh")
	if(nrow(CSccomp@ca)>1&info$taxon!="NEP"){
		  CSccomp<-alkLgthRec(CSccomp,type='stepIncr',param$stepIncr,preview=FALSE,postview=FALSE,update=TRUE)
	}
	#define the dbedis object
	dbedis<-dbeObject(species=info$species,taxon=info$species,catchCat='DIS',
			  strataDesc=myStr,methodDesc="analytical")
	#raise the discaaaaards
	listspecies<-c(unique(as.character(CSccomp@sl$spp)),unique(as.character(CLccomp@cl$taxon)))
	listspecies<-listspecies[!grepl("Argentina silus",listspecies)]
	#listspecies<-unique(CLccomp@cl$taxon)
	dbedis<-totVolume(dbedis,CSccomp,CEc,CLccomptmp,type='landings',landSpp=listspecies)
	return(dbedis)

	if(F){
	#explo to understand stuff
	cl1<-CLccomptmp@cl%>%filter(taxon!=info$species)%>%
		mutate(uu="other0")%>%#ifelse(taxon==info$species,"cible0","other0"))%>%
		group_by(time,space,technical,uu)%>%
		summarise(tot=sum(landWt))%>%ungroup()%>%
		tidyr::pivot_wider(values_from=tot,names_from=uu)
	#filter haul
	u1<-sampledFO(CSccomp,info$species,"DIS")[[1]]
	u2<-sampledFO(CSccomp,listspecies,"LAN")[[1]]
	idhaulok<-!is.na(u1)&!is.na(u2)
	idhaulok<-CSccomp@hh%>%filter(idhaulok)%>%select(trpCode,staNum)

	sl1<-left_join(idhaulok,CSccomp@sl)%>%mutate(uu=ifelse(as.character(spp)==info$species,"cible","other"))%>%
		group_by(time,space,technical,PSUid,SSUid,uu,type=substr(sort,1,3))%>%
		summarise(tot=sum(wt))%>%
		group_by(time,space,technical,PSUid)%>%
		mutate(nbhaul=n_distinct(SSUid))%>%ungroup()%>%
		filter((uu=="cible"&type=="DIS")|(uu=="other"&type=="LAN"))
	nbhaultot<-CSccomp@hh%>%group_by(time,space,technical,PSUid)%>%
		summarise(nbhaultot=n_distinct(SSUid))%>%ungroup()
	sl1<-left_join(sl1,nbhaultot)%>%mutate(tot2=tot*nbhaultot/nbhaul)
	sl1<-sl1%>%group_by(time,space,technical,uu)%>%summarise(tot=sum(tot))%>%ungroup()#*nbhaultot/nbhaul)
	sl1<-sl1%>%tidyr::pivot_wider(values_from=tot,names_from=uu)
	pipo<-full_join(sl1,cl1)%>%mutate(test=cible/other*other0)%>%left_join(dbedis@totalW$estim)%>%
		mutate(test2=(test-value)/1000)
	pipo%>%filter(!is.na(value))
	}


}
#subset a CSr according to another CSr on trpCode and staNum
subCS<-function(CSrcomp,CSr){
	idtr<-unique(CSr@tr$trpCode)
	CSrcomp@tr<-CSrcomp@tr%>%filter(trpCode%in%idtr)
	idtr<-paste(CSr@hh$trpCode,CSr@hh$staNum)%>%unique()
	CSrcomp@hh<-CSrcomp@hh%>%mutate(idtrlocal=paste(trpCode,staNum))%>%
					filter(idtrlocal%in%idtr)%>%
					select(-idtrlocal)
	idtr<-paste(CSr@sl$trpCode,CSr@sl$staNum)%>%unique()
	CSrcomp@sl<-CSrcomp@sl%>%mutate(idtrlocal=paste(trpCode,staNum))%>%
					filter(idtrlocal%in%idtr)%>%
					select(-idtrlocal)
	idtr<-paste(CSr@hl$trpCode,CSr@hl$staNum)%>%unique()
	CSrcomp@hl<-CSrcomp@hl%>%mutate(idtrlocal=paste(trpCode,staNum))%>%
					filter(idtrlocal%in%idtr)%>%
					select(-idtrlocal)
	pipo<-csData(CSrcomp@tr,CSrcomp@hh,CSrcomp@sl,CSrcomp@hl,CSr@ca)
	return(pipo)
}

aggcl<-function(CLr,taxonnew){
	if(F){
		load("./data/CLrall.rdata")
		taxon<-"pipo"
		library(dplyr)
	}
	desc<-CLr@desc
	newcl<-CLr@cl%>%group_by(landCtry,vslFlgCtry,year,quarter,month,area,rect,subRect,
			    taxon=taxonnew,landCat,commCatScl,commCat,foCatNat,foCatEu5,foCatEu6,
			    harbour,vslLenCat,unallocCatchWt,misRepCatchWt,landMult)%>%
  		   summarise(landWt=sum(landWt,na.rm=T),landValue=sum(landValue,na.rm=T))%>%
		   ungroup()%>%data.frame()
  	newcl<-newcl[,names(CLr@cl)]
  	CLr<-clData(newcl)
  	CLr@desc<-desc
	return(CLr)
}

aggcs<-function(CSr,spp){
	if(F){
		load("./data/CSrall.rdata")
		spp<-"pipo spp"
		library(dplyr)
	}
	desc<-CSr@desc
	CSr@sl$spp <- CSr@hl$spp <- spp
#hl
#c('sampType', 'landCtry', 'vslFlgCtry', 'year', 'proj', 'trpCode', 'staNum', 'spp', 
#	'catchCat', 'landCat', 'commCatScl', 'commCat', 'subSampCat', 'sex', 'lenCls')
	newhl<-CSr@hl%>%group_by(sampType,landCtry,vslFlgCtry,year,proj,trpCode,staNum,spp,
			catchCat,landCat,commCatScl,commCat,subSampCat,sex,lenCls)%>%
			summarise(lenNum=sum(lenNum,na.rm=T))%>%ungroup()%>%data.frame()
	newhl<-newhl[,names(CSr@hl)]
#sl
#c('sampType', 'landCtry', 'vslFlgCtry', 'year', 'proj', 'trpCode', 'staNum', 'spp', 
#	'catchCat', 'landCat', 'commCatScl', 'commCat', 'subSampCat')
	newsl<-CSr@sl%>%group_by(sampType,landCtry,vslFlgCtry,year,proj,trpCode,staNum,spp,
			catchCat,landCat,commCatScl,commCat,subSampCat,sex)%>%
			summarise(lenCode=paste(unique(lenCode),collapse=","),
				  wt=sum(wt,na.rm=T),subSampWt=sum(subSampWt,na.rm=T))%>%
			ungroup()%>%data.frame()

	newsl<-newsl[,names(CSr@sl)]

	CSr<-csData(CSr@tr,CSr@hh,newsl,newhl,CSr@ca,check=F)
	CSr@desc<-desc

	return(CSr)
}

#recodification metier sole-bisc following muriel spec
metsolmuriel<-function(CLr,CEr,CSr){
	if(F){
		library(dplyr)
		load("./data/CLrall.rdata")
		load("./data/CErall.rdata")
		load("./data/CSrall.rdata")
	}
	#compute length cat for CSr
	classloa<-c("u10","10-12","12-15","15-18","18-24","24-40","o40")
	tr<-CSr@tr%>%select(sampType:vslLen)%>%
			mutate(classloa=classloa[1+findInterval(vslLen,vec=c(10,12,15,18,24,40))])%>%
			select(-vslLen)
	hh<-left_join(CSr@hh,tr,by=c("sampType","landCtry","vslFlgCtry","year","proj","trpCode"))%>%
		mutate(foCatEu6=paste(foCatEu6,classloa,sep=":"))%>%select(-classloa)
	if(nrow(hh)!=nrow(CSr@hh)){stop("pb dim CSr")}
	CSr@hh<-hh
	#new metier list with size
	CLr@cl$foCatEu6<-paste(CLr@cl$foCatEu6,CLr@cl$vslLenCat,sep=":")
	CEr@ce$foCatEu6<-paste(CEr@ce$foCatEu6,CEr@ce$vslLenCat,sep=":")
	#generate metier list
	#Inshore-Gillnets: gillnet<12m
	#Offshore-Gillnets: gillnet>12m
	#Inshore-trawlers: trawler<12m
	#Offshore-trawlers: trawler>12m
	metier<-data.frame(foCatEu6=unique(c(CSr@hh$foCatEu6,CSr@hh$foCatEu6,CEr@ce$foCatEu6)),final="MIS_MIS_0_0_0_HC")
	metier$final[substr(metier$foCatEu6,1,1)%in%c("T","O")&grepl(":10-12",metier$foCatEu6)]<-"Inshore-trawlers"
	metier$final[substr(metier$foCatEu6,1,1)%in%c("T","O")&grepl(":u10",metier$foCatEu6)]<-"Inshore-trawlers"
	metier$final[substr(metier$foCatEu6,1,2)%in%c("PT")&grepl(":10-12",metier$foCatEu6)]<-"Inshore-trawlers"
	metier$final[substr(metier$foCatEu6,1,2)%in%c("PT")&grepl(":u10",metier$foCatEu6)]<-"Inshore-trawlers"
	metier$final[substr(metier$foCatEu6,1,1)%in%c("T","O")&metier$final=="MIS_MIS_0_0_0_HC"]<-"Offshore-trawlers"
	metier$final[substr(metier$foCatEu6,1,2)%in%c("PT")&metier$final=="MIS_MIS_0_0_0_HC"]<-"Offshore-trawlers"
	metier$final[substr(metier$foCatEu6,1,1)%in%c("G")&grepl(":10-12",metier$foCatEu6)]<-"Inshore-Gillnets"
	metier$final[substr(metier$foCatEu6,1,1)%in%c("G")&grepl(":u10",metier$foCatEu6)]<-"Inshore-Gillnets"
	metier$final[substr(metier$foCatEu6,1,1)%in%c("G")&metier$final=="MIS_MIS_0_0_0_HC"]<-"Offshore-Gillnets"

	return(list(CSr=CSr,CLr=CLr,CEr=CEr,metier=metier))
}

metbssmikael<-function(CLr,CEr,CSr){
	if(F){
		library(dplyr)
		load("./data/CLrall.rdata")
		load("./data/CErall.rdata")
		load("./data/CSrall.rdata")
	}
	metier<-data.frame(foCatEu6=unique(c(CSr@hh$foCatEu6,CSr@hh$foCatEu6,CEr@ce$foCatEu6)),final="MIS_MIS_0_0_0_HC")
	#ligne
	metier$final[substr(metier$foCatEu6,1,3)%in%c("LHM","LHP","LTL")]<-"LHM_DEF"
	metier$final[substr(metier$foCatEu6,1,3)%in%c("LLS","LLD","LL_")]<-"LLS_DEF"
	#filet
	metier$final[substr(metier$foCatEu6,1,3)%in%c("GTR","GTN")]<-"GTR_DEF_all_0_0_all"
	metier$final[substr(metier$foCatEu6,1,3)%in%c("GN_","GNS","GND","GNC")]<-"GNS_DEF_all_0_0_all"
	#seine
	metier$final[substr(metier$foCatEu6,1,3)%in%c("SDN","SSC","SDV")]<-"SSC_DEF_All_0_0_All"
	metier$final[substr(metier$foCatEu6,1,3)%in%c("PS_","PS1")]<-"PS_SPF_0_0_0"
	#chalut de fond
	metier$final[substr(metier$foCatEu6,1,3)%in%c("OTB","OT_","OTT","TB_","TBB","TBN","TBS")]<-"OTB_DEF"
	metier$final[substr(metier$foCatEu6,1,7)%in%c("OTB_CRU")]<-"OTB_CRU"
	#chalut pelagique
	metier$final[substr(metier$foCatEu6,1,3)%in%c("OTM")]<-"OTM_DEF"
	metier$final[substr(metier$foCatEu6,1,3)%in%c("PT_","PTB","PTM")]<-"PTM-DEF"
	
	return(list(CSr=CSr,CLr=CLr,CEr=CEr,metier=metier))
	

}


diagmet<-function(CLrtmp,CSrtmp,CErtmp){
	if(F){
		library(dplyr);library(ggplot2)
		load("testmet.rdata")
		#bidouille
		areafish<-dorawmapmet(CLrtmp,CSrtmp,CErtmp)
		#extract mesh size
		areafish<-areafish%>%#mutate(met6=gsub("","",foCatEu6))
			tidyr::separate(foCatEu6,c("gear","spp","min","max","sel"),sep="_",remove=FALSE)%>%
			mutate(max=ifelse(grepl(">=",min),gsub(">=","",min),max),
			       min=ifelse(grepl(">=",min),gsub(">=","",min),min),
			       min=as.numeric(min),max=as.numeric(max),
			       autoices="MIS_MIS_0_0_0"
			       )
		pipo<-areafish%>%select(foCatEu6,gear,spp,min,max,sel,engin,autoices)%>%distinct()
		#load("test.rdata")
		
	}
	#bidouille
	areafish<-dorawmapmet(CLrtmp,CSrtmp,CErtmp)
	#extract mesh size
	areafish<-areafish%>%#mutate(met6=gsub("","",foCatEu6))
		tidyr::separate(foCatEu6,c("gear","spp","min","max","sel"),sep="_",remove=FALSE)%>%
		mutate(max=ifelse(grepl(">=",min),gsub(">=","",min),max),
		       min=ifelse(grepl(">=",min),gsub(">=","",min),min),
		       min=as.numeric(min),max=as.numeric(max),
		       autoices="MIS_MIS_0_0_0"
		       )
	#basic group: refactoring to match intercatch datacall metier more or less
	#
	#MIS, DIV, FOO, FWR, LN, OTH, SB, SPR, REC to MIS_MIS: done intrinsically
	#DRB
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("DRB")&spp%in%c("MOL"),"DRB_MOL_0_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("DRB")&!spp%in%c("MOL"),"DRB_all_0_0_all",autoices))
	#FPO
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("FPO")&spp%in%c("CRU"),"FPO_CRU_0_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("FPO")&spp%in%c("MOL","CEP"),"FPO_MOL_0_0_0_all",autoices))
	#FYK
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("FYK"),"FYC_C",autoices))
	#GNS, GND
	#if GND/GNC then GNS (my guess...)
	areafish<-areafish%>% mutate(gear=ifelse(gear%in%c("GND","GNC","GN_"),"GNS",gear))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("GNS")&spp%in%c("CRU"),
				       "GNS_CRU_0_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("GNS")&!spp%in%c("CRU")&min>=220,
				       "GNS_DEF_>=220_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("GNS")&!spp%in%c("CRU")&min==100&max==100,
				       "GNS_DEF_>=100_0_0",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("GNS")&!spp%in%c("CRU")&120<=min&max<=219,
				       "GNS_DEF_120-219_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("GNS")&!spp%in%c("CRU")&100<=min&max<=119,
				       "GNS_DEF_100-119_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("GNS")&!spp%in%c("CRU")&80<=min&max<=99,
				       "GNS_DEF_80-99_0_0",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("GNS")&!spp%in%c("CRU")&60<=min&max<=79,
				       "GNS_DEF_60-79_0_0",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("GNS")&!spp%in%c("CRU")&10<=min&max<=30,
				       "GNS_DEF_10-30_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("GNS")&!spp%in%c("CRU")&!grepl("GN",autoices),
				       "GNS_DEF_all_0_0_all",autoices))
	#GTN
	#if GTN then GTR (my guess...)
	areafish<-areafish%>% mutate(gear=ifelse(gear%in%c("GTN"),"GTR",gear))
	#GTR
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("GTR")&spp%in%c("CRU"),
				       "GTR_CRU_0_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("GTR")&!spp%in%c("CRU")&min>=220,
				       "GTR_DEF_>=220_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("GTR")&!spp%in%c("CRU")&120<=min&max<220,
				       "GTR_DEF_120-219_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("GTR")&!spp%in%c("CRU")&100<=min&max<120,
				       "GTR_DEF_100-119_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("GTR")&!spp%in%c("CRU")&90<=min&max<=99,
				       "GTR_DEF_90-99_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("GTR")&!spp%in%c("CRU")&60<=min&max<=79,
				       "GTR_DEF_60-79_0_0",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("GTR")&!spp%in%c("CRU")&40<=min&max<=59,
				       "GTR_DEF_40-59_0_0",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("GTR")&!spp%in%c("CRU")&!grepl("GTR",autoices),
				       "GTR_DEF_all_0_0_all",autoices))
	#ligne
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("LLS","LLD","LL_","LX"),
				       "LLS_DEF",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("LLS","LLD","LL_")&spp%in%c("DWS"),
				       "LLS_DWS",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("LHP","LHM","LTL"),
				       "LHM_DEF",autoices))
	#OTB+PTB
	areafish<-areafish%>%
		mutate(gear=ifelse(gear%in%c("PTB"),"OTB",gear))
	#OTB_CRU
	#areafish<-areafish%>% mutate(gear=ifelse(gear%in%c("PTB"),"OTB",gear))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("CRU"),
				       "OTB_CRU_All_0_0_All",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("CRU")&16<=min&max<=31,
				       "OTB_CRU_16-31_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("CRU")&32<=min&max<=69,
				       "OTB_CRU_32-69_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("CRU")&70<=min&max<=99,
				       "OTB_CRU_70-99_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("CRU")&70==min&max==70,
				       "OTB_CRU_>=70_0_0",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("CRU")&100<=min&max<=119,
				       "OTB_CRU_100-119_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("CRU")&max>=120,
				       "OTB_CRU_>=120_0_0_all",autoices))
	#OTB_DWS
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("DWS"),
				       "OTB_DWS",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("DWS")&100<=min&max<=119,
				       "OTB_DWS_100-119_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("DWS")&max>=120,
				       "OTB_DWS_>=120_0_0_all",autoices))
	#OTB_DEF
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("DEF","ANA","CAT","CEP"),
				       "OTB_DEF_All_0_0_All",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("DEF","ANA","CAT","CEP")&max<=16,
				       "OTB_DEF_<16_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("DEF","ANA","CAT","CEP")&16<=min&max<=31,
				       "OTB_DEF_16-31_0_0",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("DEF","ANA","CAT","CEP")&32<=min&max<=69,
				       "OTB_DEF_32-69_0_0",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("DEF","ANA","CAT","CEP")&70<=min&max<=99,
				       "OTB_DEF_70-99_0_0",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("DEF","ANA","CAT","CEP")&70==min&max==70,
				       "OTB_DEF_>=70_0_0",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("DEF","ANA","CAT","CEP")&100<=min&max<=119,
				       "OTB_DEF_100-119_0_0",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("DEF","ANA","CAT","CEP")&max>=120,
				       "OTB_DEF_>=120_0_0",autoices))
	#OTB_SPF
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("SPF"),
				       "OTB_DEF_All_0_0_All",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("SPF")&32<=min&max<=69,
				       "OTB_SPF_32-69_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("SPF")&70<=min&max<=99,
				       "OTB_SPF_70-99_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("SPF")&max>=120,
				       "OTB_SPF_>=120_0_0_all",autoices))
	#OTB_MOL
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("MOL"),
				       "OTB_DEF_All_0_0_All",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("MOL")&32<=min&max<=69,
				       "OTB_MOL_32-69_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("MOL")&70<=min&max<=99,
				       "OTB_MOL_70-99_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("MOL")&100<=min&max<=119,
				       "OTB_MOL_100-119_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("MOL")&max>=120,
				       "OTB_MOL_>=120_0_0_all",autoices))
	#OTM and PTM
	areafish<-areafish%>%
		mutate(gear=ifelse(gear%in%c("PTM"),"OTM",gear))
	#OTM_DEF
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTM")&spp%in%c("DEF","ANA","CAT","CRU","DES","FIF","LPF","MIS","MOL"),
				       "OTM_DEF",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTM")&spp%in%c("DEF","ANA","CAT","CRU","DES","FIF","LPF","MIS","MOL")&32<=min&max<=69,
				       "OTM_DEF_32-69_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTM")&spp%in%c("DEF","ANA","CAT","CRU","DES","FIF","LPF","MIS","MOL")&70<=min&max<=99,
				       "OTM_DEF_70-99_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTM")&spp%in%c("DEF","ANA","CAT","CRU","DES","FIF","LPF","MIS","MOL")&100<=min&max<=119,
				       "OTM_DEF_100-119_0_0_all",autoices))
	#OTM_SPF
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTM")&spp%in%c("SPF"),
				       "OTM_DEF",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTM")&spp%in%c("SPF")&16<=min&max<=31,
				       "OTM_SPF_16-31_0_0",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTM")&spp%in%c("SPF")&32<=min&max<=69,
				       "OTM_SPF_32-69_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTM")&spp%in%c("SPF")&70<=min&max<=99,
				       "OTM_SPF_70-99_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTB")&spp%in%c("SPF")&max>=120,
				       "OTB_SPF_>=120_0_0_all",autoices))
	#OTM_DWS
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTM")&spp%in%c("DWS"),
				       "OTM_DWS",autoices))
	#OTM le reste
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTM")&grepl("MIS_MIS",autoices),
				       "OTM_DEF",autoices))
	#OTT_DEF
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTT")&spp%in%c("DEF","ANA","CAT","CEP","MOL"),
				       "OTT-DEF",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTT")&spp%in%c("DEF","ANA","CAT","CEP","MOL")&16<=min&max<=31,
				       "OTT_DEF_16-31_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTT")&spp%in%c("DEF","ANA","CAT","CEP","MOL")&70<=min&max<=99,
				       "OTT_DEF_70-99_0_0",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTT")&spp%in%c("DEF","ANA","CAT","CEP","MOL")&70==min&max==70,
				       "OTT_DEF_>=70_0_0",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTT")&spp%in%c("DEF","ANA","CAT","CEP","MOL")&100<=min&max<=119,
				       "OTT_DEF_100-119_0_0",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTT")&spp%in%c("DEF","ANA","CAT","CEP","MOL")&max>=120,
				       "OTT_DEF_>=120_0_0_all",autoices))
	#OTT_CRU
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTT")&spp%in%c("CRU"),
				       "OTT-CRU",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTT")&spp%in%c("CRU")&16<=max,
				       "OTT_CRU_<16_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTT")&spp%in%c("CRU")&32<=min&max<=69,
				       "OTT_CRU_32-69_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTT")&spp%in%c("CRU")&70<=min&max<=99,
				       "OTT_CRU_70-99_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTT")&spp%in%c("CRU")&70==min&max==70,
				       "OTT_CRU_>=70_0_0",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTT")&spp%in%c("CRU")&100<=min&max<=119,
				       "OTT_CRU_100-119_0_0",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTT")&spp%in%c("CRU")&max>=120,
				       "OTT_CRU_>=120_0_0_all",autoices))
	#OTT_DWS
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTT")&spp%in%c("DWS"),
				       "OTT-DWS",autoices))
	#OTT le reste
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("OTT")&grepl("MIS_MIS",autoices),
				       "OTT-DEF",autoices))
	#seine
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("SDN","PS","PS1","SV","SB","SSC")&spp%in%c("SPF"),
				       "PS_SPF_0_0_0",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("SDN","PS","PS1","SV","SB","SSC")&70<=min&max<=99,
				       "SSC_DEF_70-99_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("SDN","PS","PS1","SV","SB","SSC")&100<=min&max<=119,
				       "SSC_DEF_100-119_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("SDN","PS","PS1","SV","SB","SSC")&120<=max,
				       "SSC_DEF_>=120_0_0_all",autoices))
	#seine le reste
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("SDN","PS","PS1","SV","SB","SSC")&grepl("MIS_MIS",autoices),
				       "SSC_DEF_All_0_0_All",autoices))
	#TBB/TBS
	areafish<-areafish%>%
		mutate(gear=ifelse(gear%in%c("TBS"),"TBB",gear))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("TBB"),"TBB_DEF_all_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("TBB")&70<=min&max<=99,"TBB_DEF_70-99_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("TBB")&100<=min&max<=119,"TBB_DEF_100-119_0_0_all",autoices))
	areafish<-areafish%>%
		mutate(autoices=ifelse(gear%in%c("TBB")&120<=max,"TBB_DEF_>=120_0_0_all",autoices))


	#####################################################
	areafish<-areafish%>%select(-nb)%>%
		mutate(catchCat=ifelse(is.na(catchCat),"LAN",catchCat),nb=ifelse(is.na(nbsamp),0,nbsamp))%>%select(-nbsamp)
	areafish<-tidyr::spread(areafish,catchCat,nb,fill=0,drop=T)
	if(!any(names(areafish)%in%"DIS")){
		areafish<-areafish%>%mutate(DIS=0)
	}
	if(!any(names(areafish)%in%"LAN")){
		areafish<-areafish%>%mutate(LAN=0)
	}
	areafish<-areafish%>%transmute(area,time,foCatEu6,autoices,wlan=lan,nlan=LAN,ndis=DIS)
	

	#calcul de truc annuel
	automet<-areafish%>%mutate(year=substr(time,1,4))%>%group_by(year,autoices)%>%
		summarise(wlan=sum(wlan,na.rm=T),ndis=sum(ndis,na.rm=T),nlan=sum(nlan,na.rm=T))%>%ungroup()%>%
		arrange(desc(year),desc(wlan))%>%ungroup()%>%select(autoices,year,wlan,nlan,ndis)%>%distinct()
	automet<-automet%>%group_by(year)%>%arrange(year,wlan)%>%mutate(propwlan=cumsum(wlan)/sum(wlan))%>%ungroup()

	#a graph by autoices metier
  	areafish0<-areafish%>%mutate(tps=as.numeric(substr(time,1,4))+as.numeric(substr(time,6,6))/5,
					area=as.factor(area))%>%
		group_by(autoices,tps,area)%>%
		summarise(wlan=sum(wlan,na.rm=T),nlan=sum(nlan,na.rm=T),ndis=sum(ndis,na.rm=T))%>%
		ungroup()
	mainmet<-automet%>%filter(propwlan>.1)%>%transmute(autoices,main=1)%>%distinct()
	areafish0<-left_join(areafish0,mainmet,by="autoices")%>%mutate(main=ifelse(is.na(main),0,1))
	labelmet<-areafish0%>%select(autoices,main)%>%distinct()%>%arrange(autoices)
        idli<-areafish0%>%transmute(autoices,engin=substr(autoices,1,3))%>%
		distinct()%>%arrange(desc(autoices))%>%group_by(engin)%>%summarise(n=n())
        idli<-c(cumsum(idli$n)+0.5)

        pltrawarea<-ggplot(areafish0)+
                geom_raster(data=areafish0,aes(x=tps,y=autoices,fill=wlan),stat="identity",alpha=.75)+ 
                geom_point(data=areafish0%>%filter(nlan>0),aes(x=tps,y=autoices,size=nlan),shape=1,alpha=.3,color="blue")+ 
                geom_point(data=areafish0%>%filter(ndis>0),aes(x=tps,y=autoices,size=ndis),shape=1,alpha=.3,color="red")+ 
                scale_radius(breaks=c(1,3,10,50))+
                scale_fill_distiller(palette='Spectral',
                        name="Landing (t)",trans="log10")+
                geom_hline(yintercept=idli,color="dark grey",linetype="dashed")+
                facet_grid(~area,drop=T)+ylab("")+
		ylab("")+
                theme(legend.position="bottom",
                        axis.text.x  = element_text(angle=90, vjust=0.5, size=5),
                        axis.text.y  = element_text(angle=0, vjust=0, size=5,
					face=ifelse(labelmet$main==1,"bold","plain")),
                        strip.text.x = element_text(size=8, angle=0))
	return(list(automet=automet,pltrawarea=pltrawarea,areafish=areafish))

	

}


#a function to clean cache before a new run
clean<-function(){
	system("rm -r ./credo_cache")
	system("rm -r ./credo_files")
	system("rm *.rdata")
	system("rm *.Rdata")
	system("rm *.r")
}


#discards all methods
disallmeth<-function(listyear=2016:2014,info,CSc,CEc,CLc){
	test<-F
	if(test){
		listyear<-2016:2014
		load('info.rdata')
		load("datavalcons.Rdata")
		load("myStr.Rdata")
		library(COSTdbe)
		library(dplyr)
		library(ggplot2)

	}
	listyear<-2014:2016
	rez1<-data.frame()
	rez2<-data.frame()
	for(i in 2014:2016){
		print(i)
		rez<-compdisall(CSc,CEc,CLc,info,myStr,year=i)
		#plot(rez[[2]])
		#rezall<-rbind(rezall,rez[[1]])
		rez1<-rbind(rez1,rez[[1]])
		rez2<-rbind(rez1,data.frame(rez[[3]]))
	}

}

compdisall<-function(CSc,CEc,CLc,info,myStr,year=2016){
	if(F){
		library(dplyr)
		library(COSTdbe)
		load("datavalcons.Rdata")
		load("info.rdata")
		load("myStr.Rdata")
		year<-2004

	}
	CLc@cl$taxon<-info$species
	#annual extraction
	CSc<-subset(CSc,grepl(year,time),table="hh")
	CEc<-subset(CEc,grepl(year,time),table="ce")
	CLc<-subset(CLc,grepl(year,time),table="cl")
	#some checks to avoid computational errors
	aa<-realthing(CLc,CSc)
	bb<-aa%>%filter(!is.na(nbsamplanraw)|!is.na(nbsampdisraw))
	CSc<-subset(CSc,time%in%unique(bb$time),table="hh")
	CEc<-subset(CEc,time%in%unique(bb$time),table="ce")
	CLc<-subset(CLc,time%in%unique(bb$time),table="cl")
	CSc<-subset(CSc,space%in%unique(bb$space),table="hh")
	CEc<-subset(CEc,space%in%unique(bb$space),table="ce")
	CLc<-subset(CLc,space%in%unique(bb$space),table="cl")
	CSc<-subset(CSc,technical%in%unique(bb$technical),table="hh")
	CEc<-subset(CEc,technical%in%unique(bb$technical),table="ce")
	CLc<-subset(CLc,technical%in%unique(bb$technical),table="cl")
	

	if(file.exists(paste0("disdbeall",year,".rdata"))){
		#print("loading old file")
		load(paste0("disdbeall",year,".rdata"))
	}else{
		#print("Computing discards")
		#dis
		dbedis <- dbeObject(species=info$species,taxon=info$species,
				  catchCat='DIS',strataDesc=myStr,methodDesc='analytical')
		if(any(grepl("DIS",CSc@hl$sort))){
			dbedistrip<-totVolume(dbedis,CSc,CEc,type='trip',val='weight')
			dbedisfo<-totVolume(dbedis,CSc,CEc,type='fo',val='weight')
			dbedisfd<-totVolume(dbedis,CSc,CEc,type='fd',val='weight')
			dbedislan<-dbedis
			try( dbedislan<-totVolume(dbedis,CSc,CEc,CLc,type='landings',val='weight') ,silent=T)
		}else{
			dbedistrip<-dbedisfo<-dbedisfd<-dbedislan<-dbedis
		}
		#lan
		#print("Computing landings")
		dbelan <- dbeObject(species=info$species,taxon=info$species,
				  catchCat='LAN',strataDesc=myStr,methodDesc='analytical')
		if(any(grepl("LAN",CSc@hl$sort))){
			dbelantrip<-totVolume(dbelan,CSc,CEc,type='trip',val='weight')
			dbelanfo<-totVolume(dbelan,CSc,CEc,type='fo',val='weight')
			dbelanfd<-totVolume(dbelan,CSc,CEc,type='fd',val='weight')
			dbelanlan<-totVolume(dbelan,CSc,CEc,CLc,type='landings',val='weight')
		}else{
			dbelantrip<-dbelanfo<-dbelanfd<-dbelanlan<-dbelan
		}
		save(dbedistrip,dbedisfo,dbedisfd,dbedislan,
			dbelantrip,dbelanfo,dbelanfd,dbelanlan,file=paste0("disdbeall",year,".rdata"))
	}

	#compute real landings
	realinfo<-realthing(CLc,CSc)
	l1<-extrinfo(dbelantrip,nom="lan",meth="trip")
	l2<-extrinfo(dbelanfo,nom="lan",meth="fo")
	l3<-extrinfo(dbelanfd,nom="lan",meth="fd")
	l4<-extrinfo(dbelanlan,nom="lan",meth="lan")
	lestim<-rbind(l1,l2,l3,l4)
	d1<-extrinfo(dbedistrip,nom="dis",meth="trip")
	d2<-extrinfo(dbedisfo,nom="dis",meth="fo")
	d3<-extrinfo(dbedisfd,nom="dis",meth="fd")
	d4<-extrinfo(dbedislan,nom="dis",meth="lan")
	destim<-rbind(d1,d2,d3,d4)

	#allinfo<-full_join(full_join(realinfo,lestim))#,destim)
	allinfo<-left_join(lestim,realinfo)%>%arrange(time,space,technical,meth)#,destim)
	alldis<-left_join(destim,realinfo)
	
	
	tmp<-allinfo%>%ungroup()%>%group_by(time,space,technical)%>%mutate(test=sum(lan_val))%>%ungroup()%>%
			filter((!is.na(lan_val)&!is.na(lan_obs))&(test>0)) %>% 
			mutate(ts=paste(gsub("27.","",space),"/",substr(time,8,8)),met=sub("_","\n",sub("_","-",gsub("_0","",gsub("_all","",technical)))))

	p1<-ggplot(tmp,aes(x=meth,y=lan_val,ymin=lan_inf,ymax=lan_sup,color=meth,group=meth))+
		geom_errorbar(width=.2,size=1,alpha=.7)+geom_point()+
		geom_line(data=tmp,aes(x=meth,y=lan_obs,group="obs"),color="blue",alpha=.7)+#,linetype="dotted")+
		scale_y_continuous(trans="log10")+
		#geom_point(data=tmp,aes(x=meth,y=lan_obs,group="obs"),color="black")+
		xlab("Methods")+ylab("Landings")+
		ggtitle(paste0("Landings vs landings prediction based on auxiliary parameters ",year))+
		theme(axis.text.x = element_text(size=8, angle=90),
		      strip.text.x=element_text(size=8,angle=0),
		      legend.position="bottom")+
	facet_grid(met~ts,scale="free_y")

	#summary of estimation performance
	rez<-tmp%>%ungroup()%>%filter(meth!="lan")%>%mutate(testci=ifelse(lan_inf<=lan_obs & lan_obs<=lan_sup,1,0),
			diff=abs(lan_val-lan_obs)/lan_obs)
	rez<-rez%>%group_by(technical)%>%mutate(tottest=n())%>%group_by(technical,meth)%>%
		summarise(perftot=sum(testci)/unique(tottest),diff=median(diff))%>%ungroup()
	if(nrow(rez)>0){ rez$year<-year }

	return(list(rez,p1,tmp,alldis))

}

compdisall2<-function(CSc,CEc,CLc,info,myStr,year=2016,div="27.8.a",tech="",allspecies){
	if(F){
		library(dplyr)
		library(COSTdbe)
		load("datavalcons.Rdata")
		load("info.rdata")
		load("myStr.Rdata")
		year<-2004

	}
	CLc@cl$taxon<-info$species
	#annual extraction
	CSc<-subset(CSc,grepl(year,time),table="hh")
	CEc<-subset(CEc,grepl(year,time),table="ce")
	CLc<-subset(CLc,grepl(year,time),table="cl")
	#spatial extraction
	CSc<-subset(CSc,grepl(div,space),table="hh")
	CEc<-subset(CEc,grepl(div,space),table="ce")
	CLc<-subset(CLc,grepl(div,space),table="cl")
	#spatial extraction
	CSc<-subset(CSc,grepl(tech,technical),table="hh")
	CEc<-subset(CEc,grepl(tech,technical),table="ce")
	CLc<-subset(CLc,grepl(tech,technical),table="cl")
	#some checks to avoid computational errors
	aa<-realthing2(CLc,CSc,info$species)
	bb<-aa%>%filter(!is.na(nbsamplanraw)|!is.na(nbsampdisraw))
	CSc<-subset(CSc,time%in%unique(bb$time),table="hh")
	CEc<-subset(CEc,time%in%unique(bb$time),table="ce")
	CLc<-subset(CLc,time%in%unique(bb$time),table="cl")
	CSc<-subset(CSc,space%in%unique(bb$space),table="hh")
	CEc<-subset(CEc,space%in%unique(bb$space),table="ce")
	CLc<-subset(CLc,space%in%unique(bb$space),table="cl")
	CSc<-subset(CSc,technical%in%unique(bb$technical),table="hh")
	CEc<-subset(CEc,technical%in%unique(bb$technical),table="ce")
	CLc<-subset(CLc,technical%in%unique(bb$technical),table="cl")
	

	if(file.exists(paste0("disdbeall",info$species,year,div,tech,".rdata"))){
		#print("loading old file")
		load(paste0("disdbeall",info$species,year,div,tech,".rdata"))
	}else{
		#print("Computing discards")
		#dis
		dbedis <- dbeObject(species=info$species,taxon=info$species,
				  catchCat='DIS',strataDesc=myStr,methodDesc='analytical')
		dbedistrip<-dbedisfo<-dbedisfd<-dbedislan<-dbedis
		try(dbedistrip<-totVolume(dbedis,CSc,CEc,type='trip',val='weight'),silent=T)
		try(dbedisfo<-totVolume(dbedis,CSc,CEc,type='fo',val='weight'),silent=T)
		try(dbedisfd<-totVolume(dbedis,CSc,CEc,type='fd',val='weight'),silent=T)
		try(dbedislan<-totVolume(dbedis,CSc,CEc,CLc,type='landings',val='weight') ,silent=T)
		#lan
		#print("Computing landings")
		dbelan <- dbeObject(species=info$species,taxon=info$species,
				  catchCat='LAN',strataDesc=myStr,methodDesc='analytical')
		dbelantrip<-dbelanfo<-dbelanfd<-dbelanlan<-dbelan
		try(dbelantrip<-totVolume(dbelan,CSc,CEc,type='trip',val='weight'),silent=T)
		try(dbelanfo<-totVolume(dbelan,CSc,CEc,type='fo',val='weight'),silent=T)
		try(dbelanfd<-totVolume(dbelan,CSc,CEc,type='fd',val='weight'),silent=T)
		try(dbelanlan<-totVolume(dbelan,CSc,CEc,CLc,type='landings',val='weight'),silent=T)
		save(dbedistrip,dbedisfo,dbedisfd,dbedislan,
			dbelantrip,dbelanfo,dbelanfd,dbelanlan,file=paste0("disdbeall",info$species,year,div,tech,".rdata"))
	}

	#compute real landings
	realinfo<-realthing2(CLc,CSc,info$species)
	l1<-extrinfo(dbelantrip,nom="lan",meth="trip")
	l2<-extrinfo(dbelanfo,nom="lan",meth="fo")
	l3<-extrinfo(dbelanfd,nom="lan",meth="fd")
	l4<-extrinfo(dbelanlan,nom="lan",meth="lan")
	lestim<-rbind(l1,l2,l3,l4)
	d1<-extrinfo(dbedistrip,nom="dis",meth="trip")
	d2<-extrinfo(dbedisfo,nom="dis",meth="fo")
	d3<-extrinfo(dbedisfd,nom="dis",meth="fd")
	d4<-extrinfo(dbedislan,nom="dis",meth="lan")
	destim<-rbind(d1,d2,d3,d4)

	#allinfo<-full_join(full_join(realinfo,lestim))#,destim)
	allinfo<-left_join(lestim,realinfo)%>%arrange(time,space,technical,meth)#,destim)
	alldis<-left_join(destim,realinfo)
	
	
	tmp<-allinfo%>%ungroup()%>%group_by(time,space,technical)%>%mutate(test=sum(lan_val))%>%ungroup()%>%
			filter((!is.na(lan_val)&!is.na(lan_obs))&(test>0)) %>% 
			mutate(ts=paste(gsub("27.","",space),"/",substr(time,8,8)),met=sub("_","\n",sub("_","-",gsub("_0","",gsub("_all","",technical)))))

	p1<-ggplot(tmp,aes(x=meth,y=lan_val,ymin=lan_inf,ymax=lan_sup,color=meth,group=meth))+
		geom_errorbar(width=.2,size=1,alpha=.7)+geom_point()+
		geom_line(data=tmp,aes(x=meth,y=lan_obs,group="obs"),color="blue",alpha=.7)+#,linetype="dotted")+
		scale_y_continuous(trans="log10")+
		#geom_point(data=tmp,aes(x=meth,y=lan_obs,group="obs"),color="black")+
		xlab("Methods")+ylab("Landings")+
		ggtitle(paste0("Landings vs landings prediction based on auxiliary parameters ",year))+
		theme(axis.text.x = element_text(size=8, angle=90),
		      strip.text.x=element_text(size=8,angle=0),
		      legend.position="bottom")+
	facet_grid(met~ts,scale="free_y")

	#summary of estimation performance
	rez<-tmp%>%ungroup()%>%filter(meth!="lan")%>%mutate(testci=ifelse(lan_inf<=lan_obs & lan_obs<=lan_sup,1,0),
			diff=abs(lan_val-lan_obs)/lan_obs)
	rez<-rez%>%group_by(technical)%>%mutate(tottest=n())%>%group_by(technical,meth)%>%
		summarise(perftot=sum(testci)/unique(tottest),diff=median(diff))%>%ungroup()
	if(nrow(rez)>0){ rez$year<-year }

	return(list(rez,p1,tmp,alldis))

}



#extr info from dbe object
extrinfo<-function(dbe=dbedislan,nom="lan",meth="trip"){
  dbew<-dbe@totalWnum$ci
  eval(parse(text= paste0("dbew<-dbew%>%transmute(time,space,technical,",nom,"_val=value,",nom,"_inf=inf,",nom,"_sup=sup)")))
  eval(parse(text= paste0("dbes<-dbe@nSamp$len%>%transmute(time,space,technical,nsamp_",nom,"=value)")))
  eval(parse(text= paste0("dben<-dbe@nMeas$len%>%transmute(time,space,technical,nmeas_",nom,"=value)")))
  eval(parse(text= paste0("dbevar<-dbe@totalWvar%>%transmute(time,space,technical,var_",nom,"=value)")))
  pipo<-full_join(full_join(full_join(dbew,dbevar),dbes),dben)%>%mutate(meth=meth)
  return(pipo)
}

#extract global info from CLc and CSc
realthing<-function(CLc,CSc){
	lanreal<-CLc@cl%>%group_by(time,space,technical)%>%summarize(lan_obs=sum(landWt,na.rm=T))
	nbfishlan<-CSc@hl%>%filter(grepl("LAN",sort))%>%group_by(time,space,technical)%>%summarise(nblanraw=sum(lenNum))
	nbfishdis<-CSc@hl%>%filter(grepl("DIS",sort))%>%group_by(time,space,technical)%>%summarise(nbdisraw=sum(lenNum))
	nbsamp<-full_join(nbfishlan,nbfishdis)
	nbsamplan<-CSc@sl%>%filter(grepl("LAN",sort))%>%select(time,space,technical,trpCode,staNum)%>%
		distinct()%>%group_by(time,space,technical)%>%summarise(nbsamplanraw=n())
	nbsampdis<-CSc@sl%>%filter(grepl("DIS",sort))%>%select(time,space,technical,trpCode,staNum)%>%
		distinct()%>%group_by(time,space,technical)%>%summarise(nbsampdisraw=n())
	nbsamp<-full_join(nbsamp,full_join(nbsamplan,nbsampdis))
	nbsamp<-full_join(lanreal,nbsamp)
	return(nbsamp)
}

#extract global info from CLc and CSc
realthing2<-function(CLc,CSc,species){
	lanreal<-CLc@cl%>%group_by(time,space,technical)%>%summarize(lan_obs=sum(landWt,na.rm=T))
	nbfishlan<-CSc@hl%>%filter(grepl("LAN",sort)&spp==species)%>%group_by(time,space,technical)%>%summarise(nblanraw=sum(lenNum))
	nbfishdis<-CSc@hl%>%filter(grepl("DIS",sort)&spp==species)%>%group_by(time,space,technical)%>%summarise(nbdisraw=sum(lenNum))
	nbsamp<-full_join(nbfishlan,nbfishdis)
	nbsamplan<-CSc@sl%>%filter(grepl("LAN",sort)&spp==species)%>%select(time,space,technical,trpCode,staNum)%>%
		distinct()%>%group_by(time,space,technical)%>%summarise(nbsamplanraw=n())
	nbsampdis<-CSc@sl%>%filter(grepl("DIS",sort)&spp==species)%>%select(time,space,technical,trpCode,staNum)%>%
		distinct()%>%group_by(time,space,technical)%>%summarise(nbsampdisraw=n())
	nbsamp<-full_join(nbsamp,full_join(nbsamplan,nbsampdis))
	nbsamp<-full_join(lanreal,nbsamp)
	return(nbsamp)
}




#dis in time and landings
disexplo<-function(CSc){

restimew<-disCorrPlot(object=CSc,timeStrata=FALSE,techStrata=TRUE, aux='time',
spaceStrata=FALSE,sampPar=TRUE,reg=TRUE,show.legend=TRUE,val='weight',plot=FALSE)

restimew$meth<-"time_weight"
reslanw<-disCorrPlot(object=CSc,timeStrata=FALSE,techStrata=TRUE, aux='landings',
spaceStrata=F,sampPar=TRUE,reg=TRUE,show.legend=TRUE,val='weight',plot=F)
reslanw$meth<-"lan_weight"
p1<-ggplot(reslanw,aes(x=auxVar/1e3+1,y=disVol/1e3+1,shape=space,color=time))+
	geom_point()+facet_wrap(~technical,scales="free")+
	#geom_smooth(mapping=aes(group=1),method="loess",se=T,colour="black",alpha=.5)+
	theme(strip.text.x = element_text(size=6, angle=0))+
	ylab("Discards (kg)")+xlab("Landings (kg)")+
	scale_y_log10()+scale_x_log10()+ggtitle("Discards~landings (log +1)")
p2<-ggplot(restimew,aes(x=auxVar+1,y=disVol/1e3+1,shape=space,color=time))+
	geom_point()+facet_wrap(~technical,scales="free")+
	#geom_smooth(mapping=aes(group=1),method="loess",se=T,colour="black")+
	theme(strip.text.x = element_text(size=6, angle=0))+
	ylab("Discards (kg)")+xlab("Time (mn)")+
	scale_y_log10()+scale_x_log10()+ggtitle("Discards~time (log +1)")
return(list(p1,p2))
}
disvslan2<-function(CSc,info){
	if(F){
		library(dplyr);library(COSTcore)
		library(COSTdbe)
		library(ggplot2)
		load("test.rdata")

	}
	#myStr@timeStrata<-"year"

	#print("length")
	dbelan<-dbeObject(species=info$species,taxon=info$species,catchCat='LAN',strataDesc=myStr)
	dbelan<-RaiseLgth(dbelan,CSc,spp=info$species,taxon=info$species)
	dbedis<-dbeObject(species=info$species,taxon=info$species,catchCat='DIS',strataDesc=myStr)
	dbedis<-RaiseLgth(dbedis,CSc,spp=info$species,taxon=info$species)
	disdf<-rbind(data.frame(dbedis@lenStruc$estim,fraction='DIS'),
		  data.frame(dbelan@lenStruc$estim,fraction='LAN'))%>%
		mutate(time=gsub("20","",gsub(" ","",substr(time,1,4))))%>%
		mutate(technical=gsub("_0_0","",gsub("_all","",technical)))%>%
		mutate(length=as.numeric(as.character(length)))%>%
		mutate(value=as.numeric(as.character(value)))
	disdf<-disdf%>%group_by(time,space,technical,fraction,length)%>%summarize(value=sum(value,na.rm=T))

p1<-ggplot(disdf,aes(x=as.numeric(as.character(length)),y=value,color=time,fill=time))+
	#geom_bar(stat="identity",position="dodge",alpha=.5)+
	geom_line(alpha=.5)+
	ylab("Numbers of individuals")+xlab("length")+
	facet_grid(technical+fraction~space,scale="free")+
	theme(axis.text.x = element_text(size=8, angle=0),
	      strip.text.x=element_text(size=6,angle=0),
	      strip.text.y=element_text(size=6,angle=0))
	      #legend.position="bottom")
	p1

return(p1)

}

#dis vs lan
disvslan<-function(){
#ces commandes sont soumises à R 2.15.3 via le script 3_rejetvsretenu.r
script<-file(paste("3_rejetvsretenu.r",sep=""),open="wt")
writeLines("#utiliser R 2.15.3",con=script)
writeLines("library(COSTdbe);load('datavalcons.Rdata');load('myStr.Rdata');
	   load('info.Rdata')",con=script)
writeLines(paste0("dbelan <- dbeObject(species='",info$species,"',taxon='",info$taxon,"',
		  catchCat='LAN',strataDesc=myStr)"),con=script)
writeLines(paste0("dbedis <- dbeObject(species='",info$species,"',taxon='",info$taxon,"',
		  catchCat='DIS',strataDesc=myStr)"),con=script)
writeLines(paste0("dbelan<-RaiseLgth(dbelan,CSc,spp='",info$species,"',
		  taxon='",info$taxon,"')"),con=script)
writeLines(paste0("dbedis<-RaiseLgth(dbedis,CSc,spp='",info$species,"',
		  taxon='",info$taxon,"')"),con=script)
writeLines(paste0("disdf<-rbind(data.frame(dbedis@lenStruc$estim,fraction='DIS'),
		  data.frame(dbelan@lenStruc$estim,fraction='LAN'))"),con=script)
writeLines(paste0("rez<-tapply(disdf$value,list(disdf$length,disdf$fraction,disdf$time,disdf$technical),sum)"),con=script)
writeLines("save(disdf,rez,file='rejetvsretenu.Rdata')",con=script)
close(script)

system("wine Rscript.exe -e \"source('3_rejetvsretenu.r')\"")
load("rejetvsretenu.Rdata")
disdf<-disdf%>%mutate(time=gsub("20","",gsub(" ","",time)))
p1<-ggplot(disdf,aes(x=as.numeric(as.character(length)),y=value,fill=fraction,color=fraction))+
	geom_bar(stat="identity",position="dodge",alpha=.6)+
	ylab("Numbers of individuals")+xlab("length")+
	facet_grid(time~technical,scale="free")+
	theme(axis.text.x = element_text(size=8, angle=0),
	      strip.text.x=element_text(size=8,angle=0),
	      legend.position="bottom")
return(p1)



}



#age interpolation
interage<-function(CSc,info,param,rtp,maxage=10){
	#test area
	if(FALSE){
		source("credo_fct.R")
		require(ggplot2);require(pander)
		library(COSTdbe);library(COSTeda)
		library(dplyr)
		load("test.rdata")
		#merging test
		load("./data/CSr2016.rdata")
		load("./data/CLr2016.rdata")
		load("./data/CSrall.rdata")
		load("./data/CLrall.rdata")
		CSr<-subset(CSr,year%in%2014:2016,table="hh")
		CLr<-subset(CLr,year%in%2014:2016)
		load("myStr.Rdata")
		load("info.rdata")
		load("param.rdata")
		load("rtp.Rdata")
		CSr<-corrbase(CSr)
		rtp<-findrtp(CSr,CLr)
		CSr<-corrsampw(CSr,rtp)[[1]]
		CSr<-corrcaw(CSr,rtp)[[1]]
		CSc<-csDataCons(csDataVal(CSr),myStr)
		CLc<-clDataCons(clDataVal(CLr),myStr)
		load("datavalcons.Rdata")
		(table(CSc@ca$lenCls))
		(table(CSc@hl$lenCls))

		pipo<-interage(CSc,info,param,rtp)
		CSc2<-pipo[[1]]
		summary(CSc2@ca)
		sort(unique(CSc2@ca$lenCls))

		plot(pipo[[2]])
		plot(pipo[[3]])
		fitvb<-pipo[[4]]
		save(fitvb,file="fitvb.rdata")
		#CSc2@ca<-CSc2@ca%>%filter(is.na(technical))
		pipo<-dolen(CSc2,CLc,myStr,info,param)
		dbelan<-pipo[[1]];wdbelan<-pipo[[2]]
		aa<-wdbelan@lenStruc$estim%>%select(time,space,technical)%>%distinct()%>%mutate(tab="wbdelan")
		bb<-CSc2@ca%>%select(time,space,technical)%>%distinct()%>%mutate(tab2="CSc2")
		cc<-CSc@ca%>%select(time,space,technical)%>%distinct()%>%mutate(tab3="CSc")
		dd<-CSc@hl%>%select(time,space,technical)%>%distinct()%>%mutate(tab4="hl")
		full_join(bb,dd)%>%View()

		summary(wdbelan@lenStruc$estim)
		graphlen(dbelan,wdbelan,info,rtp,wrtp=F)

		#CSc2@ca<-CSc2@ca%>%filter(is.na(technical))
		pipo<-doage(CSc2,dbelan,wdbelan)
		dbelan<-pipo[[1]];wdbelan<-pipo[[2]]
		graphage(dbelan,wdbelan,info,rtp,checkn=0,checks=0,wrtp=T,fitvb)

	}

	#homogeneise steplenth
	#CSc<-alkLgthRec(CSc,type='stepIncr',param$stepIncr,preview=FALSE,postview=FALSE,update=TRUE)
	#graph
	p1<-ggplot(CSc@ca,aes(x=lenCls,y=indWt,type=stock))+geom_point(size=2) +facet_grid(space~time)
	#prepar suivi corr
	CSc@ca$stock<-as.character(CSc@ca$stock)
	#correction and completion age
	CSc@ca$stock[CSc@ca$age==-1] <-"corrNAage"
	CSc@ca$age[CSc@ca$age==-1]<-NA
	CSc@ca$age<-as.numeric(CSc@ca$age)
	#if(diff(range(CSc@ca$age,na.rm=T))>15){
	#	maxage<-quantile(CSc@ca$age,0.99,na.rm=T)
	#}else{
	#	maxage<-max(CSc@ca$age,na.rm=T)
	#}
	CSc@ca$age[CSc@ca$age>=maxage]<-maxage
	testnow<-!(CSc@ca$indWt!=-1 & !is.na(CSc@ca$indWt))
	CSc@ca$indWt[testnow]<-1000*rtp$a*(as.numeric(as.character(CSc@ca$lenCls[testnow]))/10)^rtp$b
	CSc@ca$stock[testnow]<-paste(CSc@ca$stock[testnow],"corrw")

	#complete missing age with data using von berta 
	catmp<-CSc@ca%>%mutate(age=ifelse(is.na(age),-1,age))
	p3<-ggplot(catmp,aes(x=age,y=lenCls))+geom_point() +facet_wrap(space~time)
	#p1
	#require(fishmethods)
	datvb<-data.frame(size=catmp$lenCls,age=catmp$age,
			trim=as.numeric(sapply (strsplit (as.character (catmp$time), " - "), FUN = function(x){x[2]})))%>%
			filter(age!=-1 & !is.na(size))%>%
			mutate(age=age+(trim-1)/4,meth="obs")%>%select(-trim)
	#simplification des données pour assurer la convergence...
	datvb<-datvb%>%group_by(age,meth)%>%summarise(size=median(size,na.rm=T))
	pltvb<-ggplot(datvb,aes(x=age,y=size))+geom_point()#+geom_line(data=datvb2,aes(x=age,y=size))
	modlm<-lm(age~size,data=datvb)
	fitvb<-nls(size~Sinf*(1-exp(-(K*(age-t0)))), data=datvb,
				control=nls.control(maxiter=10000,warnOnly=T,minFactor=1/1e100),
				algorithm="plinear",
				start=list(Sinf=max(datvb$size),K=coef(modlm)[2],t0=coef(modlm)[1]),
				#start=list(Sinf=quantile(datvb$size,0.95),K=coef(modlm)[2],t0=coef(modlm)[1]),
				trace=F)
	datvb2<-data.frame(size=predict(fitvb,newdata=data.frame(age=seq(0,maxage+1,0.1))),age=seq(0,maxage+1,0.1))
	#datvb2<-data.frame(size=predict(fitvb,newdata=data.frame(age=seq(0,100,0.1))),age=seq(0,100,0.1))
	pltvb<-ggplot(datvb,aes(x=age,y=size))+
		geom_point()+
		geom_line(data=datvb2,aes(x=age,y=size))+
		geom_point(data=CSc@ca,aes(x=age,y=lenCls),color="blue",alpha=.5)
	#pltvb
	#add missing age in the ca table
	#identification of all space time length in CSc
	#vecage<-seq(min(catmp$age[catmp$age!=-1],na.rm=T),max(catmp$age,na.rm=T),1)
	vecage<-seq(0,max(catmp$age,na.rm=T),1)
	aa<-list(time=unique(as.character(CSc@hl$time)),space=unique(as.character(CSc@hl$space)),
		   age=vecage)
	aa<-data.frame(expand.grid(aa))
	#limit them to the CSc@hl strata
	bb<-left_join(CSc@hl%>%select(space,time)%>%distinct,aa)
	nbfishplus<-10
	newind<-data.frame(do.call("rbind", replicate(nbfishplus, bb, simplify = FALSE)))
	#add trimester to age if any
	newind$trim<- as.numeric(sapply (strsplit (as.character (newind$time), " - "), FUN = function(x){x[2]}))
	newind$trim[is.na(newind$trim)]<-0
	newind$obsage<-newind$age
	newind$age<-newind$age+(newind$trim-1)/4
	newind$m2<-predict(fitvb,newind)
        statage<-catmp%>%filter(age>-1)%>%group_by(time,age)%>%summarise(m=mean(lenCls,na.rm=T),sd=sd(lenCls,na.rm=T))
	newind$s2<-mean(statage$sd,na.rm=T)
	newind<-left_join(newind,statage)%>%
		mutate(m=ifelse(is.na(m),m2,m),sd=ifelse(is.na(sd),s2,sd))
	newind$m3<-rnorm(nrow(newind),newind$m2,newind$sd)
	#convert length
	min1<-min(CSc@hl$lenCls)
	newind$m4<-round(newind$m3/param$stepIncr)*param$stepIncr
	pipo<-newind%>%transmute(PSUid=999,SSUid=999,time,space,technical=factor(NA),sort="LAN-HUC-NA",
				 sampType="V",landCtry="FRA",
				   vslFlgCtry="FRA",proj="BioPar",trpCode="",staNum="",spp=info$species,
				   sex=-1,stock="",lenCls=m4,age=trunc(age),fishId="",lenCode=paste(param$stepIncr,"mm"),
				   ageMeth=-1,plusGrp=NA,otoWt=-1,otoSide=NA,
				   indWt= 1000*rtp$a*(m3/10)^rtp$b,
				   matMeth="visual",
				   matScale="",matStage="")
	CSc@ca<-rbind(CSc@ca,pipo)
	#remove negative fishes
	CSc@ca$lenCls[CSc@ca$lenCls<=0]<-min(CSc@ca$lenCls[CSc@ca$lenCls>0])
	#interval length class
	intsl<-seq(min(CSc@hl$lenCls)%%param$stepIncr,max(c(CSc@hl$lenCls)),param$stepIncr)
	CSc@ca$lenCls<-intsl[findInterval(CSc@ca$lenCls,intsl)]
	#and just in case remove the duplicated lines
	CSc@ca<-distinct(CSc@ca)
	#rep rez by metier
	listmet<-na.omit(unique(CSc@hl$technical));nbmet<-length(listmet)
	camet<-data.frame(do.call("rbind", replicate(nbmet, CSc@ca, simplify = FALSE)))
	camet$technical<-rep(as.vector(listmet),each=nrow(CSc@ca))
	CSc@ca<-rbind(CSc@ca,camet)
	pltage<-ggplot(CSc@ca,aes(x=age,y=lenCls,shape=stock,color=gsub("20","",time),group=time))+
		geom_jitter(alpha=.6)+ 
		facet_wrap(~sub("_","\n",sub("_","-",gsub("_0_0","",technical))))
	#return stuff

	return(list(CSc,pltvb,pltage,fitvb))
}


##########################
#age interpolation
#multinomial approach with age borrowing from other quarter when needed
interage2<-function(CSc,info,param,rtp,maxage=10,p=10){
	#test area
	if(FALSE){
		source("credo_fct.R")
		require(ggplot2);require(pander)
		library(COSTdbe);library(COSTeda)
		library(dplyr)
		#data test
		load("info.rdata")
		load("param.rdata")
		load("rtp.rdata")
		print(load("datavalcons.Rdata"))
		maxage<-10
		p<-10
		#return(list(CSc,pltvb,pltage,fitvb))


	}
	#check data consistency
	pipo<-CSc@ca%>%transmute(year=substr(time,1,4),quarter=substr(time,8,8),space,sex,age,indWt,lenCls)
	ggplot(pipo,aes(x=age,y=lenCls,color=sex))+geom_point()+facet_grid(year~quarter+space)
	#si age -1 ou NA que faire ?
	#######################################################
	#adjust vb just in case
	datvb<-data.frame(size=CSc@ca$lenCls,age=CSc@ca$age,
			trim=as.numeric(sapply (strsplit (as.character (CSc@ca$time), " - "), FUN = function(x){x[2]})))%>%
			filter(age!=-1 & !is.na(size))%>%
			mutate(age=age+(trim-1)/4,meth="obs")%>%select(-trim)
	#simplification des données pour assurer la convergence...
	datvb<-datvb%>%group_by(age,meth)%>%summarise(size=median(size,na.rm=T))
	modlm<-lm(age~size,data=datvb)
	fitvb<-nls(size~Sinf*(1-exp(-(K*(age-t0)))), data=datvb,
				control=nls.control(maxiter=10000,warnOnly=T,minFactor=1/1e100),
				algorithm="plinear",
				start=list(Sinf=max(datvb$size),K=coef(modlm)[2],t0=coef(modlm)[1]),
				trace=F)
	datvb2<-data.frame(size=predict(fitvb,newdata=data.frame(age=seq(0,maxage+1,0.1))),age=seq(0,maxage+1,0.1))
	pltvb<-ggplot(datvb,aes(x=age,y=size))+
		geom_point()+
		geom_line(data=datvb2,aes(x=age,y=size))+
		geom_point(data=CSc@ca,aes(x=age,y=lenCls),color="blue",alpha=.5)
	#######################################################
	#multinom pour compléter
	strata<-unique(CSc@hl[,c("time","space")])
	#model complet si pas données alors bim
	multall <- multinom(age~lenCls,data=CSc@ca%>%filter(age>-1),trace=F)
	#plot(CSc@ca$age[CSc@ca$age>-1],predict(mult))
	for(i in 1:nrow(strata)){
		i<-1
		#filter data
		CA<-CSc@ca%>%
			filter(as.character(time)==strata$time[i]&
			       as.character(space)==strata$space[i])
		HL<-CSc@hl%>%
			filter(as.character(time)==strata$time[i]&
			       as.character(space)==strata$space[i])
		#multinom if possible
		if(nrow(CA)>1&length(unique(CA$age))>2){
			mult <- multinom(age~lenCls,data=CA%>%filter(age>-1),trace=F)
			#mult <- multinom(age~lenCls,data=data.frame(age=c(1,2,1,2,NA),lenCls=c(100,110,110,120,130)),trace=T)
		}

	}


	#homogeneise steplenth
	#CSc<-alkLgthRec(CSc,type='stepIncr',param$stepIncr,preview=FALSE,postview=FALSE,update=TRUE)
	#graph
	p1<-ggplot(CSc@ca,aes(x=lenCls,y=indWt,type=stock))+geom_point(size=2) +facet_grid(space~time)
	#prepar suivi corr
	CSc@ca$stock<-as.character(CSc@ca$stock)
	#correction and completion age
	CSc@ca$stock[CSc@ca$age==-1] <-"corrNAage"
	CSc@ca$age[CSc@ca$age==-1]<-NA
	CSc@ca$age<-as.numeric(CSc@ca$age)
	#if(diff(range(CSc@ca$age,na.rm=T))>15){
	#	maxage<-quantile(CSc@ca$age,0.99,na.rm=T)
	#}else{
	#	maxage<-max(CSc@ca$age,na.rm=T)
	#}
	CSc@ca$age[CSc@ca$age>=maxage]<-maxage
	testnow<-!(CSc@ca$indWt!=-1 & !is.na(CSc@ca$indWt))
	CSc@ca$indWt[testnow]<-1000*rtp$a*(as.numeric(as.character(CSc@ca$lenCls[testnow]))/10)^rtp$b
	CSc@ca$stock[testnow]<-paste(CSc@ca$stock[testnow],"corrw")

	#complete missing age with data using von berta 
	catmp<-CSc@ca%>%mutate(age=ifelse(is.na(age),-1,age))
	p3<-ggplot(catmp,aes(x=age,y=lenCls))+geom_point() +facet_wrap(space~time)
	#p1
	#require(fishmethods)
	datvb<-data.frame(size=catmp$lenCls,age=catmp$age,
			trim=as.numeric(sapply (strsplit (as.character (catmp$time), " - "), FUN = function(x){x[2]})))%>%
			filter(age!=-1 & !is.na(size))%>%
			mutate(age=age+(trim-1)/4,meth="obs")%>%select(-trim)
	#simplification des données pour assurer la convergence...
	datvb<-datvb%>%group_by(age,meth)%>%summarise(size=median(size,na.rm=T))
	pltvb<-ggplot(datvb,aes(x=age,y=size))+geom_point()#+geom_line(data=datvb2,aes(x=age,y=size))
	modlm<-lm(age~size,data=datvb)
	fitvb<-nls(size~Sinf*(1-exp(-(K*(age-t0)))), data=datvb,
				control=nls.control(maxiter=10000,warnOnly=T,minFactor=1/1e100),
				algorithm="plinear",
				start=list(Sinf=max(datvb$size),K=coef(modlm)[2],t0=coef(modlm)[1]),
				#start=list(Sinf=quantile(datvb$size,0.95),K=coef(modlm)[2],t0=coef(modlm)[1]),
				trace=F)
	datvb2<-data.frame(size=predict(fitvb,newdata=data.frame(age=seq(0,maxage+1,0.1))),age=seq(0,maxage+1,0.1))
	#datvb2<-data.frame(size=predict(fitvb,newdata=data.frame(age=seq(0,100,0.1))),age=seq(0,100,0.1))
	pltvb<-ggplot(datvb,aes(x=age,y=size))+
		geom_point()+
		geom_line(data=datvb2,aes(x=age,y=size))+
		geom_point(data=CSc@ca,aes(x=age,y=lenCls),color="blue",alpha=.5)
	#pltvb
	#add missing age in the ca table
	#identification of all space time length in CSc
	#vecage<-seq(min(catmp$age[catmp$age!=-1],na.rm=T),max(catmp$age,na.rm=T),1)
	vecage<-seq(0,max(catmp$age,na.rm=T),1)
	aa<-list(time=unique(as.character(CSc@hl$time)),space=unique(as.character(CSc@hl$space)),
		   age=vecage)
	aa<-data.frame(expand.grid(aa))
	#limit them to the CSc@hl strata
	bb<-left_join(CSc@hl%>%select(space,time)%>%distinct,aa)
	nbfishplus<-10
	newind<-data.frame(do.call("rbind", replicate(nbfishplus, bb, simplify = FALSE)))
	#add trimester to age if any
	newind$trim<- as.numeric(sapply (strsplit (as.character (newind$time), " - "), FUN = function(x){x[2]}))
	newind$trim[is.na(newind$trim)]<-0
	newind$obsage<-newind$age
	newind$age<-newind$age+(newind$trim-1)/4
	newind$m2<-predict(fitvb,newind)
        statage<-catmp%>%filter(age>-1)%>%group_by(time,age)%>%summarise(m=mean(lenCls,na.rm=T),sd=sd(lenCls,na.rm=T))
	newind$s2<-mean(statage$sd,na.rm=T)
	newind<-left_join(newind,statage)%>%
		mutate(m=ifelse(is.na(m),m2,m),sd=ifelse(is.na(sd),s2,sd))
	newind$m3<-rnorm(nrow(newind),newind$m2,newind$sd)
	#convert length
	min1<-min(CSc@hl$lenCls)
	newind$m4<-round(newind$m3/param$stepIncr)*param$stepIncr
	pipo<-newind%>%transmute(PSUid=999,SSUid=999,time,space,technical=factor(NA),sort="LAN-HUC-NA",
				 sampType="V",landCtry="FRA",
				   vslFlgCtry="FRA",proj="BioPar",trpCode="",staNum="",spp=info$species,
				   sex=-1,stock="",lenCls=m4,age=trunc(age),fishId="",lenCode=paste(param$stepIncr,"mm"),
				   ageMeth=-1,plusGrp=NA,otoWt=-1,otoSide=NA,
				   indWt= 1000*rtp$a*(m3/10)^rtp$b,
				   matMeth="visual",
				   matScale="",matStage="")
	CSc@ca<-rbind(CSc@ca,pipo)
	#remove negative fishes
	CSc@ca$lenCls[CSc@ca$lenCls<=0]<-min(CSc@ca$lenCls[CSc@ca$lenCls>0])
	#interval length class
	intsl<-seq(min(CSc@hl$lenCls)%%param$stepIncr,max(c(CSc@hl$lenCls)),param$stepIncr)
	CSc@ca$lenCls<-intsl[findInterval(CSc@ca$lenCls,intsl)]
	#and just in case remove the duplicated lines
	CSc@ca<-distinct(CSc@ca)
	#rep rez by metier
	listmet<-na.omit(unique(CSc@hl$technical));nbmet<-length(listmet)
	camet<-data.frame(do.call("rbind", replicate(nbmet, CSc@ca, simplify = FALSE)))
	camet$technical<-rep(as.vector(listmet),each=nrow(CSc@ca))
	CSc@ca<-rbind(CSc@ca,camet)
	pltage<-ggplot(CSc@ca,aes(x=age,y=lenCls,shape=stock,color=gsub("20","",time),group=time))+
		geom_jitter(alpha=.6)+ 
		facet_wrap(~sub("_","\n",sub("_","-",gsub("_0_0","",technical))))
	#return stuff

	return(list(CSc,pltvb,pltage,fitvb))
}


grrr<-function(){

	#complete missing age using rf
	rfFit<-randomForest::randomForest(age~lenCls,data=na.omit(CSc@ca%>%select(age,lenCls)))
	pipo<-predict(rfFit,newdata=data.frame(lenCls=CSc@ca$lenCls))
	plot(CSc@ca$lenCls,pipo)
	points(CSc@ca$lenCls,CSc@ca$age,pch=20)
	pipo<-predict(rfFit,newdata=data.frame(lenCls=CSc@ca$lenCls[is.na(CSc@ca$age)]))
	CSc@ca$stock[is.na(CSc@ca$age)]<-paste(CSc@ca$stock[is.na(CSc@ca$age)],"ageRF")
	CSc@ca$age[is.na(CSc@ca$age)]<-as.numeric(pipo)

	p1<-ggplot(CSc@ca,aes(y=lenCls,x=indWt,color=stock))+geom_point(alpha=.5) +facet_grid(space~time)
	p1
	p1<-ggplot(CSc@ca,aes(y=age,x=lenCls,color=stock))+geom_point(alpha=.5) +facet_grid(space~time)
	p1
	#rep ca by strate and metier
	#CSc@ca$indWt<-1000*rtp$a*(as.numeric(as.character(CSc@ca$lenCls))/10)^rtp$b
	CSc@ca<-rbind(CSc@ca,CSc@ca)


	return(list(CSc,p1))
}

grrr<-function(){

#*****************************************
	#correction and completion age
	testnow<-!(CSr@ca$indWt!=-1 & !is.na(CSr@ca$indWt))
	CSr@ca$indWt[testnow]<-1000*rtp$a*(as.numeric(as.character(CSr@ca$lenCls[testnow]))/10)^rtp$b
	#tabtest[7,2]<-TRUE
	CSr@ca<-CSr@ca%>%filter(!is.na(lenCls))
	CSr@ca<-CSr@ca%>%filter(lenCls>0)
	

	pipo<-CSr@ca%>%group_by(year,quarter,lenCls,age)%>%summarize(n=n())%>%ungroup()
	p1<-ggplot(pipo,aes(y=lenCls,x=age,fill=n))+geom_raster() +facet_grid(year~quarter)+
			scale_fill_distiller(palette='Spectral',name="Number of fish")
	p1<-ggplot(CSr@ca,aes(y=lenCls,x=indWt))+geom_point() +facet_grid(year~quarter)
	return(list(CSr,p1))

	#some correction
	CSc@ca<-CSc@ca%>%filter(!is.na(lenCls))
	CSc@ca<-CSc@ca%>%filter(lenCls>0)
	#complete missing age using rf
	#rfFit<-randomForest::randomForest(age~lenCls,data=na.omit(CSc@ca%>%select(age,lenCls)))
	#pipo<-predict(rfFit,newdata=data.frame(lenCls=CSc@ca$lenCls))
	#plot(CSc@ca$lenCls,pipo)
	#points(CSc@ca$lenCls,CSc@ca$age,pch=20)
	#pipo<-predict(rfFit,newdata=data.frame(lenCls=CSc@ca$lenCls[is.na(CSc@ca$age)]))
	#CSc@ca$age[is.na(CSc@ca$age)]<-as.numeric(pipo)
	#von Bertalanffy fit on ca
	#convert oldest fish age to the quantile age at 95
	CSc@ca$age<-as.numeric(CSc@ca$age)
	maxage<-quantile(CSc@ca$age,0.95)
	CSc@ca$age[CSc@ca$age>=maxage]<-maxage
	CSc@ca$age[CSc@ca$age==-1]<-NA
	CSc@ca<-CSc@ca%>%filter(!is.na(age))
	#a plot
	pipo<-CSc@ca%>%group_by(year=substr(time,1,4),quarter=substr(time,8,8),space,lenCls,age)%>%
		summarize(n=n())%>%ungroup()
	p1<-ggplot(pipo,aes(y=lenCls,x=age,fill=n))+geom_raster() +facet_grid(year~quarter+space) 
	p1

	#require(fishmethods)
	cakey<-paste(CSc@ca$time,CSc@ca$space,sep="/")
	dat<-data.frame(size=CSc@ca$lenCls,age=CSc@ca$age+(as.numeric(substr(CSc@ca$time,8,8))-1)*0.25,
					time=CSc@ca$time,space=CSc@ca$space,key=cakey)
	modlm<-lm(age~size,data=dat)
	fitvb<-nls(size~Sinf*(1-exp(-(K*(age-t0)))),
		   		     data=dat,
				     control=nls.control(maxiter=10000),
				start=list(Sinf=200,K=coef(modlm)[2],t0=coef(modlm)[1]),
				trace=FALSE)
	plot(dat[,2:1],main="age length model",xlim=c(0,maxage),ylim=c(0,max(CSc@ca$lenCls,na.rm=T)))
	lines(0:maxage,predict(fitvb,newdata=data.frame(age=0:maxage)))
	coef(fitvb)
	predict(fitvb,new=data.frame(age=c(0,0.25,0.5,0.75)))
	#list stepincr for future use
	unique(CSc@ca$lenCls)

	#check missing age strata relatively to CSc@sl
	requestlen<-CSc@hl%>%group_by(time,space,lenCls)%>%distinct()
	availlen<-CSc@ca%>%filter(!is.na(lenCls))%>%transmute(time=time,space=space,lenCls)%>%distinct()
	mislen<-anti_join(requestlen,availlen)
	if(nrow(mislen)>0){
		newind<-CSc@ca[1:nrow(mislen),]
		newind$space<-mislen$space
		newind$time<-mislen$time
		newind$lenCls<-mislen$lenCls
		newind$indWt<-1000*rtp$a*(newind$lenCls/10)^rtp$b
		newind$age<-NA
	}
	CSc@ca<-rbind(CSc@ca,newind)
	

	#check missing age strata (space and time)
	availage<-CSc@ca%>%filter(!is.na(age))%>%transmute(time=time,space=space,age=age)%>%distinct()
	#availage$age[is.na(availage$age)]<-""
	availage<-tidyr::spread(availage,"age","age",drop=F)
	eval(parse(text=paste0("availage<-tidyr::gather(availage,'age','n',",
		       paste(which(names(availage)%in%0:maxage),collapse=',')
		       ,")")))
	idmissing<-which(is.na(availage$n))
	#missing age
	misage<-CSc@ca%>%filter(is.na(age))%>%transmute(time,space,age,lenCls,indWt)




	sdlen<-CSc@ca%>%group_by(age)%>%summarise(sd=sd(lenCls,na.rm=T))%>%ungroup()%>%transmute(sd=mean(sd,na.rm=T))%>%distinct()
	nblen<-CSc@ca%>%group_by(age)%>%summarise(n=n())%>%ungroup()%>%transmute(n=mean(n))%>%distinct()
	nblen$n<-10

if(length(idmissing)>0){
	newind<-CSc@ca[1:length(idmissing),]
	newind$space<-availage$space[idmissing]
	newind$time<-availage$time[idmissing]
	newind$age<-availage$age[idmissing]
	#newind$lenCls<-predict(fitvb,new= data.frame(age=as.numeric(availage$age[idmissing])+(as.numeric(substr(availage$time[idmissing],8,8))-1)*0.25))
	newind2<-data.frame()
	#on génère 50 individus pour chaque âge manquant (rnorm avec sd obs)
	for(i in 1:nrow(newind)){
		newindtmp<-rep(newind[1,],nblen$n)
		newindtmp<-do.call("rbind", replicate(nblen$n, newind[i,], simplify = FALSE))
		newindtmp$lenCls<-rnorm(nrow(newindtmp),mean=mean(newindtmp$lenCls),sd=sdlen$sd)
		#pipo<-rnorm(nrow(newindtmp),mean=mean(newindtmp$lenCls),sd=sdlen$sd)
		#l0<-sort(unique(CSc@ca$lenCls))
		#sl0[findInterval(pipo,l0)]
		newindtmp$indWt<-1000*rtp$a*(newindtmp$lenCls/10)^rtp$b
		newind2<-rbind(newind2,newindtmp)
	}
}

newind2<-newind2[newind2$lenCls>0,]
newind2$age<-as.numeric(newind2$age)
CSc@ca<-rbind(CSc@ca,newind2)
#try(CSc<-fillALKmult(CSc,info$species,p=15),silent=T)

	CSc<-alkLgthRec(CSc,type='stepIncr',param$stepIncr,preview=FALSE,postview=FALSE,update=TRUE)

	pipo<-CSc@ca%>%group_by(year=substr(time,1,4),quarter=substr(time,8,8),lenCls,age)%>%
		summarize(n=n())%>%ungroup()
	p2<-ggplot(data.frame(pipo),aes(y=lenCls,x=age,fill=n))+geom_raster() +facet_grid(year~quarter) 
	p3<-ggplot(data.frame(pipo),aes(y=lenCls,x=age,fill=n))+geom_point() +facet_grid(year~quarter) 
#pipo<-CSc@ca%>%group_by(year=substr(time,1,4),quarter=substr(time,8,8),lenCls,age)%>%summarize(n=n())%>%ungroup()
#p2<-ggplot(pipo,aes(y=lenCls,x=age,fill=n))+geom_raster() +facet_grid(year~quarter) 
#pipo<-CSc@ca%>%mutate(year=substr(time,1,4),quarter=substr(time,8,8))
#ggplot(pipo,aes(x=age,y=lenCls,color=space))+geom_point()+facet_grid(year~quarter)

save(CSv,CSDv,CLv,CEv,CSc,CSDc,CLc,CEc,file="datavalcons2.Rdata")
return(list(CSc,p1,p2))



}


#compute age 
doage<-function(CSc,dbelan,wdbelan,param){

	#test area
	test<-F
	if(test){
		require(ggplot2);require(pander)
		library(COSTdbe);library(COSTeda)
		#CSr<-csData()
		load("datavalcons.Rdata")
		load("info.rdata")
		load("myStr.Rdata")

	}
	#print("age")
	dbelan<-RaiseAge(dbelan,CSc)
	#print("weight@age")
	if(tabtest[7,2]){
		wdbelan<-bpEstim(wdbelan,CSc,dbelan,adjust=TRUE)
	}
	return(list(dbelan=dbelan,wdbelan=wdbelan))
}
#compute dis in length and volume
dodis<-function(CSc,CLc,myStr,info,param,tabtest){
	if(F){
		require(ggplot2);require(pander)
		library(COSTdbe);library(COSTeda)
		#CSr<-csData()
		load("datavalcons.Rdata")
		load("info.rdata")
		load("myStr.Rdata")
		load("param.rdata")

	}
	#species in taxon due to COST bug
	CLc@cl$taxon<-as.character(CLc@cl$taxon)
	CLc@cl$taxon[CLc@cl$taxon==info$taxon]<-info$species
	#compute dis
	dbedis <- dbeObject(species=info$species,taxon=info$species,
			     catchCat='DIS',strataDesc=myStr,methodDesc='analytical')
	wdbedis<-dbeObject(species=info$species,taxon=info$species,
			   catchCat='DIS',strataDesc=myStr,
			   methodDesc='analytical',param='weight',adjust=F)
	dbedis<-totVolume(dbedis,CSc,CEc,CLc,type='landings')
	#trick to have bpEstim working (rep ca ntimes per technical strata)
	listmet<-unique(dbedis@lenStruc$estim$technical)
	catmp<-CSc@ca
	techtmp<-rep(listmet,each=nrow(catmp))
	length(rep(seq_len(nrow(catmp)),length(listmet)))
	catmp<-catmp[rep(seq_len(nrow(catmp)),length(listmet)),]
	catmp$technical<-techtmp
	CSc@ca<-catmp

	if(tabtest[7,2]){
		wdbedis<-bpEstim(wdbedis,CSc,dbedis,adjust=TRUE)
	}

	return(list(dbedis=dbedis,wdbedis=wdbedis))

	#dbedis@totalW$estim%>%filter(is.finite(value))

}

#compute  length
dolen<-function(CSc,CLc,myStr,info,param,tabtest){
	#test area
	if(F){
		require(ggplot2);require(pander)
		library(COSTdbe);library(COSTeda)
		#CSr<-csData()
		load("datavalcons.Rdata")
		load("info.rdata")
		load("myStr.Rdata")
		load("param.rdata")

	}
	CLc@cl$taxon<-as.character(CLc@cl$taxon)
	CLc@cl$taxon[CLc@cl$taxon==info$taxon]<-info$species
	#info$species
	#print("length")
	dbelan<-dbeObject(species=info$species,taxon=info$species,catchCat='LAN',strataDesc=myStr)
	dbelan<-RaiseLgth(dbelan,CSc,CLc,spp=info$species,taxon=info$species)
	#print("weight")
	wdbelan <-dbeObject(species=info$species,taxon=info$species,catchCat='LAN',strataDesc=myStr,
		methodDesc='analytical',param='weight',adjust=F)
	#trick to have bpEstim working (rep ca ntimes per technical strata)
	listmet<-unique(dbelan@lenStruc$estim$technical)
	catmp<-CSc@ca
	techtmp<-rep(listmet,each=nrow(catmp))
	length(rep(seq_len(nrow(catmp)),length(listmet)))
	catmp<-catmp[rep(seq_len(nrow(catmp)),length(listmet)),]
	catmp$technical<-techtmp
	CSc@ca<-catmp

	if(tabtest[7,2]){
		wdbelan<-bpEstim(wdbelan,CSc,dbelan,adjust=TRUE)
	}
	return(list(dbelan=dbelan,wdbelan=wdbelan))
}

#compute  length
dolen2<-function(CSc,CLc,myStr,info,param,allspecies){
	#test area
	test<-F
	if(test){
		require(ggplot2);require(pander)
		library(COSTdbe);library(COSTeda)
		#CSr<-csData()
		load("datavalcons.Rdata")
		load("info.rdata")
		load("myStr.Rdata")
		load("param.rdata")

	}
	CLc@cl$taxon<-info$species
	#print("length")
	dbelan<-dbeObject(species=info$species,taxon=info$species,catchCat='LAN',strataDesc=myStr)
	dbelan<-RaiseLgth(dbelan,CSc,CLc,spp=allspecies,taxon=info$species)
	#print("weight")
	wdbelan <-dbeObject(species=info$species,taxon=info$species,catchCat='LAN',strataDesc=myStr,
		methodDesc='analytical',param='weight',adjust=F)
	if(tabtest[7,2]&any(CSc@ca$spp%in%info$specie)){
		wdbelan<-bpEstim(wdbelan,CSc,dbelan,adjust=TRUE)
	}
	return(list(dbelan=dbelan,wdbelan=wdbelan))
}
#combine all info to build up check on estimation and compute sop
combineSISDCL<-function(SDtmp,SDagetmp,SI,CLc){
        #SOP + masse reelle
        SDverif<-SDtmp%>%mutate(wtot=NumberCaught*MeanWeight)%>%
                group_by(Year,Season,
                        Fleet,
                        FishingArea,
                        CatchCategory)%>%
                summarise(CATONSD=sum(wtot,na.rm=T))%>%
                ungroup()
        SDageverif<-SDagetmp%>%mutate(wtot=NumberCaught*MeanWeight)%>%
                group_by(Year,Season,
                        Fleet,
                        FishingArea,
                        CatchCategory)%>%
                summarise(CATONSDage=sum(wtot,na.rm=T))%>%
                ungroup()
        SIverif<-SI%>%select(Year,Season,Fleet,FishingArea,CatchCategory,s,n,CATON,OffLandings)
        if(param$timestratif=="Year"){
          CLcverif<-CLc@cl%>%mutate(Year=time,Season=time)
        }else{
          CLcverif<-CLc@cl%>%mutate(
                Year=sapply(strsplit(as.character(time)," - "),FUN = function(x){x[1]}),
                Season=sapply(strsplit(as.character(time)," - "),FUN = function(x){x[2]}))
        }
        CLcverif<-CLcverif%>%group_by(Year,Season,
                                      Fleet=as.character(technical),
                                      FishingArea=as.character(space),
                                      CatchCategory="L")%>%
                summarise(CATONcl=sum(landWt,na.rm=T))%>%ungroup()
        allverif<-full_join(SIverif,SDverif)%>%full_join(SDageverif)%>%full_join(CLcverif)%>%
                mutate(sopl=CATONSD/CATON,sopa=CATONSDage/CATON,diffCATONcl=CATON-CATONcl)

	return(allverif)
}
#CLc with stat rectangle
cons2rect<-function(CLc,CEc,myStr){
 #compute rect stat
        myStr0<-myStr
        myStr0@spaceStrata<-"rect"
        myStr0@timeStrata<-"month"
        CLv <- clDataVal(CLr); CLcrect <- clDataCons(CLv,myStr0)
        CEv <- ceDataVal(CEr); CEcrect <- ceDataCons(CEv,myStr0)
        lanrect<-CLcrect@cl%>%
		group_by(time,space,technical)%>%
		summarise(lan_kg=sum(landWt,na.rm=T))%>%
		ungroup()
	effrect <-CEcrect@ce%>%
		group_by(time,space,technical)%>%
		summarise_at(vars(trpNum,daysAtSea,foDur,effGtDays),sum,na.rm=T)%>%
		ungroup()
	return(list(lanrect,effrect))
}



#plot and table for SI : time series year
yearSI<-function(SIfinal){
        catches<-SIfinal%>%
		mutate(CATON=ifelse(OffLandings!=-9,OffLandings,CATON))%>%
		group_by(Year=as.numeric(Year),Fleet,FishingArea,CatchCategory)%>%
                summarise(CATON=sum(CATON/1000,na.rm=T))%>%
                ungroup()
        catchesall<-SIfinal%>%
		mutate(CATON=ifelse(OffLandings!=-9,OffLandings,CATON))%>%
		group_by(Year=as.numeric(Year),Fleet="all",FishingArea,CatchCategory)%>%
                summarise(CATON=sum(CATON/1000,na.rm=T))%>%
                ungroup()

	catchestmp<-rbind(catchesall,catches)
        pltcatches<-ggplot(data=catchestmp,aes(x=Year,y=CATON,fill=CatchCategory))+
                geom_bar(stat="identity")+
                facet_grid(Fleet~FishingArea,scale="free")+
                theme_bw()+ylab("W in tons")+
                theme(axis.text.x = element_text(size=6, angle=0),
                      axis.text.y = element_text(size=6, angle=0),
                      strip.text.x=element_text(size=6,angle=0),
                      strip.text.y=element_text(size=6,angle=0),
                      legend.position="bottom")
	catches<-catches%>%pivot_wider(names_from=CatchCategory,
				       values_from=CATON,
				       values_fill=list(CATON=0))#%>%
	#if no discards add an empty column
	if(!any(names(catches)%in%c("D"))){catches$D<-0}
	catches<-catches%>%mutate(tx=(D)/(L+D))
	return(list(catches,pltcatches))
}
yearHI<-function(HIfinal){
        effort<-HIfinal%>%group_by(Year,Fleet,FishingArea)%>%
                summarise(Effort=sum(Effort,na.rm=T))%>%
                ungroup()
        effortall<-HIfinal%>%group_by(Year,Fleet="all",FishingArea)%>%
                summarise(Effort=sum(Effort,na.rm=T))%>%
                ungroup()
	effort<-rbind(effortall,effort)
        plt<-ggplot(data=effort,aes(x=Year,y=Effort))+
                geom_bar(stat="identity")+
                facet_grid(Fleet~FishingArea,scale="free")+
                theme_bw()+
                theme(axis.text.x = element_text(size=6, angle=0),
                      axis.text.y = element_text(size=6, angle=0),
                      strip.text.x=element_text(size=6,angle=0),
                      strip.text.y=element_text(size=6,angle=0),
                      legend.position="bottom")
	return(list(effort,plt))
}
yearSD<-function(SD,param){
	if(F){
		library(dplyr)
		library(ggplot2)
		load("test.rdata")
		SD<-SDfinal

	}
 	lend<-SD%>%group_by(Sex,Year,Fleet,FishingArea,CatchCategory,
                                 AgeLength=as.numeric(AgeLength),Sex)%>%
                summarise(NumberCaught=sum(NumberCaught,na.rm=T))%>%
                ungroup()
 	lendall<-SD%>%group_by(Sex,Year,Fleet="all",FishingArea,CatchCategory,
                                 AgeLength=as.numeric(AgeLength))%>%
                summarise(NumberCaught=sum(NumberCaught,na.rm=T))%>%
                ungroup()
	lend<-rbind(lend,lendall)
        plt<-ggplot(data=lend,aes(x=AgeLength,
                              y=NumberCaught,
			      colour=CatchCategory,linetype=Sex))+
                geom_path()+
                facet_grid(Fleet~Year+FishingArea,scale="free") +
                theme_bw()+
                theme(axis.text.x = element_text(size=6, angle=0),
                      axis.text.y = element_text(size=6, angle=0),
                      strip.text.x=element_text(size=6,angle=0),
                      strip.text.y=element_text(size=6,angle=0),
                      legend.position="bottom")

        pltbis<-ggplot(data=lend,aes(x=AgeLength,
                              y=NumberCaught,colour=FishingArea,
			      linetype=CatchCategory))+
                geom_path()+
                facet_grid(Fleet+Sex~Year,scale="free") +
                theme_bw()+
                theme(axis.text.x = element_text(size=6, angle=0),
                      axis.text.y = element_text(size=6, angle=0),
                      strip.text.x=element_text(size=6,angle=0),
                      strip.text.y=element_text(size=6,angle=0),
                      legend.position="bottom")

	#cohort tracking, a try
 	lencoh<-SD%>%group_by(Year=as.numeric(Year),AgeLength=as.numeric(AgeLength),CatchCategory)%>%
                summarise(NumberCaught=sum(NumberCaught,na.rm=T))%>%
                ungroup()%>%
		mutate(Year=as.character(Year),type="D and L")
 	lencohtot<-SD%>%group_by(Year=as.numeric(Year),AgeLength=as.numeric(AgeLength),CatchCategory="all")%>%
                summarise(NumberCaught=sum(NumberCaught,na.rm=T))%>%
                ungroup()%>%
		mutate(Year=as.character(Year),type="total")
	allcohyear<-expand.grid(as.character(param$currentyear:(param$currentyear-param$nbyear)),c("D and L","total"))%>%transmute(Year=Var1,type=Var2)
	lencoh<-left_join(allcohyear,rbind(lencoh,lencohtot))#%>%filter(!is.na(CatchCategory))

        pltcoh<-ggplot()+
		geom_point(data=lencoh,aes(x=AgeLength,y=Year,size=NumberCaught,color=CatchCategory))+
		facet_wrap(~type)

	#another try (following the fao plot seen here
	##http://www.fao.org/3/T0535E/T0535E03.htm)
	#finally using geom_density_ridge
 	lencoh2<-SD%>%group_by(Sex,Year=as.character(Year),AgeLength=as.numeric(AgeLength),
				CatchCategory)%>%
                summarise(NumberCaught=sum(NumberCaught,na.rm=T))%>%
		ungroup()%>%
		group_by(Year)%>%mutate(totn=sum(NumberCaught))%>%
                ungroup()%>%mutate(NumberCaught=NumberCaught/totn)

	allcohyear<-data.frame(Year=as.character(param$currentyear:(param$currentyear-param$nbyear)))
	lencoh2<-left_join(allcohyear,lencoh2)
	pltcoh2<-ggplot(lencoh2,aes(x=AgeLength,y=Year,height=NumberCaught,
				    fill=paste(Sex,CatchCategory))) +
		ggridges::geom_density_ridges(stat="identity",alpha=.5)+
		ggridges::theme_ridges()

 	lencoh3<-SD%>%group_by(Sex,Year=as.character(Year),AgeLength=as.numeric(AgeLength))%>%
                summarise(NumberCaught=sum(NumberCaught,na.rm=T))%>%
		ungroup()%>%
		group_by(Year)%>%mutate(totn=sum(NumberCaught))%>%
                ungroup()%>%mutate(NumberCaught=NumberCaught/totn)

	allcohyear<-data.frame(Year=as.character(param$currentyear:(param$currentyear-param$nbyear)))
	lencoh3<-left_join(allcohyear,lencoh3)
	pltcoh3<-ggplot(lencoh3,aes(x=AgeLength,y=Year,height=NumberCaught,fill=Sex)) +
		ggridges::geom_density_ridges(stat="identity",alpha=.5)+
		ggridges::theme_ridges()


        return(list(plt,pltcoh,pltcoh2,pltcoh3,pltbis))
}


diagSDwl<-function(SD){
 	lenw<-SD%>%transmute(Sex,Year,Season,FishingArea,CatchCategory,
                                 AgeLength=as.numeric(AgeLength),MeanWeight)%>%
                distinct()
        pltw<-ggplot(data=lenw,aes(x=AgeLength,
                              y=MeanWeight,colour=paste(Sex,CatchCategory),shape=Season))+
                geom_point()+
                facet_grid(Year~FishingArea,scale="free") +
                theme_bw()+
                theme(axis.text.x = element_text(size=6, angle=0),
                      axis.text.y = element_text(size=6, angle=0),
                      strip.text.x=element_text(size=6,angle=0),
                      strip.text.y=element_text(size=6,angle=0),
                      legend.position="bottom")

 	agel<-SD%>%transmute(Sex,Year,Season,FishingArea,CatchCategory,
                                 AgeLength=as.numeric(AgeLength),MeanLength)%>%
                distinct()
        pltl<-ggplot(data=agel,aes(x=AgeLength,
                              y=MeanLength,colour=paste(Sex,CatchCategory),shape=Season))+
                geom_point()+
                facet_grid(Year~FishingArea,scale="free") +
                theme_bw()+
                theme(axis.text.x = element_text(size=6, angle=0),
                      axis.text.y = element_text(size=6, angle=0),
                      strip.text.x=element_text(size=6,angle=0),
                      strip.text.y=element_text(size=6,angle=0),
                      legend.position="bottom")
		pltl

        return(list(pltw=pltw,pltl=pltl))


}
vb4dis<-function(dbedis,dbelan,CSc,param,info){
	if(F){
		library(dplyr)
		library(ggplot2)
		load("test.rdata")
	}
	#identify strata where age is available 
	ageindata<-CSc@ca%>%select(time,space,age,lenCls)%>%distinct()%>%filter(age>=0)
	stratawithage<-ageindata%>%select(time,space)%>%distinct()
	#identify lenclass with no age information in dis
	lenindis<-dbedis@lenStruc$estim%>%filter(is.finite(value))%>%filter(value>0)%>%
		transmute(time,space,lenCls=as.numeric(length))%>%distinct()%>%
		semi_join(stratawithage)
	#if no data in len do nothing !
	if(nrow(lenindis)>0){
	#identify lenclass with no age information in lan 
	leninlan<-dbelan@lenStruc$estim%>%filter(is.finite(value))%>%filter(value>0)%>%
		transmute(time,space,lenCls=as.numeric(length))%>%distinct()%>%
		semi_join(stratawithage)
	agedismiss<-left_join(lenindis,ageindata)%>%mutate(type="dis")
	agelanmiss<-left_join(leninlan,ageindata)%>%mutate(type="lan")
	pipo<-rbind(agedismiss,agelanmiss)
	#plt of the initial situation
	pltinit<-ggplot()+geom_point(data=pipo%>%mutate(age=ifelse(is.na(age),-1,age)),
	       aes(x=lenCls,y=age,color=type,shape=type),alpha=.5)+
	   facet_grid(time~space,drop=T,scale="free")

	#identify size min in lan and dis for strata with age
	#sizeminhl<-CSc@hl%>%group_by(time,space,type=substr(sort,1,3))%>%
	#	summarise(minsize=min(lenCls))%>%ungroup()%>%
	#	semi_join(stratawithage)%>%tidyr::pivot_wider(names_from=type,values_from=minsize)

	#agedismiss<-left_join(lenindis,ageindata)%>%mutate(agedis=age)%>%select(-age)
	#agelanmiss<-left_join(leninlan,ageindata)%>%group_by(time,space)%>%
	#	summarise(sizelanmin=min(lenCls,na.rm=T),agelanmin=min(age,na.rm=T))%>%
	#	ungroup()
	#pipo<-left_join(agedismiss,agelanmiss)%>%
	#	mutate(flag=ifelse(lenCls<sizelanmin&is.na(agedis),"vb","no"))
		


	#vb model
        #add rtp info (NO SEX : to be adapted for sexed species)
        #need a von berta here
        agelen<-CSc@ca%>%filter(age>=0)%>%transmute(age,size=lenCls,type="obs")
        agelen0<-agelen%>%group_by(age)%>%summarise(size=median(size))%>%ungroup()%>%
                mutate(type="median")

        modlm<-lm(age~size,data=agelen0)
        fitvb<-nls(size~Sinf*(1-exp(-(K*(age-t0)))), data=agelen0, 
                   control=nls.control(maxiter=10000,warnOnly=T,minFactor=1/1e100), 
                   algorithm="plinear", 
                   start=list(Sinf=max(agelen$size),K=coef(modlm)[2],t0=coef(modlm)[1]),
                   trace=F)
        agemod<-seq(0,max(agelen0$age),0.05)
        agelen2<-data.frame(age=agemod,size=predict(fitvb,newdata=data.frame(age=agemod)),type="model")
	#compute length age for truncated length corresponding to lenCls in dbe
	#object
	lenclass<-sort(unique(agedismiss$lenCls))
	lenclass<-seq(min(lenclass),max(lenclass),by=param$stepIncr)
	agelen4dis<-agelen2%>%filter(size>=min(lenclass)&size<=max(lenclass))%>%
		mutate(lenCls=lenclass[findInterval(size,lenclass,all.inside=T)]) %>%
		group_by(lenCls)%>%summarise(agemean=round(mean(age),0))%>%
		ungroup()%>%
		transmute(age=agemean,size=lenCls,type="truncsize4dis")
        pipo<-rbind(agelen,agelen0,agelen2,agelen4dis)
        pltagesize<-ggplot(pipo,aes(x=age,y=size,color=type))+geom_point()
	#add fish for lencls without age
	newfish<-agedismiss%>%filter(is.na(age))%>%left_join(agelen4dis%>%transmute(lenCls=size,newage=age))
	#update the initial graph
	pltinit<-pltinit+geom_point(data=newfish,aes(x=lenCls,y=newage),col="black",alpha=1)

	pipo<-newfish%>%filter(!is.na(newage))%>%
		transmute(PSUid=999,SSUid=999,time=time,space=space,technical=NA,sort="LAN-HUC-NA",
				 sampType="V",landCtry="FRA",
				   vslFlgCtry="FRA",proj="BioPar",trpCode="",staNum="",spp=info$species,
				   sex=-1,stock="",lenCls=lenCls,age=newage,fishId=111,lenCode=paste(param$stepIncr,"mm"),
				   ageMeth=-1,plusGrp=NA,otoWt=-1,otoSide=NA,
				   indWt= NA,
				   matMeth="visual",
				   matScale="",matStage="")

	CSc@ca<-rbind(CSc@ca,pipo)
	}else{
		pltinit<-ggplot();pltagesize<-ggplot()
	}
	return(list(CSc,pltinit,pltagesize))

}


#add weight and length to icesdatage
addwl2age<-function(icesdat,wdbelan,wdbedis,ldbelan,ldbedis,rtp,CSc,param){
	if(F){
		library(dplyr)
		load("test.rdata")
	}
	if(param$timestratif=="Year"){
		wdbelan<-convtimedbe(wdbelan)
		wdbedis<-convtimedbe(wdbedis)
		ldbelan<-convtimedbe(ldbelan)
		ldbedis<-convtimedbe(ldbedis)
	}
    	wfishage<-rbind(wdbelan@ageStruc$estim%>%mutate(CatchCategory="L"),
                     wdbedis@ageStruc$estim%>%mutate(CatchCategory="D"))%>%
                transmute(Year=sapply(strsplit(as.character(time)," - "),FUN = function(x){x[1]}),
                          Season=sapply(strsplit(as.character(time)," - "),FUN = function(x){x[2]}),
                          Fleet=as.character(technical),
                          AgeLength=(as.character(age)),
                          FishingArea=space,CatchCategory,
                          wfish=value/1000)

    	lfishage<-rbind(ldbelan@ageStruc$estim%>%mutate(CatchCategory="L"),
                     ldbedis@ageStruc$estim%>%mutate(CatchCategory="D"))%>%
                transmute(Year=sapply(strsplit(as.character(time)," - "),FUN = function(x){x[1]}),
                          Season=sapply(strsplit(as.character(time)," - "),FUN = function(x){x[2]}),
                          Fleet=as.character(technical),
                          AgeLength=(as.character(age)),
                          FishingArea=space,CatchCategory,
                          lfish=value)

	#add to wfish and lfish
        SDagetmp<-left_join(icesdat$SDage,wfishage)%>%left_join(lfishage)
        #add rtp info (NO SEX : to be adapted for sexed species)
        #need a von berta here
        agelen<-CSc@ca%>%filter(age>=0)%>%transmute(age,size=lenCls,type="obs")
        agelen0<-agelen%>%group_by(age)%>%summarise(size=median(size))%>%ungroup()%>%
                mutate(type="median")

        modlm<-lm(age~size,data=agelen0)
        fitvb<-nls(size~Sinf*(1-exp(-(K*(age-t0)))), data=agelen0, 
                   control=nls.control(maxiter=10000,warnOnly=T,minFactor=1/1e100), 
                   algorithm="plinear", 
                   start=list(Sinf=max(agelen$size),K=coef(modlm)[2],t0=coef(modlm)[1]),
                   trace=F)
        agemod<-seq(0,max(agelen0$age),0.1)
        agelen2<-data.frame(age=agemod,size=predict(fitvb,newdata=data.frame(age=agemod)),type="model")
        pipo<-rbind(agelen,agelen0,agelen2)
        pltagesize<-ggplot(pipo,aes(x=age,y=size,color=type))+geom_point()+ylab("size (mm)")

        #add length at age
        SDagetmp$lfishvb<- predict(fitvb,newdata=data.frame(age=as.numeric(SDagetmp$AgeLength)))
        SDagetmp<-SDagetmp%>%mutate(MeanLength=ifelse(is.finite(lfish),lfish,lfishvb))

        #add weight at length using rtp if needed
        rtptmp<-rtp%>%filter(sex=="")%>%
                transmute(Season=as.character(as.numeric(quarter)),FishingArea=area,a,b)
	if(param$timestratif=="Year"){
		rtptmp<-rtptmp%>%select(-Season)
	}
	if(param$timestratif=="Year"){
		SDagetmp<-SDagetmp%>%left_join(rtptmp,by=c("FishingArea"))%>%
			mutate(wfishrtp=a*(as.numeric(MeanLength)/10)^b)%>%
			mutate(MeanWeight=ifelse(is.finite(wfish),wfish,wfishrtp))
	}else{
		SDagetmp<-SDagetmp%>%left_join(rtptmp,by=c("FishingArea","Season"))%>%
			mutate(wfishrtp=a*(as.numeric(MeanLength)/10)^b)%>%
			mutate(MeanWeight=ifelse(is.finite(wfish),wfish,wfishrtp))
	}
	#pb with MeanLength >1000 in Intercatch for mm so convert everythin to
	#cm
	SDagetmp<-SDagetmp%>%mutate(UnitMeanLength="cm",MeanLength=MeanLength/10)

        pwrtpage<-ggplot(SDagetmp,aes(x=wfishrtp,y=MeanWeight))+geom_point()
        plvbage<-ggplot(SDagetmp,aes(x=lfishvb,y=MeanLength))+geom_point()+ylab("size (cm)")
	return(list(SDage=SDagetmp,pwrtpage=pwrtpage,plvbage=plvbage,pltagesize=pltagesize))
}

#add weight to icesdat

addw2len<-function(icesdat,wdbelan,wdbedis,rtp,param){
	if(F){
		library(ggplot2);library(dplyr);library(COSTdbe)
		load("test.rdata")
	}
	if(param$timestratif=="Year"){
		wdbelan<-convtimedbe(wdbelan)
		wdbedis<-convtimedbe(wdbedis)
	}
    	wfishlen<-rbind(wdbelan@lenStruc$estim%>%mutate(CatchCategory="L"),
                     wdbedis@lenStruc$estim%>%mutate(CatchCategory="D"))%>%
                transmute(Year=sapply(strsplit(as.character(time)," - "),FUN = function(x){x[1]}),
                          Season=sapply(strsplit(as.character(time)," - "),FUN = function(x){x[2]}),
                          Fleet=as.character(technical),
                          AgeLength=(as.character(length)),
                          FishingArea=space,CatchCategory,
                          wfish=value/1000)


        #add to SD fishw
        SDtmp<-left_join(icesdat$SD,wfishlen)
        #add rtp info 
        rtptmp<-rtp%>%#filter(sex=="")%>%
                transmute(Sex=ifelse(sex=="","N",sex),
			  Season=as.character(as.numeric(quarter)),FishingArea=area,a,b)

	if(param$timestratif=="Year"){
		rtptmp<-rtptmp%>%select(-Season)
	}
	if(param$timestratif=="Year"){
		SDtmp<-SDtmp%>%left_join(rtptmp,by=c("Sex","FishingArea"))%>%
			mutate(wfishrtp=a*(as.numeric(AgeLength)/10)^b)%>%
			mutate(MeanWeight=ifelse(is.finite(wfish),wfish,wfishrtp))
	}else{
		SDtmp<-SDtmp%>%left_join(rtptmp,by=c("Sex","FishingArea","Season"))%>%
			mutate(wfishrtp=a*(as.numeric(AgeLength)/10)^b)%>%
			mutate(MeanWeight=ifelse(is.finite(wfish),wfish,wfishrtp))
	}

        pwrtp<-ggplot(SDtmp,aes(x=wfishrtp,y=MeanWeight,color=Sex))+geom_point()
	return(list(SD=SDtmp,pwrtp=pwrtp))
}


#compute and check SOP for lendata
soplen<-function(dbelan,wdbelan,rtp){
	if(F){
		require(dplyr)
		require(ggplot2)
		load("dbelan.Rdata")
		load("info.rdata")
		load("./data/CLrall.rdata")
		load("./data/CSrall.rdata")
		rtp<-findrtp(CLr,CSr)

	}
	sizefish<-dbelan@lenStruc$estim%>%
		mutate(time=as.character(time),
		       space=as.character(space),
		       technical=as.character(technical),
		       value=as.numeric(value))
	nfish<-dbelan@nMeas$len%>%transmute(time=as.character(time),
					    space=as.character(space),
					    technical=as.character(technical),
					    n=value)
	nsamp<-dbelan@nSamp$len%>%
		transmute(time=as.character(time),
			  space=as.character(space),
			  technical=as.character(technical),
			  s=value)
	wage<-wdbelan@lenStruc$estim%>%
		transmute(time=as.character(time),
			  space=as.character(space),
			  technical=as.character(technical),
			  length,
			  w=value)%>%filter(is.finite(w))
	sizefish<-left_join(sizefish,nfish,by=c("time","space","technical"))%>%
		left_join(nsamp,by=c("time","space","technical"))%>%
		left_join(wage,by=c("time","space","technical","length"))%>%
		mutate(year=substr(time,1,4),
		       quarter=gsub(" ","",(substr(time,7,8))),
		       length=as.numeric(as.character(length)))

	#add rtp info (NO SEX : to be adapted to sexed species)
	rtptmp<-rtp%>%filter(sex=="")%>%
		transmute(quarter=as.character(quarter),space=area,a,b)
	sizefish<-sizefish%>%left_join(rtptmp,by=c("quarter","space"))%>%
		mutate(wtot=value*a*(length/10)^b)




	#sop

	wtot<-sizefish%>%group_by(time,space,technical)%>%
		summarise(n=unique(n),s=unique(s),wtot=sum(wtot,na.rm=T))%>%ungroup()
	verif<-full_join(dbelan@totalW$estim,wtot)%>%ungroup()%>%
		mutate(SOP=round(wtot/value,2),value=round(value/1000,1),wtot=round(wtot/1000,1))
	wdbelan@lenStruc$estim%>%filter(time=="2009 - 1" & space=="27.7.d" & technical =="OTB_DEF_70-99_0_0")
	sizefish%>%filter(time=="2009 - 1" & space=="27.7.d" & technical =="OTB_DEF_70-99_0_0")
	sizefish%>%filter(time=="2005 - 4" & space=="27.7.d" & technical =="TBB_DEF_70-99_0_0_all")
	verif%>%filter(time=="2005 - 4" & space=="27.7.d" & technical =="TBB_DEF_70-99_0_0_all")
	dbelan@lenStruc$estim%>%filter(time=="2005 - 4" & space=="27.7.d" & technical =="TBB_DEF_70-99_0_0_all")
	CSc@hl%>%filter(time=="2005 - 4" & space=="27.7.d" & technical =="TBB_DEF_70-99_0_0_all")
	CSc@sl%>%filter(time=="2005 - 4" & space=="27.7.d" & technical =="TBB_DEF_70-99_0_0_all")

	#verif0<-verif%>%mutate(SOP=as.character(SOP))%>%select(time,space,technical,SOP)
	#sizefish<-left_join(sizefish,verif0)%>%mutate(nb=paste(s,n,sep="/"))
}


graphlen2<-function(dbelan,info){
	if(F){
		require(dplyr)
		require(ggplot2)
		load("dbelan.Rdata")
		load("info.rdata")
		rtp<-findrtp(CSr,info)
	}
	sizefish<-dbelan@lenStruc$estim%>%mutate(time=as.character(time),space=as.character(space),technical=as.character(technical))
	nfish<-dbelan@nMeas$len%>%transmute(time=as.character(time),space=as.character(space),technical=as.character(technical),n=value)
	nsamp<-dbelan@nSamp$len%>%transmute(time=as.character(time),space=as.character(space),technical=as.character(technical),s=value)
	sizefish<-left_join(left_join(sizefish,nfish),nsamp)%>%mutate(time=as.character(time),
		     space=as.character(space),
		     technical=as.character(technical),
		     year=substr(time,1,4),
		     quarter=substr(time,7,8),
		     value=as.numeric(value),
		     length=as.numeric(as.character(length)))
		plt2<-ggplot(sizefish%>%mutate(year=gsub("20","",year)),
			     #aes(x=length,y=value,color=year,linetype=space))+
			     aes(x=length,y=value,color=quarter,linetype=space))+
			ylab("Numbers of individuals")+xlab("length")+
			#facet_grid(technical~quarter,scale="free_y") +
			facet_grid(technical~year,scale="free_y") +
			geom_line(alpha=.6)+
			#geom_point(alpha=.6,size=1)+
			#geom_bar(alpha=.6,stat="identity", position=position_dodge())+
			ggtitle(paste0(info$stock))+
			theme(axis.text.x = element_text(size=6, angle=90),
			      axis.text.y = element_text(size=6, angle=0),
			      strip.text.x=element_text(size=6,angle=0),
			      strip.text.y=element_text(size=6,angle=0),
		      legend.position="bottom")
	return(plt2)

}



graphlen<-function(dbelan,wdbelan,info,rtp,checkn=50,checks=3,wrtp=T){
	if(F){
		require(dplyr)
		require(ggplot2)
		load("dbelan.Rdata")
		load("info.rdata")
		rtp<-findrtp(CSr,info)

	}

		

sizefish<-dbelan@lenStruc$estim%>%mutate(time=as.character(time),space=as.character(space),technical=as.character(technical))
nfish<-dbelan@nMeas$len%>%transmute(time=as.character(time),space=as.character(space),technical=as.character(technical),n=value)
	nsamp<-dbelan@nSamp$len%>%transmute(time=as.character(time),space=as.character(space),technical=as.character(technical),s=value)
	sizefish<-left_join(left_join(sizefish,nfish),nsamp)%>%mutate(time=as.character(time),
						     space=as.character(space),
						     technical=as.character(technical),
						     year=substr(time,1,4),
						     quarter=substr(time,7,8),
						     value=as.numeric(value),
						     length=as.numeric(as.character(length)))

	#SOP
	if(wrtp){
	#pipo<-full_join(aa,requestlen)
		sizefish$w<-sizefish$value*rtp$a*(sizefish$length/10)^rtp$b
	}else{
		#try to add w at age
		wage<-wdbelan@lenStruc$estim%>%mutate(time=as.character(time),
						      space=as.character(space),
						      technical=as.character(technical),
						      length=as.numeric(as.character(length)),w=value)%>%select(-value)
		sizefish<-left_join(sizefish,wage)%>%mutate(w=value*w/1000)
	}

	wtot<-sizefish%>%group_by(time,space,technical)%>%
		summarise(n=unique(n),s=unique(s),totw=sum(w,na.rm=T))%>%ungroup()
	verif<-full_join(dbelan@totalW$estim,wtot)%>%ungroup()%>%
		mutate(SOP=round(totw/value,2),value=round(value/1000,1),totw=round(totw/1000,1))
	verif0<-verif%>%mutate(SOP=as.character(SOP))%>%select(time,space,technical,SOP)
	sizefish<-left_join(sizefish,verif0)%>%mutate(nb=paste(s,n,sep="/"))
	sizefish<-sizefish%>%mutate(technical=paste0(substr(technical,1,7),"\n",gsub("_0","",substr(technical,8,nchar(technical)))))
	#remove strata according to the check
	sizefish<-sizefish%>%filter(n>=checkn&s>=checks)
	if(nrow(sizefish)>0){
	#graph par metier
		maxx<-max(sizefish$length)
		minx<-min(sizefish$length)
		maxy<-max(sizefish$value)
		miny<-min(sizefish$value)
		plt1<-ggplot(sizefish%>%mutate(year=gsub("20","",year)),aes(x=length,y=value))+
			ylab("Numbers of individuals")+xlab("length")+
			facet_grid(year+technical~space+quarter,scale="free_y") +
			#facet_grid(technical~space+quarter,scale="free_y") +
			geom_line(alpha=.6)+ggtitle(paste0(info$stock,"\n(nsamp/nmeas, SOP in blue)"))+
			geom_text(aes(x=maxx,y=miny,label=SOP),hjust=1,vjust=0,size=2,col="blue",alpha=.25)+
			geom_text(aes(x=minx,y=miny,label=nb),hjust=-.1,vjust=0,size=2,col="blue",alpha=.25)+
			theme(axis.text.x = element_text(size=8, angle=90),
			      axis.text.y = element_text(size=8, angle=0),
			      strip.text.x=element_text(size=10,angle=0),
			      strip.text.y=element_text(size=6,angle=0),
		      legend.position="right")
		plt2<-ggplot(sizefish%>%mutate(year=gsub("20","",year)),
			     aes(x=length,y=value,color=year,linetype=space))+
			ylab("Numbers of individuals")+xlab("length")+
			facet_grid(technical~quarter,scale="free_y") +
			#facet_grid(technical~space+quarter,scale="free_y") +
			geom_line(alpha=.6)+ggtitle(paste0(info$stock))+
			theme(axis.text.x = element_text(size=8, angle=90),
			      axis.text.y = element_text(size=8, angle=0),
			      strip.text.x=element_text(size=10,angle=0),
			      strip.text.y=element_text(size=6,angle=0),
		      legend.position="right")

	}else{
		plt1<-ggplot()+
		geom_text(aes(x=0,y=0,label=paste0("NO DATA for the tresholds ",checks,"/", checkn)))
		plt2<-ggplot()+
		geom_text(aes(x=0,y=0,label=paste0("NO DATA for the tresholds ",checks,"/", checkn)))
	}
	return(list(plt1,plt2))

}

graphage<-function(dbelan,wdbelan,info,rtp,checkn=50,checks=3,wrtp=T,fitvb){
	if(FALSE){
		load("5_agestruc.Rdata")
		
	}

	sizefish<-dbelan@ageStruc$estim%>%mutate(time=as.character(time),
						 space=as.character(space),technical=as.character(technical))
	#nfish<-dbelan@nMeas$age%>%transmute(time=as.character(time),space=as.character(space),n=value)
	#nsamp<-dbelan@nSamp$age%>%transmute(time=as.character(time),space=as.character(space),s=value)
	nfish<-dbelan@nMeas$len%>%transmute(time=as.character(time),space=as.character(space),
					    technical=as.character(technical),n=value)
	nsamp<-dbelan@nSamp$len%>%transmute(time=as.character(time),space=as.character(space),
					    technical=as.character(technical),s=value)
	sizefish<-left_join(left_join(sizefish,nfish),nsamp)%>%mutate(time=as.character(time),
						     space=as.character(space),
						     technical=as.character(technical),
						     year=substr(time,1,4),
						     quarter=substr(time,7,8),
						     value=as.numeric(value),
						     age=as.numeric(as.character(age)))
	#SOP

	if(wrtp){
	#pipo<-full_join(aa,requestlen)
		sizefish$w<-predict(fitvb,sizefish)
		sizefish$w<-sizefish$value*rtp$a*(sizefish$w/10)^rtp$b
	}else{
		#try to add w at age
		wage<-wdbelan@ageStruc$estim%>%mutate(time=as.character(time),
						      space=as.character(space),
						      technical=as.character(technical),
						      age=as.numeric(as.character(age)),w=value)%>%select(-value)
		sizefish<-left_join(sizefish,wage)%>%mutate(w=value*w/1000)

	}
	#try to add w at age
	#wage<-wdbelan@ageStruc$estim%>%mutate(age=as.numeric(as.character(age)),w=value)%>%select(-value)
	#sizefish<-left_join(sizefish,wage)

	#SOP
	wtot<-sizefish%>%group_by(time,space,technical)%>%
		summarise(n=unique(n),s=unique(s),totw=sum(w,na.rm=T))%>%ungroup()
	verif<-full_join(dbelan@totalW$estim,wtot)%>%ungroup()%>%
		mutate(SOP=round(totw/value,2),value=round(value/1000,1),totw=round(totw/1000,1))
	verif0<-verif%>%mutate(SOP=as.character(SOP))%>%select(time,space,technical,SOP)
	sizefish<-left_join(sizefish,verif0)%>%mutate(nb=paste(s,n,sep="/"))
	sizefish<-sizefish%>%mutate(technical=paste0(substr(technical,1,7),"\n",gsub("_0","",substr(technical,8,nchar(technical)))))
	#remove strata according to the check
	sizefish<-sizefish%>%filter(n>=checkn&s>=checks)
	if(nrow(sizefish)>0){
	#graph par metier
		maxx<-max(sizefish$age)
		minx<-min(sizefish$age)
		maxy<-max(sizefish$value)
		miny<-min(sizefish$value)
		#plt1<-ggplot(sizefish,aes(x=age,y=value))+
		#	ylab("Numbers of individuals")+xlab("Age")+
		#	facet_grid(year+technical~space+quarter,scale="free_y") +
		#	geom_line(alpha=.6)+ggtitle(paste0(info$stock,"\n(nsamp/nmeas, SOP in blue)"))+
		#	geom_point(alpha=.6)+
		#	geom_text(aes(x=maxx,y=miny,label=SOP),hjust=1,vjust=0,size=3,col="blue",alpha=.5)+
		#	geom_text(aes(x=minx,y=miny,label=nb),hjust=-.1,vjust=0,size=3,col="blue",alpha=.5)+
		#	theme(axis.text.x = element_text(size=8, angle=90),
		#	      axis.text.y = element_text(size=8, angle=0),
		#	      strip.text.x=element_text(size=10,angle=),
		#     legend.position="right")

		plt1<-ggplot(sizefish%>%mutate(year=gsub("20","",year)),aes(x=age,y=value))+
			ylab("Numbers of individuals")+xlab("Age")+
			facet_grid(year+technical~space+quarter,scale="free_y") +
			#facet_grid(technical~space+quarter,scale="free_y") +
			geom_line(alpha=.6)+ggtitle(paste0(info$stock,"\n(nsamp/nmeas, SOP in blue)"))+
			geom_text(aes(x=maxx,y=miny,label=SOP),hjust=1,vjust=0,size=2,col="blue",alpha=.25)+
			geom_text(aes(x=minx,y=miny,label=nb),hjust=-.1,vjust=0,size=2,col="blue",alpha=.25)+
			theme(axis.text.x = element_text(size=8, angle=90),
			      axis.text.y = element_text(size=8, angle=0),
			      strip.text.x=element_text(size=10,angle=0),
			      strip.text.y=element_text(size=6,angle=0),
		      legend.position="right")

		plt2<-ggplot(sizefish%>%mutate(year=gsub("20","",year)),
			     aes(x=age,y=value,color=year,linetype=space))+
			ylab("Numbers of individuals")+xlab("Age")+
			facet_grid(technical~quarter,scale="free_y") +
			#facet_grid(technical~space+quarter,scale="free_y") +
			geom_line(alpha=.6)+ggtitle(paste0(info$stock))+
			theme(axis.text.x = element_text(size=8, angle=90),
			      axis.text.y = element_text(size=8, angle=0),
			      strip.text.x=element_text(size=10,angle=0),
			      strip.text.y=element_text(size=6,angle=0),
		      legend.position="right")




	}else{
		plt1<-ggplot()+
		geom_text(aes(x=0,y=0,label=paste0("NO DATA for the tresholds ",checks,"/", checkn)))
		plt2<-ggplot()+
		geom_text(aes(x=0,y=0,label=paste0("NO DATA for the tresholds ",checks,"/", checkn)))
	}
	return(list(plt1,plt2))

}

#compute the delta and graph
delta<-function(CSr,myStr,info){
	#test area
	test<-F
	if(test){
		require(ggplot2);require(pander)
		library(COSTdbe);library(COSTeda)
		load("./data/CSrall.rdata")
		load("./data/CLrall.rdata")
		load("info.rdata")
		load("myStr.Rdata")
		CSr<-corrbase(CSr)
		rtp<-findrtp(CSr,CLr)
		rez<-corrsampw(CSr,rtp)
		CSr<-rez[[1]]
		table(CSr@hl$sampType,CSr@hl$year)
		CLr@cl<-CLr@cl%>%filter(year%in%2016:2014)
		CSr<-subset(CSr,year%in%2014:2016,table="hh")
	#CSr<-csData()
	}
	ca<-CSr@ca

	#raiselength annuel par metier
	CSv <- csDataVal(CSr); CSc <- csDataCons(CSv,myStr)
	#CSr@hh$foCatEu6<-CSc@hh$technical
	if(nrow(CSr@hl)>1){
		deltalan<-deltadis<-data.frame()
		if(any(CSr@hl$catchCat=='LAN')){
			deltalan<-deltCalc(CSr,myStr,species=info$species,fraction='LAN',strategy='metier')
			delta<-deltalan@outPut$SampDeltaMat%>%mutate(SampNum=as.numeric(samp))
			idsamp<-deltalan@outPut$DFsamp
			pipo<-full_join(delta,idsamp)%>%transmute(trpCode,staNum=as.numeric(staNum),delta,spp)
			deltalan<-left_join(pipo,CSc@sl%>%select(time,space,technical,trpCode,staNum,spp)%>%distinct())%>%mutate(type="LAN")
		}
		if(any(CSr@hl$catchCat=='DIS')){
			deltadis<-deltCalc(CSr,myStr,species=info$species,fraction='DIS',strategy='metier')
			delta<-deltadis@outPut$SampDeltaMat%>%mutate(SampNum=as.numeric(samp))
			idsamp<-deltadis@outPut$DFsamp
			pipo<-full_join(delta,idsamp)%>%transmute(trpCode,staNum=as.numeric(staNum),delta,spp)
			deltadis<-left_join(pipo,CSc@sl%>%select(time,space,technical,trpCode,staNum,spp)%>%distinct())%>%mutate(type="DIS")
		}

		delta<-rbind(deltalan,deltadis)%>%mutate(check=T,year=substr(time,1,4),quarter=substr(time,8,8),technical=as.character(technical))
		#outlier detection
		listtc<-(unique(delta$technical))
		listyear<-unique(delta$year)
		listtype<-unique(delta$type)
		for(i in listtc){
			for(j in listyear){
				for(k in listtype){
					#print(paste(i,j,k))
					test<-delta$technical==i & delta$year==j & delta$type==k
					if(any(test) & nrow(delta[test,])>3){
						delta$check[test]<-robustbase::adjOutlyingness(delta$delta[test],clower=0,cupper=0,alpha=.75)$nonOut
					}
				}
			}
		}
		pltdelta<-ggplot(delta%>%mutate(),aes(x=quarter,y=delta,group=check,color=check,shape=type))+
			geom_jitter(alpha=.6)+
			theme(strip.text.x = element_text(size=8, angle=0),
			      strip.text.y = element_text(size=6, angle=0),
			      )+
			facet_grid(sub("_","\n",sub("_","",gsub("_0_0","",technical)))~gsub("20","",year),scale="free_y")+xlab("")+ylab("")

		deltaid<-delta%>%
			filter(!check)%>%distinct()%>%
			select(trpCode,staNum)%>%mutate(key=paste(trpCode,staNum))%>%distinct()
		key<-paste(CSr@hh$trpCode,CSr@hh$staNum)
		CSr@hh$foVal[key%in%deltaid$key]<-"I"
		CSr<-subset(CSr,foVal=="V",table="hh",link=F)
		return(list(CSr,pltdelta))
	}
}



dorawlenhist<-function(CSr,info,param){
	#test area
	test<-F
	if(test){
		require(ggplot2);require(pander)
		library(COSTdbe);library(COSTeda)
		load("./data/CSrall.rdata")
		load("./data/CLrall.rdata")
		load("info.rdata")
		CSr<-corrbase(CSr)
		rtp<-findrtp(CSr,CLr)
		rez<-corrsampw(CSr,rtp)
		CSr<-rez[[1]]
	}
	myStr<-strIni(timeStrata='year', techStrata='foCatEu6')#,spaceStrata='area')
	#raiselength annuel par metier
	CSv <- csDataVal(CSr); CSc <- csDataCons(CSv,myStr)
	if(nrow(CSc@ca)>1){
	CSc<-alkLgthRec(CSc,type='stepIncr',param$stepIncr,preview=FALSE,postview=FALSE,update=TRUE)
	}
	dbelan <- dbeObject(species=info$species,taxon=info$taxon,
			  catchCat='LAN',strataDesc=myStr)
	dbelan<-RaiseLgth(dbelan,CSc,spp=info$species,taxon=info$taxon)
	sizefish<-dbelan@lenStruc$estim%>%transmute(year=as.numeric(substr(time,1,4)),
				 quarter=as.numeric(substr(time,8,8)),
				 space=as.character(space),
				 technical=as.character(technical),length=as.numeric(as.character(length)),
				 value=as.numeric(value))
	nfish<-dbelan@nMeas$len%>%transmute(year=as.numeric(substr(time,1,4)),
			    quarter=as.numeric(substr(time,8,8)),space=as.character(space),
			technical=as.character(technical),n=as.numeric(value))
	sizefish<-left_join(sizefish,nfish)
	sizefish<-sizefish%>%group_by(technical,year,quarter)%>%mutate(tot=sum(value,na.rm=T))%>%ungroup()%>%mutate(prop=value/tot)
	return(sizefish)
	
}

dorawmapmet<-function(CLr,CSr,CEr){
	#test area
	test<-F
	if(test){
		require(ggplot2);require(pander)
		library(COSTdbe);library(COSTeda)
		load("./data/CSrall.rdata")
		load("./data/CLrall.rdata")
		load("./data/CErall.rdata")
		load("info.rdata")
		CSr<-corrbase(CSr)
		rtp<-findrtp(CSr,CLr)
		rez<-corrsampw(CSr,rtp)
		CSr<-rez[[1]]
		CLr@cl<-CLr@cl%>%filter(year%in%2016:2014)
		CEr@ce<-CEr@ce%>%filter(year%in%2016:2014)
		CSr<-subset(CSr,year%in%2014:2016,table="hh")
		CSr<-csData()

	nbech<-CSr@sl%>%select(trpCode,staNum,catchCat,year)%>%distinct()%>%
		group_by(year,catchCat)%>%summarise(nb_sampled_haul=n())%>%ungroup()
	nbfish<-CSr@hl%>%select(lenNum,catchCat,year)%>%group_by(year,catchCat)%>%
		summarise(nb_fish=sum(lenNum,na.rm=T))%>%ungroup()

	nbsamp<-tbl_df(left_join(CSr@sl,CSr@hh)) %>% mutate(samp=paste(trpCode,staNum))%>%
		group_by(year,foCatEu6,catchCat)%>%
		summarise(nbsamp=n_distinct(samp))%>%
		ungroup()
	pipo<-full_join(nbsamp,nbech) %>%
		mutate(time=paste(year,quarter))%>%
		select(area,time,foCatEu6,catchCat,nbsamp)%>%filter(!grepl(",",area))
	}

	areafish<-tbl_df(CLr@cl)%>%
		group_by(area,year,quarter,foCatEu6)%>%
		summarise(lan=sum(landWt,na.rm=T)/1000)%>%ungroup()%>%
		group_by(foCatEu6,area)%>%mutate(tot=sum(lan,na.rm=T),prop=lan/tot)%>%
		ungroup()%>%
		mutate(time=paste(year,quarter))%>%
		select(area,time,foCatEu6,prop,lan)

	effortfish<-tbl_df(CEr@ce)%>%
		group_by(area,year,quarter,foCatEu6)%>%
		summarise(days=sum(daysAtSea,na.rm=T))%>%ungroup()%>%
		group_by(foCatEu6,area)%>%mutate(totdays=sum(days,na.rm=T),propdays=days/totdays)%>%
		ungroup()%>%
		mutate(time=paste(year,quarter))%>%
		select(area,time,foCatEu6,propdays,days)

	nbfish<-tbl_df(left_join(left_join(CSr@hl,CSr@sl),CSr@hh)) %>%
		group_by(area,year,quarter=ceiling(as.numeric(substr(date,6,7))/3),foCatEu6,catchCat)%>%
		summarise(nb=sum(lenNum,na.rm=T))%>%
		ungroup()%>%
		mutate(time=paste(year,quarter))%>%
		select(area,time,foCatEu6,catchCat,nb)%>%filter(!grepl(",",area))
	nbsamp<-tbl_df(left_join(CSr@sl,CSr@hh)) %>% mutate(samp=paste(trpCode,staNum))%>%
		group_by(area,year,quarter=ceiling(as.numeric(substr(date,6,7))/3),foCatEu6,catchCat)%>%
		summarise(nbsamp=n_distinct(samp))%>%
		ungroup() %>%
		mutate(time=paste(year,quarter))%>%
		select(area,time,foCatEu6,catchCat,nbsamp)%>%filter(!grepl(",",area))
	othermet<-tbl_df(CSr@hh)%>%
		transmute(area,year,quarter=ceiling(as.numeric(substr(date,6,7))/3),foCatEu6)%>%
		mutate(time=paste(year,quarter))%>%
		select(area,time,foCatEu6)%>%filter(!grepl(",",area))%>%distinct()


	#areafish<-full_join(full_join(full_join(areafish,nbfish),nbsamp),effortfish)%>%
	areafish<-full_join(full_join(areafish,nbfish),nbsamp)%>%
		full_join(othermet)%>%
	#areafish<-left_join(areafish,effortfish)%>%
		ungroup()%>%filter(!is.na(foCatEu6))%>%
		mutate(engin=substr(foCatEu6,1,3),area=gsub("27.","",area),
			foCatEu6=gsub(",","\n",foCatEu6))

	return(areafish)

}

dorawdisop<-function(CLr,CEr,CSr,yeartmp){
	#test area
	if(F){
		require(ggplot2);require(pander)
		library(COSTdbe);library(COSTeda)
		load("./data/CSrall.rdata")
		load("./data/CLrall.rdata")
		load("./data/CErall.rdata")
		load("info.rdata")
		yeartmp<-2016
		CSr<-corrbase(CSr)

		rtp<-findrtp(CSr,CLr)
		rez<-corrsampw(CSr,rtp)
		CSr<-rez[[1]]
		CLr@cl<-CLr@cl%>%filter(year%in%2016:2014)
		CSr<-subset(CSr,year%in%2014:2016,table="hh")
		CSr<-csData()
	}
	#test area gloub
	myStr<-strIni(timeStrata='year', techStrata='foCatEu6')#,spaceStrata='area')
	#raiselength annuel par metier
	CSv <- csDataVal(CSr); CSc <- csDataCons(CSv,myStr)
	if(nrow(CSc@ca)>1){
	CSc<-alkLgthRec(CSc,type='stepIncr',param$stepIncr,preview=FALSE,postview=FALSE,update=TRUE)
	}

	myStr<-strIni(timeStrata='quarter', techStrata='trpCode',spaceStrata='area')
	CSv <- csDataVal(CSr); CSc <- csDataCons(CSv,myStr)
	dbelan <- dbeObject(species=info$species,taxon=info$taxon,
			  catchCat='DIS',strataDesc=myStr)
	dbedis<-RaiseLgth(dbelan,CSc,spp=info$species,taxon=info$taxon)
	dbelan <- dbeObject(species=info$species,taxon=info$taxon,
			  catchCat='LAN',strataDesc=myStr)
	dbelan<-RaiseLgth(dbelan,CSc,spp=info$species,taxon=info$taxon)
	a1<-dbedis@totalW$estim%>%mutate(dis=value)%>%select(-value)
	a2<-dbelan@totalW$estim%>%mutate(lan=value)%>%select(-value)
	aa<-full_join(a1,a2)
	aa[is.na(aa)]<-0
	aa<-aa%>%mutate(tx=round(dis/(lan+dis),2))
	#add metier
	idmet<-CSr@hh%>%group_by(technical=trpCode)%>%summarise(met=paste(unique(foCatEu6),collapse=","))
	aa<-left_join(aa,idmet)

	pipo<-left_join(CSr@sl,CSr@hh)%>%group_by(foCatEu6,trpCode,catchCat)%>%summarise(wt=sum(wt,na.rm=T))
	pipo<-tidyr::spread(pipo,catchCat,wt)
	pipo[is.na(pipo)]<-0
	pipo<-pipo%>%mutate(tx=DIS/(LAN+DIS))

	CSr<-subset(CSr,year%in%2014:2016,table="hh")

	aa<-landisVol(CSr,species=info$species,fraction="DIS")

	ggplot(pipo,aes(y=foCatEu6,x=tx))+geom_point()
	dbedis <- dbeObject(species=info$species,taxon=info$species,
			  catchCat='DIS',strataDesc=myStr,methodDesc='analytical')

	#subset in time
	CSr<-eval(parse(text=paste0("subset(CSr,year%in%",yeartmp,",table='hh')")))
	CLr<-eval(parse(text=paste0("subset(CLr,year%in%",yeartmp,")")))
	CEr<-eval(parse(text=paste0("subset(CEr,year%in%",yeartmp,")")))
	myStr<-strIni(techStrata='foCatEu6',spaceStrat='area')
	CSc<-csDataCons(csDataVal(CSr),myStr)
	CEc<-ceDataCons(ceDataVal(CEr),myStr)
	CLc<-clDataCons(clDataVal(CLr),myStr)

	dbedis <- dbeObject(species=info$species,taxon=info$species,
			  catchCat='DIS',strataDesc=myStr,methodDesc='analytical')
	dbedistrip<-totVolume(dbedis,CSc,CEc,type='trip',val='weight')

	if(any(grepl("DIS",CSc@sl$sort))){
		dbedistrip<-totVolume(dbedis,CSc,CEc,type='trip',val='weight')
		dbedisfo<-totVolume(dbedis,CSc,CEc,type='fo',val='weight')
		dbedisfd<-totVolume(dbedis,CSc,CEc,type='fd',val='weight')
		dbedislan<-totVolume(dbedis,CSc,CEc,CLc,type='landings',val='weight')
	}
	distrip<-dbedistrip@totalW$estim%>%mutate(distrip=value)%>%select(-value)
	disfd<-dbedisfd@totalW$estim%>%mutate(disfd=value)%>%select(-value)
	dislan<-dbedislan@totalW$estim%>%mutate(dislan=value)%>%select(-value)
	lanall<-CLc@cl%>%group_by(time,space,technical)%>%summarise(tot=sum(landWt,na.rm=T))
	disall<-full_join(full_join(distrip,disfd),dislan)
	landisall<-left_join(disall,lanall)%>%
		mutate(txtrip=distrip/(distrip+tot),txfd=disfd/(disfd+tot),txlan=dislan/(dislan+tot)) %>%
		mutate(distrip=round(distrip/1000,2),disfd=round(disfd/1000,2),dislan=round(dislan/1000,2),
		       tot=round(tot/1000,2),txtrip=round(txtrip,2),txfd=round(txfd,2),txlan=round(txlan,2))
	return(landisall)
}

dorawdismet<-function(CLr,CEr,CSr,yeartmp){
	#test area
	if(F){
		require(ggplot2);require(pander)
		library(COSTdbe);library(COSTeda)
		load("./data/CSrall.rdata")
		load("./data/CLrall.rdata")
		load("./data/CErall.rdata")
		load("info.rdata")
		yeartmp<-2016
		CSr<-corrbase(CSr)

		rtp<-findrtp(CSr,CLr)
		rez<-corrsampw(CSr,rtp)
		CSr<-rez[[1]]
		CLr@cl<-CLr@cl%>%filter(year%in%2016:2014)
		CSr<-subset(CSr,year%in%2014:2016,table="hh")
		CSr<-csData()
	}
	#subset in time
	CSr<-eval(parse(text=paste0("subset(CSr,year%in%",yeartmp,",table='hh')")))
	CLr<-eval(parse(text=paste0("subset(CLr,year%in%",yeartmp,")")))
	CEr<-eval(parse(text=paste0("subset(CEr,year%in%",yeartmp,")")))
	myStr<-strIni(techStrata='foCatEu6',spaceStrat='area')
	CSc<-csDataCons(csDataVal(CSr),myStr)
	CEc<-ceDataCons(ceDataVal(CEr),myStr)
	CLc<-clDataCons(clDataVal(CLr),myStr)

	dbedis <- dbeObject(species=info$species,taxon=info$species,
			  catchCat='DIS',strataDesc=myStr,methodDesc='analytical')
	if(any(grepl("DIS",CSc@sl$sort))){
		dbedistrip<-totVolume(dbedis,CSc,CEc,type='trip',val='weight')
		dbedisfo<-totVolume(dbedis,CSc,CEc,type='fo',val='weight')
		dbedisfd<-totVolume(dbedis,CSc,CEc,type='fd',val='weight')
		dbedislan<-totVolume(dbedis,CSc,CEc,CLc,type='landings',val='weight')
	}
	distrip<-dbedistrip@totalW$estim%>%mutate(distrip=value)%>%select(-value)
	disfd<-dbedisfd@totalW$estim%>%mutate(disfd=value)%>%select(-value)
	dislan<-dbedislan@totalW$estim%>%mutate(dislan=value)%>%select(-value)
	lanall<-CLc@cl%>%group_by(time,space,technical)%>%summarise(tot=sum(landWt,na.rm=T))
	disall<-full_join(full_join(distrip,disfd),dislan)
	landisall<-left_join(disall,lanall)%>%
		mutate(txtrip=distrip/(distrip+tot),txfd=disfd/(disfd+tot),txlan=dislan/(dislan+tot)) %>%
		mutate(distrip=round(distrip/1000,2),disfd=round(disfd/1000,2),dislan=round(dislan/1000,2),
		       tot=round(tot/1000,2),txtrip=round(txtrip,2),txfd=round(txfd,2),txlan=round(txlan,2))
	return(landisall)
}


#basic correction
corrbase<-function(CSr){
	#validation obsmer
	CSr@hh$foVal[CSr@hh$sampType=="M"]<-"V"
	#lenCode in mm
	CSr@sl$lenCode<-"mm"
	#correction aggregation level
	CSr@hh$aggLev[CSr@hh$aggLev=="TRUE"]<-T
	#si subSampWt ou wt  = 0 alors NA
	CSr@sl$wt[CSr@sl$wt<=0]<-NA
	CSr@sl$subSampWt[CSr@sl$subSampWt<=0]<-NA
	#remove NA lenCls in CSr@hl
	CSr@hl<-CSr@hl%>%filter(is.finite(lenCls))
	return(CSr)
}
pipo<-function(){
	library(dplyr);library(ggplot2)
	load("./data/CLrall.rdata")
	load("./data/CSrall.rdata")
	rtp<-findrtp(CSr,CLr)
	
}

#correction sample weight
corrsampw<-function(CSr,rtp){
	#test area
	test<-F
	if(test){
		require(ggplot2);require(pander)
		require(dplyr)
		load("./data/CSrall.rdata")
		load("./data/CLrall.rdata")
		CSr<-corrbase(CSr)
		rtp<-findrtp(CLr,CSr)
		pipo<-corrsampw(CSr,rtp)

	}
	#preparing output for data quality in SIH
	idinvalidsih<-data.frame()

	#0. estimate subsamp weigh with rtp
	hlslhh<-left_join(left_join(CSr@hl,CSr@sl),CSr@hh)%>%mutate(quarter=ceiling(as.numeric(substr(date,6,7))/3))
	hlslhh<-left_join(hlslhh,rtp,by=c("spp","area","quarter","sex"))
	hlslhh<-hlslhh%>% mutate(fishw=1000*lenNum*(a*(lenCls/10)^b))%>%
		filter(!is.na(subSampWt))
	slcheckw<-hlslhh%>%group_by(sampType,landCtry,vslFlgCtry,year,proj,
			      trpCode,staNum,spp,catchCat,landCat,
			      commCatScl,commCat,subSampCat,sex,lenCode)%>%
		summarise(subSampWt=unique(subSampWt),wt=unique(wt),totestim=sum(fishw))%>%
		ungroup()%>%
		mutate(diffw=abs((subSampWt-totestim)/totestim),valdiff=abs(subSampWt-totestim))%>%
		mutate(testdiffw=diffw>=0.5,labeltestdiffw=ifelse(testdiffw,"w diff >=50%","w diff <50%"))

	#1. if subsampwt is NA invalid the sample
	slnow<-CSr@sl%>%filter(is.na(wt)|is.na(subSampWt))%>%mutate(comment="NA weight")

	#2. weight diff between data and rtp weight > 50% : plot and extract
	#info
	p1<-ggplot(slcheckw%>%filter(!is.na(labeltestdiffw))%>%mutate(year=sub("20","",year)),
		   aes(x=subSampWt/1000,y=totestim/1000,color=sex,size=valdiff/1000))+
			geom_point(alpha=.5)+ #geom_smooth(method="lm",se=F) +
			#geom_point(data=slcheckw,aes(x=subSampWt/1000,y=totestim/1000,shape=sex,size=diffw,color=sampType,alpha=.5))+
			facet_grid(year~spp+sampType+labeltestdiffw,scales="free")+
			xlab("Measured weight (kg)")+ ylab("Weight estimates using length-weight relationship (kg)")+
			theme(legend.position="bottom")+
			ggtitle("Samples weights observed and estimated")


	#3. observation ratio
	slobs50<-CSr@sl%>%mutate(ratio=wt/subSampWt)%>%mutate(ratio=ifelse(is.na(ratio),1,ratio))%>%filter(ratio>=50 & sampType=="S")%>%
		mutate(comment="observation ratio (wt/subSampWt) exceed 50 in sl")


	#3. gather info and save a file with the data problem
	slw50<-slcheckw%>%filter(diffw>=0.5)%>%mutate(comment="weight difference between observation and rtp calculated one > 50%")
	slwrong<-rbind(slnow,slw50[,names(slnow)],slobs50[,names(slnow)])

	#4. invalid trpCode staNum
	idwrong<-unique(paste(slwrong$trpCode,slwrong$staNum))
	CSr@hh<-CSr@hh%>%mutate(id=paste(trpCode,staNum))%>%mutate(foVal=ifelse(id%in%idwrong,"I",foVal))%>%select(-id)

	#remove invalid samples
	CSr<-subset(CSr,foVal=="V",table="hh",link=FALSE)
	#return graph and co
	write.table(slwrong,file="./ices/idinvalidsih.csv",sep=",",row.names=F,col.names=T)
	return(list(CSr,p1))

}

#correction for age and weight stuff and duplication of the age length key
corrcaw<-function(CSr,rtp){
	#test area
	if(F){
		require(ggplot2);require(pander)
		require(dplyr)
		load("./data/CSrall.rdata")
		load("./data/CLrall.rdata")
		y3<-2019:2009
		CLr<-subset(CLr,year%in%y3,table="cl")
		CSr<-subset(CSr,year%in%y3,table="hh")
		CSr<-corrbase(CSr)
		rtp<-findrtp(CLr,CSr)
		CSr<-corrsampw(CSr,rtp)
		CSr<-CSr[[1]]

		pipo<-corrcaw(CSr,rtp)


	}
	#fin test
	#length NA removed
	if(is.numeric(CSr@ca$lenCls)){
		CSr@ca<-CSr@ca[!is.na(CSr@ca$lenCls),]
		CSr@ca<-CSr@ca[CSr@ca$lenCls>0,]
	}else{
		stop("lenCls in ca table is not a numeric !!!")
	}
	#remove imprecise location (aka with ",")
	#pipo<-CSr@ca%>%filter(!grepl(",",area))
	#############################################################################
	#duplicate age info over division where sample area if possible
	divhh<-unique(CSr@hh$area)
	divhh<-divhh[!grepl(",",divhh)]
	CSrtmp<-CSr
	CSrtmp@ca<-csData()@ca[-1,]
	CSrtmp@tr<-CSrtmp@tr%>%filter(!trpCode%in%unique(CSr@ca$trpCode))
	#duplicate caage over divhh
	caage<-CSr@ca
	caagedup<-caage[rep(1:nrow(caage),length(divhh)),]
	divdup<-rep(divhh,each=nrow(caage))
	divdupnb<-rep(1:length(divhh),each=nrow(caage))
	caagedup$area<-divdup
	caagedup$trpCode<-paste0(divdup,caagedup$trpCode)
	caagedup$fishId<-as.numeric(paste0(divdupnb,caagedup$fishId))
	trcaagedup<-caagedup%>%select(sampType:trpCode)%>%distinct()%>%
		mutate(Vessel_length=NA,Vessel_power=NA,Vessel_size=NA,
		       Vessel_type=NA,Harbour=NA,No_SetsHauls_on_Trip=NA,
		       Days_at_sea=NA,Vessel_identifier=NA,sampCtry=NA,
		       Sampling_method="Observer")
	newca<-csData(tr=data.frame(trcaagedup),ca=data.frame(caagedup),check=F)
	#old vs new graph
	pltoldvsnew<-ggplot(CSr@ca%>%mutate(type="old"),
			    aes(x=age,y=lenCls,shape=area,color=paste(sex,area)))+
		geom_point(alpha=.5)+
		geom_point(data=newca@ca%>%mutate(type="new"),
			   aes(x=age,y=lenCls,shape=area,color=paste(sex,area)),alpha=.5)+
		facet_grid(year~quarter+type)
	#add duplicated age to CSr
	CSr<-rbind2(CSrtmp,newca)

	
	#add weight using rtp if needed
	pipo<-left_join(CSr@ca%>%mutate(sex=ifelse(sex%in%c("M","F"),sex,""),type="raw"),
				  rtp,by=c("spp","area","quarter","sex"))
	testnow<-!(pipo$indWt!=-1 & !is.na(pipo$indWt))
	pipo$type[testnow]<-"no data: rtp w"
	#pipo$newWt<-1000*pipo$a*(pipo$lenCls/10)^pipo$b
	pipo$indWt[testnow]<- 1000*pipo$a[testnow]*(pipo$lenCls[testnow]/10)^pipo$b[testnow]
	pipo$rtpw<- 1000*pipo$a*(pipo$lenCls/10)^pipo$b
 	#detect outliers and replace them with rtp weight
  	#checkout<-robustbase::adjOutlyingness(pipo%>%select(indWt,rtpw),coef=10,clower=0,cupper=0)
	#pipo$type[!checkout$nonOut]<-"outliers"
	pipo$check50<- (abs(pipo$indWt-pipo$rtpw)/pipo$rtpw)>0.5
	pipo$type[pipo$check50]<-"outliers 50%"
	msg<-paste0("n:",nrow(pipo),"/nb NA:",length(testnow[testnow]),
		    "/nb outliers:",nrow(pipo[pipo$check50,]))
	trucseg<-pipo[pipo$check50,]%>%mutate(x1=indWt,x2=rtpw,y1=rtpw,y2=rtpw)

	 #now correct outliers w with rtp
	 pipo$indWt[pipo$check50]<-pipo$rtpw[pipo$check50]

	#produce a general plot
	p0<-ggplot(pipo%>%mutate(quarter=as.character(quarter)),
		   aes(x=indWt,y=rtpw,color=paste(sex,type)))+#,size=lenCls/100))+
		geom_point(alpha=.4)+
		geom_curve(data=trucseg,aes(x=indWt,xend=rtpw,y=rtpw,yend=rtpw),
			     arrow = arrow(length = unit(0.01, "npc")),
			     alpha=.5)+
	   geom_point(data=trucseg,aes(x=indWt,y=rtpw),color="black",pch="+")+
	   geom_point(data=trucseg,aes(x=rtpw,y=rtpw),color="black",pch="+")+
		ggtitle(paste0("Weigths outliers detection and correction\n",msg))
		#geom_smooth(method="loess",se=F)+
		#facet_grid(year~area,scale="free")


	p1<-ggplot(pipo%>%mutate(quarter=as.character(quarter)),
		   aes(y=lenCls,x=indWt,color=paste(sex,quarter)))+
		geom_point(alpha=.7) +
		#geom_smooth(method="loess",se=F)+
		geom_point(data=pipo[testnow,],aes(y=lenCls,x=indWt),
			   color="red",shape="+",size=5) +
		facet_grid(year~area,scale="free")

	p2<-ggplot(pipo%>%mutate(quarter=as.character(quarter)),
		   aes(y=lenCls,x=age,color=paste(sex,quarter)))+
		geom_point(alpha=.7) +
		#geom_smooth(method="loess",se=F)+
		geom_point(data=pipo[testnow,],aes(y=lenCls,x=age),
			   color="red",shape="+",size=5) +
		facet_grid(year~area,scale="free")
	#rebuild a CSr with the new ca
	CSrtmp<-CSr
	CSrtmp@ca<-csData()@ca[-1,]
	CSrtmp@tr<-CSrtmp@tr%>%filter(!trpCode%in%unique(CSr@ca$trpCode))
	trcanew<-pipo[,names(csData()@ca)]%>%select(sampType:trpCode)%>%distinct()%>%
		mutate(Vessel_length=NA,Vessel_power=NA,Vessel_size=NA,
		       Vessel_type=NA,Harbour=NA,No_SetsHauls_on_Trip=NA,
		       Days_at_sea=NA,Vessel_identifier=NA,sampCtry=NA,
		       Sampling_method="Observer")
	newca<-csData(tr=data.frame(trcanew),ca=data.frame(pipo[,names(csData()@ca)]),check=F)
	CSr<-rbind2(CSrtmp,newca)

	#return(list(pipo[,names(CSr@ca)],p1,p2,p0,pltoldvsnew))
	return(list(CSr,p1,p2,p0,pltoldvsnew))
}
#get old estimates from old directories for year n-1
getoldices<-function(param,nbef=0,info,HIfinal,SIfinal,SDfinal,SDagefinal,wgname=""){
	if(wgname!=""){ info$stock<-wgname }
 	tmpoldices<- paste0("../../../analyses_stock_",param$currentyear-nbef,"/",info$wg,"/",info$stock,"/ices/")
        oldlen<-dir(path=tmpoldices,patt="length.csv",full=T)
        oldlenshort<-dir(path=tmpoldices,patt="length.csv",full=F)
	oldlen<-oldlen[grepl(paste0(tolower(info$taxon),"\\."),oldlenshort)]
        oldage<-dir(path=tmpoldices,patt="age.csv",full=T)
        oldageshort<-dir(path=tmpoldices,patt="age.csv",full=F)
	oldage<-oldage[grepl(paste0(tolower(info$taxon),"\\."),oldageshort)]
	if(length(oldlen)==0){
        	oldlen<-dir(path=tmpoldices,patt=".csv",full=T)
	}
        HIold<-HIfinal[0,]
        SIold<-SIfinal[0,]
        SDold<-SDfinal[0,]
        SDageold<-SDagefinal[0,]
        for(listfich in oldlen){
                pipo<-readLines(listfich)
                idnumHI<-names(Filter(is.numeric,HIfinal))
                HI<-data.frame(HI=pipo[grepl("HI,",pipo)])%>%
                        separate(HI,into=names(HIfinal),sep=",")%>%
                        mutate_at(idnumHI,as.numeric)
                HIold<-rbind(HIold,HI)
                idnumSI<-names(Filter(is.numeric,SIfinal))
                SI<-data.frame(SI=pipo[grepl("SI,",pipo)])%>%
                        separate(SI,into=names(SIfinal),sep=",")%>%
                        mutate_at(idnumSI,as.numeric)%>%
			mutate(CATON=ifelse(UnitCATON=="t",CATON*1000,CATON))
                SIold<-rbind(SIold,SI)
                idnumSD<-names(Filter(is.numeric,SDfinal))
                SD<-data.frame(SD=pipo[grepl("SD,",pipo)])%>%
                        separate(SD,into=names(SDfinal),sep=",")%>%
                        mutate_at(idnumSD,as.numeric)
                SDold<-rbind(SDold,SD)
        }
        for(listfich in oldage){
                pipo<-readLines(listfich)
                idnumSD<-names(Filter(is.numeric,SDagefinal))
                SDage<-data.frame(SD=pipo[grepl("SD,",pipo)])%>%
                        separate(SD,into=names(SDfinal),sep=",")%>%
                        mutate_at(idnumSD,as.numeric)
                SDageold<-rbind(SDageold,SDage)
        }
	return(list(HI=HIold,SI=SIold,SD=SDold,SDage=SDageold))


}

anaoldices<-function(SIfinal,oldestim1,oldestim2,oldestim3){

	stock0<-SIfinal%>%
		group_by(Year=as.numeric(Year),FishingArea,CatchCategory,id=param$currentyear+1)%>%
		summarise(CATON=sum(CATON,na.rm=T))%>%ungroup()
	stock1<-oldestim1$SI %>%group_by(Year=as.numeric(Year),FishingArea,CatchCategory,id=param$currentyear)%>%
		summarise(CATON=sum(CATON,na.rm=T))%>%ungroup()
	stock2<-oldestim2$SI %>%group_by(Year=as.numeric(Year),FishingArea,CatchCategory,id=param$currentyear-1)%>%
		summarise(CATON=sum(CATON,na.rm=T))%>%ungroup()
	stock3<-oldestim3$SI %>%group_by(Year=as.numeric(Year),FishingArea,CatchCategory,id=param$currentyear-2)%>%
		summarise(CATON=sum(CATON,na.rm=T))%>%ungroup()
	allversion<-rbind(stock0,stock1,stock2,stock3)%>%
		mutate(Year=as.numeric(Year),id=paste("estim",as.character(id)))%>%
		filter(Year>=param$currentyear-3)
	plt<-ggplot(data=allversion,aes(x=Year,y=CATON,fill=id))+
		geom_bar(alpha=0.7,stat="identity",position=position_dodge())+
		facet_grid(CatchCategory~FishingArea,scale="free")+
			theme(legend.position="bottom",
			axis.text.x  = element_text(angle=90, vjust=0.5, size=6),
			axis.text.y  = element_text(angle=0, vjust=0, size=6),
		      	strip.text.x = element_text(size=6, angle=0))

	tab0<-allversion%>%group_by(Year,id,CatchCategory)%>%summarise(w=sum(CATON))%>%ungroup()
	alltruc<-expand.grid((param$currentyear-3):param$currentyear,
			     paste("estim",(param$currentyear-2):(param$currentyear+1)),c("D","L"))%>%
		transmute(Year=Var1,id=Var2,CatchCategory=Var3)
	tab0<-left_join(alltruc,tab0)

	#landings
	tab1<-tab0%>%filter(CatchCategory=="L")%>%pivot_wider(values_from=w,names_from=id)
	#compute diff with last year using diagonal of numerical value of tab1
	diagL<-tab1[,grepl("estim",names(tab1))]%>%as.matrix()%>%diag()
	diagL[is.na(diagL)]<-0
	lastL<-tab1%>%pull(ncol(tab1))
	lastL[is.na(lastL)]<-0
	lastL<-lastL+1e-9
	diffL<-data.frame(V1=round(100*(lastL-diagL)/lastL,2))
	names(diffL)<-paste0("rounded diff in %\nwith ",param$currentyear+1)
	tab1<-cbind(tab1,diffL)
	names(tab1)<-sub("estim ","est",names(tab1))
	names(tab1)<-sub("CatchCategory","CatchCat",names(tab1))
	tab1bis<-t(tab1[,-1])
	colnames(tab1bis)<-tab1[,1]


	#discards
	tab2<-tab0%>%filter(CatchCategory=="D")%>%pivot_wider(values_from=w,names_from=id)
	#compute diff with last year using diagonal of numerical value of tab1
	diagD<-tab2[,grepl("estim",names(tab2))]%>%as.matrix()%>%diag()
	diagD[is.na(diagD)]<-0
	lastD<-tab2%>%pull(ncol(tab2))
	lastD[is.na(lastD)]<-0
	lastD<-lastD+1e-9
	diffD<-data.frame(V1=round(100*(lastD-diagD)/lastD,2))
	names(diffD)<-paste0("rounded diff in %\nwith ",param$currentyear+1)
	tab2<-cbind(tab2,diffD)
	names(tab2)<-sub("estim ","est",names(tab2))
	names(tab2)<-sub("CatchCategory","CatchCat",names(tab2))
	tab2bis<-t(tab2[,-1])
	colnames(tab2bis)<-tab2[,1]

	return(list(plt,tab1bis,tab2bis))
}


grapholdnewlen<-function(SDfinal,oldestim1,oldestim2,oldestim3){

	stock0<-SDfinal%>%
		group_by(Year=as.numeric(Year),FishingArea,CatchCategory,AgeLength,id=param$currentyear+1)%>%
		summarise(NumberCaught=sum(NumberCaught,na.rm=T))%>%ungroup()
	stock1<-oldestim1$SD %>%
		group_by(Year=as.numeric(Year),FishingArea,CatchCategory,AgeLength,id=param$currentyear)%>%
		summarise(NumberCaught=sum(NumberCaught,na.rm=T))%>%ungroup()
	stock2<-oldestim2$SD %>%
		group_by(Year=as.numeric(Year),FishingArea,CatchCategory,AgeLength,id=param$currentyear-1)%>%
		summarise(NumberCaught=sum(NumberCaught,na.rm=T))%>%ungroup()
	stock3<-oldestim3$SD %>%
		group_by(Year=as.numeric(Year),FishingArea,CatchCategory,AgeLength,id=param$currentyear-2)%>%
		summarise(NumberCaught=sum(NumberCaught,na.rm=T))%>%ungroup()

	allversion<-rbind(stock0,stock1,stock2,stock3)%>%
		mutate(id=paste("estim",as.character(id)),Year=as.character(Year),AgeLength=as.numeric(AgeLength))%>%
		filter(Year>=param$currentyear-3)

	if(nrow(allversion)>0){
	plt<-ggplot(data=allversion,aes(x=AgeLength,y=NumberCaught,group=id,color=id))+
		geom_line(alpha=0.7)+#,stat="identity",position=position_dodge())+
		facet_grid(Year+CatchCategory~FishingArea,scale="free")+
			theme(legend.position="bottom",
			axis.text.x  = element_text(angle=90, vjust=0.5, size=6),
			axis.text.y  = element_text(angle=0, vjust=0, size=6),
		      	strip.text.x = element_text(size=6, angle=0))
	}else{
		plt<-ggplot()
	}
	return(plt)
}

grapholdnewage<-function(SDagefinal,oldestim1,oldestim2,oldestim3){

	stock0<-SDagefinal%>%
		group_by(Year=as.numeric(Year),FishingArea,CatchCategory,AgeLength,id=param$currentyear+1)%>%
		summarise(NumberCaught=sum(NumberCaught,na.rm=T))%>%ungroup()
	stock1<-oldestim1$SDage %>%
		group_by(Year=as.numeric(Year),FishingArea,CatchCategory,AgeLength,id=param$currentyear)%>%
		summarise(NumberCaught=sum(NumberCaught,na.rm=T))%>%ungroup()
	stock2<-oldestim2$SDage %>%
		group_by(Year=as.numeric(Year),FishingArea,CatchCategory,AgeLength,id=param$currentyear-1)%>%
		summarise(NumberCaught=sum(NumberCaught,na.rm=T))%>%ungroup()
	stock3<-oldestim3$SDage %>%
		group_by(Year=as.numeric(Year),FishingArea,CatchCategory,AgeLength,id=param$currentyear-2)%>%
		summarise(NumberCaught=sum(NumberCaught,na.rm=T))%>%ungroup()

	allversion<-rbind(stock0,stock1,stock2,stock3)%>%
		mutate(id=paste("estim",as.character(id)),Year=as.character(Year),AgeLength=as.numeric(AgeLength))%>%
		filter(Year>=param$currentyear-3)

	if(nrow(allversion)>0){
	plt<-ggplot(data=allversion,aes(x=AgeLength,y=NumberCaught,group=id,color=id))+
		geom_line(alpha=0.7)+#,stat="identity",position=position_dodge())+
		facet_grid(Year+CatchCategory~FishingArea)+
			theme(legend.position="bottom",
			axis.text.x  = element_text(angle=90, vjust=0.5, size=6),
			axis.text.y  = element_text(angle=0, vjust=0, size=6),
		      	strip.text.x = element_text(size=6, angle=0))
	}else{
		plt<-ggplot()
	}
	return(plt)
}
#compare old CL/CE with the new one
 getoldcl<-function(CLr,info,param,wgname=""){
	 if(F){
		 library(ggplot2);library(dplyr)
		 load("test.rdata")
	 }
	if(wgname!=""){ info$stock<-wgname }
         #landings
         s1<-CLr@cl%>%group_by(year,area)%>%summarise(w=sum(landWt))%>%ungroup()%>%
                 mutate(type=paste0("raw data ",param$currentyear+1))
         oldfich<-paste0("../../../analyses_stock_",param$currentyear-0,"/",info$wg,"/",info$stock,"/data/CLrall.rdata")
	 if(file.exists(oldfich)){ load(oldfich) }else{ CLr<-clData() }
         s2<-CLr@cl%>%group_by(year,area)%>%summarise(w=sum(landWt))%>%ungroup()%>%
                 mutate(type=paste0("raw data ",param$currentyear-0))
         oldfich<-paste0("../../../analyses_stock_",param$currentyear-1,"/",info$wg,"/",info$stock,"/data/CLrall.rdata")
	 if(file.exists(oldfich)){ load(oldfich) }else{ CLr<-clData() }
         s3<-CLr@cl%>%group_by(year,area)%>%summarise(w=sum(landWt))%>%ungroup()%>%
                 mutate(type=paste0("raw data ",param$currentyear-1))
         oldfich<-paste0("../../../analyses_stock_",param$currentyear-2,"/",info$wg,"/",info$stock,"/data/CLrall.rdata")
	 if(file.exists(oldfich)){ load(oldfich) }else{ CLr<-clData() }
         s4<-CLr@cl%>%group_by(year,area)%>%summarise(w=sum(landWt))%>%ungroup()%>%
                 mutate(type=paste0("raw data ",param$currentyear-2))
         s123lan<-rbind(s1,s2,s3,s4)%>%transmute(year,area,name="landings",value=w,type)
         #effort
        s1<-CEr@ce%>%group_by(year,area)%>%
         summarise_at(vars(trpNum,foNum,foDur,daysAtSea,effKwDays),sum)%>%ungroup()%>%
         pivot_longer(trpNum:effKwDays)%>%
                mutate(type=paste0("raw data ",param$currentyear+1))
         oldfich<-paste0("../../../analyses_stock_",param$currentyear-0,"/",info$wg,"/",info$stock,"/data/CErall.rdata")
	 if(file.exists(oldfich)){ load(oldfich) }else{ CEr<-ceData() }
        s2<-CEr@ce%>%group_by(year,area)%>%
         summarise_at(vars(trpNum,foNum,foDur,daysAtSea,effKwDays),sum)%>%ungroup()%>%
         pivot_longer(trpNum:effKwDays)%>%
                mutate(type=paste0("raw data ",param$currentyear-0))
         oldfich<-paste0("../../../analyses_stock_",param$currentyear-1,"/",info$wg,"/",info$stock,"/data/CErall.rdata")
	 if(file.exists(oldfich)){ load(oldfich) }else{ CEr<-ceData() }
        s3<-CEr@ce%>%group_by(year,area)%>%
         summarise_at(vars(trpNum,foNum,foDur,daysAtSea,effKwDays),sum)%>%ungroup()%>%
         pivot_longer(trpNum:effKwDays)%>%
                mutate(type=paste0("raw data ",param$currentyear-1))
         oldfich<-paste0("../../../analyses_stock_",param$currentyear-2,"/",info$wg,"/",info$stock,"/data/CErall.rdata")
	 if(file.exists(oldfich)){ load(oldfich) }else{ CEr<-ceData() }
        s4<-CEr@ce%>%group_by(year,area)%>%
         summarise_at(vars(trpNum,foNum,foDur,daysAtSea,effKwDays),sum)%>%ungroup()%>%
         pivot_longer(trpNum:effKwDays)%>%
                mutate(type=paste0("raw data ",param$currentyear-2))

 s12<-rbind(rbind(s1,s2,s3,s4)%>%filter(!name%in%c("foDur","foNum")),s123lan)
        return(s12)
 }






findprecatch<-function(CLr,CSr){
	if(F){
		library(dplyr)
		load("./data/CLrall.rdata")
		load("./data/CSrall.rdata")
	}
	if(file.exists("../../precatch.csv")){
		precatch<-read.csv("../../precatch.csv")
		precatch$catch<-apply(precatch[,6:7],1,sum,na.rm=T)
		precatch<-precatch%>%transmute(year=Year,
					       spp=Species.Latin.Name,
					       area=tolower(gsub("_",".",Area)),
					       country=Country,
					       lan=catch)%>%filter(country=="FR")
		#clean area precatch
		unique(precatch$area)
		precatch$area[substr(precatch$area,1,7)=="27.10.a"]<-"27.10.a"
		precatch$area[substr(precatch$area,1,7)=="27.12.a"]<-"27.12.a"
		precatch$area[substr(precatch$area,1,7)=="27.14.a"]<-"27.14.a"
		precatch$area[substr(precatch$area,1,7)=="27.14.b"]<-"27.14.b"
		precatch$area[substr(precatch$area,1,6)=="27.2.a"]<-"27.2.a"
		precatch$area[substr(precatch$area,1,6)=="27.2.b"]<-"27.2.b"
		precatch$area[substr(precatch$area,1,6)=="27.3.a"]<-"27.3.a"
		precatch$area[substr(precatch$area,1,6)=="27.3.b"]<-"27.3.b"
		precatch$area[substr(precatch$area,1,6)=="27.3.c"]<-"27.3.c"
		precatch$area[substr(precatch$area,1,6)=="27.3.d"]<-"27.3.d"
		precatch$area[substr(precatch$area,1,6)=="27.5.a"]<-"27.5.a"
		precatch$area[substr(precatch$area,1,6)=="27.5.b"]<-"27.5.b"
		precatch$area[substr(precatch$area,1,6)=="27.6.b"]<-"27.6.b"
		precatch$area[substr(precatch$area,1,6)=="27.7.b"]<-"27.7.b"
		precatch$area[substr(precatch$area,1,6)=="27.7.c"]<-"27.7.c"
		precatch$area[substr(precatch$area,1,6)=="27.7.k"]<-"27.7.k"
		precatch$area[substr(precatch$area,1,6)=="27.7.j"]<-"27.7.j"
		precatch$area[substr(precatch$area,1,6)=="27.8.d"]<-"27.8.d"
		precatch$area[substr(precatch$area,1,6)=="27.8.e"]<-"27.8.e"
		precatch$area[substr(precatch$area,1,6)=="27.9.b"]<-"27.9.b"
		#clean spp
		precatch$spp[precatch$spp=="Beryx"]<-"Beryx spp"

		#species id
		spp0<-unique(CSr@sl$spp)
		area1<-unique(unlist(strsplit(paste0(sort(unique(CLr@cl$area)),collapse=","),",")))
		area2<-unique(unlist(strsplit(paste0(sort(unique(CSr@hh$area)),collapse=","),",")))
		area0<-sort(unique(c(area1,area2)))

		#test spp
		#if(any(spp0%in%unique(precatch$spp))){
			pretot<-precatch%>%filter(spp%in%spp0)
		#}else{
		#	pretot<-data.frame()
		#	for(i in 1:length(spp0)){
		#		pretot0<-precatch%>%filter(spp0[i])
		#	}
		#	pretot<-rbind(pretot,pretot0)
		#}
		pretot<-precatch%>% filter(area%in%area0) %>%filter(spp%in%spp0)
		
		if(nrow(pretot)==0){pretot<-0}else{pretot<-sum(pretot$lan,na.rm=T)}
	}else{
		pretot<-0
	}
	return(pretot)

}

#find rtp
findrtp<-function(CLr,CSr){
	if(F){
	load("./data/CLrall.rdata")
	load("./data/CSrall.rdata")
	}
	rtpspp<-read.csv("/home/moi/ifremer/data/refTables/ISIH-19926-taille_poids-2020.txt",
			 sep=";",stringsAsFactors=FALSE)
	#convert space in ICES code
	fct1<-function(a="001A00"){
		if(nchar(a)==6){
			a1<-as.numeric(substr(a,1,3))#001
			a2<-tolower(gsub("0","",substr(a,4,6)))#A
			rez<-paste0("27.",a1,".",a2)
		}else{
			rez<-a
		}
		return(rez)
	}
	pipo<-data.frame(area=unique(rtpspp$SECT_COD))
	pipo$areanew<-apply(pipo,1,fct1)
	pipo1<-pipo%>%filter(substr(areanew,nchar(areanew),nchar(areanew))==".")
	pipo2<-pipo%>%filter(substr(areanew,nchar(areanew),nchar(areanew))!=".")
	#duplicate line if area finish by a dot (to complete a,b,c...)
	npipo1<-nrow(pipo1)
	pipo1<-pipo1[rep(1:npipo1,each=26),]
	pipo1$areanew<-paste0(pipo1$areanew,rep(letters,npipo1))
	pipo<-rbind(pipo1,pipo2)
	rtpspp<-left_join(rtpspp,pipo,by=c("SECT_COD"="area"))

	#convert month info into quarter if available
	rtpspp<-rtpspp%>%mutate(timeint=paste(MOIS_DEB,MOIS_FIN,sep="-"))
	rtpspp$quarter<-NA
	rtpspp$quarter[rtpspp$timeint=="1-3"]<-1
	rtpspp$quarter[rtpspp$timeint=="4-6"]<-2
	rtpspp$quarter[rtpspp$timeint=="7-9"]<-3
	rtpspp$quarter[rtpspp$timeint=="10-12"]<-4
	#spp with quarterly rtp
	rtpspp1<-rtpspp%>%filter(!is.na(quarter))
	#spp with yearly rtp
	rtpspp2<-rtpspp%>%filter(is.na(quarter))
	#then repeat single value 4 times if timeint is 1-12
	rtpspp2<-rtpspp2[rep(1:nrow(rtpspp2),each=4),]%>%mutate(quarter=rep(1:4,nrow(rtpspp2)))
	#final
	rtpspp<-rbind(rtpspp1,rtpspp2)
	#recode SEXE
	rtpspp$sex<-""
	rtpspp$sex[rtpspp$SEXE=="Male"]<-"M"
	rtpspp$sex[rtpspp$SEXE=="Femelle"]<-"F"
	#recode a and b
	rtpspp$a<-rtpspp$RTP_COEF_A_CM
	rtpspp$b<-rtpspp$RTP_COEF_B
	#recode nom
	rtpspp$spp<-rtpspp$NOM_VALIDE

	#now create a rtp dataframe with area and quarter as it should be
	#according to CLr
	spacespp1<-unique(sort(unlist(strsplit(paste(unique(CLr@cl$area),collapse=","),","))))
	spacespp2<-unique(sort(unlist(strsplit(paste(unique(CSr@hh$area),collapse=","),","))))
	spacespp<-unique(sort(c(spacespp1,spacespp2)))
	listspp<-unique(CSr@sl$spp)
	rtpfinal<-data.frame(expand.grid(spp=listspp,area=spacespp,quarter=1:4,sex=c("","F","M"),stringsAsFactors=F))
	#filter ab stuff for nephrops
	#if(unique(CSr@sl$spp)=="Nephrops norvegicus"){
	if("Nephrops norvegicus"%in%unlist(strsplit(CSr@sl$spp,","))){
		rtpspp<-rtpspp[grepl("LENGTH_CARAPACE",rtpspp$TMESURE_COD),]
		rtpspp$a<-10*rtpspp$RTP_COEF_A_CM*10^rtpspp$RTP_COEF_B
	}

	#rtp
	#rtpfinal$spp<-"Dicentrarchus labrax"
	#rtpfinal$spp<-"Mullus surmuletus"
	#rtpfinal$spp<-"Raja clavata"
	#rtpfinal$spp<-"Solea solea"
	pipo<-left_join(rtpfinal,
			rtpspp%>%transmute(spp,area=areanew,quarter,sex,a,b),
			by=c("spp","area","quarter","sex"))%>%
		group_by(spp,area,quarter,sex)%>%summarise(a=median(a,na.rm=T),b=median(b,na.rm=T))%>%
		ungroup()
	#if F or M or no sex empty use the median one 
	pipo<-pipo%>%group_by(spp,area,quarter)%>%mutate(amed=median(a,na.rm=T),bmed=median(b,na.rm=T),
							     a=ifelse(is.na(a),amed,a),
							     b=ifelse(is.na(b),bmed,b))%>%
			ungroup()%>%select(-amed,-bmed)
	#if no data then use median a and b for the species 
	rtpmed<-rtpspp[rtpspp$spp%in%unique(CSr@sl$spp),]
	amed<-median(rtpmed$a,na.rm=T)
	bmed<-median(rtpmed$b,na.rm=T)
	pipo$a[is.na(pipo$a)]<-amed
	pipo$b[is.na(pipo$b)]<-bmed
	#then if no data use generic formulas
	pipo$a[is.na(pipo$a)]<-0.01/1000
	pipo$b[is.na(pipo$b)]<-3
	return(pipo)

	
}

#simple stock definition
defstock<-function(stock,wg,year,CLr,CSr){
	if(F){
		library(dplyr)
		stock<-"alf.27.nea"
		wg<-"WGDEEP"
		load("./data/CLrall.rdata")
		load("./data/CSrall.rdata")
		year<-2017
		defstock(stock,wg,year,CLr,CSr)
	}
	#load fao ref table
	reffao<-read.csv2("/home/moi/ifremer/data/refTables/ISIH-19926-espece_fao-2019.txt",fileEncoding="ISO-8859-1")
	#check info using wgparam and areaciem
	load("/home/moi/ifremer/analyses/areaciem2017.rdata")
	load("/home/moi/ifremer/analyses/wgparam2017.rdata")
	areaciem<-areaciem%>%filter(newstock==stock)
	wgparam<-wgparam%>%filter(newstock==stock)
	#function unit
	if(!grepl(".fu.",tolower(stock))){
		listarea<-paste(gsub("27.","",sort(unique(areaciem$area[areaciem$type=="Div"]))),collapse=",")
	}else{
		listarea<-paste(gsub("27.","",sort(unique(areaciem$area[areaciem$type=="StatRec"]))),collapse=",")
	}
	commonname0<-paste(reffao[reffao$ESPF_COD==wgparam$fao,6:5],collapse="/")
	info0<-data.frame(wg=wgparam$wg,
			  stock=wgparam$newstock,
			  commonname=commonname0,
			  species=wgparam$spp,
			  taxon=wgparam$fao,
			  area=listarea)

	commonname<-paste(reffao[reffao$ESPF_COD==unique(CLr@cl$taxon),6:5],collapse="/")
	wyear<-CLr@cl%>%group_by(year)%>%summarise(w=sum(landWt,na.rm=T)/1000)
	#area with fu stuff
		if(!grepl(".fu.",tolower(stock))){
		 	area=paste(eval(parse(text=paste0("unique(sort(c('",paste(gsub("27.","",gsub(",","','",
							     sort(unique(c(CSr@hh$area,CLr@cl$area))))),
							   collapse="','"),"')))"))),collapse=",")
		}else{
		 	area=paste(eval(parse(text=paste0("unique(sort(c('",paste(gsub("27.","",gsub(",","','",
							     sort(unique(c(CSr@hh$rect,CLr@cl$rect))))),
							   collapse="','"),"')))"))),collapse=",")
		}

	precatch<-findprecatch(CLr,CSr)
	info<-data.frame(wg=wg,
		 stock=stock,
		 commonname,
		 species=paste(sort(unique(CSr@sl$spp)),collapse=","),
		 taxon=paste(sort(unique(CLr@cl$taxon)),collapse=","),
		 area=sort(unique(area)),
		 #area=paste(na.omit(unique(c(CSr@hh$area,CLr@cl$area))),collapse=","),
		 year=paste(gsub("20","",sort(unique(CLr@cl$year))),collapse=","))
	info$precatch<-precatch
	eval(parse(text= paste0("info$catch_t_",year,"<-as.numeric(wyear%>%filter(year==",year,")%>%select(w))")))
	eval(parse(text= paste0("info$catch_t_",year-1,"<-as.numeric(wyear%>%filter(year==",year-1,")%>%select(w))")))
	eval(parse(text= paste0("info$catch_t_",year-2,"<-as.numeric(wyear%>%filter(year==",year-2,")%>%select(w))")))
	info$sacrois<-"3.3.8"
	#check some stuff and replace
	if(info$species=="NA"){info$species<-info0$species}
	if(info$taxon==""){info$taxon<-info0$taxon}
	if(info$area=="NA"){info$area<-info0$area}
	if(grepl("character",info$commonname)){info$commonname<-info0$commonname}

	if(info$year!=""){
		t1<-info[,c("stock","wg","species","taxon","area")]
		t2<-info0[,c("stock","wg","species","taxon","area")]
	#	if(!isTRUE(all.equal(t1,t2))){
	#		print("Stock definition")
	#		print(all.equal(t1,t2))
	#		print(t(rbind(t1,t2)))
	#	}
	}


	return(list(info,info0))

}


#test data availability against parameters
teststock<-function(info,CLr,CSr,param){
	if(F){

		library(dplyr)
		load("./data/CLrall.rdata")
		load("./data/CSrall.rdata")
		load("test.rdata")

	}
	#allspecies vector if more than 1
	allspecies<-unlist(strsplit(info$species,","))

	#data test
	clcur<-CLr@cl%>%filter(year==param$currentyear)
	#eval(parse(text=paste0("cscur<-subset(CSr,year==",param$currentyear,",table='sl')")))
	#cscur<-subset(CSr,year==param$currentyear,table="sl")
	nbech<-CSr@sl%>%filter(spp%in%allspecies)%>%
		select(trpCode,staNum,catchCat,year)%>%distinct()%>%
		group_by(year,catchCat)%>%summarise(nb_sampled_haul=n())%>%ungroup()
	nbfish<-CSr@hl%>%filter(spp%in%unlist(strsplit(info$species,split=",")))%>%
		select(lenNum,catchCat,year)%>%group_by(year,catchCat)%>%
		summarise(nb_fish=sum(lenNum,na.rm=T))%>%ungroup()

	nbCL<-nrow(clcur)
	testCL<-nrow(clcur)>0 & any(!is.na(clcur),na.rm=T)
	#nbECH<-nbech%>%filter(year==param$currentyear&catchCat=="DIS")%>%select(nb_sampled_haul)%>%as.numeric(na.rm=T)
	#nbECH<-ifelse(is.na(nbECH),0,nbECH)
	nbECH<-nbech%>%filter(year==param$currentyear&catchCat=="LAN")%>%select(nb_sampled_haul)%>%as.numeric(na.rm=T)
	nbECH<-ifelse(is.na(nbECH),0,nbECH)
	testECH<-nbECH>=param$seuilsample
	nbDIS<-nbech%>%filter(year==param$currentyear&catchCat=="DIS")%>%select(nb_sampled_haul)%>%as.numeric(na.rm=T)
	nbDIS<-ifelse(is.na(nbDIS),0,nbDIS)
	testDIS<-nbDIS>=param$seuilsampledis
	nbfishLAN<-nbfish%>%filter(year==param$currentyear&catchCat=="LAN")%>%select(nb_fish)%>%as.numeric(na.rm=T)
	nbfishLAN<-ifelse(is.na(nbfishLAN),0,nbfishLAN)
	testfishLAN<-nbfishLAN>=param$seuilfish & testECH
	nbfishDIS<-nbfish%>%filter(year==param$currentyear&catchCat=="DIS")%>%select(nb_fish)%>%as.numeric(na.rm=T)
	nbfishDIS<-ifelse(is.na(nbfishDIS),0,nbfishDIS)
	testfishDIS<-nbfishDIS>=param$seuilfish & testDIS
	nbage<-as.numeric(CSr@ca%>%filter(year==param$currentyear & age!=-1 & spp%in%allspecies)%>%filter(!is.na(age))%>%summarise(nb=n()))
	testage<-nbage>=param$seuilage & (testfishLAN|testfishDIS)
	nbmasse<-as.numeric(CSr@ca%>%filter(year==param$currentyear & indWt!=-1 & spp%in%allspecies)%>%filter(!is.na(age))%>%summarise(nb=n()))
	testmasse<-nbmasse>=param$seuilweight & (testfishLAN|testfishDIS)

	tabtest<-data.frame(Type=c("LAN","LAN samples","DIS samples","sizeLAN","sizeDIS","age","weight"),
		Results=c(testCL,testECH,testDIS,testfishLAN,testfishDIS,testage,testmasse),
		Tresholds=c(0,param$seuilsample,param$seuilsampledis,param$seuilfish,param$seuilfish,param$seuilage,param$seuilweight),
		Nb_currentyear=c(paste(nbCL,"rows in CL"),
		     paste(nbECH,"LAN samples"),
		     paste(nbDIS,"DIS samples"),
		     paste(nbfishLAN,"fishes in LAN"),
		     paste(nbfishDIS,"fishes in DIS"),
		     paste(nbage,"ages in CA"),
		     paste(nbmasse,"weight in CA")))
	#tabtest for a period of interest using param info

	#data test
	listyear<-param$currentyear:(param$currentyear-param$nbyear)
	clcur<-CLr@cl%>%filter(year%in%listyear)
	#eval(parse(text=paste0("cscur<-subset(CSr,year==",param$currentyear,",table='sl')")))
	#cscur<-subset(CSr,year==param$currentyear,table="sl")
	nbech<-CSr@sl%>%filter(spp%in%allspecies)%>%
		select(trpCode,staNum,catchCat,year)%>%distinct()%>%
		group_by(year,catchCat)%>%summarise(nb_sampled_haul=n())%>%ungroup()
	#nbfish<-CSr@hl%>%filter(spp%in%info$species)%>%
	nbfish<-CSr@hl%>%filter(spp%in%unlist(strsplit(info$species,split=",")))%>%
		select(lenNum,catchCat,year)%>%group_by(year,catchCat)%>%
		summarise(nb_fish=sum(lenNum,na.rm=T))%>%ungroup()

	nbCL<-nrow(clcur)
	testCL<-nrow(clcur)>0 & any(!is.na(clcur),na.rm=T)
	#nbECH<-nbech%>%filter(year%in%listyear&catchCat=="DIS")%>%pull(nb_sampled_haul)%>%as.numeric(na.rm=T)%>%sum
	#nbECH<-ifelse(is.na(nbECH),0,nbECH)
	nbECH<-nbech%>%filter(year%in%listyear&catchCat=="LAN")%>%pull(nb_sampled_haul)%>%as.numeric(na.rm=T)%>%sum
	nbECH<-ifelse(is.na(nbECH),0,nbECH)
	testECH<-nbECH>=param$seuilsample
	nbDIS<-nbech%>%filter(year%in%listyear&catchCat=="DIS")%>%pull(nb_sampled_haul)%>%as.numeric(na.rm=T)%>%sum
	nbDIS<-ifelse(is.na(nbDIS),0,nbDIS)
	testDIS<-nbDIS>=param$seuilsampledis
	nbfishLAN<-nbfish%>%filter(year%in%listyear&catchCat=="LAN")%>%pull(nb_fish)%>%as.numeric(na.rm=T)%>%sum
	nbfishLAN<-ifelse(is.na(nbfishLAN),0,nbfishLAN)
	testfishLAN<-nbfishLAN>=param$seuilfish & testECH
	nbfishDIS<-nbfish%>%filter(year%in%listyear&catchCat=="DIS")%>%pull(nb_fish)%>%as.numeric(na.rm=T)%>%sum
	nbfishDIS<-ifelse(is.na(nbfishDIS),0,nbfishDIS)
	testfishDIS<-nbfishDIS>=param$seuilfish & testDIS
	nbage<-as.numeric(CSr@ca%>%filter(year%in%listyear& age!=-1 & spp%in%allspecies)%>%filter(!is.na(age))%>%summarise(nb=n()))
	testage<-nbage>=param$seuilage & (testfishLAN|testfishDIS)
	nbmasse<-as.numeric(CSr@ca%>%filter(year%in%listyear& indWt!=-1 & spp%in%allspecies)%>%filter(!is.na(age))%>%summarise(nb=n()))
	testmasse<-nbmasse>=param$seuilweight & (testfishLAN|testfishDIS)

	tabtest2<-data.frame(Type=c("LAN","LAN samples","DIS samples","sizeLAN","sizeDIS","age","weight"),
		Results=c(testCL,testECH,testDIS,testfishLAN,testfishDIS,testage,testmasse),
		Tresholds=c(0,param$seuilsample,param$seuilsampledis,param$seuilfish,param$seuilfish,param$seuilage,param$seuilweight),
		Nb_allyear=c(paste(nbCL,"rows in CL"),
		     paste(nbECH,"LAN samples"),
		     paste(nbDIS,"DIS samples"),
		     paste(nbfishLAN,"fishes in LAN"),
		     paste(nbfishDIS,"fishes in DIS"),
		     paste(nbage,"ages in CA"),
		     paste(nbmasse,"weight in CA")))

	tabtest$Nb_allyear<-tabtest2$Nb_allyear
	tabtest$Results<-tabtest2$Results
	names(tabtest)[4]<-paste0("Nb_",param$currentyear)
	names(tabtest)[5]<-paste0("Nb_",param$currentyear,"to",param$currentyear-param$nbyear)


	return(tabtest)
}



#test data availability against parameters
teststockold<-function(CLr,CEr,CSr,param){

	#data test
	clcur<-CLr@cl%>%filter(year==param$currentyear)
	cecur<-CEr@ce%>%filter(year==param$currentyear)
	#eval(parse(text=paste0("cscur<-subset(CSr,year==",param$currentyear,",table='sl')")))
	#cscur<-subset(CSr,year==param$currentyear,table="sl")
	nbech<-CSr@sl%>%select(trpCode,staNum,catchCat,year)%>%distinct()%>%
		group_by(year,catchCat)%>%summarise(nb_sampled_haul=n())%>%ungroup()
	nbfish<-CSr@hl%>%select(lenNum,catchCat,year)%>%group_by(year,catchCat)%>%
		summarise(nb_fish=sum(lenNum,na.rm=T))%>%ungroup()

	nbCL<-nrow(clcur)
	testCL<-nrow(clcur)>0 & any(!is.na(clcur),na.rm=T)
	nbECH<-nbech%>%filter(year==param$currentyear&catchCat=="DIS")%>%select(nb_sampled_haul)%>%as.numeric(na.rm=T)
	nbECH<-ifelse(is.na(nbECH),0,nbECH)
	nbECH<-nbECH+nbech%>%filter(year==param$currentyear&catchCat=="LAN")%>%select(nb_sampled_haul)%>%as.numeric(na.rm=T)
	nbECH<-ifelse(is.na(nbECH),0,nbECH)
	testECH<-nbECH>=param$seuilsample
	nbDIS<-nbech%>%filter(year==param$currentyear&catchCat=="DIS")%>%select(nb_sampled_haul)%>%as.numeric(na.rm=T)
	nbDIS<-ifelse(is.na(nbDIS),0,nbDIS)
	testDIS<-nbDIS>=param$seuilsampledis
	nbfishLAN<-nbfish%>%filter(year==param$currentyear&catchCat=="LAN")%>%select(nb_fish)%>%as.numeric(na.rm=T)
	nbfishLAN<-ifelse(is.na(nbfishLAN),0,nbfishLAN)
	testfishLAN<-nbfishLAN>=param$seuilfish
	nbfishDIS<-nbfish%>%filter(year==param$currentyear&catchCat=="DIS")%>%select(nb_fish)%>%as.numeric(na.rm=T)
	nbfishDIS<-ifelse(is.na(nbfishDIS),0,nbfishDIS)
	testfishDIS<-nbfishDIS>=param$seuilfish
	nbage<-as.numeric(CSr@ca%>%filter(year==param$currentyear & age!=-1)%>%filter(!is.na(age))%>%summarise(nb=n()))
	testage<-nbage>=param$seuilage
	nbmasse<-as.numeric(CSr@ca%>%filter(year==param$currentyear & indWt!=-1)%>%filter(!is.na(age))%>%summarise(nb=n()))
	testmasse<-nbmasse>=param$seuilweight

	tabtest<-data.frame(Type=c("LAN","samples","DIS","sizeLAN","sizeDIS","age","weight"),
		Results=c(testCL,testECH,testDIS,testfishLAN,testfishDIS,testage,testmasse),
		Tresholds=c(0,param$seuilsample,param$seuilsampledis,param$seuilfish,param$seuilfish,param$seuilage,param$seuilweight),
		Nb=c(paste(nbCL,"rows in CL"),
		     paste(nbECH,"samples"),
		     paste(nbDIS,"DIS samples"),
		     paste(nbfishLAN,"fishes in LAN"),
		     paste(nbfishDIS,"fishes in DIS"),
		     paste(nbage,"ages in CA"),
		     paste(nbmasse,"weight in CA")))
	return(tabtest)
}


#plot landings ts and co
tscatches<-function(CLrall,info,param){
	#test
 	ts1<-tbl_df(CLrall)%>%group_by(year,foCatEu6)%>%summarise(w=sum(landWt,na.rm=T)/1000)%>%
			group_by(year)%>%mutate(wyear=sum(w),wperc=w/wyear,test1=wperc>=0.1,test2=wperc>=0.01)%>%
			mutate(metier=ifelse(test2,foCatEu6,"<1%"),
			       condition=ifelse(test1,">10%","<10%"))
	ts1<-ts1%>%mutate(metier=sub("_","",sub("_","-",sub("_","",gsub("_0","",metier)))))
	ts1plt<-ggplot(ts1,aes(x=year,y=w,colour=metier,fill=metier))+
		geom_bar(stat="identity")+
		facet_grid(condition~.,scale="free")+
		xlab("")+ylab("Landings (t)")+
		ggtitle(paste0("French landings for ",info$stock,"\nin ",gsub("27.","",info$area)))
	#par division CIEM 
	ts2<-tbl_df(CLrall)%>%group_by(year,area)%>%summarise(w=sum(landWt,na.rm=T)/1000)%>%
			group_by(year)%>%mutate(wyear=sum(w),wperc=w/wyear,test1=wperc>=0.1)%>%
			mutate(condition=ifelse(test1,">10%","<10%"))
	ts2plt<-ggplot(ts2,aes(x=year,y=w,colour=area,fill=area))+
		geom_bar(stat="identity")+
		facet_grid(condition~.,scale="free")+
		xlab("")+ylab("Landings (t)")+
		ggtitle(paste0("French landings for ",info$stock,"\nin ",gsub("27.","",info$area)))
	return(list(ts1plt,ts2plt))
}

#plot landings ts and co
tseffort<-function(CErall,info,param){
	#test
 	ts1<-tbl_df(CErall)%>%group_by(year,foCatEu6)%>%summarise(fd=sum(daysAtSea,na.rm=T))%>%
			group_by(year)%>%mutate(fdyear=sum(fd),fdperc=fd/fdyear,test1=fdperc>=0.1,test2=fdperc>=0.01)%>%
			mutate(metier=ifelse(test2,foCatEu6,"<1%"),
			       condition=ifelse(test1,">10%","<10%"))
	ts1<-ts1%>%mutate(metier=sub("_","",sub("_","-",sub("_","",gsub("_0","",metier)))))
	ts1plt<-ggplot(ts1,aes(x=year,y=fd,colour=metier,fill=metier))+
		geom_bar(stat="identity")+
		facet_grid(condition~.,scale="free")+
		xlab("")+ylab("Days at sea")+
		ggtitle(paste0("French days at sea","\nin ",gsub("27.","",info$area)))
	#par division CIEM 
	ts2<-tbl_df(CErall)%>%group_by(year,area)%>%summarise(fd=sum(daysAtSea,na.rm=T))%>%
			group_by(year)%>%mutate(fdyear=sum(fd),fdperc=fd/fdyear,test1=fdperc>=0.1)%>%
			mutate(condition=ifelse(test1,">10%","<10%"))
	ts2plt<-ggplot(ts2,aes(x=year,y=fd,colour=area,fill=area))+
		geom_bar(stat="identity")+
		facet_grid(condition~.,scale="free")+
		xlab("")+ylab("Days at sea")+
		ggtitle(paste0("French days at sea","\nin ",gsub("27.","",info$area)))
	return(list(ts1plt,ts2plt))
}
#map lanrect and effrect
mapconsrect<-function(lanrect,effrect,param,info){
	if(F){
		load("test.rdata")
		library(dplyr)
		library(ggplot2)
		library(sf)
	}

	#geography
	load("/home/moi/ifremer/data/refTables/geosihsextant/rect.rdata")
	load("/home/moi/ifremer/data/refTables/geosihsextant/div.rdata")
	#prepa data
	#rect
	listrect<-unique(unique(lanrect$space))
	#listrect<-unique(unique(effrect$space),unique(lanrect$space))
	rect<-rect%>%filter(ICESNAME%in%listrect)%>%
			transmute(rect=ICESNAME,div=F_DIVISION,label=ICESNAME,geometry)
	rectxy<-st_coordinates(st_centroid(rect))
	rect$x<-round(rectxy[,1],3)
	rect$y<-round(rectxy[,2],3)
	#div
	listdiv<-unique(rect$div)#sort(unique(c(CLr@cl$area,CSr@hh$area)))
	div<-div%>%
		transmute(div=F_DIVISION,label=ETIQUETTE,geometry)%>%
		filter(div%in%listdiv)
	divxy<-st_coordinates(st_centroid(div))
	div$x<-divxy[,1]
	div$y<-divxy[,2]
	#rect<-rect%>%filter(F_DIVISION=="27.7.d")
	rx<-range(rect$x)
	ry<-range(rect$y)

	#mapbase<-ggplot()+geom_sf(data=div,fill=NA,colour="black",lwd=2)+
	#	    geom_sf(data=rect,fill=NA,alpha=.2)+
	#		borders("world",xlim=rx,ylim=ry,fill="grey",alpha=.5)+
	#		    geom_sf_label(data=div,aes(label=div))+
	#			coord_sf(xlim=rx,ylim=ry)+
	#			    theme_bw()+ggtitle("")+xlab("Longitude")+ylab("Latitude")
	#lan year
	maplan<-lanrect%>%mutate(year=as.numeric(substr(time,1,4)),rect=as.character(space))%>%
		filter(year>=(param$currentyear-param$nbyear))%>%
		group_by(year,rect,technical)%>%summarise(lan_kg=sum(lan_kg,na.rm=T)/1000)%>%ungroup()
	maplan<-maplan%>%left_join(rect,by=c("rect"))
	pltlanyear<-ggplot()+#geom_sf(data=div,fill=NA,colour="black",lwd=2)+
		geom_raster(data=maplan,aes(x=x,y=y,fill=lan_kg),stat="identity",alpha=1)+
		scale_fill_distiller(palette='Spectral',name="Landings (kg)")+
		geom_sf(data=div,fill=NA,colour="black",lwd=.3,alpha=.5)+
		#geom_sf_text(data=div,aes(label=gsub("\\.","",gsub("27.","",div))),size=3,color="grey")+
		#borders("world",xlim=rx,ylim=ry,fill="grey",colour=NA,alpha=.5)+
		borders("world",fill="grey",colour=NA,alpha=.5)+
		facet_grid(technical~year,drop=FALSE)+
		coord_sf(xlim=rx,ylim=ry)+
		theme_bw()+
		ggtitle(paste0("French landings for ",info$stock))+
		xlab("Longitude")+ylab("Latitude")+
		theme(legend.position="bottom",
		      strip.text.y=element_text(size=6,angle=0))


	#eff year
	mapeff<-effrect%>%mutate(year=as.numeric(substr(time,1,4)),rect=as.character(space))%>%
		filter(year>=(param$currentyear-param$nbyear))%>%
		group_by(year,rect,technical)%>%summarise(GtDays=sum(effGtDays,na.rm=T)/1000)%>%ungroup()
	mapeff<-mapeff%>%left_join(rect,by=c("rect"))
	plteffyear<-ggplot()+#geom_sf(data=div,fill=NA,colour="black",lwd=2)+
		geom_raster(data=mapeff,aes(x=x,y=y,fill=GtDays),stat="identity",alpha=1)+
		scale_fill_distiller(palette='Spectral',name="Effort\n(GtDays)")+
		geom_sf(data=div,fill=NA,colour="black",lwd=.3,alpha=.5)+
		#geom_sf_text(data=div,aes(label=gsub("\\.","",gsub("27.","",div))),size=3,color="grey")+
		#borders("world",xlim=rx,ylim=ry,fill="grey",colour=NA,alpha=.5)+
		borders("world",fill="grey",colour=NA,alpha=.5)+
		facet_grid(technical~year,drop=FALSE)+
		coord_sf(xlim=rx,ylim=ry)+
		theme_bw()+
		ggtitle(paste0("French effort for ",info$stock))+
		xlab("Longitude")+ylab("Latitude")+
		theme(legend.position="bottom",
		      strip.text.y=element_text(size=6,angle=0))

	return(list(pltlanyear=pltlanyear,plteffyear=plteffyear))




}

#map raw landings by rectangle
maplan2<-function(CLr,CSr,info,tabtest,param,nbyear=5){
	if(F){
		load("./data/CLrall.rdata")
		load("./data/CSrall.rdata")
		load("test.rdata")
		nbyear<-9
		library(dplyr)
		library(sf)

	}
	#geography
	load("/home/moi/ifremer/data/refTables/geosihsextant/rect.rdata")
	load("/home/moi/ifremer/data/refTables/geosihsextant/div.rdata")
	#prepa data
	y3<-param$currentyear:(param$currentyear-nbyear)
	listdiv<-sort(unique(c(CLr@cl$area[CLr@cl$year%in%y3],CSr@hh$area[CSr@hh$year%in%y3])))
	div<-div%>%
		transmute(div=F_DIVISION,label=ETIQUETTE,geometry)%>%
		filter(div%in%listdiv)
	divxy<-st_coordinates(st_centroid(div))
	div$x<-divxy[,1]
	div$y<-divxy[,2]
	#rect<-rect%>%filter(F_DIVISION=="27.7.d")
	#listrect<-sort(unique(c(CLr@cl$rect,CSr@hh$rect)))
	listrect<-sort(unique(c(CLr@cl$rect[CLr@cl$year%in%y3],CSr@hh$rect[CSr@hh$year%in%y3])))
	rect<-rect%>%filter(F_DIVISION%in%listdiv)%>%
			transmute(rect=ICESNAME,div=F_DIVISION,label=ICESNAME,geometry)%>%
			filter(rect%in%listrect)
	rectxy<-st_coordinates(st_centroid(rect))
	rect$x<-round(rectxy[,1],3)
	rect$y<-round(rectxy[,2],3)

	rx<-range(rect$x)
	ry<-range(rect$y)

	#mapbase<-ggplot()+geom_sf(data=div,fill=NA,colour="black",lwd=2)+
	#	    geom_sf(data=rect,fill=NA,alpha=.2)+
	#		borders("world",xlim=rx,ylim=ry,fill="grey",alpha=.5)+
	#		    geom_sf_label(data=div,aes(label=div))+
	#			coord_sf(xlim=rx,ylim=ry)+
	#			    theme_bw()+ggtitle("")+xlab("Longitude")+ylab("Latitude")
	#lan year
	maplan<-CLr@cl%>%filter(year>=(param$currentyear-nbyear))%>%
		group_by(year,rect)%>%summarise(w=sum(landWt,na.rm=T)/1000)%>%ungroup()
	maplan<-maplan%>%left_join(rect,by=c("rect"))
	pltlanyear<-ggplot()+#geom_sf(data=div,fill=NA,colour="black",lwd=2)+
		geom_raster(data=maplan,aes(x=x,y=y,fill=w),stat="identity",alpha=1)+
		scale_fill_distiller(palette='Spectral',name="Landings (t)")+
		geom_sf(data=div,fill=NA,colour="black",lwd=.3,alpha=.5)+
		geom_sf_text(data=div,aes(label=gsub("\\.","",gsub("27.","",div))),size=3,color="grey")+
		#borders("world",xlim=rx,ylim=ry,fill="grey",colour=NA,alpha=.5)+
		borders("world",fill="grey",colour=NA,alpha=.5)+
		facet_wrap(~year,ncol=4,drop=FALSE)+
		coord_sf(xlim=rx,ylim=ry)+
		theme_bw()+
		ggtitle(paste0("French landings for ",info$stock,"\n","in ",info$area))+
		xlab("Longitude")+ylab("Latitude")#+theme(legend.position="bottom")

	#lan quarter 
	maplan<-CLr@cl%>%filter(year>=(param$currentyear-nbyear))%>%
		group_by(year,quarter,rect)%>%summarise(w=sum(landWt,na.rm=T)/1000)%>%ungroup()
	maplan<-maplan%>%left_join(rect,by=c("rect"))
	pltlanquarter<-ggplot()+#geom_sf(data=div,fill=NA,colour="black",lwd=2)+
		geom_raster(data=maplan,aes(x=x,y=y,fill=w),stat="identity",alpha=1)+
		scale_fill_distiller(palette='Spectral',name="Landings (t)")+
		geom_sf(data=div,fill=NA,colour="black",lwd=.3,alpha=.5)+
		geom_sf_text(data=div,aes(label=gsub("\\.","",gsub("27.","",div))),size=3,color="grey")+
		#borders("world",xlim=rx,ylim=ry,fill="grey",colour=NA,alpha=.5)+
		borders("world",fill="grey",colour=NA,alpha=.5)+
		facet_grid(year~quarter,drop=FALSE)+
		coord_sf(xlim=rx,ylim=ry)+
		theme_bw()+
		ggtitle(paste0("French landings for ",info$stock,"\n","in ",info$area))+
		xlab("Longitude")+ylab("Latitude")
	#sample
	mapsamp<-left_join(CSr@sl,CSr@hh)%>%filter(year>=(param$currentyear-nbyear))%>%
		mutate(idsamp=paste(trpCode,staNum))%>%
		group_by(year,rect)%>%summarise(n=n_distinct(idsamp))%>%ungroup()
	mapsamp<-mapsamp%>%left_join(rect,by=c("rect"))
	pltsampyear<-ggplot()+#geom_sf(data=div,fill=NA,colour="black",lwd=2)+
		geom_raster(data=mapsamp,aes(x=x,y=y,fill=n),stat="identity",alpha=1)+
		scale_fill_distiller(palette='Spectral',name="Number of sample")+
		geom_sf(data=div,fill=NA,colour="black",lwd=.3,alpha=.5)+
		geom_sf_text(data=div,aes(label=gsub("\\.","",gsub("27.","",div))),size=3,color="grey")+
		#borders("world",xlim=rx,ylim=ry,fill="grey",colour=NA,alpha=.5)+
		borders("world",fill="grey",colour=NA,alpha=.5)+
		facet_wrap(~year,ncol=4,drop=FALSE)+
		coord_sf(xlim=rx,ylim=ry)+
		theme_bw()+
		ggtitle(paste0("Sampling numbers for ",info$stock,"\n","in ",info$area))+
		xlab("Longitude")+ylab("Latitude")

	mapsamp<-left_join(CSr@sl,CSr@hh)%>%filter(year>=(param$currentyear-nbyear))%>%
		mutate(idsamp=paste(trpCode,staNum),quarter=ceiling(as.numeric(substr(date,6,7))/3))%>%
		group_by(year,quarter,rect)%>%summarise(n=n_distinct(idsamp))%>%ungroup()
	mapsamp<-mapsamp%>%left_join(rect,by=c("rect"))
	pltsampquarter<-ggplot()+#geom_sf(data=div,fill=NA,colour="black",lwd=2)+
		geom_raster(data=mapsamp,aes(x=x,y=y,fill=n),stat="identity",alpha=1)+
		scale_fill_distiller(palette='Spectral',name="Number of sample")+
		geom_sf(data=div,fill=NA,colour="black",lwd=.3,alpha=.5)+
		geom_sf_text(data=div,aes(label=gsub("\\.","",gsub("27.","",div))),size=3,color="grey")+
		#borders("world",xlim=rx,ylim=ry,fill="grey",colour=NA,alpha=.5)+
		borders("world",fill="grey",colour=NA,alpha=.5)+
		facet_grid(year~quarter,drop=FALSE)+
		coord_sf(xlim=rx,ylim=ry)+
		theme_bw()+
		ggtitle(paste0("Sampling numbers for ",info$stock,"\n","in ",info$area))+
		xlab("Longitude")+ylab("Latitude")


	return(list(pltlanyear=pltlanyear,pltlanquarter=pltlanquarter,
		    pltsampyear=pltsampyear,pltsampquarter=pltsampquarter))




}

#map raw landings by rectangle
maplan<-function(CLr,CSr,info,tabtest,param,nbyear=5){
	#try ices_area stuff
		#mapdiv<-rgdal::readOGR(dsn = "/home/moi/ifremer/data/refTables/ICES_areas",layer="ICES_Areas_20160601_dense")
		#mapdiv2<-fortify(mapdiv)
		#icesdiv<-left_join(mapdiv2,mapdiv@data%>%mutate(id=as.character(rownames(mapdiv@data))),by=c("id"="id"))
		#icesdiv<-icesdiv%>%mutate(area=paste(Major_FA,SubArea,Division,sep="."))
		#iddiv<-cbind(mapdiv@data,data.frame(x=sp::coordinates(mapdiv)[,1],y=sp::coordinates(mapdiv)[,2]))
		#iddiv<-iddiv%>%mutate(area=paste(Major_FA,SubArea,Division,sep="."))
		#save(icesdiv,iddiv,file="/home/moi/ifremer/data/refTables/ICES_areas/icesdiv.rdata")
	load(file="/home/moi/ifremer/data/refTables/ICES_areas/icesdiv.rdata")
	#summarise data: year/trim
	map1<-CLr%>%filter(year>=(param$currentyear-nbyear))%>%group_by(year,area)%>%summarise(w=sum(landWt,na.rm=T)/1000)%>%ungroup()
	icesdivtmp<-inner_join(icesdiv,map1,by=c("area"="area"))
	iddivtmp<-iddiv%>%filter(area%in%unique(icesdivtmp$area))%>%group_by(area)%>%summarise(x=mean(x),y=mean(y))%>%
		ungroup()%>%mutate(shortarea=gsub("\\.","",gsub("27.","",area)))
	map1<-left_join(map1,iddivtmp)
	map1moy<-map1%>%group_by(area,x,y)%>%summarise(w=mean(w))
	map2<-CLr%>%filter(year>=(param$currentyear-nbyear))%>%group_by(year,rect)%>%summarise(w=sum(landWt,na.rm=T)/1000)%>%ungroup()
	map2$lon<-map2$lat<-map2$lonc<-map2$latc<-NA
	for(i in 1:nrow(map2)){
		if(nchar(map2$rect[i])==4){
		map2$lon[i]<-DATRAS::icesSquare2coord(map2$rect[i],"midpoint")$lon
		map2$lat[i]<-DATRAS::icesSquare2coord(map2$rect[i],"midpoint")$lat
		map2$lonc[i]<-DATRAS::icesSquare2coord(map2$rect[i],"corner")$lon
		map2$latc[i]<-DATRAS::icesSquare2coord(map2$rect[i],"corner")$lat
		}
	}
	rangex<-c(min(map2$lon,na.rm=T)-.5,max(map2$lon,na.rm=T)+.5)
	rangey<-c(min(map2$lat,na.rm=T)-.5,max(map2$lat,na.rm=T)+.5)
	#poly map
		map<-ggplot()+theme_bw()+
			theme(panel.grid.minor.y= element_blank(),
		              panel.grid.minor.x = element_blank())+
			#geom_polygon(data=icesdivtmp,aes(x=long,y=lat,group=group),fill=NA,border="black")+
			geom_path(data=icesdivtmp,aes(x=long,y=lat,group=group),color="dark grey")+#,fill=NA,border="black")+
			geom_raster(data=map2,aes(x=lon,y=lat,fill=w),stat="identity",alpha=.75)+ 
			scale_fill_distiller(palette='Spectral',name="Landings (t)")+
			geom_point(data=map1,aes(x=x,y=y,size=w),alpha=.25)+
			geom_point(data=map1moy,aes(x=x,y=y,size=w),alpha=.5,col="red",fill=NA,shape=1)+
			scale_size_continuous(name="Landings (t, by area)")+
			#geom_polygon(data=coast_map,aes(x=long,y=lat,group=group),fill="grey")+#coord_fixed(1)
			#geom_vline(xintercept=seq(-12, 30, by=1),col="light grey")+
			#geom_hline(yintercept=seq(39,90, by=0.5),col="light grey")+
			borders("world",xlim=rangex,ylim=rangey,fill="light grey",colour="light grey")+
			geom_text(data=iddivtmp,aes(x=x,y=y+1.5,label=shortarea))+
			coord_quickmap(xlim=range(icesdivtmp$long),ylim=range(icesdivtmp$lat))+
			facet_wrap(~year)+xlab("")+ylab("")+
			ggtitle(paste0("French landings for ",info$stock,"\nin ",gsub("27.","",info$area)))

	return(map)

}


#map raw landings by rectangle
maplannico<-function(CLr,CSr,info,tabtest,param){
	#try ices_area stuff
		#mapdiv<-rgdal::readOGR(dsn = "/home/moi/ifremer/data/refTables/ICES_areas",layer="ICES_Areas_20160601_dense")
		#mapdiv2<-fortify(mapdiv)
		#icesdiv<-left_join(mapdiv2,mapdiv@data%>%mutate(id=as.character(rownames(mapdiv@data))),by=c("id"="id"))
		#icesdiv<-icesdiv%>%mutate(area=paste(Major_FA,SubArea,Division,sep="."))
		#iddiv<-cbind(mapdiv@data,data.frame(x=sp::coordinates(mapdiv)[,1],y=sp::coordinates(mapdiv)[,2]))
		#iddiv<-iddiv%>%mutate(area=paste(Major_FA,SubArea,Division,sep="."))
		#save(icesdiv,iddiv,file="/home/moi/ifremer/data/refTables/ICES_areas/icesdiv.rdata")
	if(F){
	load("./data/CSrall.rdata")
	load("./data/CLrall.rdata")
	load("info.rdata")
	load("tabtest.rdata")
	load("param.rdata")
	library(ggplot2)
	library(dplyr)
	library(maps)
	library(mapdata)
	CLr<-CLr@cl%>%filter(harbour=="XSM")

	}
	load(file="/home/moi/ifremer/data/refTables/ICES_areas/icesdiv.rdata")
	#summarise data: year/trim
	map1<-CLr%>%filter(year>=(param$currentyear-5))%>%group_by(year,area,quarter)%>%summarise(w=sum(landWt,na.rm=T)/1000)%>%ungroup()
	icesdivtmp<-inner_join(icesdiv,map1,by=c("area"="area"))
	map2<-CLr%>%filter(year>=(param$currentyear-5))%>%group_by(year,quarter,rect)%>%summarise(w=sum(landWt,na.rm=T)/1000)%>%ungroup()
	map2$lon<-map2$lat<-map2$lonc<-map2$latc<-NA
	for(i in 1:nrow(map2)){
		if(nchar(map2$rect[i])==4){
		map2$lon[i]<-DATRAS::icesSquare2coord(map2$rect[i],"midpoint")$lon
		map2$lat[i]<-DATRAS::icesSquare2coord(map2$rect[i],"midpoint")$lat
		map2$lonc[i]<-DATRAS::icesSquare2coord(map2$rect[i],"corner")$lon
		map2$latc[i]<-DATRAS::icesSquare2coord(map2$rect[i],"corner")$lat
		}
	}
	rangex<-c(min(map2$lon,na.rm=T)-.0,max(map2$lon,na.rm=T)+.0)
	rangey<-c(min(map2$lat,na.rm=T)-.0,max(map2$lat,na.rm=T)+.0)
	rangex<-c(-17.5,0)
	rangey<-c(45,62)

	#poly map
	map<-ggplot()+theme_bw()+
			theme(panel.grid.minor.y= element_blank(),
		              panel.grid.minor.x = element_blank())+
			#geom_polygon(data=icesdivtmp,aes(x=long,y=lat,group=group),fill=NA,border="black")+
			#geom_path(data=icesdivtmp,aes(x=long,y=lat,group=group),color="dark grey")+#,fill=NA,border="black")+
			geom_raster(data=map2,aes(x=lon,y=lat,fill=w),stat="identity",alpha=.75)+ 
			scale_fill_distiller(palette='Spectral',name="Landings (t)")+
			#scale_size_continuous(name="Landings (t, by area)")+
			#geom_polygon(data=coast_map,aes(x=long,y=lat,group=group),fill="grey")+#coord_fixed(1)
			#geom_vline(xintercept=seq(-12, 30, by=1),col="light grey")+
			#geom_hline(yintercept=seq(39,90, by=0.5),col="light grey")+
			borders("worldHires",xlim=rangex,ylim=rangey,fill="light grey",colour="light grey")+
			#geom_text(data=iddivtmp,aes(x=x,y=y+1.5,label=shortarea))+
			#coord_quickmap(xlim=range(icesdivtmp$long),ylim=range(icesdivtmp$lat))+
			coord_quickmap(xlim=rangex,ylim=rangey)+#range(icesdivtmp$long),ylim=range(icesdivtmp$lat))+
			facet_grid(year~quarter)+xlab("")+ylab("")+
			ggtitle(paste0("French landings for ",info$stock))#,"\nin ",gsub("27.","",info$area)))
			map
	ggsave(file="map.pdf")
	ggsave(file="map.png")

	return(map)
	#extract sacrois

}

test<-function(){
print(load("./data/CLrall.rdata"))
print(load("./data/CSrall.rdata"))
library(dplyr)
library(ggplot2)
library(pander)
load("info.rdata")
load("param.rdata")
load("tabtest.rdata")
load("datavalcons.Rdata")
CLc<-CLc@cl
CLv<-CLv@cl
CEc<-CEc@ce
exporticessimple(CLc,CEc,CLv,param,info,excel=TRUE)

}

#make empty HI and SD in case of no data
makeemptyIC<-function(param,info,year){
        load("/home/moi/ifremer/analyses/wgparam2017.rdata")
        aa<-wgparam%>%filter(newstock==info$stock)
        load("/home/moi/ifremer/analyses/areaciem2017.rdata")
        bb<-areaciem%>%filter(newstock==info$stock)%>%filter(type=="Div")
        typearea<-"Div"
        if(nrow(bb)==0){
                        bb<-areaciem%>%filter(newstock==info$stock)%>%filter(type=="SubDiv")
                        listspbb<-data.frame(old=bb$area,new=NA)
                        listspbb$new<-listspbb$old
                        typearea<-"SubDiv"
        }else{
                        #table de correspondance des zones FAO vs ICES
                        listspbb<-data.frame(old=bb$area,new=NA)
                        listspbb$new<-listspbb$old
        }
        HI<-data.frame("HI","FR",year,"Year",year,Fleet="MIS_MIS_0_0_0",typearea,FishingArea=listspbb$new,NA,"kwd",-9,NA)
        names(HI) <- c("RecordType","Country","Year","SeasonType","Season","Fleet","AreaType","FishingArea","DepthRange",
                         "UnitEffort"," Effort","AreaQualifier")
        SI<-data.frame("SI","FR",year,"Year",year,Fleet="MIS_MIS_0_0_0",typearea,FishingArea=listspbb$new,NA,aa$fao,NA,"L","R",NA,"H",NA,NA,"t",0,-9,-9,NA,NA,NA)
        names(SI) <- c("RecordType","Country","Year","SeasonType","Season","Fleet","AreaType","FishingArea","DepthRange","Species","Stock",
                        "CatchCategory","ReportingCategory","DataToFrom","Usage","SamplesOrigin",
                        "QualityFlag","UnitCATON","CATON","OffLandings","varCATON","InfoFleet",
                        "InfoStockCoordinator","InfoGeneral")
	return(list(HI,SI))
}
#export intercatch
exporticesempty<-function(param,info){
        load("/home/moi/ifremer/analyses/wgparam2017.rdata")
        aa<-wgparam%>%filter(newstock==info$stock)
        load("/home/moi/ifremer/analyses/areaciem2017.rdata")
        bb<-areaciem%>%filter(newstock==info$stock)%>%filter(type=="Div")
        typearea<-"Div"
        if(nrow(bb)==0){
                        bb<-areaciem%>%filter(newstock==info$stock)%>%filter(type=="SubDiv")
                        listspbb<-data.frame(old=bb$area,new=NA)
                        listspbb$new<-listspbb$old
                        typearea<-"SubDiv"
        }else{
                        #table de correspondance des zones FAO vs ICES
                        listspbb<-data.frame(old=bb$area,new=NA)
                        listspbb$new<-listspbb$old
        }

        #pander(listspbb[,c("old","new")],style="simple")
        #création des tables HI et SI
        HI<-data.frame("HI","FR",param$currentyear,"Year",param$currentyear,Fleet="MIS_MIS_0_0_0",typearea,FishingArea=listspbb$new,NA,"kwd",-9,NA)
        names(HI) <- c("RecordType","Country","Year","SeasonType","Season","Fleet","AreaType","FishingArea","DepthRange",
                         "UnitEffort"," Effort","AreaQualifier")
        SI<-data.frame("SI","FR",param$currentyear,"Year",param$currentyear,Fleet="MIS_MIS_0_0_0",typearea,FishingArea=listspbb$new,NA,aa$fao,NA,"L","R",NA,"H",NA,NA,"t",0,-9,-9,NA,NA,NA)
        names(SI) <- c("RecordType","Country","Year","SeasonType","Season","Fleet","AreaType","FishingArea","DepthRange","Species","Stock",
                        "CatchCategory","ReportingCategory","DataToFrom","Usage","SamplesOrigin",
                        "QualityFlag","UnitCATON","CATON","OffLandings","varCATON","InfoFleet",
                        "InfoStockCoordinator","InfoGeneral")
        nomfichlen<-paste0("./ices/FRA",param$currentyear,gsub("-","",info$stock),"_length.csv")
        COSTdbe::makeICfile(list(HI=HI,SI=SI), filename=nomfichlen,append=F)
        #export ices  excel too
        wb<-createWorkbook()
        addWorksheet(wb,"SI");writeData(wb,"SI",SI)
        addWorksheet(wb,"HI");writeData(wb,"HI",HI)
        accfile<-paste0("./ices/FRA_",param$currentyear,"_",info$wg,"_",info$stock,".xlsx")
        saveWorkbook(wb,file=accfile,overwrite=T)
        #check
        test<-SI%>%group_by(technical=as.character(Fleet),space=as.character(FishingArea))%>%summarise(totwexport=sum(CATON))
        return(list(HI,SI,test))


}


#export intercatch
exporticessimple<-function(yeartmp,CLc,CEc,CLv,param,info,excel=FALSE){
	CLc<-CLc%>%filter(grepl(yeartmp,time))
	CEc<-CEc%>%filter(grepl(yeartmp,time))
	CLv<-CLv%>%filter(grepl(yeartmp,year))
	CL_met2<-CLc
	CE_met2<-CEc
	CL2<-aggregate(list(Landings=CL_met2$landWt),
	       	       list(Year=substr(CL_met2$time,1,4),Quarter=substr(CL_met2$time,8,8),
			Area=CL_met2$space,Metier=CL_met2$technical),FUN=sum,na.rm=TRUE)
	CL2$Landings<-CL2$Landings/1000 #conversion en tonne
	CE2<-aggregate(list(kwd=CE_met2$effKwDays),
	       	       list(Year=substr(CE_met2$time,1,4),Quarter=substr(CE_met2$time,8,8),
			Area=CE_met2$space,Metier=CE_met2$technical),FUN=sum,na.rm=TRUE)
	CL_CE_met<-merge(CL2,CE2,all.x=TRUE,all.y=FALSE)
	#generate intercatch file
	CatchCategory<-"L";Country<-"FR"
	Species<-unique(CLc$taxon);ReportingCategory = "R"
	#SeasonType<-"Year"
	SeasonType<-param$timestratif
	#add Div or SubArea for AreaTYpe
	checkletter<-function(a){rez<-FALSE;for(i in 1:26){rez<-rez|grepl(letters[i],a)};return(rez)}
	testlett<-checkletter(CL_CE_met$Area)
	AreaType<-ifelse(testlett,"Div","SubArea")
	#Year<-"2014"
	FishingArea = CL_CE_met$Area
	Fleet=CL_CE_met$Metier
	if(nrow(CE2)>0){
		Effort = CL_CE_met$kwd
	}else{
		Effort<--1
	}
	if(SeasonType=="Year"){Season<-CL_CE_met$Year}
	if(SeasonType=="Quarter"){Season = CL_CE_met$Quarter}
	#Season = CL_CE_met$Year
	CATON = CL_CE_met$Landings
	HI <- data.frame(RecordType = rep("HI",nrow(CL_CE_met)), Country=rep(Country,nrow(CL_CE_met)), 
		 		 Year = CL_CE_met$Year, SeasonType = rep(SeasonType,nrow(CL_CE_met)), 
				 Season = Season, Fleet = Fleet, AreaType = AreaType, FishingArea = FishingArea, 
				 DepthRange = NA,  UnitEffort = rep("kwd",nrow(CL_CE_met)), Effort = Effort, 
				AreaQualifier = NA, stringsAsFactors = FALSE)
	SI <- data.frame(RecordType = "SI", Country, Year =CL_CE_met$Year, SeasonType =rep(SeasonType,nrow(CL_CE_met)),
		 		 Season=Season, Fleet, AreaType, FishingArea , DepthRange = NA, Species, Stock = NA, 
				 CatchCategory, ReportingCategory, DataToFrom = "NA", Usage = "H", SamplesOrigin = "NA", 
				QualityFlag = NA, UnitCATON = "t", CATON = CATON, OffLandings = -9, varCATON = -9, InfoFleet = NA, 
				InfoStockCoordinator = NA, InfoGeneral = NA, stringsAsFactors = FALSE)
	COSTdbe::makeICfile(list(HI=HI,SI=SI), filename=paste0("./ices/FRA",yeartmp,gsub("-","",info$stock),".csv"), append=FALSE)

	if(excel){
		require(openxlsx)
		require(tidyr)
		source("./data2icesxls.R")
		#compute area official stuff
		#CLc and CLv have the same order in line
		clrect<-CLv%>%dplyr::select(foCatEu6,rect,area,quarter,landWt)%>%
			mutate(final=CLc$technical,quarter=paste0("q",quarter))%>%filter(!is.na(final))%>%
			group_by(rect,area,quarter)%>%summarise(tot=sum(landWt,na.rm=T)/1000)%>%
				mutate(quarter=factor(quarter,c("q1","q2","q3","q4")))%>%
				ungroup()

				datclrect<-tidyr::spread(clrect,quarter,tot,drop=F)%>%transmute(rect,q1,q2,q3,q4,area)
				dim(datclrect)
				icesrect<-COSTeda::ICESAreaRects%>%transmute(lat,lon,rect=gsub(" ","",as.character(StatRect)),areabis=division,subdiv=subdivision)%>%distinct()
				datclrect<-left_join(datclrect,icesrect)%>%transmute(rect,q1,q2,q3,q4,aa="",lat,lon,area)
				SIrect<-data.frame(datclrect)
				#kg to t
				#save(HI,SI,SIrect,info,file="testexport.Rdata")
				accfile<-paste0("./ices/FRA_",yeartmp,"_",info$wg,"_",info$stock,".xlsx")
				data2icesxls(info=info,param=param,HI=HI,SI=SI,SIrect=SIrect,SDlength="",SDage="",
					     nomfich=accfile)
				#system(paste0("mv ",accfile,"./ices/",accfile))
				#system(paste0("mv FRA_2016_",info$stock,"simple.xlsx ./ices/FRA_2016_",info$stock,"simple.xlsx"))
	}

			#check
			test1<-SI%>%group_by(technical=as.character(Fleet),space=as.character(FishingArea),
					     time=as.character(paste(Year,Season,sep=" - ")))%>%summarise(tot=sum(CATON))
			test2<-CLc%>%group_by(technical=as.character(technical),space=as.character(space),
					      time=as.character(time))%>%summarise(totinit=sum(landWt)/1000)
			if(param$timestratif=="Year"){test2<-test2%>%mutate(time=paste(time,time,sep=" - "))}

			test<-full_join(test1,test2)
			testmet<-test%>%group_by(technical)%>%summarise(tot=sum(tot,na.rm=T),totinit=sum(totinit,na.rm=T))%>%
				transmute(technical,space="all",time="all",tot,totinit)
			test<-rbind(data.frame(test),data.frame(testmet))
			test<-rbind(data.frame(test),
				    aa<-data.frame(technical="all",space="all",time="all",
					tot=sum(SI$CATON),totinit=sum(CLc$landWt)/1000)
				    )
			#test<-test%>%mutate(diff=totwinit-totwexport)
			return(list(test=test,HI=HI,SI=SI))

	
}

#export ices length
exporticesonlydis<-function(year,CEc,dbedis,wdbedis,param,info,rtp,fitvb,checkn=100,checks=1,checksdis=3,wrtp=T,excel=F,age=T){
	if(F){
#		library(COSTdbe)
#		library(dplyr)
#		library(ggplot2)
#		source("credo_fct.R")
#		load('info.rdata')
#		load('param.rdata')
#		load("datavalcons.Rdata")
#		load("myStr.Rdata")
#		load("dbelan.Rdata")
#		dbelan<-dbeObject(species=info$species,taxon=info$taxon,catchCat="LAN",strataDesc=myStr)
#		wdbelan<-dbeObject(species=info$species,taxon=info$taxon,catchCat="LAN",strataDesc=myStr,
#				   methodDesc="analytical",param="weight")
#		dbedis<-dbeObject(species=info$species,taxon=info$taxon,catchCat="DIS",strataDesc=myStr)
#		wdbedis<-dbeObject(species=info$species,taxon=info$taxon,catchCat="DIS",strataDesc=myStr,
#				   methodDesc="analytical",param="weight")
#		load("dbedis.Rdata")
#		load("wdbedis.Rdata")
#		#load("5_agestruc.Rdata")
#		year<-2016
#		rtp<-findrtp(CSc)
#		load("fitvb.rdata")
#		checkn=50
#		checks=3
		library(dplyr);library(ggplot2);library(COSTdbe)
		load("test.rdata")
		year<-2017;checks<-3;checkn<-50;checksdis<-3

	}
	#gloub
	#if annual data add year - year to time if needed
	if(param$timestratif=="Year"){
		#CEc, CLc
		CEc@ce<-CEc@ce%>%mutate(time=paste(time,time,sep=" - "))
		convtimedbe<-function(dbelan){
		for (slot in slotNames(dbelan)){ 
			slottmp<-eval(parse(text=paste0("dbelan@",slot)))
			idslot<-names(slottmp)
			for(i in idslot){
				eval(parse(text=paste0("datatmp<-dbelan@",slot,"$",i)))
				if(class(datatmp)=="data.frame"){
					if(any(names(datatmp)%in%"time")){
						datatmp<-datatmp%>%mutate(time=paste(time,time,sep=" - ")) 
						eval(parse(text=paste0("dbelan@",slot,"$",i,"<-datatmp")))
					}
				}
			}
		} 
		return(dbelan)
		}
		dbedis<-convtimedbe(dbedis)
		wdbedis<-convtimedbe(wdbedis)
	}

	######################################################
	#filter data for a year
	#CSc<-subset(CSc,grepl(year,time),table="hh")
	CEc<-subset(CEc,grepl(year,time),table="ce")
	######################################################
	#ices file generation
	disices<-makeICdf(dbedis,Country="FR",Species=info$taxon,ReportingCategory="R",
				  SamplesOrigin="M",CANUMtype="length",Usage="H",unitMeanWeight="kg")
	######################################################
	#add nbsamp and nbmeas to SD
	disices <-addsampmeas(disices,dbedis)
	#filter
	disices$SD<- disices$SD%>%filter(NumSamplesLngt>=checksdis & NumLngtMeas>=checkn) %>%filter(Year%in%year)
	disices$SI<-disices$SI%>%filter(s>=checksdis) %>%filter(Year%in%year)
	disices$SI<-disices$SI%>%filter(!grepl("MIS_MIS",Fleet)) %>%filter(Year%in%year)
	disices$HI<-disices$HI%>%filter(!grepl("MIS_MIS",Fleet)) %>%filter(Year%in%year)

	######################################################
	#add effort
	disices$HI$UnitEffort<-"kWd"
	effort<-CEc@ce%>%group_by(time,space,technical)%>%summarise(value=sum(effKwDays,na.rm=T))%>%
		ungroup%>%mutate(Year=substr(as.character(time),1,4),Season=substr(as.character(time),8,8))%>%
		transmute(Year,Season,Fleet=technical,FishingArea=space,value)%>%
		ungroup()
	if(param$timestratif=="Year"){effort<-effort%>%mutate(Season=Year)}
	effort<-left_join(disices$HI,effort)
	effort$Effort<-effort$value
	effort<-effort[,names(disices$HI)]
	disices$HI<-effort
	disices$HI<-disices$HI%>%mutate(Effort=ifelse(is.na(Effort),0,Effort))
	disices$SI<-disices$SI%>%mutate(CATON=ifelse(CATON==-9,0,CATON))
	######################################################
	#correction division
	disices$SI<-corrarea(disices$SI)
	disices$HI<-corrarea(disices$HI)
	disices$SD<-corrarea(disices$SD)

	######################################################
	#SD,SI,HI various correction
	disices$SD$NumberCaught<-as.numeric(disices$SD$NumberCaught)
	disices$SD$varNumLanded<-as.numeric(disices$SD$varNumLanded)
	disices$SI$varCATON<--9
	disices$SI$OffLandings<--9
	if(nrow(disices$SD)>0){ disices$SD$CANUMtype<-"Lngt" }
	if(nrow(disices$SD)>0){ disices$SD$varNumLanded<--9 }
	#disices$SI<-disices$SI[disices$SI$CATON>=0,]
	#disices$SD<-disices$SD[disices$SD$NumberCaught>=0,]
	#remove inf and -9 if needed

	disices$SI<-disices$SI[disices$SI$CATON!=Inf,]
	disices$SD<-disices$SD[disices$SD$NumberCaught>0,]

	######################################################
	#weight
	#add meanWeight via wdbelan+completion using rtp
	wfish<-wdbedis@lenStruc$estim%>%transmute(
	Year=sapply(strsplit(as.character(time)," - "),FUN = function(x){x[1]}),
	Season=sapply(strsplit(as.character(time)," - "),FUN = function(x){x[2]}),
						Fleet=as.character(technical),
						AgeLength=as.numeric(as.character(length)),
						FishingArea=space,CatchCategory="D",
					wfish=value)%>%group_by(Year,Season,Fleet,FishingArea,CatchCategory,AgeLength)%>%
		summarise(wfish=mean(wfish,na.rm=T))%>%ungroup()

	#use rtp to complete empty weights
	wfish$wfish[is.na(wfish$wfish)]<-1000*rtp$a*(wfish$AgeLength[is.na(wfish$wfish)]/10)^rtp$b
	wfish<-wfish%>%mutate(AgeLength=as.character(AgeLength))
	wfish<-distinct(wfish)
	#pb de jointure ici
	pipo<-left_join(disices$SD,wfish)%>%mutate(MeanWeight=wfish/1000) %>%select(-wfish)
	if(nrow(pipo)==nrow(disices$SD)){ disices$SD<-pipo}else{stop("merging pb")}
	disices$SD<-left_join(disices$SD,wfish)%>%mutate(MeanWeight=wfish/1000)%>%select(-wfish)
	#complete with rtp if needed
	disices$SD$MeanWeight[is.na(disices$SD$MeanWeight)]<-(1)*rtp$a*(as.numeric(disices$SD$AgeLength[is.na(disices$SD$MeanWeight)])/10)^rtp$b
	if(any(is.na(disices$SD$MeanWeight))){ stop("w pb")}
	######################################################
	#verif et correction SOP
	disices$SD$w<-as.numeric(disices$SD$NumberCaught)*disices$SD$MeanWeight
	verif<-disices$SD%>%group_by(Season,Year,Fleet,FishingArea,CatchCategory)%>%summarise(totw=sum(w,na.rm=T))%>%ungroup()
	verif<-left_join(verif,disices$SI) %>%select(Season,Year,Fleet,FishingArea,CatchCategory,CATON,totw)%>%ungroup()%>%
		mutate(sop=totw/CATON)%>%mutate(sop=ifelse(is.na(sop),1,sop),facteur=ifelse(.8<=sop&sop<=1.2,sop,0.8))
	verif$facteur[verif$sop>1.2]<-1.2
	if(any(is.na(disices$SD$MeanWeight))){ stop("w pb 1")}

	#correction
	disices$SD<-left_join(disices$SD,verif)%>%mutate(MeanWeight=MeanWeight/sop*facteur)
	disices$SD$w<-as.numeric(disices$SD$NumberCaught)*disices$SD$MeanWeight
	verif<-disices$SD%>%group_by(Season,Year,Fleet,FishingArea,CatchCategory)%>%summarise(totw=sum(w,na.rm=T))%>%ungroup()
	verif<-left_join(verif,disices$SI)%>%select(Season,Year,Fleet,FishingArea,CatchCategory,CATON,totw)%>%
		mutate(sop=totw/CATON)
	if(any(is.na(disices$SD$MeanWeight))){ stop("w pb 2")}

	#file generation 
	pipo<-makeICdf(dbedis,Country="FR",Species=info$taxon,ReportingCategory="R",
		       SamplesOrigin="M",CANUMtype="length",Usage="H",unitMeanWeight="kg")
	HI<-disices$HI[,names(pipo$HI)]
	SI<-disices$SI[,names(pipo$SI)]
	SD<-disices$SD[,names(pipo$SD)]
	makeICfile(list(HI=HI,SI=SI,SD=SD), filename=paste0("./ices/FRA",year,gsub("-","",info$stock),"_length.csv"), append=FALSE)
	#on sauve tout
	save(HI,SI,SD,file="HISISD.Rdata")
}
#convert CL to intercatch SI with HI :catchtype : B for BMS, R for logdis
convCLc2ices<-function(CLc,CEc,myStr,param,info,catchtype="B"){
	if(F){
#		checks=3
		library(dplyr);library(ggplot2);library(COSTdbe)
		load("test.rdata")
		load("test2.rdata")
		CLc<-bmsCLc
		myStr<-myStr2
		catchtype<-"B"

	}
	#if annual data add year - year to time if needed
	if(param$timestratif=="Year"){
		#CEc, CLc
		CEc@ce<-CEc@ce%>%mutate(time=paste(time,time,sep=" - "))
		CLc@cl<-CLc@cl%>%mutate(time=paste(time,time,sep=" - "))
	}

	######################################################
	#generate empty dbelan
	dbelan<-dbeObject(species=info$species,taxon=info$species,catchCat='LAN',
			  strataDesc=myStr)
	#add mis_mis and other metier in totalW of dbelan 
	landings<-CLc@cl%>%group_by(time,space,technical)%>%
		summarise(value=sum(landWt,na.rm=T))%>%
		ungroup()
	dbelan@totalW$estim<-landings
	######################################################
	#SI HI SD generation
	lanices<-makeICdf(dbelan,Country="FR",Species=info$taxon,ReportingCategory="R",
				  SamplesOrigin="M",CANUMtype="length",Usage="H",unitMeanWeight="kg")
	lanices$SI$CatchCategory<-catchtype
	######################################################
	#add effort
	lanices$HI$UnitEffort<-"kWd"
	effort<-CEc@ce%>%group_by(time,space,technical)%>%summarise(value=sum(effKwDays,na.rm=T))%>%
		ungroup%>%mutate(Year=substr(as.character(time),1,4),Season=substr(as.character(time),8,8))%>%
		transmute(Year,Season,Fleet=technical,FishingArea=space,value)%>%
		ungroup()
	if(param$timestratif=="Year"){effort<-effort%>%mutate(Season=Year)}
	effort<-merge(lanices$HI,effort,all.x=T)
	effort$Effort<-effort$value
	effort<-effort[,names(lanices$HI)]
	lanices$HI<-effort

	######################################################
	#SD,SI,HI various correction
	lanices$SI<-corrarea(lanices$SI)
	lanices$HI<-corrarea(lanices$HI)
	#put value in OffLandings
	lanices$SI$OffLandings<-lanices$SI$CATON
	#-9 not accepted... so 0
	lanices$SI$CATON<-0
	lanices$SI$varCATON<--9
	#put 0 in HI if no info
	lanices$HI<-lanices$HI%>%mutate(Effort=ifelse(is.na(Effort),0,Effort))
	######################################################
	#filter empty line
	lanices$SI<-lanices$SI[!is.na(lanices$SI$FishingArea),]
	lanices$HI<-lanices$HI[!is.na(lanices$HI$FishingArea),]


	lanices<-list(HI=lanices$HI,SI=lanices$SI)
	return(lanices)

}
#add year in year if stratif is year

		convtimedbe<-function(dbelan){
		for (slot in slotNames(dbelan)){ 
			slottmp<-eval(parse(text=paste0("dbelan@",slot)))
			idslot<-names(slottmp)
			for(i in idslot){
				eval(parse(text=paste0("datatmp<-dbelan@",slot,"$",i)))
				if(class(datatmp)=="data.frame"){
					if(any(names(datatmp)%in%"time")){
						datatmp<-datatmp%>%mutate(time=paste(time,time,sep=" - ")) 
						eval(parse(text=paste0("dbelan@",slot,"$",i,"<-datatmp")))
					}
				}
			}
		} 
		return(dbelan)
		}
#convert COST object in SI, HI, SD stuff
convert2ices<-function(CLc,CEc,dbelan,dbedis,param,info){
	if(F){
#		checks=3
		library(dplyr);library(ggplot2);library(COSTdbe)
		load("test.rdata")

	}
	#gloub
	#if annual data add year - year to time if needed
	if(param$timestratif=="Year"){
		#CEc, CLc
		CEc@ce<-CEc@ce%>%mutate(time=paste(time,time,sep=" - "))
		CLc@cl<-CLc@cl%>%mutate(time=paste(time,time,sep=" - "))
		dbelan<-convtimedbe(dbelan)
		dbedis<-convtimedbe(dbedis)
		wdbelan<-convtimedbe(wdbelan)
		wdbedis<-convtimedbe(wdbedis)
	}

	######################################################
	#add mis_mis and other metier in totalW of dbelan 
	landings<-CLc@cl%>%group_by(time,space,technical)%>%
		summarise(valuereal=sum(landWt,na.rm=T))%>%
		ungroup()
	dbelan@totalW$estim<-merge(dbelan@totalW$estim,landings,all=TRUE)
	dbelan@totalW$estim$value[is.na(dbelan@totalW$estim$value)]<-dbelan@totalW$estim$valuereal[is.na(dbelan@totalW$estim$value)]
	dbelan@totalW$estim<-na.omit(dbelan@totalW$estim)
	######################################################
	#SI HI SD generation
	lanices<-makeICdf(dbelan,Country="FR",Species=info$taxon,ReportingCategory="R",
				  SamplesOrigin="M",CANUMtype="length",Usage="H",unitMeanWeight="kg")
	######################################################
	#add nbsamp and nbmeas to SD
	lanices <-addsampmeas(lanices,dbelan)
	lanices$SD<-lanices$SD%>%filter(is.finite(NumberCaught))%>%filter(NumberCaught>0)

	#filter
	#lanices$SD<-lanices$SD%>%filter(NumSamplesLngt>=checks & NumLngtMeas>=checkn)
	#add discards if any
	disices<-makeICdf(dbedis,Country="FR",Species=info$taxon,ReportingCategory="R",SamplesOrigin="M",
			  	CANUMtype="length",Usage="H",unitMeanWeight="kg")
	disices <-addsampmeas(disices,dbedis)
	#remove Inf and 0 value
	disices$SD<-disices$SD%>%filter(is.finite(NumberCaught))%>%filter(NumberCaught>0)
	disices$SI<-disices$SI%>%filter(is.finite(CATON))%>%filter(CATON>=0)#%>%filter(NumberCaught>0)
	#dbedis@lenStruc$estim%>%filter(technical=="GTR_DEF_120-219_0_0_all" & time =="2016 - 2")
	#dbedis@totalW$estim%>%filter(technical=="GTR_DEF_120-219_0_0_all" & time =="2016 - 2")
	#CSc@hl%>%filter(technical=="GTR_DEF_120-219_0_0_all" & time =="2016 - 2")
	#CSc@sl%>%filter(technical=="GTR_DEF_120-219_0_0_all" & time =="2016 - 2")
	#CSc@hh%>%filter(technical=="GTR_DEF_120-219_0_0_all" & time =="2016 - 2")
	
	
	#remove MIS_MIS in discards
	#disices$SI<-disices$SI%>%filter(!grepl("MIS_MIS",Fleet))
	######################################################
	#bind SI and SI	AND filter only on the considered year
	lanices$SD<-rbind(lanices$SD,disices$SD)%>%#filter(Year%in%year)
		filter(!is.na(Fleet))
	lanices$SI<-rbind(lanices$SI,disices$SI)%>%#filter(Year%in%year)
		filter(!is.na(Fleet))
	lanices$HI<-rbind(lanices$HI,disices$HI)%>%distinct()%>%#filter(Year%in%year)
		filter(!is.na(Fleet))
	######################################################
	#add effort
	lanices$HI$UnitEffort<-"kWd"
	effort<-CEc@ce%>%group_by(time,space,technical)%>%summarise(value=sum(effKwDays,na.rm=T))%>%
		ungroup%>%mutate(Year=substr(as.character(time),1,4),Season=substr(as.character(time),8,8))%>%
		transmute(Year,Season,Fleet=technical,FishingArea=space,value)%>%
		ungroup()
	if(param$timestratif=="Year"){effort<-effort%>%mutate(Season=Year)}
	effort<-merge(lanices$HI,effort,all.x=T)
	effort$Effort<-effort$value
	effort<-effort[,names(lanices$HI)]

	lanices$HI<-effort
	lanices$HI<-lanices$HI%>%mutate(Effort=ifelse(is.na(Effort),0,Effort))
	lanices$SI<-lanices$SI%>%mutate(CATON=ifelse(CATON==-9,0,CATON))
	######################################################
	#correction division
	lanices$SI<-corrarea(lanices$SI)
	lanices$HI<-corrarea(lanices$HI)
	lanices$SD<-corrarea(lanices$SD)
	######################################################
	#SD,SI,HI various correction
	lanices$SD$NumberCaught<-as.numeric(lanices$SD$NumberCaught)
	lanices$SD$varNumLanded<-as.numeric(lanices$SD$varNumLanded)
	lanices$SI$varCATON<--9
	lanices$SI$OffLandings<--9
	if(nrow(lanices$SD)>0){ lanices$SD$CANUMtype<-"Lngt" }
	if(nrow(lanices$SD)>0){ lanices$SD$varNumLanded<--9 }
	if(nrow(lanices$SD)>0){ lanices$SD<-lanices$SD%>%filter(NumberCaught>0)}
	lanices$SI$CATON[lanices$SI$CATON==Inf]<-0

	#age
	#ices file generation
	ageices<-makeICdf(dbelan,Country="FR",Species=info$taxon,ReportingCategory="R",
				  SamplesOrigin="M",CANUMtype="age",Usage="H",unitMeanWeight="kg")
	#add nb of sampling and length for filtering
	ageices<-addsampmeas(ageices,dbelan)
	#add discards if any
	disices<-makeICdf(dbedis,Country="FR",Species=info$taxon,ReportingCategory="R",SamplesOrigin="M",
			  	CANUMtype="age",Usage="H",unitMeanWeight="kg")
	#add nb of sampling and length for filtering
	disices<-addsampmeas(disices,dbedis)
	######################################################
	#bind SI and SI	AND filter only on the considered year
	ageices$SD<-rbind(ageices$SD,disices$SD)%>%#filter(Year%in%year)
		filter(!is.na(Fleet))

	######################################################
	#correction division
	ageices$SD<-corrarea(ageices$SD)
	######################################################
	#SD,SI,HI various correction
	ageices$SD$NumberCaught<-as.numeric(ageices$SD$NumberCaught)
	ageices$SD$varNumLanded<-as.numeric(ageices$SD$varNumLanded)
	if(nrow(ageices$SD)>0){ ageices$SD$varNumLanded<--9 }
	#disices$SI<-disices$SI[disices$SI$CATON>=0,]
	lanices$SDage<-ageices$SD

	return(lanices)


}

trucapres<-function(){
#exporticeslength<-function(year,CLc,CLv,CEc,dbelan,dbedis,wdbelan,wdbedis,param,info,rtp,fitvb,checkn=100,checks=1,checksdis=3,wrtp=T,excel=F,age=T,dbelanage,wdbelange,dbedisage,wdbedisage){

	######################################################
	#weight
	#add meanWeight via wdbelan+completion using rtp
	wfishl<-wdbelan@lenStruc$estim%>%transmute(
	Year=sapply(strsplit(as.character(time)," - "),FUN = function(x){x[1]}),
	Season=sapply(strsplit(as.character(time)," - "),FUN = function(x){x[2]}),
						Fleet=as.character(technical),
						AgeLength=as.numeric(as.character(length)),
						FishingArea=space,CatchCategory="L",
					wfish=value)%>%group_by(Year,Season,Fleet,FishingArea,CatchCategory,AgeLength)%>%
		summarise(wfish=mean(wfish,na.rm=T))%>%ungroup()
	wfishd<-wdbedis@lenStruc$estim%>%transmute(
	Year=sapply(strsplit(as.character(time)," - "),FUN = function(x){x[1]}),
	Season=sapply(strsplit(as.character(time)," - "),FUN = function(x){x[2]}),
						Fleet=as.character(technical),
						AgeLength=as.numeric(as.character(length)),
						FishingArea=space,CatchCategory="D",
					wfish=value)%>%group_by(Year,Season,Fleet,FishingArea,CatchCategory,AgeLength)%>%
		summarise(wfish=mean(wfish,na.rm=T))%>%ungroup()
	wfish<-rbind(wfishl,wfishd)%>%distinct()
	#use rtp to complete empty weights
	wfish$wfish[is.na(wfish$wfish)]<-1000*rtp$a*(wfish$AgeLength[is.na(wfish$wfish)]/10)^rtp$b
	wfish<-wfish%>%mutate(AgeLength=as.character(AgeLength))
	wfish<-distinct(wfish)
	#pb de jointure ici
	pipo<-left_join(lanices$SD,wfish)%>%mutate(MeanWeight=wfish/1000) %>%select(-wfish)
	if(nrow(pipo)==nrow(lanices$SD)){ lanices$SD<-pipo}else{stop("merging pb")}
	lanices$SD<-left_join(lanices$SD,wfish)%>%mutate(MeanWeight=wfish/1000)%>%select(-wfish)
	#complete with rtp if needed
	lanices$SD$MeanWeight[is.na(lanices$SD$MeanWeight)]<-(1)*rtp$a*(as.numeric(lanices$SD$AgeLength[is.na(lanices$SD$MeanWeight)])/10)^rtp$b
	if(any(is.na(lanices$SD$MeanWeight))){ stop("w pb")}
	######################################################
	#verif et correction SOP
	lanices$SD$w<-as.numeric(lanices$SD$NumberCaught)*lanices$SD$MeanWeight
	verif<-lanices$SD%>%group_by(Season,Year,Fleet,FishingArea,CatchCategory)%>%summarise(totw=sum(w,na.rm=T))%>%ungroup()
	verif<-left_join(verif,lanices$SI) %>%select(Season,Year,Fleet,FishingArea,CatchCategory,CATON,totw)%>%ungroup()%>%
		mutate(sop=totw/CATON)%>%mutate(sop=ifelse(is.na(sop),1,sop),facteur=ifelse(.8<=sop&sop<=1.2,sop,0.8))
	verif$facteur[verif$sop>1.2]<-1.2
	if(any(is.na(lanices$SD$MeanWeight))){ stop("w pb 1")}

	#correction
	lanices$SD<-left_join(lanices$SD,verif)%>%mutate(MeanWeight=MeanWeight/sop*facteur)
	lanices$SD$w<-as.numeric(lanices$SD$NumberCaught)*lanices$SD$MeanWeight
	verif<-lanices$SD%>%group_by(Season,Year,Fleet,FishingArea,CatchCategory)%>%summarise(totw=sum(w,na.rm=T))%>%ungroup()
	verif<-left_join(verif,lanices$SI)%>%select(Season,Year,Fleet,FishingArea,CatchCategory,CATON,totw)%>%
		mutate(sop=totw/CATON)
	if(any(is.na(lanices$SD$MeanWeight))){ stop("w pb 2")}

	#adding stratum in intercatch file
	if(F){
		metconv<-dbelan@strataDesc@tcRec
		metconv<-data.frame(Fleet=metconv$to,from=metconv$from)%>%
			group_by(Fleet)%>%summarise(from=paste(gsub(",","/",unique(from)),collapse="/"))%>%
			ungroup()%>%mutate(from=ifelse(grepl("MIS_MIS",from),"other",gsub("_","",gsub("_0","",from))))
		lanices$SI<-left_join(lanices$SI,metconv)%>%mutate(InfoStockCoordinator=substr(from,1,249))%>%select(-from)
	}
	

	#file generation 
	pipo<-makeICdf(dbelan,Country="FR",Species=info$taxon,ReportingCategory="R",
		       SamplesOrigin="M",CANUMtype="length",Usage="H",unitMeanWeight="kg")
	HI<-lanices$HI[,names(pipo$HI)]
	SI<-lanices$SI[,names(pipo$SI)]
	SD<-lanices$SD[,names(pipo$SD)]
	makeICfile(list(HI=HI,SI=SI,SD=SD), filename=paste0("./ices/FRA",year,gsub("-","",info$stock),"_length.csv"), append=FALSE)
	#on sauve tout
	save(HI,SI,SD,file="HISISD.Rdata")
	if(age){
		SDage<-exportage(year,dbelanage,dbedisage,SI,wdbelanage,wdbedisage,param,rtp,fitvb,checks,checksdis,checkn)
		verifage<-SDage[[2]]
		SDage<-SDage[[1]]
		SDage<-SDage%>%filter(Year%in%year)
		makeICfile(list(HI=HI,SI=SI,SD=SDage), filename=paste0("./ices/FRA",year,gsub("-","",info$stock),"_age.csv"), append=FALSE)
		#on sauve tout
		save(HI,SI,SD,SDage,file="HISISD.Rdata")
	}

if(excel){
	require(openxlsx)
	require(tidyr)
	source("./data2icesxls.R")
	#compute area official stuff
	#CLc and CLv have the same order in line
	clrect<-CLv@cl%>%dplyr::select(foCatEu6,rect,area,quarter,landWt)%>%
		mutate(final=CLc@cl$technical,quarter=paste0("q",quarter))%>%filter(!is.na(final))%>%
		group_by(rect,area,quarter)%>%summarise(tot=sum(landWt,na.rm=T)/1000)%>%
			mutate(quarter=factor(quarter,c("q1","q2","q3","q4")))%>%
			ungroup()

			datclrect<-tidyr::spread(clrect,quarter,tot,drop=F)%>%transmute(rect,q1,q2,q3,q4,area)
			dim(datclrect)
			icesrect<-COSTeda::ICESAreaRects%>%transmute(lat,lon,rect=gsub(" ","",as.character(StatRect)),areabis=division,subdiv=subdivision)%>%distinct()
			datclrect<-left_join(datclrect,icesrect)%>%transmute(rect,q1,q2,q3,q4,aa="",lat,lon,area)
			SIrect<-data.frame(datclrect)
			#kg to t
			#save(HI,SI,SIrect,info,file="testexport.Rdata")
			SItmp<-SI%>%mutate(CATON=CATON/1000)
			accfile<-paste0("./ices/FRA_",year,"_",info$wg,"_",info$stock,".xlsx")
			if(nrow(SD)>0){
				if(age){
					data2icesxls(info=info,param=param,HI=HI,SI=SItmp,
						     SIrect=SIrect,SDlength=SD,SDage=SDage,nomfich=accfile)
				}else{
					data2icesxls(info=info,param=param,HI=HI,SI=SItmp,
						     SIrect=SIrect,SDlength=SD,SDage="",nomfich=accfile)
				}
			}else{
				if(age){
					data2icesxls(info=info,param=param,HI=HI,SI=SItmp,
						     SIrect=SIrect,SDlength="",SDage=SDage,nomfich=accfile)
				}else{
					data2icesxls(info=info,param=param,HI=HI,SI=SItmp,
						     SIrect=SIrect,SDlength="",SDage="",nomfich=accfile)
				}
			}
			#system(paste0("mv ",accfile,"./ices/",accfile))
			#system(paste0("mv FRA_2016_",info$stock,"simple.xlsx ./ices/FRA_2016_",info$stock,"simple.xlsx"))
}

#some graph as output
#pander::pander(tidyr::unite(verif,"Season/Fleet/Area",Season,Fleet,FishingArea,CatchCategory,sep="/"),style="simple",split.tables=Inf)

plt1len<-ggplot(SD,aes(x=as.numeric(as.character(AgeLength)),y=NumberCaught,fill=Season,color=Season))+
	ylab("Numbers of individuals")+xlab("Length")+
	facet_grid(Fleet+CatchCategory~Year+FishingArea,scale="free") +
	geom_line(alpha=.6)+#+geom_point(alpha=.6)
	ggtitle("Final length distribution")+
	theme(axis.text.x = element_text(size=8, angle=90),
	axis.text.y = element_text(size=10, angle=0),
	strip.text.x=element_text(size=8,angle=90),
strip.text.y=element_text(size=8,angle=0),
legend.position="bottom")

veriflen<-verif
#graph w@length
	plt1wlen<-ggplot(SD,aes(x=as.numeric(as.character(AgeLength)),y=MeanWeight,fill=Season,color=Season))+
		ylab("Individual weight (kg)")+xlab("Length (mm)")+
		facet_grid(Fleet+CatchCategory~Year+FishingArea)+#,scale="free") +
		geom_line(alpha=.6)+#+geom_point(alpha=.6)
		geom_point(alpha=.6)+#+geom_point(alpha=.6)
		ggtitle("Individual weigth at length")+
		theme(axis.text.x = element_text(size=8, angle=90),
		axis.text.y = element_text(size=10, angle=0),
		strip.text.x=element_text(size=8,angle=90),
	strip.text.y=element_text(size=8,angle=0),
	legend.position="bottom")


if(age){
	plt1age<-ggplot(SDage,aes(x=as.numeric(as.character(AgeLength)),y=NumberCaught,fill=Season,color=Season))+
		ylab("Numbers of individuals")+xlab("Age")+
		facet_grid(Fleet+CatchCategory~Year+FishingArea,scale="free") +
		geom_line(alpha=.6)+#+geom_point(alpha=.6)
		geom_point(alpha=.6)+#+geom_point(alpha=.6)
		ggtitle("Final age distribution")+
		theme(axis.text.x = element_text(size=8, angle=90),
		axis.text.y = element_text(size=10, angle=0),
		strip.text.x=element_text(size=8,angle=90),
	strip.text.y=element_text(size=8,angle=0),
	legend.position="bottom")
	plt1wage<-ggplot(SDage,aes(x=as.numeric(as.character(AgeLength)),y=MeanWeight,fill=Season,color=Season))+
		ylab("Individual weight (kg)")+xlab("Age (year)")+
		facet_grid(Fleet+CatchCategory~Year+FishingArea)+#,scale="free") +
		geom_line(alpha=.6)+#+geom_point(alpha=.6)
		geom_point(alpha=.6)+#+geom_point(alpha=.6)
		ggtitle("Individual weigth at age")+
		theme(axis.text.x = element_text(size=8, angle=90),
		axis.text.y = element_text(size=10, angle=0),
		strip.text.x=element_text(size=8,angle=90),
	strip.text.y=element_text(size=8,angle=0),
	legend.position="bottom")

	return(list(HI=HI,SI=SI,SD=SD,SDage=SDage,plt1len=plt1len,plt1age=plt1age,veriflen=veriflen,verifage=verifage,
		    plt1wlen=plt1wlen,plt1wage=plt1wage))
}else{
	return(list(HI=HI,SI=SI,SD=SD,plt1len=plt1len,veriflen=veriflen,plt1wlen=plt1wlen))
}


}
#function for exporticeslen
addsampmeas<-function(lanices,dbelan){
	#add nb of sampling and length for filtering
	nsamp<-dbelan@nSamp$len%>%transmute(time,FishingArea=space,Fleet=technical,s=value)
	nsamp$Season = sapply (strsplit (as.character (nsamp$time), " - "), FUN = function(x){x[2]}) 
	nsamp$Year = sapply (strsplit (as.character (nsamp$time), " - "), FUN = function(x){x[1]}) 
	nsamp<-nsamp%>%select(-time)
	nlen<-dbelan@nMeas$len%>%transmute(time,FishingArea=space,Fleet=technical,n=value)
	nlen$Season = sapply (strsplit (as.character (nlen$time), " - "), FUN = function(x){x[2]}) 
	nlen$Year = sapply (strsplit (as.character (nlen$time), " - "), FUN = function(x){x[1]}) 
	nlen<-nlen%>%select(-time)
	lanices$SD<-left_join(lanices$SD,nsamp)%>%mutate(s=ifelse(is.na(s),0,s),NumSamplesLngt=s)%>%select(-s)
	lanices$SD<-left_join(lanices$SD,nlen)%>%mutate(n=ifelse(is.na(n),0,n),NumLngtMeas=n)%>%select(-n)
	lanices$SD<-lanices$SD%>%filter(!grepl("MIS_MIS",Fleet))
	lanices$SI<-left_join(lanices$SI,nsamp)%>%mutate(s=ifelse(is.na(s),0,s))
	lanices$SI<-left_join(lanices$SI,nlen)%>%mutate(n=ifelse(is.na(n),0,n))
	lanices$SI<-lanices$SI%>%mutate(InfoGeneral=paste0("nsamp:",s,"/nlength:",n))
	return(lanices)
}
#add Div or SubArea for AreaTYpe
checkletter<-function(a){rez<-FALSE;for(i in 1:26){rez<-rez|grepl(letters[i],a)};return(rez)}
corrarea<-function(SI){
	testlett<-checkletter(SI$FishingArea)
	SI$AreaType<-ifelse(testlett,"Div","SubArea")
	testFU<-grepl("FU",SI$FishingArea)
	SI$AreaType[testFU]<-"Div"
	testFU<-grepl("FU.22",SI$FishingArea)
	SI$AreaType[testFU]<-"SubDiv"
	return(SI)
}




exportage<-function(year,dbelan,dbedis,SI,wdbelan,wdbedis,param,rtp,fitvb,checks,checksdis,checkn){
	#ices file generation
	ageices<-makeICdf(dbelan,Country="FR",Species=info$taxon,ReportingCategory="R",
				  SamplesOrigin="M",CANUMtype="age",Usage="H",unitMeanWeight="kg")
	#add nb of sampling and length for filtering
	ageices<-addsampmeas(ageices,dbelan)
	#filter
	ageices$SD<-ageices$SD%>%filter(NumSamplesLngt>=checks & NumLngtMeas>=checkn)
	#on vire MIS_MIS aussi au cas où
	ageices$SD<-ageices$SD%>%filter(!grepl("MIS_MIS",Fleet))

	#add discards if any
	disices<-makeICdf(dbedis,Country="FR",Species=info$taxon,ReportingCategory="R",SamplesOrigin="M",
			  	CANUMtype="age",Usage="H",unitMeanWeight="kg")
	#add nb of sampling and length for filtering
	disices<-addsampmeas(disices,dbedis)
	#remove strata according to checks and checks in SI and SD
	disices$SD<-disices$SD%>%filter(NumSamplesLngt>=checksdis & NumLngtMeas>=checkn)
	######################################################
	#bind SI and SI	AND filter only on the considered year
	ageices$SD<-rbind(ageices$SD,disices$SD)%>%filter(Year%in%year)
	######################################################
	#correction division
	ageices$SD<-corrarea(ageices$SD)

	######################################################
	#SD,SI,HI various correction
	ageices$SD$NumberCaught<-as.numeric(ageices$SD$NumberCaught)
	ageices$SD$varNumLanded<-as.numeric(ageices$SD$varNumLanded)
	if(nrow(ageices$SD)>0){ ageices$SD$varNumLanded<--9 }
	#disices$SI<-disices$SI[disices$SI$CATON>=0,]
	#disices$SD<-disices$SD[disices$SD$NumberCaught>=0,]

	#weight
	#add meanWeight via wdbelan+completion using rtp
	wfishl<-wdbelan@ageStruc$estim%>%transmute(
	Year=sapply(strsplit(as.character(time)," - "),FUN = function(x){x[1]}),
	Season=sapply(strsplit(as.character(time)," - "),FUN = function(x){x[2]}),
						Fleet=as.character(technical),
						AgeLength=as.numeric(as.character(age)),
						FishingArea=space,CatchCategory="L",
					wfish=value)%>%group_by(Year,Season,Fleet,FishingArea,CatchCategory,AgeLength)%>%
		summarise(wfish=mean(wfish,na.rm=T))%>%ungroup()
	wfishd<-wdbedis@ageStruc$estim%>%transmute(
	Year=sapply(strsplit(as.character(time)," - "),FUN = function(x){x[1]}),
	Season=sapply(strsplit(as.character(time)," - "),FUN = function(x){x[2]}),
						Fleet=as.character(technical),
						AgeLength=as.numeric(as.character(age)),
						FishingArea=space,CatchCategory="D",
					wfish=value)%>%group_by(Year,Season,Fleet,FishingArea,CatchCategory,AgeLength)%>%
		summarise(wfish=mean(wfish,na.rm=T))%>%ungroup()
	wfish<-rbind(wfishl,wfishd)%>%distinct()

	wfish$wfish[is.na(wfish$wfish)]<-
		1000*rtp$a*(predict(fitvb,data.frame(age=wfish$AgeLength[is.na(wfish$wfish)]))/10)^rtp$b

	ggplot(wfish,aes(x=AgeLength,y=wfish,color=Season))+geom_point()+facet_grid(Fleet~Year)
	wfish<-wfish%>%mutate(AgeLength=as.character(AgeLength))
	wfish<-distinct(wfish)

	#pb de jointure ici
	pipo<-left_join(ageices$SD,wfish)%>%mutate(MeanWeight=wfish/1000)%>%select(-wfish)
	if(nrow(pipo)==nrow(ageices$SD)){ ageices$SD<-pipo}else{stop("merging pb")}

	ageices$SD<-left_join(ageices$SD,wfish)%>%mutate(MeanWeight=wfish/1000)%>%select(-wfish)
	ageices$SD$MeanWeight[is.na(ageices$SD$MeanWeight)]<-(1/1000)*rtp$a*(as.numeric(ageices$SD$AgeLength)[is.na(ageices$SD$MeanWeight)]/10)^rtp$b
	#ggplot(ageices$SD,aes(x=as.numeric(AgeLength),y=MeanWeight,color=Season))+geom_point()+facet_wrap(Fleet~Year)

	#ageices$SD$MeanWeight<-rtp$a*(as.numeric(ageices$SD$AgeLength)/10)^rtp$b
	#verif et correction SOP

	ageices$SD$w<-as.numeric(ageices$SD$NumberCaught)*ageices$SD$MeanWeight
	verif<-ageices$SD%>%group_by(Season,Year,Fleet,FishingArea,CatchCategory)%>%summarise(totw=sum(w,na.rm=T))%>%ungroup()
	verif<-left_join(verif,SI) %>%select(Season,Year,Fleet,FishingArea,CatchCategory,CATON,totw)%>%
		mutate(sop=totw/CATON,facteur=ifelse(.8<=sop&sop<=1.2,sop,0.8))
	verif$facteur[verif$sop>1.2]<-1.2

	ageices$SD<-left_join(ageices$SD,verif)%>%mutate(MeanWeight=MeanWeight/sop*facteur)
	ageices$SD$w<-as.numeric(ageices$SD$NumberCaught)*ageices$SD$MeanWeight
	verif<-ageices$SD%>%group_by(Season,Year,Fleet,FishingArea,CatchCategory)%>%summarise(totw=sum(w,na.rm=T))%>%ungroup()
	verif<-left_join(verif,SI)%>%select(Season,Year,Fleet,FishingArea,CatchCategory,CATON,totw)%>%
		mutate(sop=totw/CATON)

	#pander::pander(tidyr::unite(verif,"Season/Fleet/Area",Season,Fleet,FishingArea,CatchCategory,sep="/"),style="simple",split.tables=Inf)

	plt1<-ggplot(ageices$SD,aes(x=as.numeric(as.character(AgeLength)),y=NumberCaught,fill=Season,color=Season))+
		ylab("Numbers of individuals")+xlab("length")+
		facet_grid(Fleet+CatchCategory~Year,scale="free") +
		geom_line(alpha=.6)+#+geom_point(alpha=.6)
		ggtitle("Final landings")+
		theme(axis.text.x = element_text(size=8, angle=90),
		axis.text.y = element_text(size=10, angle=0),
		strip.text.x=element_text(size=8,angle=90),
	strip.text.y=element_text(size=8,angle=0),
	legend.position="bottom")

	#disices$SI<-disices$SI[disices$SI$CATON>=0,]
	#disices$SD<-disices$SD[disices$SD$NumberCaught>=0,]
	#SI<-SI[!is.na(SI$Year),]
	#SI$CATON[SI$CATON==Inf]<-0
	pipo<-makeICdf(dbelan,Country="FR",Species=info$taxon,ReportingCategory="R",
		       SamplesOrigin="M",CANUMtype="age",Usage="H",unitMeanWeight="kg")
	#add meanLength using fitvb
						     
	ageices$SD$MeanLength[is.na(ageices$SD$MeanLength)]<-round(predict(fitvb,data.frame(age=as.numeric(ageices$SD$AgeLength[is.na(ageices$SD$MeanLength)])))/10,2)
	if(nrow(ageices$SD)>0){
	ageices$SD$UnitMeanLength<-"cm"
	}
	#ageices$SD$MeanLength<-sprintf("%09.0f",round(predict(fitvb,data.frame(age=as.numeric(ageices$SD$AgeLength))),0))
							     

	SDage<-ageices$SD[,names(pipo$SD)]
	return(list(SDage,verif,plt1))
}


