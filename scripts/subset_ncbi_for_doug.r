
#read from NCBI pull
gb<-read.csv(file="~/Desktop/ants_for_doug/Formicidae.combined.acc.nos.csv",stringsAsFactors=FALSE)

#Read FL list
fl<-read.csv(file="~/Desktop/ants_for_doug/FLantsMattNelson-2_mpn_mods.csv",stringsAsFactors=FALSE)

#modify Pheidole morrisi to be Pheidole morrisii as it is in NCBI
fl$species[fl$species %in% "Pheidole morrisi"]<-"Pheidole morrisii"
fl$CurrentName[fl$CurrentName %in% "Pheidole morrisi"]<-"Pheidole morrisii"

#get rid of sp.'s
fl<-fl[1:214,]

#Check for which FL species absent from NCBI pull
missing<-fl$CurrentName[!fl$CurrentName %in% gb$Taxon]
missing

tax<-read.csv(file="~/Desktop/ants_for_doug/Bolton_List_of_Valid_Species_27_July_2018.csv",header=TRUE,stringsAsFactors=FALSE)
nrow(tax)
nrow(fl)

invalid<-fl$CurrentName[!fl$CurrentName %in% tax$TaxonName]
invalid


fl.gb<-gb[gb$Taxon %in% fl$CurrentName,]
write.csv(fl.gb,file="~/Desktop/ants_for_doug/NCBI_FL_ANTS.csv",row.names=FALSE)



loci<-c("nuSSU", "nuLSU", "AbdA", "COI", "LR", "Wg", "EF1aF1", "EF1aF2", "ArgK", "CAD", "Top1", "Ubx")

for(x in 1:length(loci)){
	gg<-as.data.frame(as.matrix(fl.gb[,loci[x]][!is.na(fl.gb[,loci[x]])]))
	write.table(gg,file=paste("~/Desktop/ants_for_doug/fl_ants_first_try/",loci[x],".acc.nos.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE)
}


ga<-read.csv(file="~/Desktop/ants_for_doug/GA_ants_mods.csv",stringsAsFactors=FALSE)
ga$SPECIES[ga$SPECIES %in% "morrisi"]<-"morrisii"
ga$CurrentName<-paste(ga$GENUS,ga$SPECIES,sep=" ")
ga<-ga[1:211,]
nrow(ga)
tax<-read.csv(file="~/Desktop/ants_for_doug/Bolton_List_of_Valid_Species_27_July_2018.csv",header=TRUE,stringsAsFactors=FALSE)
nrow(tax)

invalid<-ga$CurrentName[!ga$CurrentName %in% tax$TaxonName]
invalid

#read from NCBI pull
gb<-read.csv(file="~/Desktop/ants_for_doug/Formicidae.combined.acc.nos.csv",stringsAsFactors=FALSE)

#Check for which GA species absent from NCBI pull
present.ga<-sort(ga$CurrentName[ga$CurrentName %in% gb$Taxon])
present.ga

#Check for which GA species absent from NCBI pull
missing.ga<-sort(ga$CurrentName[!ga$CurrentName %in% gb$Taxon])
missing.ga






#make data frame with genus, species and higher-level bolton taxonomy and presence/absence from FL/GA

#Read FL list
fl<-read.csv(file="~/Desktop/ants_for_doug/FLantsMattNelson-2_mpn_mods.csv",stringsAsFactors=FALSE)

#modify Pheidole morrisi to be Pheidole morrisii as it is in NCBI
fl$species[fl$species %in% "Pheidole morrisi"]<-"Pheidole morrisii"
fl$CurrentName[fl$CurrentName %in% "Pheidole morrisi"]<-"Pheidole morrisii"

#get rid of sp.'s
fl<-fl[1:214,]
fl$Taxon<-fl$CurrentName
fl$Genus<-NA
fl$Species<-NA
fl$FL<-1
for(x in 1:nrow(fl)){
	fl$Genus[x]<-strsplit(fl$Taxon[x]," ")[[1]][1]
	fl$Species[x]<-strsplit(fl$Taxon[x]," ")[[1]][2]
}
fl<-fl[,c("Taxon","Genus","Species","FL")]


ga<-read.csv(file="~/Desktop/ants_for_doug/GA_ants_mods.csv",stringsAsFactors=FALSE)
ga$SPECIES[ga$SPECIES %in% "morrisi"]<-"morrisii"
ga$CurrentName<-paste(ga$GENUS,ga$SPECIES,sep=" ")
ga<-ga[1:211,]
ga$Taxon<-ga$CurrentName
ga$Genus<-NA
ga$Species<-NA
ga$GA<-1
for(x in 1:nrow(ga)){
	ga$Genus[x]<-strsplit(ga$Taxon[x]," ")[[1]][1]
	ga$Species[x]<-strsplit(ga$Taxon[x]," ")[[1]][2]
}
ga<-ga[,c("Taxon","Genus","Species","GA")]

require(plyr)
comb<-join(fl,ga,by="Taxon",type="full",match="all")
comb$FL[is.na(comb$FL)]<-0
comb$GA[is.na(comb$GA)]<-0

#add taxonomy from Bolton
tax<-read.csv(file="~/Desktop/ants_for_doug/Bolton_List_of_Valid_Species_27_July_2018.csv",header=TRUE,stringsAsFactors=FALSE)

invalid<-comb$Taxon[!comb$Taxon %in% tax$TaxonName]
invalid

colnames(tax)[1]<-"Taxon"
comb.tax<-join(comb,tax,by="Taxon",type="left",match="all")
comb.tax<-comb.tax[,c(6,7,2,3,1,4:5,9,11:22)]
comb.tax.ord<-comb.tax[with(comb.tax,order(Subfamily,Tribe,Genus,Species)),]


#read from NCBI pull
gb<-read.csv(file="~/Desktop/ants_for_doug/Formicidae.combined.acc.nos.csv",stringsAsFactors=FALSE)
gb<-gb[,c(1:14)]
comb.tax.ord.gb<-join(comb.tax.ord,gb,by="Taxon",type="left",match="all")
comb.tax.ord.gb$GB<-0
comb.tax.ord.gb$GB[!is.na(comb.tax.ord.gb$NCBI.ID)]<-1
write.csv(comb.tax.ord.gb,file="~/Desktop/ants_for_doug/FL_&_GA_Combined_Species_w_GB_Info.csv",row.names=FALSE)
present<-comb.tax.ord.gb$Taxon[comb.tax.ord.gb$GB==1]
missing<-comb.tax.ord.gb$Taxon[comb.tax.ord.gb$GB==0]


tots<-table(comb.tax.ord.gb$Genus)
tots.pres<-table(comb.tax.ord.gb$Genus[comb.tax.ord.gb$GB==1])
tots.miss<-table(comb.tax.ord.gb$Genus[comb.tax.ord.gb$GB==0])

gen.sum<-as.data.frame(tots)
colnames(gen.sum)<-c("Genus","NoSpecies")
gen.sum[,c("Present","Missing")]<-0

for(x in 1:nrow(gen.sum)){
	if(gen.sum$Genus[x] %in% names(tots.pres)){
		gen.sum$Present[x]<-tots.pres[names(tots.pres)==gen.sum$Genus[x]][[1]]
	}
	if(gen.sum$Genus[x] %in% names(tots.miss)){	
		gen.sum$Missing[x]<-tots.miss[names(tots.miss)==gen.sum$Genus[x]][[1]]
	}
}

gen.sum$PropMissing<-round(gen.sum$Missing/gen.sum$NoSpecies,2)
gen.sum

gen.sum[gen.sum$PropMissing==1,]
gen.sum[gen.sum$PropMissing>0.33 & gen.sum$PropMissing<1,]
write.csv(gen.sum,file="~/Desktop/ants_for_doug/FL_&_GA_Combined_Genera_w_GB_Info.csv",row.names=FALSE)
