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
#Dorymymrex grandulus present twice
ga<-ga[!duplicated(ga$Taxon),]


require(plyr)
comb<-join(fl,ga,by="Taxon",type="full",match="all")
comb$FL[is.na(comb$FL)]<-0
comb$GA[is.na(comb$GA)]<-0

#add in Doug's notes about which taxa are important
dbn<-read.csv(file="~/Desktop/ants_for_doug/DBB_on_FL_&_GA_Combined_Species_w_GB_Info.csv",stringsAsFactors=FALSE)
#Dorymymrex grandulus present twice
dbn<-dbn[!duplicated(dbn$Taxon),]

comb$DBB_Notes<-NA
for(x in 1:nrow(comb)){
	comb$DBB_Notes[x]<-dbn$Notes[dbn$Taxon %in% comb$Taxon[x]]
	#print(paste(comb$Taxon[x], "______",dbn$Notes[dbn$Taxon %in% comb$Taxon[x]]))
}
comb$MPN_Notes<-NA

#And add in alternates, which were selected for those deemed of importance (except Leptogenys, which was not important, but included anyways)
alts<-read.csv(file="~/Desktop/ants_for_doug/alternates_for_use.csv",stringsAsFactors=FALSE)
alts$Genus<-NA
alts$Species<-NA
alts$MPN_Notes<-NA
alts$FL<-0
alts$GA<-0
for(x in 1:nrow(alts)){
	alts$Genus[x]<-strsplit(alts$Alternate[x]," ")[[1]][1]
	alts$Species[x]<-strsplit(alts$Alternate[x]," ")[[1]][2]
	alts$MPN_Notes[x]<-paste("Alternate for ",alts$Taxon[x],sep="")
}

alts.to.add<-alts[!is.na(alts$Alternate),]
alts.to.add<-alts.to.add[,c("Alternate","Genus","Species","FL","GA","DBB_Notes","MPN_Notes")]
colnames(alts.to.add)<-c("Taxon","Genus","Species","FL","GA","DBB_Notes","MPN_Notes")

combnew<-rbind(comb,alts.to.add)
write.csv(combnew,file="FL_GA_list_w_alts.csv",row.names=FALSE)


#read in list of FL and GA taxa with alternates
comb<-read.csv(file="FL_GA_list_w_alts.csv",stringsAsFactors=FALSE)
#add taxonomy from Bolton
tax<-read.csv(file="~/Desktop/ants_for_doug/Bolton_List_of_Valid_Species_27_July_2018.csv",header=TRUE,stringsAsFactors=FALSE)

invalid<-comb$Taxon[!comb$Taxon %in% tax$TaxonName]
invalid

colnames(tax)[1]<-"Taxon"
comb.tax<-join(comb,tax,by="Taxon",type="left",match="all")
comb.tax<-comb.tax[,c(8,9,2,3,1,4:5,11,13:24,6:7)]
comb.tax.ord<-comb.tax[with(comb.tax,order(Subfamily,Tribe,Genus,Species)),]

#read from NCBI pull
gb<-read.csv(file="~/Desktop/ants_for_doug/Formicidae.combined.acc.nos.csv",stringsAsFactors=FALSE)
gb<-gb[,c(1:14)]
comb.tax.ord.gb<-join(comb.tax.ord,gb,by="Taxon",type="left",match="all")
comb.tax.ord.gb$GB<-0
comb.tax.ord.gb$GB[!is.na(comb.tax.ord.gb$NCBI.ID)]<-1
write.csv(comb.tax.ord.gb,file="~/Desktop/ants_for_doug/FL_&_GA_Combined_Species_w_GB_Info_w_Alts.csv",row.names=FALSE)
present<-comb.tax.ord.gb$Taxon[comb.tax.ord.gb$GB==1]
missing<-comb.tax.ord.gb$Taxon[comb.tax.ord.gb$GB==0]










#17Nov2019
#Subsequently added in more accessions that are shown in colors in FL_&_GA_Combined_Species_w_GB_Info_w_Alts_17Nov19.xlsx

#save individual accession number files
gb<-read.csv(file="/Users/matthewnelsen/Documents/papers_reviews/papers/ants_for_doug/FL_&_GA_Combined_Species_w_GB_Info_w_Alts_17Nov19.csv",stringsAsFactors=FALSE)
loci<-c("nuSSU", "nuLSU", "AbdA", "COI", "LR", "Wg", "EF1aF1", "EF1aF2", "ArgK", "CAD", "Top1", "Ubx")
for(x in 1:length(loci)){
	gg<-as.data.frame(as.matrix(gb[,loci[x]][!is.na(gb[,loci[x]])]))
	write.table(gg,file=paste("/Users/matthewnelsen/Documents/papers_reviews/papers/ants_for_doug/17nov2019_1/Formicidae/",loci[x],".acc.nos.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE)
}

#use 3_multi_seq_fetcher_17nov19.py to retrieve sequences
#use 4_multi_parser_revised_17nov19.py to parse sequence files
#use 5_name_changer_17nov19.r
#use 6_blasting_ants_tax_hybrid_family2_17nov19.py
#use 7_make_out_files_17nov19.r
#use 8_add.outs.to.master_17nov19.r
#use 9_concatenate_outs_and_reals_17nov19.txt
#use 10_initial_align_fullalign_einsi_nogappy12_defaultmacse.r
#but note that in 11, i chose to not exclude the really short ones and made another loop to remove some that aligned strangely
#use 11_profile_aligning_fullalign_einse_nogappy6_short_altorg_sparseseqs2_defaultmacse_correctedlength.r 
#concatenate in mesquite
