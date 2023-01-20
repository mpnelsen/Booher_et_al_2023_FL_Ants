master<-read.csv(file="/Users/matthewnelsen/Documents/papers_reviews/papers/ants_for_doug/17nov2019_1/Formicidae/combined.updated.updated.pruned.csv",stringsAsFactors=FALSE)
#master<-master[,2:28]
new.adds<-read.csv(file="/Users/matthewnelsen/Documents/papers_reviews/papers/ants_for_doug/FL_&_GA_Combined_Species_w_GB_Info_w_Alts_17Nov19.csv",stringsAsFactors=FALSE)
new.adds<-new.adds[new.adds$Taxon %in% c("Pheidole_dentigula","Myrmica_pinetorum","Crematogaster_ashmeadi"),]
master<-merge(master,new.adds,all=TRUE)
outs<-read.csv(file="/Users/matthewnelsen/Documents/papers_reviews/papers/ant_plant_interactions/ants_june_to_august_2016/antcat_distribution_4_june_2016/outgroup_things/standards.w.taxonomy.dupaccsrem.alltax.wouts.outsonly.csv",stringsAsFactors=FALSE)
outs<-outs[,c(1:10,12:18)]
combo<-merge(master,outs,all=TRUE)
for(v in 1:nrow(combo)){
	if(combo$Tribe[v]=="Incertae_Sedis"){
		combo$Tribe[v]<-combo$Genus[v]
	}
}
combo$Order<-"Hymenoptera"
for(x in 1:nrow(combo)){
	if(is.na(combo$Tribe[x]) | combo$Tribe[x] %in% ""){
		combo$Tribe[x]<-combo$Genus[x]
	}
}
write.csv(combo,file="/Users/matthewnelsen/Documents/papers_reviews/papers/ants_for_doug/17nov2019_1/Formicidae/combined.updated.updated.pruned.wouts.csv",row.names=FALSE)
write.csv(combo,file="/Users/matthewnelsen/Documents/papers_reviews/papers/ants_for_doug/17nov2019_1/Formicidae/alignments/Formicidae/combined.updated.updated.pruned.wouts.csv",row.names=FALSE)

#Pheidole_dentigula DB1AC_ConsensusCO1 DB1AC_CONSENSO.fas -> Pheidole_dentigula_CO1.fasta
#Myrmica_pinetorum DB4AC_Consensus DB4AC_Consensus.fas -> Myrmica_pinetorum_CO1.fasta
#Crematogaster_ashmeadi DB2AC_LCO_CO1 DB2AC_LCO_MPNtrim.fas -> Crematogaster_ashmeadi_CO1.fasta