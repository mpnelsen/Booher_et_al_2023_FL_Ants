#from 4_name.changer.regular.seqs.r

require(seqinr)
loci<-c("nuSSU","nuLSU","AbdA","COI","LR","Wg","EF1aF1","EF1aF2","ArgK","CAD","Top1","Ubx")
dat.type<-c("DNA","DNA","DNA","DNA","DNA","DNA","DNA","DNA","DNA","DNA","DNA","DNA")
seq.refs<-read.csv(file="/Users/matthewnelsen/Documents/papers_reviews/papers/ants_for_doug/FL_&_GA_Combined_Species_w_GB_Info_w_Alts_17Nov19.csv",stringsAsFactors=FALSE)
#add these in
#Pheidole_dentigula DB1AC_ConsensusCO1 DB1AC_CONSENSO.fas -> Pheidole_dentigula_CO1.fasta
#Myrmica_pinetorum DB4AC_Consensus DB4AC_Consensus.fas -> Myrmica_pinetorum_CO1.fasta
#Crematogaster_ashmeadi DB2AC_LCO_CO1 DB2AC_LCO_MPNtrim.fas -> Crematogaster_ashmeadi_CO1.fasta
seq.refs$GB[seq.refs$Taxon %in% c("Pheidole dentigula","Myrmica pinetorum","Crematogaster ashmeadi")]<-1
seq.refs$COI[seq.refs$Taxon %in% c("Pheidole dentigula","Myrmica pinetorum","Crematogaster ashmeadi")]<-"Manual"
write.csv(seq.refs,file="/Users/matthewnelsen/Documents/papers_reviews/papers/ants_for_doug/FL_&_GA_Combined_Species_w_GB_Info_w_Alts_17Nov19_extras.csv",row.names=FALSE)

#read in revised
seq.refs<-read.csv(file="/Users/matthewnelsen/Documents/papers_reviews/papers/ants_for_doug/FL_&_GA_Combined_Species_w_GB_Info_w_Alts_17Nov19_extras.csv",stringsAsFactors=FALSE)

path<-"/Users/matthewnelsen/Documents/papers_reviews/papers/ants_for_doug/17nov2019_1/Formicidae"

for(locus in loci[3:12]){
	align<-read.fasta(file=paste(path,"/",locus,"/",locus,"_CDS_Parsed.fasta",sep=""),seqtype="DNA",forceDNAtolower=FALSE)
	if(locus=="COI"){
		#add in extras
		g.align<-read.fasta(file="/Users/matthewnelsen/Documents/papers_reviews/papers/ants_for_doug/fwdco1sequences/Pheidole_dentigula_CO1.fasta",seqtype="DNA",forceDNAtolower=FALSE)
		combo.align<-c(g.align,align)
		g.align<-read.fasta(file="/Users/matthewnelsen/Documents/papers_reviews/papers/ants_for_doug/fwdco1sequences/Myrmica_pinetorum_CO1.fasta",seqtype="DNA",forceDNAtolower=FALSE)
		combo.alignb<-c(g.align,combo.align)
		g.align<-read.fasta(file="/Users/matthewnelsen/Documents/papers_reviews/papers/ants_for_doug/fwdco1sequences/Crematogaster_ashmeadi_CO1.fasta",seqtype="DNA",forceDNAtolower=FALSE)
		combo.alignc<-c(g.align,combo.alignb)
		align<-combo.alignc
	}
	write.fasta(sequences=align,names=names(align),file.out=paste(path,"/",locus,"/",locus,"_Parsed.fasta",sep=""),nbchar=1000000)
}


for(locus in loci){
	require(ape)
	new.align<-read.fasta(file=paste(path,"/",locus,"/",locus,"_Parsed.fasta",sep=""),seqtype=dat.type[loci==locus],forceDNAtolower=FALSE)
	names(new.align)<-gsub("_"," ",names(new.align))
	for(p in 1:length(names(new.align))){
		#print(seq.refs$Taxon[seq.refs[,locus] %in% names(new.align)[p]])
		#print(length(new.align))
		#print(names(new.align)[duplicated(names(new.align))])
		#print(names(new.align)[p])
		if(!names(new.align)[p] %in% c("Pheidole dentigula","Myrmica pinetorum","Crematogaster ashmeadi")){
			names(new.align)<-gsub("[.].*","",names(new.align))
			names(new.align)[p]<-seq.refs$Taxon[seq.refs[,locus] %in% names(new.align)[p]]
		}
	}
	write.fasta(sequences=new.align,names=gsub(" ","_",names(new.align)),file=paste(path,"/",locus,"/",locus,"_Renamed_Parsed.fasta",sep=""),nbchar=1000000)
	names.w.accs<-seq.refs$Taxon[!is.na(seq.refs[,locus])]
	names.missing<-names.w.accs[!names.w.accs %in% names(new.align)]
	print(locus)
	print(names.missing)
	#print(seq.refs[,locus][names.missing %in% seq.refs$Taxon])
	seq.refs[seq.refs$Taxon %in% names.missing,locus]<-NA
	#line below will remove taxa lacking any ribosomal sequences
	#seq.refs<-seq.refs[rowSums(is.na(seq.refs[,loci]))!=length(loci),]
}


seq.refs$Taxon[is.na(seq.refs$nuSSU) & is.na(seq.refs$nuLSU) & is.na(seq.refs$AbdA) & is.na(seq.refs$Wg) & is.na(seq.refs$LR) & is.na(seq.refs$COI) & is.na(seq.refs$EF1aF1) & is.na(seq.refs$EF1aF2) & is.na(seq.refs$ArgK) & is.na(seq.refs$CAD) & is.na(seq.refs$Top1) & is.na(seq.refs$Ubx)]
seq.refs$Family<-"Formicidae"
write.csv(seq.refs,file=paste(path,"/","combined.updated.pruned.csv",sep=""),row.names=FALSE)

