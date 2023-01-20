require(seqinr)
loci<-c("nuSSU","nuLSU","AbdA","COI","LR","Wg","EF1aF1","EF1aF2","ArgK","CAD","Top1","Ubx")
in.path<-"/Users/matthewnelsen/Documents/papers_reviews/papers/ant_plant_interactions/ants_june_to_august_2016/antcat_distribution_4_june_2016/outgroup_things/"
out.path<-"/Users/matthewnelsen/Documents/papers_reviews/papers/ants_for_doug/17nov2019_1/Formicidae/"
for(locus in loci){
	align<-read.fasta(file=paste(in.path,locus,"/",locus,"_Parsed_Renamed.fasta",sep=""),seqtype="DNA",forceDNAtolower=FALSE)
	out.align<-align[!grepl("Formicidae",names(align))]
	names(out.align)<-gsub("XM_[0-9]*","XM",names(out.align))
	names(out.align)<-gsub("NM_[0-9]*","NM",names(out.align))	
	for(q in 1:length(out.align)){
		names(out.align)[q]<-paste(strsplit(names(out.align)[q],"_")[[1]][4:(length(strsplit(names(out.align)[q],"_")[[1]])-1)],collapse="_")
	}	
	write.fasta(sequences=out.align,names=names(out.align),file.out=paste(out.path,locus,"_OUTS_ONLY_Parsed_Renamed.fasta",sep=""),nbchar=1000000)
}
