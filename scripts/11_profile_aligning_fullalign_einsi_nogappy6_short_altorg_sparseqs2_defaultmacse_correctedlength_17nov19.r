profile.align.tree.list.maker<-function(to.make,levels.for.tree,clade.level,locus,input.accfile){
	require(plyr)
	require(ape)
	input.accfile<-read.csv(input.accfile,stringsAsFactors=FALSE)
	mylist<-subset(input.accfile,select=c(levels.for.tree,locus))
	newlist<-mylist[!is.na(mylist[,locus]),]
	newlist<-as.data.frame(unclass(newlist))
	seqs.per.grp<-as.data.frame(table(newlist[clade.level]))
	colnames(seqs.per.grp)<-c(as.character(clade.level),"Abund")
	newlist<-join(newlist,seqs.per.grp,by=as.character(clade.level))
	newlist.b<-newlist
	newlist.b[,locus]<-NULL
	tree.dat<-unique(newlist.b)
	tax.hier<-paste("~",gsub(" ","/",paste(levels.for.tree,collapse=" ")),sep="")
	tr<-as.phylo(as.formula(tax.hier),data=tree.dat)
	tr.ord<-reorder(tr,"postorder")
	org<-as.data.frame(tr.ord$edge)
	colnames(org)<-c("From","To")
	nn<-unique(org$From)
	nn.fr<-plyr::count(org$From)
	colnames(nn.fr)<-c("Node.No","Desc")
	nn.fr<-nn.fr[match(nn.fr$Node.No,nn),]
	org$To.Name<-tr.ord$tip.label[tr.ord$edge[,2]]
	org$To.Name[is.na(org$To.Name)]<-org$To[is.na(org$To.Name)]
	for(tax in 1:nrow(org)){
		org$Abund[tax]<-newlist$Abund[newlist[clade.level]==org$To.Name[tax]][1]
	}
	if(to.make=="tree"){
		return(tr.ord)
	}	
	if(to.make=="large.summary"){
		return(org)
	}	
	if(to.make=="par.desc"){
		return(nn.fr)
	}
}


profile.aligner.merge<-function(clade,clade.level,locus,organization,node.nos.fr,folderpath,path.to.mafft,dat.type,threads.mafft,tax.stand=NULL){
	require(ape)
	require(seqinr)
	for(z in 1:nrow(node.nos.fr)){
		nds<-organization$To.Name[organization$From==node.nos.fr$Node.No[z]]
		#cat.files<-as.vector(paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",nds,"_",locus,"_aligned.fasta",sep=""))
		cat.files<-as.character(paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",nds,"_",locus,"_aligned.fasta",sep="", collapse=" "))
		cat.outfile<-as.character(paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",node.nos.fr$Node.No[z],"_",locus,"_cat_aligned.fasta",sep=""))
		cat.commands<-as.character(paste("cat",paste(cat.files,sep=" "),">",cat.outfile,sep=" "))
		system(cat.commands,intern=TRUE)
		ruby.table.maker.file<-as.character(paste(folderpath,"makemergetable.rb",sep=""))
		sub.msa.table<-as.character(paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",node.nos.fr$Node.No[z],"_",locus,"_submsa_table",sep=""))
		ruby.commands<-as.character(paste("ruby",ruby.table.maker.file, cat.files, ">",sub.msa.table,sep=" "))
		system(ruby.commands,intern=TRUE)
		output.file<-as.character(paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",node.nos.fr$Node.No[z],"_",locus,"_aligned.fasta",sep=""))	
		#mafft.commands<-paste(path.to.mafft,"--localpair --maxiterate 1000 --preservecase --merge", sub.msa.table, cat.outfile, ">",output.file, sep=" ")
		#make something here to get number of seqs going into alignment
		#print(cat.files)
		#this is problematic if both files are joined together like this "file1 file2" instead of "file1" "file2"
		#ntax<-sum(sapply(cat.files,function(zaz) length(read.fasta(file=zaz,seqtype=dat.type[ind.loc],forceDNAtolower=FALSE))))
		#so instead count number in cat.outfile
		if(dat.type!="AA"){
			ntax<-length(read.fasta(file=cat.outfile,seqtype="DNA",forceDNAtolower=FALSE))	
		}
		if(dat.type=="AA"){
			ntax<-length(read.fasta(file=cat.outfile,seqtype="AA",forceDNAtolower=FALSE))				
		}
		#then make small decision tree for mafft commands - if less than 200 seqs, use G-INS-i (and --leavegappyregion), but if over, use FFT-NS-i (w/o --leavegappyregion)
		#--leavegappyregion only works well in small alignments (Katoh & Standley 2016), so remove from large aligns
		#does not seem wise to set unalignlevel (VSM) higher than 0.8 (Katoh & Standley 2016).
		if(ntax<200){
			#mafft.commands<-paste(path.to.mafft,"--globalpair --maxiterate 1000 --thread", threads, "--preservecase --merge", sub.msa.table, cat.outfile, ">",output.file, sep=" ")
			#did G-INS-i instead
			#mafft.commands<-paste(path.to.mafft,"--globalpair --maxiterate 1000 --thread", threads, "--preservecase --unalignlevel 0.8 --leavegappyregion --merge", sub.msa.table, cat.outfile, ">", output.file, sep=" ")	
			#did E-INS-i instead
			mafft.commands<-paste(path.to.mafft,"--genafpair --maxiterate 1000 --thread", threads.mafft, "--preservecase --merge", sub.msa.table, cat.outfile, ">", output.file, sep=" ")	
		}
		if(ntax>=200){
			#did FFT-NS-i...cannot add unalignlevel as this is only available when using G-INS-i
			mafft.commands<-paste(path.to.mafft,"--retree 2 --maxiterate 1000 --thread", threads.mafft, "--preservecase --merge", sub.msa.table, cat.outfile, ">",output.file, sep=" ")	
		}
		system(mafft.commands,intern=FALSE)
		Sys.sleep(1)
		#Refine alignment w MACSE
		#if(dat.type=="codon"){
		#	output.file<-as.character(paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",node.nos.fr$Node.No[z],"_",locus,"_aligned.fasta",sep=""))	
		#	output.file.unrefined<-as.character(paste(strsplit(output.file,split=".fasta")[[1]],"_unrefined.fasta",sep=""))
		#	copy.commands.first<-paste("cp", output.file, output.file.unrefined, sep=" ")
		#	system(copy.commands.first)
		#	commands.to.macse<-paste("java -jar -Xmx1200m", path.to.macse,"-prog refineAlignment -align", output.file.unrefined, "-gc_def", translation.cat, "-out_AA", paste(output.file,"_AA.fasta",sep=""), "-out_NT", paste(output.file,"_NT.fasta",sep=""), "-optim 1 -ext_gap_ratio 0.0001 -gap_op 1", sep=" ")
		#	system(commands.to.macse)
		#	copy.commands<-paste("cp", as.character(paste(output.file,"_NT.fasta",sep="")), output.file, sep=" ")
		#	system(copy.commands)
		#	commands.to.perl<-paste("perl -i -pe 's/[!]/-/g'", output.file, sep=" ")
		#	system(commands.to.perl)
		#	Sys.sleep(1)
		#}	
	}
	rename.commands<-paste("cp", paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",node.nos.fr$Node.No[nrow(node.nos.fr)],"_",locus,"_aligned.fasta",sep=""), paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",locus,"_aligned.fasta",sep=""), sep=" ")
	system(rename.commands,intern=TRUE)
	#Remove taxon standards if included previously
	if(!is.null(tax.stand)){	
		if(dat.type!="AA"){
			mat<-read.fasta(file=paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",locus,"_aligned.fasta",sep=""),seqtype="DNA",forceDNAtolower=FALSE)
			mat<-mat[-grep("Standard_.*",names(mat))]
			write.fasta(sequences=mat,names=names(mat),file=paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",locus,"_aligned.fasta",sep=""),nbchar=1000000)
		}
		if(dat.type=="AA"){
			mat<-read.fasta(file=paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",locus,"_aligned.fasta",sep=""),seqtype="AA",forceDNAtolower=FALSE)		
			mat<-mat[-grep("Standard_.*",names(mat))]
			write.fasta(sequences=mat,names=names(mat),file=paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",locus,"_aligned.fasta",sep=""),nbchar=1000000)
		}
	}
}

macse.refine<-function(folderpath,clade,locus,clade.level,path.to.macse,translation.cat){
	#Refine alignment w MACSE
	output.file<-as.character(paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",locus,"_aligned.fasta",sep=""))	
	output.file.unrefined<-as.character(paste(strsplit(output.file,split=".fasta")[[1]],"_unrefined.fasta",sep=""))
	copy.commands.first<-paste("cp", output.file, output.file.unrefined, sep=" ")
	system(copy.commands.first)
	commands.to.macse<-paste("java -jar -Xmx1200m", path.to.macse,"-prog refineAlignment -align", output.file, "-gc_def", translation.cat, "-out_AA", paste(output.file,"_AA.fasta",sep=""), "-out_NT", paste(output.file,"_NT.fasta",sep=""), "-optim 1", sep=" ")
	system(commands.to.macse)
	copy.commands<-paste("cp", as.character(paste(output.file,"_NT.fasta",sep="")), output.file, sep=" ")
	system(copy.commands)
	commands.to.perl<-paste("perl -i -pe 's/[!]/-/g'", output.file, sep=" ")
	system(commands.to.perl)
	Sys.sleep(1)	
}

gblocks.remove<-function(dat.type,folderpath,clade,locus,clade.level){
	#added condition that gblocks only used if alignment specified
	if(dat.type!="AA"){
		full.align<-read.fasta(file=paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",locus,"_aligned.fasta",sep=""),seqtype="DNA",forceDNAtolower=FALSE)
	}
	if(dat.type=="AA"){
		full.align<-read.fasta(file=paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",locus,"_aligned.fasta",sep=""),seqtype="AA",forceDNAtolower=FALSE)
	}
	ntax<-length(full.align)
	if(dat.type=="DNA"){
		t.var<-as.character(paste("-t=","d",sep=""))
	}
	if(dat.type=="AA"){
		t.var<-as.character(paste("-t=","p",sep=""))
	}
	if(dat.type=="codon"){
		t.var<-as.character(paste("-t=","c",sep=""))
	}
	subdat<-ceiling((ntax/2)+0.5)
	b2.var<-as.character(paste("-b2=",subdat,sep=""))
	b4.var<-as.character(paste("-b4=","5",sep=""))	
	b5.var<-as.character(paste("-b5=","h",sep=""))
	file.extension<-as.character(paste("-e=","-gb5",sep=""))	
	gblocks.commands<-paste(path.to.gblocks, paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",locus,"_aligned.fasta",sep=""), t.var, b2.var, b4.var, b5.var, file.extension, sep=" ")
	print(gblocks.commands)
	system(gblocks.commands,intern=TRUE)
	Sys.sleep(1)
	b4.var<-as.character(paste("-b4=","2",sep=""))
	file.extension<-as.character(paste("-e=","-gb2",sep=""))	
	gblocks.commands<-paste(path.to.gblocks, paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",locus,"_aligned.fasta",sep=""), t.var, b2.var, b4.var, b5.var, file.extension, sep=" ")
	cat(gblocks.commands)
	system(gblocks.commands,intern=TRUE)
	Sys.sleep(1)
	#to.rename<-read.fasta(file=as.character(paste(output.file,"-gb",sep="")),seqtype=dat.type[ind.loc],forceDNAtolower=FALSE)
	#write.fasta(sequences=to.rename,names=names(to.rename),file=output.file,nbchar=1000000)
	#rename.commands<-paste("cp", paste(output.file,"-gb",sep=""), output.file, sep=" ")
	#system(rename.commands,intern=TRUE)
	#Sys.sleep(1)
}

profile.to.end.single<-function(levels.for.tree=tree.levs,locus=locus,input.accfile=input.accfile,tax.stand=tax.stand,clade=clade,clade.level=clade.level,organization=ind.org,node.nos.fr=ind.node.no.fr,folderpath=folderpath,path.to.mafft=path.to.mafft, threads.mafft=threads.mafft,dat.type=dat.type,translation.cat=translation.cat,seq.refs=seq.refs){
	ind.org<-profile.align.tree.list.maker(to.make="large.summary",levels.for.tree=tree.levs,clade.level=clade.level,locus=locus,input.accfile=input.accfile)
	ind.node.no.fr<-profile.align.tree.list.maker(to.make="par.desc",levels.for.tree=tree.levs,clade.level=clade.level,locus=locus,input.accfile=input.accfile)
	profile.aligner.merge(clade=clade,clade.level=clade.level,locus=locus,dat.type=dat.type,organization=ind.org,node.nos.fr=ind.node.no.fr,folderpath=folderpath,path.to.mafft=path.to.mafft,threads.mafft=threads.mafft,tax.stand=tax.stand)
}

gaps.to.qs<-function (x, leading=TRUE, trailing=TRUE){
	#this is modified from the trimSpace function of seqinr 3.1.3
	require(seqinr)
    for(q in 1:length(x)){
    	if(leading){
        	leading.pattern<-regmatches(x[[q]][1],regexpr("^[-]*",x[[q]][1]))
        	lead.no<-nchar(leading.pattern)
        	rep.pattern<-gsub("-","?",leading.pattern)
        	x[[q]][1]<-sub(pattern=paste("^[-]{",lead.no,"}",sep=""), replacement=rep.pattern, x=x[[q]][1], perl=TRUE)
    	}
    	if(trailing){
        	trailing.pattern<-regmatches(x[[q]][1],regexpr("[-]*$",x[[q]][1]))
        	trail.no<-nchar(trailing.pattern)
 			rep.pattern<-gsub("-","?",trailing.pattern)
        	x[[q]][1]<-sub(pattern=paste("[-]{",trail.no,"}","$",sep=""), replacement=rep.pattern, x=x[[q]][1], perl=TRUE)
    	}
    }
    return(x)
}

loci<-c("nuSSU","nuLSU","AbdA","COI","LR","Wg","EF1aF1","EF1aF2","ArgK","CAD","Top1","Ubx")
dat.type<-c("DNA","DNA","codon","codon","codon","codon","codon","codon","codon","codon","codon","codon")
translation.cat<-c(0,0,1,5,1,1,1,1,1,1,1,1)
folderpath<-"/Users/matthewnelsen/Documents/papers_reviews/papers/ants_for_doug/17nov2019_1/Formicidae/alignments/"
clade="Formicidae"
clade.level<-"Tribe"
#mafft v7.305b 2016/Aug/16
path.to.mafft<-"/Users/matthewnelsen/mafft"
path.to.macse<-"/Applications/macse_v1.2.jar"
tree.levs<-c("Family","Subfamily","Tribe")
path.to.gblocks<-"/Applications/Gblocks_0.91b/Gblocks"
threads<-2
#path.to.raxml<-"/farmshare/user_data/mpnelsen/standard-RAxML/raxmlHPC-PTHREADS-SSE3"
input.accfile<-as.character(paste(folderpath,clade,"/","combined.updated.updated.pruned.wouts.csv",sep=""))
seq.refs<-read.csv(file=as.character(input.accfile),stringsAsFactors=FALSE)
#potential.outgroups<-c(seq.refs$Taxon[seq.refs$Family=="Bethylidae"])
taxon.standards<-c("Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera")
out.accfile<-as.character(paste(folderpath,clade,"/","combined.updated.updated.pruned.wouts.gblocksremoved.csv",sep=""))
seq.type<-c("DNA","DNA","DNA","DNA","DNA","DNA","DNA","DNA","DNA","DNA","DNA","DNA")

if(any(dat.type!="codon")){
	dat.type.non.codon<-dat.type[dat.type!="codon"]
	loci.non.codon<-loci[dat.type!="codon"]
	if(!missing(taxon.standards)){
		taxon.standards.non.codon<-taxon.standards[dat.type!="codon"]
	}
	if(missing(taxon.standards)){
		taxon.standards.non.codon<-NULL		
	}
	for(ind.loc in 1:length(loci.non.codon)){
		profile.to.end.single(levels.for.tree=tree.levs,locus=loci.non.codon[ind.loc],tax.stand=taxon.standards.non.codon[ind.loc],input.accfile=input.accfile,clade=clade,clade.level=clade.level,organization=ind.org,node.nos.fr=ind.node.no.fr,folderpath=folderpath,path.to.mafft=path.to.mafft, threads.mafft=threads,dat.type=dat.type.non.codon[ind.loc],seq.refs=seq.refs)
		#replace empty w question marks before gblocks
		f.al<-read.fasta(file=paste(folderpath,clade,"/",loci.non.codon[ind.loc],"/",clade.level,"_Matrices","/",loci.non.codon[ind.loc],"_aligned.fasta",sep=""),seqtype=dat.type.non.codon[ind.loc],forceDNAtolower=FALSE,as.string=TRUE)
		f.al.q<-gaps.to.qs(f.al)
		write.fasta(sequences=f.al.q,names=names(f.al.q),file.out=paste(folderpath,clade,"/",loci.non.codon[ind.loc],"/",clade.level,"_Matrices","/",loci.non.codon[ind.loc],"_aligned.fasta",sep=""),nbchar=1000000)
		gblocks.remove(dat.type=dat.type.non.codon[ind.loc],folderpath=folderpath,clade=clade,locus=loci.non.codon[ind.loc],clade.level=clade.level)
	}
}

if(any(dat.type=="codon")){
	require(doParallel)
	require(seqinr)
	dat.type.codon<-dat.type[dat.type=="codon"]
	loci.codon<-loci[dat.type=="codon"]
	translation.cat.codon<-translation.cat[dat.type=="codon"]
	if(!missing(taxon.standards)){
		taxon.standards.codon<-taxon.standards[dat.type=="codon"]
	}
	if(missing(taxon.standards)){
		taxon.standards.codon<-NULL		
	}
	for(qv in 1:length(loci.codon)){
		profile.to.end.single(levels.for.tree=tree.levs,locus=loci.codon[qv],tax.stand=taxon.standards.codon[qv],input.accfile=input.accfile,clade=clade,clade.level=clade.level,organization=ind.org,node.nos.fr=ind.node.no.fr,folderpath=folderpath,path.to.mafft=path.to.mafft, threads.mafft=threads,dat.type=dat.type.codon[qv],seq.refs=seq.refs)
	}
	if(length(loci.codon)<threads){
		cl<-makeCluster(length(loci.codon))
	}
	if(length(loci.codon)>=threads){
		cl<-makeCluster(threads)
	}
	length(cl)
	registerDoParallel(cl)
	foreach(qa=1:length(loci.codon)) %dopar% macse.refine(folderpath=folderpath,clade=clade,locus=loci.codon[qa],clade.level=clade.level,path.to.macse=path.to.macse,translation.cat=translation.cat.codon[qa])
	for(zz in 1:length(loci.codon)){
		#replace empty w question marks
		f.al<-read.fasta(file=paste(folderpath,clade,"/",loci.codon[zz],"/",clade.level,"_Matrices","/",loci.codon[zz],"_aligned.fasta",sep=""),seqtype="DNA",forceDNAtolower=FALSE,as.string=TRUE)
		f.al.q<-gaps.to.qs(f.al)
		write.fasta(sequences=f.al.q,names=names(f.al.q),file.out=paste(folderpath,clade,"/",loci.codon[zz],"/",clade.level,"_Matrices","/",loci.codon[zz],"_aligned.fasta",sep=""),nbchar=1000000)
		gblocks.remove(dat.type=dat.type.codon[zz],folderpath=folderpath,clade=clade,locus=loci.codon[zz],clade.level=clade.level)
	}	
}
stopCluster(cl)

#this will then remove sequences from alignments with small proportion of sequence remaining after gblocks removal and update spreadsheet
for (l in 1:length(loci)){
	require(seqinr)
	require(plyr)
	if(seq.type[l]=="AA"){
		big.matrix<-read.fasta(file=paste(folderpath,clade,"/",loci[l],"/",clade.level,"_Matrices","/",loci[l],"_aligned.fasta-gb2",sep=""),seqtype="AA",forceDNAtolower=FALSE)
	}
	if(seq.type[l]!="AA"){
		big.matrix<-read.fasta(file=paste(folderpath,clade,"/",loci[l],"/",clade.level,"_Matrices","/",loci[l],"_aligned.fasta-gb2",sep=""),seqtype="DNA",forceDNAtolower=FALSE)
	}
	#this extracts the proportion of the sequence that is NOT alphabetic
	len.mat<-data.frame(matrix(nrow=length(big.matrix),ncol=4))
	colnames(len.mat)<-c("Taxon","Length","NumberMissing","PropMissing")
	for(qq in 1:length(big.matrix)){
		ga<-plyr::count(big.matrix[qq],1)
		len.mat[qq,1]<-colnames(ga)[1]
		len.mat[qq,2]<-sum(ga$freq)
		len.mat[qq,3]<-sum(ga$freq[!ga[,1] %in% LETTERS])
		len.mat[qq,4]<-sum(ga$freq[!ga[,1] %in% LETTERS])/sum(ga$freq)
	}
	#sparse<-len.mat$Taxon[len.mat$PropMissing>0.25]
	#big.matrix<-big.matrix[!names(big.matrix) %in% sparse]
	write.fasta(sequences=big.matrix,names=names(big.matrix),file.out=paste(folderpath,clade,"/",loci[l],"/",clade.level,"_Matrices","/",loci[l],"_aligned.fasta-gb2_removed.length.fasta",sep=""))
	write.csv(seq.refs[seq.refs$Taxon %in% sparse,c("Taxon",loci[l])],file=paste(folderpath,clade,"/",loci[l],"/",clade.level,"_Matrices","/",loci[l],"_removed.length.csv",sep=""),quote=FALSE,row.names=FALSE)
	seq.refs[seq.refs$Taxon %in% sparse,loci[l]]<-NA
	seq.refs<-seq.refs[rowSums(is.na(seq.refs[,loci]))!=length(loci),]
}
#write.csv(seq.refs,file=out.accfile,quote=FALSE,row.names=FALSE)
write.csv(seq.refs,file=out.accfile,row.names=FALSE)













#####################
#####################
#####################
#####################








#I don't like the removal of the short sequences, but did notice some accessions that did not align well and would like to remove them.
removes<-read.csv(file="~/Documents/papers_reviews/papers/ants_for_doug/24sep2018_flga_ants/Formicidae/alignments/Formicidae/accessions_to_remove.csv",stringsAsFactors=FALSE)

#first, remove from spreadsheet that includes short seqs
folderpath<-"/Users/matthewnelsen/Documents/papers_reviews/papers/ants_for_doug/17nov2019_1/Formicidae/alignments/"
clade="Formicidae"
input.accfile<-as.character(paste(folderpath,clade,"/","combined.updated.updated.pruned.wouts.csv",sep=""))
seq.refs<-read.csv(file=as.character(input.accfile),stringsAsFactors=FALSE)
for(r in 1:nrow(removes)){
	seq.refs[seq.refs$Taxon %in% removes$Taxon[r],removes$Locus[r]]<-NA
}
out.accfile<-as.character(paste(folderpath,clade,"/","combined.updated.updated.pruned.wouts.badalignsremoved.csv",sep=""))
write.csv(seq.refs,file=out.accfile,quote=FALSE,row.names=FALSE)

loci<-c("nuSSU","nuLSU","AbdA","COI","LR","Wg","EF1aF1","EF1aF2","ArgK","CAD","Top1","Ubx")
dat.type<-c("DNA","DNA","codon","codon","codon","codon","codon","codon","codon","codon","codon","codon")
seq.type<-c("DNA","DNA","DNA","DNA","DNA","DNA","DNA","DNA","DNA","DNA","DNA","DNA")
#this will then remove sequences from alignments with small proportion of sequence remaining after gblocks removal and update spreadsheet
for (l in 1:length(loci)){
	require(seqinr)
	require(plyr)
	if(seq.type[l]=="AA"){
		big.matrix<-read.fasta(file=paste(folderpath,clade,"/",loci[l],"/",clade.level,"_Matrices","/",loci[l],"_aligned.fasta-gb2",sep=""),seqtype="AA",forceDNAtolower=FALSE)
	}
	if(seq.type[l]!="AA"){
		big.matrix<-read.fasta(file=paste(folderpath,clade,"/",loci[l],"/",clade.level,"_Matrices","/",loci[l],"_aligned.fasta-gb2",sep=""),seqtype="DNA",forceDNAtolower=FALSE)
	}
	#remove taxa from individual alignments that are wonky
	if(loci[l] %in% removes$Locus){
		drop.tax<-removes$Taxon[removes$Locus %in% loci[l]]
		big.matrix<-big.matrix[!names(big.matrix) %in% drop.tax]
	}
	write.fasta(sequences=big.matrix,names=names(big.matrix),file.out=paste(folderpath,clade,"/",loci[l],"/",clade.level,"_Matrices","/",loci[l],"_aligned.fasta-gb2_badalignsremoved.fasta",sep=""))
}
