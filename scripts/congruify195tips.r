library(phytools)
library(geiger)
library(plyr)
std.tree<-read.tree(file="/home/mpnelsen/timescaling_doug/Dryad_Supplementary_File_7_ML_TREE_treepl_185.tre")
ml.tree<-read.tree(file="/home/mpnelsen/timescaling_doug/RAxML_bestTree_195tips.result")
outtree.name<-"/home/mpnelsen/timescaling_doug/RAxML_bestTree.result_treepl_prime_195tips.tre"
treepl.filename<-"/home/mpnelsen/timescaling_doug/RAxML_bestTree.result_primefile_195tips.txt"
ml.tree.path<-"/home/mpnelsen/timescaling_doug/RAxML_bestTree_195tips.result"

initial.congruify.to.treepl<-function(reftree=NULL,target.tree=NULL,numsites=NULL,penalty=NULL,n.threads=1,outtree.name=NULL,treepl.filename=NULL,ml.tree.path=NULL){
		new.vector<-c()
		res<-congruify.phylo(reftree,target.tree)
		new.vector[1]<-paste("treefile = ",ml.tree.path,sep="")
		new.vector[2]<-paste("numsites = ",numsites,sep="")
		new.vector[nrow(res$calibrations)*3+3]<-paste("thorough")
		new.vector[nrow(res$calibrations)*3+4]<-paste("prime")
		new.vector[nrow(res$calibrations)*3+5]<-penalty
		new.vector[nrow(res$calibrations)*3+6]<-paste("nthreads =", n.threads ,sep="")
		new.vector[nrow(res$calibrations)*3+7]<-paste("outfile =", outtree.name ,sep="")
		comp<-matchNodes(reftree,res$reference,method="descendants")
		#if below is false, stop and print error, if true, go on
		if(!isTRUE(all(!is.na(comp)) && all(comp[,"tr1"]==comp[,"tr2"]))){
			print("PROBLEM w TREE")
		}
		if(isTRUE(all(!is.na(comp)) && all(comp[,"tr1"]==comp[,"tr2"]))){
			res$calibrations[c("MaxAge","MinAge")]<-sapply(res$calibrations[c("MaxAge","MinAge")],as.numeric)
			res$calibrations$TargetNodeNo<-NA
			for (i in 1:nrow(res$calibrations)){
				new.vector[i*3+0]<-paste("mrca","=",paste("Calibration_",i,sep=""), res$calibrations$taxonA[i],res$calibrations$taxonB[i],sep=" ")
				new.vector[i*3+1]<-paste("min","=", paste("Calibration_",i,sep=""), round(res$calibrations$MinAge[i],2),sep=" ")
				new.vector[i*3+2]<-paste("max","=", paste("Calibration_",i,sep=""), round(res$calibrations$MaxAge[i],2),sep=" ")
			}	
		}
		write(new.vector,treepl.filename)
}

initial.congruify.to.treepl(reftree=std.tree,target.tree=ml.tree,numsites=9196,penalty="log_pen",n.threads=6,outtree.name=outtree.name,treepl.filename=treepl.filename,ml.tree.path=ml.tree.path)




library(phyloch)
library(phytools)
library(geiger)
std.tree<-read.tree(file="/home/mpnelsen/timescaling_doug/Dryad_Supplementary_File_7_ML_TREE_treepl_185.tre")
ml.tree<-read.tree(file="/home/mpnelsen/timescaling_doug/RAxML_bestTree_195tips.result")
ml.tree.path<-"/home/mpnelsen/timescaling_doug/RAxML_bestTree_195tips.result"
outtree.name<-"/home/mpnelsen/timescaling_doug/RAxML_bestTree.result_treepl_CV_195tips.tre"

congruify.to.treepl.cv.maker<-function(filepath=NULL,reftree=NULL,target.tree=NULL,cviter=100,cvstart=0.000001,cvstop=1000,cvmultstep=0.1,cviterno=2,randomcviterno=10,numsites=NULL,penalty=NULL,n.threads=1,ml.tree.path=NULL,outtree.name=NULL){
		new.vector<-c()
		res<-congruify.phylo(reftree,target.tree)
		new.vector[1]<-paste("treefile = ",ml.tree.path,sep="")
		new.vector[2]<-paste("numsites = ",numsites,sep="")
		new.vector[nrow(res$calibrations)*3+3]<-paste("thorough")
		new.vector[nrow(res$calibrations)*3+4]<-penalty
		new.vector[nrow(res$calibrations)*3+5]<-paste("nthreads = ", n.threads ,sep="")
		new.vector[nrow(res$calibrations)*3+6]<-paste("cviter = ",cviterno,sep="")
		new.vector[nrow(res$calibrations)*3+7]<-paste("cvstart = ",cvstart,sep="")
		new.vector[nrow(res$calibrations)*3+8]<-paste("cvstop = ",cvstop,sep="")
		new.vector[nrow(res$calibrations)*3+9]<-paste("cvmultstep = ",cvmultstep,sep="")
		new.vector[nrow(res$calibrations)*3+10]<-paste("randomcv",sep="")
		new.vector[nrow(res$calibrations)*3+11]<-paste("randomcviter = ",randomcviterno,sep="")
		new.vector[nrow(res$calibrations)*3+12]<-paste("outfile = ",outtree.name,sep="")
		new.vector[nrow(res$calibrations)*3+13]<-paste("cvoutfile = ",outtree.name,"_cv.outfile.txt",sep="")
		new.vector[nrow(res$calibrations)*3+14]<-"opt = 2"
		new.vector[nrow(res$calibrations)*3+15]<-"moredetail"		
		new.vector[nrow(res$calibrations)*3+16]<-"optad = 2"
		new.vector[nrow(res$calibrations)*3+17]<-"moredetailad"
		new.vector[nrow(res$calibrations)*3+18]<-"optcvad = 5"
		comp<-matchNodes(reftree,res$reference,method="descendants")
		#if below is false, stop and print error, if true, go on
		if(!isTRUE(all(!is.na(comp)) && all(comp[,"tr1"]==comp[,"tr2"]))){
			print("PROBLEM w TREE")
		}
		if(isTRUE(all(!is.na(comp)) && all(comp[,"tr1"]==comp[,"tr2"]))){
			res$calibrations[c("MaxAge","MinAge")]<-sapply(res$calibrations[c("MaxAge","MinAge")],as.numeric)
			res$calibrations$TargetNodeNo<-NA
			for (i in 1:nrow(res$calibrations)){
				new.vector[i*3+0]<-paste("mrca","=",paste("Calibration_",i,sep=""), res$calibrations$taxonA[i],res$calibrations$taxonB[i],sep=" ")
				new.vector[i*3+1]<-paste("min","=", paste("Calibration_",i,sep=""), round(res$calibrations$MinAge[i],2),sep=" ")
				new.vector[i*3+2]<-paste("max","=", paste("Calibration_",i,sep=""), round(res$calibrations$MaxAge[i],2),sep=" ")
			}	
		}
		write(new.vector,paste("cv_configfile_","RAxML_bestTree_195tips.result",sep=""))
}



congruify.to.treepl.cv.maker(filepath="/home/mpnelsen/timescaling_doug/",reftree=std.tree,target.tree=ml.tree,cviter=100,cvstart=0.000001,cvstop=1000,cvmultstep=0.1,cviterno=2,randomcviterno=10,numsites=9196,penalty="log_pen",n.threads=6,outtree.name=outtree.name,ml.tree.path=ml.tree.path)




