#graft doug's strumigenys tree on


require(phytools)
tr1<-read.tree(file="/Users/matthewnelsen/Documents/papers_reviews/papers/ants_for_doug/17nov2019_1/Formicidae/alignments/Formicidae/doug_partition_finder/2020_second_run_new_partitions_195tips_2/RAxML_bestTree.result_treepl_CV_195tips.tre")
tr2<-read.tree(file="/Users/matthewnelsen/Documents/papers_reviews/papers/ants_for_doug/StrumigenysTree.nexus.tre")
#need to do this: http://www.phytools.org/static.help/paste.tree.html
tr2b<-tr2
tr2b$root.edge<-0


#find strumigenys
str<-tr1$tip.label[grep("Strumigenys",tr1$tip.label)]

tr1b<-drop.tip(tr1,str[1:(length(str)-1)])
tr1b$tip.label[tr1b$tip.label %in% "Strumigenys_membranifera"]<-"NA"
plot(tr1b,cex=0.4)
edgelabels(tr1b$edge.length, col="black", cex=0.5)
tr1b$edge.length



ht<-max(nodeHeights(tr2))

mod<-paste.tree(tr1b,tr2b)
plot(mod,cex=0.5)
edgelabels(mod$edge.length, col="black", cex=0.5)

tr1b$edge.length[!tr1b$edge.length %in% mod$edge.length]

#81.870001 #171
mod$edge.length[mod$edge.length>81.8 & mod$edge.length<81.9]


mod$edge.length[mod$edge.length==81.870001]<-81.870001-ht


plot(mod,cex=0.4)
edgelabels(mod$edge.length, col="black", cex=0.5)

is.ultrametric(mod)

plot(ladderize(mod,FALSE),cex=0.7)

mod.lad<-ladderize(mod,FALSE)
write.tree(mod.lad,file="/Users/matthewnelsen/Documents/papers_reviews/papers/ants_for_doug/17nov2019_1/Formicidae/alignments/Formicidae/doug_partition_finder/2020_second_run_new_partitions_195tips_2/RAxML_bestTree.result_treepl_CV_195tips_STRUMIGENYSADDED.tre")

mod.lad$tip.label
outs<-c("Apis_mellifera","Pristocera_MAD01","Chyphotes_mellipes","Chyphotes_sp","Philanthus_sp","Dasymutilla_aureola","Odontophotopsis_sp","Aporus_niger","Evagetes_sp","Sapyga_pumila","Scolia_verticalis","Chalybion_californicum","Aglyptacros_cf_sulcatus","Metapolybia_cingulata","Mischocyttarus_flavitarsis","Vespula_sp")
mod.lad.in<-drop.tip(mod.lad,outs)
plot(mod.lad.in,cex=0.5)
axisPhylo()

write.tree(mod.lad.in,file="/Users/matthewnelsen/Documents/papers_reviews/papers/ants_for_doug/17nov2019_1/Formicidae/alignments/Formicidae/doug_partition_finder/2020_second_run_new_partitions_195tips_2/RAxML_bestTree.result_treepl_CV_195tips_STRUMIGENYSADDED_INGROUP.tre")
