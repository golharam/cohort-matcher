install.packages("igraph")
library(igraph)
library(reshape2)
library(plyr)
#Clean
rm(list = ls())
snpt<-150
snpot<-0.7
#Read the edges that have the SNPs overlap and SNP num
cm<-read.table("meltedResults.txt",sep="\t",head=TRUE)
cm$weighted<-cm$Fraction_Match*log(cm$SNPs_Compared)
#get an idea about the set
cm.sorted<-cm[with(cm,order(Fraction_Match)),]
cm.sorted$color<-'black'
cm.sorted$color[cm.sorted$SNPs_Compared<100]<-'red'
cm.sorted$color[cm.sorted$SNPs_Compared<150&cm.sorted$SNPs_Compared>99]<-'orange'
png("snpdistr.png")
plot(cm.sorted$Fraction_Match,col=cm.sorted$color,main="SNP
overlap, red is <100 SNPs and orange <150 SNPs")
abline(0.7,0,col="red")
abline(0.6,0,col="orange")
dev.off()
png("weightDistr.png")
cm.sorted<-cm[with(cm,order(weighted)),]
cm.sorted$color<-'black'
cm.sorted$color[cm.sorted$SNPs_Compared<100]<-'red'
cm.sorted$color[cm.sorted$SNPs_Compared<150&cm.sorted$SNPs_Compared>99]<-'orange'
plot(cm.sorted$weighted,col=cm.sorted$color,main="SNP
overlap, red is <100 SNPs and orange <150 SNPs")
abline(0.7,0,col="red")
abline(0.6,0,col="orange")
dev.off()
#Filter, should be args?
#cm.graph<-cm[cm$Fraction_Match>snpot&cm$SNPs_Compared>snpt,c('Sample1','Sample2')]
cm.graph<-cm[cm$weighted>4,c('Sample1','Sample2')]
#Get the graph
cm.g<-graph_from_edgelist(as.matrix(cm.graph),directed=FALSE)
#Look only at the cliques with size >2 and do not report
sub-cliques
cm.maxcliques<-max_cliques(cm.g,min=2)
#Kludgy, but works
capture.output(cm.maxcliques,file="maxcliques.txt")

#Define how to get patients per clique
checkIfClean<-function(samples) {
patientPerClique<-unique(patients[patients$Sample %in%
samples,'Patient'])
if (length(patientPerClique)>1) { cliqt<-"mixed"} else
{cliqt<-"clean"}
return(c(cliqt,as.character(patientPerClique)))
}

#Get the patient info
patients<-read.table("patientSheet.txt",sep="\t",head=TRUE)
patients.cnts<-count(patients,c('Patient'))
paste("We have
",length(patients.cnts[patients.cnts$freq==1,c('Patient')]),
" patients with single sample")
#get the patients per clique
cm.named<-lapply(seq_along(cm.maxcliques),function(i) {
paste(checkIfClean(attributes(cm.maxcliques[[i]])$names)) })

storePdata<-function(i) {
pd<-checkIfClean(attributes(cm.maxcliques[[i]])$names);
mattr1<-attributes(cm.maxcliques[[i]])
mattr1$type <- pd[1];
mattr1$patients <- tail(pd,-1);
return(mattr1);
}

#Set the patient vector and type as attributes on a clique
l2<-lapply(seq_along(cm.maxcliques),function (i) {
attributes(cm.maxcliques[[i]])<-storePdata(i) })

#get the clean cliques to check them for overlaps
getClean<-function(i) {
return(c(i,l2[[i]]$type,l2[[i]]$patients))
}

clean.cliques<-do.call(rbind, lapply(seq_along(l2),
function(i)if (l2[[i]]$type=='clean') {c(l2[[i]]$patients,i)}))

#here are the cliques we should combine
cliques.grp<-dcast(as.data.frame(clean.cliques),V1 ~ V2)
#Now we should combine if a row has more than one non-NA value
cl.tomerge<-cliques.grp[apply(cliques.grp,1,function(x)
sum(!is.na(x))>2),]
mergeCliques<-function (mergeRow) {
clnums<-as.numeric(tail(mergeRow[!is.na(mergeRow)]))
allnodes<-do.call('rbind',lapply(cm.maxcliques[as.numeric(tail(clnums,-1))],function(x)
attributes(x)$names))
allnodes<-unique(as.vector(allnodes))
return(c(mergeRow[1],allnodes))
}
cm.mlist<-apply(cl.tomerge,1,mergeCliques)#Merged cliques
#We now have list of merged clean cliques, lets get the
simple cliques as well
cl.nomerge<-cliques.grp[apply(cliques.grp,1,function(x)
sum(!is.na(x))==2),]
cm.nomlist<-apply(cl.nomerge,1,mergeCliques)#Merged cliques
paste("We have ",length(cm.nomlist)+length(cm.mlist)," clean
cliques")
#Finally list of mixed(dirty cliques)
mixed.cliques<- lapply(seq_along(l2), function(i)if
(l2[[i]]$type=='mixed') {c(l2[[i]]$patients,i)})
mixed.cliques<-mixed.cliques[-which(sapply(mixed.cliques,
is.null))]
paste("We have also ",length(mixed.cliques)," cases with
potential swaps")

#Mixed maxcliques
mixed.maxcliques<-cm.maxcliques[as.numeric(unlist(lapply(mixed.cliques,function(x)
{tail(x,1)})))]
write(unlist(lapply(seq_along(mixed.maxcliques),function(i){
paste(unlist(attributes(mixed.maxcliques[[i]])$names),collapse="
")})),file="mixed.cliques.txt")

#OK samples
clean.samples<-sort(unlist(lapply(cm.nomlist,function(x)
tail(x,-1))))
clean.samples<-append(clean.samples,sort(unlist(lapply(cm.mlist,function(x)
tail(x,-1)))))
paste("We have ",length(clean.samples)," clean samples")
write(clean.samples,file="clean.samples.txt")
