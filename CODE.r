list.files()->fls
read.table(fls[1],header=T,sep="\t")->g3
read.table(fls[2],header=T,sep="\t")->g4
read.table(fls[3],header=T,sep="\t")->g42
read.table(fls[4],header=T,sep="\t")->g5
read.table(fls[5],header=T,sep="\t")->tcga
read.table("StageIsamples.txt",header=T,sep="\t")->s1
g3<-g3[g3$sampleID %in% s1$sampleID,]
g4<-g4[g4$sampleID %in% s1$sampleID,]
g42<-g42[g42$sampleID %in% s1$sampleID,]
g5<-g5[g5$sampleID %in% s1$sampleID,]
tcga<-tcga[tcga$sampleID %in% s1$sampleID,]

genes<-c("RBL1","MCM4","MCM2","CDC27","CCNE1","CDKN1A","WEE1")
indexs<-function(x,y){
	idx<-1:3
	for(i in 1:length(y)){
		idx<-c(idx,which(x[i]==colnames(y)))
	}
	return(idx)
}
cox.analysis<-function(genes,coeff,dataset){
	coeff<-as.numeric(coeff[,1])
	samples<-dataset$sampleID
	dataset<-dataset[,c(ncol(dataset),1:(ncol(dataset)-1))]
	genes<-c(genes, "OS", "OS_IND")
	idx<-c()
	for(i in 1:length(genes)){
		idx<-c(idx,which(colnames(dataset)==genes[i]))
	}
	workplace<-dataset[,idx]
	OS<-workplace$OS;OS_IND<-workplace$OS_IND;workplace$OS<-NULL;workplace$OS_IND<-NULL
	riskscore<-c()
	for(i in 1:nrow(workplace)){
		riskscore<-c(riskscore,sum(workplace[i,]*coeff))
	}
	medians <- median(as.numeric(riskscore))
	medianB <- as.numeric(as.numeric(riskscore)>=medians)
	# pValue
	ss <- survdiff(Surv(OS, OS_IND) ~ medianB)
	p_os <- 1-pchisq(ss$chisq, 1)
	sf <- survfit(Surv(OS, OS_IND) ~ medianB)
	workplace$sampleID<-samples
	workplace<-workplace[,c(ncol(workplace),1:(ncol(workplace)-1))]
	data.write<-cbind(workplace,OS, OS_IND, riskscore, medianB, stringsAsFactors=FALSE)
	data.write<-as.data.frame(data.write, stringsAsFactors=FALSE)
	return(data.write)
}

g3s<-g3[,indexs(genes,g3)]
g4s<-g4[,indexs(genes,g4)]
g42s<-g42[,indexs(genes,g42)]
g5s<-g5[,indexs(genes,g5)]
tcgas<-tcga[,indexs(genes,tcga)]
g3s$OS<-g3s$OS/30
tcgaw$OS<-tcgaw$OS*12
coxm<-coxph(Surv(OS, OS_IND) ~ ., tcgas[,-1])
b<-as.data.frame(summary(coxm)$coefficient, stringsAsFactors=FALSE)
par(mfrow=c(2,3))
tcgaw<-cox.analysis(genes,b,tcgas)
g3w<-cox.analysis(genes,b,g3s)
g4w<-cox.analysis(genes,b,g4s)
g42w<-cox.analysis(genes,b,g42s)
g5w<-cox.analysis(genes,b,g5s)

setwd('../biomarker/NSCLC/')
write.table(tcgaw,"NSCLC_SevenGenes.txt",sep="\t",quote=F,row.names=F)
write.table(g3w,"GSE30219_SevenGenes.txt",sep="\t",quote=F,row.names=F)
write.table(g4w,"GSE41271_SevenGenes.txt",sep="\t",quote=F,row.names=F)
write.table(g42w,"GSE42127_SevenGenes.txt",sep="\t",quote=F,row.names=F)
write.table(g5w,"GSE50081_SevenGenes.txt",sep="\t",quote=F,row.names=F)

par(mar=c(7,4,4,4))
plot(c(),xlim=c(0,12),ylim=c(3,5.5),xlab=NA,ylab="Risk score",axes=F)
boxplot(m$riskscore~m$age,add=T,col=2:3,at=1:2,names=NA,pch=16,axes=F)
boxplot(m$riskscore~m$gender,add=T,col=2:3,at=4:5,names=NA,pch=16,axes=F)
boxplot(m$riskscore~m$iDimenson,add=T,col=2:3,at=7:8,names=NA,pch=16,axes=F)
boxplot(m$riskscore~m$Stage,add=T,col=2:3,at=10:11,names=NA,pch=16,axes=F)
axis(1,at=c(1,2,4,5,7,8,10,11),NA)
axis(2)
text(c(1,2,4,5,7,8,10,11),2.7,c("Age<60","Age>=60","Female","Male","Diameter>=1cm","Diameter<1cm","Stage 1-2","Stage 3-4"),xpd=T,pos=1,srt=45,adj=.5)
text(1.5,5.4,"p>0.01",col=2)
text(4.5,5.4,"p>0.01",col=2)
text(7.5,5.4,"p>0.01",col=2)
text(10.5,5.4,"p>0.01",col=2)

 



plot(dss,mark.time=T,lwd=2,col=1:2,xlab="Disease-free survival months",ylab="Survival rate",xlim=c(0,100))
legend("bottomleft",c("Low-risk","High-risk"),col=1:2,lwd=2)
text(20,.4,"p=1e-5")


plot(pfs,mark.time=T,lwd=2,col=1:2,xlab="Progression-free survival months",ylab="Survival rate",xlim=c(0,100))
legend("bottomleft",c("Low-risk","High-risk"),col=1:2,lwd=2)
text(80,.9,"p=4e-6")













