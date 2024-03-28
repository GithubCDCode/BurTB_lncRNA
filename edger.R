library("edgeR")#3.24.3   
library('gplots')
foldchange=1
padj=0.05
rt=read.csv("gene_count_matrix.csv",header = TRUE,row.names = 1,check.names = FALSE)
exp=as.matrix(rt)
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames = dimnames)
#data=data[rowMeans(data)>1,]
data=data[rowSums(data)>1,]

group <- read.csv("ground.txt",header=TRUE,row.names = 1,check.names = FALSE)
group <- group[,1]
design <- model.matrix(~group)

y <- DGEList(counts=data,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y,pair = c("CR","XR"))
topTags(et)
summary(de <- decideTestsDGE(et))
#ordered_tags <- topTags(et, n=100000)
ordered_tags <- topTags(et,n=NULL,sort.by='none')
allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff
newData=y$pseudo.counts

write.csv(diff, "CRNF.csv")
diffSig = diff[(diff$FDR < padj & (diff$logFC>foldchange | diff$logFC<(-foldchange))),]
write.table(diffSig, file="CRNFdiffSig.xls",sep="\t",quote=F)
write.csv(diffSig, "CRNFdiffSig.csv")
diffUp = diff[(diff$FDR < padj & (diff$logFC>foldchange)),]
write.table(diffUp, file="CRNFup.xls",sep="\t",quote=F)
write.csv(diffUp, "CRNFdiffUp.csv")
diffDown = diff[(diff$FDR < padj & (diff$logFC<(-foldchange))),]
write.table(diffDown, file="CRNFdown.xls",sep="\t",quote=F)
write.csv(diffDown, "CRNFdiffDown.csv")
