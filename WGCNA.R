rm(list=ls())
options(stringsAsFactors = F)
setwd('C:/Myproject/20230504_HW_lncRNA/')

exprset <- read.table('genes.readcount.txt',sep = '\t',header = T,row.names = 1,check.names = F)
brucella <- c(1,2,3,4,5,9,10)
exprset <- exprset[,brucella]
pdata <- c('control',rep('brucella',4),'control','control')
library(edgeR)
dge <- DGEList(counts = exprset,group = pdata)
dge <- calcNormFactors(dge)
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)
et <- exactTest(dge)
topTags(et)
order_tags <- topTags(et,n=100000)
order_tags <- data.frame(order_tags)
write.csv(order_tags,file="controlvsbrucella.mRNA.gene_level.deg.csv",row.names=T)

library(ggplot2)
data <- read.table('controlvsbrucella.mRNA.gene_level.deg.csv',sep = ',',header = T)
data$log2FoldChange
ggplot(data,aes(log2FoldChange,-log10(pvalue))) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed',color = '#999999') +
  geom_vline(xintercept = c(-1,1),linetype = 'dashed',color = '#999999')+
  geom_point(aes(size = -log10(pvalue),color = -log10(pvalue))) +
  scale_color_gradientn(values = seq(0,1,0.2),colors = c('#39489f','#39bbec','#f9ed36',"#f38466",'#b81f25')) +
  scale_size_continuous(range = c(1,3)) +
  theme_bw() +
  theme(panel.grid = element_blank())

exprset <- read.table('lnc_transcripts.readcount',sep = '\t',header = T,row.names = 1,check.names = F)
brucella <- c(1,2,3,4,5,9,10)
exprset <- exprset[,brucella]
pdata <- c('control',rep('brucella',4),'control','control')
library(edgeR)
dge <- DGEList(counts = exprset,group = pdata)
dge <- calcNormFactors(dge)
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)
et <- exactTest(dge)
topTags(et)
order_tags <- topTags(et,n=100000)
order_tags <- data.frame(order_tags)
write.csv(order_tags,file="controlvsbrucella.lncRNA.transcript_level.deg.csv",row.names=T)

library(ggplot2)
data <- read.table('controlvsbrucella.lncRNA.transcript_level.deg.csv',sep = ',',header = T)
data$log2FoldChange
ggplot(data,aes(log2FoldChange,-log10(pvalue))) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed',color = '#999999') +
  geom_vline(xintercept = c(-1,1),linetype = 'dashed',color = '#999999')+
  geom_point(aes(size = -log10(pvalue),color = -log10(pvalue))) +
  scale_color_gradientn(values = seq(0,1,0.2),colors = c('#39489f','#39bbec','#f9ed36',"#f38466",'#b81f25')) +
  scale_size_continuous(range = c(1,3)) +
  theme_bw() +
  theme(panel.grid = element_blank())

exprset <- read.table('genes.readcount.txt',sep = '\t',header = T,row.names = 1,check.names = F)
brucella <- c(1,6,7,8,9,10)
exprset <- exprset[,TB]
pdata <- c('control',rep('TB',3),'control','control')
library(edgeR)
dge <- DGEList(counts = exprset,group = pdata)
dge <- calcNormFactors(dge)
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)
et <- exactTest(dge)
topTags(et)
order_tags <- topTags(et,n=100000)
order_tags <- data.frame(order_tags)
write.csv(order_tags,file="controlvsTB.mRNA.gene_level.deg.csv",row.names=T)

library(ggplot2)
data <- read.table('controlvsTB.mRNA.gene_level.deg.csv',sep = ',',header = T)
data$log2FoldChange
ggplot(data,aes(log2FoldChange,-log10(pvalue))) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed',color = '#999999') +
  geom_vline(xintercept = c(-1,1),linetype = 'dashed',color = '#999999')+
  geom_point(aes(size = -log10(pvalue),color = -log10(pvalue))) +
  scale_color_gradientn(values = seq(0,1,0.2),colors = c('#39489f','#39bbec','#f9ed36',"#f38466",'#b81f25')) +
  scale_size_continuous(range = c(1,3)) +
  theme_bw() +
  theme(panel.grid = element_blank())

exprset <- read.table('lnc_transcripts.readcount',sep = '\t',header = T,row.names = 1,check.names = F)
brucella <- c(1,6,7,8,9,10)
exprset <- exprset[,TB]
pdata <- c('control',rep('TB',3),'control','control')
library(edgeR)
dge <- DGEList(counts = exprset,group = pdata)
dge <- calcNormFactors(dge)
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)
et <- exactTest(dge)
topTags(et)
order_tags <- topTags(et,n=100000)
order_tags <- data.frame(order_tags)
write.csv(order_tags,file="controlvsTB.lncRNA.transcript_level.deg.csv",row.names=T)

library(ggplot2)
data <- read.table('controlvsTB.lncRNA.transcript_level.deg.csv',sep = ',',header = T)
data$log2FoldChange
ggplot(data,aes(log2FoldChange,-log10(pvalue))) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed',color = '#999999') +
  geom_vline(xintercept = c(-1,1),linetype = 'dashed',color = '#999999')+
  geom_point(aes(size = -log10(pvalue),color = -log10(pvalue))) +
  scale_color_gradientn(values = seq(0,1,0.2),colors = c('#39489f','#39bbec','#f9ed36',"#f38466",'#b81f25')) +
  scale_size_continuous(range = c(1,3)) +
  theme_bw() +
  theme(panel.grid = element_blank())

exprset <- read.table('genes.readcount.txt',sep = '\t',header = T,row.names = 1,check.names = F)
TB_brucella <- c(2,3,4,5,6,7,8)
exprset <- exprset[,TB_brucella]
pdata <- c(rep('brucella',4),rep('TB',3))
library(edgeR)
dge <- DGEList(counts = exprset,group = pdata)
dge <- calcNormFactors(dge)
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)
et <- exactTest(dge)
topTags(et)
order_tags <- topTags(et,n=100000)
order_tags <- data.frame(order_tags)
write.csv(order_tags,file="TBvsbrucella.mRNA.gene_level.deg.csv",row.names=T)

library(ggplot2)
data <- read.table('TBvsbrucella.mRNA.gene_level.deg.csv',sep = ',',header = T)
data$log2FoldChange
ggplot(data,aes(log2FoldChange,-log10(pvalue))) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed',color = '#999999') +
  geom_vline(xintercept = c(-1,1),linetype = 'dashed',color = '#999999')+
  geom_point(aes(size = -log10(pvalue),color = -log10(pvalue))) +
  scale_color_gradientn(values = seq(0,1,0.2),colors = c('#39489f','#39bbec','#f9ed36',"#f38466",'#b81f25')) +
  scale_size_continuous(range = c(1,3)) +
  theme_bw() +
  theme(panel.grid = element_blank())

exprset <- read.table('lnc_transcripts.readcount.txt',sep = '\t',header = T,row.names = 1,check.names = F)
TB_brucella <- c(2,3,4,5,6,7,8)
exprset <- exprset[,TB_brucella]
pdata <- c(rep('brucella',4),rep('TB',3))
library(edgeR)
dge <- DGEList(counts = exprset,group = pdata)
dge <- calcNormFactors(dge)
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)
et <- exactTest(dge)
topTags(et)
order_tags <- topTags(et,n=100000)
order_tags <- data.frame(order_tags)
write.csv(order_tags,file="TBvsbrucella.lncRNA.transcript_level.deg.csv",row.names=T)

library(ggplot2)
data <- read.table('TBvsbrucella.lncRNA.transcript_level.deg.csv',sep = ',',header = T)
data$log2FoldChange
ggplot(data,aes(log2FoldChange,-log10(pvalue))) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed',color = '#999999') +
  geom_vline(xintercept = c(-1,1),linetype = 'dashed',color = '#999999')+
  geom_point(aes(size = -log10(pvalue),color = -log10(pvalue))) +
  scale_color_gradientn(values = seq(0,1,0.2),colors = c('#39489f','#39bbec','#f9ed36',"#f38466",'#b81f25')) +
  scale_size_continuous(range = c(1,3)) +
  theme_bw() +
  theme(panel.grid = element_blank())


library(clusterProfiler)
library(org.Hs.eg.db)
Enrichment <- read.table("jVenn.csv",sep = ",",stringsAsFactors = F,header = T)
gene <- as.character(Enrichment$CvsB.CvsT.TvsB)
ego_CC <- enrichGO(gene = gene$ENTREZID,OrgDb=org.Hs.eg.db,ont = "CC",keyType = "ENTREZID",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 1,qvalueCutoff = 1,readable = TRUE)
CC = as.data.frame(ego_CC @ result)
ego_BP <- enrichGO(gene = gene$ENTREZID,OrgDb=org.Hs.eg.db,ont = "BP",keyType = "ENTREZID",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 1,qvalueCutoff = 1,readable = TRUE)
BP = as.data.frame(ego_BP @ result)
ego_MF <- enrichGO(gene = gene$ENTREZID,OrgDb=org.Hs.eg.db,ont = "MF",keyType = "ENTREZID",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 1,qvalueCutoff = 1,readable = TRUE)
MF = as.data.frame(ego_MF @ result)
KEGG <- enrichKEGG(gene = gene$ENTREZID, organism = "hsa", keyType = "kegg",
                   pvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = 1,
                   qvalueCutoff = 1, use_internal_data = FALSE)
KEGG = setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
KEGG_Pathway <- as.data.frame(KEGG @ result)
Enrichment_all=rbind(CC,BP,MF,KEGG_Pathway)
write.table(Enrichment_all,"Enrichment_all.txt",row.names = F,quote = F,sep="\t")

library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
library(dplyr)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(org.Hs.eg.db)
library(patchwork)
library(WGCNA)
library(GSEABase)
library(randomcoloR)
traitData <- read.table('traitData.txt',sep = '\t',header = T,check.names = F)
datExpr0 <- read.table('select_WGCAN.txt',sep = '\t',header = T,row.names = 1,check.names = F)
for (df in colnames(traitData)) {
  traitData[,df]=traitData[,df]/max(traitData[,df])
  print(sd(traitData[,df]))
}
max(traitData)
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
allTraits = traitData
dim(allTraits)
names(allTraits)
# Form a data frame analogous to expression data that will hold the clinical traits.
fpkmSamples = rownames(datExpr0)
traitSamples =rownames(allTraits)
traitRows = match(fpkmSamples, traitSamples)
datTraits = allTraits[traitRows,]
rownames(datTraits) 
collectGarbage()
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
#sizeGrWindow(12,12)
pdf(file=paste0(afdir,"Sample dendrogram and trait heatmap.pdf"),width=12,height=11)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap", 
                    marAll = c(1, 12, 3, 1))
dev.off()
enableWGCNAThreads()
# Choose a set of soft-thresholding powers
powers = c(1:30)
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
# Plot the results:
#sizeGrWindow(9, 5)
pdf(file=paste0(afdir,"Scale independence.pdf"),width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


######chose the softPower
softPower =sft$powerEstimate
#显示软阈值
print(softPower)

adjacency = adjacency(datExpr0, power = softPower)

##### Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)

#sizeGrWindow(12,9)
pdf(file=paste0(afdir,"Gene clustering on TOM-based dissimilarity.pdf"),width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 50
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
pdf(file=paste0(afdir,"Dynamic Tree Cut.pdf"),width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
#sizeGrWindow(7, 6)
pdf(file=paste0(afdir,"Clustering of module eigengenes.pdf"),width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25######模块剪切高度
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

# Call an automatic merging function
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

#sizeGrWindow(12, 9)
pdf(file=paste0(afdir,"merged dynamic.pdf"), width = 9, height = 6.5)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
# Save module colors and labels for use in subsequent parts
#save(MEs, TOM, dissTOM,  moduleLabels, moduleColors, geneTree, sft, file = "networkConstruction-stepByStep.RData")

##############################relate modules to external clinical triats######################################
# Define numbers of genes and samples
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#sizeGrWindow(10,6)
pdf(file=paste0(afdir,"Module-trait relationships.pdf"),width=7,height=7.5)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(12, 8, 3, 3))

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


######## Define variable weight containing all column of datTraits

###MM and GS

modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

for (trait in traitNames){
  traitColumn=match(trait,traitNames)
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      #sizeGrWindow(7, 7)
      pdf(file=paste(afdir,"/9_", trait, "_", module,"_Module membership vs gene significance.pdf",sep=""),width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
    }
  }
}
#####
names(datExpr0)
probes = names(datExpr0)
#################export GS and MM############### 
geneInfo0 = data.frame(probes= probes,moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = paste0(afdir,"/10_GS_and_MM.xls"),sep="\t",row.names=F)