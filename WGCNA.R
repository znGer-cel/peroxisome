#load package#
library(WGCNA)
library(flashClust)
library(iterators)
#load data#
library(readxl)
femData <- read_excel("VST.xlsx")
femData<-aggregate(.~GeneSymbol,mean,data=femData)
rownames(femData) <- femData$GeneSymbol
femData <- femData[ , -c(1)]
#datExpr
datExpr = as.data.frame(t(femData))
#check data
gsg = goodSamplesGenes(datExpr, verbose = 3) 
gsg$allOK 
if (!gsg$allOK)
{if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}
#Check outliers
sampleTree = hclust(dist(datExpr), method = "average")
sizeGrWindow(12,9) 
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 400, col = "red")
clust = cutreeStatic(sampleTree, cutHeight = 400, minSize = 10)
table(clust)   
keepSamples = (clust==1)  
datExpr = datExpr[keepSamples, ]  
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#load Traits#
allTraits = read.csv("Traits.csv", sep=',',header= 1);
dim(allTraits)  
names(allTraits)
femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$sample);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];
collectGarbage()
sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE) 
plotDendroAndColors(sampleTree2, traitColors,groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
#save#
save(datExpr, datTraits, file = "FemaleLiver-01-dataInput.RData")
###Creat WGCNA network###
library(WGCNA)
options(stringsAsFactors = FALSE);
lnames = load(file = "FemaleLiver-01-dataInput.RData");
lnames
###chose soft power###
powers = c(seq(1, 30, by = 1)) 
sft <- pickSoftThreshold(datExpr,powerVector = powers,networkType = "signed",
                         RsquaredCut = 0.80,verbose = 5)
sft$powerEstimate
par(mfrow = c(1,2));
cex1 = 0.85;
#Scale Free Topology Model Fit, signed R^2
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red", pos=4);
abline(h=0.80,col="red")
#Mean Connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.85,col="red")
##chose power##
net = blockwiseModules(datExpr, power = 18,
                       TOMType = "signed",networkType = "signed", minModuleSize = 100,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)
table(net$colors)
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "FemaleLiver-02-networkConstruction-auto.RData")
###heatmap###
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(datTraits),
               yLabels = names(MEs),ySymbols = names(MEs),
               colorLabels = FALSE,colors = greenWhiteRed(50),
               textMatrix = textMatrix,setStdMargins = FALSE,
               cex.text = 0.7,cex.lab = 1,zlim = c(-1,1),
               main = paste("Module-trait relationships"))
position = as.data.frame(datTraits$TIME);
names(position) = "time";
#result#
geneModule <- data.frame(GeneSymbol = colnames(datExpr),Module = moduleColors)
write.csv(geneModule, file = "WGCNA_gene_module.csv", row.names = FALSE)

## Repeated random subsampling stability analysis
set.seed(123)
## 
moduleColors <- labels2colors(net$colors)
names(moduleColors) <- colnames(datExpr)
nIter <- 10
sampleFraction <- 0.8
moduleLevels <- sort(unique(moduleColors))
bootstrapRetention <- matrix(NA,nrow = length(moduleLevels),ncol = nIter,
  dimnames = list(moduleLevels,paste0("Iter_",1:nIter)))
for(i in 1:nIter){
  cat("Bootstrap:",i,"\n")
  sample_idx <- sample(1:nrow(datExpr),floor(sampleFraction*nrow(datExpr)),replace = FALSE)
  dat_boot <- datExpr[sample_idx,]
  net_boot <- blockwiseModules(dat_boot,power = 18,TOMType = "signed",
                               networkType = "signed",minModuleSize = 100,
                               reassignThreshold = 0,mergeCutHeight = 0.25,
                               numericLabels = TRUE,pamRespectsDendro = FALSE,verbose = 0)
  bootColors <- labels2colors(net_boot$colors)
  names(bootColors) <- colnames(dat_boot)
  for(m in moduleLevels){
    genes1 <- names(moduleColors)[moduleColors==m]
    score <- numeric()
    for(b in unique(bootColors)){
      genes2 <- names(bootColors)[bootColors==b]
      jaccard <-
        length(intersect(genes1,genes2))/
        length(union(genes1,genes2))
      score <- c(score,jaccard)
    }
    bootstrapRetention[m,i] <- max(score)
  }
}
bootstrap_summary <- data.frame(Module = rownames(bootstrapRetention),
                                Mean_Jaccard = rowMeans(bootstrapRetention),
                                SD = apply(bootstrapRetention,1,sd),
                                Lower95 = apply(bootstrapRetention,1,quantile,0.025),
                                Upper95 = apply(bootstrapRetention,1,quantile,0.975))
print(bootstrap_summary)
write.csv(bootstrap_summary,"Bootstrap_Module_Stability.csv",row.names = FALSE)
## Leave-one-out stability
leaveRetention <- matrix(NA,nrow = length(moduleLevels),ncol = nrow(datExpr),
  dimnames = list(moduleLevels,rownames(datExpr)))
for(i in 1:nrow(datExpr)){cat("Leave sample:",rownames(datExpr)[i],"\n")
  dat_leave <- datExpr[-i,]
  net_leave <- blockwiseModules(dat_leave,power = 18,TOMType = "signed",
                                networkType = "signed",minModuleSize = 100,
                                reassignThreshold = 0,mergeCutHeight = 0.25,
                                numericLabels = TRUE,pamRespectsDendro = FALSE,verbose = 0)
  leaveColors <- labels2colors(net_leave$colors)
  names(leaveColors) <- colnames(dat_leave)
  for(m in moduleLevels){
    genes1 <- names(moduleColors)[moduleColors==m]
    score <- numeric()
    for(b in unique(leaveColors)){
      genes2 <- names(leaveColors)[leaveColors==b]
      jaccard <-
        length(intersect(genes1,genes2))/
        length(union(genes1,genes2))
      score <- c(score,jaccard)}
    leaveRetention[m,i] <- max(score)}}
leave_summary <- data.frame(
  Module = rownames(leaveRetention),
  Mean_Jaccard = rowMeans(leaveRetention),
  SD = apply(leaveRetention,1,sd),
  Lower95 = apply(leaveRetention,1,quantile,0.025),
  Upper95 = apply(leaveRetention,1,quantile,0.975))
print(leave_summary)
write.csv(leave_summary,"LeaveOneOut_Module_Stability.csv",row.names = FALSE)
## Merge stability results
stability_summary <- merge(bootstrap_summary,leave_summary,by = "Module",
  suffixes = c("_Bootstrap","_LeaveOneOut"))
write.csv(stability_summary,"WGCNA_Module_Stability_Summary.csv",row.names = FALSE)
print(stability_summary)