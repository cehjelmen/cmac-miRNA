countData <- read.csv("C:/Users/Carl Hjelmen/Desktop/cmacdiff.csv", header=T, row.names=1)
barplot(colSums(countData)*1e-6, names=colnames(countData), ylab="Library size (millions)")
library(DESeq2)
#period_only <- "Con_Early_Pf", "Con_Early_Pf", "Con_Early_Pf", "Con_Mid_Pup_2", "Con_Mid_Pup_2", "Con_Mid_Pup_2", "Early_PF", "Early_PF", "Early_PF", "Early_Pup", "Early_Pup", "Early_Pup", "Fast_Early_Pf", "Fast_Early_Pf", "Fast_Early_Pf", "Fast_Mid_Pup_2", "Fast_Mid_Pup_2", "Fast_Mid_Pup_2", "Feed_3rd", "Feed_3rd", "Feed_3rd", "Late_PF", "Late_PF", "Late_PF", "Late_Pup", "Late_Pup", "Late_Pup", "Mid_Pup_1", "Mid_Pup_1", "Mid_Pup_1", "Mid_Pup_2", "Mid_Pup_2", "Mid_Pup_2", "Slow_Early_Pf", "Slow_Early_Pf", "Slow_Early_Pf", "Slow_Mid_Pup_2", "Slow_Mid_Pup_2", "Slow_Mid_Pup_2"
#colData <- DataFrame(condition=factor(c(period_only)))
colData <- DataFrame(condition=factor(c("Con_Early_Pf", "Con_Early_Pf", "Con_Early_Pf", "Con_Mid_Pup_2", "Con_Mid_Pup_2", "Con_Mid_Pup_2", "Early_PF", "Early_PF", "Early_PF", "Early_Pup", "Early_Pup", "Early_Pup", "Fast_Early_Pf", "Fast_Early_Pf", "Fast_Early_Pf", "Fast_Mid_Pup_2", "Fast_Mid_Pup_2", "Fast_Mid_Pup_2", "Feed_3rd", "Feed_3rd", "Feed_3rd", "Late_PF", "Late_PF", "Late_PF", "Late_Pup", "Late_Pup", "Late_Pup", "Mid_Pup_1", "Mid_Pup_1", "Mid_Pup_1", "Mid_Pup_2", "Mid_Pup_2", "Mid_Pup_2", "Slow_Early_Pf", "Slow_Early_Pf", "Slow_Early_Pf", "Slow_Mid_Pup_2", "Slow_Mid_Pup_2", "Slow_Mid_Pup_2")))
dds <- DESeqDataSetFromMatrix(countData, colData, formula(~ condition))
dds <- DESeq(dds)
plotMA(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
head(resOrdered)
sig <- resOrdered[!is.na(resOrdered$padj) & resOrdered$padj<0.10 & abs(resOrdered$log2FoldChange)>=1,]
head(sig)
selected <- rownames(sig);selected
library("RColorBrewer")
library("gplots")
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100) ## hmcol <- heat.colors
heatmap.2( log2(counts(dds,normalized=TRUE)[rownames(dds) %in% selected,]+1), col = hmcol, scale="row", trace="none", margin=c(4,6), cexRow=0.5, cexCol=0.5, keysize=1 )
write.csv(counts(dds,normalized=TRUE), "C:/Users/Carl Hjelmen/Desktop/normal_data.csv")