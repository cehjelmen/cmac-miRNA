### Working on compiling all of the Rscript data for the C. macellaria miRNA paper
###Carl Hjelmen
###11-20-21

#setwd to R_files folder


###Running DESeq2
####DESeq2 on wtpuparial data####
#read in data
countData <- read.csv("DESeq2/wt_pupae/wt_pupae_counts.csv", header=T, row.names=1)
barplot(colSums(countData)*1e-6, names=colnames(countData), ylab="Library size (millions)")
library(DESeq2)
colData <- DataFrame(condition=factor(c("Early_Pup", "Early_Pup", "Early_Pup", "Mid_Pup_1", "Mid_Pup_1", "Mid_Pup_1", "Mid_Pup_2", "Mid_Pup_2", "Mid_Pup_2","Late_Pup", "Late_Pup", "Late_Pup")))
dds <- DESeqDataSetFromMatrix(countData, colData, formula(~ condition))
dds <- DESeq(dds)
plotMA(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
#ouput resOrdered info
write.csv(resOrdered, "DESeq2/wt_pupae/puparial_wt_resOrdered.csv")
head(resOrdered)
sig <- resOrdered[!is.na(resOrdered$padj) & resOrdered$padj<0.10 & abs(resOrdered$log2FoldChange)>=1,]
head(sig)
#ouput sig info
write.csv(sig, "DESeq2/wt_pupae/cmacpuparial_wt_sig.csv")
selected <- rownames(sig);selected
library("RColorBrewer")
library("gplots")
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100) ## hmcol <- heat.colors
heatmap.2( log2(counts(dds,normalized=TRUE)[rownames(dds) %in% selected,]+1), col = hmcol, scale="row", trace="none", margin=c(8,8), cexRow=0.7, cexCol=0.75, keysize=1 )
#output the normalized data for future analysis
write.csv(counts(dds,normalized=TRUE), "DESeq2/wt_pupae/pupaecounts_normal_data.csv")

####DESeq2 on wt larval data####
countData <- read.csv("DESeq2/wt_larval/wt_larval_counts.csv", header=T, row.names=1)
barplot(colSums(countData)*1e-6, names=colnames(countData), ylab="Library size (millions)")
library(DESeq2)
colData <- DataFrame(condition=factor(c("Feed_3rd", "Feed_3rd", "Feed_3rd", "Early_PF", "Early_PF", "Early_PF", "Late_PF", "Late_PF", "Late_PF")))
dds <- DESeqDataSetFromMatrix(countData, colData, formula(~ condition))
dds <- DESeq(dds)
plotMA(dds)
res <- results(dds)
#write output table for res
write.table(res, "DESeq2/wt_larval/wt_larval_res.txt")
#order the results by significance
resOrdered <- res[order(res$padj),]
head(resOrdered)
#ouput resOrdered info
write.csv(resOrdered, "DESeq2/wt_larval/wt_larval_resOrdered.csv")
#what is the data that is just significant
sig <- resOrdered[!is.na(resOrdered$padj) & resOrdered$padj<0.10 & abs(resOrdered$log2FoldChange)>=1,]
head(sig)
#ouput sig info
write.csv(sig, "DESeq2/wt_larval/wt_larval_sig.csv")
selected <- rownames(sig);selected
library("RColorBrewer")
library("gplots")
#make heatmap
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100) ## hmcol <- heat.colors
heatmap.2( log2(counts(dds,normalized=TRUE)[rownames(dds) %in% selected,]+1), col = hmcol, scale="row", trace="none", margin=c(4,6), cexRow=0.8, cexCol=0.8, keysize=1 )
#output the normalized data for future analysis
write.csv(counts(dds,normalized=TRUE), "DESeq2/wt_larval/wt_larval_counts_normal_data.csv")

####DESeq2 Selection Early Postfeeding####
countData <- read.csv("DESeq2/sel_epf/sel_epf_counts.csv", header=T, row.names=1)
barplot(colSums(countData)*1e-6, names=colnames(countData), ylab="Library size (millions)")
library(DESeq2)
colData <- DataFrame(condition=factor(c("Con_Early_Pf", "Con_Early_Pf", "Con_Early_Pf", "Fast_Early_Pf", "Fast_Early_Pf", "Fast_Early_Pf", "Slow_Early_Pf", "Slow_Early_Pf", "Slow_Early_Pf")))
dds <- DESeqDataSetFromMatrix(countData, colData, formula(~ condition))
dds <- DESeq(dds)
plotMA(dds)
res <- results(dds)
#output results
write.table(res, "DESeq2/sel_epf/sel_epf_res.txt")

resOrdered <- res[order(res$padj),]
head(resOrdered)

#ouput resOrdered info
write.csv(resOrdered, "DESeq2/sel_epf/sel_epf_resOrdered.csv")

#what is the data that is just significant
sig <- resOrdered[!is.na(resOrdered$padj) & resOrdered$padj<0.10 & abs(resOrdered$log2FoldChange)>=1,]
head(sig)
#output sig info
write.csv(sig, "DESeq2/sel_epf/sel_epf_sig.csv")
selected <- rownames(sig);selected
library("RColorBrewer")
library("gplots")
#make heatmap
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100) ## hmcol <- heat.colors
heatmap.2( log2(counts(dds,normalized=TRUE)[rownames(dds) %in% selected,]+1), col = hmcol, scale="row", trace="none", margin=c(6,8), cexRow=0.8, cexCol=0.8, keysize=1 )
#output the normalized data for future analysis
write.csv(counts(dds,normalized=TRUE), "DESeq2/sel_epf/sel_epf_counts_normal_data.csv")

####DESeq2 Selection MP2####
countData <- read.csv("DESeq2/sel_mp2/sel_mp2_counts.csv", header=T, row.names=1)
barplot(colSums(countData)*1e-6, names=colnames(countData), ylab="Library size (millions)")
library(DESeq2)
colData <- DataFrame(condition=factor(c( "Con_Mid_Pup_2", "Con_Mid_Pup_2", "Con_Mid_Pup_2", "Fast_Mid_Pup_2", "Fast_Mid_Pup_2", "Fast_Mid_Pup_2", "Slow_Mid_Pup_2", "Slow_Mid_Pup_2", "Slow_Mid_Pup_2")))
dds <- DESeqDataSetFromMatrix(countData, colData, formula(~ condition))
dds <- DESeq(dds)
plotMA(dds)
#overall results
res <- results(dds)
#outputresults
write.table(res, "DESeq2/sel_mp2/sel_mp2_res.txt")
resOrdered <- res[order(res$padj),]
head(resOrdered)
#output resOrdered
write.csv(resOrdered, "DESeq2/sel_mp2/sel_mp2_resOrdered.csv")

#what is the data that is just significant
sig <- resOrdered[!is.na(resOrdered$padj) & resOrdered$padj<0.10 & abs(resOrdered$log2FoldChange)>=1,]
head(sig)
#output sig
write.csv(sig, "DESeq2/sel_mp2/sel_mp2_sig.csv")
selected <- rownames(sig);selected
library("RColorBrewer")
library("gplots")
#make heatmap
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100) ## hmcol <- heat.colors
heatmap.2( log2(counts(dds,normalized=TRUE)[rownames(dds) %in% selected,]+1), col = hmcol, scale="row", trace="none", margin=c(6,6), cexRow=0.8, cexCol=0.8, keysize=1 )
#output normalized data
write.csv(counts(dds,normalized=TRUE), "DESeq2/sel_mp2/sel_mp2_counts_normal_data.csv")

####DESeq2 Sex Specific Data####
countData <- read.csv("DESeq2/sex/sex_counts.csv", header=T, row.names=1)
barplot(colSums(countData)*1e-6, names=colnames(countData), ylab="Library size (millions)")
library(DESeq2)
colData <- DataFrame(condition=factor(c("Female_25C","Female_25C","Female_25C","Female_25C","Female_25C","Female_25C","Female_25C","Female_25C","Female_25C","Male_25C","Male_25C","Male_25C","Male_25C","Male_25C","Male_25C","Male_25C","Male_25C")))
dds <- DESeqDataSetFromMatrix(countData, colData, formula(~ condition))
dds <- DESeq(dds)
plotMA(dds)
res <- results(dds)
#outputresults
write.table(res, "DESeq2/sex/early_intrapuparial_sex_res.txt")
resOrdered <- res[order(res$padj),]
head(resOrdered)
#output resOrdered
write.csv(resOrdered, "DESeq2/sex/early_intrapuparial_sex_resOrdered.csv")

sig <- resOrdered[!is.na(resOrdered$padj) & resOrdered$padj<0.10 & abs(resOrdered$log2FoldChange)>=1,]
head(sig)
#output sig
write.csv(sig, "DESeq2/sex/early_intrapuparial_sex_sig.csv")
selected <- rownames(sig);selected
library("RColorBrewer")
library("gplots")
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100) ## hmcol <- heat.colors
heatmap.2( log2(counts(dds,normalized=TRUE)[rownames(dds) %in% selected,]+1), dendrogram="none", col = hmcol, scale="row", trace="none", margin=c(4,6), cexRow=0.5, cexCol=0.5, keysize=1 )
#output normalized data
write.csv(counts(dds,normalized=TRUE), "DESeq2/sex/sex_counts_normal_data.csv")









####







###PCA
#significantly differentially expressed PCA data

####Larval samples####
#datafile made from getting normalized expression information for
#the 10 miRNA identified as signifnicantly differentially expressed in
#wild-type larval samples
#file includes normalized expression information for development-time-selected
#samples as well, as well as selection regime

library(ggplot2)
library(ggbiplot)
larval_pca<-read.csv("pca/larvae_pca_sig.csv", as.is=TRUE)
#structure of data
#str(larval_pca)
#change to factors
larval_pca$Stage<-as.factor(larval_pca$Stage)
larval_pca$Grouping<-as.factor(larval_pca$Grouping)
larval_pca$Type<-as.factor(larval_pca$Type)

#confirm that levels are right for factors
levels(larval_pca$Stage)


larval_pca.active<-larval_pca[1:9,4:13]
larval_pca.active.stage<-(larval_pca[1:9,3])
larval_pca.active.groups<-larval_pca[1:9, 2]

#str(larval_pca.active.stage)
#ordered the stages chronologically
larval_pca.active.stage<-ordered(larval_pca.active.stage, levels=c("Feeding 3rd", "Early Postfeeding", "Late Postfeeding"))

#ordered stages chronologically with selection info
larval_pca.stage<-(larval_pca[,2])
larval_pca.stage<-ordered(larval_pca.stage, levels=c("Feeding 3rd", "Early Postfeeding", "Late Postfeeding", "Control EPF", "Fast EPF", "Slow EPF"))

#input the information on selection regime
larval_pca.selinfo<-(larval_pca[,1])
#compute pca
#run PCA
larval_pca.pca<-prcomp(larval_pca.active, center=TRUE, scale.=TRUE)
#summary of results
summary(larval_pca.pca)

library(ggbiplot)
library(viridis)
#making plot


#projection of selection samples onto PCA
larval_pca.selection<-larval_pca[10:18, 4:13]
larval_pca.selection.stage<-(larval_pca[10:18, 2])
  
ls.sc<-scale(larval_pca.selection, center=larval_pca.pca$center)
ls.pred<-ls.sc %*% larval_pca.pca$rotation
  
larval_pca.plusproj.pca <- larval_pca.pca
larval_pca.plusproj.pca$x <- rbind(larval_pca.plusproj.pca$x, ls.pred)

##making larvae plot
larv<-ggbiplot(larval_pca.plusproj.pca, ellipse=TRUE, obs.scale = 1, var.scale=1, groups=larval_pca.stage, var.axes = FALSE)+
  scale_color_viridis_d(name="Values", end=0.9)+geom_point(aes(shape=larval_pca.selinfo, fill=larval_pca.stage, size=3))+
  scale_shape_manual(name= "Values", values=c(24,22,25,21))+  scale_fill_viridis_d(end=0.9)+
  ggtitle("Larval development")+theme_minimal()+ theme(legend.position="right")+
  guides(fill="none", size="none", shape=guide_legend(override.aes=list(size=5)), color=guide_legend(override.aes = list(size=5)))+
  theme(legend.title = element_blank())+
  theme(axis.title = element_text(size=14, face="bold"), axis.text = element_text(face="bold", size=10), title=element_text(face="bold", size=16))+
  theme(legend.text = element_text(size=12))+xlim(c(-5, 7.5))+ylim(c(-5,3))
larv



####Pupae samples####
#datafile made from getting normalized expression information for
#the 17 miRNA identified as signifnicantly differentially expressed in
#wild-type intrapuparial samples
#file includes normalized expression information for development-time-selected
#samples as well, as well as selection regime
library(ggplot2)
library(ggbiplot)
pupae_pca<-read.csv("pca/pupae_pca_sig.csv", as.is=TRUE)
#structure of data
#str(larval_pca)
#change to factors
pupae_pca$Stage<-as.factor(pupae_pca$Stage)
pupae_pca$Grouping<-as.factor(pupae_pca$Grouping)
pupae_pca$Type<-as.factor(pupae_pca$Type)

#confirm that levels are right for factors
levels(pupae_pca$Stage)


pupae_pca.active<-pupae_pca[1:12, 4:20]
pupae_pca.active.stage<-(pupae_pca[1:12,3])
pupae_pca.active.groups<-pupae_pca[1:12, 2]

#str(pupae_pca.active.stage)
#ordered the stages chronologically
pupae_pca.active.stage<-ordered(pupae_pca.active.stage, levels=c("Early Puparial", "Mid Puparial 1", "Mid Puparial 2", "Late Puparial"))

#ordered stages chronologically with selection info
pupae_pca.stage<-(pupae_pca[,2])
pupae_pca.stage<-ordered(pupae_pca.stage, levels=c("Early", "Mid 1", "Mid 2", "Late", "Control Mid", "Fast Mid", "Slow Mid"))

#input the information on selection regime
pupae_pca.selinfo<-(pupae_pca[,1])
#compute pca
#run PCA
pupae_pca.pca<-prcomp(pupae_pca.active, center=TRUE, scale.=TRUE)
#summary of results
summary(pupae_pca.pca)

library(ggbiplot)
library(viridis)
#making plot


#projection of selection samples onto PCA
pupae_pca.selection<-pupae_pca[13:21, 4:20]
pupae_pca.selection.stage<-(pupae_pca[13:21, 2])

s.sc<-scale(pupae_pca.selection, center=pupae_pca.pca$center)
s.pred<-s.sc %*% pupae_pca.pca$rotation

pupae_pca.plusproj.pca <- pupae_pca.pca
pupae_pca.plusproj.pca$x <- rbind(pupae_pca.plusproj.pca$x, s.pred)

##making pupae plot
pupae<-ggbiplot(pupae_pca.plusproj.pca, ellipse=TRUE, obs.scale = 1, var.scale=1, groups=pupae_pca.stage, var.axes = FALSE)+
  scale_color_viridis_d(name="Values", end=0.9)+geom_point(aes(shape=pupae_pca.selinfo, fill=pupae_pca.stage, size=3))+
  scale_shape_manual(name= "Values", values=c(24,22,25,21))+  scale_fill_viridis_d(end=0.9)+
  ggtitle("Intrapuparial development")+theme_minimal()+ theme(legend.position="right")+
  guides(fill="none", size="none", shape=guide_legend(override.aes=list(size=5)), color=guide_legend(override.aes = list(size=5)))+
  theme(legend.title = element_blank())+
  theme(axis.title = element_text(size=14, face="bold"), axis.text = element_text(face="bold", size=10), title=element_text(face="bold", size=16))+
  theme(legend.text = element_text(size=12))+xlim(c(-5, 7.5))+ylim(c(-5,3))
pupae


####making a combination plot for PCA of larvae and puparial####
#you should run both the larval and pupae samples pca scripts before this one
library(ggpubr)
ggarrange(larv, pupae, labels=c("A", "B"), ncol=1, nrow=2, font.label=list(size=25))


####Making pca will all data####
###first plot is including vertebrate data
#read in data
#this is all normalized expression data for all stages
all_data_pca<-read.csv("pca/all_normalized_stages_pca.csv", header=T, as.is=TRUE)
#read in active form of data, the normal developmen data
cmac.active<-all_data_pca[1:21, 4:222]
cmac.active.stage<-(all_data_pca[1:21,3])
cmac.active.stage<-factor(cmac.active.stage, levels=c("Feeding 3rd", "Early Postfeeding", "Late Postfeeding", "Early Puparial", "Mid Puparial 1", "Mid Puparial 2", "Late Puparial"))

cmac.stage<-(all_data_pca[,3])

cmac.stage

#compute pca
res.pca<-prcomp(cmac.active, center=TRUE, scale.=TRUE)
summary(res.pca)

library(ggbiplot)
library(viridis)
all_pca<-ggbiplot(res.pca, ellipse=TRUE, obs.scale = 1, var.scale=1, groups=cmac.active.stage, var.axes = FALSE)+
  scale_fill_viridis_d(end=1)+scale_color_viridis_d(begin=0,end=1)+
  geom_point(shape=21,aes(fill=cmac.active.stage),colour="black", size=4)+
  theme(axis.title = element_text(face="bold", size=14), axis.text=element_text(face="bold", size=10, colour="black"))+
  theme(legend.title=element_blank(), legend.text =element_text(face="bold", size=9), legend.position = "bottom" )+
  guides(color="none")+xlim(-17, 15)+ylim(-10,15)
all_pca

####pca with no vertebrate mirna####

#read in data
#file was made by removing all mirna which were identified as vertebrate

no_vert_data<-read.csv("pca/all_normalized_stages_pca_no_vert.csv", as.is=TRUE)

#read in active form of data, the normal developmen data

nv.cmac.active<-no_vert_data[1:21, 3:163]
nv.cmac.active.stage<-(no_vert_data[1:21,2])
str(nv.cmac.active.stage)
nv.cmac.active.stage<-factor(nv.cmac.active.stage, levels=c("Feeding 3rd", "Early Postfeeding", "Late Postfeeding", "Early Puparial", "Mid Puparial 1", "Mid Puparial 2", "Late Puparial"))

nv.cmac.stage<-(no_vert_data[,2])

nv.cmac.stage
#read in other selection data

#nv.cmac.selection<-no_vert_data[22:39, 3:163]
#nv.cmac.selection.stage<-(no_vert_data[22:39,2])

#compute pca

nv.res.pca<-prcomp(nv.cmac.active, center=TRUE, scale.=TRUE)

summary(nv.res.pca)

library(ggbiplot)

library(viridis)
novert<-ggbiplot(nv.res.pca, ellipse=TRUE, obs.scale = 1, var.scale=1, groups=nv.cmac.active.stage, var.axes = FALSE)+
  scale_fill_viridis_d(end=1)+scale_color_viridis_d(begin=0,end=1)+
  geom_point(shape=21,aes(fill=nv.cmac.active.stage),colour="black", size=4)+
  theme(axis.title = element_text(face="bold", size=14), axis.text=element_text(face="bold", size=10, colour="black"))+
  theme(legend.title=element_blank(), legend.text =element_text(face="bold", size=9), legend.position = "bottom" )+
  guides(color="none")+xlim(-17, 15)+ylim(-10,15)
novert



####making a combination plot for PCA of all stages with and without vert data####
#Be sure to run the pca script from both groups in the above
library(ggpubr)
ggarrange(all_pca, novert, labels=c("A", "B"), ncol=2, nrow=1, font.label=list(size=25), common.legend = TRUE, legend="bottom")




###Random forest

####using all data together####
#training models with no vertebrate miRNA included
library(randomForest)

percent.dev<-read.csv("rf/all_dat/cmacnovertwtrf.csv")
str(percent.dev)
stagedata<-read.csv("rf/all_dat/cmacnovertwt.csv")
seldata<-read.csv("rf/all_dat/cmacselectionrf.csv")

cmacstage.rf<-randomForest(Stage~., data=stagedata, ntree=5000)
cmacstage.rf


cmacper.rf<-randomForest(DevPer ~., data=percent.dev, ntree=5000)
cmacper.rf

#see how many trees you should use, look at plot to see where error drops
plot(cmacstage.rf)
plot(cmacper.rf)

print(cmacstage.rf)
print(cmacper.rf)

varImpPlot(cmacstage.rf, sort=TRUE, n.var=50, type=NULL, class=NULL, scale=TRUE, main="Random Forest of Cmac miRNA Stage"  )
varImpPlot(cmacper.rf, sort=TRUE, n.var=50, type=NULL, class=NULL, scale=TRUE, main="Random Forest of Cmac miRNA Percentage no vertebrate"  )


importance(cmacstage.rf)
importance(cmacper.rf)

write.csv(importance(cmacstage.rf), "rf/all_dat/rfmeandecreaseginistage.csv")
write.csv(importance(cmacper.rf), "rf/all_dat/rfmeandecreaseginiper.csv")

dev.per.output<-predict(cmacper.rf, seldata, type="response", norm.votes = TRUE)
stage.output<-predict(cmacstage.rf, seldata, type="response", norm.votes = TRUE)

write.csv(dev.peroutput, "rf/all_dat/randomforestselectionpredictiondevpercent.csv")
write.csv(stage.output, "rf/all_dat/randomforestselectionpredictionstage.csv")



#with sig DE data
####larval samples####
#prediction proportion of development with data from only miRNA found to be significantly differentially
#expressed by DESeq2 in larvae
library(randomForest)
larvae_rf<-read.csv("rf/larval/larvae_rf.csv")
larvae_sel_rf<-read.csv("rf/larval/larvae_sel_rf.csv")

cmacper.rf<-randomForest(Proportion ~., data=larvae_rf, ntree=5000)
cmacper.rf

#see how many trees you should use, look at plot to see where error drops
plot(cmacper.rf)
print(cmacper.rf)

varImpPlot(cmacper.rf, sort=TRUE, n.var=10, type=NULL, class=NULL, scale=TRUE, main="Random Forest of Cmac miRNA Percentage Puparial"  )
importance(cmacper.rf)

write.csv(importance(cmacper.rf), "rf/larval/rfmeandecreaseginiperlarvae.csv")

output<-predict(cmacper.rf, larvae_sel_rf, type="response", norm.votes = TRUE)
?predict.randomForest
write.csv(output, "rf/larval/randomforestselectionpredictiondevpercentlarvae.csv")


####intrapuparial samples####
#files are  made from normalized expression information for the miRNA found to be
#significanlty differentially expressed in intrapuparial samples by DESeq2
library(randomForest)

pupae_rf<-read.csv("rf/pupae/pupae_rf.csv")
str(pupae_rf)
pupae_sel_rf<-read.csv("rf/pupae/pupae_sel_rf.csv")

cmac_pup_per.rf<-randomForest(Proportion ~., data=pupae_rf, ntree=5000)
cmac_pup_per.rf

#see how many trees you should use, look at plot to see where error drops
plot(cmac_pup_per.rf)
print(cmac_pup_per.rf)

varImpPlot(cmac_pup_per.rf, sort=TRUE, n.var=17, type=NULL, class=NULL, scale=TRUE, main="Random Forest of Cmac miRNA Percentage Puparial"  )

importance(cmac_pup_per.rf)
write.csv(importance(cmac_pup_per.rf), "rf/pupae/rfmeandecreaseginiper.csv")


output<-predict(cmac_pup_per.rf, pupae_sel_rf, type="response", norm.votes = TRUE)
write.csv(output, "rf/pupae/randomforestselectionpredictiondevpercent.csv")





####data specific to validation####
pup<-read.csv("rf/qPCR_dat/pcr_val_pup.csv")
str(pup)
selpup<-read.csv("rf/qPCR_dat/sel_pcr_val_pup.csv")

pcr.cmac.per.rf<-randomForest(Proportion ~., data=pup, ntree=5000)
print(pcr.cmac.per.rf)
plot(pcr.cmac.per.rf)

varImpPlot(pcr.cmac.per.rf, sort=TRUE, type=NULL, class=NULL, scale=TRUE,  n.var=17, main="Random Forest of Cmac miRNA Percentage")
importance(pcr.cmac.per.rf)
write.csv(importance(pcr.cmac.per.rf), "rf/qPCR_dat/rfmeandecreaseginiper.csv")

output<-predict(pcr.cmac.per.rf, selpup, type="response", norm.votes = TRUE)
output
write.csv(output, "rf/qPCR_dat/randomforestselectionpredictiondevpercent.csv")





###Making boxplots for qPCR and RNA Seq data

####Scripts for boxplots####
qpcr<-read.csv("boxplots/qpcr.csv", as.is=T)
rnaseq<-read.csv("boxplots/rnaseq.csv", as.is=T)

str(qpcr)
qpcr$Stage<-factor(qpcr$Stage, levels=c("Early Puparial", "Mid Puparial 1", "Mid Puparial 2", "Late Puparial"))

str(rnaseq)
rnaseq$Stage<-factor(rnaseq$Stage, levels=c("Feeding 3rd Instar", "Early Postfeeding", "Late Postfeeding", "Early Puparial", "Mid Puparial 1", "Mid Puparial 2", "Late Puparial"))



library(ggplot2)
####qpcr data boxplots####
m957q<-ggplot(qpcr, aes(Stage, miR957))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("miR957 qPCR")+ylab("Fold Change Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
m957q
mlet7q<-ggplot(qpcr, aes(Stage, let7))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("let7 qPCR")+ylab("Fold Change Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mlet7q
m184q<-ggplot(qpcr, aes(Stage, miR184))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("miR184 qPCR")+ylab("Fold Change Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12), axis.title = element_text( size=16), title = element_text(size=18))
m184q
m277q<-ggplot(qpcr, aes(Stage, miR277))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("miR277 qPCR")+ylab("Fold Change Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
m277q
m92bq<-ggplot(qpcr, aes(Stage, miR92b))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("miR92b qPCR")+ylab("Fold Change Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
m92bq  
m8q<-ggplot(qpcr, aes(Stage, miR8))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("miR8 qPCR")+ylab("Fold Change Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
m8q  
m31aq<-ggplot(qpcr, aes(Stage, miR31a))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("miR31a qPCR")+ylab("Fold Change Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
m31aq
m276bq<-ggplot(qpcr, aes(Stage, miR276b))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("miR276b qPCR")+ylab("Fold Change Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
m276bq  
mbantamq<-ggplot(qpcr, aes(Stage, bantam))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("bantam qPCR")+ylab("Fold Change Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mbantamq
m277oq<-ggplot(qpcr, aes(Stage, miR277o))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("miR277 qPCR")+ylab("Fold Change Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
m277oq
m317q<-ggplot(qpcr, aes(Stage, miR317))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("miR317 qPCR")+ylab("Fold Change Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
m317q
m305q<-ggplot(qpcr, aes(Stage, miR305))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("miR305 qPCR")+ylab("Fold Change Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
m305q
m10aq<-ggplot(qpcr, aes(Stage, miR10a))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("miR10a qPCR")+ylab("Fold Change Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
m10aq
m10q<-ggplot(qpcr, aes(Stage, miR10))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("miR10 qPCR")+ylab("Fold Change Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
m10q
m1q<-ggplot(qpcr, aes(Stage, miR1))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("miR1 qPCR")+ylab("Fold Change Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
m1q
m22q<-ggplot(qpcr, aes(Stage, miR22))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("miR22 qPCR")+ylab("Fold Change Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
m22q


####rnaseq data boxplots####

mir8<-ggplot(rnaseq, aes(Stage, Dana.miR.8.RA))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("Dana-miR-8-RA RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mir8

mir10<-ggplot(rnaseq, aes(Stage, Dmel.mir.10.RA))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("Dmel-miR-10-RA RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mir10

mir276a<-ggplot(rnaseq, aes(Stage, Dana.miR.276a.RB))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("Dana-miR-276a-RB RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mir276a

mir184<-ggplot(rnaseq, aes(Stage, Dana.miR.184.RB))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("Dana-miR-184-RB RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mir184

mir31a<-ggplot(rnaseq, aes(Stage, Dana.miR.31a.RA))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("Dana-miR-31a-RA RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mir31a

miranabantam<-ggplot(rnaseq, aes(Stage, Dana.bantam.RA))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("Dana-bantam-RA RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
miranabantam

mirdana10<-ggplot(rnaseq, aes(Stage, Dana.miR.10.RA))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("Dana-miR-10-RA RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirdana10

mirdana92b<-ggplot(rnaseq, aes(Stage, Dana.miR.92b.RA))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("Dana-miR-92b-RA RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirdana92b

mirmel10<-ggplot(rnaseq, aes(Stage, Dmel.mir.10.RB))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("Dmel-mir-10-RB RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirmel10

mir317<-ggplot(rnaseq, aes(Stage, ame.miR.317))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("ame-miR-317 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mir317

mirhme8<-ggplot(rnaseq, aes(Stage, hme.miR.8))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("hme-miR-8 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirhme8

mirdvir8.5<-ggplot(rnaseq, aes(Stage, dvi.miR.8.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("dvi-miR-8-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirdvir8.5

mirdvir305<-ggplot(rnaseq, aes(Stage, Dvir.miR.305.RA))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("Dvir-miR-305-RA RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirdvir305

miraae305<-ggplot(rnaseq, aes(Stage, aae.miR.305.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("aae-miR-305-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
miraae305

mirdana276b<-ggplot(rnaseq, aes(Stage, Dana.miR.276b.RB))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("Dana-miR-276b-RB RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirdana276b

mirdana277<-ggplot(rnaseq, aes(Stage, Dana.miR.277.RA))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("Dana-miR-277-RA RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirdana277

miraaebantam<-ggplot(rnaseq, aes(Stage, aae.bantam.3p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("aae-bantam-3p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
miraaebantam

mirxtrmir22<-ggplot(rnaseq, aes(Stage, xtr.miR.22.3p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("xtr-miR-22-3p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirxtrmir22

mirhme277<-ggplot(rnaseq, aes(Stage, hme.miR.277))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("hme-miR-277 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirhme277

mirbfl184<-ggplot(rnaseq, aes(Stage, bfl.miR.184.3p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("bfl-miR-184-3p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirbfl184

mirhme305<-ggplot(rnaseq, aes(Stage, hme.miR.305))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("hme-miR-305 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirhme305

mirdanalet7<-ggplot(rnaseq, aes(Stage, Dana.let.7.RA))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("Dana-let7-RA RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirdanalet7

miraae276<-ggplot(rnaseq, aes(Stage, aae.miR.276.3p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("aae-miR-276-3p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
miraae276

mirhmelet7<-ggplot(rnaseq, aes(Stage, hme.let.7))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("hme-let7 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirhmelet7

mircqubantam<-ggplot(rnaseq, aes(Stage, cqu.bantam.3p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("cqu-bantam-3p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mircqubantam

mirdana305<-ggplot(rnaseq, aes(Stage, Dana.miR.305.RA))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("Dana-miR-305-RA RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirdana305

mirisc10<-ggplot(rnaseq, aes(Stage, isc.miR.10_))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("isc-miR-10 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirisc10

mirdana1<-ggplot(rnaseq, aes(Stage, Dana.miR.1.RA))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("Dana-miR-1-RA RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirdana1

mirmel957<-ggplot(rnaseq, aes(Stage, Dmel.mir.957.RA))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("Dmel-mir-957-RA RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirmel957

mirtur317<-ggplot(rnaseq, aes(Stage, tur.miR.317.3p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("tur-miR-317-3p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirtur317

mirtur10<-ggplot(rnaseq, aes(Stage, tur.miR.10.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("tur-miR-10-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirtur10

miraae317<-ggplot(rnaseq, aes(Stage, aae.miR.317))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("aae-miR-317 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
miraae317

mirmse317<-ggplot(rnaseq, aes(Stage, mse.miR.317))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("mse-miR-317 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirmse317


mirlmi8<-ggplot(rnaseq, aes(Stage, lmi.miR.8.3p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("lmi-miR-8-3p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirlmi8

mirbta22<-ggplot(rnaseq, aes(Stage, bta.miR.22.3p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("bta-miR-22-3p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirbta22

mirmselet7a<-ggplot(rnaseq, aes(Stage, mse.let.7a))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("mse-let7a RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirmselet7a

mirchi22<-ggplot(rnaseq, aes(Stage, chi.miR.22.3p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("chi-miR-22-3p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirchi22

mirapi8<-ggplot(rnaseq, aes(Stage, api.miR.8))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("api-miR-8 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirapi8

miraga92b<-ggplot(rnaseq, aes(Stage, aga.miR.92b))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("aga-miR-92b RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
miraga92b

####mammal data boxplots####
mirbta2904<-ggplot(rnaseq, aes(Stage, bta.miR.2904))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("bta-miR-2904 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirbta2904

mirhsa141<-ggplot(rnaseq, aes(Stage, hsa.miR.141.3p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("hsa-miR-141-3p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirhsa141

mirhsa25<-ggplot(rnaseq, aes(Stage, hsa.miR.25.3p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("hsa-miR-25-3p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirhsa25

mirggo148a<-ggplot(rnaseq, aes(Stage, ggo.miR.148a))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("ggo-miR-148a RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirggo148a


mirchi21<-ggplot(rnaseq, aes(Stage, chi.miR.21.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("chi-miR-21-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirchi21

miroar21<-ggplot(rnaseq, aes(Stage, oar.miR.21))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("oar-miR-21 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
miroar21

mirptr26b<-ggplot(rnaseq, aes(Stage, ptr.miR.26b))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("ptr-miR-26b RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirptr26b

mirbta21<-ggplot(rnaseq, aes(Stage, bta.miR.21.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("bta-miR-21-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirbta21

mirgga30c<-ggplot(rnaseq, aes(Stage, gga.miR.30c.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("gga-miR-30c-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirgga30c

mirmmulet7i<-ggplot(rnaseq, aes(Stage, mmu.let.7i.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("mmu-let-7i-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirmmulet7i

mirchi99a<-ggplot(rnaseq, aes(Stage, chi.miR.99a.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("chi-miR-99a-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirchi99a

mirmml130a<-ggplot(rnaseq, aes(Stage, mml.miR.130a.3p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("mml-miR-130a-3p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirmml130a

mircgr192<-ggplot(rnaseq, aes(Stage, cgr.miR.192))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("cgr-miR-192 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mircgr192

mirgga194<-ggplot(rnaseq, aes(Stage, gga.miR.194))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("gga-miR-194 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirgga194

mirggolet7f<-ggplot(rnaseq, aes(Stage, ggo.let.7f))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("ggo-let7f RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirggolet7f

mirgga30a<-ggplot(rnaseq, aes(Stage, gga.miR.30a.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("gga-miR-30a-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirgga30a

mirmdo33<-ggplot(rnaseq, aes(Stage, mdo.miR.22.3p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("mdo-miR-22-3p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirmdo33

mirmdo186<-ggplot(rnaseq, aes(Stage, mdo.miR.186.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("mdo-miR-186-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirmdo186

mircfa486<-ggplot(rnaseq, aes(Stage, cfa.miR.486))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("cfa-miR-486 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mircfa486

miraca101<-ggplot(rnaseq, aes(Stage, aca.miR.101.3p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("aca-miR-101-3p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
miraca101

mircfa199<-ggplot(rnaseq, aes(Stage, cfa.miR.199))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("cfa-miR-199 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mircfa199

miracalet7f<-ggplot(rnaseq, aes(Stage, aca.let.7f.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("aca-let7f-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
miracalet7f

miroar30a<-ggplot(rnaseq, aes(Stage, oar.miR.30a.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("oar-miR-30a-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
miroar30a

mirchi26b<-ggplot(rnaseq, aes(Stage, chi.miR.26b.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("chi-miR-26b-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirchi26b

mirmmu99b<-ggplot(rnaseq, aes(Stage, mmu.miR.99b.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("mmu-miR-99b-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirmmu99b

mirhsa1246<-ggplot(rnaseq, aes(Stage, hsa.miR.1246))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("hsa-miR-1246 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirhsa1246


mirhsa151a<-ggplot(rnaseq, aes(Stage, hsa.miR.151a.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("hsa-miR-151a-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirhsa151a

mirtgu142<-ggplot(rnaseq, aes(Stage, tgu.miR.142.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("tgu-miR-142-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirtgu142

mirxtr10b<-ggplot(rnaseq, aes(Stage, xtr.miR.10b))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("xtr-miR-10b RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirxtr10b

mirmdo26<-ggplot(rnaseq, aes(Stage, mdo.miR.26.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("mdo-miR-26-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirmdo26

mirola192<-ggplot(rnaseq, aes(Stage, ola.miR.192.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("ola-miR-192-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirola192

mirssa30c<-ggplot(rnaseq, aes(Stage, ssa.miR.30c.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("ssa-miR-30c-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirssa30c

mirmne184<-ggplot(rnaseq, aes(Stage, mne.miR.184))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("mne-miR-184 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirmne184

mircfa125a<-ggplot(rnaseq, aes(Stage, cfa.miR.125a))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("cfa-miR-125a RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mircfa125a

mirchi22<-ggplot(rnaseq, aes(Stage, chi.miR.22.3p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("chi.miR.22.3p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirchi22

mirbta191<-ggplot(rnaseq, aes(Stage, bta.miR.191))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("bta-miR-191 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirbta191

mirchi100<-ggplot(rnaseq, aes(Stage, chi.miR.100.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("chi-miR-100-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirchi100

mirchi125b<-ggplot(rnaseq, aes(Stage, chi.miR.125b.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("chi-miR-125b-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirchi125b

mirppy92<-ggplot(rnaseq, aes(Stage, ppy.miR.92))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("ppy-miR-92 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirppy92

mirssc21<-ggplot(rnaseq, aes(Stage, ssc.miR.21))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("ssc-miR-21 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirssc21

mireca423<-ggplot(rnaseq, aes(Stage, eca.miR.423.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("eca-miR-423-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mireca423

mirbta30a<-ggplot(rnaseq, aes(Stage, bta.miR.30a.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("bta-miR-30a-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirbta30a

mirbta22<-ggplot(rnaseq, aes(Stage, bta.miR.22.3p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("bta-miR-22-3p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirbta22

mircgr1260<-ggplot(rnaseq, aes(Stage, cgr.miR.1260))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("cgr-miR-1260 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mircgr1260

mirhsa186<-ggplot(rnaseq, aes(Stage, hsa.miR.186.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("hsa-miR-186-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirhsa186

mirmmu27b<-ggplot(rnaseq, aes(Stage, mmu.miR.27b.3p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("mmu-miR-27b-3p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirmmu27b

mireca126<-ggplot(rnaseq, aes(Stage, eca.miR.126.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("eca-miR-126-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mireca126

mirhsa143<-ggplot(rnaseq, aes(Stage, hsa.miR.143.3p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("hsa-miR-143-3p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirhsa143

miraca9<-ggplot(rnaseq, aes(Stage, aca.miR.9.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("aca-miR-9-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
miraca9

miroan133<-ggplot(rnaseq, aes(Stage, oan.miR.133a.3p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("oan-miR-133a-3p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
miroan133

mirppa16<-ggplot(rnaseq, aes(Stage, ppa.miR.16))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("ppa-miR-16 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirppa16

miroha125b<-ggplot(rnaseq, aes(Stage, oha.miR.125b.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("oha-miR-125b-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
miroha125b

mirccr143<-ggplot(rnaseq, aes(Stage, ccr.miR.143))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("ccr-miR-143 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirccr143

mirtch30a<-ggplot(rnaseq, aes(Stage, tch.miR.30a.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("tch-miR-30a-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirtch30a

mirtch26a<-ggplot(rnaseq, aes(Stage, tch.miR.26a.5p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("tch-miR-26a-5p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirtch26a

mirssc192<-ggplot(rnaseq, aes(Stage, ssc.miR.192))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("ssc-miR-192 RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirssc192

miroan148<-ggplot(rnaseq, aes(Stage, oan.miR.148.3p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("oan-miR-148-3p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
miroan148

mirxtr22<-ggplot(rnaseq, aes(Stage, xtr.miR.22.3p))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("xtr.miR.22.3p RNAseq")+ylab("Normalized Expression")+xlab("")+
  theme(axis.text = element_text(colour="black", size=12, angle=45, hjust=0.9), axis.title = element_text( size=16), title = element_text(size=18))
mirxtr22

####multiple boxplots for RNAseq qPCR comparisons####
library(ggpubr)
ggarrange(mbantamq, miranabantam, miraaebantam, mircqubantam, ncol=2, nrow=2, labels="AUTO", font.label=list(size=24, face="bold"))

ggarrange(mbantamq, miranabantam, ncol=2, nrow=1, labels="AUTO", font.label=list(size=24, face="bold"))

ggarrange(mlet7q, mirdanalet7, mirhmelet7, mirmselet7a, ncol=2, nrow=2, labels="AUTO", font.label=list(size=24, face="bold"))

ggarrange(m277q, m277oq, mirdana277, mirhme277, ncol=2, nrow=2, labels="AUTO", font.label=list(size=24, face="bold"))

ggarrange(m92bq, mirdana92b, miraga92b, ncol=3, nrow=1, labels="AUTO", font.label=list(size=24, face="bold"))

ggarrange(m92bq,miraga92b, ncol=2, nrow=1, labels="AUTO", font.label=list(size=24, face="bold"))

ggarrange(m10q, mir10, mirdana10, mirmel10, ncol=2, nrow=2, labels="AUTO", font.label=list(size=24, face="bold"))

ggarrange(m317q, mir317, ncol=2, nrow=1, labels="AUTO", font.label=list(size=24, face="bold"))

ggarrange(m1q, mirdana1, ncol=2, nrow=1, labels="AUTO", font.label=list(size=24, face="bold"))

ggarrange(m8q, mir8, mirdvir8.5,  ncol=3, nrow=1, labels="AUTO", font.label=list(size=24, face="bold"))

ggarrange(m31aq, mir31a, ncol=2, nrow=1, labels="AUTO", font.label=list(size=24, face="bold"))

ggarrange(m184q, mir184, mirbfl184,  ncol=3, nrow=1, labels="AUTO", font.label=list(size=24, face="bold"))

ggarrange(m276bq, mirdana276b, ncol=2, nrow=1, labels="AUTO", font.label=list(size=24, face="bold"))

ggarrange(m305q, mirdvir305, miraae305, mirhme305, mirdana305, ncol=3, nrow=2, labels="AUTO", font.label=list(size=24, face="bold"))

ggarrange(m957q, mirmel957, ncol=2, nrow=1, labels="AUTO", font.label=list(size=24, face="bold"))











###boxplots for vertebrate data specifically
####all vertebrate data together####
library(reshape2)
library(ggplot2)
vert_data<-read.csv("boxplots/vertdata.csv", as.is=T)
str(vert_data)

vert_data$Stage<-factor(vert_data$Stage, levels=c("Feeding 3rd Instar", "Early Postfeeding", "Late Postfeeding", "Early Puparial", "Mid Puparial 1", "Mid Puparial 2", "Late Puparial"))
x<-melt(vert_data, id="Stage" )
head(x)
str(x)

ggplot(data=x, aes(x=Stage, y=value))+geom_boxplot(fill="slateblue", alpha=0.2)+ylim(0,5)+ggtitle("Vertebrate miRNA")+ylab("Normalized Expression")  

####bos taurus plot####
library(ggplot2)
library(reshape2)
btadata<-read.csv("boxplots/btadata.csv", as.is=T)
btadata$Stage<-factor(btadata$Stage, levels=c("Feeding 3rd Instar", "Early Postfeeding", "Late Postfeeding", "Early Puparial", "Mid Puparial 1", "Mid Puparial 2", "Late Puparial"))
str(btadata)

y<-melt(btadata, id="Stage" )
ggplot(data=y, aes(x=Stage, y=value))+geom_boxplot(fill="slateblue", alpha=0.2)+ggtitle("Bos taurus miRNA")+ylab("Normalized Expression")


