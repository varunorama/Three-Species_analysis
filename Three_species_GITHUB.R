###############################################################################################
library(readr)
library(DESeq2)
library(HTSFilter)
library(data.table)
d = "F:/LabWork/Three-Species/"
setwd(d)

three_species_annot = fread("Trinotate_annotation/AMEX_ThreeSpecies_associations_wTrinotate_allSPECIES.txt")

amex.eData = read.table("Ax_RSEM.gene.counts.matrix", row.names = 1, header = TRUE, sep='\t')

amex.eData = amex.eData[which(rownames(amex.eData)%in%three_species_annot$AxREFv2),]

#amex.eData = amex.eData[,c(-3,-8)]

amex.pData = data.frame("samp.names" = c("AMEX_D0_1","AMEX_D0_2","AMEX_D0_3",
                                         "AMEX_D1_1","AMEX_D1_2","AMEX_D1_3",
                                         "AMEX_D2_1","AMEX_D2_2","AMEX_D2_3",
                                         "AMEX_D3_1","AMEX_D3_2","AMEX_D3_3",
                                         "AMEX_D4_1","AMEX_D4_2","AMEX_D4_3",
                                         "AMEX_D5_1","AMEX_D5_2","AMEX_D5_3"),
                        "dpa" = as.factor(c(rep(0,3),rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3))))
rownames(amex.pData) = colnames(amex.eData)

#amex.dds <- amex.dds[rowSums(fpm(amex.dds)>0)>=18]
register(bpstart(SnowParam(6)))
##### AMEX - generate DESeq set  #######
amex.dds <- DESeqDataSetFromMatrix(countData = round(amex.eData[,c(1,2,4,6,7,9,10,12,13,14,17,18)]), colData = amex.pData[c(1,2,4,6,7,9,10,12,13,14,17,18),], design = ~dpa)
amex_dds2 <- DESeq(amex.dds, test=c("LRT"), reduced =~1, parallel = TRUE)
#amex_filtered <- HTSFilter(amex_dds2, s.len=100, normalization = "DESeq", parallel = TRUE)$filteredData
#amex_LRT = results(amex_dds2, independentFiltering=TRUE)
#hist(amex_LRT$pvalue)
#summary(amex_LRT)

###
amex.dds_1v0 <- DESeqDataSetFromMatrix(countData = round(amex.eData[,c(1,2,4,6)]), colData = amex.pData[c(1,2,4,6),], design = ~dpa)
amex.dds_1v0 <- amex.dds_1v0[rowSums(fpm(amex.dds_1v0)>0)>=4]
amex_dds2_D1vs0 <- DESeq(amex.dds_1v0, test=c("Wald"), parallel = TRUE)
amex_D1vD0 = results(amex_dds2_D1vs0, contrast = c("dpa","1","0"), independentFiltering=TRUE)
hist(amex_D1vD0$pvalue)
summary(amex_D1vD0, alpha = 0.05)

amex.dds_2v0 <- DESeqDataSetFromMatrix(countData = round(amex.eData[,c(c(1,2,7,9))]), colData = amex.pData[c(c(1,2,7,9)),], design = ~dpa)
amex.dds_2v0 <- amex.dds_2v0[rowSums(fpm(amex.dds_2v0)>0)>=4]
amex_dds2_D2vs0 <- DESeq(amex.dds_2v0, test=c("Wald"), parallel = TRUE)
amex_D2vD0 = results(amex_dds2_D2vs0, contrast = c("dpa","2","0"), independentFiltering=TRUE)
hist(amex_D2vD0$pvalue)
summary(amex_D2vD0, alpha = 0.05)

amex.dds_3v0 <- DESeqDataSetFromMatrix(countData = round(amex.eData[,c(1:2,10,12)]), colData = amex.pData[c(1:2,10,12),], design = ~dpa)
amex.dds_3v0 <- amex.dds_3v0[rowSums(fpm(amex.dds_3v0)>0)>=4]
amex_dds2_D3vs0 <- DESeq(amex.dds_3v0, test=c("Wald"), parallel = TRUE)
amex_D3vD0 = results(amex_dds2_D3vs0, contrast = c("dpa","3","0"), independentFiltering=TRUE)
hist(amex_D3vD0$pvalue)
summary(amex_D3vD0, alpha = 0.05)

amex.dds_4v0 <- DESeqDataSetFromMatrix(countData = round(amex.eData[,c(1:2,13:14)]), colData = amex.pData[c(1:2,13:14),], design = ~dpa)
amex.dds_4v0 <- amex.dds_4v0[rowSums(fpm(amex.dds_4v0)>0)>=4]
amex_dds2_D4vs0 <- DESeq(amex.dds_4v0, test=c("Wald"), parallel = TRUE)
amex_D4vD0 = results(amex_dds2_D4vs0, contrast = c("dpa","4","0"), independentFiltering=TRUE)
hist(amex_D4vD0$pvalue)
summary(amex_D4vD0, alpha = 0.05)

amex.dds_5v0 <- DESeqDataSetFromMatrix(countData = round(amex.eData[,c(1:2,17:18)]), colData = amex.pData[c(1:2,17:18),], design = ~dpa)
amex.dds_5v0 <- amex.dds_5v0[rowSums(fpm(amex.dds_5v0)>0)>=4]
amex_dds2_D5vs0 <- DESeq(amex.dds_5v0, test=c("Wald"), parallel = TRUE)
amex_D5vD0 = results(amex_dds2_D5vs0, contrast = c("dpa","5","0"), independentFiltering=TRUE)
hist(amex_D5vD0$pvalue)
summary(amex_D5vD0, alpha = 0.05)
###
amex_D1vD0_df = data.frame(amex_D1vD0)
amex_D1vD0_df$txid = rownames(amex_D1vD0_df)
amex_D2vD0_df = data.frame(amex_D2vD0)
amex_D2vD0_df$txid = rownames(amex_D2vD0_df)
amex_D3vD0_df = data.frame(amex_D3vD0)
amex_D3vD0_df$txid = rownames(amex_D3vD0_df)
amex_D4vD0_df = data.frame(amex_D4vD0)
amex_D4vD0_df$txid = rownames(amex_D4vD0_df)
amex_D5vD0_df = data.frame(amex_D5vD0)
amex_D5vD0_df$txid = rownames(amex_D5vD0_df)
###

amex_independent_tests = union(union(union(union(union(rownames(amex_D1vD0[which(amex_D1vD0$padj < 0.05),]),
                                                       rownames(amex_D2vD0[which(amex_D2vD0$padj < 0.05),])),
                                                 rownames(amex_D3vD0[which(amex_D3vD0$padj < 0.05),])),
                                           rownames(amex_D4vD0[which(amex_D4vD0$padj < 0.05),])),
                                     rownames(amex_D4vD0[which(amex_D4vD0$padj < 0.05),])),
                               rownames(amex_D5vD0[which(amex_D5vD0$padj < 0.05),]))

amex_DEG = data.frame(trans = amex_independent_tests)
rownames(amex_DEG) = amex_DEG[,1]
amex_DEG$D1v0_FC = amex_D1vD0_df[match(rownames(amex_DEG),rownames(amex_D1vD0_df)),c(2)]
amex_DEG$D2v0_FC = amex_D2vD0_df[match(rownames(amex_DEG),rownames(amex_D2vD0_df)),c(2)]
amex_DEG$D3v0_FC = amex_D3vD0_df[match(rownames(amex_DEG),rownames(amex_D3vD0_df)),c(2)]
amex_DEG$D4v0_FC = amex_D4vD0_df[match(rownames(amex_DEG),rownames(amex_D4vD0_df)),c(2)]
amex_DEG$D5v0_FC = amex_D5vD0_df[match(rownames(amex_DEG),rownames(amex_D5vD0_df)),c(2)]
amex_DEG$D1v0_FDR = amex_D1vD0_df[match(rownames(amex_DEG),rownames(amex_D1vD0_df)),c(6)]
amex_DEG$D2v0_FDR = amex_D2vD0_df[match(rownames(amex_DEG),rownames(amex_D2vD0_df)),c(6)]
amex_DEG$D3v0_FDR = amex_D3vD0_df[match(rownames(amex_DEG),rownames(amex_D3vD0_df)),c(6)]
amex_DEG$D4v0_FDR = amex_D4vD0_df[match(rownames(amex_DEG),rownames(amex_D4vD0_df)),c(6)]
amex_DEG$D5v0_FDR = amex_D5vD0_df[match(rownames(amex_DEG),rownames(amex_D5vD0_df)),c(6)]

write.table(amex_DEG, "DEG_list/Amex_DEG_fulllist.txt",row.names = TRUE, col.names = TRUE, sep='\t')
pearson_microarray<-cor(fpm(amex_dds2)[rownames(counts(amex_dds2, norm=TRUE))%in%amex_independent_tests,],method="spearman")
diff_pearson<-1-pearson_microarray
dist_pearson_micro<-as.dist(diff_pearson)
pearson_micro_clust<-hclust(dist_pearson_micro,method="average")
plot(pearson_micro_clust)

amex_fpm_norm = data.frame(fpm(amex_dds2)[rownames(counts(amex_dds2, norm=TRUE))%in%amex_independent_tests,])
amex_fpm_norm$D0_avg = rowMeans(amex_fpm_norm[,c(1:2)])
amex_fpm_norm$D1_avg = rowMeans(amex_fpm_norm[,c(3:4)])
amex_fpm_norm$D2_avg = rowMeans(amex_fpm_norm[,c(5:6)])
amex_fpm_norm$D3_avg = rowMeans(amex_fpm_norm[,c(7:8)])
amex_fpm_norm$D4_avg = rowMeans(amex_fpm_norm[,c(9:10)])
amex_fpm_norm$D5_avg = rowMeans(amex_fpm_norm[,c(11:12)])

library(pvclust)
y = pvclust(data=amex_fpm_norm[,c(13:18)],method.hclust="complete",method.dist="correlation",nboot=1000, parallel = TRUE)
tiff(filename = "figures/AMEX_TimepointCluster_Bootstrap.tiff", width = 6, height = 4, units = "in", res = 300)
plot(y)
pvrect(y, alpha=0.94)
dev.off()

#plotMDS(fpm(amex_dds2)[rownames(counts(amex_dds2, norm=TRUE))%in%amex_independent_tests,])

nrow(amex_D1vD0[which(amex_D1vD0$padj < 0.05 & amex_D1vD0$log2FoldChange > 0),])
nrow(amex_D2vD0[which(amex_D2vD0$padj < 0.05& amex_D2vD0$log2FoldChange > 0),])
nrow(amex_D3vD0[which(amex_D3vD0$padj < 0.05& amex_D3vD0$log2FoldChange > 0),])
nrow(amex_D4vD0[which(amex_D4vD0$padj < 0.05& amex_D4vD0$log2FoldChange > 0),])
nrow(amex_D5vD0[which(amex_D5vD0$padj < 0.05& amex_D5vD0$log2FoldChange > 0),])

nrow(amex_D1vD0[which(amex_D1vD0$padj < 0.05& amex_D1vD0$log2FoldChange < 0),])
nrow(amex_D2vD0[which(amex_D2vD0$padj < 0.05& amex_D2vD0$log2FoldChange < 0),])
nrow(amex_D3vD0[which(amex_D3vD0$padj < 0.05& amex_D3vD0$log2FoldChange < 0),])
nrow(amex_D4vD0[which(amex_D4vD0$padj < 0.05& amex_D4vD0$log2FoldChange < 0),])
nrow(amex_D5vD0[which(amex_D5vD0$padj < 0.05& amex_D5vD0$log2FoldChange < 0),])


write.table(amex_D1vD0_df, "amex_D1vD0_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(amex_D2vD0_df, "amex_D2vD0_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(amex_D3vD0_df, "amex_D3vD0_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(amex_D4vD0_df, "amex_D4vD0_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(amex_D5vD0_df, "amex_D5vD0_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')

write.table(amex_D1vD0_df[which(amex_D1vD0_df$padj<0.05),], "amex_D1vD0_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(amex_D2vD0_df[which(amex_D2vD0_df$padj<0.05),], "amex_D2vD0_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(amex_D3vD0_df[which(amex_D3vD0_df$padj<0.05),], "amex_D3vD0_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(amex_D4vD0_df[which(amex_D4vD0_df$padj<0.05),], "amex_D4vD0_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(amex_D5vD0_df[which(amex_D5vD0_df$padj<0.05),], "amex_D5vD0_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')


###########


am1 = rownames(amex_D1vD0_df[which(amex_D1vD0_df$padj<0.05),])
am2 = unique(c(rownames(amex_D2vD0_df[which(amex_D2vD0_df$padj<0.05),]),am1))
am3 = unique(c(rownames(amex_D3vD0_df[which(amex_D3vD0_df$padj<0.05),]),am2))
am4 = unique(c(rownames(amex_D4vD0_df[which(amex_D4vD0_df$padj<0.05),]),am3))
am5 = unique(c(rownames(amex_D5vD0_df[which(amex_D5vD0_df$padj<0.05),]),am4))

c(length(am1), length(am2), length(am3), length(am4), length(am5))

#save.image(file = "ThreeSpecies.RData")
########################################
########## AAND ########################
########################################

library(readr)
library(DESeq2)
library(HTSFilter)
d = "F:/LabWork/Three-Species/"
setwd(d)

and.eData = read.table("And_RSEM.gene.counts.matrix", row.names = 1, header = TRUE, sep='\t')

and.eData = and.eData[which(rownames(and.eData)%in%three_species_annot$AndREFv2),]

and.pData = data.frame("samp.names" = c("And_D0_1","And_D0_2","And_D0_3",
                                        "And_D1_1","And_D1_2","And_D1_3",
                                        "And_D2_1","And_D2_2","And_D2_3",
                                        "And_D3_1","And_D3_2","And_D3_3",
                                        "And_D4_1","And_D4_2","And_D4_3",
                                        "And_D5_1","And_D5_2","And_D5_3"),
                       "dpa" = as.factor(c(rep(0,3),rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3))))
rownames(and.pData) = colnames(and.eData)
and.dds <- DESeqDataSetFromMatrix(countData = round(and.eData[,c(1:3,4,6,7:9,10,12,13,15,16,17)]), colData = and.pData[c(1:3,4,6,7:9,10,12,13,15,16,17),], design = ~dpa)
#and.dds <- and.dds[rowSums(fpm(and.dds)>0)>=9]
register(bpstart(SnowParam(6)))
and_dds <- DESeq(and.dds, test=c("LRT"), reduced =~1, parallel = TRUE)
##### AAND - run comparisons  #######
#and_LRT = results(and_dds, independentFiltering=TRUE)
#hist(and_LRT$pvalue)
#summary(and_LRT, alpha = 0.05)

###
and.dds_1v0 <- DESeqDataSetFromMatrix(countData = round(and.eData[,c(1:3,4,6)]), colData = and.pData[c(1:3,4,6),], design = ~dpa)
and.dds_1v0 <- and.dds_1v0[rowSums(fpm(and.dds_1v0)>0)>=5]
and_dds_D1vs0 <- DESeq(and.dds_1v0, test=c("Wald"), parallel = TRUE)
and_D1vD0 = results(and_dds_D1vs0, contrast = c("dpa","1","0"), independentFiltering=TRUE)
hist(and_D1vD0$pvalue)
summary(and_D1vD0, alpha = 0.05)

and.dds_2v0 <- DESeqDataSetFromMatrix(countData = round(and.eData[,c(1:3,7:9)]), colData = and.pData[c(1:3,7:9),], design = ~dpa)
and.dds_2v0 <- and.dds_2v0[rowSums(fpm(and.dds_2v0)>0)>=6]
and_dds_D2vs0 <- DESeq(and.dds_2v0, test=c("Wald"), parallel = TRUE)
and_D2vD0 = results(and_dds_D2vs0, contrast = c("dpa","2","0"), independentFiltering=TRUE)
hist(and_D2vD0$pvalue)
summary(and_D2vD0, alpha = 0.05)

and.dds_3v0 <- DESeqDataSetFromMatrix(countData = round(and.eData[,c(1:3,10,12)]), colData = and.pData[c(1:3,10,12),], design = ~dpa)
and.dds_3v0 <- and.dds_3v0[rowSums(fpm(and.dds_3v0)>0)>=5]
and_dds_D3vs0 <- DESeq(and.dds_3v0, test=c("Wald"), parallel = TRUE)
and_D3vD0 = results(and_dds_D3vs0, contrast = c("dpa","3","0"), independentFiltering=TRUE)
hist(and_D3vD0$pvalue)
summary(and_D3vD0, alpha = 0.05)


and.dds_4v0 <- DESeqDataSetFromMatrix(countData = round(and.eData[,c(1:3,13,15)]), colData = and.pData[c(1:3,13,15),], design = ~dpa)
and.dds_4v0 <- and.dds_4v0[rowSums(fpm(and.dds_4v0)>0)>=5]
and_dds_D4vs0 <- DESeq(and.dds_4v0, test=c("Wald"), parallel = TRUE)
and_D4vD0 = results(and_dds_D4vs0, contrast = c("dpa","4","0"), independentFiltering=TRUE)
hist(and_D4vD0$pvalue)
summary(and_D4vD0, alpha = 0.05)


and.dds_5v0 <- DESeqDataSetFromMatrix(countData = round(and.eData[,c(1:3,16,17)]), colData = and.pData[c(1:3,16,17),], design = ~dpa)
and.dds_5v0 <- and.dds_5v0[rowSums(fpm(and.dds_5v0)>0)>=5]
and_dds_D5vs0 <- DESeq(and.dds_5v0, test=c("Wald"), parallel = TRUE)
and_D5vD0 = results(and_dds_D5vs0, contrast = c("dpa","5","0"), independentFiltering=TRUE)
hist(and_D5vD0$pvalue)
summary(and_D5vD0, alpha = 0.05)


and_independent_tests = union(union(union(union(union(rownames(and_D1vD0[which(and_D1vD0$padj < 0.05),]),
                                                      rownames(and_D2vD0[which(and_D2vD0$padj < 0.05),])),
                                                rownames(and_D3vD0[which(and_D3vD0$padj < 0.05),])),
                                          rownames(and_D4vD0[which(and_D4vD0$padj < 0.05),])),
                                    rownames(and_D4vD0[which(and_D4vD0$padj < 0.05),])),
                              rownames(and_D5vD0[which(and_D5vD0$padj < 0.05),]))


and_D1vD0_df = data.frame(and_D1vD0)
and_D1vD0_df$txid = rownames(and_D1vD0_df)
and_D2vD0_df = data.frame(and_D2vD0)
and_D2vD0_df$txid = rownames(and_D2vD0_df)
and_D3vD0_df = data.frame(and_D3vD0)
and_D3vD0_df$txid = rownames(and_D3vD0_df)
and_D4vD0_df = data.frame(and_D4vD0)
and_D4vD0_df$txid = rownames(and_D4vD0_df)
and_D5vD0_df = data.frame(and_D5vD0)
and_D5vD0_df$txid = rownames(and_D5vD0_df)

and_DEG = data.frame(trans = and_independent_tests)
rownames(and_DEG) = and_DEG[,1]
and_DEG$D1v0_FC = and_D1vD0_df[match(rownames(and_DEG),rownames(and_D1vD0_df)),c(2)]
and_DEG$D2v0_FC = and_D2vD0_df[match(rownames(and_DEG),rownames(and_D2vD0_df)),c(2)]
and_DEG$D3v0_FC = and_D3vD0_df[match(rownames(and_DEG),rownames(and_D3vD0_df)),c(2)]
and_DEG$D4v0_FC = and_D4vD0_df[match(rownames(and_DEG),rownames(and_D4vD0_df)),c(2)]
and_DEG$D5v0_FC = and_D5vD0_df[match(rownames(and_DEG),rownames(and_D5vD0_df)),c(2)]
and_DEG$D1v0_FDR = and_D1vD0_df[match(rownames(and_DEG),rownames(and_D1vD0_df)),c(6)]
and_DEG$D2v0_FDR = and_D2vD0_df[match(rownames(and_DEG),rownames(and_D2vD0_df)),c(6)]
and_DEG$D3v0_FDR = and_D3vD0_df[match(rownames(and_DEG),rownames(and_D3vD0_df)),c(6)]
and_DEG$D4v0_FDR = and_D4vD0_df[match(rownames(and_DEG),rownames(and_D4vD0_df)),c(6)]
and_DEG$D5v0_FDR = and_D5vD0_df[match(rownames(and_DEG),rownames(and_D5vD0_df)),c(6)]

nrow(and_D1vD0[which(and_D1vD0$padj < 0.05 & and_D1vD0$log2FoldChange > 0),])
nrow(and_D2vD0[which(and_D2vD0$padj < 0.05& and_D2vD0$log2FoldChange > 0),])
nrow(and_D3vD0[which(and_D3vD0$padj < 0.05& and_D3vD0$log2FoldChange > 0),])
nrow(and_D4vD0[which(and_D4vD0$padj < 0.05& and_D4vD0$log2FoldChange > 0),])
nrow(and_D5vD0[which(and_D5vD0$padj < 0.05& and_D5vD0$log2FoldChange > 0),])

nrow(and_D1vD0[which(and_D1vD0$padj < 0.05& and_D1vD0$log2FoldChange < 0),])
nrow(and_D2vD0[which(and_D2vD0$padj < 0.05& and_D2vD0$log2FoldChange < 0),])
nrow(and_D3vD0[which(and_D3vD0$padj < 0.05& and_D3vD0$log2FoldChange < 0),])
nrow(and_D4vD0[which(and_D4vD0$padj < 0.05& and_D4vD0$log2FoldChange < 0),])
nrow(and_D5vD0[which(and_D5vD0$padj < 0.05& and_D5vD0$log2FoldChange < 0),])

write.table(and_D1vD0_df, "and_D1vD0_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(and_D2vD0_df, "and_D2vD0_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(and_D3vD0_df, "and_D3vD0_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(and_D4vD0_df, "and_D4vD0_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(and_D5vD0_df, "and_D5vD0_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')

write.table(and_D1vD0_df[which(and_D1vD0_df$padj<0.05),], "and_D1vD0_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(and_D2vD0_df[which(and_D2vD0_df$padj<0.05),], "and_D2vD0_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(and_D3vD0_df[which(and_D3vD0_df$padj<0.05),], "and_D3vD0_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(and_D4vD0_df[which(and_D4vD0_df$padj<0.05),], "and_D4vD0_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(and_D5vD0_df[which(and_D5vD0_df$padj<0.05),], "and_D5vD0_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')


ad1 = rownames(and_D1vD0_df[which(and_D1vD0_df$padj<0.05),])
ad2 = unique(c(rownames(and_D2vD0_df[which(and_D2vD0_df$padj<0.05),]),ad1))
ad3 = unique(c(rownames(and_D3vD0_df[which(and_D3vD0_df$padj<0.05),]),ad2))
ad4 = unique(c(rownames(and_D4vD0_df[which(and_D4vD0_df$padj<0.05),]),ad3))
ad5 = unique(c(rownames(and_D5vD0_df[which(and_D5vD0_df$padj<0.05),]),ad4))

c(length(ad1), length(ad2), length(ad3), length(ad4), length(ad5))

write.table(and_DEG, "DEG_list/and_DEG_fulllist.txt",row.names = TRUE, col.names = TRUE, sep='\t')

pearson_microarray<-cor(fpm(and_dds)[rownames(counts(and_dds, norm=TRUE))%in%and_independent_tests,],method="spearman")
diff_pearson<-1-pearson_microarray
dist_pearson_micro<-as.dist(diff_pearson)
pearson_micro_clust<-hclust(dist_pearson_micro,method="average")
plot(pearson_micro_clust)

and_fpm_norm = data.frame(fpm(and_dds)[rownames(counts(and_dds, norm=TRUE))%in%and_independent_tests,])
and_fpm_norm$D0_avg = rowMeans(and_fpm_norm[,c(1:3)])
and_fpm_norm$D1_avg = rowMeans(and_fpm_norm[,c(4:5)])
and_fpm_norm$D2_avg = rowMeans(and_fpm_norm[,c(6:8)])
and_fpm_norm$D3_avg = rowMeans(and_fpm_norm[,c(9:10)])
and_fpm_norm$D4_avg = rowMeans(and_fpm_norm[,c(11:12)])
and_fpm_norm$D5_avg = rowMeans(and_fpm_norm[,c(13:14)])


library(pvclust)
y = pvclust(data=and_fpm_norm[,c(15:20)],method.hclust="complete",method.dist="correlation",nboot=1000, parallel = TRUE)
tiff(filename = "figures/AAND_TimepointCluster_Bootstrap.tiff", width = 6, height = 4, units = "in", res = 300)
plot(y)
pvrect(y, alpha=0.94)
dev.off()



########################################
########## AMAC ########################
########################################

#load("ThreeSpecies.RData")
mac.eData = read.table("Mac_RSEM.gene.counts.matrix", row.names = 1, header = TRUE, sep='\t')

mac.eData = mac.eData[which(rownames(mac.eData)%in%three_species_annot$MacREFv2),]

#mac.eData = mac.eData[,c(-6)]
mac.pData = data.frame("samp.names" = c("Mac_D0_1","Mac_D0_2","Mac_D0_3",
                                        "Mac_D1_1","Mac_D1_2","Mac_D1_3",
                                        "Mac_D2_1","Mac_D2_2","Mac_D2_3",
                                        "Mac_D3_1","Mac_D3_2","Mac_D3_3",
                                        "Mac_D4_1","Mac_D4_2","Mac_D4_3",
                                        "Mac_D5_1","Mac_D5_2","Mac_D5_3"),
                       "dpa" = as.factor(c(rep(0,3),rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3))))
rownames(mac.pData) = colnames(mac.eData)
mac.dds <- DESeqDataSetFromMatrix(countData = round(mac.eData[,c(1:3,4,5,7:8,10:11,13:14,16:17)]), colData = mac.pData[c(1:3,4,5,7:8,10:11,13:14,16:17),], design = ~dpa)
#mac.dds <- mac.dds[rowSums(fpm(mac.dds)>0)>=9]
register(bpstart(SnowParam(6)))
mac_dds <- DESeq(mac.dds, test=c("LRT"), reduced =~1, parallel = TRUE )

##### AMAC - generate DESeq set  #######
mac_LRT = results(mac_dds, independentFiltering=TRUE)
hist(mac_LRT$pvalue)
summary(mac_LRT)


### comparison
mac.dds_1v0 <- DESeqDataSetFromMatrix(countData = round(mac.eData[,c(1:5)]), colData = mac.pData[c(1:5),], design = ~dpa)
mac.dds_1v0 <- mac.dds_1v0[rowSums(fpm(mac.dds_1v0)>0)>=5]
mac_dds_D1vs0 <- DESeq(mac.dds_1v0, test=c("Wald"), parallel = TRUE)
mac_D1vD0 = results(mac_dds_D1vs0, contrast = c("dpa","1","0"), independentFiltering=TRUE)
hist(mac_D1vD0$pvalue)
summary(mac_D1vD0, alpha = 0.05)

mac.dds_2v0 <- DESeqDataSetFromMatrix(countData = round(mac.eData[,c(1:3,7:8)]), colData = mac.pData[c(1:3,7:8),], design = ~dpa)
mac.dds_2v0 <- mac.dds_2v0[rowSums(fpm(mac.dds_2v0)>0)>=5]
mac_dds_D2vs0 <- DESeq(mac.dds_2v0, test=c("Wald"), parallel = TRUE)
mac_D2vD0 = results(mac_dds_D2vs0, contrast = c("dpa","2","0"), independentFiltering=TRUE)
hist(mac_D2vD0$pvalue)
summary(mac_D2vD0, alpha = 0.05)

mac.dds_3v0 <- DESeqDataSetFromMatrix(countData = round(mac.eData[,c(1:3,10:11)]), colData = mac.pData[c(1:3,10:11),], design = ~dpa)
mac.dds_3v0 <- mac.dds_3v0[rowSums(fpm(mac.dds_3v0)>0)>=5]
mac_dds_D3vs0 <- DESeq(mac.dds_3v0, test=c("Wald"), parallel = TRUE)
mac_D3vD0 = results(mac_dds_D3vs0, contrast = c("dpa","3","0"), independentFiltering=TRUE)
hist(mac_D3vD0$pvalue)
summary(mac_D3vD0, alpha = 0.05)

mac.dds_4v0 <- DESeqDataSetFromMatrix(countData = round(mac.eData[,c(1:3,13:14)]), colData = mac.pData[c(1:3,13:14),], design = ~dpa)
mac.dds_4v0 <- mac.dds_4v0[rowSums(fpm(mac.dds_4v0)>0)>=5]
mac_dds_D4vs0 <- DESeq(mac.dds_4v0, test=c("Wald"), parallel = TRUE)
mac_D4vD0 = results(mac_dds_D4vs0, contrast = c("dpa","4","0"), independentFiltering=TRUE)
hist(mac_D4vD0$pvalue)
summary(mac_D4vD0, alpha = 0.05)

mac.dds_5v0 <- DESeqDataSetFromMatrix(countData = round(mac.eData[,c(1:3,16:17)]), colData = mac.pData[c(1:3,16:17),], design = ~dpa)
mac.dds_5v0 <- mac.dds_5v0[rowSums(fpm(mac.dds_5v0)>0)>=5]
mac_dds_D5vs0 <- DESeq(mac.dds_5v0, test=c("Wald"), parallel = TRUE)
mac_D5vD0 = results(mac_dds_D5vs0, contrast = c("dpa","5","0"), independentFiltering=TRUE)
hist(mac_D5vD0$pvalue)
summary(mac_D5vD0, alpha = 0.05)

mac_independent_tests = union(union(union(union(union(rownames(mac_D1vD0[which(mac_D1vD0$padj < 0.05),]),
                                                      rownames(mac_D2vD0[which(mac_D2vD0$padj < 0.05),])),
                                                rownames(mac_D3vD0[which(mac_D3vD0$padj < 0.05),])),
                                          rownames(mac_D4vD0[which(mac_D4vD0$padj < 0.05),])),
                                    rownames(mac_D4vD0[which(mac_D4vD0$padj < 0.05),])),
                              rownames(mac_D5vD0[which(mac_D5vD0$padj < 0.05),]))

mac_D1vD0_df = data.frame(mac_D1vD0)
mac_D1vD0_df$txid = rownames(mac_D1vD0_df)
mac_D2vD0_df = data.frame(mac_D2vD0)
mac_D2vD0_df$txid = rownames(mac_D2vD0_df)
mac_D3vD0_df = data.frame(mac_D3vD0)
mac_D3vD0_df$txid = rownames(mac_D3vD0_df)
mac_D4vD0_df = data.frame(mac_D4vD0)
mac_D4vD0_df$txid = rownames(mac_D4vD0_df)
mac_D5vD0_df = data.frame(mac_D5vD0)
mac_D5vD0_df$txid = rownames(mac_D5vD0_df)

mac_DEG = data.frame(trans = mac_independent_tests)
rownames(mac_DEG) = mac_DEG[,1]
mac_DEG$D1v0_FC = mac_D1vD0_df[match(rownames(mac_DEG),rownames(mac_D1vD0_df)),c(2)]
mac_DEG$D2v0_FC = mac_D2vD0_df[match(rownames(mac_DEG),rownames(mac_D2vD0_df)),c(2)]
mac_DEG$D3v0_FC = mac_D3vD0_df[match(rownames(mac_DEG),rownames(mac_D3vD0_df)),c(2)]
mac_DEG$D4v0_FC = mac_D4vD0_df[match(rownames(mac_DEG),rownames(mac_D4vD0_df)),c(2)]
mac_DEG$D5v0_FC = mac_D5vD0_df[match(rownames(mac_DEG),rownames(mac_D5vD0_df)),c(2)]
mac_DEG$D1v0_FDR = mac_D1vD0_df[match(rownames(mac_DEG),rownames(mac_D1vD0_df)),c(6)]
mac_DEG$D2v0_FDR = mac_D2vD0_df[match(rownames(mac_DEG),rownames(mac_D2vD0_df)),c(6)]
mac_DEG$D3v0_FDR = mac_D3vD0_df[match(rownames(mac_DEG),rownames(mac_D3vD0_df)),c(6)]
mac_DEG$D4v0_FDR = mac_D4vD0_df[match(rownames(mac_DEG),rownames(mac_D4vD0_df)),c(6)]
mac_DEG$D5v0_FDR = mac_D5vD0_df[match(rownames(mac_DEG),rownames(mac_D5vD0_df)),c(6)]

write.table(mac_DEG, "DEG_list/mac_DEG_fulllist.txt",row.names = TRUE, col.names = TRUE, sep='\t')

pearson_microarray<-cor(fpm(mac_dds)[rownames(counts(mac_dds, norm=TRUE))%in%mac_independent_tests,],method="spearman")
diff_pearson<-1-pearson_microarray
dist_pearson_micro<-as.dist(diff_pearson)
pearson_micro_clust<-hclust(dist_pearson_micro,method="average")
plot(pearson_micro_clust)

mac_fpm_norm = data.frame(fpm(mac_dds)[rownames(counts(mac_dds, norm=TRUE))%in%mac_independent_tests,])
mac_fpm_norm$D0_avg = rowMeans(mac_fpm_norm[,c(1:3)])
mac_fpm_norm$D1_avg = rowMeans(mac_fpm_norm[,c(4:5)])
mac_fpm_norm$D2_avg = rowMeans(mac_fpm_norm[,c(6:7)])
mac_fpm_norm$D3_avg = rowMeans(mac_fpm_norm[,c(8:9)])
mac_fpm_norm$D4_avg = rowMeans(mac_fpm_norm[,c(10:11)])
mac_fpm_norm$D5_avg = rowMeans(mac_fpm_norm[,c(12:13)])


library(pvclust)
y = pvclust(data=mac_fpm_norm[,c(14:19)],method.hclust="complete",method.dist="correlation",nboot=1000, parallel = TRUE)
tiff(filename = "figures/AMAC_TimepointCluster_Bootstrap.tiff", width = 6, height = 4, units = "in", res = 300)
plot(y)
pvrect(y, alpha=0.94)
dev.off()

nrow(mac_D1vD0[which(mac_D1vD0$padj < 0.05 & mac_D1vD0$log2FoldChange > 0),])
nrow(mac_D2vD0[which(mac_D2vD0$padj < 0.05& mac_D2vD0$log2FoldChange > 0),])
nrow(mac_D3vD0[which(mac_D3vD0$padj < 0.05& mac_D3vD0$log2FoldChange > 0),])
nrow(mac_D4vD0[which(mac_D4vD0$padj < 0.05& mac_D4vD0$log2FoldChange > 0),])
nrow(mac_D5vD0[which(mac_D5vD0$padj < 0.05& mac_D5vD0$log2FoldChange > 0),])

nrow(mac_D1vD0[which(mac_D1vD0$padj < 0.05 & mac_D1vD0$log2FoldChange < 0),])
nrow(mac_D2vD0[which(mac_D2vD0$padj < 0.05& mac_D2vD0$log2FoldChange < 0),])
nrow(mac_D3vD0[which(mac_D3vD0$padj < 0.05& mac_D3vD0$log2FoldChange < 0),])
nrow(mac_D4vD0[which(mac_D4vD0$padj < 0.05& mac_D4vD0$log2FoldChange < 0),])
nrow(mac_D5vD0[which(mac_D5vD0$padj < 0.05& mac_D5vD0$log2FoldChange < 0),])

write.table(mac_D1vD0_df, "mac_D1vD0_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(mac_D2vD0_df, "mac_D2vD0_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(mac_D3vD0_df, "mac_D3vD0_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(mac_D4vD0_df, "mac_D4vD0_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(mac_D5vD0_df, "mac_D5vD0_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')

write.table(mac_D1vD0_df[which(mac_D1vD0_df$padj<0.05),], "mac_D1vD0_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(mac_D2vD0_df[which(mac_D2vD0_df$padj<0.05),], "mac_D2vD0_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(mac_D3vD0_df[which(mac_D3vD0_df$padj<0.05),], "mac_D3vD0_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(mac_D4vD0_df[which(mac_D4vD0_df$padj<0.05),], "mac_D4vD0_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(mac_D5vD0_df[which(mac_D5vD0_df$padj<0.05),], "mac_D5vD0_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')

mac1 = rownames(mac_D1vD0_df[which(mac_D1vD0_df$padj<0.05),])
mac2 = unique(c(rownames(mac_D2vD0_df[which(mac_D2vD0_df$padj<0.05),]),mac1))
mac3 = unique(c(rownames(mac_D3vD0_df[which(mac_D3vD0_df$padj<0.05),]),mac2))
mac4 = unique(c(rownames(mac_D4vD0_df[which(mac_D4vD0_df$padj<0.05),]),mac3))
mac5 = unique(c(rownames(mac_D5vD0_df[which(mac_D5vD0_df$padj<0.05),]),mac4))

c(length(mac1), length(mac2), length(mac3), length(mac4), length(mac5))

#################################################################################################

amex.dds_1v0 <- DESeqDataSetFromMatrix(countData = round(amex.eData[,c(1,2,4,6)]), colData = amex.pData[c(1,2,4,6),], design = ~dpa)
amex.dds_1v0 <- amex.dds_1v0[rowSums(fpm(amex.dds_1v0)>0)>=4]
amex_dds2_D1vs0 <- DESeq(amex.dds_1v0, test=c("Wald"), parallel = TRUE)
amex_D1vD0 = results(amex_dds2_D1vs0, contrast = c("dpa","1","0"), independentFiltering=TRUE)
hist(amex_D1vD0$pvalue)
summary(amex_D1vD0, alpha = 0.05)

pearson_microarray<-cor(log2(counts(amex_dds2_D1vs0,norm=TRUE)+1),method="pearson")
diff_pearson<-1-pearson_microarray
dist_pearson_micro<-as.dist(diff_pearson)
pearson_micro_clust<-hclust(dist_pearson_micro,method="ward.D")
plot(pearson_micro_clust)

amex.dds_2v0 <- DESeqDataSetFromMatrix(countData = round(amex.eData[,c(4,6,7,9)]), colData = amex.pData[c(4,6,7,9),], design = ~dpa)
amex.dds_2v0 <- amex.dds_2v0[rowSums(fpm(amex.dds_2v0)>0)>=4]
amex_dds2_D2vs0 <- DESeq(amex.dds_2v0, test=c("Wald"), parallel = TRUE)
amex_D2vD0 = results(amex_dds2_D2vs0, contrast = c("dpa","2","1"), independentFiltering=TRUE)
hist(amex_D2vD0$pvalue)
summary(amex_D2vD0, alpha = 0.05)

pearson_microarray<-cor(log2(counts(amex_dds2_D2vs0,norm=TRUE)+1),method="pearson")
diff_pearson<-1-pearson_microarray
dist_pearson_micro<-as.dist(diff_pearson)
pearson_micro_clust<-hclust(dist_pearson_micro,method="complete")
plot(pearson_micro_clust)

amex.dds_3v0 <- DESeqDataSetFromMatrix(countData = round(amex.eData[,c(7,9,10,12)]), colData = amex.pData[c(7,9,10,12),], design = ~dpa)
amex.dds_3v0 <- amex.dds_3v0[rowSums(fpm(amex.dds_3v0)>0)>=4]
amex_dds2_D3vs0 <- DESeq(amex.dds_3v0, test=c("Wald"), parallel = TRUE)
amex_D3vD0 = results(amex_dds2_D3vs0, contrast = c("dpa","3","2"), independentFiltering=TRUE)
hist(amex_D3vD0$pvalue)
summary(amex_D3vD0, alpha = 0.05)

pearson_microarray<-cor(log2(counts(amex_dds2_D3vs0,norm=TRUE)+1),method="pearson")
diff_pearson<-1-pearson_microarray
dist_pearson_micro<-as.dist(diff_pearson)
pearson_micro_clust<-hclust(dist_pearson_micro,method="ward.D")
plot(pearson_micro_clust)

amex.dds_4v0 <- DESeqDataSetFromMatrix(countData = round(amex.eData[,c(10,12,13:14)]), colData = amex.pData[c(10,12,13:14),], design = ~dpa)
amex.dds_4v0 <- amex.dds_4v0[rowSums(fpm(amex.dds_4v0)>0)>=4]
amex_dds2_D4vs0 <- DESeq(amex.dds_4v0, test=c("Wald"), parallel = TRUE)
amex_D4vD0 = results(amex_dds2_D4vs0, contrast = c("dpa","4","3"), independentFiltering=TRUE)
hist(amex_D4vD0$pvalue)
summary(amex_D4vD0, alpha = 0.05)

pearson_microarray<-cor(log2(counts(amex_dds2_D4vs0,norm=TRUE)+1),method="pearson")
diff_pearson<-1-pearson_microarray
dist_pearson_micro<-as.dist(diff_pearson)
pearson_micro_clust<-hclust(dist_pearson_micro,method="complete")
plot(pearson_micro_clust)

amex.dds_5v0 <- DESeqDataSetFromMatrix(countData = round(amex.eData[,c(13:14,17:18)]), colData = amex.pData[c(13:14,17:18),], design = ~dpa)
amex.dds_5v0 <- amex.dds_5v0[rowSums(fpm(amex.dds_5v0)>0)>=4]
amex_dds2_D5vs0 <- DESeq(amex.dds_5v0, test=c("Wald"), parallel = TRUE)
amex_D5vD0 = results(amex_dds2_D5vs0, contrast = c("dpa","5","4"), independentFiltering=TRUE)
hist(amex_D5vD0$pvalue)
summary(amex_D5vD0, alpha = 0.05)

pearson_microarray<-cor(log2(counts(amex_dds2_D5vs0,norm=TRUE)+1),method="pearson")
diff_pearson<-1-pearson_microarray
dist_pearson_micro<-as.dist(diff_pearson)
pearson_micro_clust<-hclust(dist_pearson_micro,method="complete")
plot(pearson_micro_clust)


###
amex_D1vD0_df = data.frame(amex_D1vD0)
amex_D1vD0_df$txid = rownames(amex_D1vD0_df)
amex_D1vD0_df$geneID = amex_D1vD0_df[which(rownames(amex_D1vD0_df)%in%three_species_annot$AxREFv2),]
amex_D2vD0_df = data.frame(amex_D2vD0)
amex_D2vD0_df$txid = rownames(amex_D2vD0_df)
amex_D2vD0_df$geneID = amex_D2vD0_df[which(rownames(amex_D2vD0_df)%in%three_species_annot$AxREFv2),]
amex_D3vD0_df = data.frame(amex_D3vD0)
amex_D3vD0_df$txid = rownames(amex_D3vD0_df)
amex_D3vD0_df$geneID = amex_D3vD0_df[which(rownames(amex_D3vD0_df)%in%three_species_annot$AxREFv2),]
amex_D4vD0_df = data.frame(amex_D4vD0)
amex_D4vD0_df$txid = rownames(amex_D4vD0_df)
amex_D4vD0_df$geneID = amex_D4vD0_df[which(rownames(amex_D4vD0_df)%in%three_species_annot$AxREFv2),]
amex_D5vD0_df = data.frame(amex_D5vD0)
amex_D5vD0_df$txid = rownames(amex_D5vD0_df)
amex_D5vD0_df$geneID = amex_D5vD0_df[which(rownames(amex_D5vD0_df)%in%three_species_annot$AxREFv2),]
###

amex_independent_tests = union(union(union(union(union(rownames(amex_D1vD0[which(amex_D1vD0$padj < 0.05),]),
                                                       rownames(amex_D2vD0[which(amex_D2vD0$padj < 0.05),])),
                                                 rownames(amex_D3vD0[which(amex_D3vD0$padj < 0.05),])),
                                           rownames(amex_D4vD0[which(amex_D4vD0$padj < 0.05),])),
                                     rownames(amex_D4vD0[which(amex_D4vD0$padj < 0.05),])),
                               rownames(amex_D5vD0[which(amex_D5vD0$padj < 0.05),]))

amex_DEG = data.frame(trans = amex_independent_tests)
rownames(amex_DEG) = amex_DEG[,1]
amex_DEG$D1v0_FC = amex_D1vD0_df[match(rownames(amex_DEG),rownames(amex_D1vD0_df)),c(2)]
amex_DEG$D2v0_FC = amex_D2vD0_df[match(rownames(amex_DEG),rownames(amex_D2vD0_df)),c(2)]
amex_DEG$D3v0_FC = amex_D3vD0_df[match(rownames(amex_DEG),rownames(amex_D3vD0_df)),c(2)]
amex_DEG$D4v0_FC = amex_D4vD0_df[match(rownames(amex_DEG),rownames(amex_D4vD0_df)),c(2)]
amex_DEG$D5v0_FC = amex_D5vD0_df[match(rownames(amex_DEG),rownames(amex_D5vD0_df)),c(2)]
amex_DEG$D1v0_FDR = amex_D1vD0_df[match(rownames(amex_DEG),rownames(amex_D1vD0_df)),c(6)]
amex_DEG$D2v0_FDR = amex_D2vD0_df[match(rownames(amex_DEG),rownames(amex_D2vD0_df)),c(6)]
amex_DEG$D3v0_FDR = amex_D3vD0_df[match(rownames(amex_DEG),rownames(amex_D3vD0_df)),c(6)]
amex_DEG$D4v0_FDR = amex_D4vD0_df[match(rownames(amex_DEG),rownames(amex_D4vD0_df)),c(6)]
amex_DEG$D5v0_FDR = amex_D5vD0_df[match(rownames(amex_DEG),rownames(amex_D5vD0_df)),c(6)]
#amex_DEG$geneID = amex_DEG[which(rownames(amex_DEG)%in%three_species_annot$AxREFv2),]
write.table(amex_DEG, "DEG_list/Amex_DEG_fulllist_adjacent.txt",row.names = TRUE, col.names = TRUE, sep='\t')
pearson_microarray<-cor(fpm(amex_dds2)[rownames(counts(amex_dds2, norm=TRUE))%in%amex_independent_tests,],method="spearman")
#pearson_microarray<-cor(fpm(amex_dds2),method="spearman")
diff_pearson<-1-pearson_microarray
dist_pearson_micro<-as.dist(diff_pearson)
pearson_micro_clust<-hclust(dist_pearson_micro,method="complete")
plot(pearson_micro_clust)

amex_fpm_norm = data.frame(fpm(amex_dds2)[rownames(counts(amex_dds2, norm=TRUE))%in%amex_independent_tests,])
amex_fpm_norm$D0_avg = rowMeans(amex_fpm_norm[,c(1:2)])
amex_fpm_norm$D1_avg = rowMeans(amex_fpm_norm[,c(3:4)])
amex_fpm_norm$D2_avg = rowMeans(amex_fpm_norm[,c(5:6)])
amex_fpm_norm$D3_avg = rowMeans(amex_fpm_norm[,c(7:8)])
amex_fpm_norm$D4_avg = rowMeans(amex_fpm_norm[,c(9:10)])
amex_fpm_norm$D5_avg = rowMeans(amex_fpm_norm[,c(11:12)])

library(pvclust)
y = pvclust(data=amex_fpm_norm[,c(13:18)],method.hclust="complete",method.dist="correlation",nboot=1000, parallel = TRUE)
tiff(filename = "figures/AMEX_TimepointCluster_adjacent_Bootstrap.tiff", width = 6, height = 4, units = "in", res = 300)
plot(y)
pvrect(y, alpha=0.94)
dev.off()

nrow(amex_D1vD0[which(amex_D1vD0$padj < 0.05 & amex_D1vD0$log2FoldChange > 0),])
nrow(amex_D2vD0[which(amex_D2vD0$padj < 0.05 & amex_D2vD0$log2FoldChange > 0),])
nrow(amex_D3vD0[which(amex_D3vD0$padj < 0.05 & amex_D3vD0$log2FoldChange > 0),])
nrow(amex_D4vD0[which(amex_D4vD0$padj < 0.05 & amex_D4vD0$log2FoldChange > 0),])
nrow(amex_D5vD0[which(amex_D5vD0$padj < 0.05 & amex_D5vD0$log2FoldChange > 0),])

nrow(amex_D1vD0[which(amex_D1vD0$padj < 0.05 & amex_D1vD0$log2FoldChange < 0),])
nrow(amex_D2vD0[which(amex_D2vD0$padj < 0.05 & amex_D2vD0$log2FoldChange < 0),])
nrow(amex_D3vD0[which(amex_D3vD0$padj < 0.05 & amex_D3vD0$log2FoldChange < 0),])
nrow(amex_D4vD0[which(amex_D4vD0$padj < 0.05 & amex_D4vD0$log2FoldChange < 0),])
nrow(amex_D5vD0[which(amex_D5vD0$padj < 0.05 & amex_D5vD0$log2FoldChange < 0),])

write.table(amex_D1vD0_df, "amex_D1vD0_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(amex_D2vD0_df, "amex_D2vD1_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(amex_D3vD0_df, "amex_D3vD2_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(amex_D4vD0_df, "amex_D4vD3_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(amex_D5vD0_df, "amex_D5vD4_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')

write.table(amex_D1vD0_df[which(amex_D1vD0_df$padj<0.05),], "amex_D1vD0_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(amex_D2vD0_df[which(amex_D2vD0_df$padj<0.05),], "amex_D2vD1_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(amex_D3vD0_df[which(amex_D3vD0_df$padj<0.05),], "amex_D3vD2_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(amex_D4vD0_df[which(amex_D4vD0_df$padj<0.05),], "amex_D4vD3_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(amex_D5vD0_df[which(amex_D5vD0_df$padj<0.05),], "amex_D5vD4_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')

am1 = rownames(amex_D1vD0_df[which(amex_D1vD0_df$padj<0.05),])
am2 = unique(c(rownames(amex_D2vD0_df[which(amex_D2vD0_df$padj<0.05),]),am1))
am3 = unique(c(rownames(amex_D3vD0_df[which(amex_D3vD0_df$padj<0.05),]),am2))
am4 = unique(c(rownames(amex_D4vD0_df[which(amex_D4vD0_df$padj<0.05),]),am3))
am5 = unique(c(rownames(amex_D5vD0_df[which(amex_D5vD0_df$padj<0.05),]),am4))

##### AAND - run comparisons  #######
and.dds <- DESeqDataSetFromMatrix(countData = round(and.eData[,c(1:3,4,6,7:9,10,12,13,15,16,17)]), colData = and.pData[c(1:3,4,6,7:9,10,12,13,15,16,17),], design = ~dpa)
#and.dds <- and.dds[rowSums(fpm(and.dds)>0)>=15]
and_dds <- DESeq(and.dds, test=c("LRT"), reduced =~1, parallel = TRUE)
#and_LRT = results(and_dds, independentFiltering=TRUE)
#hist(and_LRT$pvalue)
#summary(and_LRT, alpha = 0.05)

pearson_microarray<-cor(log2(counts(and_dds,norm=TRUE)+1),method="spearman")
diff_pearson<-1-pearson_microarray
dist_pearson_micro<-as.dist(diff_pearson)
pearson_micro_clust<-hclust(dist_pearson_micro,method="average")
plot(pearson_micro_clust)

###
and.dds_1v0 <- DESeqDataSetFromMatrix(countData = round(and.eData[,c(1:3,4,6)]), colData = and.pData[c(1:3,4,6),], design = ~dpa)
and.dds_1v0 <- and.dds_1v0[rowSums(fpm(and.dds_1v0)>0)>=5]
and_dds_D1vs0 <- DESeq(and.dds_1v0, test=c("Wald"), parallel = TRUE)
and_D1vD0 = results(and_dds_D1vs0, contrast = c("dpa","1","0"), independentFiltering=TRUE)
hist(and_D1vD0$pvalue)
summary(and_D1vD0, alpha = 0.05)

pearson_microarray<-cor(log2(counts(and_dds_D1vs0,norm=TRUE)+1),method="spearman")
diff_pearson<-1-pearson_microarray
dist_pearson_micro<-as.dist(diff_pearson)
pearson_micro_clust<-hclust(dist_pearson_micro,method="average")
plot(pearson_micro_clust)

and.dds_2v0 <- DESeqDataSetFromMatrix(countData = round(and.eData[,c(4,6,7:9)]), colData = and.pData[c(4,6,7:9),], design = ~dpa)
and.dds_2v0 <- and.dds_2v0[rowSums(fpm(and.dds_2v0)>0)>=5]
and_dds_D2vs0 <- DESeq(and.dds_2v0, test=c("Wald"), parallel = TRUE)
and_D2vD0 = results(and_dds_D2vs0, contrast = c("dpa","2","1"), independentFiltering=TRUE)
hist(and_D2vD0$pvalue)
summary(and_D2vD0, alpha = 0.05)

pearson_microarray<-cor(log2(counts(and_dds_D2vs0,norm=TRUE)+1),method="spearman")
diff_pearson<-1-pearson_microarray
dist_pearson_micro<-as.dist(diff_pearson)
pearson_micro_clust<-hclust(dist_pearson_micro,method="average")
plot(pearson_micro_clust)


and.dds_3v0 <- DESeqDataSetFromMatrix(countData = round(and.eData[,c(7:9,10,12)]), colData = and.pData[c(7:9,10,12),], design = ~dpa)
and.dds_3v0 <- and.dds_3v0[rowSums(fpm(and.dds_3v0)>0)>=5]
and_dds_D3vs0 <- DESeq(and.dds_3v0, test=c("Wald"), parallel = TRUE)
and_D3vD0 = results(and_dds_D3vs0, contrast = c("dpa","3","2"), independentFiltering=TRUE)
hist(and_D3vD0$pvalue)
summary(and_D3vD0, alpha = 0.05)

pearson_microarray<-cor(log2(counts(and_dds_D3vs0,norm=TRUE)+1),method="spearman")
diff_pearson<-1-pearson_microarray
dist_pearson_micro<-as.dist(diff_pearson)
pearson_micro_clust<-hclust(dist_pearson_micro,method="average")
plot(pearson_micro_clust)

and.dds_4v0 <- DESeqDataSetFromMatrix(countData = round(and.eData[,c(10,12,13,15)]), colData = and.pData[c(10,12,13,15),], design = ~dpa)
and.dds_4v0 <- and.dds_4v0[rowSums(fpm(and.dds_4v0)>0)>=4]
and_dds_D4vs0 <- DESeq(and.dds_4v0, test=c("Wald"), parallel = TRUE)
and_D4vD0 = results(and_dds_D4vs0, contrast = c("dpa","4","3"), independentFiltering=TRUE)
hist(and_D4vD0$pvalue)
summary(and_D4vD0, alpha = 0.05)

pearson_microarray<-cor(log2(counts(and_dds_D4vs0,norm=TRUE)+1),method="spearman")
diff_pearson<-1-pearson_microarray
dist_pearson_micro<-as.dist(diff_pearson)
pearson_micro_clust<-hclust(dist_pearson_micro,method="complete")
plot(pearson_micro_clust)

and.dds_5v0 <- DESeqDataSetFromMatrix(countData = round(and.eData[,c(13,15,16:17)]), colData = and.pData[c(13,15,16:17),], design = ~dpa)
and.dds_5v0 <- and.dds_5v0[rowSums(fpm(and.dds_5v0)>0)>=4]
and_dds_D5vs0 <- DESeq(and.dds_5v0, test=c("Wald"), parallel = TRUE)
and_D5vD0 = results(and_dds_D5vs0, contrast = c("dpa","5","4"), independentFiltering=TRUE)
hist(and_D5vD0$pvalue)
summary(and_D5vD0, alpha = 0.05)

pearson_microarray<-cor(log2(counts(and_dds_D5vs0,norm=TRUE)+1),method="pearson")
diff_pearson<-1-pearson_microarray
dist_pearson_micro<-as.dist(diff_pearson)
pearson_micro_clust<-hclust(dist_pearson_micro,method="complete")
plot(pearson_micro_clust)

####

and_independent_tests = union(union(union(union(union(rownames(and_D1vD0[which(and_D1vD0$padj < 0.05),]),
                                                      rownames(and_D2vD0[which(and_D2vD0$padj < 0.05),])),
                                                rownames(and_D3vD0[which(and_D3vD0$padj < 0.05),])),
                                          rownames(and_D4vD0[which(and_D4vD0$padj < 0.05),])),
                                    rownames(and_D4vD0[which(and_D4vD0$padj < 0.05),])),
                              rownames(and_D5vD0[which(and_D5vD0$padj < 0.05),]))


and_D1vD0_df = data.frame(and_D1vD0)
and_D1vD0_df$txid = rownames(and_D1vD0_df)
and_D2vD0_df = data.frame(and_D2vD0)
and_D2vD0_df$txid = rownames(and_D2vD0_df)
and_D3vD0_df = data.frame(and_D3vD0)
and_D3vD0_df$txid = rownames(and_D3vD0_df)
and_D4vD0_df = data.frame(and_D4vD0)
and_D4vD0_df$txid = rownames(and_D4vD0_df)
and_D5vD0_df = data.frame(and_D5vD0)
and_D5vD0_df$txid = rownames(and_D5vD0_df)

and_DEG = data.frame(trans = and_independent_tests)
rownames(and_DEG) = and_DEG[,1]
and_DEG$D1v0_FC = and_D1vD0_df[match(rownames(and_DEG),rownames(and_D1vD0_df)),c(2)]
and_DEG$D2v0_FC = and_D2vD0_df[match(rownames(and_DEG),rownames(and_D2vD0_df)),c(2)]
and_DEG$D3v0_FC = and_D3vD0_df[match(rownames(and_DEG),rownames(and_D3vD0_df)),c(2)]
and_DEG$D4v0_FC = and_D4vD0_df[match(rownames(and_DEG),rownames(and_D4vD0_df)),c(2)]
and_DEG$D5v0_FC = and_D5vD0_df[match(rownames(and_DEG),rownames(and_D5vD0_df)),c(2)]
and_DEG$D1v0_FDR = and_D1vD0_df[match(rownames(and_DEG),rownames(and_D1vD0_df)),c(6)]
and_DEG$D2v0_FDR = and_D2vD0_df[match(rownames(and_DEG),rownames(and_D2vD0_df)),c(6)]
and_DEG$D3v0_FDR = and_D3vD0_df[match(rownames(and_DEG),rownames(and_D3vD0_df)),c(6)]
and_DEG$D4v0_FDR = and_D4vD0_df[match(rownames(and_DEG),rownames(and_D4vD0_df)),c(6)]
and_DEG$D5v0_FDR = and_D5vD0_df[match(rownames(and_DEG),rownames(and_D5vD0_df)),c(6)]

write.table(and_DEG, "DEG_list/and_DEG_fulllist_adjacent.txt",row.names = TRUE, col.names = TRUE, sep='\t')


library(pvclust)
y = pvclust(data=amby_DEGs_Rv0[,c(6:20)],method.hclust="complete",method.dist="correlation",nboot=10000, parallel = TRUE)
tiff(filename = "figures/Ambystoma_DEG_Rv0Comparisons.tiff", width = 6, height = 8, units = "in", res = 300)
plot(y)
pvrect(y, alpha=0.94)
dev.off()


nrow(and_D1vD0[which(and_D1vD0$padj < 0.05 & and_D1vD0$log2FoldChange > 0),])
nrow(and_D2vD0[which(and_D2vD0$padj < 0.05 & and_D2vD0$log2FoldChange > 0),])
nrow(and_D3vD0[which(and_D3vD0$padj < 0.05& and_D3vD0$log2FoldChange > 0),])
nrow(and_D4vD0[which(and_D4vD0$padj < 0.05& and_D4vD0$log2FoldChange > 0),])
nrow(and_D5vD0[which(and_D5vD0$padj < 0.05& and_D5vD0$log2FoldChange > 0),])

nrow(and_D1vD0[which(and_D1vD0$padj < 0.05& and_D1vD0$log2FoldChange < 0),])
nrow(and_D2vD0[which(and_D2vD0$padj < 0.05& and_D2vD0$log2FoldChange < 0),])
nrow(and_D3vD0[which(and_D3vD0$padj < 0.05& and_D3vD0$log2FoldChange < 0),])
nrow(and_D4vD0[which(and_D4vD0$padj < 0.05& and_D4vD0$log2FoldChange < 0),])
nrow(and_D5vD0[which(and_D5vD0$padj < 0.05& and_D5vD0$log2FoldChange < 0),])

write.table(and_D1vD0_df, "and_D1vD0_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(and_D2vD0_df, "and_D2vD1_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(and_D3vD0_df, "and_D3vD2_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(and_D4vD0_df, "and_D4vD3_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(and_D5vD0_df, "and_D5vD4_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')

write.table(and_D1vD0_df[which(and_D1vD0_df$padj<0.05),], "and_D1vD0_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(and_D2vD0_df[which(and_D2vD0_df$padj<0.05),], "and_D2vD1_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(and_D3vD0_df[which(and_D3vD0_df$padj<0.05),], "and_D3vD2_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(and_D4vD0_df[which(and_D4vD0_df$padj<0.05),], "and_D4vD3_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(and_D5vD0_df[which(and_D5vD0_df$padj<0.05),], "and_D5vD4_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')

and1 = rownames(and_D1vD0_df[which(and_D1vD0_df$padj<0.05),])
and2 = unique(c(rownames(and_D2vD0_df[which(and_D2vD0_df$padj<0.05),]),and1))
and3 = unique(c(rownames(and_D3vD0_df[which(and_D3vD0_df$padj<0.05),]),and2))
and4 = unique(c(rownames(and_D4vD0_df[which(and_D4vD0_df$padj<0.05),]),and3))
and5 = unique(c(rownames(and_D5vD0_df[which(and_D5vD0_df$padj<0.05),]),and4))

##### AMAC - generate DESeq set  #######
#mac_LRT = results(mac_dds, independentFiltering=TRUE)
#hist(mac_LRT$pvalue)
#summary(mac_LRT)


### comparison
mac.dds_1v0 <- DESeqDataSetFromMatrix(countData = round(mac.eData[,c(1:3,4:5)]), colData = mac.pData[c(1:3,4:5),], design = ~dpa)
mac.dds_1v0 <- mac.dds_1v0[rowSums(fpm(mac.dds_1v0)>0)>=5]
mac_dds_D1vs0 <- DESeq(mac.dds_1v0, test=c("Wald"), parallel = TRUE)
mac_D1vD0 = results(mac_dds_D1vs0, contrast = c("dpa","1","0"), independentFiltering=TRUE)
hist(mac_D1vD0$pvalue)
summary(mac_D1vD0, alpha = 0.05)

pearson_microarray<-cor(log2(counts(mac_dds_D1vs0,norm=TRUE)+1),method="spearman")
diff_pearson<-1-pearson_microarray
dist_pearson_micro<-as.dist(diff_pearson)
pearson_micro_clust<-hclust(dist_pearson_micro,method="complete")
plot(pearson_micro_clust)

mac.dds_2v0 <- DESeqDataSetFromMatrix(countData = round(mac.eData[,c(4:5,7:8)]), colData = mac.pData[c(4:5,7:8),], design = ~dpa)
mac.dds_2v0 <- mac.dds_2v0[rowSums(fpm(mac.dds_2v0)>0)>=4]
mac_dds_D2vs0 <- DESeq(mac.dds_2v0, test=c("Wald"), parallel = TRUE)
mac_D2vD0 = results(mac_dds_D2vs0, contrast = c("dpa","2","1"), independentFiltering=TRUE)
hist(mac_D2vD0$pvalue)
summary(mac_D2vD0, alpha = 0.05)

pearson_microarray<-cor(log2(counts(mac_dds_D2vs0,norm=TRUE)+1),method="spearman")
diff_pearson<-1-pearson_microarray
dist_pearson_micro<-as.dist(diff_pearson)
pearson_micro_clust<-hclust(dist_pearson_micro,method="complete")
plot(pearson_micro_clust)

mac.dds_3v0 <- DESeqDataSetFromMatrix(countData = round(mac.eData[,c(7:8,10:11)]), colData = mac.pData[c(7:8,10:11),], design = ~dpa)
mac.dds_3v0 <- mac.dds_3v0[rowSums(fpm(mac.dds_3v0)>0)>=4]
mac_dds_D3vs0 <- DESeq(mac.dds_3v0, test=c("Wald"), parallel = TRUE)
mac_D3vD0 = results(mac_dds_D3vs0, contrast = c("dpa","3","2"), independentFiltering=TRUE)
hist(mac_D3vD0$pvalue)
summary(mac_D3vD0, alpha = 0.05)

pearson_microarray<-cor(log2(counts(mac_dds_D3vs0,norm=TRUE)+1),method="spearman")
diff_pearson<-1-pearson_microarray
dist_pearson_micro<-as.dist(diff_pearson)
pearson_micro_clust<-hclust(dist_pearson_micro,method="complete")
plot(pearson_micro_clust)

mac.dds_4v0 <- DESeqDataSetFromMatrix(countData = round(mac.eData[,c(10:11,13:14)]), colData = mac.pData[c(10:11,13:14),], design = ~dpa)
mac.dds_4v0 <- mac.dds_4v0[rowSums(fpm(mac.dds_4v0)>0)>=4]
mac_dds_D4vs0 <- DESeq(mac.dds_4v0, test=c("Wald"), parallel = TRUE)
mac_D4vD0 = results(mac_dds_D4vs0, contrast = c("dpa","4","3"), independentFiltering=TRUE)
hist(mac_D4vD0$pvalue)
summary(mac_D4vD0, alpha = 0.05)

pearson_microarray<-cor(log2(counts(mac_dds_D4vs0,norm=TRUE)+1),method="spearman")
diff_pearson<-1-pearson_microarray
dist_pearson_micro<-as.dist(diff_pearson)
pearson_micro_clust<-hclust(dist_pearson_micro,method="ward.D")
plot(pearson_micro_clust)

mac.dds_5v0 <- DESeqDataSetFromMatrix(countData = round(mac.eData[,c(13:14,16:17)]), colData = mac.pData[c(13:14,16:17),], design = ~dpa)
mac.dds_5v0 <- mac.dds_5v0[rowSums(fpm(mac.dds_5v0)>0)>=4]
mac_dds_D5vs0 <- DESeq(mac.dds_5v0, test=c("Wald"), parallel = TRUE)
mac_D5vD0 = results(mac_dds_D5vs0, contrast = c("dpa","5","4"), independentFiltering=TRUE)
hist(mac_D5vD0$pvalue)
summary(mac_D5vD0, alpha = 0.05)

pearson_microarray<-cor(log2(counts(mac_dds_D5vs0,norm=TRUE)+1),method="spearman")
diff_pearson<-1-pearson_microarray
dist_pearson_micro<-as.dist(diff_pearson)
pearson_micro_clust<-hclust(dist_pearson_micro,method="ward.D")
plot(pearson_micro_clust)


mac_independent_tests = union(union(union(union(union(rownames(mac_D1vD0[which(mac_D1vD0$padj < 0.05),]),
                                                      rownames(mac_D2vD0[which(mac_D2vD0$padj < 0.05),])),
                                                rownames(mac_D3vD0[which(mac_D3vD0$padj < 0.05),])),
                                          rownames(mac_D4vD0[which(mac_D4vD0$padj < 0.05),])),
                                    rownames(mac_D4vD0[which(mac_D4vD0$padj < 0.05),])),
                              rownames(mac_D5vD0[which(mac_D5vD0$padj < 0.05),]))

mac_D1vD0_df = data.frame(mac_D1vD0)
mac_D1vD0_df$txid = rownames(mac_D1vD0_df)
mac_D2vD0_df = data.frame(mac_D2vD0)
mac_D2vD0_df$txid = rownames(mac_D2vD0_df)
mac_D3vD0_df = data.frame(mac_D3vD0)
mac_D3vD0_df$txid = rownames(mac_D3vD0_df)
mac_D4vD0_df = data.frame(mac_D4vD0)
mac_D4vD0_df$txid = rownames(mac_D4vD0_df)
mac_D5vD0_df = data.frame(mac_D5vD0)
mac_D5vD0_df$txid = rownames(mac_D5vD0_df)

mac_DEG = data.frame(trans = mac_independent_tests)
rownames(mac_DEG) = mac_DEG[,1]
mac_DEG$D1v0_FC = mac_D1vD0_df[match(rownames(mac_DEG),rownames(mac_D1vD0_df)),c(2)]
mac_DEG$D2v0_FC = mac_D2vD0_df[match(rownames(mac_DEG),rownames(mac_D2vD0_df)),c(2)]
mac_DEG$D3v0_FC = mac_D3vD0_df[match(rownames(mac_DEG),rownames(mac_D3vD0_df)),c(2)]
mac_DEG$D4v0_FC = mac_D4vD0_df[match(rownames(mac_DEG),rownames(mac_D4vD0_df)),c(2)]
mac_DEG$D5v0_FC = mac_D5vD0_df[match(rownames(mac_DEG),rownames(mac_D5vD0_df)),c(2)]
mac_DEG$D1v0_FDR = mac_D1vD0_df[match(rownames(mac_DEG),rownames(mac_D1vD0_df)),c(6)]
mac_DEG$D2v0_FDR = mac_D2vD0_df[match(rownames(mac_DEG),rownames(mac_D2vD0_df)),c(6)]
mac_DEG$D3v0_FDR = mac_D3vD0_df[match(rownames(mac_DEG),rownames(mac_D3vD0_df)),c(6)]
mac_DEG$D4v0_FDR = mac_D4vD0_df[match(rownames(mac_DEG),rownames(mac_D4vD0_df)),c(6)]
mac_DEG$D5v0_FDR = mac_D5vD0_df[match(rownames(mac_DEG),rownames(mac_D5vD0_df)),c(6)]

write.table(mac_DEG, "DEG_list/mac_DEG_fulllist_adjacent.txt",row.names = TRUE, col.names = TRUE, sep='\t')

pearson_microarray<-cor(fpm(mac_dds)[rownames(counts(mac_dds, norm=TRUE))%in%mac_independent_tests,],method="spearman")
#pearson_microarray<-cor(fpm(amex_dds2),method="spearman")
diff_pearson<-1-pearson_microarray
dist_pearson_micro<-as.dist(diff_pearson)
pearson_micro_clust<-hclust(dist_pearson_micro,method="complete")
plot(pearson_micro_clust)

library(pvclust)
y = pvclust(data=mac_fpm_norm[,c(14:19)],method.hclust="complete",method.dist="correlation",nboot=1000, parallel = TRUE)
tiff(filename = "figures/AMAC_TimepointCluster_adjacent_Bootstrap.tiff", width = 6, height = 4, units = "in", res = 300)
plot(y)
pvrect(y, alpha=0.94)
dev.off()



nrow(mac_D1vD0[which(mac_D1vD0$padj < 0.05 & mac_D1vD0$log2FoldChange > 0),])
nrow(mac_D2vD0[which(mac_D2vD0$padj < 0.05 & mac_D2vD0$log2FoldChange > 0),])
nrow(mac_D3vD0[which(mac_D3vD0$padj < 0.05& mac_D3vD0$log2FoldChange > 0),])
nrow(mac_D4vD0[which(mac_D4vD0$padj < 0.05& mac_D4vD0$log2FoldChange > 0),])
nrow(mac_D5vD0[which(mac_D5vD0$padj < 0.05& mac_D5vD0$log2FoldChange > 0),])

nrow(mac_D1vD0[which(mac_D1vD0$padj < 0.05& mac_D1vD0$log2FoldChange < 0),])
nrow(mac_D2vD0[which(mac_D2vD0$padj < 0.05& mac_D2vD0$log2FoldChange < 0),])
nrow(mac_D3vD0[which(mac_D3vD0$padj < 0.05& mac_D3vD0$log2FoldChange < 0),])
nrow(mac_D4vD0[which(mac_D4vD0$padj < 0.05& mac_D4vD0$log2FoldChange < 0),])
nrow(mac_D5vD0[which(mac_D5vD0$padj < 0.05& mac_D5vD0$log2FoldChange < 0),])

write.table(mac_D1vD0_df, "mac_D1vD0_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(mac_D2vD0_df, "mac_D2vD1_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(mac_D3vD0_df, "mac_D3vD2_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(mac_D4vD0_df, "mac_D4vD3_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(mac_D5vD0_df, "mac_D5vD4_all.txt",row.names = TRUE, col.names = TRUE, sep='\t')

write.table(mac_D1vD0_df[which(mac_D1vD0_df$padj<0.05),], "mac_D1vD0_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(mac_D2vD0_df[which(mac_D2vD0_df$padj<0.05),], "mac_D2vD1_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(mac_D3vD0_df[which(mac_D3vD0_df$padj<0.05),], "mac_D3vD2_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(mac_D4vD0_df[which(mac_D4vD0_df$padj<0.05),], "mac_D4vD3_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')
write.table(mac_D5vD0_df[which(mac_D5vD0_df$padj<0.05),], "mac_D5vD4_DEG.txt",row.names = TRUE, col.names = TRUE, sep='\t')

mac1 = rownames(mac_D1vD0_df[which(mac_D1vD0_df$padj<0.05),])
mac2 = unique(c(rownames(mac_D2vD0_df[which(mac_D2vD0_df$padj<0.05),]),mac1))
mac3 = unique(c(rownames(mac_D3vD0_df[which(mac_D3vD0_df$padj<0.05),]),mac2))
mac4 = unique(c(rownames(mac_D4vD0_df[which(mac_D4vD0_df$padj<0.05),]),mac3))
mac5 = unique(c(rownames(mac_D5vD0_df[which(mac_D5vD0_df$padj<0.05),]),mac4))
