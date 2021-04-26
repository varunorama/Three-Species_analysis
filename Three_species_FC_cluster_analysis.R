############### ANALYZE ADJACENT FOLD CHANGES FIRST ######################
setwd("F:/LabWork/Three-Species/")
library(data.table)
amex_D1vD0 = read.table("amex_D1vD0_all.txt",row.name = 1, header = TRUE, sep='\t')
amex_D2vD1 = read.table("amex_D2vD1_all.txt",row.name = 1, header = TRUE, sep='\t')
amex_D3vD2 = read.table("amex_D3vD2_all.txt",row.name = 1, header = TRUE, sep='\t')
amex_D4vD3 = read.table("amex_D4vD3_all.txt",row.name = 1, header = TRUE, sep='\t')
amex_D5vD4 = read.table("amex_D5vD4_all.txt",row.name = 1, header = TRUE, sep='\t')

and_D1vD0 = read.table("and_D1vD0_all.txt",row.name = 1, header = TRUE, sep='\t')
and_D2vD1 = read.table("and_D2vD1_all.txt",row.name = 1, header = TRUE, sep='\t')
and_D3vD2 = read.table("and_D3vD2_all.txt",row.name = 1, header = TRUE, sep='\t')
and_D4vD3 = read.table("and_D4vD3_all.txt",row.name = 1, header = TRUE, sep='\t')
and_D5vD4 = read.table("and_D5vD4_all.txt",row.name = 1, header = TRUE, sep='\t')

mac_D1vD0 = read.table("mac_D1vD0_all.txt",row.name = 1, header = TRUE, sep='\t')
mac_D2vD1 = read.table("mac_D2vD1_all.txt",row.name = 1, header = TRUE, sep='\t')
mac_D3vD2 = read.table("mac_D3vD2_all.txt",row.name = 1, header = TRUE, sep='\t')
mac_D4vD3 = read.table("mac_D4vD3_all.txt",row.name = 1, header = TRUE, sep='\t')
mac_D5vD4 = read.table("mac_D5vD4_all.txt",row.name = 1, header = TRUE, sep='\t')

amby_DEGs_adjacent = read.table("DEG_genes_union_across_all_species_Adjacent.txt",header = TRUE, sep='\t')
amby_DEGs_adjacent$amex_D1vD0 = amex_D1vD0[match(amby_DEGs_adjacent$AMEX_trans,rownames(amex_D1vD0)),2]
amby_DEGs_adjacent$amex_D2vD1 = amex_D2vD1[match(amby_DEGs_adjacent$AMEX_trans,rownames(amex_D2vD1)),2]
amby_DEGs_adjacent$amex_D3vD2 = amex_D3vD2[match(amby_DEGs_adjacent$AMEX_trans,rownames(amex_D3vD2)),2]
amby_DEGs_adjacent$amex_D4vD3 = amex_D4vD3[match(amby_DEGs_adjacent$AMEX_trans,rownames(amex_D4vD3)),2]
amby_DEGs_adjacent$amex_D5vD4 = amex_D5vD4[match(amby_DEGs_adjacent$AMEX_trans,rownames(amex_D5vD4)),2]

amby_DEGs_adjacent$and_D1vD0 = and_D1vD0[match(amby_DEGs_adjacent$AND_trans,rownames(and_D1vD0)),2]
amby_DEGs_adjacent$and_D2vD1 = and_D2vD1[match(amby_DEGs_adjacent$AND_trans,rownames(and_D2vD1)),2]
amby_DEGs_adjacent$and_D3vD2 = and_D3vD2[match(amby_DEGs_adjacent$AND_trans,rownames(and_D3vD2)),2]
amby_DEGs_adjacent$and_D4vD3 = and_D4vD3[match(amby_DEGs_adjacent$AND_trans,rownames(and_D4vD3)),2]
amby_DEGs_adjacent$and_D5vD4 = and_D5vD4[match(amby_DEGs_adjacent$AND_trans,rownames(and_D5vD4)),2]

amby_DEGs_adjacent$mac_D1vD0 = mac_D1vD0[match(amby_DEGs_adjacent$MAC_trans,rownames(mac_D1vD0)),2]
amby_DEGs_adjacent$mac_D2vD1 = mac_D2vD1[match(amby_DEGs_adjacent$MAC_trans,rownames(mac_D2vD1)),2]
amby_DEGs_adjacent$mac_D3vD2 = mac_D3vD2[match(amby_DEGs_adjacent$MAC_trans,rownames(mac_D3vD2)),2]
amby_DEGs_adjacent$mac_D4vD3 = mac_D4vD3[match(amby_DEGs_adjacent$MAC_trans,rownames(mac_D4vD3)),2]
amby_DEGs_adjacent$mac_D5vD4 = mac_D5vD4[match(amby_DEGs_adjacent$MAC_trans,rownames(mac_D5vD4)),2]

amby_DEGs_adjacent[is.na(amby_DEGs_adjacent)] <- 0

library(pvclust)
y = pvclust(data=amby_DEGs_adjacent[,c(6:20)],method.hclust="complete",method.dist="correlation",nboot=10000, parallel = TRUE)
tiff(filename = "figures/Ambystoma_DEG_adjacentComparisons_all.tiff", width = 6, height = 8, units = "in", res = 300)
plot(y)
pvrect(y, alpha=0.94)
dev.off()

library(corrplot)
M <- cor(amby_DEGs_adjacent[,c(6:20)], method = "spearman")
tiff(filename = "figures/Ambystoma_DEG_adjacentComparisons_corrplot.tiff", width = 8, height = 8, units = "in", res = 300)
corrplot(M,type = "upper",method = "number", tl.cex=1, tl.srt = 45,tl.col = "black", number.cex = .7)
dev.off()
cor.test(amby_DEGs_adjacent$amex_D1vD0, amby_DEGs_adjacent$amex_D2vD1, method = "spearman")
cor.test(amby_DEGs_adjacent$amex_D2vD1, amby_DEGs_adjacent$amex_D3vD2, method = "spearman")
cor.test(amby_DEGs_adjacent$amex_D3vD2, amby_DEGs_adjacent$amex_D4vD3, method = "spearman")
cor.test(amby_DEGs_adjacent$amex_D4vD3, amby_DEGs_adjacent$amex_D5vD4, method = "spearman")

############################### REGENERATION VS 0 ##################################################

amex_D1vD0 = read.table("amex_D1vD0_all.txt",row.name = 1, header = TRUE, sep='\t')
amex_D2vD0 = read.table("amex_D2vD0_all.txt",row.name = 1, header = TRUE, sep='\t')
amex_D3vD0 = read.table("amex_D3vD0_all.txt",row.name = 1, header = TRUE, sep='\t')
amex_D4vD0 = read.table("amex_D4vD0_all.txt",row.name = 1, header = TRUE, sep='\t')
amex_D5vD0 = read.table("amex_D5vD0_all.txt",row.name = 1, header = TRUE, sep='\t')

and_D1vD0 = read.table("and_D1vD0_all.txt",row.name = 1, header = TRUE, sep='\t')
and_D2vD0 = read.table("and_D2vD0_all.txt",row.name = 1, header = TRUE, sep='\t')
and_D3vD0 = read.table("and_D3vD0_all.txt",row.name = 1, header = TRUE, sep='\t')
and_D4vD0 = read.table("and_D4vD0_all.txt",row.name = 1, header = TRUE, sep='\t')
and_D5vD0 = read.table("and_D5vD0_all.txt",row.name = 1, header = TRUE, sep='\t')

mac_D1vD0 = read.table("mac_D1vD0_all.txt",row.name = 1, header = TRUE, sep='\t')
mac_D2vD0 = read.table("mac_D2vD0_all.txt",row.name = 1, header = TRUE, sep='\t')
mac_D3vD0 = read.table("mac_D3vD0_all.txt",row.name = 1, header = TRUE, sep='\t')
mac_D4vD0 = read.table("mac_D4vD0_all.txt",row.name = 1, header = TRUE, sep='\t')
mac_D5vD0 = read.table("mac_D5vD0_all.txt",row.name = 1, header = TRUE, sep='\t')

amby_DEGs_Rv0 = read.table("DEG_genes_union_across_all_species_RegVD0.txt",header = TRUE, sep='\t')
amby_DEGs_Rv0$amex_D1vD0 = amex_D1vD0[match(amby_DEGs_Rv0$AMEX_trans,rownames(amex_D1vD0)),2]
amby_DEGs_Rv0$amex_D2vD0 = amex_D2vD0[match(amby_DEGs_Rv0$AMEX_trans,rownames(amex_D2vD0)),2]
amby_DEGs_Rv0$amex_D3vD0 = amex_D3vD0[match(amby_DEGs_Rv0$AMEX_trans,rownames(amex_D3vD0)),2]
amby_DEGs_Rv0$amex_D4vD0 = amex_D4vD0[match(amby_DEGs_Rv0$AMEX_trans,rownames(amex_D4vD0)),2]
amby_DEGs_Rv0$amex_D5vD0 = amex_D5vD0[match(amby_DEGs_Rv0$AMEX_trans,rownames(amex_D5vD0)),2]

amby_DEGs_Rv0$and_D1vD0 = and_D1vD0[match(amby_DEGs_Rv0$AND_trans,rownames(and_D1vD0)),2]
amby_DEGs_Rv0$and_D2vD0 = and_D2vD0[match(amby_DEGs_Rv0$AND_trans,rownames(and_D2vD0)),2]
amby_DEGs_Rv0$and_D3vD0 = and_D3vD0[match(amby_DEGs_Rv0$AND_trans,rownames(and_D3vD0)),2]
amby_DEGs_Rv0$and_D4vD0 = and_D4vD0[match(amby_DEGs_Rv0$AND_trans,rownames(and_D4vD0)),2]
amby_DEGs_Rv0$and_D5vD0 = and_D5vD0[match(amby_DEGs_Rv0$AND_trans,rownames(and_D5vD0)),2]

amby_DEGs_Rv0$mac_D1vD0 = mac_D1vD0[match(amby_DEGs_Rv0$AMAC_trans,rownames(mac_D1vD0)),2]
amby_DEGs_Rv0$mac_D2vD0 = mac_D2vD0[match(amby_DEGs_Rv0$AMAC_trans,rownames(mac_D2vD0)),2]
amby_DEGs_Rv0$mac_D3vD0 = mac_D3vD0[match(amby_DEGs_Rv0$AMAC_trans,rownames(mac_D3vD0)),2]
amby_DEGs_Rv0$mac_D4vD0 = mac_D4vD0[match(amby_DEGs_Rv0$AMAC_trans,rownames(mac_D4vD0)),2]
amby_DEGs_Rv0$mac_D5vD0 = mac_D5vD0[match(amby_DEGs_Rv0$AMAC_trans,rownames(mac_D5vD0)),2]

amby_DEGs_Rv0[is.na(amby_DEGs_Rv0)] <- 0

library(pvclust)
y = pvclust(data=amby_DEGs_Rv0[,c(6:20)],method.hclust="complete",method.dist="correlation",nboot=10000, parallel = TRUE)
tiff(filename = "figures/Ambystoma_DEG_Rv0Comparisons_all.tiff", width = 6, height = 8, units = "in", res = 300)
plot(y)
pvrect(y, alpha=0.94)
dev.off()

library(corrplot)
M <- cor(amby_DEGs_Rv0[,c(6:20)], method = "spearman")
tiff(filename = "figures/Ambystoma_DEG_Rv0Comparisons_corrplot.tiff", width = 8, height = 8, units = "in", res = 300)
corrplot(M,type = "upper",method = "pie", tl.cex=1, tl.srt = 45,tl.col = "black", number.cex = .7)
dev.off()
