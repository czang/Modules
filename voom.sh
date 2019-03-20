#!/bin/sh

INPUT=$1
TUMORDATA=/data5/home/smei/TCGA_data/RNAseq_raw_count/raw/${INPUT}.new.rpkm.tumor.rds
NORMALDATA=/data5/home/smei/TCGA_data/RNAseq_raw_count/raw/normal/${INPUT}.new.rpkm.normal.rds

R --no-save -q << EOF


library(limma)

Tumor <- readRDS("$TUMORDATA")
Normal <- readRDS("$NORMALDATA")
tumor.count <- length(Tumor[1,])
normal.count <- length(Normal[1,])

data=cbind(Tumor,Normal)

design <- model.matrix(~0+factor(c(rep("0",tumor.count),rep("1",normal.count))))
colnames(design) <- c("Tumor", "Normal")
v <- voom(data,design)
fit <- lmFit(v,design)

contrast.mat <- makeContrasts(Tumor-Normal, levels=design)
fit2 <- contrasts.fit(fit, contrast.mat)
fit2 <- eBayes(fit2)
limma.RES <- topTable(fit2, number=nrow(data))
limma.RES_UP_FC2<-subset(limma.RES,adj.P.Val<=0.05&logFC>=1)
#limma.RES_DOWN_FC2<-subset(limma.RES,adj.P.Val<=0.05&logFC<=(-1))
write.table(limma.RES_UP_FC2,"${INPUT}_up_fc2_fdr0.05.txt",quote=F,row.names=F,sep='\t')


EOF

