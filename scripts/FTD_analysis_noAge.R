library("DESeq2")

setwd("/home/liz/Collab/Bruno_FTD/analysis")

source( '/home/liz/RNAseq/MiSeq/Decon/code/load.salmon.R' )
source( '/home/liz/RNAseq/MiSeq/Decon/code/extract.expression.R' )

pheno.HiSeq = read.csv('/40/Cruchaga_Data/bulkRNASeq/201703_MendelianVsSporadics/03.-Phenotype/archive/parietal/Pheno_11_04_2018v2.csv',header=T)
pheno.HiSeq$AGE_at_death= as.numeric(pheno.HiSeq$AGE_at_death)
pheno.HiSeq$GENDER = as.factor(pheno.HiSeq$GENDER)
pheno.HiSeq$QC.RIN= as.numeric(pheno.HiSeq$QC.RIN)
pheno.FTD.sporadic = pheno.HiSeq[pheno.HiSeq$condition == "FTD" & pheno.HiSeq$Primary_Group %in% c("FTD_Case") & pheno.HiSeq$brain_region %in% c("parietal","insular"),]
pheno.CO = pheno.HiSeq[pheno.HiSeq$condition == "CO" & pheno.HiSeq$DementedControl == 0 & pheno.HiSeq$brain_region %in% c("parietal","insular"),]
pheno.FTD.carrier = pheno.HiSeq[pheno.HiSeq$condition == "FTD" & pheno.HiSeq$Primary_Group %in% c("FTD_Carrier") & pheno.HiSeq$brain_region %in% c("parietal","insular"),]

# compare FTD and CO age #
pheno.FTD = rbind(pheno.FTD.carrier,pheno.FTD.sporadic)
shapiro.test(pheno.FTD$AGE_at_death) # p-value = 0.4869 so normally distributed
shapiro.test(pheno.CO$AGE_at_death) # p-value = 0.747 so normally distributed

t.test(pheno.FTD$AGE_at_death, pheno.CO$AGE_at_death) #  p-value = 0.0001352 age with sig difference

HiSeq = loadTable(path="/40/Cruchaga_Data/bulkRNASeq/201703_MendelianVsSporadics/02.-ProcessedData/04.-Salmon/Hg19/parietal")
HiSeq.NumReads = HiSeq$NumReads
ERCCrow = unique(grep(paste(c("DQ","EF"),collapse="|"), 
                       rownames(HiSeq.NumReads), value=TRUE))
HiSeq.gene = HiSeq.NumReads[!is.element(rownames(HiSeq.NumReads),ERCCrow),]
dim(HiSeq.gene)
HiSeq.gene.int =matrix(as.integer(unlist(HiSeq.gene)),nrow(HiSeq.gene),ncol(HiSeq.gene))
colnames(HiSeq.gene.int) = colnames(HiSeq.gene)
geneID = read.csv("/40/AD/Expression/2017_03_MendelianVsSporadics/06.-References/spike-in/gencode_v19_annotation_gene_ID.txt", header = F, sep = "\t")
row.match = match(rownames(HiSeq.gene), geneID$V1)
rownames(HiSeq.gene.int) = geneID$V2[row.match]


## DESeq Analysis ##
# C9ORF72 N = 5 #
pheno.HiSeq.select = rbind(pheno.FTD.carrier[pheno.FTD.carrier$Primary_gene.subgroup == "C9ORF72",], pheno.CO)
HiSeq.gene.select = HiSeq.gene.int[,pheno.HiSeq.select$ID_RNAseq]

dds.HiSeq <- DESeqDataSetFromMatrix(countData = HiSeq.gene.select, colData = pheno.HiSeq.select, design = ~ condition + QC.RIN + GENDER)
cts <- counts(dds.HiSeq)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds.HiSeq <- estimateSizeFactors(dds.HiSeq, geoMeans=geoMeans)
dds.HiSeq <- DESeq(dds.HiSeq)
res.HiSeq <- results(dds.HiSeq,contrast = c("condition", "FTD","CO"))
res.HiSeq.sig <- res.HiSeq[which(res.HiSeq$padj < 0.05),]
write.csv(res.HiSeq,"C9ORF72vsCO_Sex_RIN_DESeq_res.csv",quote = F)
write.csv(res.HiSeq.sig,"C9ORF72vsCO_Sex_RIN_DESeq_sig_res.csv",quote = F)
rm(pheno.HiSeq.select,HiSeq.gene.select,dds.HiSeq,cts,geoMeans,res.HiSeq,res.HiSeq.sig)

# GRN N = 5 #
pheno.HiSeq.select = rbind(pheno.FTD.carrier[pheno.FTD.carrier$Primary_gene.subgroup == "GRN",], pheno.CO)
HiSeq.gene.select = HiSeq.gene.int[,pheno.HiSeq.select$ID_RNAseq]

dds.HiSeq <- DESeqDataSetFromMatrix(countData = HiSeq.gene.select, colData = pheno.HiSeq.select, design = ~ condition + QC.RIN + GENDER)
cts <- counts(dds.HiSeq)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds.HiSeq <- estimateSizeFactors(dds.HiSeq, geoMeans=geoMeans)
dds.HiSeq <- DESeq(dds.HiSeq)
res.HiSeq <- results(dds.HiSeq,contrast = c("condition", "FTD","CO"))
res.HiSeq.sig <- res.HiSeq[which(res.HiSeq$padj < 0.05),]
write.csv(res.HiSeq,"GRNvsCO_Sex_RIN_DESeq_res.csv",quote = F)
write.csv(res.HiSeq.sig,"GRNvsCO_Sex_RIN_DESeq_sig_res.csv",quote = F)
rm(pheno.HiSeq.select,HiSeq.gene.select,dds.HiSeq,cts,geoMeans,res.HiSeq,res.HiSeq.sig)


# MAPT N = 3 becasue the sample size is too small; try both add covar and without any covar #
# with covar #
pheno.HiSeq.select = rbind(pheno.FTD.carrier[pheno.FTD.carrier$Primary_gene.subgroup == "MAPT",], pheno.CO)
HiSeq.gene.select = HiSeq.gene.int[,pheno.HiSeq.select$ID_RNAseq]

dds.HiSeq <- DESeqDataSetFromMatrix(countData = HiSeq.gene.select, colData = pheno.HiSeq.select, design = ~ condition + QC.RIN + GENDER)
cts <- counts(dds.HiSeq)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds.HiSeq <- estimateSizeFactors(dds.HiSeq, geoMeans=geoMeans)
dds.HiSeq <- DESeq(dds.HiSeq)
res.HiSeq <- results(dds.HiSeq,contrast = c("condition", "FTD","CO"))
res.HiSeq.sig <- res.HiSeq[which(res.HiSeq$padj < 0.05),]
write.csv(res.HiSeq,"MAPTvsCO_Sex_RIN_DESeq_res.csv",quote = F)
write.csv(res.HiSeq.sig,"MAPTvsCO_Sex_RIN_DESeq_sig_res.csv",quote = F)
rm(pheno.HiSeq.select,HiSeq.gene.select,dds.HiSeq,cts,geoMeans,res.HiSeq,res.HiSeq.sig)

# Sporadic N = 8 #
pheno.HiSeq.select = rbind(pheno.FTD.sporadic, pheno.CO)
HiSeq.gene.select = HiSeq.gene.int[,pheno.HiSeq.select$ID_RNAseq]

dds.HiSeq <- DESeqDataSetFromMatrix(countData = HiSeq.gene.select, colData = pheno.HiSeq.select, design = ~ condition + QC.RIN + GENDER)
cts <- counts(dds.HiSeq)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds.HiSeq <- estimateSizeFactors(dds.HiSeq, geoMeans=geoMeans)
dds.HiSeq <- DESeq(dds.HiSeq)
res.HiSeq <- results(dds.HiSeq,contrast = c("condition", "FTD","CO"))
res.HiSeq.sig <- res.HiSeq[which(res.HiSeq$padj < 0.05),]
write.csv(res.HiSeq,"FTDsporadicvsCO_Sex_RIN_DESeq_res.csv",quote = F)
write.csv(res.HiSeq.sig,"FTDsporadicvsCO_Sex_RIN_DESeq_sig_res.csv",quote = F)
rm(pheno.HiSeq.select,HiSeq.gene.select,dds.HiSeq,cts,geoMeans,res.HiSeq,res.HiSeq.sig)
