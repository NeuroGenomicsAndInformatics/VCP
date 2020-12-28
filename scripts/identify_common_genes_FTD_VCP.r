setwd("/home/eteleeb/projects/Bruno_FTD")
library(readxl)
libray(gplots)
library(RColorBrewer)

## MAKE DIRECTORIES 
#dir.create('/home/eteleeb/projects/Bruno_FTD/de_results/Figures')

## function to plot the effect 
plot_effect <- function(FTD.res, VCP.res, name, c) {
  
  FTD.res.sig <- FTD.res[!is.na(FTD.res$pvalue) & FTD.res$pvalue < 0.05, ]
  if (c %in% c('OneMonth', 'TwoMonth')) {
    vcp <- as.data.frame(VCP.res[toupper(VCP.res$`Gene Symbol`) %in% toupper(FTD.res.sig$GeneSymbol),  c('Gene Symbol', 'log(2)FC')])
  } else {
    vcp <- as.data.frame(VCP.res[toupper(VCP.res$`Gene Symbol`) %in% toupper(FTD.res.sig$GeneSymbol),  c('Gene Symbol', 'avg_logFC')])
  }
  colnames(vcp) <- c('GeneSymbol', 'log(2)FC')
  
  
  vcp$GeneSymbol <- toupper(vcp$GeneSymbol)
  FTD.res.sig <- merge(FTD.res.sig, vcp, sort =F)
  FTD.res.sig$is.sameDir <- 'Y'
  FTD.res.sig[(FTD.res.sig$log2FoldChange > 0 & FTD.res.sig$`log(2)FC` < 0) |
            (FTD.res.sig$log2FoldChange < 0 & FTD.res.sig$`log(2)FC` > 0) , 'is.sameDir'] <- 'N'
  
  plots <- list()
  for (d in c('N', 'Y')) {
   
    if (d == "N") {
      full.name <- paste0(name, '\nAll overlapped genes')
      r2 <- lm(FTD.res.sig$log2FoldChange ~ FTD.res.sig$`log(2)FC`, df = FTD.res.sig)
      r2 <- signif(summary(r2)$r.squared, digits = 2)
      n <- nrow(FTD.res.sig)    
    } else {
      full.name <- paste0(name, '\nSame direction overlapped genes')
      FTD.res.sig <- FTD.res.sig[FTD.res.sig$is.sameDir == "Y", ]
      r2 <- lm(FTD.res.sig$log2FoldChange ~ FTD.res.sig$`log(2)FC`, df = FTD.res.sig)
      r2 <- signif(summary(r2)$r.squared, digits = 2)
      n <- nrow(FTD.res.sig)
    }
    
    pp = ggplot(FTD.res.sig, aes(log2FoldChange, `log(2)FC`)) + geom_point(color="orange2") + theme_bw()
    pp = pp + labs(x='FTD effect', y = 'VCP effect')
    pp = pp + ggtitle(paste0(full.name, ' (n=', n, ')'))
    pp = pp + geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x, size=0.5)
    pp = pp + theme(axis.text.x=element_text(size=10, color="black"),
                  axis.text.y=element_text(size=10, color="black"), legend.title=element_blank(),
                  axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
                  plot.title = element_text(size = 13, hjust=0.5, color="black", face="bold"),
                  legend.position="none",  
                  panel.border = element_rect(linetype='solid', color='black'))
    pp = pp + annotate("text", x = min(FTD.res.sig$log2FoldChange)+0.3, y=max(FTD.res.sig$`log(2)FC`), label = paste0('R2 = ', r2), fontface="bold", size=3)
    
    plots[[d]] <- pp    
  
  }
  
  pdf(paste0('/home/eteleeb/projects/Bruno_FTD/de_results/Weihlexcitatoryneuron_cluster/Figures/', name, '.pdf'), width = 10, height = 5)
  grid.arrange(grobs=plots, ncol=2, nrwo=1)
  dev.off()
  
}

###################################################

## function to compute the significance of overlap 
compute_sig_ov <- function(genes, FTD.sig.genes.res, all.genes, ov.res) {
  pvals <- list()
  FTD.sig.genes <- rownames(FTD.sig.genes.res)
  for (d in c('Y', 'N')) {
    
    if (d == 'Y') {
      idx <- 'SAME'
      ov <- length(ov.res$GeneSymbol[ov.res$is.sameDir == d])  
    } else {
      idx <- 'ALL'
      ov <- length(intersect(genes, FTD.sig.genes))
    }
    
    ss <- length(all.genes) - (length(genes) - ov) - (length(FTD.sig.genes) - ov) - ov
    test.mat <-matrix(c(ov, length(genes) - ov,
                        length(FTD.sig.genes) - ov,  ss), nrow = 2,
                      dimnames = list(Weihl = c("yes", "no"),
                                      FTD = c("yes", "no")))
    test.pval <- fisher.test(test.mat, alternative = "two.sided")$p.value
    pvals[[idx]] <- signif(test.pval, digits = 4)
  }
  return (pvals)
}
#######

## Function to plot VennDiagram for the overlap 
plot_venn <- function(all.sig, vcp.genes, name) {
  
  ov <- intersect(rownames(all.sig), vcp.genes)
  write.table(ov, file=paste0('/home/eteleeb/projects/Bruno_FTD/de_results/WeihlVCPproteomic/overlapped_genes/', name, '_overlapped_genes.tsv'), quote = F, sep="\t", row.names = F, col.names = F)
              
  p = draw.pairwise.venn (length(rownames(all.sig)), length(vcp.genes), length(ov),
                          scaled = F,
                          category = c(paste0('All-FTD-DEGs\n(n=',length(rownames(all.sig)), ')\n'),
                                       paste0('VCP\n(n=', length(vcp.genes),')\n' )), alpha = rep(0.3, 2),
                          fill = c("skyblue2", 'orange3'), cat.pos = c(0.02,0.02), cat.cex = 1, cex=1.5, lty =0
  )


  pdf(paste0('/home/eteleeb/projects/Bruno_FTD/de_results/WeihlVCPproteomic/Figures/', name, '_overlap.pdf'), width=5, height=5)
  #grid.draw(p, top ="test")
  grid.arrange(gTree(children=p), top=paste0("Comparison: ", name))
  dev.off()
}

## extract data 
res.files <- list.files(path = "Zeran_analysis/analysis", pattern = "*res.csv")
res.files <- res.files[!grepl('_sig', res.files)]

comp = c('C9ORF72vsCO_Sex', 'C9ORF72vsCO_Age', 'FTDsporadicvsCO_Sex', 'FTDsporadicvsCO_Age', 'GRNvsCO_Sex', 'GRNvsCO_Age', 'MAPTvsCO_Sex', 'MAPTvsCO_Age', 'MAPTvsCO_DESeq')

##########################################################################
##### Extract DE information for genes in file WeihlVCPproteomic.xlsx #########
##########################################################################
ov.sig1 <- NULL
for (c in c('OneMonth', 'TwoMonth')) {
  cat('Running for: ', c, '\n')
  if (c == "OneMonth") {
    VCP_proteomic <- read_excel("ChrisData/WeihlVCPproteomic.xlsx", sheet = 1)
  } else {
    VCP_proteomic <- read_excel("ChrisData/WeihlVCPproteomic.xlsx", sheet = 2)
  }
  ## loop through files 
  for (k in 1:length(comp)) {
    cmp = comp[k] 
    file.name = res.files[grepl(cmp, res.files)]
    outFile <- gsub('.csv', '', file.name)
    de.res <- read.csv(paste0('Zeran_analysis/analysis/', file.name), sep=",", stringsAsFactors = F, row.names = 1)
    res <- de.res[toupper(rownames(de.res)) %in% toupper(VCP_proteomic$`Gene Symbol`) ,]
    res$GeneSymbol <- rownames(res)
    res <- res[, c('GeneSymbol','baseMean',	'log2FoldChange',	'lfcSE','stat',	'pvalue',	'padj')]
    ## write DE information for all genes
    res.dir = paste0("de_results/WeihlVCPproteomic/", c)
    dir.create(res.dir)
    write.table(res, file=paste0(res.dir, '/',c,'_', outFile, '.tsv'), sep="\t", quote = F, row.names = F)
    
    ## plot the effect between genes 
    plot_effect (res, VCP_proteomic, paste0(c, '_', cmp), c )
    
    ## extract sig results 
    res.sig <- res[!is.na(res$pvalue) & res$pvalue < 0.05, ]  ## using pvalue 
    #res.sig <- res[!is.na(res$padj) & res$padj < 0.05, ]  ## using FDR 
    
    if (nrow(res.sig) > 0 ) {
      ## check if the effect in the same direction or not
      is.sameDir <- as.data.frame(VCP_proteomic[toupper(VCP_proteomic$`Gene Symbol`) %in% toupper(res.sig$GeneSymbol),  c('Gene Symbol', 'log(2)FC')])
      colnames(is.sameDir) <- c('GeneSymbol', 'log(2)FC')
      is.sameDir$GeneSymbol <- toupper(is.sameDir$GeneSymbol)
      res.sig <- merge(res.sig, is.sameDir, sort =F)
      res.sig$is.sameDir <- 'Y'
      res.sig[(res.sig$log2FoldChange > 0 & res.sig$`log(2)FC` < 0) |
                (res.sig$log2FoldChange < 0 & res.sig$`log(2)FC` > 0) , 'is.sameDir'] <- 'N'
      res.sig$`log(2)FC` <- NULL
      
      ## write DE information for all significant genes
      write.table(res.sig, file=paste0(res.dir, '/',c,'_', outFile, '_sig.tsv'), sep="\t", quote = F, row.names = T)
      
      ## plot the venndiagram 
      plot_venn(de.res[!is.na(de.res$pvalue) & de.res$pvalue < 0.05, ], toupper(VCP_proteomic$`Gene Symbol`), paste0(c, '_', cmp))
  
    } 
    ## check the significane of overlap 
    fisher.p <- compute_sig_ov(toupper(VCP_proteomic$`Gene Symbol`), de.res[!is.na(de.res$pvalue) & de.res$pvalue < 0.05, ], rownames(de.res),  res.sig)
    #fisher.p <- compute_sig_ov(toupper(VCP_proteomic$`Gene Symbol`), de.res[!is.na(de.res$padj) & de.res$padj < 0.05, ], rownames(de.res),  res.sig)
    d <- data.frame(data = c, comparison = cmp, sig_ov_all = fisher.p[['ALL']], sig_ov_same_dir = fisher.p[['SAME']] )
    ov.sig1 <- rbind(ov.sig1, d)
    
  }
  
}

## write the significance of overlap results 
write.table(ov.sig1, file=paste0('de_results/Weihlexcitatoryneuron_cluster/significance_of_overlap_PValue.tsv'), sep="\t", quote = F, row.names = F)


##########################################################################
##### Extract DE information for genes in file WeihlVCPdata.xlsx #########
##########################################################################
ov.sig2 <- NULL
for (c in c('cluster10_1m', 'cluster10_2m', 'cluster25_1m', 'cluster25_2m')) {
  cat('Running for: ', c)
  if (c == "cluster10_1m") {
    cc <- read_excel("WeihlVCPdata.xlsx", sheet = 1)
  } else if (c == "cluster10_2m") {
    cc <- read_excel("WeihlVCPdata.xlsx", sheet = 2)
  } else if (c == "cluster25_1m") {
    cc <- read_excel("WeihlVCPdata.xlsx", sheet = 3)
  } else if (c == "cluster25_2m") {
    cc <- read_excel("WeihlVCPdata.xlsx", sheet = 4)
  }
  
  colnames(cc) <- gsub("...1", "Gene Symbol", colnames(cc))
  
  ## loop through files 
  for (k in 1:length(comp)) {
    cmp = comp[k] 
    file.name = res.files[grepl(cmp, res.files)]
    outFile <- gsub('.csv', '', file.name)
    de.res <- read.csv(paste0('Zeran_analysis/analysis/', file.name), sep=",", stringsAsFactors = F, row.names = 1)
    res <- de.res[toupper(rownames(de.res)) %in% toupper(cc$`Gene Symbol`) ,]
    res$GeneSymbol <- rownames(res)
    res <- res[, c('GeneSymbol','baseMean',	'log2FoldChange',	'lfcSE','stat',	'pvalue',	'padj')]
    ## write DE information for all genes
    #write.table(res, file=paste0('de_results/WeihlVCPdata/', c, '/',c,'_', outFile, '.tsv'), sep="\t", quote = F, row.names = F)
    
    ## plot the effect between genes 
    #plot_effect (res, cc, paste0(c, '_', cmp), c )
    
    ## extract sig results 
    #res.sig <- res[!is.na(res$pvalue) & res$pvalue < 0.05, ] ## using PValue
    res.sig <- res[!is.na(res$padj) & res$padj < 0.05, ] ## using FDR
    
    if (nrow(res.sig) > 0 ) {
      ## check if the effect in the same direction or not
      is.sameDir <- as.data.frame(cc[toupper(cc$`Gene Symbol`) %in% toupper(res.sig$GeneSymbol),  c('Gene Symbol', 'avg_logFC')])
      colnames(is.sameDir) <- c('GeneSymbol', 'log(2)FC')
      is.sameDir$GeneSymbol <- toupper(is.sameDir$GeneSymbol)
      res.sig <- merge(res.sig, is.sameDir, sort =F)
      res.sig$is.sameDir <- 'Y'
      res.sig[(res.sig$log2FoldChange > 0 & res.sig$`log(2)FC` < 0) |
                (res.sig$log2FoldChange < 0 & res.sig$`log(2)FC` > 0) , 'is.sameDir'] <- 'N'
      res.sig$`log(2)FC` <- NULL
      
      ## write DE information for all significant genes
      #write.table(res.sig, file=paste0('de_results/WeihlVCPdata/', c, '/',c, '_', outFile, '_sig.tsv'), sep="\t", quote = F, row.names = T)
      
      ## plot the venndiagram 
      #plot_venn(de.res[!is.na(de.res$pvalue) & de.res$pvalue < 0.05, ], toupper(cc$`Gene Symbol`), paste0(c, '_', cmp))
      
    }     
    
    ## check the significane of overlap 
    #fisher.p <- compute_sig_ov(toupper(cc$`Gene Symbol`), de.res[!is.na(de.res$pvalue) & de.res$pvalue < 0.05, ], rownames(de.res),  res.sig)
    fisher.p <- compute_sig_ov(toupper(cc$`Gene Symbol`), de.res[!is.na(de.res$padj) & de.res$padj < 0.05, ], rownames(de.res),  res.sig)
    d <- data.frame(data = c, comparison = cmp, sig_ov_all = fisher.p[['ALL']], sig_ov_same_dir = fisher.p[['SAME']] )
    ov.sig2 <- rbind(ov.sig2, d)
    
  }
  
}

## write the significance of overlap results 
write.table(ov.sig2, file=paste0('de_results/WeihlVCPdata/significance_of_overlap_adjusted_PValue.tsv'), sep="\t", quote = F, row.names = F)


### comparison between first list and second list 
for (c in c('OneMonth', 'TwoMonth')) {
  for (k in 1:length(comp)) {
    cmp = comp[k]
    f1 = list.files(paste0('de_results/WeihlVCPproteomic/',c))[grepl(cmp, list.files(paste0('de_results/WeihlVCPproteomic/',c)))]
    f1 = f1[grepl('_sig', f1)]
    f2 = list.files(paste0('de_results/Weihlexcitatoryneuron_cluster/',c))[grepl(cmp, list.files(paste0('de_results/Weihlexcitatoryneuron_cluster/',c)))]
    f2 = f2[grepl('_sig', f2)]
    if( length(f1) !=1 | length(f2) !=1) { next }
    
    res1 <- read.table(paste0('de_results/WeihlVCPproteomic/',c,'/', f1), header =T, stringsAsFactors = F, check.names = F)
    res2 <- read.table(paste0('de_results/Weihlexcitatoryneuron_cluster/',c,'/', f2), header =T, stringsAsFactors = F, check.names = F)
    
    ## plot overlap 
    res12.ov <- intersect(res1$GeneSymbol, res2$GeneSymbol) 
    
    p = draw.pairwise.venn (nrow(res1), nrow(res2), length(res12.ov), 
                            scaled = F, 
                            category = c(paste0('VCPproteomic_res1\n(n=',nrow(res1), ')\n'), 
                                         paste0('VCPproteomic_res2\n(n=', nrow(res2),')\n' )), alpha = rep(0.3, 2),
                            fill = c("skyblue2", 'orange3'), cat.pos = c(0.02,0.02), cat.cex = 1, cex=1.5, lty =0
    )
    
    
    pdf(paste0('/home/eteleeb/projects/Bruno_FTD/de_results/Weihlexcitatoryneuron_cluster/comparison/',c,'_res1_res2_overlap_', cmp, '.pdf'), width=5, height=5)
    #grid.draw(p, top ="test")
    grid.arrange(gTree(children=p), top=paste0("VCPproteomic_res1_res2_",cmp))
    dev.off()
    
    ## write results 
    d <- list(specific2res1 = res1$GeneSymbol[!res1$GeneSymbol %in% res12.ov], specific2res2 = res2$GeneSymbol[!res2$GeneSymbol %in% res12.ov], res12.overalp=res12.ov)
    dd <- data.frame(matrix(ncol=0, nrow=max(sapply(d,length))) , stringsAsFactors = FALSE)  ## for gene names results 
    
    for (j in 1:3) {
      col.value = as.data.frame((d[j]))
      diff <- (max(sapply(d,length)) - nrow(col.value))
      diff= as.data.frame(rep('', diff))
      colnames(diff) = colnames(col.value)
      col.value <- rbind(col.value, diff)
      dd <- cbind(dd, col.value)
      
    }
    
    write.table(dd, file = paste0("/home/eteleeb/projects/Bruno_FTD/de_results/Weihlexcitatoryneuron_cluster/comparison/",c,"_res1_res2_overlap_", cmp,".tsv"), sep="\t", row.names = F, quote = F )
    
  }
  
}

####################################################################################
### plot the common genes between VCP one/two month an snRNAs and proteomics 
####################################################################################
source( '/home/liz/RNAseq/MiSeq/Decon/code/load.salmon.R' )
source( '/home/liz/RNAseq/MiSeq/Decon/code/extract.expression.R' )
#"MAPTvsCO_Sex"        "MAPTvsCO_Age"        "MAPTvsCO_DESeq" 
cmp = "C9ORF72vsCO_Sex"
res.type = "TwoMonth"
dir.create(paste0('/home/eteleeb/projects/Bruno_FTD/de_results/exp_heatmaps/', cmp))

prt.res <-  read.table(paste0('de_results/WeihlVCPproteomic/', res.type, '/PValue/', res.type, '_', cmp,'_RIN_DESeq_res_sig.tsv'), header =T, stringsAsFactors = F)
sc.res <-  read.table(paste0('de_results/Weihlexcitatoryneuron_cluster/', res.type, '/PValue/', res.type, '_', cmp,'_RIN_DESeq_res_sig.tsv'), header =T, stringsAsFactors = F)

# plot VennDiagram 
shared.genes <- intersect(prt.res$GeneSymbol, sc.res$GeneSymbol)
p = draw.pairwise.venn (nrow(prt.res), nrow(sc.res), length(shared.genes), 
                        scaled = F, 
                        category = c(paste0('VCP_proteomic\n(n=',nrow(prt.res), ')\n'), 
                                     paste0('VCP_scRNA\n(n=', nrow(sc.res),')\n' )), alpha = rep(0.3, 2),
                        fill = c("skyblue2", 'orange3'), cat.pos = c(0.02,0.02), cat.cex = 1, cex=1.5, lty =0
)

pdf(paste0('/home/eteleeb/projects/Bruno_FTD/de_results/exp_heatmaps/',cmp,'/',res.type, '_shared_genes_VCPproteomic_snRNA_',cmp,'.pdf'), width=5, height=5)
#grid.draw(p, top ="test")
grid.arrange(gTree(children=p), top=paste0("Comparison: ", cmp))
dev.off()


## plot the heatmap
pheno.HiSeq = read.csv('/40/Cruchaga_Data/bulkRNASeq/201703_MendelianVsSporadics/03.-Phenotype/archive/parietal/Pheno_11_04_2018v2.csv',header=T)
pheno.HiSeq$AGE_at_death= as.numeric(pheno.HiSeq$AGE_at_death)
pheno.HiSeq$GENDER = as.factor(pheno.HiSeq$GENDER)
pheno.HiSeq$QC.RIN= as.numeric(pheno.HiSeq$QC.RIN)
pheno.FTD.sporadic = pheno.HiSeq[pheno.HiSeq$condition == "FTD" & pheno.HiSeq$Primary_Group %in% c("FTD_Case") & pheno.HiSeq$brain_region %in% c("parietal","insular"),]
pheno.CO = pheno.HiSeq[pheno.HiSeq$condition == "CO" & pheno.HiSeq$DementedControl == 0 & pheno.HiSeq$brain_region %in% c("parietal","insular"),]
pheno.FTD.carrier = pheno.HiSeq[pheno.HiSeq$condition == "FTD" & pheno.HiSeq$Primary_Group %in% c("FTD_Carrier") & pheno.HiSeq$brain_region %in% c("parietal","insular"),]

# compare FTD and CO age #
pheno.FTD = rbind(pheno.FTD.carrier,pheno.FTD.sporadic)

HiSeq = loadTable(path="/40/Cruchaga_Data/bulkRNASeq/201703_MendelianVsSporadics/02.-ProcessedData/04.-Salmon/Hg19/parietal")
HiSeq.NumReads = HiSeq$tpm
ERCCrow = unique(grep(paste(c("DQ","EF"),collapse="|"),
                      rownames(HiSeq.NumReads), value=TRUE))
HiSeq.gene = HiSeq.NumReads[!is.element(rownames(HiSeq.NumReads),ERCCrow),]
dim(HiSeq.gene)
HiSeq.gene.int =matrix(unlist(HiSeq.gene),nrow(HiSeq.gene),ncol(HiSeq.gene))
colnames(HiSeq.gene.int) = colnames(HiSeq.gene)
geneID = read.csv("/40/AD/Expression/2017_03_MendelianVsSporadics/06.-References/spike-in/gencode_v19_annotation_gene_ID.txt", header = F, sep = "\t")
row.match = match(rownames(HiSeq.gene), geneID$V1)
rownames(HiSeq.gene.int) = geneID$V2[row.match]

## for all genes 
VCP_res <- read_excel("ChrisData/Weihlexcitatoryneuron_cluster.xlsx", sheet = 1)  ### for OneMonth
#VCP_res <- read_excel("ChrisData/Weihlexcitatoryneuron_cluster.xlsx", sheet = 1)  ### for TwoMonth

# GRN N = 5 #
#pheno.HiSeq.select = rbind(pheno.FTD.carrier[pheno.FTD.carrier$Primary_gene.subgroup == "GRN",], pheno.CO)
# MAPT N = 3 with covar 
pheno.HiSeq.select = do.call('rbind', 
                             list(#pheno.FTD.sporadic, 
                                  pheno.FTD.carrier[pheno.FTD.carrier$Primary_gene.subgroup == "C9ORF72",],
                                  pheno.CO, 
                                  pheno.FTD.carrier[pheno.FTD.carrier$Primary_gene.subgroup == "GRN",])
                                  #pheno.FTD.carrier[pheno.FTD.carrier$Primary_gene.subgroup == "MAPT",])
                             )

##########################################################################################################
########################################## Plot heatmaps ##########################################
##########################################################################################################
# extract shared genes
s="WeihlVCPproteomic"
name="VCP Proteomics"
s="Weihlexcitatoryneuron_cluster"
name = "VCP Transcriptomics"
res.type="TwoMonth"
dir.create(paste0('/home/eteleeb/projects/Bruno_FTD/de_results/manuscript_figures/',res.type))

## read and organize the data 
grn.shared.genes <- read.table(paste0('de_results/',s,'/',res.type,'/PValue/',res.type,'_GRNvsCO_Sex_RIN_DESeq_res_sig.tsv'), header =T)
c9.shared.genes <- read.table(paste0('de_results/',s,'/',res.type,'/PValue/',res.type,'_C9ORF72vsCO_Sex_RIN_DESeq_res_sig.tsv'), header =T)
if (res.type =="TwoMonth" & s=="WeihlVCPproteomic") {
  c9.shared.genes <- c9.shared.genes[c9.shared.genes$GeneSymbol !="TSC22D2", ]
} else if (res.type =="TwoMonth" & s=="Weihlexcitatoryneuron_cluster")  {
  c9.shared.genes <- c9.shared.genes[c9.shared.genes$GeneSymbol !="INO80", ]
}
shared.genes <- rbind(c9.shared.genes, grn.shared.genes)
shared.genes$data.set <- 'B'
shared.genes[shared.genes$GeneSymbol %in% grn.shared.genes$GeneSymbol & !shared.genes$GeneSymbol %in% c9.shared.genes$GeneSymbol, 'data.set'] <- 'G'
shared.genes[!shared.genes$GeneSymbol %in% grn.shared.genes$GeneSymbol & shared.genes$GeneSymbol %in% c9.shared.genes$GeneSymbol, 'data.set'] <- 'C'
shared.genes$direction <- 'dn'
shared.genes[shared.genes$log2FoldChange > 0, 'direction'] <- 'up'
#shared.genes <- shared.genes[order(shared.genes$log2FoldChange), c('GeneSymbol', 'data.set', 'log2FoldChange')]
shared.genes <- shared.genes[!duplicated(shared.genes$GeneSymbol), ]

## extract up and down and order the data 
shared.genes.up.g = shared.genes[shared.genes$direction=="up" & shared.genes$data.set=="G",] 
shared.genes.up.b = shared.genes[shared.genes$direction=="up" & shared.genes$data.set=="B",]
shared.genes.up.c = shared.genes[shared.genes$direction=="up" & shared.genes$data.set=="C",]
shared.genes.up.g= shared.genes.up.g[order(shared.genes.up.g$log2FoldChange, decreasing = T), ]
shared.genes.up.b= shared.genes.up.b[order(shared.genes.up.b$log2FoldChange, decreasing = T), ]
shared.genes.up.c= shared.genes.up.c[order(shared.genes.up.c$log2FoldChange, decreasing = T), ]
shared.genes.up = do.call('rbind', list(shared.genes.up.g, shared.genes.up.b, shared.genes.up.c))

shared.genes.dn.g = shared.genes[shared.genes$direction=="dn" & shared.genes$data.set=="G",] 
shared.genes.dn.b = shared.genes[shared.genes$direction=="dn" & shared.genes$data.set=="B",]
shared.genes.dn.c = shared.genes[shared.genes$direction=="dn" & shared.genes$data.set=="C",]
shared.genes.dn.g= shared.genes.dn.g[order(shared.genes.dn.g$log2FoldChange, decreasing = F), ]
shared.genes.dn.b= shared.genes.dn.b[order(shared.genes.dn.b$log2FoldChange, decreasing = F), ]
shared.genes.dn.c= shared.genes.dn.c[order(shared.genes.dn.c$log2FoldChange, decreasing = F), ]
shared.genes.dn = do.call('rbind', list(shared.genes.dn.g, shared.genes.dn.b, shared.genes.dn.c))
# merge all 
shared.genes = rbind(shared.genes.dn, shared.genes.up)

### plot side bar 
shared.genes$GeneSymbol <- with(shared.genes, factor(shared.genes$GeneSymbol, levels=rev(unique(shared.genes$GeneSymbol))))
side.bar <- melt(shared.genes[,c('GeneSymbol', 'data.set')])
colnames(side.bar) <- c('GeneSymbol', 'data.set')
side.bar$data.set <- factor(side.bar$data.set, levels = c('G','B', 'C'))
pdf(paste0('/home/eteleeb/projects/Bruno_FTD/de_results/manuscript_figures/',res.type, '/',s,'_',res.type, '_side_bar.pdf'), width = 3.5, height = 28, useDingbats=F)
p <- ggplot(side.bar, aes(x='', y=GeneSymbol)) + geom_point(size=3.2, aes(color=data.set), stroke=0.3)
p <- p + labs(x='', y='', title='') + theme_bw(base_size=12)
p <- p + scale_x_discrete(position = "top", labels = c('GRN', 'Both', 'C9ORF72'))
p <- p + theme(axis.text.x=element_blank(),
               axis.text.x.top = element_text(vjust = 0.5, face="bold"),
               axis.text.y=element_text(size=5, vjust=0.5),
               axis.ticks=element_blank(), axis.line = element_blank(),
               panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
               panel.background=element_blank(), legend.position = "right")
p +  scale_colour_viridis_d("Detected in:", labels=c('GRN','GRN/C9ORF72', 'C9ORF72'), guide= guide_legend(override.aes = list(size = 5)))
dev.off()

## plot the heatmap 
## extract expression data 
HiSeq.gene.select = as.data.frame(HiSeq.gene.int[,pheno.HiSeq.select$ID_RNAseq])
#HiSeq.gene.select = HiSeq.gene.select[rownames(HiSeq.gene.select) %in% shared.genes$GeneSymbol, ]
HiSeq.gene.select$GeneSymbol <- rownames(HiSeq.gene.select)
rownames(HiSeq.gene.select) <- NULL
HiSeq.gene.select = merge(shared.genes[,c('GeneSymbol', 'data.set')], HiSeq.gene.select, sort =F)
rownames(HiSeq.gene.select) <- HiSeq.gene.select$GeneSymbol
HiSeq.gene.select$GeneSymbol <- NULL
HiSeq.gene.select$data.set <- NULL
sidebar.col <- c(rep("#2c7fb8" ,table(shared.genes$log2FoldChange>0)[[1]]), rep("#C84540", table(shared.genes$log2FoldChange>0)[[2]]))


## prepare the data 
exp <-log2(HiSeq.gene.select+1)

colBreaks = seq(-2, 2, 0.01)
#col.pan <- colorpanel(length(colBreaks)-1, "blue", "white", "red")
#col.pan <- colorRampPalette(c("red","white","darkgreen"))(length(colBreaks)-1)
#col.pan <- rev(brewer.pal(length(colBreaks)-1,"Spectral"))
col.pan = magma(length(colBreaks)-1)
## make the ColSideColors
myColSidebar = c(#rep('#2ca25f', nrow(pheno.FTD.sporadic)), 
                 rep('#FDE725FF', nrow(pheno.FTD.carrier[pheno.FTD.carrier$Primary_gene.subgroup == "C9ORF72",])), 
                 rep('#d95f0e', nrow(pheno.CO)), 
                 rep('#440154FF', nrow(pheno.FTD.carrier[pheno.FTD.carrier$Primary_gene.subgroup == "GRN",])) 
                 #rep('black', nrow(pheno.FTD.carrier[pheno.FTD.carrier$Primary_gene.subgroup == "MAPT",]))
                 )

pdf(paste0('/home/eteleeb/projects/Bruno_FTD/de_results/manuscript_figures/',res.type, '/', s,'_',res.type, '_exp_heatmap.pdf'), width = 8, height = 25)
par(oma=c(0, 0, 0, 5))
par(cex.main=0.75)
heatmap.2(as.matrix(exp), 
          dendrogram = "none", 
          Colv = F, 
          Rowv = F, 
          cexRow=0.7, 
          labCol="",
          labRow = as.expression(lapply(rownames(exp), function(a) bquote(bold(.(a))))), 
          #colRow = ifelse(rownames(exp) %in% grn.shared.genes$GeneSymbol & !rownames(exp) %in% c9.shared.genes$GeneSymbol,"", 
          #                 ifelse(!rownames(exp) %in% grn.shared.genes$GeneSymbol & rownames(exp) %in% c9.shared.genes$GeneSymbol, "#FDE725FF", "#35B779FF")),
          scale = "row", 
          #main = paste0("Expression profiles of genes overlapped between FTD and VCP\nfor GRN and C9ORF72 carriers (n=", nrow(exp), ")\n(", name," - ",res.type,")"),
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(128)),
          #breaks = colBreaks, 
          #key = TRUE, keysize = 2,
          #key.par=list(mar=c(0,5,25,10)), 
          density.info = "none", 
          key.title = NA, 
          key.xlab = "Expression",
          trace = "none",
          ColSideColors = myColSidebar, 
          RowSideColors = sidebar.col,
          #sepwidth=0.3, 
          #sepcolor = 'blue',
          margin=c(15, 10), 
          lhei=c(2,45), lwid=c(2,5)
)
col.pos <- structure(list(x = c(0.260, 0.510, 0.800), y = 0.992), .Names = c("x", "y"))
col.labels <- c('C9ORF72 (n=5)', 'CONTROL (n=16)', 'GRN (n=5)')
text(col.pos$x[1], col.pos$y[1], labels=col.labels[1], srt=0, xpd=TRUE, adj = 0, cex=0.5, font=2)
text(col.pos$x[2], col.pos$y[1], labels=col.labels[2], srt=0, xpd=TRUE, adj = 0, cex=0.5, font=2)
text(col.pos$x[3], col.pos$y[1], labels=col.labels[3], srt=0, xpd=TRUE, adj = 0, cex=0.5, font=2)
#legend (0.02, 0.988, c("C9ORF72","Control", "GRN"), col=c("#43a2ca","#31a354", "red"), pch= c(15,15), cex=0.5, xpd = T, pt.cex = 0.8)
dev.off()


####################################################################
###################### Plot the effect (logFC) #####################  
####################################################################
# extract logFC from VCP data 
if (res.type == "OneMonth") {
  vcp.data <- as.data.frame(read_excel(paste0("ChrisData/", s, ".xlsx"), sheet = 1))
} else {
  vcp.data <- read_excel(paste0("ChrisData/",s ,".xlsx"), sheet = 2)
}
vcp.data <- vcp.data[, c('Gene Symbol', 'log(2)FC')]
colnames(vcp.data) <- c('GeneSymbol', 'VCP.logFC')
vcp.data$GeneSymbol <- toupper(vcp.data$GeneSymbol)
p.shared.genes <- merge(shared.genes, vcp.data, sort = F)

p.shared.genes$data.set <- factor(p.shared.genes$data.set, levels = c('G','B', 'C'))
t.shared.genes$data.set <- factor(t.shared.genes$data.set, levels = c('G','B', 'C'))

name="VCP Proteomics"
s1="WeihlVCPproteomic"
pval.res <- read.table(paste0('de_results/',s1, '/', s1, '_significance_of_overlap_PValue.tsv'), header =T, stringsAsFactors = F)
mytable1 <- cbind(sites=c("C9ORF72","GRN"),pval.res[pval.res$data==res.type & grepl('Sex',pval.res$comparison) & grepl('GRN|C9ORF72', pval.res$comparison),c('sig_ov_all', 'sig_ov_same_dir')])
colnames(mytable1) <- c('Carrier', 'Overlap', 'Directionality')
mytable1$Overlap <- signif(mytable1$Overlap, digits = 2)
mytable1$Directionality <- signif(mytable1$Directionality, digits = 2)
colnames(mytable1)[colnames(mytable1)=="Overlap"] <- paste0("Overlap (", nrow(p.shared.genes),")")
colnames(mytable1)[colnames(mytable1)=="Directionality"] <- paste0("Directionality (", nrow(p.shared.genes[p.shared.genes$is.sameDir=="Y",]),")")
mytable1 <- ggtexttable(mytable1, rows = NULL, theme = ttheme("mBlue", base_size = 7, padding = unit(c(2, 2),"mm")))
mytable1 <- tab_add_title(mytable1, 'Significance of overlap:',  padding = unit(0.5, "line"), face="bold", size = 7.2)

mytitle = paste0(name," - ",res.type," (n=",nrow(p.shared.genes),")") 
pp <- ggplot(p.shared.genes, aes(x=log2FoldChange, VCP.logFC)) + geom_point(size=2.5, aes(alpha=factor(ifelse(p.shared.genes$is.sameDir=="N", 0.3, 1)), color=data.set)) + theme_bw()
pp <- pp + labs(x='Human RNA-Seq', y = 'VCP proteomics')
pp <- pp + ggtitle(mytitle)
pp <- pp + geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x, size=1, linetype="dashed")
pp <- pp + theme(axis.text.x=element_text(size=10, color="black"),
                axis.text.y=element_text(size=10, color="black"), 
                legend.title=element_text(size=12, face="bold"),
                axis.title.x=element_text(size=12, face="bold"), axis.title.y=element_text(size=12, face="bold"),
                plot.title = element_text(size = 12, hjust=0.5, color="black", face="bold"),
                legend.position="right", legend.key.size = unit(1,"line"),  
                panel.border = element_rect(linetype='solid', color='black'))
#pp <- pp + scale_colour_viridis_d(name="Same direction", direction = -1, option = "D", labels= c('No', 'Yes'), guide= guide_legend(override.aes = list(size = 5)))
pp <- pp + scale_color_manual(name="Detected in:", values=c('G'='#440154FF','B'='#35B779FF', 'C'='#FDE725FF'), labels=c('G'='GRN','B'='GRN/C9ORF72', 'C'='C9ORF72'), guide= guide_legend(override.aes = list(size = 5)))
pp <- pp + annotation_custom(ggplotGrob(mytable1), xmin= -5, xmax=max(p.shared.genes$log2FoldChange), ymin=3.6, ymax=max(p.shared.genes$VCP.logFC)+0.9)
pp <- pp + scale_alpha_discrete(name ="Directionality", labels=c('NO', 'YES'), guide= guide_legend(override.aes = list(size = 5)))
#pp <- pp + scale_shape_manual(values=c('G'=19, 'B'=17, 'C'=15), labels = c('G'='GRN', 'B'='Both', 'C'='C9ORF72'), guide= guide_legend(override.aes = list(size = 4, color="gray80")))
#pp <- pp + annotation_custom(tableGrob(mytable1, rows=NULL, theme = ttheme_default(base_size = 7, padding = unit(c(2, 2),"mm"), 
#                                                                                   colhead=list(fg_params=list(col="navyblue", fontface="bold"))
#                                                                     )), xmin=0.2, xmax=max(p.shared.genes$log2FoldChange), ymin=2, ymax=max(p.shared.genes$VCP.logFC)+0.9)
pp

## Nueron
name2 = "VCP Transcriptomics"
s2="Weihlexcitatoryneuron_cluster"
pval.res2 <- read.table(paste0('de_results/',s2,'/significance_of_overlap_PValue.tsv'), header =T, stringsAsFactors = F)
mytable2 <- cbind(sites=c("C9ORF72","GRN"),pval.res2[pval.res2$data==res.type & grepl('Sex',pval.res2$comparison) & grepl('GRN|C9ORF72', pval.res2$comparison),c('sig_ov_all', 'sig_ov_same_dir')])
colnames(mytable2) <- c('Carrier', 'Overlap', 'Directionality')
mytable2$Overlap <- signif(mytable2$Overlap, digits = 2)
mytable2$Directionality <- signif(mytable2$Directionality, digits = 2)
colnames(mytable2)[colnames(mytable2)=="Overlap"] <- paste0("Overlap (", nrow(t.shared.genes),")")
colnames(mytable2)[colnames(mytable2)=="Directionality"] <- paste0("Directionality (", nrow(t.shared.genes[t.shared.genes$is.sameDir=="Y",]),")")
mytable2 <- ggtexttable(mytable2, rows = NULL, theme = ttheme("mBlue", base_size = 7, padding = unit(c(2, 2),"mm"))) 
mytable2 <- tab_add_title(mytable2, 'Significance of overlap:', padding = unit(0.5, "line"), face="bold", size = 7.2)

mytitle2 = paste0(name2," - ",res.type," (n=",nrow(t.shared.genes),")") 
pp2 <- ggplot(t.shared.genes, aes(x=log2FoldChange, VCP.logFC)) + geom_point(size=2.5, aes(alpha=factor(ifelse(t.shared.genes$is.sameDir=="N", 0.3, 1)), color=data.set)) + theme_bw()
pp2 <- pp2 + labs(x='Human RNA-Seq', y = 'VCP transcriptomics')
pp2 <- pp2 + ggtitle(mytitle2)
pp2 <- pp2 + geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x, size=1, linetype="dashed")
pp2 <- pp2 + theme(axis.text.x=element_text(size=10, color="black"),
                 axis.text.y=element_text(size=10, color="black"), 
                 axis.title.x=element_text(size=12, face="bold"), axis.title.y=element_text(size=12, face="bold"),
                 plot.title = element_text(size = 12, hjust=0.5, color="black", face="bold"),
                 legend.position="none", legend.key.size = unit(1,"line"),  
                 panel.border = element_rect(linetype='solid', color='black'))
#pp2 <- pp2 + scale_colour_viridis_d(name="Same direction", direction = -1, option = "D", labels= c('No', 'Yes'), guide= guide_legend(override.aes = list(size = 5)))
#pp2 <- pp2 + scale_shape_manual(values=c('G'=19, 'B'=17, 'C'=15), labels = c('G'='GRN', 'B'='Both', 'C'='C9ORF72'), guide= guide_legend(override.aes = list(size = 4, color="gray80")))
pp2 <- pp2 + scale_color_manual(values=c('G'='#440154FF','B'='#35B779FF', 'C'='#FDE725FF'), labels=c('G'='GRN','B'='GRN/C9ORF72', 'C'='C9ORF72'))
pp2 <- pp2 + annotation_custom(ggplotGrob(mytable2), xmin=-0.5, xmax=0, ymin=0.1, ymax=max(t.shared.genes$VCP.logFC)+1)
#pp2 <- pp2 + annotation_custom(ggplotGrob(mytable2, rows=NULL, theme = ttheme_default(base_size = 7, padding = unit(c(2, 2),"mm"), 
#                                                                                     colhead=list(fg_params=list(col="navyblue", fontface=4L))
#                                                                       )), xmin=-0.5, xmax=0, ymin=0.4, ymax=max(t.shared.genes$VCP.logFC)+1)
pp2

#main.title = textGrob('The correlation of the effect between FTD and VCP data for GRN and C9ORF72 carriers\n', gp=gpar(fontsize=14, fontface="bold"))
pdf(paste0('/home/eteleeb/projects/Bruno_FTD/de_results/manuscript_figures/',res.type, '/', res.type, '_effect_size.pdf'), width = 12, height = 5, useDingbats = F, title="")
#grid.arrange(pp, pp2, ncol=2)
grid.arrange(grobs=list(pp, pp2), ncol=2, widths= c(5,4))  
dev.off()

#pp = pp + annotate("text", x = min(FTD.res.sig$log2FoldChange)+0.3, y=max(FTD.res.sig$`log(2)FC`), label = paste0('R2 = ', r2), fontface="bold", size=3)
#"Detected in:", labels=c('GRN', 'C9PRF72C','Both')


##############################################################################
################### generate supplementary tables ################### 
##############################################################################
grn.samples = as.character(pheno.FTD.carrier[pheno.FTD.carrier$Primary_gene.subgroup == "GRN", "ID_RNAseq"]) 
c9.samples = as.character(pheno.FTD.carrier[pheno.FTD.carrier$Primary_gene.subgroup == "C9ORF72", "ID_RNAseq"]) 
co.samples = as.character(pheno.CO$ID_RNAseq)

# compute expression mean 
exp$GRN.Mean.TPM = rowMeans(exp[, grn.samples])
exp$C9ORF72.Mean.TPM = rowMeans(exp[, c9.samples])
exp$Control.Mean.TPM = rowMeans(exp[, co.samples])

mean.exp <- exp[, c('GRN.Mean.TPM', 'C9ORF72.Mean.TPM', 'Control.Mean.TPM')]
mean.exp$GeneSymbol = rownames(mean.exp)
rownames(mean.exp) = NULL

## merge all 
shared.genes <- merge(shared.genes, mean.exp, sort =F)

## add metadata 
#annot <- read.table('annotation.tsv', header =F, check.names = F,  stringsAsFactors = F)
#annot <- annot[annot$Gene !="ENSG00000278217.1", ]
#colnames(annot) <- c('Gene', 'GeneSymbol', 'Gene.locus','Gene.biotype')

final.shared.genes <- merge(annot, shared.genes, sort =F)
final.shared.genes <- final.shared.genes[, c('Gene', 'GeneSymbol', 'Gene.locus','Gene.biotype','data.set','GRN.Mean.TPM','C9ORF72.Mean.TPM','Control.Mean.TPM')]
colnames(final.shared.genes) <- c('Gene', 'Gene.name', 'Gene.locus','Gene.biotype', 'Detected.In','GRN.Mean.TPM','C9ORF72.Mean.TPM','Control.Mean.TPM')
final.shared.genes$Detected.In[final.shared.genes$Detected.In=="G"] <- "GRN"
final.shared.genes$Detected.In[final.shared.genes$Detected.In=="C"] <- "C9ORF72"
final.shared.genes$Detected.In[final.shared.genes$Detected.In=="B"] <- "GRN/C9ORF72"

## extract logFC , PVAlues 
grn.shared.genes.all <- read.table(paste0('de_results/',s,'/',res.type,'/PValue/',res.type,'_GRNvsCO_Sex_RIN_DESeq_res.tsv'), header =T)
grn.shared.genes.all <- grn.shared.genes.all[, c('GeneSymbol', 'log2FoldChange', 'pvalue')]
colnames(grn.shared.genes.all) <- c('Gene.name', 'GRNvsCO.LogFC', 'GRNvsCO.PValue')
grn.shared.genes.all$GRNvsCO.Direction <- 'Up'
grn.shared.genes.all[grn.shared.genes.all$GRNvsCO.LogFC < 0, 'GRNvsCO.Direction'] <- 'Down'

c9.shared.genes.all <- read.table(paste0('de_results/',s,'/',res.type,'/PValue/',res.type,'_C9ORF72vsCO_Sex_RIN_DESeq_res.tsv'), header =T)
c9.shared.genes.all <- c9.shared.genes.all[, c('GeneSymbol', 'log2FoldChange', 'pvalue')]
colnames(c9.shared.genes.all) <- c('Gene.name', 'C9ORF72vsCO.LogFC', 'C9ORF72vsCO.PValue')
c9.shared.genes.all$C9ORF72vsCO.Direction <- 'Up'
c9.shared.genes.all[c9.shared.genes.all$C9ORF72vsCO.LogFC < 0, 'C9ORF72vsCO.Direction'] <- 'Down'

## merge all 
logFC_PValue_grn_c9 <- merge(grn.shared.genes.all, c9.shared.genes.all)
final.shared.genes <- merge(final.shared.genes, logFC_PValue_grn_c9, sort =F)
final.shared.genes <- final.shared.genes[, c('Gene', 'Gene.name', 'Gene.locus','Gene.biotype', 'Detected.In',
                                             'GRNvsCO.LogFC', 'GRNvsCO.PValue', 'GRNvsCO.Direction', 
                                             'C9ORF72vsCO.LogFC', 'C9ORF72vsCO.PValue', 'C9ORF72vsCO.Direction',
                                             'GRN.Mean.TPM','C9ORF72.Mean.TPM','Control.Mean.TPM')]

write.table(final.shared.genes, file=paste0("de_results/manuscript_figures/Final/supp_table_for_", s, '_', res.type, '.tsv'), sep="\t", quote = F, row.names = F)


