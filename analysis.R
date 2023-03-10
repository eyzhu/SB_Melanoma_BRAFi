# install devtools if necessary
library(DESeq2)
library(fgsea)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(limma)
library(topGO)
library(ggVennDiagram)
library(ggvenn)
library(ggforce)
library(reshape2)
library(fgsea)
library(ggpubr)
library("EnsDb.Hsapiens.v86")
library("ensembldb")
library("tximport")
library("gtools")

#### Load VAV1 data ########
# Generate transcript abundance files using kallisto, version 0.46.0 
# cmd: "kallisto quant -i GRCH38.86.idx -t 3 -o output fastq1 fastq2

# Annotate transcripts
# edb = EnsDb.Hsapiens.v86
# df = transcripts(edb, columns=c("tx_id","gene_name"), return.type="data.frame")
# tx2gene = df
# 
# kallFs = file.path(Samples,"output","abundance.tsv"), samples is the path to the samples below
# 
# names(kallFs) = c("A375_E_VEH_1", "A375_E_VEH_2", "A375_E_VEH_3",
#                   "A375_VAV1_VEH_1", "A375_VAV1_VEH_2", "A375_VAV1_VEH_3",
#                   "A375_VAV1_SRCi_1", "A375_VAV1_SRCi_2", "A375_VAV1_SRCi_3"
# )
#
# txi <- tximport(kallFs, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion=TRUE)
#############################

sampleLable = colnames(txi$abundance)
subtype = sub("_[1-9]$", "", sampleLable)

sampInfo = subtype
sampInfo = factor(sampInfo)

sampInfo = data.frame(group=sampInfo)
rownames(sampInfo) = sampleLable

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = sampInfo,
                                   design = ~ group)

dds <- estimateSizeFactors(ddsTxi)
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
dds <- dds[idx,]
dds <- DESeq(dds)

ntd = normTransform(dds)
normExpression =  assay(ntd)

colname = colnames(normExpression)[1]

A375_VAV1 = data.frame(results(dds, contrast=c("group","A375_VAV1_VEH", "A375_E_VEH")))
A375_VAV1_SRCi = data.frame(results(dds, contrast=c("group","A375_VAV1_SRCi", "A375_VAV1_VEH")))
A375_VAV1_sig = A375_VAV1[A375_VAV1[,"padj"]<0.01,]

A375_VAV1_SRCi_sig = A375_VAV1_SRCi[A375_VAV1_SRCi[,"padj"]<0.01,]
A375_VAV1_SRCi_sig = A375_VAV1_SRCi[!is.na(A375_VAV1_SRCi[,2]),]

orderedVAV1 =  A375_VAV1_sig[order(-A375_VAV1_sig[,2]), ]
v1Names = rownames(head(orderedVAV1[-1,],100))

write.csv(A375_VAV1, file="a375_vav1.csv")
write.csv(A375_VAV1_SRCi, file="a375_vav1_src.csv")

normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}

scaled = t(as.data.frame(apply(normExpression_R1KD, 1, normalize)))

## violin plot
inv = read.csv("data/Invasive_genes2.csv", stringsAsFactors = FALSE, header=FALSE)
inv2 = unlist(inv)
inv2 = inv2[inv2 %in% rownames(scaled)]

a375_expr = melt(scaled[inv2, ])
a375_expr = cbind(a375_expr, rep=rep(1, nrow(a375_expr)))
a375_expr[grep("VAV1_VEH", a375_expr[,2]), 4] = 2 
a375_expr[grep("VAV1_SRCi", a375_expr[,2]), 4] = 3
a375_expr[,4] = factor(a375_expr[,4])

a375_expr = cbind(a375_expr, rep=rep(1, nrow(a375_expr)))
a375_expr[grep("_2", a375_expr[,2]), 5] = 2
a375_expr[grep("_3", a375_expr[,2]), 5] = 3
a375_expr[,5] = factor(a375_expr[,5])

colnames(a375_expr) = c("gene", "expName", "expression", "group", "rep")

## Figure 4E
p <- ggplot(a375_expr, aes(fill=rep, x=group, y=expression)) + geom_violin() +
  geom_sina(alpha=0.5) + theme_classic() + theme(text = element_text(size = 20)) +
  scale_fill_manual(values=c("aliceblue", "aliceblue", "aliceblue"))
###

a375_ranks = A375_VAV1[,"log2FoldChange"]
names(a375_ranks) = rownames(A375_VAV1)
a375_ranks = a375_ranks[order(a375_ranks)]

invList = list(invasive=unlist(inv))

fgseaRes <- fgsea(pathways = invList, stats = a375_ranks)

## Figure 4A subpanel
pdf(file="vav1_enrich.pdf", width=3.5, height=3)
  plotEnrichment(invList$invasive, a375_ranks, gseaParam = 1, ticksSize = 0.3) +  
    theme(axis.title = element_text(size = 15), 
          axis.text = element_text(size = 15), legend.text = element_text(size = 20), 
          legend.title = element_text(size = 25))
dev.off()

######### Load Rac1_P29S data ########
# Annotate transcripts
# edb = EnsDb.Hsapiens.v86
# df = transcripts(edb, columns=c("tx_id","gene_name"), return.type="data.frame")
# tx2gene = df
# 
# kallFs = file.path(Samples,"output","abundance.tsv"), samples is the path to the samples below
# 
# names(kallFs) = c("A375_EM_1", "A375_EM_2", "A375_EM_3",
#                   "A375_P29S_1", "A375_P29S_2", "A375_P29S_3"
# )
#
# txi <- tximport(kallFs, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion=TRUE)
#######################################

SampleLable = colnames(txi$abundance)
subtype = sub("_[1-9]$", "", sampleLable)

sampInfo = subtype
sampInfo = factor(sampInfo)

sampInfo = data.frame(group=sampInfo)
rownames(sampInfo) = sampleLable

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = sampInfo,
                                   design = ~ group)

dds <- estimateSizeFactors(ddsTxi)
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
dds <- dds[idx,]
dds <- DESeq(dds)

ntd = normTransform(dds)
normExpressionP29S =  assay(ntd)
# normExpressionP29SSub = normExpressionP29S[,grep("A375", colnames(normExpressionP29S))]
# scaledP29S = t(as.data.frame(apply(normExpressionP29SSub, 1, normalize)))

A375_R1_OE = data.frame(results(dds, contrast=c("group","A375_P29S", "A375_EM")))
sigA375_P29S = A375_R1_OE
sigA375 = sigA375_P29S[sigA375_P29S[,"padj"]<0.01,]
sigA375 = sigA375[!is.na(sigA375[,2]),]
inv = read.csv("invasive_genes2.csv", stringsAsFactors = FALSE, header=FALSE)
invList = list(invasive=unlist(inv))

a375_ranks = A375_R1_OE[,"log2FoldChange"]
names(a375_ranks) = rownames(A375_R1_OE)
a375_ranks = a375_ranks[order(a375_ranks)]

fgseaRes <- fgsea(pathways = invList,
                  stats    = a375_ranks)

# Figure 4B embedded
pdf(file="P29S_enrich.pdf", width=3.5, height=3)
  plotEnrichment(invList$invasive, a375_ranks, gseaParam = 1, ticksSize = 0.3) +  
  theme(axis.title = element_text(size = 15), 
        axis.text = element_text(size = 15), legend.text = element_text(size = 20), 
        legend.title = element_text(size = 25))
dev.off()

# Figure 4D analysis
a375_ranks = A375_VAV1_SRCi[,"log2FoldChange"]
names(a375_ranks) = rownames(A375_VAV1_SRCi)
a375_ranks = a375_ranks[order(a375_ranks)]

vav1U = vav1U[vav1U!="VAV1"]

v1Ulist = list(vav1=vav1U)

fgseaRes <- fgsea(pathways = v1Ulist,
                  stats    = a375_ranks)
# Figure 4D
pdf(file="VAV1_Srci_enrich.pdf", width=3.5, height=3)
plotEnrichment(v1Ulist$vav1, a375_ranks, gseaParam = 1, ticksSize = 0.3) +  
  theme(axis.title = element_text(size = 15), 
        axis.text = element_text(size = 15), legend.text = element_text(size = 20), 
        legend.title = element_text(size = 25))
dev.off()

# Figure 4F
qq = read.csv("data/sigs_vav1.csv", header=TRUE, stringsAsFactors = FALSE)
qq2 = data.frame(unlist(qq[1:5,1]), unlist(qq[1:5,4]))

qq2[,2] = -log10(qq2[,2])
colnames(qq2) = c("x", "y")

# Figure 4F top left
pdf("sigsVav1.pdf", 9,5)

ggplot(qq2, aes(x = reorder(x, y), y = y)) +
  geom_segment(aes(x = reorder(x, y),
                   xend = reorder(x, y),
                   y = 0, yend = y),
               color = "black", lwd = 3) +
  geom_point(size = 8, pch = 21, bg = 4, col = 1) +
  xlab("Group") +
  ylab("") +
  coord_flip() + theme(text = element_text(size = 45), axis.title.x = element_blank(),
                       axis.title.y = element_blank()) + ylim(0,15)

dev.off()

qq = read.csv("data/sigs_p29s.csv", header=TRUE, stringsAsFactors = FALSE)
qq2 = data.frame(unlist(qq[1:5,1]), unlist(qq[1:5,4]))

qq2[,2] = -log10(qq2[,2])
colnames(qq2) = c("x", "y")

# Figure 4F top right
pdf("sigsP29S.pdf", 7,5)

ggplot(qq2, aes(x = reorder(x, y), y = y)) +
  geom_segment(aes(x = reorder(x, y),
                   xend = reorder(x, y),
                   y = 0, yend = y),
               color = "black", lwd = 3) +
  geom_point(size = 8, pch = 21, bg = 4, col = 1) +
  xlab("Group") +
  ylab("") +
  coord_flip() + theme(text = element_text(size = 45), axis.title.x = element_blank(),
                       axis.title.y = element_blank()) + ylim(0,15)

dev.off()

# Figure 4F bottom left
qq = read.csv("data/sigs_undiff.csv", header=TRUE, stringsAsFactors = FALSE)
qq2 = data.frame(unlist(qq[1:5,1]), unlist(qq[1:5,4]))

qq2[,2] = -log10(qq2[,2])
colnames(qq2) = c("x", "y")

pdf("sigsUndiff.pdf", 7.5,5)

ggplot(qq2, aes(x = reorder(x, y), y = y)) +
  geom_segment(aes(x = reorder(x, y),
                   xend = reorder(x, y),
                   y = 0, yend = y),
               color = "black", lwd = 3) +
  geom_point(size = 8, pch = 21, bg = 4, col = 1) +
  xlab("Group") +
  ylab("") +
  coord_flip() + theme(text = element_text(size = 45), axis.title.x = element_blank(),
                       axis.title.y = element_blank()) + ylim(0,40)

dev.off()

# Figure 4F bottom right
qq = read.csv("data/sigs_teads.csv", header=TRUE, stringsAsFactors = FALSE)
qq2 = data.frame(unlist(qq[1:5,1]), unlist(qq[1:5,4]))

qq2[,2] = -log10(qq2[,2])
colnames(qq2) = c("x", "y")

pdf("sigsTeads.pdf", 7.9,5)

ggplot(qq2, aes(x = reorder(x, y), y = y)) +
  geom_segment(aes(x = reorder(x, y),
                   xend = reorder(x, y),
                   y = 0, yend = y),
               color = "black", lwd = 3) +
  geom_point(size = 8, pch = 21, bg = 4, col = 1) +
  xlab("Group") +
  ylab("") +
  coord_flip() + theme(text = element_text(size = 45), axis.title.x = element_blank(),
                       axis.title.y = element_blank()) + ylim(0,40)

dev.off()

# Figure 4B
pathways.hallmark <- gmtPathways("h.all.v6.2.symbols.gmt.txt")

a375_ranks = A375_VAV1[,"log2FoldChange"]
names(a375_ranks) = rownames(A375_VAV1)
# a375_ranks = a375_ranks[abs(a375_ranks)>0.5]

fgseaRes <- fgsea(pathways = pathways.hallmark,
                  stats = a375_ranks)

fgseaResFiltVAV1 = fgseaRes[which(abs(fgseaRes$padj) <0.1), ]
# fgseaResFilt = data.frame(pathway=fgseaResFilt$pathway, Padj=fgseaResFilt$padj, NES=fgseaResFilt$NES)

pdf("hallmark_pathways_V1.pdf", width=6, height=3)
  ggplot(fgseaResFilt, aes(reorder(pathway, -NES), NES)) + geom_col() +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score") +
    theme_minimal()
dev.off()

#### 
a375_ranks = A375_R1_OE[,"log2FoldChange"]
names(a375_ranks) = rownames(A375_R1_OE)
# a375_ranks = a375_ranks[abs(a375_ranks)>0.5]
fgseaRes <- fgsea(pathways = pathways.hallmark,
                  stats = a375_ranks)

fgseaResFiltR1 = fgseaRes[which(abs(fgseaRes$padj) <0.1), ]
# fgseaResFilt = data.frame(pathway=fgseaResFilt$pathway, Padj=fgseaResFilt$padj, NES=fgseaResFilt$NES)

pdf("hallmark_pathways_P29S.pdf", width=6, height=3)
  ggplot(fgseaResFiltR1, aes(reorder(pathway, -NES), NES)) + geom_col() +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score") +
    theme_minimal()
dev.off()

### combined VAV1/Rac1 enrichment
fgseaResFiltV1 = data.frame(fgseaResFiltVAV1, condition=rep("VAV1", nrow(fgseaResFiltVAV1)))
fgseaResFiltR1 = data.frame(fgseaResFiltR1, condition=rep("Rac1P29S", nrow(fgseaResFiltR1)))

fgseaFinal = rbind(fgseaResFiltV1, fgseaResFiltR1)

fgseaFinal[,"condition"] = as.factor(fgseaFinal[,"condition"])
fgseaFinal = fgseaFinal[, c("pathway", "NES", "condition")]

allpaths = unique(fgseaFinal[,"pathway"])

vav1Paths = fgseaFinal[fgseaFinal[,"condition"]=="VAV1",1]
rac1Paths = fgseaFinal[fgseaFinal[,"condition"]=="Rac1P29S",1]

addVav1 = allpaths[!(allpaths %in% vav1Paths)]
addRac1 = allpaths[!(allpaths %in% rac1Paths)]

addVav1Data = data.frame(pathway=addVav1, NES=rep(0.02, length(addVav1)), condition=rep("VAV1", length(addVav1)))
addRac1Data = data.frame(pathway=addRac1, NES=rep(0.02, length(addRac1)), condition=rep("Rac1P29S", length(addRac1)))

fgseaFinal2 = rbind(fgseaFinal, addVav1Data, addRac1Data)
fgseaFinal2[,3] = factor(fgseaFinal2[,3],levels=c("VAV1","Rac1P29S"))

# Figure 4F
pdf("hallmark_pathways_Vav1andP29S.pdf", width=8, height=3)
  ggplot(fgseaFinal2, aes( x=reorder(pathway, -NES), y=NES, fill=condition)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score") + scale_fill_manual(values=c("#56B4E9","#E69F00")) +
    theme_minimal()
dev.off()

## Figure 4A left volcano for VAV1
qq = data.frame(log2FoldChange=as.numeric(A375_VAV1[,2]), pvalue=as.numeric(A375_VAV1[,6]))
rownames(qq) = rownames(A375_VAV1)

keyvals.shape = rep(1, nrow(qq))
keyvals.shape[rownames(qq) %in% inv2] = 17
names(keyvals.shape)[keyvals.shape==1] = "zz"
names(keyvals.shape)[keyvals.shape==17] = "Tsoi Un-Differentiated"

pdf(file="vav1_OE.pdf", width=13, height=10)
  EnhancedVolcano(qq,
                lab = rownames(qq),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'VAV1 Overexpression',
                FCcutoff = 2,
                pointSize = 3.0,
                labSize = 4.0,
                col=c('grey', 'grey', 'grey', "dodgerblue3"),
                colCustom = NULL,
                shapeCustom = keyvals.shape,
                pCutoff = 1e-10,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                arrowheads = FALSE)
dev.off()

## ## Figure 4A (right) volcano for Rac1P29S 
qq = data.frame(log2FoldChange=as.numeric(A375_R1_OE[,2]), pvalue=as.numeric(A375_R1_OE[,6]))
rownames(qq) = rownames(A375_R1_OE)

keyvals.shape = rep(1, nrow(qq))
keyvals.shape[rownames(qq) %in% inv2] = 17
names(keyvals.shape)[keyvals.shape==1] = "zz"
names(keyvals.shape)[keyvals.shape==17] = "Tsoi Un-Differentiated"

## P29S OE volcano
pdf(file="P29S_OE.pdf", width=13, height=10)
EnhancedVolcano(qq,
                lab = rownames(qq),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'VAV1 Overexpression',
                FCcutoff = 2,
                pointSize = 3.0,
                labSize = 4.0,
                col=c('grey', 'grey', 'grey', "#E69F00"),
                colCustom = NULL,
                shapeCustom = keyvals.shape,
                pCutoff = 1e-10,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                arrowheads = FALSE)
dev.off()

###### Figure 4C vav1/rac1 overlap
sigP29S2 = sigA375_P29S[abs(sigA375_P29S[,"log2FoldChange"]) > 2,]
sigVAV1 = A375_VAV1_sig[ abs(A375_VAV1_sig[,"log2FoldChange"]) > 2, ]

qq = list(RAC1P29S=rownames(sigP29S2), VAV1=rownames(sigVAV1), Tsoi=inv2)

pdf(file="venn.pdf", width = 7, height = 4)
  ggvenn(qq, columns = c("VAV1", "RAC1P29S"), fill_color = c("#56B4E9", "#EFC000FF"), auto_scale = TRUE )
dev.off()

###### Figure 4G
load("data/tcga_expr.Rdata")
load("data/clust_assign.Rdata")

wantGenes = c("CSF1R", "PDGFRB", "ABCB1", "ABCG2")

InvWant = normExpressionTCGA_ALL[wantGenes, names(clusts2[clusts2=="Invasive"]) ]
ProfWant = normExpressionTCGA_ALL[wantGenes, names(clusts2[clusts2=="Proliferative"]) ]

ii = melt(InvWant)
ii = data.frame(ii, group="invasive")
pp = melt(ProfWant)
pp = data.frame(pp, group="proliferative")

data = rbind(ii,pp)
data[,4] = factor(data[,4])

pdf("TCGA_cands.pdf", 10,5)
ggplot(data, aes(fill=group, y=log2(value), x=Var1))  +
  geom_violin(width=0.9, position=position_dodge(0.85)) +
  geom_boxplot(width=0.2, outlier.shape = NA, position=position_dodge(0.85)) + scale_fill_brewer(palette="RdBu") + theme_minimal() +
  theme(text = element_text(size = 20)) + stat_compare_means(aes(group = group), method = "t.test")
dev.off()

##### Figure S4
load("data/tcga_expr.Rdata")
load("data/clust_assign.Rdata")

wantGenes = c("CSF1", "IL34", "PDGFA", "PDGFB", "PDGFC", "PDGFD")

InvWant = normExpressionTCGA_ALL[wantGenes, names(clusts2[clusts2=="Invasive"]) ]
ProfWant = normExpressionTCGA_ALL[wantGenes, names(clusts2[clusts2=="Proliferative"]) ]

ii = melt(InvWant)
ii = data.frame(ii, group="invasive")
pp = melt(ProfWant)
pp = data.frame(pp, group="proliferative")

data = rbind(ii,pp)
data[,4] = factor(data[,4])

pdf("TCGA_ligands.pdf", 10,5)
ggplot(data, aes(fill=group, y=log2(value), x=Var1))  +
  geom_violin(width=0.9, position=position_dodge(0.85)) +
  geom_boxplot(width=0.2, outlier.shape = NA, position=position_dodge(0.85)) + scale_fill_brewer(palette="RdBu") + theme_minimal() +
  theme(text = element_text(size = 20)) + stat_compare_means(aes(group = group), method = "t.test")
dev.off()

########## Fig 4H
v1data = data.frame(A375_VAV1[wantGenes,c(2,6)], gene = wantGenes)
v1data = data.frame( v1data, bool = v1data[,"padj"] < 0.01)
v1data = data.frame(v1data, condition=rep("VAV1", 4))

r1data = data.frame(A375_R1_OE[wantGenes,c(2,6)], gene = wantGenes)
r1data = data.frame(r1data, bool= v1data[,"padj"] < 0.01)
r1data = data.frame(r1data, condition=rep("RAC1", 4))

dat = rbind(v1data,r1data)
dat[,"condition"] = factor(dat[,"condition"], levels = c("VAV1", "RAC1"))
# dat[,"padj"] = -log10(dat[,"padj"])
dat[,"bool"] = factor(dat[,"bool"], levels = c("TRUE", "FALSE") )

pdf("vav1OECands.pdf",5.5,3)
ggplot(dat,aes(x=gene,y=log2FoldChange,fill=bool)) + 
  geom_col(position="dodge",width=0.4) +
  coord_flip() + scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  facet_grid(condition~.)+
  theme(strip.text.y = element_text(angle = 0))  + 
  theme(text = element_text(size = 30), axis.title.x=element_blank(), axis.title.y=element_blank()) +
  ylim(0,2.5) + theme_classic()
dev.off()

