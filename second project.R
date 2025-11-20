#PACKAGES INSTALLATION
install.packages(c(
  "tidyverse",
  "matrixStats",
  "cowplot",
  "tibble",
  "ggplot2",
  "plotly",
  "RColorBrewer",
  "gplots",
  "gameofthrones",
  "d3heatmap",
  "gprofiler2",
  "readxl"
))
install.packages("BiocManager")
BiocManager::install(c(
  "edgeR",
  "limma",
  "biomaRt",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "clusterProfiler",
  "enrichplot",
  "GO.db",
  "BiocParallel"
))
a
BiocManager::install("org.Hs.eg.db")
a
#import the Packages
library(tidyverse)
library(dplyr) 
library(matrixStats) 
library(cowplot) 
library(tibble) 
library(ggplot2) 
library(plotly) 
library(RcolorBrewer)
library(gplots) 
library(gameofthrones) 
library(d3heatmap) 
library(gprofiler2) 
library(readxl) 
library(edgeR)
library(limma)
library(biomaRt)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(GO.db)


#Setting Working Directory
setwd("C:\\Users\\HP\\R content")
rnadata <- read.delim("GSE287775_All_gene_raw_counts.tsv")
head(rnadata)

#LOADING DATA

rnadata$Control1 <- rnadata$sample.MM.1S_WT_GFP_1
rnadata$Control2 <- rnadata$sample.MM.1S_WT_GFP_2
rnadata$Control3 <- rnadata$sample.MM.1S_WT_GFP_3
rnadata$Treatment1 <- rnadata$sample.MM.1S_E9_GFP_1
rnadata$Treatment2 <- rnadata$sample.MM.1S_E9_GFP_2
rnadata$Treatment3 <- rnadata$sample.MM.1S_E9_GFP_3
View(rnadata)
#Changing col name from ensembl_gene_id to gene_id
colnames(rnadata)[colnames(rnadata)=="ensembl_gene_id"] <- "gene_id"
#Select needed data
rnadata2 <- select(rnadata,gene_id,Control1,Control2,Control3,
                   Treatment1,Treatment2,Treatment3)

View(rnadata2)
#preparing count matrix dataframe
mycounts <- as.data.frame(rnadata2)

#Setting Column Order
column_order <- c("Control1", "Control2","Control3",
                  "Treatment1","Treatment2","Treatment3")
#Making row.names unique and setting gene_id as rownames
rownames(mycounts) <- make.unique(mycounts$gene_id)
#Setting GeneID as rownames
#rownames(mycounts) <- mycounts$gene_id
View(rownames(mycounts))
#Removing GeneID column and setting values as row names
mycounts2 <- mycounts[,-1]

View(mycounts2)

#Defining Experimental Groups
group <- factor(c("Control","Control","Control",
                  "Treatment","Treatment","Treatment"))
levels(group)
nlevels(group)
#Creating the experimental design data
targets <- data.frame(row.names =colnames(mycounts2),group)
#SET THE ROWNAMES OF THE targets dataframe as the ROWNAMES
sampleLabels <- rownames(targets)
view(sampleLabels)

#making DGEList object
myDGEList <- DGEList(mycounts2)


#Generating Count Per Million (CPM) values
cpm <- cpm(myDGEList)

#Confirm cpm values calculation
colSums(cpm)

log2.cpm <- cpm(myDGEList,log=TRUE)

#Converting data matrix to dataframe
log2.cpm.df <- as_tibble(log2.cpm,rownames="gene_id")
head(log2.cpm.df)
colnames(log2.cpm.df) <- c("gene_id",sampleLabels)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df,cols = -1,
                                  names_to = "samples",values_to = "expression")

#Plotting violin plot for all samples
plot1 <- ggplot(log2.cpm.df.pivot)+aes(x=samples,y=expression,fill = samples)+
  geom_violin(trim = FALSE,show.legend = FALSE)+
  stat_summary(fun = "median",geom = "point",shape=95,size=10,
               coor="black",show.legend = FALSE)+
  labs(y="log2 expression",x="Sample",
       title="Log2 Counts per Million (CPM) of All Samples",
       subtitle = "unfiltered,non-normalized")+
  theme_bw()+
  coord_flip()


ggplot(log2.cpm.df.pivot)+aes(x=samples,y=expression,fill = samples)+
  geom_violin(trim = FALSE,show.legend = FALSE)+
  stat_summary(fun = "median",geom = "point",shape=95,size=10,coor="black",
               show.legend = FALSE)+
  labs(y="log2 expression",x="Sample",
       title="Log2 Counts per Million (CPM) of All Samples",
       subtitle = "unfiltered,non-normalized")+
  theme_bw()+
  coord_flip()
#filtering 
cpm <- cpm(myDGEList)
filt_threshold <- rowSums(cpm>1)>=2
myDGEList.filtered <-myDGEList[filt_threshold,]  
#myDGEList.filtered <- myDGEList #not filtered for GSEA
log2.cpm.filtered <- cpm(myDGEList.filtered,log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered,rownames="gene_id")
View(log2.cpm.filtered.df)
View(colnames(log2.cpm.filtered.df))
colnames(log2.cpm.filtered.df) <- c("gene_id",sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df,
                                           cols = -1,
                                           names_to = "samples",
                                           values_to = "expression")
View(log2.cpm.filtered.df.pivot)
plot2 <- ggplot(log2.cpm.filtered.df.pivot)+aes(x=samples,
                                                y=expression,fill = samples)+
  geom_violin(trim = FALSE,show.legend = FALSE)+
  stat_summary(fun = "median",geom = "point",shape=95,size=10,coor="black",
               show.legend = FALSE)+
  labs(y="log2 expression",x="Sample",
       title="Log2 Counts per Million (CPM) of All Samples",
       subtitle = "filtered,non-normalized")+
  theme_bw()+
  coord_flip()
plot2

#Normalization
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered,method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm,log=TRUE)
View(log2.cpm.filtered.norm)
log2.cpm.filtered.norm.df <- as.tibble(log2.cpm.filtered.norm,rownames="gene_id")

colnames(log2.cpm.filtered.norm.df) <- c('gene_id',sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df,
                                                cols = -1,names_to = "samples",
                                                values_to = "expression")
plot3 <- ggplot(log2.cpm.filtered.norm.df.pivot)+
  aes(x=samples,y=expression,fill = samples)+
  geom_violin(trim = FALSE,show.legend = FALSE)+
  stat_summary(fun = "median",
               geom = "point",
               shape=95,
               size=10,
               coor="black",
               show.legend = FALSE)+
  labs(y="log2 expression",
       x="Sample",
       title="Log2 Counts per Million (CPM) of All Samples",
       subtitle = "filtered,TMM normalized")+
  theme_bw()+
  coord_flip()
plot3
norm_plot <- plot_grid(plot1,plot2,plot3,labels = c("A","B","C"),
                       label_size = 12)
View(norm_plot)
ggsave("norm_plot.png",norm_plot,width = 10,height = 9,dpi = 300)
ggsave("plot_unfil_fil_norm.png",norm_plot,width = 10,height = 9,dpi = 300)

#hierarchical clustering
#setting cohort and group

group <- targets$group
group <- factor(group)

#hierarchical clustering
png("Hierarchial clustering.png",width = 400,height = 500)

View(log2.cpm.filtered.norm)
t(log2.cpm.filtered.norm)
distance <- dist(t(log2.cpm.filtered.norm),
                 method = "euclidean")
View(distance)
clusters <- hclust(distance,method = "complete")
View(clusters) 
plot(clusters,labels=sampleLabels)
plot(clusters,labels=c("Control1","Control2","Control3",
                       "Treatment1","Treatment2","Treatment3"))
dev.off()  

#Principal Component Analysis
pca.res <- prcomp(t(log2.cpm.filtered.norm),
                  scale. = F,retx = T)
pca.res$x
View(pca.res$x)
pc.var <- pca.res$sdev^2
pc.per <- round(pc.var/sum(pc.var)*100,1)

#visualizing PCA result
pca.res.df <- as.tibble(pca.res$x)
png("PCA.png",width = 400,height = 500)
ggplot(pca.res.df)+
  aes(x=PC1,y=PC2,labels=sampleLabels,color=group)+
  geom_point(size=4)+
  #geom_label() +
  xlab(paste0("PC1 (",pc.per[1],"%",")"))+
  ylab(paste0("PC2 (",pc.per[2],"%",")"))+
  labs(title = "Principal Component Analysis")+
  theme_bw()

dev.off()

#Differential Gene Expression Analysis
#Matrix Design
group <- factor(targets$group)
design <- model.matrix(~0+group)
colnames(design) <- levels(group)

#Getting the Mean Variance Relationships with Voom
voom_result <- voom(myDGEList.filtered.norm,design)

# Fitting Linear Model to Data
model_fit <-lmFit(voom_result,design)
View(model_fit)

#Contrast Matrix
contrast_matrix <- makeContrasts(infection=Treatment - Control,levels = design)

#Extracting Linear Model Fit
linear_model <- contrasts.fit(model_fit,contrast_matrix)

#Getting BAYESIAN STAT For Linear Fit Model
bay_stat <- eBayes(linear_model)

#Differential Expressed Genes
mygenes <- topTable(bay_stat,adjust="BH",coef = 1,number = Inf,sort.by = "P")
View(mygenes)

#Extract Differentially Expressed Genes
diff_exp_genes <- mygenes[(mygenes$adj.P.Val< 0.05 &
                             mygenes$logFC > 1) | (mygenes$adj.P.Val < 0.05 &
                                                     mygenes$logFC < -1),]
mygenes[(mygenes$adj.P.Val< 0.01 & mygenes$logFC > 1) | (mygenes$adj.P.Val <
                                                           0.01 & mygenes$logFC <
                                                           -1),]

#Extracting Ensemble_ids
ensemble_ids <- gsub("\\..*","",row.names(diff_exp_genes))
#add the modified Ensemble_IDs as a new column in the diff_exp_genes dataframe
diff_exp_genes$modified_ensembl <- ensemble_ids

#Connect to Ensembl's Homo sapiens dataset using the biomaRt package

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
View(ensembl)
#Ensembl Query to get homo sapiens gene official names
result <- getBM(mart = ensembl,
                attributes = c("ensembl_gene_id","hgnc_symbol"),
                filters = "ensembl_gene_id",
                values = ensemble_ids)
#merge on 'ensebl_gene_id'colum and keep all rows from diff_exp_genes dataframe
annotated_deg <- merge(diff_exp_genes, result, by.x ="modified_ensembl",
                       by.y ="ensembl_gene_id", all.x=TRUE)
#Replace rows where hgnc_symbol is NA with modified_ensembl values
annotated_deg$hgnc_symbol[is.na(annotated_deg$hgnc_symbol)] <- annotated_deg$
  modified_ensembl[is.na(annotated_deg$hgnc_symbol)] 
#Check Duplicated gene
#dup_gene_ids <- annotated_degunf$hgnc_symbol[duplicated(annotated_degunf$
hgnc_symbol|duplicated(annotated_degunf$hgnc_symbol,fromLast = TRUE)
#Remove all the rows with duplicate genes
#annotated_deg <- annotated_degunf[!annotated_degunf$hgnc_symbol %in%
dup_gene_ids,]
#set annotated_deg row names as hgnc_symbol
#row.names(annotated_deg) <- annotated_deg$hgnc_symbol
#set annotated_deg rownames as hgnc symbol

annotated_deg$hgnc_symbol <- as.character(annotated_deg$hgnc_symbol)
annotated_deg$unique_symbol <- make.unique(annotated_deg$hgnc_symbol)
rownames(annotated_deg) <- annotated_deg$unique_symbol

#row.names(annotated_deg) <- annotated_deg$hgnc_symbol
#removing the hgnc_symbol and modified_emsembl_id
colnames(annotated_deg) %in% c("modified_ensembl",
                               "hgnc_symbol","unique_symbol")
annotated_deg <- annotated_deg[,!(colnames(annotated_deg) %in%
                                    c("modified_ensembl",
                                      "hgnc_symbol","unique_symbol"))]

#extract upregulated & downregulated genes
upregulated_genes <-annotated_deg[annotated_deg$adj.P.Val< 0.05 &
                                    annotated_deg$logFC > 1, ]
downregulated_genes <-annotated_deg[annotated_deg$adj.P.Val< 0.05 &
                                      annotated_deg$logFC < -1, ]
uupregulated <- upregulated_genes %>% 
  arrange(desc(logFC))
View(uupregulated)
head(uupregulated,10)

#write total, upregulated &downregulated degs to csv
write.csv(annotated_deg, file = "total_degs.csv")
write.csv(annotated_deg, file = "upregulated_degs.csv")
write.csv(annotated_deg, file = "downregulated_degs.csv")

#Making Volcano Plot For Differentially Expressed Genes

#convert to a tibble
degs <- mygenes %>% 
  as_tibble(rownames = "gene_id")
#degs <- (as_tibble(mygenes,rownames = "GeneID"))
#convert to volcano plot
png("volcanoplots.png",width = 500,height = 500)

ggplot(degs) + aes(y = -log10(adj.P.Val), x = logFC)+
  geom_point(aes(colour = ifelse(adj.P.Val<0.05 & logFC > 1, "upregulated",
                                 ifelse(adj.P.Val<0.05 & logFC < -1,
                                        "downregulated","Not Significant"))),
             size = 2) + geom_hline(yintercept =  -log10(0.05),
                                    linetype = "longdash", colour = "grey",
                                    linewidth = 1)+
  geom_vline(xintercept = 1, linetype = "longdash", colour = "grey",
             linewidth = 1)+
  geom_vline(xintercept = -1, linetype = "longdash", colour = "grey",
             linewidth = 1)+
  labs(title = "Multiple Myeloma vs Control",
       subtitle = "741 upregulated genes and 682 downregulated genes")+ 
  theme_bw()+scale_colour_manual(values = c("upregulated"= "red",
                                            "downregulated"="blue",
                                            "Not Significant"="grey"))+
  theme(legend.title = element_blank())
dev.off()

#plot points with conditional coloring based on thresholds
#ggplot(degs) + aes(y = -log10(adj.P.Val), x = logFC)+
  #geom_point(aes(colour = ifelse(adj.P.Val<0.05 & logFC > 1, "upregulated",
                                 #ifelse(adj.P.Val<0.05 & logFC < -1,
                                        "#downregulated","Not Significant"))),
             #setting p-value threshold:p<0.05
             #geom_hline(yintercept =  -log10(0.05),
                        #linetype = "longdash", colour = "grey",
                        #linewidth = 1)
             
             #Setting fold change threshold: |log2FC|
             #geom_vline(xintercept = 1, linetype = "longdash", colour = "grey",
                        #linewidth = 1)+
               #geom_vline(xintercept = -1, linetype = "longdash", colour = "grey",
                          #linewidth = 1)+
               #Adding titles and captions
               #labs(title = "Prostate Cancer vs Control",
                    #subtitle = "4 upregulated genes and 3 downregulated genes")  
             #Apply theme
             #theme_bw()
             #Customize the colours
            # scale_colour_manual(values = c("upregulated"= "red",
                                            "downregulated"="blue",
                                            "Not Significant"="grey"))
             #Remove legend title
            # theme(legend.title = element_blank())
             
#GENE ONTOLOGY
human_mart <- useMart(biomart = "ensembl",
                                   host = "ensembl.org",
                                   dataset = "hsapiens_gene_ensembl",)
annot_diff <- getBM(mart = human_mart,
                                 attributes = c("hgnc_symbol","entrezgene_id","description"),
                                 filters = "hgnc_symbol",
                                 values = row.names(annotated_deg))
annot_diff$entrezgene_id <- as.character(annot_diff$entrezgene_id)
#Biological Process Analysis
ora_analysis_BP <- enrichGO(gene = annot_diff$entrezgene_id,
                                         OrgDb = org.Hs.eg.db,
                                         keyType = "ENTREZID",
                                         ont = "BP",
                                         pAdjustMethod = "BH",
                                         qvalueCutoff = 0.05,
                                         readable = FALSE,
                                         pool = FALSE)
ora_analysis_BP_table <- as.data.frame(ora_analysis_BP)
View(ora_analysis_BP_table)
ora_analysis_BP_final <- clusterProfiler::simplify(ora_analysis_BP)
write_delim(
               x = as.data.frame(ora_analysis_BP@result),
               path = "Biological process.csv",
               delim = ","
             )
#dot plot
png("BP dotplot.png",width = 470,height = 550)
dotplot(ora_analysis_BP_final,showCategory = 10)
dev.off()
#barplot
png("BP barplot.png",width = 470,height = 550)
barplot(ora_analysis_BP_final,showCategory = 10)
dev.off()
#Cellular component
ora_analysis_CC <- enrichGO(gene = annot_diff$entrezgene_id,
                                         OrgDb = org.Hs.eg.db,
                                         keyType = "ENTREZID",
                                         ont = "CC",
                                         pAdjustMethod = "BH",
                                         qvalueCutoff = 0.05,
                                         readable = FALSE,
                                         pool = FALSE)
ora_analysis_CC_table <- as.data.frame(ora_analysis_CC)
View(ora_analysis_CC_table)
ora_analysis_CC_final <- clusterProfiler::simplify(ora_analysis_CC)
write_delim(
            x = as.data.frame(ora_analysis_CC@result),
            path = "Cellular component.csv",
            delim = ","
             )
#dotplot
png("CC dotplot.png",width = 470,height = 550)
dotplot(ora_analysis_CC_final,showCategory = 10)
dev.off()
#barplot
png("CC barplot.png",width = 470,height = 550)
barplot(ora_analysis_CC_final,showCategory = 10)
dev.off()
#molecular Function
ora_analysis_MF <- enrichGO(gene = annot_diff$entrezgene_id,
                                         OrgDb = org.Hs.eg.db,
                                         keyType = "ENTREZID",
                                         ont = "MF",
                                         pAdjustMethod = "BH",
                                         qvalueCutoff = 0.05,
                                         readable = FALSE,
                                         pool = FALSE)
ora_analysis_MF_table <- as.data.frame(ora_analysis_MF)
View(ora_analysis_MF_table)
ora_analysis_MF_final <- clusterProfiler::simplify(ora_analysis_MF)
write_delim(
               x = as.data.frame(ora_analysis_MF@result),
               path = "Molecular Function.csv",
               delim = ","
             )
#dotplot
png("MF dotplot.png",width = 470,height = 550)
dotplot(ora_analysis_MF_final,showCategory = 10)
dev.off()
#barplot
png("MF barplot.png",width = 470,height = 550)
barplot(ora_analysis_MF_final,showCategory = 10)
dev.off()
#KEGG analysis
ora_analysis_KEGG <- enrichKEGG(gene = annot_diff$entrezgene_id,
                                             organism ='hsa', 
                                             pAdjustMethod = "BH",
                                             qvalueCutoff = 0.05)
ora_analysis_KEGG_table <- as.data.frame(ora_analysis_KEGG)
write.csv(ora_analysis_KEGG_table, file = "KEGG pathway.csv")
png("KEGG dotplot.png",width = 470,height = 550)
dotplot(ora_analysis_KEGG,showCategory = 10)
        dev.off()
#barplot
png("KEGG barplot.png",width = 470,height = 550)
barplot(ora_analysis_KEGG,showCategory = 10)
             dev.off()
#Gene Set Enrichment Analysis(GSEA)
#data processing for GSEA
myDGEList.unfiltered <- myDGEList #not filtered for GSEA
log2.cpm.unfiltered <- cpm(myDGEList.unfiltered,log=TRUE)
log2.cpm.unfiltered.df <- as_tibble(log2.cpm.unfiltered,rownames="gene_id")
View(colnames(log2.cpm.unfiltered.df))
colnames(log2.cpm.unfiltered.df) <- c("gene_id",sampleLabels)
log2.cpm.unfiltered.df.pivot <- pivot_longer(log2.cpm.unfiltered.df,
                                           cols = -1,
                                           names_to = "samples",
                                           values_to = "expression")
View(log2.cpm.unfiltered.df.pivot)


#Normalization
myDGEList.unfiltered.norm <- calcNormFactors(myDGEList.unfiltered,method = "TMM")
log2.cpm.unfiltered.norm <- cpm(myDGEList.unfiltered.norm,log=TRUE)
log2.cpm.unfiltered.norm.df <- as.tibble(log2.cpm.unfiltered.norm,rownames="gene_id")


colnames(log2.cpm.unfiltered.norm.df) <- c('gene_id',sampleLabels)
log2.cpm.unfiltered.norm.df.pivot <- pivot_longer(log2.cpm.unfiltered.norm.df,
                                                cols = -1,names_to = "samples",
                                                values_to = "expression")

#Getting the Mean Variance Relationships with Voom
voom_unfresult <- voom(myDGEList.unfiltered.norm,design)

# Fitting Linear Model to Data
model_fitunf <-lmFit(voom_unfresult,design)

View(model_fitunf)

#Contrast Matrix
contrast_matrixunf <- makeContrasts(infection=Treatment - Control,
                                    levels = design)

#Extracting Linear Model Fit
linear_modelunf <- contrasts.fit(model_fitunf,contrast_matrixunf)

#Getting BAYESIAN STAT For Linear Fit Model
bay_statunf <- eBayes(linear_modelunf)
#Differential Expressed Genes
mygenesunf <- topTable(bay_statunf,adjust="BH",coef = 1,
                       number = Inf,sort.by = "P")
 
#Extracting Ensemble_ids by removing version number
ensemble_idsunf <- gsub("\\..*","",row.names(mygenesunf))
#add the modified Ensemble_IDs as a new column in the diff_exp_genes dataframe
mygenesunf$modified_ensembl <- ensemble_idsunf

#Connect to Ensembl's Homo sapiens dataset using the biomaRt package

ensemblunf <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
View(ensembl)
#Ensembl Query to get homo sapiens gene official names
resultunf <- getBM(mart = ensemblunf,
                attributes = c("ensembl_gene_id","hgnc_symbol"),
                filters = "ensembl_gene_id",
                values = ensemble_idsunf)
#merge on 'ensebl_gene_id'colum and keep all rows from diff_exp_genes dataframe
annotated_degunf <- merge(mygenesunf, resultunf, by.x ="modified_ensembl",
                       by.y ="ensembl_gene_id", all.x=TRUE)
#set annotated_deg rownames as hgnc symbol
#Replace rows where hgnc_symbol is NA with modified_ensembl values
annotated_degunf$hgnc_symbol[is.na(annotated_degunf$
                                     hgnc_symbol)] <- annotated_degunf$
  modified_ensembl[is.na(annotated_degunf$hgnc_symbol)] 
#Check Duplicated gene
dup_gene_ids <- annotated_degunf$hgnc_symbol[duplicated(annotated_degunf$
                                                          hgnc_symbol)|
                                               duplicated(annotated_degunf$
                                                            hgnc_symbol,
                                                          fromLast = TRUE)]
#Remove all the rows with duplicate genes
annotated_degunf <- annotated_degunf[!annotated_degunf$hgnc_symbol %in%
                                       dup_gene_ids,]
#set annotted_degunf row names as hgnc_symbol
row.names(annotated_degunf) <- annotated_degunf$hgnc_symbol
#removing the hgnc_symbol and modified_emsembl_id
colnames(annotated_degunf) %in% c("modified_ensembl",
                               "hgnc_symbol")
annotated_degunf <- annotated_degunf[,!(colnames(annotated_degunf) %in%
                                    c("modified_ensembl","hgnc_symbol"))]

#Gene Set Enrichment Analysis
data <- annotated_degunf %>% 
  dplyr::arrange(desc(logFC))
gene_list <- data$logFC
View(gene_list)
names(gene_list) <- row.names(data)
head(gene_list)

gse <- gseGO(
  geneList = gene_list,
  ont = "All",
  keyType = "SYMBOL",
  nPerm = 1000,
  minGSSize = 3,
  maxGSSize = 100,
  pvalueCutoff = 0.05,
  verbose = TRUE,
  OrgDb = org.Hs.eg.db,
  pAdjustMethod = "none"
)            

#to filter by NES
gse@result <- gse@result %>% 
  arrange(desc(NES))
view(summary(gse))
gse_results <- gse@result
#save gse result to file as CSV
write.csv(gse_results,file = "GSEA_results.csv")

png("gse_dotplot.png",width = 500, height = 470)
dotplot(gse,showCategory = 10,split = ".sign") +
  facet_grid(.~.sign) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text.y = element_text(size = 7))
dev.off()

#plot specific GSE result
png("gsel.png",width =500,height = 400)
#Find the index,NES, and p.adjust values for the desired gene set
gene_set_index <- which(gse@result$Description=="outer kinetochore")
gene_set_NES <- gse@result$NES[gene_set_index]
gene_set_padj <- gse@result$p.adjust[gene_set_index]
#Create the plot with NES and p.adjust included in the title
gseaplot2(
  gse,
  geneSetID = gene_set_index,
  title = paste(
    "outer kinetochore",
    "\nNES =",round(gene_set_NES,2),
    ", p.adjust =", formatC(gene_set_padj,format = "e",digits = 2)
    
  )
)
dev.off()

png("gse2.png",width =500,height = 400)
#Find the index,NES, and p.adjust values for the desired gene set
gene_set_index <- which(gse@result$Description=="spindle assembly checkpoint signaling")
gene_set_NES <- gse@result$NES[gene_set_index]
gene_set_padj <- gse@result$p.adjust[gene_set_index]
#Create the plot with NES and p.adjust included in the title
gseaplot2(
  gse,
  geneSetID = gene_set_index,
  title = paste(
    "spindle assembly checkpoint signaling",
    "\nNES =",round(gene_set_NES,2),
    ", p.adjust =", formatC(gene_set_padj,format = "e",digits = 2)
    
  )
)
dev.off()

png("gse3.png",width =500,height = 400)
#Find the index,NES, and p.adjust values for the desired gene set
gene_set_index <- which(gse@result$Description=="regulation of mitotic sister chromatid separation")
gene_set_NES <- gse@result$NES[gene_set_index]
gene_set_padj <- gse@result$p.adjust[gene_set_index]
#Create the plot with NES and p.adjust included in the title
gseaplot2(
  gse,
  geneSetID = gene_set_index,
  title = paste(
    "regulation of mitotic sister chromatid separation",
    "\nNES =",round(gene_set_NES,2),
    ", p.adjust =", formatC(gene_set_padj,format = "e",digits = 2)
    
  )
)
dev.off()

png("gse4.png",width = 500,height = 400)
gene_set_index <- which(gse@result$Description=="")
gene_set_NES <- gse@result$NES[gene_set_index]
gene_set_padj <- gse@result$p.adjust[gene_set_index]
gseaplot2(
  gse,
  geneSetID = gene_set_index,
  title = paste(
    "regulation of mitotic sister chromatid separation",
    "\nNES =",round(gene_set_NES,2),
    ", p.adjust =", formatC(gene_set_padj,format = "e",digits = 2)
    
  )
)
dev.off()