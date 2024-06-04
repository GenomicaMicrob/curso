# ----- Volcano plot -> https://biostatsquid.com/volcano-plots-r-tutorial/
library(ggplot2)
library(ggrepel)
library("DESeq2")

dge <- read.csv("DGE.csv") # Import values as data frame
dge$diffexpressed <- "NO" # Add a new column for the DGE values. NO by default.
dge$diffexpressed[dge$log2FoldChange > 0.6 & dge$pvalue < 0.05] <- "UP" # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
dge$diffexpressed[dge$log2FoldChange < -0.6 & dge$pvalue < 0.05] <- "DOWN" # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
dge$delabel <- ifelse(dge$X %in% head(dge[order(dge$padj), "X"], 10), dge$X, NA) # Create a new column "delabel" to dge, that will contain the name of the top 10 differentially expressed genes (NA in case they are not)

theme_set(theme_classic(base_size = 20) +
 theme(axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(0.7), color = 'black'),
 axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(0.7), color = 'black'), plot.title = element_text(hjust = 0.5) # Biostatsquid theme
 ))

volcano <- ggplot(dge, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label = delabel)) +
   geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
   geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
   geom_point(size = 2) +
   scale_color_manual(values = c("#00AFBB", "grey", "#db487e"), labels = c("Downregulated", "Not signifficant", "Upregulated")) +
   labs(y= "Confidence (-Log10 adjusted p values)", x = "Log2 Ratio (Fold Change)")+
   ggtitle("Differential expression between treatments.", subtitle = output.file) +
   geom_text_repel(max.overlaps = Inf) # VolcanoPlot, needs chanes to scale_color and to title
volcano
