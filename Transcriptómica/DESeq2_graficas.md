### Gráficas de DGE
Con los datos obtenidos del [análisis con DESeq2](Transcriptómica/DESeq2.md), podemos generar varias gráficas.

***
#### Gráfica PCoA
Primeros debemos transformar los datos logarítmicamente:
```bash
rld <- rlog(dds, blind=TRUE)
```
Ahora creemos la gráfica:
```bash
plotPCA(rld, intgroup="Treatment")
```
***
#### Heatmap
**Nota.** Asegurarnos de haber creado el objeto *rld* como se especifica arriba.

Primero debemos cargar la librería y hacer unos análisis previos:
```bash
library(pheatmap)
```
```bash
rld_mat <- assay(rld) # crear una matriz
rld_cor <- cor(rld_mat) # calcular correlaciones pairwise
```
Grafiquemos:
```bash
pheatmap(rld_cor, fontsize = 15)
```
la gráfica podemos embellecerla con otros colores:
```bash
library(RColorBrewer)
heat.colors <- brewer.pal(4, "Greys")
```
```bash
pheatmap(rld_cor, color = heat.colors, border_color="darkgrey", fontsize = 14, fontsize_row = 14, height=20)
```
***
#### Gráficas tipo volcán
Proceso tomado de: https://biostatsquid.com/volcano-plots-r-tutorial/

Podemos cargar en RStudio el script [Volcano_plot.R](Transcriptómica/Volcano_plot.R) para correrlo desde ahí.

```bash
dge <- read.csv("DGE.csv")

dge$diffexpressed <- "NO" # Add a new column for the DGE values. NO by default.

dge$diffexpressed[dge$log2FoldChange > 0.6 & dge$pvalue < 0.05] <- "UP" # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"

dge$diffexpressed[dge$log2FoldChange < -0.6 & dge$pvalue < 0.05] <- "DOWN" # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"

dge$delabel <- ifelse(dge$X %in% head(dge[order(dge$padj), "X"], 10), dge$X, NA) # Create a new column "delabel" to dge, that will contain the name of the top 10 differentially expressed genes (NA in case they are not)
```
Creemos un tema para la gráfica:
```bash
theme_set(theme_classic(base_size = 20) +
 theme(axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(0.7), color = 'black'),
 axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(0.7), color = 'black'), plot.title = element_text(hjust = 0.5) # Biostatsquid theme
 ))
```
Creemos la gráfica:
```bash
ggplot(dge, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label = delabel)) +
   geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
   geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
   geom_point(size = 2) +
   scale_color_manual(values = c("#00AFBB", "grey", "#db487e"), labels = c("Downregulated", "Not signifficant", "Upregulated")) +
   labs(y= "Confidence (-Log10 adjusted p values)", x = "Log2 Ratio (Fold Change)")+
   ggtitle("Differential expression between treatments.") +
   geom_text_repel(max.overlaps = Inf)
```

Shrinkage of effect size (LFC) helps to remove the low count genes (by shrinking towards zero)
```bash
resLFC <- lfcShrink(dds, coef="Treatment_Seawater_vs_TSB", type="apeglm")
par(mfrow = c(1, 2))
plotMA(resLFC, main="Shrinkage of LFCs", ylim=c(-4,4))
plotMA(res, main="No shrinkage of LFCs", ylim=c(-4,4))
```
***
