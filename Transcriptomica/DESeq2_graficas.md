### Gráficas de DGE
Con los datos obtenidos del [análisis con DESeq2](DESeq2.md), podemos generar varias gráficas.

***
#### Gráfica de componentes principales (PCoA)
Primeros debemos transformar los datos en la matriz `dds` logarítmicamente:
```bash
rld <- rlog(dds, blind=TRUE)
```
Ahora creemos la gráfica separando por tratamiento:
```bash
plotPCA(rld, intgroup="Treatment")
```
***
#### Heatmap
**Nota.** Asegurarnos de haber creado el objeto `rld` como se especifica arriba.

Primero debemos cargar la librería y hacer unos análisis previos:
```bash
library(pheatmap)
```
Crear una matriz

```bash
rld_mat <- assay(rld)
```
Calcular correlaciones pairwise:

```bash
rld_cor <- cor(rld_mat)
```
Grafiquemos:
```bash
pheatmap(rld_cor, fontsize = 15)
```
La gráfica podemos mejorarla para destacar mejor las similitudes:
```bash
library(RColorBrewer)
heat.colors <- brewer.pal(4, "Greys")
```
Si no está instalada la librería pheatmap, ahy que instalarla en R.
```bash
library(pheatmap)
pheatmap(rld_cor, color = heat.colors, border_color="darkgrey", fontsize = 14, fontsize_row = 14, height=20)
```
***
#### Gráficas tipo volcán
Proceso tomado de: https://biostatsquid.com/volcano-plots-r-tutorial/

Podemos cargar en RStudio el script [Volcano_plot.R](scripts/Volcano_plot.R) para correrlo desde ahí.

```bash
dge <- read.csv("DGE.csv")
```
Añadir una nueva columna para los valores de DGE:
```bash
dge$diffexpressed <- "NO"
```
Marcar como "UP" si log2Foldchange > 0.6 and pvalue < 0.05
```bash
dge$diffexpressed[dge$log2FoldChange > 0.6 & dge$pvalue < 0.05] <- "UP"
```
Marcar como "DOWN" si log2Foldchange < -0.6 and pvalue < 0.05
```bash
dge$diffexpressed[dge$log2FoldChange < -0.6 & dge$pvalue < 0.05] <- "DOWN"
```
Añadir una nueva columna "delabel" que contendrá el nombre de los 10 primeros DGE genes (NA en caso contrario).
```bash
dge$delabel <- ifelse(dge$X %in% head(dge[order(dge$padj), "X"], 10), dge$X, NA)
```

Creemos un tema para la gráfica:
```bash
library(ggplot2)
library(ggrepel)
```
```bash
theme_set(theme_classic(base_size = 20) +
 theme(axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(0.7), color = 'black'),
 axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(0.7), color = 'black'), plot.title = element_text(hjust = 0.5)
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
resLFC <- lfcShrink(dds, coef="Treatment_MM9_vs_TSB", type="apeglm")
par(mfrow = c(1, 2))
plotMA(resLFC, main="Shrinkage of LFCs", ylim=c(-4,4))
plotMA(res, main="No shrinkage of LFCs", ylim=c(-4,4))
```
***
