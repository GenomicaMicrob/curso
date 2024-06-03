### Análisis de Genes Diferencialemnte Expresados (DEG) con [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) en RStudio

Procedimiento tomado de: https://www.reneshbedre.com/blog/deseq2.html
***
A partir del análisis de transcritos realizado [anteriormente con FADU](FADU.md), podemos calcular los DEG a partir del archivo `counts.csv` generado. DESeq2 es un **paquete de R** y por lo tanto, si no lo tenemos instalado en R, hay que hacerlo en **RStudio**.

Creemos un nuevo proyecto de **RStudio** en la carpeta donde tengamos los datos de FADU e instalemos *DESeq2*.

 ```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```
Una vez instalado carguemos los paquetes necesarios, si alguno no esta, hay que instalarlo desde el mismo RStudio:

```
library("DESeq2")
library(ggplot2)
library(ggrepel)
```
***
Importemos los archivos y creemos las matrices de datos:

**Data**

```bash
count_matrix <- read.csv("counts.csv", row.names = 1)
```
Si el archivo esta **delimitado por tabuladores**, usar entonces `read.delim` en lugar de `read.csv`

```bash
count_matrix <- round(count_matrix, 0) # round numbers
```
```bash
head(count_matrix, 2) # view first two rows
```
**Metadata**

El archivo de metadatos (`metadata.csv`) tiene los siguientes datos delimitados por comas:
 ```
Sample,Treatment,Replicate
MM9r1,MM9,1
MM9r2,MM9,2
MM9r3,MM9,3
TSBr1,TSB,1
TSBr2,TSB,2
TSBr3,TSB,3
```
Importémoslo a una matriz de datos:
```bash
coldata <- read.csv("metadata.csv", header = TRUE, row.names = 1)
```
```bash
coldata <- data.frame(coldata) # conversión a dataframe
```
```bash
coldata$Treatment <- as.factor(coldata$Treatment)
```
```bash
coldata # veamos los metadatos
```

Revisemos si tenemos el **mismo nombre y número** de las columnas (muestras) en los datos que en las lineas de los metadatos:

```bash
all(rownames(coldata) %in% colnames(count_matrix))
```
```bash
all(rownames(coldata) == colnames(count_matrix))
```
**Ambos** deben ser verdadero (`TRUE`). Si alguno de los resultados anteriores es `FALSE`,  debemos ver porque hay una discrepancia; normalmente los metadatos no están bien capturados.
***
### Análisis DGE
```bash
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~ Treatment)
```
Construcción del `DESeqDataSet` para el análisis de DGE; se filtran (quitan) genes con valores menores a 10 lecturas (se puede cambiar).

```bash
dds <- dds[rowSums(counts(dds)) >= 10,]
```
*Pre-filtering helps to remove genes that have very few mapped reads, reduces memory, and increases the speed of the DESeq2 analysis*

Seleccionar cual es el valor para referencia (medio TSB) contra el cual se compararán los datos.
```bash
dds$Treatment <- relevel(dds$Treatment, ref = "TSB")
```
*REFERENCE. Select the reference (ref parameter) level for condition comparisons. The comparisons of other conditions will be compared against this reference i.e, the log2 fold changes will be calculated based on ref value (infected/control). If this parameter is not set, comparisons will be based on alphabetical order of the levels*

**Análisis de expresión diferencial**
```bash
dds <- DESeq(dds)
```
**Resultados**
```bash
resultsNames(dds)
```
Mandemos los resultados a una nueva tabla (`res`); por default hay una filtración de genes con baja expresión (se puede quitar el filtro con `independentFiltering=FALSE`); podemos añadir el nombre del análisis también:
```bash
res <- results(dds, name="Treatment_MM9_vs_TSB")
res
```
Hagamos un ordenamiento de los resultados por el valor de *p* ajsutado (método de Benjamini-Hochberg FDR):
```bash
res[order(res$padj),]
```
*Order gene expression table by adjusted p value (Benjamini-Hochberg FDR method). Note: You may get some genes with p value set to NA. This is due to all samples have zero counts for a gene or there is extreme outlier count for a gene or that gene is subjected to independent filtering by DESeq2.*

Salvemos los resultados DGE a una tabla (`DGE.csv):
```bash
write.csv(as.data.frame(res[order(res$padj),] ), file = "DGE.csv")
```
Obtengamos un resumen de DGE con valores de *p* ajustados y solo los que son igual o menores a 0.05:

```bash
summary(results(dds, alpha=0.05))
```
Ahora obtengamos los valores normalizados:
```bash
normalized_counts <- counts(dds, normalized=TRUE)
head(normalized_counts)
```
Listo, hemos terminado.
***
### Gráficas tipo volcán
Proceso tomado de: https://biostatsquid.com/volcano-plots-r-tutorial/
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
