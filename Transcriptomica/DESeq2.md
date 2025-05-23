### Análisis de Genes Diferencialmente Expresados (DEG) con [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) en RStudio

Procedimiento tomado de: https://www.reneshbedre.com/blog/deseq2.html
***
A partir del análisis de transcritos realizado [anteriormente con FADU](Transcriptomica.md), podemos calcular los DEG a partir del archivo `counts.csv` generado. DESeq2 es un **paquete de R** y por lo tanto, si no lo tenemos instalado en R, hay que hacerlo en **RStudio**.

**Importante!** Creemos un nuevo proyecto de RStudio **en la carpeta** donde tengamos los datos de FADU e instalemos *DESeq2*.

**Nota.** DESeq2 tarda mucho en instalarse.

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

**Importante**. Si anteriormente no habíamos redondeado los valores, debemos hacerlo ahora:

```bash
count_matrix <- round(count_matrix, 0) # redondeo de decimales
```
Revisemos que los datos están redondeados sin cifras decimales.

```bash
head(count_matrix, 2)
```
**Metadata**

El archivo de metadatos [metadata.csv](data/metadata.csv) tiene los siguientes datos delimitados por comas:
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
Conversión a dataframe

```bash
coldata <- data.frame(coldata)
```
Conversión de la columna `Treatment` a un factor

```bash
coldata$Treatment <- as.factor(coldata$Treatment)
```
Veamos los metadatos

```bash
coldata
```
**Validación**

Revisemos si tenemos el **mismo nombre y número** de las columnas (muestras) en los datos que en las lineas de los metadatos:

```bash
all(rownames(coldata) %in% colnames(count_matrix))
```
```bash
all(rownames(coldata) == colnames(count_matrix))
```
**Ambos** deben ser verdadero (`TRUE`)! Si alguno de los resultados anteriores es `FALSE`,  debemos ver porque hay una discrepancia; normalmente los metadatos no están bien capturados. No podremos continuar si hay discrepancias.
***

### Análisis DGE
**Tip.** Podemos importar el script [DESeq2_script.R](scripts/DESeq2_script.R) a RStudio y desde ahí ir corriendo todos los pasos.

```bash
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~ Treatment)
```
Construcción del `DESeqDataSet` para el análisis de DGE; se filtran (quitan) genes con valores menores a 10 lecturas (se puede cambiar).

```bash
dds <- dds[rowSums(counts(dds)) >= 10,]
```
*Pre-filtering helps to remove genes that have very few mapped reads, reduces memory, and increases the speed of the DESeq2 analysis*

Seleccionar cual es el valor para referencia (medio TSB) contra el cual se compararán los datos:
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
```
Hagamos un **ordenamiento** de los resultados por el valor de *p* ajustado (método de Benjamini-Hochberg FDR):

```bash
res[order(res$padj),]
```
*Order gene expression table by adjusted p value (Benjamini-Hochberg FDR method). Note: You may get some genes with p value set to NA. This is due to all samples have zero counts for a gene or there is extreme outlier count for a gene or that gene is subjected to independent filtering by DESeq2.*

Salvemos los resultados DGE a una tabla (`DGE.csv`) en la carpeta:

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
Listo, hemos terminado, ahora podemos **graficar** los datos siguiendo la guía [DESeq2_graficas](DESeq2_graficas.md) y analizar los valores obtenidos en el archivo `DGE.csv` generado.
***
