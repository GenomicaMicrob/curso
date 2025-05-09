## Gráficas  para 16S

El script mg_classifer nos da ya varios archivos con los que podemos obtener gráficas con R de manera sencilla. Podemos generar las siguientes gráficas básicas para un análisis de microbiota por amplicones:

1. Curvas de rarefacción
2. Gráfica de barras apiladas a nivel phylum y familia
3. Heatmaps
4. Gráfica de componentes principales (PCoA)
5. Boxplot
***
### Curvas de rarefacción

Este tipo de gráficas nos permite ver si el esfuerzo de secuenciación (muestreo) fue suficiente para tener una buena representación de las OTUs en nuestra muestra, si se llega a una asíntota en la gráfica.

Para esta gráfica usaremos en RStudio la paquetería [ampvis2](https://kasperskytte.github.io/ampvis2/articles/ampvis2.html) y el archivo `/OTU_tables/qiime_ampvis2/OTU_table.tsv`

Creado el proyecto carguemos la librería de ampvis2

```bash
library(ampvis2)
```
Importemos ahora la **tabla de OTUs** y los **metadatos** a RStudio:

```bash
df <- read.delim("OTU_tables/qiime_ampvis2/OTU_table.tsv", , comment.char="#")
```

```bash
md <-
```

ya importados ambos archivos podemos generar las curvas con ampvis2:

```bash

```
