# Clasificación taxonómica de metagenomas
Mediante el programa [Focus](https://github.com/metageni/FOCUS) podemos hacer una **rápida clasificación taxonómica** de muestras metagenómicas. La base de datos con la que se clasifica no es muy extensa y tampoco muy actualizada pero nos permite hacer una clasificación aceptable.
***
#### Archivo y scripts necesarios para este ejercicio

- [dietas.V2.tar.gz](https://figshare.com/s/4e700c8c9ce853e74827)
- `multiple_fastq_merger-converter`
- `focus2ampvis`
- [Focus](https://github.com/metageni/FOCUS)
***

Usaremos el set de datos incluido en la imagen virtual **MGlinux18.2** y que se encuentra en el archivo comprimido `dietas.v2.tar.gz` en `Documents`. Para descomprimirlo ejecutemos el comando:

```bash
cd Documents
```
```bash
tar xzf dietas.v2.tar.gz
```
Tendremos una carpeta llamada Dietas en `Documents` con cinco metagenomas en la subcarpeta `seqs`.

Primero debemos activar el ambiente `conda` en donde esta instalado `Focus`:

```bash
conda activate focus
```

Vayamos a la carpeta con las secuencias que están ambos archivos `_R1` y `_R2` en las cinco subcarpetas:

```bash
cd seqs
```

Focus puede trabajar con archivos fastq pero es mejor primero **ensamblar ambos archivos** y de una vez convertirlos a formato fasta para una más ágil clasificación. Esto lo podemos hacer con `flash` fácilmente, entremos a la carpeta de una muestra `cd C08` y luego ejecutemos `flash`:

```bash
flash -t 2 C08.R1.fastq C08.R1.fastq --output-prefix=C08
```

flash nos genera varios archivos, uno con las secuencias ensambladas (`C08.extendedFrags.fastq`) y otros dos con las secuencias no ensambladas (`C08.notCombined_1.fastq C08.notCombined_2.fastq`), como todas las secuencias son útiles, **podemos unirlas** en un solo archivo con `cat`:

```bash
cat C08.extendedFrags.fastq C08.notCombined_1.fastq C08.notCombined_2.fastq > C08.fq
```

Ahora podemos **convertir este archivo a formato fasta** con `awk`:

```bash
awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' C08.fq > C08.fasta
```

Borremos los archivo innecesarios:

```bash
rm C08.extendedFrags.fastq C08.histogram C08.notCombined_2.fastq C08.hist C08.notCombined_1.fastq
```

Ahora tenemos que hacer lo mismo para las otras muestras o bien, correr un script que tenemos para procesar todas las muestras automáticamente. Para esto tenemos que estar en la carpeta `seqs` y desde allí correr el script `multiple_fastq_merger-converter`

```bash
multiple_fastq_merger-converter
```

al final tendremos una nueva carpeta llamada fastas con todos los archivos en formato fasta de las cinco muestras.

Ahora si ya podemos correr focus utilizando la carpeta fastas como entrada:

```bash
focus -q fastas -o focus -t 2
```

Focus creará una nueva carpeta de salida (`focus`) con las clasificaciones a varios niveles taxonómicos, podemos abrir estas carpetas para ver los resultados y analizarlos con alguna hoja de cálculo (Excel, gnumeric, etc.) o bien podemos convertir la hoja con todos los niveles taxonómicos para analizarla con [ampvis2](https://github.com/KasperSkytte/ampvis2) en [RStudio](https://posit.co/products/open-source/rstudio/). Esto lo podemos ejecutar con un script que tenemos para ello:
```bash
cd focus
```
```bash
focus2ampvis
```

Este script nos generará un archivo (`OTU_table.tsv`) con el formato adecuado para analizarlo con `ampvis2`, el proceso es igual que si fuera una análisis 16S, ver el proceso [aquí](https://sites.google.com/ciad.mx/cursobioinfomicrob/visualización?authuser=0).
***
