## Análisis de amplicones 16S region V4.

### Limpieza y ensamble

Lo primero que tenemos que hacer es limpiar y ensamblas las secuencias, para esto usaremos el pipeline del laboratorio [mg_pipeline](about:blank), pero antes de eso hay que descomprimir las secuencias que se encuentran en `Documents`:

```bash
$ tar xzf MiSeqSOP.tar.gz
```
lo que creará una nueva carpeta llamada `MiSeqSOP` con una subcarpeta con las secuencias, información del estudio y una tabla de metadatos. Para limpiar las secuencias, tenemos que entrar a la carpeta seqs y ahí ejecutar el script `pair-end_cleaner`:

```bash
$ cd MiSeqSOP/seqs
```

```bash
$ pair-end_cleaner
```

Nos saldrá un menú con la siguiente información:

```latex
--- pair-end_cleaner v1.0.1 -------------------------
Script to unzip, clean, assemble, convert, and rename
Illumina pair-end fastq files in all 19 subdirectories.

Select a region for

 16S
  1 V3
  2 V4
  3 V3-V4

 18S
  4 V9

  x exit
```

Debemos elegir la región del gen que se secuenció, para este ejemplo la **opción 2**, region V4 del gen 16S.

El script limpiará las secuencias para eliminar las de baja calidad, inferior a Q20, quitar las bases indefinidas n y eliminar las secuencias cortas con [cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html). Posteriormente ensamblará las secuencias pareadas con [PEAR](https://cme.h-its.org/exelixis/web/software/pear/doc.html) en una solo consenso y la convertirá de fastq a fasta. Creará un archivo por para cada muestra con terminación .fna en una nueva carpeta. Esto puede tomar un par de minutos.

### Eliminación de quimeras

A las secuencias limpias y ensambladas habrá que eliminarle las secuencias quiméricas generadas probablemente durante la PCR, para esto tenemos otro script basado en [VSEARCH](https://github.com/torognes/vsearch) que utiliza un base de datos para este fin.

Hay que entrar a la carpeta generada en el paso anterior y que empieza con assembled.fecha_hora y ejecutar el script:

```bash
$ chimera_detector *.fna
```

Después de varios minutos (15 a 20), se creará una nueva subcarpeta llamada chimera_detector.fecha_hora con las secuencias libre de quimeras con terminación *.fasta

### Clasificación de secuencias
Por último paso, tenemos que clasificar las secuencias para asignarles una taxonomía, esto lo haremos con el script mg_classifier:

```bash
$ mg_classifier *.fasta
```

Nos saldrá un menu para seleccionar la base de datos a usar

```latex
Select a database:
16S_rRNA
   1 EzBioCloud-LGM 2023
   2 SILVA v132
18S_rRNA
   3 PR2

   x exit
```

En la imagen virtual `MGlinux18.2`, por cuestiones de espacio, no tenemos todas las bases de datos, elijamos la **número 1**. Después de pocos minutos, tendremos una nueva subcarpeta con los resultados. El documento principal es una tabla de OTUs en la que está cuantas secuencias de cada OTU (líneas) se asignaron a cada muestra (columnas) analizada.



## FOCUS

Archivo y script necesarios para este ejercicio:

```text
  - Dietas.tar.gz
  - multiple_fastq_merger-converter
  - focus2ampvis
```

Mediante el programa [Focus](https://github.com/metageni/FOCUS) podemos hacer una **rápida clasificación taxonómica** de muestras metagenómicas. La base de datos con la que se clasifica no es muy extensa y tampoco muy actualizada pero nos permite hacer una clasificación aceptable.

Usaremos el set de datos incluido en la imagen virtual **MGlinux18.2** y que se encuentra en el archivo comprimido `Dietas.tar.gz` en `Documents`. Para descomprimirlo ejecutemos el comando:

```bash
$ cd Documents
$ tar xzf Dietas.tar.gz
```

Tendremos una carpeta llamada Dietas en Documents con cinco metagenomas en la subcarpeta seqs.

Primero debemos activar el ambiente conda en donde esta instalado Focus:

```bash
$ conda activate focus
```

Vayamos a la carpeta con las secuencias que están ambos archivos _R1 y _R2 en las cinco subcarpetas:

```bash
$ cd seqs
```

Focus puede trabajar con archivos fastq pero es mejor primero **ensamblar ambos archivos** y de una vez convertirlos a formato fasta para una más ágil clasificación. Esto lo podemos hacer con flash fácilmente, entremos a la carpeta de una muestra:

```bash
$ cd C08
$ flash -t 2 C08.R1.fastq C08.R1.fastq --output-prefix=C08
```

flash nos genera varios archivos, uno con las secuencias ensambladas (`C08.extendedFrags.fastq`) y otros dos con las secuencias no ensambladas (`C08.notCombined_1.fastq C08.notCombined_2.fastq`), como todas las secuencias son útiles, **podemos unirlas** en un solo archivo:

```bash
$ cat C08.extendedFrags.fastq C08.notCombined_1.fastq C08.notCombined_2.fastq > C08.fq
```

Ahora podemos **convertir este archivo a formato fasta**:

```bash
$ awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' C08.fq > C08.fasta
```

Borremos los archivo innecesarios ya:

```bash
$ rm C08.extendedFrags.fastq C08.histogram C08.notCombined_2.fastq C08.hist C08.notCombined_1.fastq
```

Ahora tenemos que hacer lo mismo para las otras muestras o bien, correr un script que tenemos para procesar todas las muestras automáticamente. Para esto tenemos que estar en la carpeta seqs y desde allí correr el script multiple_fastq_merger-converter

```bash
$ multiple_fastq_merger-converter
```

al final tendremos una nueva carpeta llamada fastas con todos los archivos en formato fasta de las cinco muestras.

Ahora si ya podemos correr focus utilizando la carpeta fastas como entrada:

```bash
(focus)$ focus -q fastas -o focus -t 2
```

Focus creará una nueva carpeta de salida (focus) con las clasificaciones a varios niveles taxonómicos, podemos abrir estas carpetas para ver los resultados y analizarlos con alguna hoja de cálculo (Excel, gnumeric, etc) o bien podemos convertir la hoja con todos los niveles taxonómicos para analizarla con [ampvis2](https://github.com/KasperSkytte/ampvis2) en [RStudio](https://posit.co/products/open-source/rstudio/). Esto lo podemos ejecutar con un scritp fácilmente:

```bash
(focus)$ cd focus
(focus)$ focus2ampvis
```

Este script nos generará un archivo (`OTU_table.tsv`) con el formato adecuado para analizarlo con ampvis2, el proceso es igual que si fuera una análisis 16S, ver el proceso [aquí](https://sites.google.com/ciad.mx/cursobioinfomicrob/visualización?authuser=0).

## SUPERFOCUS

También podemos **clasificar funcionalmente** con [SUPERFOCUS](https://github.com/metageni/SUPER-FOCUS), pero debido al tamaño de la base de datos y la intensidad computacional, no es factible hacerlo en una máquina virtual. Pero el proceso es igual a FOCUS y el resultado es una hoja tipo excel que podemos descargar de [aquí](https://docs.google.com/spreadsheets/d/1n8Eb1wdwlaIVkuR7mlNEe3i2_mZbXsAY/edit?usp=sharing&ouid=105167676353150462869&rtpof=true&sd=true).
