## Análisis de amplicones 16S region V4

### Limpieza y ensamble

Lo primero que tenemos que hacer es limpiar y ensamblas las secuencias, para esto usaremos el pipeline del laboratorio [mg_pipeline](https://github.com/GenomicaMicrob/metagenomic_pipeline), pero antes de eso hay que descomprimir las secuencias que se encuentran en `Documents`:

```bash
tar xzf MiSeqSOP.tar.gz
```
lo que creará una nueva carpeta llamada `MiSeqSOP` con una subcarpeta con las secuencias, información del estudio y una tabla de metadatos. Para limpiar las secuencias, tenemos que entrar a la carpeta seqs y ahí ejecutar el script `pair-end_cleaner`:

```bash
cd MiSeqSOP/seqs
```

***
Ahora, dependiendo de que imagen virtual estemos usando podremos usar alguno de dos scripts para limpienza:
- MGlinux18.3 --> `pair-end_cleaner`
- MGlinux22 (Apple silicon) --> `fastp_pair-end_cleaner`

#### MGlinux18
```bash
pair-end_cleaner
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

#### MGlinux en Apple silicon (M1-M4)
Ya que no tenemos el programa PEAR para ARM64 (Apple silicon), debemos usar otro script que nos limpie las secuencias en subdirectorios:

```bash
fastp_pair-end_cleaner
```
Nos creará una nueva carpeta `cleaned...` con las secuencias limpias en el subdirectorio `fastas_cleaned`, entrando a éste subdirectorio podemos proceder con la eliminación de secuencias quiméricas.

***
### Eliminación de quimeras

A las secuencias limpias y ensambladas habrá que eliminarle las secuencias quiméricas generadas probablemente durante la PCR, para esto tenemos otro script basado en [VSEARCH](https://github.com/torognes/vsearch) que utiliza un base de datos para este fin.

Hay que entrar a la carpeta generada en el paso anterior y que empieza con `assembled.fecha_hora` y ejecutar el script:

```bash
chimera_detector *.fna
```

Después de varios minutos (15 a 20), se creará una nueva subcarpeta llamada `chimera_detector.fecha_hora` con las secuencias libre de quimeras con terminación `.fasta`
***

### Clasificación de secuencias
Por último paso, tenemos que clasificar las secuencias para asignarles una taxonomía, esto lo haremos con el script `mg_classifier`:

```bash
mg_classifier *.fasta
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

En la imagen virtual **MGlinux18.2**, por cuestiones de espacio, no tenemos todas las bases de datos, elijamos la **número 1**. Después de pocos minutos, tendremos una nueva subcarpeta con los resultados. El documento principal es una tabla de OTUs en la que está cuantas secuencias de cada OTU (líneas) se asignaron a cada muestra (columnas) analizada.

#### Ejemplo de resultados

```latex
- Main table with OTUs found per sample: otus.tsv
- Log of the bioinformatic processes done: mg_classifier.log
OTU_tables subidrectory
   - Excel-type table with OTUs order from most abundant: otus.xls
   - Table with the number of sequences per taxon: OTUs_summary.tsv
   - Table with genera found per sample: genus.tsv
   - Table with families found per sample: family.tsv
   - Table with phyla found per sample: phylum.tsv
   - File with the samples first and the taxonomy, useful to make KRONA charts: samples-tax.tsv
   - Files with core taxa: core_phylum.txt and core_family.txt
   STAMP subidrectory
     - Table of results to be analyzed with STAMP
       Linux file: bacteria.spf
       Windows file: bacteria.win.spf
   QIIME subidrectory
     - OTU table useful for qiime, ampvis2, etc.: OTU_table.tsv and otus.txt
   Phyloseq subidrectory
     - Taxonomy table: taxonomy_table.tsv
     - Table of OTUs: otus_table.tsv
Diversity subdirectory
   - Alpha diversity indexes: alpha.txt
   - Selected alpha diversity indexes: alpha_indexes.tsv
   - Beta diversity indexes: beta matrix (.txt) and trees (.tree)
   - Rarefaction analyses of samples: rare.txt
Plot subdirectory
   - File for plotting the 15 most abundant families: families-plot.csv
   - File for plotting all the phyla: phyla_plot.csv
   - File for plotting only the 15 most abundant phyla: phyla15_plot.csv
```
***
