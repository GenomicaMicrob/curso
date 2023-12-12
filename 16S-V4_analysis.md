## Análisis de amplicones 16S region V4

### Limpieza y ensamble

Lo primero que tenemos que hacer es limpiar y ensamblas las secuencias, para esto usaremos el pipeline del laboratorio [mg_pipeline](about:blank), pero antes de eso hay que descomprimir las secuencias que se encuentran en `Documents`:

```bash
tar xzf MiSeqSOP.tar.gz
```
lo que creará una nueva carpeta llamada `MiSeqSOP` con una subcarpeta con las secuencias, información del estudio y una tabla de metadatos. Para limpiar las secuencias, tenemos que entrar a la carpeta seqs y ahí ejecutar el script `pair-end_cleaner`:

```bash
cd MiSeqSOP/seqs
```

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

### Eliminación de quimeras

A las secuencias limpias y ensambladas habrá que eliminarle las secuencias quiméricas generadas probablemente durante la PCR, para esto tenemos otro script basado en [VSEARCH](https://github.com/torognes/vsearch) que utiliza un base de datos para este fin.

Hay que entrar a la carpeta generada en el paso anterior y que empieza con assembled.fecha_hora y ejecutar el script:

```bash
chimera_detector *.fna
```

Después de varios minutos (15 a 20), se creará una nueva subcarpeta llamada chimera_detector.fecha_hora con las secuencias libre de quimeras con terminación *.fasta

### Clasificación de secuencias
Por último paso, tenemos que clasificar las secuencias para asignarles una taxonomía, esto lo haremos con el script mg_classifier:

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
