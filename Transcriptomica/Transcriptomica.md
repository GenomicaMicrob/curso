# Análisis transcriptómico de una cepa bacteriana
Haremos un análisis transcriptómico de la expresión de una cepa bacteriana de *Vibrio parahaemolyticus* cultivada en dos tipos de medios, uno rico en nutrientes (TSB) y otro mínimo (MM9).

#### Archivo y scripts necesarios para este ejercicio

- [transcritos.tar.gz (176.3 MB)](https://drive.google.com/file/d/1hHg_qHgwy7fPg5aPLimKPfNfZpRiBG85/view?usp=sharing)
- [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) v2.5.2
- [samtools](https://www.htslib.org/) 1.18
- [FADU](https://github.com/IGS/FADU) v1.9.0


El set de datos tiene secuencias de transcritos (ADN) de tres réplicas para cada uno de los medios TSB y MM9 (recortados a un máximo de 500,000 secuencias por archivo pareado); en total 12 archivos tipo fasta pareado todos en subcarpetas. Éstos se mapearán sólo al cromosoma II (por simplicidad) de la cepa M0904 de *Vibrio parahaemolyticus*.
En el archivo comprimido también viene el archivo fasta del cromosoma II (`M0904_ChII.fasta`) y un archivo gff con las anotaciones de las CDS encontradas en ese cromosoma (`M0904_ChII.gff`).
***

### Creación de genoma de referencia

Primero debemos crear una base de datos del cromosoma II a partir de la secuencia completa del mismo, esto lo haremos con `bowtie2`.

Entremos a la carpeta `transcritos` y empecemos.

```bash
bowtie2-build M0904_ChII.fasta refgenome
```
Tendremos así varios archivos con terminación `.bt2`.
***
### Mapeo de muestras a la referencia

Ahora podremos mapear los archivos fasta de cada una de las muestras al genoma con esta base de datos, tenemos que entrar a cada uno de los subdirectorio para **mapear cada muestra**, por ejemplo, para la muestra `MM9r1`:

```bash
cd MM9r1
```

```bash
bowtie2 --no-unal --threads 4 -x ../refgenome -f -1 MM9r1.R1.fasta -2 MM9r1.R2.fasta -S MM9r1.sam
```
Ahora debemos hacer unas conversiones de formato con samtools:
```bash
samtools view --threads 4 -bSo MM9r1.bam MM9r1.sam
```
```bash
samtools sort --threads 4 MM9r1.bam -o MM9r1.sorted.bam
```
```bash
samtools index MM9r1.sorted.bam
```
Repetir lo mismo para cada muestra, o bien correr el script [mapper4FADU.sh](scripts/mapper4FADU.sh)

***
### Análisis de transcritos
**Importante!** Si no tenemos instalado FADU, por favor seguir las instrucciones para su instalación [aquí](Instalación_FADU.md).
***
Ahora ya podremos hacer el análisis de los transcritos con FADU, para esto necesitamos el archivo GFF que tiene esta estructura; más info del formato [aquí](https://www.biobam.com/differences-between-gtf-and-gff-files-in-genomic-data-analysis/#:~:text=The%20General%20Feature%20Format%20(GFF,Sequence%20Ontology%20Project%20(v3))).

```
M0904_ChII	Geneious	CDS	1151	3121	.	+	0	ID=M0904_ChII-1;Name=DUF3346 domain-containing protein
M0904_ChII	Geneious	CDS	3410	3877	.	-	0	ID=M0904_ChII-2;Name=hypothetical protein
M0904_ChII	Geneious	CDS	3410	3877	.	-	0	ID=M0904_ChII-3;Name=transcriptional regulator
```
***
Corramos el análisis con la **primer muestra**, tener en cuenta en donde estamos situados, idealmente en la carpeta *transcritos*.

```bash
/opt/FADU-1.9.0/fadu.jl -M -p -g M0904_ChII.gff -b MM9r1.sorted.bam -o . -f "CDS" -a "ID"
```
**Nota**. Revisar si FADU está en `/opt/FADU-1.9.0/fadu.jl` o en `/opt/FADU-1.9.1/fadu.jl` y correr el comando anterior acorde a la versión.

```
-M, --remove_multimapped If enabled, remove any reads or fragments that
                        are mapped to multiple regions of the genome
-p, --keep_only_proper_pairs
-g, --gff3_file /path/to/annotation.gff3
-b, --bam_file /path/to/file.bam
-f, --feature_type FEATURE_TYPE Which GFF3 feature type (column 3)
-a, --attribute_type ATTRIBUTE_TYPE   Which GFF3 feature type (column 9)
```
El resultado de `FADU` se encuentra en un archivo delimitado por tabuladores en el **directorio superior** llamado `MM9r1.sorted.counts.txt`

```bash
cd ..
```
El archivo contiene 5 columnas:

```
featureID       uniq_len num_alignments counts   tpm
M0904_ChII-1    1971     7.0	          6.66   407.96
M0904_ChII-10   210      9.0	          5.95   3416.15
M0904_ChII-100  1137     0.0	          0.00   0.00
```

1. El CDS (featureID)
2. Número de bases que solo alinean con la CDS (uniq_len)
3. El número de alineamientos (num_alignments)
4. El conteo de transcritos (counts)
5. Los transcritos por cada millón de bases (tpm)

Para análisis posteriores solo necesitaremos el número de transcritos y los tpm en dos archivos diferentes, por lo que procederemos a extraer los datos del archivo MM9r1.sorted.counts.txt:
```bash
cut -f1,4 MM9r1.sorted.counts.txt | sed "s/counts/MM9r1/" > MM9r1.counts
```
Y ahora para los TPMs:
```bash
cut -f1,5 MM9r1.sorted.counts.txt | sed "s/tpm/MM9r1/" > MM9r1.tpm
```
Los archivos intermedios los podemos borrar:
```bash
rm *.sam *.bam *.bai
```
#### Hagamos lo mismo para cada una de las muestras.
O bien, podemos crear un sencillo script para que entre en cada subdirectorio y correr FADU, ver este ejemplo: **[fadu_script.sh](scripts/fadu_script.sh)**.

***
**Nota.** Como se ve en el resultado, la columna `featureID` contiene solo un número consecutivo del gen pero no información de quién es. Podemos sustituir esa primer columna con el nombre del gen. Esta información está en el archivo `.gff` en la última columna (ver arriba). Por ejemplo, para el primer gen (primera linea del `.gff`) `ID=M0904_ChII-1;Name=DUF3346 domain-containing protein` sustituiremos `M0904_ChII-1` por `1-DUF3346 domain-containing protein` eliminando de paso `Name=` que ya no lo necesitamos.

Los featuresID de las 3 primeras líneas del archivo quedarían así:

| Original | Nuevo nombre |
| --- | --- |
|ID=M0904_ChII-1;Name=DUF3346 domain-containing protein|ID=1-DUF3346 domain-containing protein|
|ID=M0904_ChII-2;Name=hypothetical protein|ID=2-hypothetical protein|
|ID=M0904_ChII-3;Name=transcriptional regulator|ID=3-transcriptional regulator|

Es importante mantener un numero consecutivo (`1-`) antes del nombre por si ha varios genes con el mismo nombre (p.ej. *hypothetical protein*).

Esto lo podemos hacer automáticamente con el script en python [gff_ID_renaming.py](scripts/gff_ID_renaming.py) que nos generará un nuevo archivo con los nombre ya correctos.

```bash
gff_ID_renaming.py M0904_ChII.gff
```
El archivo `M0904_ChII_mod.gff` generado lo podemos usar al principio con FADU para generar todo el análisis de nuevo (sorry) usando el script **[fadu_script.sh](scripts/fadu_script.sh)**:
```bash
fadu_script.sh M0904_ChII_mod.gff
```
***

### Compilación de datos
Teniendo todos los archivos de salida de FADU para cada muestra, podemos unirlos (compilarlos); vayamos a la carpeta superior donde están los archivos `*_counts`

```bash
awk_compiler *.counts > counts_raw.temp
```
**[awk_compiler](scripts/awk_compiler.md)** es un script para compilar tablas de datos con la misma información pero diferentes datos.

#### Redondeo de decimales
Para los siguientes análisis con DESeq2 debemos redondear los números:
```bash
awk 'NR==1 {print; next} {printf "%s", $1; for(i=2;i<=NF;i++) printf " %d", int($i+0.5); print ""}' counts_raw.temp > counts.temp
```
#### Elimnación de CDS no transcritas
Probablemente tendremos genes con cero transcritos en todas las muestras, por lo que debemos borrar las líneas que tengan solo ceros.

1. obtener solo la primer línea del archivos:
```bash
head -1 counts.temp > header.temp
```
2. eliminar lineas con solo ceros:
```bash
sed 's/ /_/g' counts.temp | awk 'NR > 1{s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' | sed 's/_/ /g' > clean_counts.temp
```
3. Unir headers con tabla y convertirla a delimitada por comas (.csv):
```bash
cat header.temp clean_counts.temp | sed 's/\t/,/g' > counts.csv
```

Podemos borrar los temporales
```bash
rm *.temp
```
***
Podemos ahora proceder a calcular con el archivo `counts.csv` los genes diferencialmente expresados con [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) en RStudio siguiendo [esta guía](DESeq2.md).
***
