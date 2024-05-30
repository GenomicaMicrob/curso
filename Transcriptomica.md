# Análisis transcriptómico de una cepa bacteriana
Haremos un análisis transcriptómico de la expresión de una cepa bacteriana de *Vibrio parahaemolyticus* cultivada en dos tipos de medios, uno rico en nutrientes (TSB) y otro mínimo (MM9).

#### Archivo y scripts necesarios para este ejercicio

- transcritos.tar.gz (176.3 MB)
- [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) v2.5.2
- [samtools](https://www.htslib.org/) 1.18
- [FADU](https://github.com/IGS/FADU) v1.9.0

Primero descarguemos el **set de datos** del servidor **Biobacter** a la carpeta de `Documents`, entremos a ella y descarguemos:

```bash
scp usuario@187.141.151.196:/raid1/datasets/transcritos.tar.gz .
```
decomprimamos el archivos

```bash
tar xzf transcritos.tar.gz
```

El set de datos tiene secuencias de transcritos (ADN) de tres réplicas para cada uno de los medios TSB y MM9 (recortados a un máximo de 500,000 secuencias por archivo pareado); en total 12 archivos tipo fasta pareado todos en subcarpetas. Éstos se mapearán sólo al cromosoma II (por simplicidad) de la cepa M0904 de *Vibrio parahaemolyticus*.
En el archivo comprimido también viene el archivo fasta del cromosoma II (`M0904_ChII.fasta`) y un archivo gff con las anotaciones de las CDS encontradas en ese cromosoma (`M0904_ChII.gff`).
***
### Creación de genoma de referencia

Primero debemos crear una base de datos del cromosoma II a partir de la secuencia completa del mismo, esto lo haremos con bowtie2.

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
Repetir lo mismo para cada muestra.

***
### Análisis de transcritos
**Importante!** Si no tenemos instalado FADU, por favor seguir las instrucciones para su instalación [aquí](/Users/bruno/Documents/GitHub/curso/FADU.md).
***
Ahora ya podremos hacer el análisis de los transcritos con FADU, para esto necesitamos el archivo GFF que tiene esta estructura; más info del formato [aquí](https://www.biobam.com/differences-between-gtf-and-gff-files-in-genomic-data-analysis/#:~:text=The%20General%20Feature%20Format%20(GFF,Sequence%20Ontology%20Project%20(v3).

```
M0904_ChII	Geneious	CDS	1151	3121	.	+	0	ID=M0904_ChII-1;Name=DUF3346 domain-containing protein
M0904_ChII	Geneious	CDS	3410	3877	.	-	0	ID=M0904_ChII-2;Name=hypothetical protein
M0904_ChII	Geneious	CDS	3410	3877	.	-	0	ID=M0904_ChII-3;Name=transcriptional regulator
```
***
Corramos el análisis con la **primer muestra**:

```bash
/opt/FADU-1.9.0/fadu.jl -M -p -g ../M0904_ChII.gff -b MM9r1.sorted.bam -o ../ -f "CDS" -a "ID"
```
```
-M, --remove_multimapped If enabled, remove any reads or fragments that
                        are mapped to multiple regions of the genome
-p, --keep_only_proper_pairs
-g, --gff3_file /path/to/annotation.gff3
-b, --bam_file /path/to/file.bam
-f, --feature_type FEATURE_TYPE Which GFF3 feature type (column 3)
-a, --attribute_type ATTRIBUTE_TYPE   Which GFF3 feature type (column 9)
```
El resultado de `FADU` se encuentra en un archivo delimitado por tabuladores en el directorio superior llamado `MM9r1.sorted.counts.txt` y contiene cinco columnas:

1. El CDS (featureID)
2. Número de bases que solo alinean con la CDS (uniq_len)
3. El número de alineamientos (num_alignments)
4. El conteo de transcritos (counts)
5. Los transcritos por cada millón de kilobases (tpm)

```
featureID       uniq_len num_alignments counts   tpm
M0904_ChII-1    1971     7.0	          6.66   407.96
M0904_ChII-10   210      9.0	          5.95   3416.15
M0904_ChII-100  1137     0.0	          0.00   0.00
```
Para análisis posteriores solo necesitaremos el número de transcritos y los tpm:
```bash
cut -f1,4 ../MM9r1.sorted.counts.txt | sed "s/counts/MM9r1/" > ../MM9r1.counts
cut -f1,5 ../MM9r1.sorted.counts.txt | sed "s/tpm/MM9r1/" > ../MM9r1.tpm
```
Ya que el programa para calcular los genes diferencialmente expresados no acepta decimales, redondeemos el archivo `counts`:
```bash
sed -i 's/\b0\.0\b/0/g; s/\b0\.00\b/0/g' ../MM9r1.counts
```

Los archivos intermedios los podemos borrar:
```bash
rm *.sam *.bam *.bai
```

#### Hagamos lo mismo para cada una de las muestras.
O bien, podemos crear un sencillo script para que entre en cada subdirectorio y correr FADU, ver este ejemplo: **[fadu_script.sh](fadu_script.sh)**.
***

### Compilación de datos
Teniendo todos los archivos de salida de FADU para cada muestra, podemos unirlos (compilarlos); vayamos a la carpeta superior donde están los archivos `*_counts.tsv`

```bash
awk_compiler *.counts > counts.temp
```
**[awk_compiler](awk_compiler.md)** es un script para compilar tablas de datos con la misma información pero diferentes datos.

Al terminar de compilar, tendremos probablemente genes con cero transcritos, por lo que debemos borrar las líneas que tengan solo ceros.

```bash
head -1 counts.temp > header.temp # obtener solo la primer línea del archivos
```
```bash
sed 's/ /_/g' counts.temp | awk 'NR > 1{s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' | sed 's/_/ /g' > clean_counts.temp # elimina las lineas que tengan puros ceros
```
```bash
cat header.temp clean_counts.temp | sed 's/\t/,/g' > counts.csv # une los archivos
```
Podemos borrar los temporales
```bash
rm *.temp
```
Podemos ahora proceder a calcular con el archivo `counts.csv` los genes diferencialmente expresados con [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) en RStudio siguiendo [esta guía](DESeq2.md).
***
