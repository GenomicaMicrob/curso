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

Ahora podremos mapear los archivos fasta de cada una de las muestras al genoma con esta base de datos, tenemos que entrar a cada uno de los subdirectorio para mapear cada muestra, por ejemplo, para la muestra `MM9r1`:

```bash
cd fastas/MM9r1
```

```bash
bowtie2 --no-unal --threads 4 -x ../../refgenome -f -1 MM9r1.R1.fasta -2 MM9r1.R2.fasta -S MM9r1.sam
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
***
### Análisis de transcritos
**Importante!** Si no tenemos instalado FADU, por favor seguir las instrucciones para su instalación [aquí](/Users/bruno/Documents/GitHub/curso/FADU.md).
***
Ahora ya podremos hacer el análisis de los transcritos con FADU:
```bash
/opt/FADU-1.9.0/fadu.jl -M -p -g ../../M0904_ChII.gff -b MM9r1.sorted.bam -o ../ -f "CDS" -a "ID"
```
El resultado de `FADU` se encuentra en un archivo delimitado por tabuladores en el directorio superior llamado `MM9r1.sorted.counts.txt` y contiene cinco columnas:

1. El CDS (featureID)
2. Número de bases que solo alinean con la CDS (uniq_len)
3. El número de alineamientos (num_alignments)
4. El conteo de transcritos (counts)
5. Los transcritos por cada millón de kilobases (tpm)

```
featureID	uniq_len	num_alignments	counts	tpm
M0904_ChII-1	1971	7.0	6.66	407.96
M0904_ChII-10	210	9.0	5.95	3416.15
M0904_ChII-100	1137	0.0	0.00	0.00
```
Los archivos intermedios los podemos borrar:
```bash
rm *.sam *.bam *.bai
```

#### Hagamos lo mismo para cada una de las muestras.
***
