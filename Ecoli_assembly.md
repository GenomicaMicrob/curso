# Ensamble de un genoma Illumina

Se tienen lecturas WGS de la cepa K12 de *E. coli* secuenciada con la plataforma Illumina (2x150 pb, 388,895 secuencias), estas lecturas se tienen que limpiar y realizar un ensamble con ellas.
***
#### Archivos y programas
Los programas ya se encuentran instalados en la imagen virtual **MGlinux18.2**.
- [ecoli-K12_Illumina.tar.gz](https://drive.google.com/file/d/1NOcflmwa6ioLDOjFVhhl5TdhJbIgpBO1/view?usp=sharing)
- [SPAdes](https://cab.spbu.ru/software/spades/)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- [Quast](https://cab.spbu.ru/software/quast/)
***
### Procedimiento
- Limpieza de lecturas
- Ensamble denovo con SPAdes
- Realizar un scaffolding
- Validación

Los datos para estos análisis de pueden bajar de aquí: [ecoli-K12_Illumina.tar.gz](https://drive.google.com/file/d/1NOcflmwa6ioLDOjFVhhl5TdhJbIgpBO1/view?usp=sharing) este archivo comprimido tiene tres archivos, `SSR6436961_1.fastq`, `SSR6436961_2.fastq` y `metadata.txt`.

### Limpieza
Revisar la calidad de las secuencias con `FastQC` como se muestra en la página de limpieza de secuencias Illumina

Con los datos obtenidos del análisis con `FastQC`, hacer es limpiar las secuencias con `cutadapt`:

```bash
cutadapt -u 15 -U 15 -a CAAGCAGAAGACGGCATACGAGAT -a CTGTCTCTTATACACATCT -A AATGATACGGCGACCACCGAGATCTACAC -A CTGTCTCTTATACACATCT --times 2 -q 30,30 --trim-n -o ecoli.R1.fastq -p ecoli.R2.fastq SRR6436961_1.fastq SRR6436961_2.fastq
```
Este proceso genera dos archivos ya limpios (`ecoli.R1.fastq` y `ecoli.R2.fastq`), podemos volver a correr `FastQC` para verificar la calidad de las secuencias.

### Ensamble de las secuencias
Existen varios programas para ensamblar secuencias en contigs, algunos son más adecuados para cierto tipo de plataforma de secuenciación y otros se pueden adaptar a casi cualquier tipo. Dos cosas son importantes para una buena cobertura, **suficiente cantidad de lecturas** (secuencias) y que sean de **buena calidad**. La primera nos la da generalmente el equipo de secuenciación pero la segunda no siempre.

#### Spades
Probaremos también ensamblar el mismo archivo con SPADES para comparar los resultados, link al manual de SPADES. Debido a que este programa consume muchos recursos y tarda tiempo en ejecutarse se puede correr en el servidor Biobacter o bien en la imagen virtual **MGlinux18.2** pero obviamente tardará más.

Para este programa solo necesitamos el archivo tipo fastq y ejecutar el siguiente comando:
```bash
spades -k 21,33,55,77 --careful -1 ecoli.R1.fastq -2 ecoli.R2.fastq -o spades
```
El resultados aparece en el directorio spades como `contigs.fasta` y como `scaffolds.fasta`. El primero son los contigs y el segundo un *scaffold* del primero.

Veamos un resumen de estos archivos; entremos al directorio donde `spades` creó los ensambles:
```bash
cd spades
```
y veamos un resumen y estadísticas de los contigs:
```bash
basic-stats contigs.fasta
```
***
### Validación del ensamble
El o los ensambles obtenidos podemos compararlos y evaluarlos con `Quast`;
hagamos la evaluación:
```bash
quast contigs.fasta scaffolds.fasta
```
Quast genera los resultados de la evaluación en un nuevo directorio (`quast_results`) y dentro del cual hay un archivo html para ser visualizado por cualquier navegador.
***
