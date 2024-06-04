# Ensamble de un genoma Illumina

Se tienen lecturas WGS de la cepa K12 de *E. coli* secuenciada con la plataforma Illumina (2x150 pb, 388,895 secuencias), estas lecturas se tienen que limpiar y realizar un ensamble con ellas.
***
#### Archivos y programas
Los programas y set de datos ya se encuentran instalados en la imagen virtual **MGlinux18.2**.
- [ecoli-K12_Illumina.tar.gz](https://drive.google.com/file/d/1NOcflmwa6ioLDOjFVhhl5TdhJbIgpBO1/view?usp=sharing)
- [SPAdes](https://cab.spbu.ru/software/spades/)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- [Quast](https://cab.spbu.ru/software/quast/)
***
### Procedimiento
- Limpieza de lecturas
- Ensamble denovo con SPAdes
- Validación

Los datos para estos análisis de pueden bajar de aquí: [ecoli-K12_Illumina.tar.gz](https://drive.google.com/file/d/1NOcflmwa6ioLDOjFVhhl5TdhJbIgpBO1/view?usp=sharing) este archivo comprimido tiene tres archivos, `SSR6436961_1.fastq`, `SSR6436961_2.fastq` y `metadata.txt`.
***

### Limpieza
Revisar la calidad de las secuencias con `FastQC` como se muestra en la página de [limpieza de secuencias Illumina](https://bioinformatica.ciad.mx/home/preparaci%C3%B3n-secuencias/limpieza-de-lecturas/illumina)

Con los datos obtenidos del análisis con `FastQC`, hacer es limpiar las secuencias con [cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html#):

```bash
cutadapt -u 15 -U 15 -a CAAGCAGAAGACGGCATACGAGAT -a CTGTCTCTTATACACATCT -A AATGATACGGCGACCACCGAGATCTACAC -A CTGTCTCTTATACACATCT --times 2 -q 30,30 --trim-n -o ecoli.R1.fastq -p ecoli.R2.fastq --json=ecoli.cutadapt.json SRR6436961_1.fastq SRR6436961_2.fastq
```
Qué le estamos pidiendo a `cutadapt` que haga? Veamos:
- `-u 15` le decimos que corte 15 bases al principio de cada secuencia.
- `-U 15` ahora que elimine 15 bases al final de cada secuencia.
- Con `-a -A` le pedimos que busque los adaptadores (secuencias) y los borre o parte de ellos. `-a` es para el archivo R1 y `-A` para el R2.
- `--times 2` busque los adaptadores dos veces para asegurarnos eliminarlos.
- `-q 30,30` Eliminar las bases al final y principio de la secuencia que tengan una calidad inferior a Q30 (1/1,000).
- `--trim-n` Eliminar las bases ambiguas: N.
- `-o`y `-p` son los nombres con los que queremos nombrar a los archivos de salida.
- `--json` creación de un reporte en formato `.json` que puede ser leído por `MultiQC`.
- Y por último, se ponen los archivos de entrada, los que queremos limpiar.

Este proceso genera dos archivos ya limpios (`ecoli.R1.fastq` y `ecoli.R2.fastq`), podemos volver a correr `FastQC` para verificar la calidad de las secuencias.
***
### Ensamble de las secuencias
Existen varios programas para ensamblar secuencias en contigs, algunos son más adecuados para cierto tipo de plataforma de secuenciación y otros se pueden adaptar a casi cualquier tipo. Dos cosas son importantes para una buena cobertura, **suficiente cantidad de lecturas** (secuencias) y que sean de **buena calidad**. La primera nos la da generalmente el equipo de secuenciación pero la segunda no siempre.

#### Spades
Probaremos ensamblar los archivos limpios con [SPAdes](https://cab.spbu.ru/software/spades/). Debido a que este programa consume muchos recursos y tarda tiempo en ejecutarse se puede correr en el servidor Biobacter o bien en la imagen virtual **MGlinux18.2** pero obviamente tardará más.

```bash
spades -k 21,33,55,77 --careful -1 ecoli.R1.fastq -2 ecoli.R2.fastq -o spades
```
El resultados aparece en el directorio spades como `contigs.fasta` y como `scaffolds.fasta`. El primero son los contigs y el segundo un *scaffold* del primero.

Veamos un resumen de estos archivos; entremos al directorio donde `spades` creó los ensambles:
```bash
cd spades
```
y veamos un resumen y estadísticas de los contigs con el script `basic-stats`:
```bash
basic-stats contigs.fasta
```
Veamos ahora el scaffolding:
```bash
basic-stats scaffolds.fasta
```
***
### Validación del ensamble
El o los ensambles obtenidos podemos compararlos y evaluarlos con [Quast](https://cab.spbu.ru/software/quast/);
hagamos la evaluación:
```bash
quast contigs.fasta scaffolds.fasta
```
Quast genera los resultados de la evaluación en un nuevo directorio (`quast_results`) y dentro del cual hay un archivo `html` para ser visualizado por cualquier navegador.
***
### Anotación del ensamble
Como ultimo paso del proceso, podemos anotar el genoma obtenido del mejor ensamble realizado anteriormente. Este proceso implica obtener CDS y RNAs y compararlos con bases de datos por lo que no podemos hacerlo fácil y rápidamente en una imagen virtual.
Si podemos hacerlo en el servido *Biobacter* que tiene el software necesario y las bases de datos instaladas.

El programa que usaremos es [PROKKA](https://github.com/tseemann/prokka) y, como para muchos programas, el nombre de las secuencias (headers) debe ser corto y sencillo. Para recortar los nombre de las secuencias usemos `awk`:
```bash
awk '/^>/{print ">'K12'_"++i; next}{print}' scaffolds.fasta > K12.fasta
```
Ahora podemos subir (con `scp`) el ensamble renombrado a nuestro `home` en el servidor y allí correr la anotación:
```bash
scp K12.fasta usuario@187.141.151.196:
```
**NOTA**: cambiar `usuario` por el login name de cada uno.

Entrar al servidor *biobacter*:
```bash
ssh usuario@187.141.151.196
```

Ya en el servidor *Biobacter* podremos correr PROKKA pero antes debemos activar el ambiente `conda` apropiado:
```bash
conda activate prokka
```

```bash
prokka --genus Escherichia --species coli --strain K12 --cpus 4 K12.fasta
```
PROKKA genera una carpeta con archivos en diferentes formatos.
***
