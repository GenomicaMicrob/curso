# SARS-CoV2

Ensamble del genoma de una muestra de SARS-CoV2. Usaremos las secuencias de una muestra de la variante *alpha* para este ejercicio.

#### Entrar a la carpeta Alpha
	cd /home/user/Documents/SARSCoV2/Alpha/
#### Limpieza de secuencias crudas
Primero es necesario limpiar las secuencias crudas fastq, principalmente eliminar bases de baja calidad (`-q 30`) y aquellas menores a 20 bases (`-l 20`) y muy importante eliminar las secuencias de los primer usados, en este caso los primers se encuentran en el archivo `Artic.v3.fasta`. Esto lo haremos con `fastp`.

	fastp -i Alpha_S001_L001_R1_001.fastq.gz -o alpha.R1.fq -I Alpha_S001_L001_R2_001.fastq.gz -O alpha.R2.fq -q 30 -l 20 --low_complexity_filter --adapter_fasta /opt/COVID/Artic.v3.fasta --thread 2 -h alpha-fastp.html -R "alpha fastp trimming report"
#### Copiar referencia
Para ensamblar las secuencias necesitaremos hacerlo usando el genoma del virus original, copiémoslo al directorio de trabajo en el que estamos, notar que el comando termina con `.`

	cp /opt/COVID/REF/NC_045512.fasta .
#### Creacion de indices
Usaremos este genoma de referencia para crear un índice con `bwa`:

	bwa index NC_045512.fasta
#### Mapeo de secuencias
Ahora mapeemos las secuencias a la referencia y hagamos algunos formateos:

	bwa mem -t 4 NC_045512.fasta alpha.R1.fq alpha.R2.fq | samtools sort -O BAM - > assembly.sorted.bam
	samtools index assembly.sorted.bam
	samtools faidx NC_045512.fasta
#### Secuencia consenso
Teniendo ya mapeadas las secuencias calculemos una secuencia consenso de nuestro nuevo genoma alpha:

	samtools mpileup -d 50000 --reference NC_045512.fasta -a -Q 30 assembly.sorted.bam | ivar consensus -t 0 -m 2 -n N -p alpha-consensus.fa

Así tenemos ya nuestra secuencia del genoma de la muestra llamado `alpha-consensus.fa`
#### Cobertura
Podemos también calcular la cobertura a todo lo largo del genoma que tuvimos con las secuencias así como aquellas zonas donde hubo una baja profundidad de secuenciación.

	bedtools genomecov -bga -ibam assembly.sorted.bam | grep -w '0$\|1$\|2$\|3$' > alpha-lowcoverage.bed
	bedtools genomecov -d -split -ibam assembly.sorted.bam | cut -f2,3 | sed 's/\t/,/' | sed '1 i\position,coverage' > alpha-coverage.csv
#### Variaciones
Por último podemos obtener las variaciones de nuestra muestra contra el genoma de referencia.

	samtools mpileup -A -d 0 --reference NC_045512.fasta -B -Q 0 assembly.sorted.bam | ivar variants -p variants -q 20 -t 0.25 -r NC_045512.fasta -g /opt/COVID/REF/NC_045512.gff
#### Limpieza de archivos
	rm *.fasta.* fastp.json *.txt
Para hacer un **análisis completo online** del genoma obtenido podemos subir el archivo `alpha-consensus.fa` a [Nextclade](https://clades.nextstrain.org/) 


## Graficas
Para ver una gráfica de la cobertura podemos importar el archivo `assembly.sorted.bam` a **R**. 
Abramos `RStudio`, carguemos ggplot e importemos el archivo.

	library(ggplot2)
	df <- alpha.coverage
	ggplot(data = df, aes(position, coverage)) +
	ggtitle("alpha") + scale_y_continuous(trans = 'log10', name = "Profundidad (X)") + scale_x_continuous(name = "Genoma SARS-CoV-2") + geom_line(colour="grey50") + geom_area(fill="grey", alpha=0.6) + geom_hline(yintercept = 10, colour = "grey60", size = 0.5
