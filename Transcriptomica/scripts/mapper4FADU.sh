#!/bin/bash
# Script para mapear secuencias
# USO mapper4FADU.sh
#----- help and usage -----
display_help(){
echo -e "
____________________ mapper4FADU ____________________________________

A simple script to map clean sequences in subdirectories for use later
with FADU.
Sequences must be clean and in fastq format (.R1.fastq and .R2.fastq).
A reference genome/contig must first be referenced with bowtie2-build
and named "refgenome".

USAGE: $NAME

_______________________________________________________________________________
"
} # -h is typed, display_help
if [ "$1" == "-h" ]
then
	display_help
	exit 1
fi
# ----- script -----
for D in ./*; do
	if [ -d "$D" ]; then
		cd "$D"
		NAME=$(basename -s .R1.fasta *.R1.fasta)
		R1=$(ls *.R1.fasta)
		R2=$(ls *.R2.fasta)
		echo "
Procesando muestra $NAME"
		bowtie2 --no-unal --threads 4 -x ../refgenome -f -1 $R1 -2 $R2 -S "$NAME.sam"
 		echo "
Convirtiendo archivos sam a bam"
		samtools view -F 4 -@ 4 -bS "$NAME.sam" > "$NAME.bam"
		samtools sort --threads 4 "$NAME.bam" -o $NAME.sorted.bam
		samtools index "$NAME.sorted.bam"
		mv "$NAME.sorted.bam" "$NAME.sorted.bam.bai" ../
		echo "
Listo
-------------------------"
		cd ..
	fi
done
