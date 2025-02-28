#!/bin/bash
# Script para mapear secuencias
# USO mapper.sh

CONDA="/home/user/miniconda/etc/profile.d/conda.sh"
ANVIO="/home/user/miniconda/envs/anvio-8"  # conda environment
source "$CONDA"  # Shell for CONDA
conda activate "$ANVIO"  # Activation of anvio in CONDA

for D in ./*; do
	if [ -d "$D" ]; then
		cd "$D"
		NAME=$(basename -s .R1.fastq *.R1.fastq)
    R1=$(ls *.R1.fastq)
    R2=$(ls *.R2.fastq)
		echo "
Procesando muestra $NAME"
  bowtie2 --threads 2 -x ../contigs -1 "$R1" -2 "$R2" -S "$NAME.sam"
	samtools view -F 4 -bS "$NAME.sam" > "$NAME.raw"
	anvi-init-bam "$NAME.raw" -o "$NAME.bam"
mv *.bam *.bai ../
		echo "
Listo
-------------------------"
		cd ..
	fi
done
