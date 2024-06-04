#!/bin/bash
# Script para clasificar transcritos en subdiretorios con FADU
# USO fadu.sh referencia.gff

REFGFF=$(echo $1) # salvar el archivo gff a una variable
for D in ./*; do
	if [ -d "$D" ]; then
		cd "$D"
		NAME=$(basename -s .R1.fasta *.R1.fasta)
		echo "
Procesando muestra $NAME"

		/opt/FADU-1.9.0/fadu.jl -M -p -g ../$REFGFF -b $NAME.sorted.bam -o ../ -f "CDS" -a "ID"
		cut -f1,4 ../$NAME.sorted.counts.txt | sed "s/counts/$NAME/" > ../$NAME.counts.tsv # cortar solo las columnas 1 y 4 (Counts)
		cut -f1,5 ../$NAME.sorted.counts.txt | sed "s/tpm/$NAME/" > ../$NAME.tpm # cortar solo las columnas 1 y 5 (TPM)
		sed -i 's/\b0\.0\b/0/g; s/\b0\.00\b/0/g' ../$NAME.counts.tsv # eliminar decimales ya que DESeq2 no lo acepata
		rm *.sam *.bam *.bai # borrar archivos ya usados

		echo "
Listo
-------------------------"
		cd ..
	fi
done
