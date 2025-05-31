##!/bin/bash
# Script para clasificar transcritos en subdiretorios con FADU
# USO fadu.sh referencia.gff

# Function to display help
function show_help {
    echo "Uso: fadu.sh referencia.gff"
    echo "Este script clasifica transcritos en subdiretorios con FADU."
    echo "Argumentos:"
    echo "  referencia.gff  Archivo GFF de referencia."
    exit 0
}

# Check if help is requested
if [[ $1 == "-h" ]]; then
    show_help
fi

# Check if the reference GFF is provided
if [[ -z $1 ]]; then
    echo "Error: No se proporcionó el archivo de referencia GFF."
    show_help
fi

REFGFF=$(echo $1) # salvar el archivo gff a una variable
for D in ./*; do
    if [ -d "$D" ]; then
        cd "$D"

        # Check for both .fasta and .fastq files
        if ls *.R1.fasta 1> /dev/null 2>&1; then
            NAME=$(basename -s .R1.fasta *.R1.fasta)
        elif ls *.R1.fastq 1> /dev/null 2>&1; then
            NAME=$(basename -s .R1.fastq *.R1.fastq)
        else
            echo "No se encontraron archivos .R1.fasta o .R1.fastq en $D"
            cd ..
            continue
        fi
				# Check if the sorted BAM file is present
        if [ ! -f "$NAME.sorted.bam" ]; then
            echo "Archivo $NAME.sorted.bam no encontrado en $D"
            cd ..
            continue
        fi

        echo "
Procesando muestra $NAME"

        /opt/FADU-1.9.0/fadu.jl -M -p -g ../$REFGFF -b $NAME.sorted.bam -o . -f "CDS" -a "ID"
        cut -f1,4 $NAME.sorted.counts.txt | sed "s/counts/$NAME/" > ../$NAME.counts # cortar solo las columnas 1 y 4 (Counts)
        cut -f1,5 $NAME.sorted.counts.txt | sed "s/tpm/$NAME/" > ../$NAME.tpm # cortar solo las columnas 1 y 5 (TPM)
        #rm *.sam *.bam *.bai # borrar archivos ya usados

        echo "
Lista muestra $NAME
-------------------------"
        cd ..
    fi
done
