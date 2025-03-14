# Limpieza de secuencias 16S
El primer paso para el análisis de secuencias de amplicones del gen ribosomal 16S o 18S es la limpieza de las secuencias. Usaremos un set de datos de 16S región V4 obtenidos de heces de ratas.

#### Archivos y programas
Los programas y set de datos ya se encuentran instalados en la imagen virtual **MGlinux18.2**.

- [MiSeqSOP.tar.gz](https://drive.google.com/file/d/1qmUd0AYSAjj2ND15Bfi2koTAU5GWRk1n/view?usp=drive_link). Set de datos
- [fastp](https://github.com/OpenGene/fastp)
- [flash](https://github.com/ebiggers/flash)

***
### Procedimiento

Descomprimir el archivo MiSeqSOP.tar.gz

```bash
tar xzf MiSeqSOP.tar.gz
```
#### Limpieza de secuencias pareadas

Se obtendrán 19 carpetas con un par se secuencias pareadas cada una. Entrar a la primer carpeta `F3D1` y realizar la limpieza de las secuencias:

```bash
fastp -w 2 -q 30 -i F3D1_S189_L001_R1_001.fastq.gz -I F3D1_S189_L001_R2_001.fastq.gz -o F3D1.R1.fastq.gz -O F3D1.R2.fastq.gz
```
Estamos usando dos núcleos (`-w 2`) y eliminando secuencias menores a una calidad de Q30 (`-q 30`) y salvándolas a dos archivos.

#### Ensamble de las secuencias

Ahora podemos ensamblar las secuencias pareadas con `flash`.

```bash
flash -t 2 -m 8 -o F3D1 F3D1.R1.fastq.gz F3D1.R2.fastq.gz
```
Aquí estamos usando 2 núcleos (`-t 2`) para ensamblar secuencias pareadas ya limpias y que tengan un traslape mínimo de 8 bases (`-m 8`) y tengan como nombre F3D1 (`-o F3D1`). Las secuencias se salvarán en el directorio actual.

#### Conversión a fasta

Por último podemos convertir las secuencias ensambladas a formato fasta:

```bash
sed '1~4s/^@/>/p;2~4p' F3D1.extendedFrags.fastq > F3D1.fasta
```
***
Este proceso habría que hacerlo para todas las secuencias en todas las carpetas, lo que sería muy tedioso. Podemos usar un simple script con [chatGPT](https://chatgpt.com/) para que procese todas las secuencias automáticamente usando el siguiente prompt:

`Crea un script muy sencillo en bash que entre a subdirectorios, obtenga el basename de cada subdirectorio, ejecute una limpieza de las dos secuencias fastq.gz con fastp (-w 2 -q 30) y luego las ensamble con flash (-t 2 -m 8) y por último convierta de fastq a fasta. Las secuencias originales tienen por ejemplo, para subdirectorio F3D1, el formato  F3D1_S189_L001_R1_001.fastq.gz F3D1_S189_L001_R2_001.fastq.gz y todas las secuencias de salida lleven el basename obtenido del directorio (F3D1 para el ejemplo).`

ChatGPT nos creará un script como el siguiente:

```bash
#!/bin/bash

# Verifica que fastp y flash estén instalados
command -v fastp >/dev/null 2>&1 || { echo >&2 "fastp no está instalado. Saliendo."; exit 1; }
command -v flash >/dev/null 2>&1 || { echo >&2 "flash no está instalado. Saliendo."; exit 1; }

# Iterar sobre cada subdirectorio
for dir in */; do
    # Obtener el nombre base del subdirectorio
    basename=$(basename "$dir")

    # Definir nombres de archivos
    read1="$dir/${basename}_S*_L001_R1_001.fastq.gz"
    read2="$dir/${basename}_S*_L001_R2_001.fastq.gz"

    # Verificar si los archivos existen
    if ls $read1 1> /dev/null 2>&1 && ls $read2 1> /dev/null 2>&1; then
        echo "Procesando: $basename"

        # Limpieza con fastp
        fastp -i $read1 -I $read2 -o "$dir/${basename}_clean_R1.fastq" -O "$dir/${basename}_clean_R2.fastq" -w 2 -q 30

        # Ensamblado con flash
        flash "$dir/${basename}_clean_R1.fastq" "$dir/${basename}_clean_R2.fastq" -t 2 -m 8 -d "$dir"

        # Convertir fastq a fasta
        awk 'NR%4==1 {print ">" substr($0,2)} NR%4==2 {print}' "$dir/out.extendedFrags.fastq" > "$dir/${basename}.fasta"

        echo "Finalizado: $basename"
    else
        echo "Archivos no encontrados en $dir, omitiendo."
    fi
done
```
Copiémoslo y salvémoslo en nuestra carpeta para ejecutarlo; podemos llamarlo `fastp_cleaner.sh` y para hacerlo ejecutable:

```bash
chmod +x fastp_cleaner.sh
```
Ahora si podemos correrlo solo llamándolo desde la terminal.

```bash
./fastp_cleaner.sh
```
***

Una vez limpias las secuencias podemos proceder a buscar y eliminar [secuencias quiméricas](Quimeras.md).
