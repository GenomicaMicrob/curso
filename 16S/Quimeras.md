## Secuencias quiméricas

Durante el proceso de PCR para amplificar el gen 16S, como en todas las PCRs, suceden errores que introducen secuencias quiméricas que son una combinación de dos moléculas para formar una tercera híbrida.

Para identificar estas secuencias quiméricas, es necesario comparar las secuencias limpias contra una base de datos y aquellas que difieran grandemente, se pueden considerar quimeras.

### Procedimiento

Usaremos el programa [vsearch](https://github.com/torognes/vsearch) para comparar nuestras secuencias limpias contra la base de datos [SILVA](https://www.arb-silva.de/).

Primero necesitamos mover todas las secuencias fasta a una nueva carpeta:
```bash
mkdir fastas
```
Buscaremos todas las secuencias `.fasta` y movámoslas a la carpeta:

```bash
find . -name "*.fasta" -exec mv {} fastas/ \;
```
Entremos a la carpeta:
```bash
cd fastas
```
Ahora si podemos buscar quimeras en una de los archivos y salvémoslas a un nuevo archivo:

```bash
vsearch --uchime_ref F3D1.fasta --db reference.db --threads 4 --nonchimeras F3D1.fna
```
Tendríamos que hacer lo mismo para todas los archivos, o bien hacerlo con un `for loop`; este proceso puede tomar mucho tiempo por lo que podemos correrlo en una sesión de `screen`.

```bash
for f in *.fasta; do vsearch --uchime_ref $f --db reference.db --threads 4 --nonchimeras $f.fna
```
Podemos renombrar el nombre de los archivos de salida para simplificarlos:

```bash
rename 's/.fasta//' *.fna
```
***

Tenemos así secuencias libres de quimeras por lo que podemos proceder a [clasificarlas taxonómicamente](Clasificacion.md).
