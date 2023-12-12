# Pangenoma con Anvi'o

Anvi'o también puede realizar filogenómica y pangenomas de una manera sencilla pero muy potente. Una descripción más detallada se encuentra en [Phylogenomics](https://merenlab.org/2017/06/07/phylogenomics/).

Primero tenemos que tener obviamente los genomas a analizar en formato fasta, usaremos el set de datos que contiene tres genomas de *E. coli*; los cuales los podemos decargar de [aquí](https://drive.google.com/file/d/1OSoJIfb7kkdGx4rJrfHHKYqy-d_Ucejw/view?usp=share_link), descomprimir y usar los archivos `.fasta`.

### Preparación de archivos
Es importante que los archivos fasta cumplan los siguientes requisitos:
- La terminación del archivo **debe ser** `.fa`.
- Los nombre **NO deben empezar con un número**.
- De ser necesario, solo pueden llevar guion bajo en el nombre, **no puntos ni guión medio**.
- Los nombres de las secuencias deben ser cortos, sin espacios y solo con caracteres alfanuméricos y cuando mucho guión bajo.
- Las bases en las secuencias **solo pueden ser A, C, G, T y N**.

### Análisis
Generar una base de datos para cada genoma, como tenemos varios genomas a procesar podemos hacer un `for loop`. El script `anvi-script-FASTA-to-contigs-db` checa los archivos fasta para ver si cumplen las características antes señaladas y si no los reformatea (`anvi-script-reformat-fasta`), crea la base de datos (`anvi-gen-contigs-database`) y busca genes housekeeping (`anvi-run-hmms`).

Este proceso puede durar horas y por lo tanto no sería raro que la terminal se desconectara antes de terminar, para evitar esto hay varias opciones, una buena es usar el comando de linux `screen`.

Primero hay que crear una nueva "sesión" con screen que la llamaremos `contigs2dbs`, pero puede ser cualquier nombre útil:

```bash
screen -S contigs2dbs
```

Una vez que se activa esta sesión podemos ya ejecutar el comando anterior pero teniendo cuidado de activar anvio primero:

```bash
conda activate anvio-8
```
#### Creación base de datos

```bash
for f in *.fa; do anvi-script-FASTA-to-contigs-db $f; done
```
Una vez corriendo podemos "desconectarnos" de la sesión de **screen** y volver a la sesión de terminal anterior presionando al mismo tiempo las teclas `Ctrl a d` :

```bash
Ctrl a d
```
Para volverse a conectar a la sesión y ver como va el análisis:

```bash
screen -r contigs2dbs
```
Si solo hemos hecho una sesión entonces no hace falta poner el nombre de la sesión, solo `screen -r`.

Al terminar el paso anterior, que puede durar horas, habrá que generar un listado de los genomas para usar con los siguientes pasos de anvio.
***
#### Listado de genomas

Generar un archivo con los nombres de los genomas y la ruta a los archivos de la base de datos a la que pertenecen

```bash
ls -1 *.fa > names.tmp
```
```bash
ls -1 *.db > dbs.tmp
```
```bash
paste names.tmp dbs.tmp > genome.list
```
```bash
rm *.tmp
```
```bash
sed -i 's/.fa//g' genome.list
```
```bash
sed -i '1i name\tcontigs_db_path' genome.list
```

El archivo de salida (`genome.list`) debe tener la siguiente estructura:

| name | contigs_db_path |
| --- | --- |
| CP001368_0157 | CP001368_0157.db |
| DH10B | DH10B.db |
| DH1 | DH1.db |

La primer columna (`name`) podemos cambiarla al nombre que deseemos, pero la segunda columna no pues es la ruta al archivo. Para los nombres de la primer columna **NO usar nombres que empiecen con número!** La separación entre columnas es con tabuladores.
***
#### Funciones
**Nota.** Esta opción solo esta disponible en el servidor **Biobacter** y no en la imagen virtual **MGlinux18.2**

Opcionalmente podemos añadirle funciones a los genes de cada genoma, para lo cual usamos la base de datos COGS de NCBI, además de otras anotaciones.

```bash
for g in *.db; do \
anvi-run-hmms -c $g --num-threads 4 \
anvi-run-ncbi-cogs -c $g --num-threads 4 \
anvi-scan-trnas -c $g --num-threads 4 \
anvi-run-scg-taxonomy -c $g --num-threads 4 \
done
```
***
#### Single Copy Genes
Obtener los *Single Copy Genes* (SCG) para el análisis filogenético, checar si anvio está activado aún, si no, volver a activarlo (`conda activate anvio-8`).

```bash
anvi-get-sequences-for-hmm-hits --external-genomes genome.list -o concatenated-proteins.fa --hmm-source Bacteria_71 --return-best-hit --get-aa-sequences --concatenate
```
#### Arbol filogenético
Crear un árbol filogenético con los SCG

```bash
anvi-gen-phylogenomic-tree -f concatenated-proteins.fa -o phylogenomic-tree.txt
```
#### Pangenoma
Para obtener el pangenoma:
```bash
anvi-gen-genomes-storage -e genome.list -o PANGENOME-GENOMES.db
```
```bash
anvi-pan-genome -g PANGENOME-GENOMES.db -n PANGENOME -T 2
```
Si se colaron archivos con nombre no aptos, ver arriba, aquí es donde se botará el proceso con un error. El segundo comando puede tardar bastante en completarse.

Opcionalmente (no disponible en la imagen virtual) podemos hacer un análisis de **Average Nucleotide Identitity** (ANI)

```bash
anvi-compute-genome-similarity -e genome.list --program pyANI -o ANI -p PANGENOME/PANGENOME-PAN.db -T 8
```
#### Visualización
Para visualizar el **árbol filogenético**:
```bash
anvi-interactive -t phylogenomic-tree.txt -p profile.db --title 'tree 1' --manual
```

El árbol se puede visualizar en cualquier programa que lea dendrogramas en formato newick

Para visualizar el **pangenoma**:
```bash
anvi-display-pan -p PANGENOME/PANGENOME-PAN.db -g PANGENOME-GENOMES.db
```
***
