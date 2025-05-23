### Clasificación de secuencias
Por último paso, tenemos que clasificar las secuencias para asignarles una taxonomía. Un programa muy útil para esto es usar [vsearch](https://github.com/torognes/vsearch) con la función `--cluster_fast`.

Entremos en la carpeta en donde hayamos obtenido las secuencias libres de quimeras y allí podemos ejecutar vsearch:
```bash
vsearch --cluster_fast Mock.fasta --id 0.97 --threads 2 --sizeout --centroids Mock.centroid --fasta_width
```

Así obtendremos un nuevo archivo (`Mock.centroid`) que tendrá las secuencias centroides (secuencias representativas de un cluster de secuencias agrupadas por tener una identidad superior al 97%). Este archivo fasta tendrá las secuencias con sus nombres que tendrá el nombre del archivo (`Mock`), un número (`_1`) y la cantidad de secuencias de ese cluster (`size=303`), p. ej.:
`Mock_1;size=303`

Teniendo ya las secuencias centroides representativas podemos clasificarlas usando una base de datos por medio de un blast:

```bash
vsearch --usearch_global Mock.centroid --threads 2 --db /opt/mg_pipeline/databases/EzBioCloud --id 0.50 --userout Mock.blast --userfields query+target+id
```

Tendremos ahora el archivo `Mock.blast` con tres columnas:
1. El nombre de la secuencia centroide
2. La taxonomía más cercana
3. El porcentaje de similitud

Tendríamos que hacer esto mismo para todas las muestras, lo que podemos hacerlo con un `for loop` o bien utilizar el script [mg_classifier](https://github.com/GenomicaMicrob/mg_classifier).
***
## mg_classifier

```bash
mg_classifier *.fasta
```

Nos saldrá un menu para seleccionar la **base de datos** a usar:

```latex
Select a database:
16S_rRNA
   1 EzBioCloud
   2 SILVA v132
18S_rRNA
   3 PR2 v4.2

   x exit
```

En la imagen virtual **MGlinux18.2**, por cuestiones de espacio, no tenemos todas las bases de datos, elijamos la **número 1**. Después de pocos minutos, tendremos una nueva subcarpeta con los resultados. El documento principal es una tabla de OTUs en la que está cuantas secuencias de cada OTU (líneas) se asignaron a cada muestra (columnas) analizada.
***

El **resultado** es una nueva carpeta con varios archivos de salida, pero el principal archivo es la tabla de OTUs `otus.tsv`. Además tenemos tres carpetas, una con cálculos de índices de diversidad (`diversity`), tablas de OTUs en diferentes formatos (`OTU_tables`) y archivos para graficar con R (`plots`).

```latex
mg_classifier.EzBioCloud.2025-05-06_13:13/
├── diversity
│   ├── alpha_indexes.tsv
│   ├── alpha.txt
│   ├── beta.bray_curtis.tree
│   ├── beta.bray_curtis.txt
│   └── rare.txt
├── mg_classifier.log
├── mg_classifier.report
├── mg_classifier.report.md
├── OTUs_summary.tsv
├── otus.tsv
├── OTU_tables
│   ├── core_families.txt
│   ├── core_phyla.txt
│   ├── family.tsv
│   ├── genus.tsv
│   ├── otus.xls
│   ├── phyloseq
│   │   ├── otus_table.tsv
│   │   └── taxonomy_table.tsv
│   ├── phylum.tsv
│   ├── qiime_ampvis2
│   │   ├── otus.txt
│   │   └── OTU_table.tsv
│   ├── samples-tax.tsv
│   └── STAMP
│       ├── bacteria.spf
│       └── bacteria.win.spf
└── plots
    ├── families_plot.csv
    ├── phyla15_plot.csv
    └── phyla_plot.csv

6 directories, 26 files
```
***
#### Graficado
Teniendo ya las secuencias clasificadas, quizá lo primero que tenemos que hacer es graficar estos resultados. Veamos como hacer esto en la siguiente [página](https://bioinformatica.ciad.mx/home/metagenomica/visualizacion/ampvis2).
***
FIN
