# Metagenómica shotgun con Anvi'o

Haremos un análisis metageómico completo de unas muestras de heces de camarón con la plataforma [Anvi'o](https://anvio.org/).

### Ensamble de contigs
Tenemos primero que crear un listado de las muestras que usaremos para cada uno de las pareadas; para esto tenemos que tener todas los archivos en una carpeta o dentro de subdirectorios pero que los archivos pareados tengan algo que los distinga y que sea del mismo tipo siempre. Por ejemplo: `Sample01.R1.fastq` y `Sample01.R2.fastq`. 

Con los siguientes comandos podemos crear una variable para que automáticamente genere un listado para cada archivo pareado separado por comas que se encuentren en subcarpetas:

	R1=$(find . -name '*.R1.*' | tr '\n' ',' | sed s'/.$//; s/\.\///g')
	
	R2=$(find . -name '*.R2.*' | tr '\n' ',' | sed s'/.$//; s/\.\///g')

Y entonces estas variables ya las podemos usar en el comando de `megahit`:

	megahit -1 $R1 -2 $R2 --min-contig-len 1000 -m 0.8 -t 2 -o megahit
Este análisis puede tardar varios minutos dependiendo del número de núcleos que le hayamos asignado; nos generará una carpeta llamada megahit con archivos, el ensamble estará en `final.contigs.fa`, podemos inspeccionarlo y ver que tiene 2,962 contigs y una N50 = 1,614 pb. Podemos ver también que el encabezado (heading) de los contigs es así:

	grep ">" megahit/final.contigs.fa | head -3
>k99_69 flag=1 multi=3.0000 len=1393
>k99_70 flag=1 multi=2.0000 len=1036
>k99_77 flag=1 multi=3.0000 len=1587

Este tipo de encabezado no le gusta a **Anvio**, por lo que tenemos que cambiarlo a algo mas sencillo, sin espacios, corto y con números consecutivos:

	awk '/^>/{print ">'contig'_"++i; next}{print}' megahit/final.contigs.fa > contigs.fa

Ahora tendremos los encabezados mas sencillos:

	grep ">" contigs.fa | head -3

>contig_1
>contig_2
>contig_3

Importante, megahit crea muchos archivos intermedios grandes que solo nos ocupan espacio, y no tenemos mucho en la imagen virtual, borrémoslos, pero antes asegurémonos que tenemos el archivo contigs.fa (generado anteriormente) en la carpeta. Como comprobación de que tanto espacio ocupa, chequemos el tamaño con `du`:

	du -h megahit/

>280M  megahit/intermediate_contigs
>285M megahit/

Ahora si borrémoslo:

	rm -fr megahit/
### Mapeo de contigs

Teniendo un archivo de contigs formado con todos los metagenomas, debemos ahora mapear cada metagenoma a ese esos contigs, esto nos dirá cuantas y cuales secuencias de cada metagenoma contribuyeron a cada contig, primero tenemos que crear un índice:

	bowtie2-build contigs.fa contigs

Ahora si mapeamos los metagenomas usando ese índice (contigs):

	bowtie2 --threads 2 -x contigs -f sample1.fna -S sample1.sam

Esto tenemos que hacer para cada uno de las cinco metagenomas:

	bowtie2 --threads 2 -x contigs -1 C08/C08.R1.fastq -2 C08/C08.R2.fastq -S C08.sam
	bowtie2 --threads 2 -x contigs -1 P08/P08.R1.fastq -2 P08/P08.R2.fastq -S P08.sam
	bowtie2 --threads 2 -x contigs -1 P18/P18.R1.fastq -2 P18/P18.R2.fastq -S P18.sam
	bowtie2 --threads 2 -x contigs -1 S04/S04.R1.fastq -2 S04/S04.R2.fastq -S S04.sam
	bowtie2 --threads 2 -x contigs -1 S19/S19.R1.fastq -2 S19/S19.R2.fastq -S S19.sam

Para lo que sigue, debemos activar el ambiente `Anvio`:

	conda activate anvio-8

Nótese que ahora aparece (`anvio-8`) antes del prompt `$`.

Convertir los archivos generados a formato sam:

	for f in *.sam; do samtools view -F 4 -bS $f > $f.raw; done

Ahora debemos inicializar los archivos generado con Anvio:

	for f in *.raw; do anvi-init-bam $f -o $f.bam; done

Limpieza de archivo ya no necesarios:

	rm *.raw *.sam

Como los nombres de los archivos van creciendo, es conveniente acortarlos con `rename`:

	rename 's/.sam.raw.bam/.bam/' *.bam

	rename 's/.sam.raw.bam.bai/.bai/' *.bai
	
Ahora tendremos nombre mas cortos para los archivos `.bam` y `.bai`
### Creación de base de datos
	anvi-gen-contigs-database -f contigs.fa -o contigs.db

Posteriormente es buena idea buscar hidden Markov models (hmm) en nuestros contigs, éstos son Single Copy Genes para bacterias y nos ayudarán a saber si un genoma bacteriano encontrado está "completo" o no:

	anvi-run-hmms -c contigs.db -I Bacteria_71 -T 2

También podemos extraer las secuencias de los genes existentes en la base de datos para una posterior clasificación taxonómica:

	anvi-get-sequences-for-gene-calls -c contigs.db -o gene-calls.fa
### Perfil para cada muestra
Teniendo ya todos los datos, debemos crear un "profile" para cada muestra, tener cuidado por que el subdirectorio va  a ser borrado y crear uno nuevo, pero de todos modos ya no necesitamos los archivos `fastq` y solo nos quitan nuestro preciado espacio.

	anvi-profile -i C08.bam -c contigs.db --output-dir C08 --sample-name C08 -W -T 2

Esto hay que hacerlo para cada una de las cinco muestras por separado:

	anvi-profile -i P08.bam -c contigs.db --output-dir P08 --sample-name P08 -W -T 2
	anvi-profile -i P18.bam -c contigs.db --output-dir P18 --sample-name P18 -W -T 2
	anvi-profile -i S04.bam -c contigs.db --output-dir S04 --sample-name S04 -W -T 2
	anvi-profile -i S19.bam -c contigs.db --output-dir S19 --sample-name S19 -W -T 2

Luego tenemos que unir todos los cinco profiles:

	anvi-merge */PROFILE.db -o SAMPLES-MERGED -c contigs.db --sample-name DIETAS

##### Limpieza

Ya no necesitamos los archivos `.bam` ni `.bai`, ni los creados durante la generación de índices, podemos eliminarlos:

	rm *.bai *.bam contigs.1.bt2 contigs.2.bt2 contigs.3.bt2 contigs.4.bt2 contigs.rev.*
	
***

#### Taxonomía
Como no corrimos clasificación taxonómica ni funcional, por falta de espacio en la imagen virtual, no podremos obtener información al respecto en la visualización.
En caso que hallamos corrido la clasificación taxonómica aparte y tengamos el resultado, podemos importarla a la base de datos.
En el set de datos *Dietas.tar.gz* tenemos los archivos pre generados en una carpeta que contiene cuatro archivos:

>`gene-calls.fa` Archivo fasta de los genes encontrados en los contigs.
>`gene-calls.tax` Taxonomía asociada a los genes anotados.
>`gene-calls.out` Anotación de los genes sin taxonomía.
>`gene-calls.krona` Archivo para generación de gráfica Krona.

El archivo `gene-calls.tax`  hay que incorporarlo a la base de datos contigs.db:

	anvi-import-taxonomy-for-genes -c contigs.db -i gene-calls.tax -p kaiju --just-do-it

Podemos generar una gráfica Krona con el archivo kaiju.krona

	ktImportText -o kaiju_gene_taxonomy.html kaiju.krona
#### Descripción
Tenemos también un archivo (`Dietas_descripcion.md`) con una breve descripción de los datos y métodos, podemos incluirlos en la base de datos para que sea desplegado en la visualización:

	anvi-update-db-description --description Dietas_descripcion.md SAMPLES-MERGED/PROFILE.db

#### Metadatos
Por último podemos importar los metadatos a Anvio:

	anvi-import-misc-data metadata.txt -p SAMPLES-MERGED/PROFILE.db --target-data-table layers

***

### Visualización de resultados

Por último podemos ver el resultado del análisis en una gráfica interactiva en el navegador:

	anvi-interactive -c contigs.db -p SAMPLES-MERGED/PROFILE.db --taxonomic-level t_species

***
