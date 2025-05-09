# Curso Bioinformática Microbiana

*Bruno Gómez-Gil*

Instrucciones para los diferentes ejercicios en el curso de [Bioinformática microbiana](https://bioinformatica.ciad.mx/) del CIAD Mazatlán.

Estas instrucciones están dirigidos a realizarse con la imagen virtual [MGlinux18.2](https://bioinformatica.ciad.mx/programas/virtualizacion/mglinux) pues allí se tienen ya los sets de datos, programas y scripts para ello.
***
### Ejercicios

Los análisis que tenemos descritos aquí son y próximamente se agregarán nuevos:
#### Genómica
- [Ensamble de un genoma de SARS-CoV2](Genomica/SARS-CoV_analysis.md)
- [Ensamble de un genoma de *E. coli*](Genomica/Ecoli_assembly.md)
- [Pangenoma de *E.coli*](Genomica/Pangenoma_Ecoli.md)

#### 16S
- [Limpieza de secuencias](16S/Limpieza.md)
- [Eliminación de secuencias quiméricas](16S/Quimeras.md)
- [Clasificación taxonómica](16S/Clasificacion.md)

#### Metagenómica
- [Análisis metagenómico con Anvio](Metagenómica/Shotgun_analysis.md)
- [Análisis de amplicones 16-V4](Metagenómica/16S-V4_analysis.md)
- [Clasificación taxonómica de metagenomas](Metagenómica/Taxonomic_metaclassification.md)
- [Clasificación funcional con Superfocus](Metagenómica/Functional_classification.md)

#### Transcriptómica
- [Análisis transcriptómico de *Vibrio parahaemolyticus*](Transcriptomica/Transcriptomica.md)
- [Estimación de genes diferencialmente expresados con DESeq2](Transcriptomica/DESeq2.md)
***

### Datasets
Los **set de datos** comprimidos también se pueden descargar de los siguientes repositorios:

| Dataset | Tamaño (MB) | Repositorios | Descripción |
| --- | ---: | :---: | --- |
| dietas.v2.tar.gz | 148.6 | [Figshare](https://figshare.com/s/4e700c8c9ce853e74827) [GDrive](https://drive.google.com/file/d/1FRdvMIERJJHSeo15n6vYlTtMvv3G16bZ/view?usp=sharing) | Cinco metagenomas de heces de camarón recortados a 100,000 secuencias pareadas y diversos análisis ya pre hechos. |
| Infant-Gut-Tutorial_v8.tar.gz | 118.9 | [GDrive](https://drive.google.com/file/d/1w3Ie2eOZaIbTffg8moX-iJwMyztPgA6i/view?usp=sharing) | Análisis con Anvio de muestras metagenómicas de heces humanas |
| MiSeqSOP.tar.gz | 33.3 | [Figshare](https://figshare.com/s/f21ca7e71285396a1020) [GDrive](https://drive.google.com/file/d/1qmUd0AYSAjj2ND15Bfi2koTAU5GWRk1n/view?usp=drive_link) | Muestras de heces de rata a diferentes tiempos, región V4 del gen ribosomal 16S, secuencias pareadas 2x150. |
| SARSCoV2.tar.gz | 96.9 | [Figshare](https://figshare.com/s/766f0052088f2dab119c) [GDrive](https://drive.google.com/file/d/1-FMZkRvzgubmlyi2qP6J8AQJm5HeYDRC/view?usp=sharing) | Tres muestras de secuencias de SARS-CoV2 (alfa, delta y omicrón) para ensamble. |
| ecoli-K12_Illumina.tar.gz | 81.5 | [GDrive](https://drive.google.com/file/d/1NOcflmwa6ioLDOjFVhhl5TdhJbIgpBO1/view?usp=sharing) | Secuencias pareadas de *E. coli* K12. |
| Francisella_pan.tar.gz | 3.1 | [Figshare](https://figshare.com/ndownloader/files/42746980) [GDrive](https://drive.google.com/file/d/1pBQLuCk5O9m2wfG8QdJ6EoLKd0Yc-03X/view?usp=drive_link) | Seis genomas de *Francisella tularensis* para realizar un análisis de pangenoma de la especie. |
| ecoli.tar.gz | 149.5 | [GDrive](https://drive.google.com/file/d/1OSoJIfb7kkdGx4rJrfHHKYqy-d_Ucejw/view?usp=drive_link) | Tres genomas de *E. coli* análisis pangenómicos y secuencias Illumina y Ion Torrent para ensamble. |
| transcritos.tar.gz | 176.3 | [GDrive](https://drive.google.com/file/d/1hHg_qHgwy7fPg5aPLimKPfNfZpRiBG85/view?usp=sharing) | Seis secuencias de transcritos de *Vibrio parahaemolyticus* crecido en dos medios de cultivo.

GDrive = Google Drive

Los set de datos también se pueden descargar desde el servidor bioinformático **Biobacter** del CIAD; se encuentran en la siguiente ruta:
`/raid1/datasets/`
***
