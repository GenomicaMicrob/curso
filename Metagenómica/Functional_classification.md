# Clasificación funcional con Superfocus

#### Archivos y programas
Los programas y set de datos ya se encuentran instalados en la imagen virtual **MGlinux18.2** pero la última versión de los datos se pueden bajar de aquí:

- [dietas.v2.tar.gz](https://drive.google.com/file/d/1FRdvMIERJJHSeo15n6vYlTtMvv3G16bZ/view?usp=sharing)
- [SUPERfocus](https://github.com/metageni/SUPER-FOCUS)
- [Krona](https://github.com/marbl/Krona)

***
Podemos **clasificar funcionalmente** metagenomas con [SUPERFOCUS](https://github.com/metageni/SUPER-FOCUS), pero debido al tamaño de la base de datos y la intensidad computacional, no es factible hacerlo en una máquina virtual. El proceso es igual a [FOCUS](Metagenómica/Taxonomic_metaclassification.md) y el resultado es una hoja tipo excel que podemos ver en la subcarpeta `cooked/superfocus/` en el set de datos `dietas.v2.tar.gz`. Este set de datos podemos bajarlo de [GDrive](https://drive.google.com/file/d/1FRdvMIERJJHSeo15n6vYlTtMvv3G16bZ/view?usp=sharing).

Superfocus clasifica cada secuencia siguiendo el esquema de [SEED](https://theseed.org/wiki/Home_of_the_SEED) que está basado en cuatro niveles jerárquicos. Superfocus genera varios archivos pero el principal es `*_All_levels.csv` que tiene tanto el conteo de secuencias como la abundancia relativa para cada muestra. Es conveniente separar estos dos grupos de valores y dejar así dos archivos; éstos son los que ya están en el set de datos (`dietas_numbers.txt` y `dietas_rel_abund.txt`).

***
### Arreglo del archivo `dietas_numbers.txt`

Vamos a graficar estos datos con R para ver las funciones clasificadas en las cinco muestras. Pero antes hay que modificar los nombres de las columnas del archivo dietas_numbers.txt para simplicifarlos y evitar problemas con los espacios. Actualmente están de la siguiente forma:

```
Subsystem Level 1	Subsystem Level 2	Subsystem Level 3	SEED Function	P08	S19	P18	C08	S04

Amino Acids and Derivatives	-	Amino acid racemase	2-methylaconitate_cis-trans_isomerase	0	0.00607096	0	0.000801796	0.004145507

Amino Acids and Derivatives	-	Amino acid racemase	2-methylaconitate_isomerase	0	0.001856409	0.001834223	0	0.000414551
```

Ejecutemos un comando para cambiar los nombres de las columnas:

 ```bash
sed -i 's/Subsystem Level 1/Level1/; s/Subsystem Level 2/Level2/; s/Subsystem Level 3/Level3/; s/SEED Function/SEED_Function/' dietas_numbers.txt
 ```

Ahora tendremos los nombre sin espacios y mas simples:

```
Level1  Level2 Level3 SEED_Function	P08	S19	P18	C08	S04
```
***
### Graficos de barras

Para graficar los datos podemos usar R usando el script [dietas_graficos.R](scripts/dietas_graficos.R) desde `RStudio` como se explica [Aqui](Graficacion_superfocus.md).

***
### Krona
También podemos generar una gráfica tipo [Krona](https://github.com/marbl/Krona) fácilmente; primero tenemos que extraer los datos para cada una de las muestras (con `cut`) y cambiar el formato (primero los datos y después los niveles, con `awk`):

```bash
cut -f1-5 dietas_numbers.txt | awk -v OFS='\t' 'BEGIN {FS="\t"}; {print $5, $1, $2, $3, $4}' | sed 1,4d > P08.krona
```
teniendo así un archivo para krona, generemos la gráfica interactiva krona:

```bash
ktImportText -o P08.html P08.krona
```
Este gráfico podemos visualizarlo con cualquier navegador. Para crear las gráficas de las demás muestras, solo tenemos que cortar las columnas correspondientes del archivo `dietas_numbers.txt`, p. ej. para la muestra S19: `cut -f1-4,6 dietas_numbers.txt ...` el resto del comando es igual al anterior, salvo el archivo de salida obviamente.
***

![Krona](plots/P08_krona.png)
***
