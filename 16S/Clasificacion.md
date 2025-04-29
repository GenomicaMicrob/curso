### Clasificación de secuencias
Por último paso, tenemos que clasificar las secuencias para asignarles una taxonomía, esto lo haremos con el script `mg_classifier`:

```bash
mg_classifier *.fasta
```

Nos saldrá un menu para seleccionar la base de datos a usar

```latex
Select a database:
16S_rRNA
   1 EzBioCloud-LGM 2023
   2 SILVA v132
18S_rRNA
   3 PR2

   x exit
```

En la imagen virtual **MGlinux18.2**, por cuestiones de espacio, no tenemos todas las bases de datos, elijamos la **número 1**. Después de pocos minutos, tendremos una nueva subcarpeta con los resultados. El documento principal es una tabla de OTUs en la que está cuantas secuencias de cada OTU (líneas) se asignaron a cada muestra (columnas) analizada.
***
