### Instalación de fadu

Fadu es un programa para analizar transcriptómica de procariontes únicamente. Si se tiene la imagen virtual MGlinux 18.2 o anteriores, es necesario instalarlo en ella.
***
### Instalación del lenguaje Julia
Ya que FADU corre en Julia, seránecesario descargar el lenguaje e instalarlo:

```bash
wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.3-linux-x86_64.tar.gz
```
Una vez descargado hay que descomprimirlo

```bash
tar xzf julia-1.10.3-linux-x86_64.tar.gz
```
así tendremos una nueva carpeta que la podemos mover al directorio `/opt/`

```bash
sudo mv julia-1.10.3/ /opt/
```
Recordar que al usar `sudo` nos preguntará el password: `user01`

Borremos el comprimido para hacer espacio

```bash
rm julia-1.10.3-linux-x86_64.tar.gz
```
Ahora es bueno hacer un link virtual al executable de julia para no tener que poner la ruta cada vez que queramos ejecutarlo:

```bash
sudo ln -s /opt/julia-1.10.3/bin/julia /usr/local/bin/julia
```
Para comprobar que todo funciona, cerremos la terminal, abramos una nueva y escribamos `julia` y debe observarse julia.

Es necesario instalar algunos paquetes para `julia`; una vez ya en el ambiente julia debemos entrar al manejador de paquetes presionando la tecla `]` y veremos que el prompt cambia a `(@v1.10) pkg>` ahora si instalaremos los paquetes necesarios en conjunto:

```bash
add ArgParse Logging BGZFStreams.jl GenomicFeatures.jl GFF3 Indexes.jl StructArrays.jl XAM.jl BED.jl
```
Podemos salir del `package manager` y luego de julia presionando `Ctrl + c` y luego `exit()`

***
### Descarga e instalación de FADU

Descarguemos la versión 1.9.0 de FADU
```bash
wget https://github.com/IGS/FADU/archive/refs/tags/v1.9.0.tar.gz
```
Descomprimamos el archivo directo a `/opt/`:

```bash
sudo tar xzf v1.9.0.tar.gz -C /opt/
```
Ya descargado y descomprimido, necesitamos cambiar al dueño del paquete para poderlo ejecutar ya que esta como `root` por default:

```bash
sudo chown -R user /opt/FADU/1.9.0/*
```
Hagámoslo ejecutable
```bash
sudo chmod +x /opt/FADU/1.9.0/fadu.jl
```
***
