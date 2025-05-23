#!/usr/bin/env python3
#Script para sustituir el ID en un archivo gff por el nombre del gen y ponerle un num consecutivo
# USO: gff_ID_renaming.py archivo.gff
import re
import argparse
import os

# Configurar el analizador de argumentos
parser = argparse.ArgumentParser(description='Modificar tabla GFF.')
parser.add_argument('input_file', type=str, help='Archivo de entrada')
args = parser.parse_args()

# Generar el nombre del archivo de salida
input_filename = args.input_file
file_root, file_ext = os.path.splitext(input_filename)
output_filename = f"{file_root}_mod{file_ext}"

# Leer el contenido del archivo de entrada
with open(input_filename, 'r') as file:
    lines = file.readlines()

# Inicializar el contador para el ID
contador = 1

# Inicializar una lista para almacenar las líneas modificadas
modified_lines = []

# Recorrer cada línea en el archivo
for line in lines:
    # Usar una expresión regular para encontrar y reemplazar el ID y eliminar Name
    modified_line = re.sub(r'ID=.*?;Name=', f'ID={contador}-', line)
    # Incrementar el contador
    contador += 1
    # Añadir la línea modificada a la lista
    modified_lines.append(modified_line)

# Escribir las líneas modificadas en el archivo de salida
with open(output_filename, 'w') as file:
    file.writelines(modified_lines)

print(f'Tabla modificada guardada en {output_filename}')
