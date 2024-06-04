### awk_compiler
es un programa para compilar tablas con la misma información pero diferentes datos.

Uso: `awk_compiler *.tsv`

```bash
#!/usr/bin/awk -f
# Set field separator to tab (\t)
BEGIN { FS = "\t" }
# Remove carriage return characters if file is in DOS format (CRLF)
{ sub(/\r/, "") }
# Increment indx by 1 (starts at 0) everytime a new file is processed
!f[FILENAME]++ { ++indx }
# Add a key ($1) to keys array every time it is first encountered
!a[$1]++ { keys[i++] = $1 }
# Store the 2nd and 3rd field to b matrix
{ b[$1, indx] = $2 }
# This block runs after all files are processed
END {
    # Traverse the keys in order
    for (i = 0; i in keys; ++i) {
        key = keys[i]
        # Print key
        printf "%s", key
        # Print columns from every file in order
        for (j = 1; j <= indx; ++j) {
            v = b[key, j]
            printf "\t%s", length(v) ? v : "0"
        }
        # End the line with a newline
        print ""
    }
}
```
***
