# Script para analizar genes diferencialmente expresados con R en RStudio
# Se necesita tener ya importados el conteo de transcritos (archivo counts.csv como matriz count_matrix) y los metadatos (coldata).
# Proceso tomado de: https://biostatsquid.com/volcano-plots-r-tutorial/

library("DESeq2")
reference <- readline(prompt="Tratamiento de referencia: ") # Mencionar cual será el tratamiento de referencia

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~ Treatment) # change Treatment if necessary
dds <- dds[rowSums(counts(dds)) >= 10,] # FILTER. remove the genes which have < 10 reads
dds$Treatment <- relevel(dds$Treatment, ref = reference)
dds <- DESeq(dds)# Perform differential gene expression analysis
resultsNames(dds) # see all comparisons (here there is only one)
res <- results(dds) # get gene expression table at this step independent filtering is applied by default to remove low count genes
res
res[order(res$padj),] # Order gene expression table by adjusted p value (Benjamini-Hochberg FDR method).
write.csv(as.data.frame(res[order(res$padj),] ), file = "DGE.csv") # Export differential gene expression analysis table to CSV file
summary(results(dds, alpha=0.05)) # Get summary of differential gene expression with adjusted p value cut-off at 0.05.
normalized_counts <- counts(dds, normalized=TRUE) # Get normalized counts.
head(normalized_counts)
# FIN
