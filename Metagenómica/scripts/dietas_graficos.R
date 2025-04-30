# dietas_graficos.R

# Cargar librerías
library(ggplot2)
library(reshape2)
library(scales)

# Leer datos
data <- read.delim("dietas_numbers.txt", sep = "\t", header = TRUE, check.names = FALSE)

# Detectar columnas numéricas, excluyendo columnas de metadatos
samples <- setdiff(names(Filter(is.numeric, data)), c("Level1"))

# --- Función: Agrupar y preparar datos sumados ---
preparar_datos_sumados <- function(data, samples) {
  df_sum <- aggregate(data[samples], by = list(Level1 = data$Level1), FUN = sum, na.rm = TRUE)
  df_long <- melt(df_sum, id.vars = "Level1", variable.name = "sample", value.name = "count")
  
  totals <- aggregate(count ~ Level1, data = df_long, sum)
  orden_niveles <- totals$Level1[order(totals$count)]
  df_long$Level1 <- factor(df_long$Level1, levels = orden_niveles)
  return(df_long)
}

# --- Función: Graficar distribución absoluta ---
graficar_distribucion <- function(df_long) {
  ggplot(df_long, aes(x = Level1, y = count, fill = sample)) +
    geom_col(position = "stack", width = 0.85) +
    coord_flip() +
    labs(title = "Distribución por categoría (orden ascendente)",
         x = "Categoría principal", y = "Conteo absoluto") +
    scale_fill_brewer(palette = "Set3", name = "Muestra") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 10))
}

# --- Función: Preparar datos porcentuales con top 15 + "Otras" ---
preparar_datos_porcentuales <- function(df_long, top_n = 15) {
  df_long <- df_long[!is.na(df_long$Level1), ]
  total_por_sample <- aggregate(count ~ sample, data = df_long, sum)
  df_perc <- merge(df_long, total_por_sample, by = "sample")
  df_perc$percentage <- df_perc$count.x / df_perc$count.y * 100

  # Calcular suma total por categoría
  suma_por_categoria <- aggregate(percentage ~ Level1, data = df_perc, sum)
  top_categorias <- head(suma_por_categoria$Level1[order(-suma_por_categoria$percentage)], top_n)

  # Reemplazar categorías no top por "Otras"
  df_perc$Level1 <- as.character(df_perc$Level1)
  df_perc$Level1[!(df_perc$Level1 %in% top_categorias)] <- "Otras"

  # Reagrupar porcentajes tras recodificar "Otras"
  df_perc <- aggregate(percentage ~ sample + Level1, data = df_perc, sum)

  # Ordenar niveles: menor arriba, mayor abajo, "Otras" al final
  orden_final <- aggregate(percentage ~ Level1, data = df_perc, sum)
  orden_final <- orden_final$Level1[order(orden_final$percentage)]
  orden_final <- c(setdiff(orden_final, "Otras"), "Otras")
  df_perc$Level1 <- factor(df_perc$Level1, levels = orden_final)

  return(df_perc)
}

# --- Función: Graficar distribución porcentual ---
graficar_porcentual <- function(df_perc) {
  n_cat <- length(unique(df_perc$Level1))
  paleta <- colorRampPalette(c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", 
                                "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", 
                                "#bcbd22", "#17becf"))(n_cat)

  ggplot(df_perc, aes(x = sample, y = percentage, fill = Level1)) +
    geom_col(position = "stack", width = 0.85) +
    labs(title = "Composición porcentual por muestra", x = "Muestra", y = "Porcentaje", fill = "Categorías Funcionales") +
    scale_y_continuous(labels = percent_format(scale = 1)) +
    scale_fill_manual(values = paleta) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.5, "cm")) +
    guides(fill = guide_legend(ncol = 1))
}

# --- Ejecutar flujo ---
df_long <- preparar_datos_sumados(data, samples)
p1 <- graficar_distribucion(df_long)
print(p1)

porcentual_df <- preparar_datos_porcentuales(df_long, top_n = 15)
p2 <- graficar_porcentual(porcentual_df)
print(p2)

# --- Guardar opcional ---
# ggsave("grafica1_orden_ascendente.png", p1, width = 10, height = 7)
# ggsave("grafica2_todas_categorias.png", p2, width = 10, height = 6, dpi = 300)
