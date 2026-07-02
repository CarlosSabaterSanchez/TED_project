library(readxl)
library(dplyr)
library(writexl)
library(stringr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(vegan)
library(pheatmap)

##BSG HEATMAP CON ANOTACION COLUMNAS 

# ---- DATOS AZUCARES ----
datosbsg <- read_excel("BSG_2025.xlsx")
datosbsg <- as.data.frame(datosbsg)
rownames(datosbsg) <- datosbsg$Isolate
datosbsg$Isolate <- NULL

# ---- NORMALIZAR AZUCARES ----
data_scaled_bsg <- scale(as.matrix(datosbsg))
data_scaled_bsg

#------	DATOS CONTEOS Y PH ------
counts_ph_bsg <- read_excel(
  "COUNTS-PH.xlsx",
  sheet = "BSG"
)

counts_ph_bsg <- as.data.frame(counts_ph_bsg)

# Renombrar columnas (para evitar problemas con espacios)
colnames(counts_ph_bsg) <- c(
  "Isolate",
  "T0_cfu",
  "T24_cfu",
  "T48_cfu",
  "T0_pH",
  "T24_pH",
  "T48_pH"
)

# ==================================
# 3. CALCULAR VARIABLES
# ==================================
str(counts_ph_bsg)

counts_ph_bsg <- counts_ph_bsg %>%
  mutate(
    T0_cfu  = as.numeric(T0_cfu),
    T24_cfu = as.numeric(T24_cfu),
    T48_cfu = as.numeric(T48_cfu),
    
    T0_pH  = as.numeric(T0_pH),
    T24_pH = as.numeric(T24_pH),
    T48_pH = as.numeric(T48_pH)
  )

counts_ph_bsg

counts_ph_bsg <- counts_ph_bsg %>%
  mutate(
    Caida_pH = T0_pH - T48_pH,
    Incremento_cfu = T48_cfu - T0_cfu
  )

summary(counts_ph_bsg)

# ==================================
# 4. ORDENAR SEGÚN HEATMAP
# ==================================

# # Asegurar mismo orden que azúcares
# counts_ph_bsg <- counts_ph_bsg[
#   match(rownames(data_scaled_bsg), counts_ph_bsg$Isolate),
# ]
# counts_ph_bsg

# ==================================
# 4. ANOTACIONES
# ==================================

annotation_col_bsg <- data.frame(
  Caida_pH = counts_ph_bsg$Caida_pH,
  Incremento_cfu = counts_ph_bsg$Incremento_cfu
)

annotation_col_bsg

rownames(annotation_col_bsg) <- counts_ph_bsg$Isolate

# Colores
ann_colors <- list(
  Caida_pH = colorRampPalette(c("white", "blue"))(100),
  Incremento_cfu = colorRampPalette(c("white", "darkgreen"))(100)
)

# ==================================
# 6. HEATMAP FINAL (con orden correcto)
# ==================================

# Dibujar una vez para obtener el orden del clustering
p_tmp <- pheatmap(
  t(data_scaled_bsg),
  
  scale = "none",
  
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  
  color = colorRampPalette(
    c("navy", "white", "firebrick3")
  )(100),
  
  fontsize_row = 10,
  fontsize_col = 10,
  
  annotation_col = annotation_col_bsg,
  annotation_colors = ann_colors,
  
  silent = TRUE   # importante para que no lo imprima aún
)

# Extraer orden real de columnas tras clustering
orden_cluster <- p_tmp$tree_col$order

# Reordenar annotation_col según clustering
annotation_col <- annotation_col_bsg[
  colnames(t(data_scaled_bsg))[orden_cluster],
]

# Ahora sí generar heatmap definitivo
svg("heatmap_azucares-bsg_and_counts_ph_PRUEBA.svg",
    width = 8,
    height = 8)

pheatmap(
  t(data_scaled_bsg),
  
  scale = "none",
  
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  
  color = colorRampPalette(
    c("navy", "white", "firebrick3")
  )(100),
  
  annotation_col = annotation_col_bsg,
  annotation_colors = ann_colors,
  
  main = "Normalized sugar with pH and counts changes",
  
  fontsize_row = 10,
  fontsize_col = 10,
  
  border_color = NA
)

dev.off()
#-------------------------------------------
#BSG NMDS
#-------------------------------------------

# ---- LEER ARCHIVOS GENES ----
pectin    <- read_excel("Pectin.xlsx")
cellulose <- read_excel("Cellulose.xlsx")
xylan     <- read_excel("Xylan.xlsx")
glucan    <- read_excel("Glucan.xlsx")

# ---- COMBINAR ARCHIVOS ----
datos_completos <- pectin %>%
  inner_join(cellulose, by = "ID") %>%
  inner_join(xylan,   by = "ID") %>%
  inner_join(glucan,  by = "ID")
  
  # ---- DISTANCIAS Y NMDS ----
dist_matrixbsg <- vegdist(data_scaled_bsg, method = "euclidean")
set.seed(123)
nmds_euclidean_bsg <- metaMDS(dist_matrixbsg, k = 2, trymax = 100)
nmds_euclidean_bsg # stress 0.15686 - débil, probamos a ver si k = 3 mejora, aunque representemos sólo
## k = 2 se podria citar que k = 3 mejora patron con su estadística

nmds_euclidean_bsg3 <- metaMDS(dist_matrixbsg, k = 3, trymax = 100)
nmds_euclidean_bsg3

# ---- COORDENADAS DE SITES ----
sitesbsg <- as.data.frame(scores(nmds_euclidean_bsg, display = "sites"))
sitesbsg$Isolate <- rownames(sitesbsg)

# ---- VARIABLES CUALITATIVAS ----
# Usamos toda la tabla de genes binarios como variables cualitativas para envfit
# Excluimos la columna ID para usar solo los datos binarios
vars_cual <- datos_completos %>% select(-ID)

# Filtrar variables constantes (desviación estándar = 0)
vars_cual2 <- vars_cual[, apply(vars_cual, 2, function(x) length(unique(x)) > 1)]

# Ajustar envfit
fit_braybsg <- envfit(nmds_euclidean_bsg, vars_cual2, permutations = 999, na.rm = TRUE)

fit_scoresbsg <- as.data.frame(scores(fit_braybsg, "vectors")) * ordiArrowMul(fit_braybsg)
fit_scoresbsg$Variable <- rownames(fit_scoresbsg)
fit_scoresbsg$r2   <- fit_braybsg$vectors$r
fit_scoresbsg$pval <- fit_braybsg$vectors$pvals

# ---- VECTORES SIGNIFICATIVOS (p < 0.1) ----
fit_scores_sigbsg <- subset(fit_scoresbsg, pval < 0.1)
fit_scores_sigbsg

arrow_factor <- 1.5   # prueba entre 1.2 y 2

fit_scores_sigbsg <- fit_scores_sigbsg %>%
  mutate(
    NMDS1_s = NMDS1 * arrow_factor,
    NMDS2_s = NMDS2 * arrow_factor
  )

fit_scores_sigbsg

bsg_scaled <- ggplot(sitesbsg, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 2.5, color = "steelblue", alpha = 0.7) +
  geom_text_repel(aes(label = Isolate), size = 4, max.overlaps = Inf,
                  box.padding = 0.6, point.padding = 0.4,
                  segment.size = 0.3, force = 2) +
  # vectores solo si hay significativos
  {if(nrow(fit_scores_sigbsg) > 0)
    geom_segment(data = fit_scores_sigbsg,
                 aes(x = 0, y = 0, xend = NMDS1_s, yend = NMDS2_s),
                 arrow = arrow(length = unit(0.15, "cm")),
                 color = "red", inherit.aes = FALSE)
  } +
  {if(nrow(fit_scores_sigbsg) > 0)
    geom_text_repel(data = fit_scores_sigbsg,
                    aes(x = NMDS1_s, y = NMDS2_s, label = Variable),
                    size = 4, color = "red",
                    max.overlaps = Inf, box.padding = 0.6, point.padding = 0.4,
                    force = 2, inherit.aes = FALSE)
  } +
  coord_equal(
    xlim = c(min(sitesbsg$NMDS1) - 0.2, max(sitesbsg$NMDS1) + 0.2),
    ylim = c(min(sitesbsg$NMDS2) - 0.2, max(sitesbsg$NMDS2) + 0.2),
    expand = FALSE
  ) +
  theme_minimal(base_size = 12) +
  labs(title = "NMDS (CAZY PROFILES bsg)")



bsg_scaled

#------------------------------------------------
##BAB HEATMAP CON ANOTACION COLUMNAS 
#-----------------------------------------------
datos <- read_excel("BAB_2025.xlsx")
datos <- as.data.frame(datos)
rownames(datos) <- datos$Isolate
datos$Isolate <- NULL

# ---- NORMALIZAR AZUCARES ----
data_scaled <- scale(as.matrix(datos))
data_scaled

# ---- DATOS PH Y CRECIMIENTO

counts_ph <- read_excel(
  "COUNTS-PH.xlsx",
  sheet = "BAB"
)

counts_ph <- as.data.frame(counts_ph)

# Renombrar columnas (para evitar problemas con espacios)
colnames(counts_ph) <- c(
  "Isolate",
  "T0_cfu",
  "T24_cfu",
  "T48_cfu",
  "T0_pH",
  "T24_pH",
  "T48_pH"
)

# ------ CALCULAR VARIABLES

str(counts_ph)

counts_ph <- counts_ph %>%
  mutate(
    T0_cfu  = as.numeric(T0_cfu),
    T24_cfu = as.numeric(T24_cfu),
    T48_cfu = as.numeric(T48_cfu),
    
    T0_pH  = as.numeric(T0_pH),
    T24_pH = as.numeric(T24_pH),
    T48_pH = as.numeric(T48_pH)
  )

counts_ph <- counts_ph %>%
  mutate(
    pH_reduction = T0_pH - T48_pH,
    log_cfuml_increase = T48_cfu - T0_cfu
  )

summary(counts_ph)

# ----- ANOTACIONES

annotation_col <- data.frame(
  pH_reduction = counts_ph$pH_reduction,
  log_cfuml_increase = counts_ph$log_cfuml_increase
)

annotation_col

rownames(annotation_col) <- counts_ph$Isolate

# Colores
ann_colors <-
  list(
  pH_reduction = colorRampPalette(c("white", "blue"))(100),
  log_cfuml_increase = colorRampPalette(c("white", "darkgreen"))(100)
)


# ----- HEATMAP FINAL (con orden correcto)

# 1 Dibujar una vez para obtener el orden del clustering
p_tmp <- pheatmap(
  t(data_scaled),
  
  scale = "none",
  
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  
  color = colorRampPalette(
    c("navy", "white", "firebrick3")
  )(100),
  
  fontsize_row = 10,
  fontsize_col = 10,
  
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  
  silent = TRUE   # importante para que no lo imprima aún
)

# 2️ Extraer orden real de columnas tras clustering
orden_cluster <- p_tmp$tree_col$order
orden_cluster

# 3 Reordenar annotation_col según clustering
annotation_col <- annotation_col[
  colnames(t(data_scaled))[orden_cluster],
]

annotation_col

# Ahora sí generar heatmap definitivo
svg("heatmap_azucares-bab_and_counts_ph_130226.svg",
    width = 8,
    height = 8)

pheatmap(
  t(data_scaled),
  
  scale = "none",
  
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  
  color = colorRampPalette(
    c("navy", "white", "firebrick3")
  )(100),
  
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  
  main = "Normalized sugar with pH and counts changes",
  
  fontsize_row = 10,
  fontsize_col = 10,
  
  border_color = NA
)

dev.off()

##---------------------------------
##BAB NMDS
##---------------------------------
# ---- DISTANCIAS Y NMDS ----
dist_matrix <- vegdist(data_scaled, method = "euclidean")
set.seed(123)
nmds_euclidean <- metaMDS(dist_matrix, k = 2, trymax = 100)
nmds_euclidean

# ---- COORDENADAS DE SITES ----
sites <- as.data.frame(scores(nmds_euclidean, display = "sites"))
sites$Isolate <- rownames(sites)

# ---- VARIABLES CUALITATIVAS ----
# Usamos toda la tabla de genes binarios como variables cualitativas para envfit
# Excluimos la columna ID para usar solo los datos binarios
vars_cual <- datos_completos %>% select(-ID)

# Filtrar variables constantes (desviación estándar = 0)
vars_cual2 <- vars_cual[, apply(vars_cual, 2, function(x) length(unique(x)) > 1)]

# Ajustar envfit
fit_bray <- envfit(nmds_euclidean, vars_cual2, permutations = 999, na.rm = TRUE)

fit_scores <- as.data.frame(scores(fit_bray, "vectors")) * ordiArrowMul(fit_bray)
fit_scores$Variable <- rownames(fit_scores)
fit_scores$r2   <- fit_bray$vectors$r
fit_scores$pval <- fit_bray$vectors$pvals

# ---- VECTORES SIGNIFICATIVOS (p < 0.1) ----
fit_scores_sig <- subset(fit_scores, pval < 0.1)

# ---- PLOT NMDS ----
bab <- ggplot(sites, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 2.5, color = "steelblue", alpha = 0.7) +
  geom_text_repel(aes(label = Isolate), size = 3, max.overlaps = Inf,
                  box.padding = 0.6, point.padding = 0.4,
                  segment.size = 0.3, force = 2) +
  # vectores solo si hay significativos
  {if(nrow(fit_scores_sig) > 0)
    geom_segment(data = fit_scores_sig,
                 aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                 arrow = arrow(length = unit(0.15, "cm")),
                 color = "red", inherit.aes = FALSE)
  } +
  {if(nrow(fit_scores_sig) > 0)
    geom_text_repel(data = fit_scores_sig,
                    aes(x = NMDS1, y = NMDS2, label = Variable),
                    size = 3, color = "red",
                    max.overlaps = Inf, box.padding = 0.6, point.padding = 0.4,
                    force = 2, inherit.aes = FALSE)
  } +
  coord_equal(
    xlim = c(min(sites$NMDS1) - 0.2, max(sites$NMDS1) + 0.2),
    ylim = c(min(sites$NMDS2) - 0.2, max(sites$NMDS2) + 0.2),
    expand = FALSE
  ) +
  theme_minimal(base_size = 12) +
  labs(title = "NMDS (CAZY PROFILES - Significant vectors only)")

bab

