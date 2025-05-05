# -------------------------------------------------------------
# Análisis transcriptómico de célula única (scRNA-seq) en Trucha:
# Integración de muestras y análisis de expresión génica
# -------------------------------------------------------------

# Cargar bibliotecas
library(Seurat)
library(ggplot2)
library(dplyr)
library(scDblFinder)
library(clustree)

set.seed(48729351)

# Funciones
qc_vlnplot <- function(objeto, features, group.by = "orig.ident", pt.size = 0.1, ncol = 3) {
  plots <- VlnPlot(objeto, features = features, group.by = group.by, pt.size = pt.size, ncol = ncol)
  titulos <- c("Genes por célula", "UMIs por célula", "Contenido mitocondrial")
  y_labels <- c("Número de genes", "UMIs totales", "% genes mitocondriales")

  for (i in seq_along(features)) {
    plots[[i]] <- plots[[i]] +
      labs(title = titulos[i], y = y_labels[i]) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  return(plots)
}

# 1. Cargar datos --------------------------------------------------------------
sp1 <- Read10X(data.dir = "Datos/Trucha/SP1_filtered_feature_bc_matrix")
sp2 <- Read10X(data.dir = "Datos/Trucha/SP2_filtered_feature_bc_matrix")
sp3 <- Read10X(data.dir = "Datos/Trucha/SP3_filtered_feature_bc_matrix")

# 2. Crear objetos Seurat y QC -------------------------------------------------
muestra1 <- CreateSeuratObject(sp1, project = "SP1", min.cells = 3, min.features = 200)
muestra2 <- CreateSeuratObject(sp2, project = "SP2", min.cells = 3, min.features = 200)
muestra3 <- CreateSeuratObject(sp3, project = "SP3", min.cells = 3, min.features = 200)

# Calcular porcentaje mitocondrial (MT-)
muestra1[["percent.mt"]] <- PercentageFeatureSet(muestra1, pattern = "^MT-")
muestra2[["percent.mt"]] <- PercentageFeatureSet(muestra2, pattern = "^MT-")
muestra3[["percent.mt"]] <- PercentageFeatureSet(muestra3, pattern = "^MT-")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Violin plot de nFeature_RNA y nCount_RNA para las tres muestras
combinado <- merge(muestra1, y = c(muestra2, muestra3), add.cell.ids = c("SP1", "SP2", "SP3"))
pre_limpieza <- qc_vlnplot(combinado, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

# Filtrar células. Conservar:
# - percent.mt <= 20%
# - nFeature_RNA < 4500
muestra1 <- subset(muestra1, subset = percent.mt <= 20 & nFeature_RNA < 4500)
muestra2 <- subset(muestra2, subset = percent.mt <= 20 & nFeature_RNA < 4500)
muestra3 <- subset(muestra3, subset = percent.mt <= 20 & nFeature_RNA < 4500)

# Convertir a SingleCellExperiment y correr scDblFinder
muestra1.sce <- as.SingleCellExperiment(muestra1)
muestra1.sce <- scDblFinder(muestra1.sce)
muestra1 <- AddMetaData(muestra1, colData(muestra1.sce)[, "scDblFinder.class"], col.name = "scDblFinder.class")

muestra2.sce <- as.SingleCellExperiment(muestra2)
muestra2.sce <- scDblFinder(muestra2.sce)
muestra2 <- AddMetaData(muestra2, colData(muestra2.sce)[, "scDblFinder.class"], col.name = "scDblFinder.class")

muestra3.sce <- as.SingleCellExperiment(muestra3)
muestra3.sce <- scDblFinder(muestra3.sce)
muestra3 <- AddMetaData(muestra3, colData(muestra3.sce)[, "scDblFinder.class"], col.name = "scDblFinder.class")

# Filtrar células dobletes (conservar singletes)
muestra1 <- subset(muestra1, scDblFinder.class == "singlet")
muestra2 <- subset(muestra2, scDblFinder.class == "singlet")
muestra3 <- subset(muestra3, scDblFinder.class == "singlet")

combinado <- merge(muestra1, y = c(muestra2, muestra3), add.cell.ids = c("SP1", "SP2", "SP3"))
post_limpieza <- qc_vlnplot(combinado, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

cat("muestra1:", nrow(muestra1), "genes,", ncol(muestra1), "células\n")
cat("muestra2:", nrow(muestra2), "genes,", ncol(muestra2), "células\n")
cat("muestra3:", nrow(muestra3), "genes,", ncol(muestra3), "células\n")

ggsave("Resultados/Trucha_MSP/vln_sin_QC_trucha.png", plot = pre_limpieza, width = 20, height = 8, dpi = 500)
ggsave("Resultados/Trucha_MSP/vln_con_QC_trucha.png", plot = post_limpieza, width = 20, height = 8, dpi = 500)
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Juntar muestras en una lista
lista_muestras <- list(muestra1, muestra2, muestra3)

# 3. Normalizar e integrar muestras --------------------------------------------
# A. Normalizar cada muestra y escalar datos
lista_muestras <- lapply(lista_muestras, SCTransform)

# B. Seleccionar genes para integrar
genes <- SelectIntegrationFeatures(lista_muestras, nfeatures = 3500)

# C. Preparar datos para integración
lista_muestras <- PrepSCTIntegration(
  object.list = lista_muestras,
  anchor.features = genes
)

# D. Encontrar anclas de integración
anclas <- FindIntegrationAnchors(
  object.list = lista_muestras,
  normalization.method = "SCT",
  anchor.features = genes
)

# E. Integrar
datos_integrados <- IntegrateData(anclas, normalization.method = "SCT")

# 4. Análisis básico -----------------------------------------------------------

# A. Reducción de dimensionalidad
# PCA
DefaultAssay(datos_integrados) <- "integrated"

datos_integrados <- RunPCA(datos_integrados)

# Gráfico de codo para elegir el número de dimensiones
elbow_plot <- ElbowPlot(datos_integrados, ndims = 30)

# UMAP y TSNE
datos_integrados <- RunUMAP(datos_integrados, dims = 1:30)
datos_integrados <- RunTSNE(datos_integrados, dims = 1:30)

# B. Clustering
datos_integrados <- FindNeighbors(datos_integrados, dims = 1:30)

# Comparación de resoluciones de clúster
datos_integrados <- FindClusters(datos_integrados, resolution = 1.2)
datos_integrados$resolution_1.2 <- datos_integrados$seurat_clusters

datos_integrados <- FindClusters(datos_integrados, resolution = 0.4)
datos_integrados$resolution_0.4 <- datos_integrados$seurat_clusters

datos_integrados <- FindClusters(datos_integrados, resolution = 0.6)
datos_integrados$resolution_0.6 <- datos_integrados$seurat_clusters

# 5. Identificación de marcadores diferenciales --------------------------------
# Análisis de Expresión Diferencial (Genes sobreexpresados en cada clúster vs el resto)
DefaultAssay(datos_integrados) <- "RNA"
datos_integrados <- JoinLayers(datos_integrados)
datos_integrados <- NormalizeData(datos_integrados)

marcadores <- FindAllMarkers(
  object = datos_integrados,
  assay = "RNA",
  slot = "data",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Seleccionar los 20 marcadores más significativos por p-valor
top_marcadores <- marcadores %>%
  dplyr::group_by(cluster) %>% # Agrupar por clúster
  dplyr::arrange(p_val_adj) %>% # Ordenar por p_val_adj
  dplyr::slice(1:20) %>% # Seleccionar las primeras 20 filas
  dplyr::ungroup() # Desagrupar

# Seleccionar 3 marcadores diferenciales más importantes de cada clúster para heatmap
top3_marcadores <- top_marcadores %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice(1:3) %>%
  dplyr::ungroup()

# 6. Generación de resultados --------------------------------------------------
DefaultAssay(datos_integrados) <- "RNA"
# - Heatmap de los 3 marcadores diferenciales más significativos por clúster
heatmap_marcadores <- DoHeatmap(
  datos_integrados,
  features = top3_marcadores$gene,
  size = 6,
  assay = "integrated",
  angle = 55
) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 14, face = "bold"),
    legend.key.size = unit(0.7, "cm")
  )

# Asignación de clústeres a tipo celular (tras AED)
new.cluster.ids <- c(
  "0. B tempranas (IgM⁺)",
  "1. Megacariocitos I",
  "2. Megacariocitos II",
  "3. T CD8⁺ efector migratorio (IL21R⁺)",
  "4. B plasmáticas (Igλ⁺ Igκ⁺)",
  "5. B naïve (IgM⁺)",
  "6. T CD4⁺ memoria central (CCR7⁺ CD28⁺)",
  "7. T CD8⁺ efector/memoria central (FOXO1⁺ SATB1⁺ CCR7⁺)",
  "8. T CD8⁺ efector/proliferativo (MKI67⁺)",
  "9. T CD8⁺ efector terminal (PRF1⁺ GZMB⁺)",
  "10. Eritrocitos",
  "11. Neutrófilos",
  "12. NCC",
  "13. Macrófagos I",
  "14. Macrófagos II",
  "15. cDCs",
  "16. B plasmáticas maduras (Igκ⁺ CXCR4⁺)",
  "17. NK-Like"
)

names(new.cluster.ids) <- levels(datos_integrados)
datos_integrados <- RenameIdents(datos_integrados, new.cluster.ids)

# Gráficos UMAP
UMAP_Muestras <- DimPlot(
  datos_integrados,
  group.by = "orig.ident",
  reduction = "umap",
  label.size = 3.5,
  pt.size = 0.4
) + labs(title = "Muestras (UMAP)") +
  NoLegend()

UMAP_Clusters <- DimPlot(
  datos_integrados,
  reduction = "umap",
  label = TRUE,
  label.size = 3.5,
  pt.size = 0.4
) +
  labs(title = "Clústeres (UMAP)") +
  theme(
    axis.title = element_text(size = 20), axis.text = element_text(size = 18),
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 4))) +
  NoLegend()

# Gráficos t-SNE
TSNE_Muestras <- DimPlot(
  datos_integrados,
  group.by = "orig.ident",
  reduction = "tsne",
  label.size = 3.5,
  pt.size = 0.4
) + labs(title = "Muestras (t-SNE)")

TSNE_Clusters <- DimPlot(
  datos_integrados,
  reduction = "tsne",
  label = TRUE,
  label.size = 3.5,
  pt.size = 0.4
) +
  labs(title = "Clústeres (t-SNE)") +
  theme(
    axis.title = element_text(size = 20), axis.text = element_text(size = 18),
    legend.text = element_text(size = 12)
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 4)))

UMAP_Clusters + TSNE_Clusters

# Otros gráficos.
# - DotPlot del gen más expresado por cada clúster
genes_clasicos <- unique(c(
  "pax-5", "blnk", "LOC110538697",
  "hyal3", "thbs1b",
  "LOC110491676", "f13a1b",
  "LOC110530298", "dock10", "il7r",
  "LOC110487484", "LOC118936491", "LOC110496825",
  "LOC110533327", "LOC110486747",
  "LOC110534451", "LOC110534358", "LOC118941410",
  "LOC110537729", "LOC118964336", "tcf7",
  "sh2d1ab", "s100w", "LOC110531827",
  "il2rb", "LOC110509811", "LOC110508453",
  "LOC110489254", "LOC110538445",
  "LOC100136017", "mmp9", "lect2",
  "si:dkey-9i23.4", "LOC110534434", "LOC110531658",
  "LOC100136950", "mpeg1.1", "plxdc2",
  "LOC100136950", "csf3r", "lyz2",
  "LOC110535338", "flt3", "si:dkey-88e18.2",
  "zgc:152968", "LOC110500099", "LOC110499928",
  "LOC110536507", "LOC110498134", "LOC110491495"
))

genes_nuevos_linfocitos <- unique(c(
  # Linfocitos B
  "LOC110485215", # CL0
  "LOC118938876", # CL0
  "LOC110533327", # CL0, CL5
  "LOC110495722", # CL4
  "LOC110500016", # CL4
  "LOC110533868", # CL4
  "LOC110527864", # CL5
  "si:dkey-24p1.1", # CL5
  "LOC118942906", # CL16
  "LOC110526114", # CL16
  # Linfocitos T
  "si:ch211-67e16.3", # CL3
  "LOC110528322", # CL3
  "LOC110534358", # CL6
  "LOC110496534", # CL6
  "LOC118937335", # CL7
  "LOC110529458", # CL7
  "s100w", # CL8
  "LOC118947720", # CL8
  "LOC110502101", # CL9
  "nitr2" # CL9
))

dotplot_clasico_trucha <- DotPlot(datos_integrados, features = genes_clasicos, dot.scale = 4) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  labs(title = "Marcadores clásicos de cada población", x = "Genes", y = "Clústers")

# Subconjunto de linfocitos
b_t_clusters <- c(
  "0. B tempranas (IgM⁺)",
  "3. T CD8⁺ efector migratorio (IL21R⁺)",
  "4. B plasmáticas (Igλ⁺ Igκ⁺)",
  "5. B naïve (IgM⁺)",
  "6. T CD4⁺ memoria central (CCR7⁺ CD28⁺)",
  "7. T CD8⁺ efector/memoria central (FOXO1⁺ SATB1⁺ CCR7⁺)",
  "8. T CD8⁺ efector/proliferativo (MKI67⁺)",
  "9. T CD8⁺ efector terminal (PRF1⁺ GZMB⁺)",
  "16. B plasmáticas maduras (Igκ⁺ CXCR4⁺)"
)
linfocitos_subset <- subset(datos_integrados, idents = b_t_clusters)

dotplot_nuevos_trucha <- DotPlot(linfocitos_subset, features = genes_nuevos_linfocitos, dot.scale = 4, assay = "RNA") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  labs(title = "Gen más expresado por clúster", x = "Genes", y = "Clústers")

# - RidgePlot de genes más característicos
ridgeplot_trucha <- RidgePlot(datos_integrados, features = c("cd3e", "hyal3", "LOC100136017"), ncol = 3)

# - FeaturePlot de genes más característicos y no identificados
# Linf B
# cl0_feature_plot <- FeaturePlot(datos_integrados, features = c("LOC118936869", "pax-5"))
# cl4_feature_plot <- FeaturePlot(datos_integrados, features = c("LOC110505920", "LOC110495722", "LOC110500016"))
# cl5_feature_plot <- FeaturePlot(datos_integrados, features = c("LOC110533327", "LOC110527864"))
# cl16_feature_plot <- FeaturePlot(datos_integrados, features = c("LOC110500099", "LOC110499553"))

# Linf T
# cl3_feature_plot <- FeaturePlot(datos_integrados, features = c("LOC110535393", "LOC110505032"))
# cl6_feature_plot <- FeaturePlot(datos_integrados, features = c("LOC110496534", "LOC110534358", "LOC110516883"))
# cl7_feature_plot <- FeaturePlot(datos_integrados, features = c("LOC110537729", "tcf7"))
# cl8_feature_plot <- FeaturePlot(datos_integrados, features = c("s100w", "LOC118936983"))
# cl9_feature_plot <- FeaturePlot(datos_integrados, features = c("LOC110509811"))

# - ClusterTree y matriz de distancias
colnames(datos_integrados@meta.data) <- sub("^resolution_", "res.", colnames(datos_integrados@meta.data))
cluster_tree <- clustree(datos_integrados@meta.data, prefix = "res.")

# 7. Guardar los resultados ---------------------------------------------------
# Guardar en txt los marcadores más representativos (p_val) de cada clúster
top_marcadores %>%
  mutate(across(where(is.numeric), ~ round(., 5))) %>% # Redondear columnas numéricas
  write.table(file = "Resultados/Trucha_MSP/top_marcadores_p_val_adj.txt", sep = "\t", quote = FALSE, row.names = FALSE, dec = ",")

# Guardar en txt los marcadores más expresados(pct1-pct2) de cada clúster
top_marcadores <- marcadores %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(desc(pct.1 - pct.2)) %>%
  dplyr::slice(1:20) %>%
  dplyr::ungroup()

top_marcadores %>%
  mutate(across(where(is.numeric), ~ round(., 5))) %>%
  write.table(file = "Resultados/Trucha_MSP/top_marcadores_pct_diff.txt", sep = "\t", quote = FALSE, row.names = FALSE, dec = ",")

# Guardar gráficos
ggsave("Resultados/Trucha_MSP/UMAP_TSNE_clusters_trucha.png", plot = UMAP_Clusters + TSNE_Clusters, width = 20, height = 8, dpi = 500)
# ggsave("Resultados/Trucha_MSP/UMAP_TSNE_muestras_trucha.png", plot = UMAP_Muestras + TSNE_Muestras, width = 16, height = 8, dpi = 500)
ggsave("Resultados/Trucha_MSP/heatmap_trucha.png", plot = heatmap_marcadores, width = 20, height = 8, dpi = 500)
ggsave("Resultados/Trucha_MSP/dotplot_nuevos_trucha.jpg", plot = dotplot_nuevos_trucha, width = 14, height = 8, dpi = 500)
ggsave("Resultados/Trucha_MSP/dotplot_clasico_trucha.jpg", plot = dotplot_clasico_trucha, width = 14, height = 8, dpi = 500)
ggsave("Resultados/Trucha_MSP/clustertree_trucha.png", plot = cluster_tree, width = 15, height = 7, dpi = 500)
ggsave("Resultados/Trucha_MSP/elbowplot_trucha.jpg", plot = elbow_plot, width = 14, height = 8, dpi = 500)
ggsave("Resultados/Trucha_MSP/featureplot_trucha.png", plot = feature_plot, width = 14, height = 8, dpi = 500)
ggsave("Resultados/Trucha_MSP/ridgeplot_trucha.png", plot = ridgeplot_trucha, width = 20, height = 8, dpi = 500)
