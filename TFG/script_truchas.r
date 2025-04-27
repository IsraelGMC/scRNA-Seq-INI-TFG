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
  titulos <- c("Genes por cél.", "UMIs por cél.", "Contenido mitocondrial")
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
# VizDimLoadings(datos_integrados, dims = 1:2, reduction = "pca") # Genes que más contribuyen a las componentes principales
elbow_plot <- ElbowPlot(datos_integrados, ndims = 30)

# UMAP y TSNE
datos_integrados <- RunUMAP(datos_integrados, dims = 1:30)
datos_integrados <- RunTSNE(datos_integrados, dims = 1:30)

# B. Clustering
datos_integrados <- FindNeighbors(datos_integrados, dims = 1:30)

# Se hacen 3 resoluciones para ver cuál es más adecuada
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
# Asignación de clústeres a tipo celular (reemplaza números)
new.cluster.ids <- c(
  "1. B tempranas",
  "2. Megacariocitos I",
  "3. Megacariocitos II",
  "4. T TCR-alfa",
  "5. B maduras",
  "6. B plasmablastos",
  "7. T CD4+ memoria",
  "8. T transición memoria/Treg",
  "9. T CD8+ prolif.",
  "10. T CD8+ dif.",
  "11. Eritrocitos",
  "12. Neutrófilos",
  "13. NCC",
  "14. Macrófagos I",
  "15. Macrófagos II",
  "16. cDCs",
  "17. B plasmáticas",
  "18. NK-like"
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
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 4)))

UMAP_Clusters + TSNE_Clusters
UMAP_Muestras + TSNE_Muestras

# Otros gráficos. DEBE USAR ASSAY INTEGRATED
DefaultAssay(datos_integrados) <- "integrated"

# Heatmap de los 3 marcadores diferenciales más significativos por clúster
heatmap_marcadores <- DoHeatmap(datos_integrados, features = top3_marcadores$gene, size = 3, angle = 55) +
  theme(axis.text.y = element_text(size = 6))

# Visualización de la expresión génica en clusterss
# - Marcadores representativos de cada clúster
features_trucha <- unique(c(
  "LOC110485215", # Cluster 0
  "hyal3", # Cluster 1
  "thbs1b", # Cluster 2
  "LOC110530298", # Cluster 3
  "LOC110487484", # Cluster 4
  "LOC110533327", # Cluster 5
  "LOC110534451", # Cluster 6
  "LOC110537729", # Cluster 7
  "sh2d1ab", # Cluster 8
  "il2rb", # Cluster 9
  "LOC110489254", # Cluster 10
  "LOC100136017", # Cluster 11
  "si:dkey-9i23.4", # Cluster 12
  "LOC100136950", # Cluster 13 y 14
  "LOC110535338", # Cluster 15
  "zgc:152968", # Cluster 16
  "LOC110536507" # Cluster 17
))

# - DotPlot del gen más expresado por cada clúster
dotplot_trucha <- DotPlot(datos_integrados, features = features_trucha, dot.scale = 4) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  labs(title = "Gen más expresado por clúster", x = "Genes", y = "Clústers")

# - RidgePlot de genes más característicos
ridgeplot_trucha <- RidgePlot(datos_integrados, features = c("cd3e", "hyal3", "LOC100136017"), ncol = 3)

# - FeaturePlot de genes más característicos
feature_plot <- FeaturePlot(datos_integrados, features = c("LOC100136251", "f13a1b"))

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
ggsave("Resultados/Trucha_MSP/vln_sin_QC_trucha.png", plot = pre_limpieza, width = 20, height = 8, dpi = 500)
ggsave("Resultados/Trucha_MSP/vln_con_QC_trucha.png", plot = post_limpieza, width = 20, height = 8, dpi = 500)
ggsave("Resultados/Trucha_MSP/UMAP_TSNE_clusters_trucha.png", plot = UMAP_Clusters + TSNE_Clusters, width = 20, height = 8, dpi = 500)
ggsave("Resultados/Trucha_MSP/UMAP_TSNE_muestras_trucha.png", plot = UMAP_Muestras + TSNE_Muestras, width = 16, height = 8, dpi = 500)
ggsave("Resultados/Trucha_MSP/heatmap_trucha.png", plot = heatmap_marcadores, width = 20, height = 8, dpi = 500)
ggsave("Resultados/Trucha_MSP/dotplot_trucha.png", plot = dotplot_trucha, width = 14, height = 8, dpi = 500)
ggsave("Resultados/Trucha_MSP/clustertree_trucha.png", plot = cluster_tree, width = 15, height = 7, dpi = 500)
ggsave("Resultados/Trucha_MSP/elbowplot_trucha.png", plot = elbow_plot, width = 14, height = 8, dpi = 500)
ggsave("Resultados/Trucha_MSP/featureplot_trucha.png", plot = feature_plot, width = 14, height = 8, dpi = 500)
ggsave("Resultados/Trucha_MSP/ridgeplot_trucha.png", plot = ridgeplot_trucha, width = 20, height = 8, dpi = 500)
