# -------------------------------------------------------------
# Análisis de transcriptómica de célula única (scRNA-seq) - Pulpo
# -------------------------------------------------------------

# 0. Cargar bibliotecas
library(Seurat)
library(ggplot2)
library(dplyr)
library(scDblFinder)

# 1. Cargar datos y crear objeto Seurat ------------------------
pulpo_data <- Read10X(data.dir = "Datos/Pulpo/filtered_feature_bc_matrix")
pulpo_obj <- CreateSeuratObject(
  counts = pulpo_data,
  project = "Pulpo_scRNAseq",
  min.cells = 3,
  min.features = 200
)

# 2. Control de calidad (QC) ------------------------------------

# Calcular porcentaje de genes mitocondriales
pulpo_obj[["percent.mt"]] <- PercentageFeatureSet(pulpo_obj, pattern = "^MT-")

# Visualizar métricas básicas de QC
VlnPlot(pulpo_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filtrar células basadas en umbrales de calidad
pulpo_obj <- subset(pulpo_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Detección de dobles usando scDblFinder
pulpo_sce <- as.SingleCellExperiment(pulpo_obj)
pulpo_sce <- scDblFinder(pulpo_sce)
pulpo_obj <- AddMetaData(pulpo_obj, colData(pulpo_sce)[, "scDblFinder.class"], col.name = "scDblFinder.class")

# Filtrar dobles detectados
pulpo_obj <- subset(pulpo_obj, scDblFinder.class == "singlet")

# Visualización posterior al filtrado
VlnPlot(pulpo_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 3. Normalización y selección de genes variables --------------

# Normalizar datos (LogNormalize, escala de 10,000)
pulpo_obj <- NormalizeData(pulpo_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Identificar los 2500 genes más variables
pulpo_obj <- FindVariableFeatures(pulpo_obj, selection.method = "vst", nfeatures = 2500)

# 4. Análisis exploratorio básico --------------------------------

# Escalar datos
all_genes <- rownames(pulpo_obj)
pulpo_obj <- ScaleData(pulpo_obj, features = all_genes)

# Reducción de dimensionalidad
pulpo_obj <- RunPCA(pulpo_obj)
VizDimLoadings(pulpo_obj, dims = 1:2, reduction = "pca")
pulpo_obj <- RunUMAP(pulpo_obj, dims = 1:30)
pulpo_obj <- RunTSNE(pulpo_obj, dims = 1:30)

# Agrupamiento (clustering)
pulpo_obj <- FindNeighbors(pulpo_obj, dims = 1:30)
pulpo_obj <- FindClusters(pulpo_obj, resolution = 0.5)

# 5. Identificación de marcadores --------------------------------

# Buscar genes marcadores de cada clúster
marcadores <- FindAllMarkers(
  pulpo_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Seleccionar los 10 mejores marcadores por clúster (por p-valor)
top_marcadores <- marcadores %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(p_val_adj) %>%
  dplyr::slice(1:10) %>%
  dplyr::ungroup()

# 6. Visualización de resultados -------------------------------

# Gráficos de clústeres en UMAP y t-SNE
UMAP_Clusters <- DimPlot(pulpo_obj, label = TRUE, reduction = "umap", pt.size = 0.5) +
  ggtitle("Clústeres de células (UMAP)") +
  theme(
    axis.title = element_text(size = 20), axis.text = element_text(size = 18),
    legend.title = element_text(size = 22), legend.text = element_text(size = 18)
  )

TSNE_Clusters <- DimPlot(pulpo_obj, label = TRUE, reduction = "tsne", pt.size = 0.5) +
  ggtitle("Clústeres de células (t-SNE)") +
  theme(
    axis.title = element_text(size = 20), axis.text = element_text(size = 18),
    legend.title = element_text(size = 22), legend.text = element_text(size = 18)
  )

# Visualización conjunta
UMAP_Clusters + TSNE_Clusters

# Heatmap de los 2 mejores marcadores por clúster
top2_marcadores <- top_marcadores %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice(1:2) %>%
  dplyr::ungroup()

subset_clusters <- subset(pulpo_obj, idents = 3:18)

heatmap_marcadores <- DoHeatmap(
  subset_clusters,
  features = top2_marcadores$gene,
  size = 6,
  angle = 55
) +
  theme(
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.key.size = unit(0.8, "cm")
  )

heatmap_marcadores

# DotPlot de genes seleccionados
features_pulpo <- unique(c(
  "OCTVUL-1NCB119692", "OCTVUL-1B007148", "OCTVUL-1B029596",
  "OCTVUL-1B008669", "OCTVUL-1B017546", "OCTVUL-1B030290",
  "OCTVUL-1NCB047077", "OCTVUL-1B026250", "OCTVUL-1B015187",
  "OCTVUL-1B028951", "OCTVUL-1B001659", "OCTVUL-1B019827",
  "OCTVUL-1B030144", "OCTVUL-1B029021", "OCTVUL-1B022559",
  "OCTVUL-1B019811", "OCTVUL-1B029953", "OCTVUL-1B003402",
  "OCTVUL-1B012201", "OCTVUL-1B016834", "OCTVUL-1B023464",
  "OCTVUL-1B031523"
))

dotplot_pulpo <- DotPlot(pulpo_obj, features = features_pulpo, dot.scale = 4) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  labs(title = "Expresión de genes representativos por clúster", x = "Genes", y = "Clústeres")

# 7. Guardado de resultados -------------------------------------

# Guardar tablas de marcadores
top_marcadores %>%
  mutate(across(where(is.numeric), ~ round(., 5))) %>%
  write.table(file = "Resultados/Pulpo/top_marcadores_p_val_adj.txt", sep = "\t", quote = FALSE, row.names = FALSE, dec = ",")

top_marcadores_pct <- marcadores %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(desc(pct.1 - pct.2)) %>%
  dplyr::slice(1:20) %>%
  dplyr::ungroup()

top_marcadores_pct %>%
  mutate(across(where(is.numeric), ~ round(., 5))) %>%
  write.table(file = "Resultados/Pulpo/top_marcadores_pct_diff.txt", sep = "\t", quote = FALSE, row.names = FALSE, dec = ",")

# Guardar gráficos
ggsave("Resultados/Pulpo/umap_pulpo.png", plot = UMAP_Clusters, width = 8, height = 7, dpi = 500)
ggsave("Resultados/Pulpo/tsne_pulpo.png", plot = TSNE_Clusters, width = 8, height = 7, dpi = 500)
ggsave("Resultados/Pulpo/heatmap_pulpo.png", plot = heatmap_marcadores, width = 24, height = 10, dpi = 500)
# ggsave("Resultados/Pulpo/dotplot_pulpo.png", plot = dotplot_pulpo, width = 14, height = 6, dpi = 500)
