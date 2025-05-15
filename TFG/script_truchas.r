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
      theme(
        plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.title.x = element_text(size = 24, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(size = 24, face = "bold", margin = margin(r = 10)),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 20),
        axis.text.y = element_text(size = 20)
      )
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

#>>>>>>>>>>>>>>>>>>>>>>>>>>
# QUITAR
library(Seurat)
library(ggplot2)
library(dplyr)
library(scDblFinder)
library(clustree)
#<<<<<<<<<<<<<<<<<<<<<<<<<<

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
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 15, face = "bold"),
    legend.key.size = unit(0.8, "cm")
  )

# - RidgePlots de genes más característicos RESULTADOS
library(patchwork)

make_ridge_panel <- function(datos, genes, ncol = length(genes)) {
  plots <- RidgePlot(datos, features = genes, combine = FALSE)
  plots <- lapply(seq_along(plots), function(i) {
    p <- plots[[i]] +
      theme(
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, face = "bold")
      )
    # Mostrar eje Y solo en el primer gráfico
    if (i == 1) {
      p <- p + theme(axis.title.y = element_text(size = 30, face = "bold", hjust = 0.5))
    } else {
      p <- p + theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank()
      )
    }
    # Mostrar eje X solo en el segundo gráfico
    if (i == 2) {
      p <- p + theme(axis.title.x = element_text(size = 30, face = "bold", hjust = 0.64))
    } else {
      p <- p + theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank()
      )
    }
    # Eliminar leyenda salvo en el último (opcional)
    if (i != length(plots)) {
      p <- p + theme(legend.position = "none")
    } else {
      p <- p + theme(
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.key.size = unit(0.9, "cm")
      )
    }

    return(p)
  })
  wrap_plots(plots, ncol = ncol)
}

# Panel Linf B RESULTADOS
ridgeplot_B_trucha <- make_ridge_panel(datos_integrados,
  c("pax-5", "blnk", "LOC110537828"),
  ncol = 3
)

ridgeplot_B_trucha
# Panel Linf T RESULTADOS
ridgeplot_T_trucha <- make_ridge_panel(datos_integrados,
  c("LOC110530298", "cd3e", "tcf7"),
  ncol = 3
)

ridgeplot_B_nuevos_trucha <- make_ridge_panel(datos_integrados,
  c("LOC110533327", "si:dkey-24p1.1", "LOC110533868"),
  ncol = 3
)

# Asignación de clústeres a tipo celular (tras AED)
nuevos_id_clusteres <- c(
  "0. B naïve IgD⁺",
  "1. Megacariocitos I",
  "2. Megacariocitos II",
  "3. T CD4⁺ Th17",
  "4. B maduras Igλ⁺/Igκ⁺",
  "5. B naïve IgT⁺",
  "6. T CD4⁺ Th0 naïve",
  "7. T CD4⁺ Th0 memoria",
  "8. T CD8⁺ Tc efector proliferativo",
  "9. T CD8⁺ Tc efector terminal",
  "10. Eritrocitos",
  "11. Neutrófilos",
  "12. NCC",
  "13. Macrófagos I",
  "14. Macrófagos II",
  "15. cDCs",
  "16. B plasmáticas Igκ⁺/CXCR4⁺",
  "17. NK-Like"
)

names(nuevos_id_clusteres) <- levels(datos_integrados)
datos_integrados <- RenameIdents(datos_integrados, nuevos_id_clusteres)

# Gráficos UMAP-TSNE. Células por muestra
UMAP_Muestras <- DimPlot(
  datos_integrados,
  group.by = "orig.ident",
  reduction = "umap",
  label.size = 3.5,
  pt.size = 0.4
) +
  labs(title = "Muestras (UMAP)") +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18)
  ) +
  NoLegend()

TSNE_Muestras <- DimPlot(
  datos_integrados,
  group.by = "orig.ident",
  reduction = "tsne",
  label.size = 3.5,
  pt.size = 0.4
) +
  labs(title = "Muestras (t-SNE)") +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18)
  )

# Gráficos UMAP-TSNE. Clústeres
UMAP_Clusters <- DimPlot(
  datos_integrados,
  reduction = "umap",
  label = TRUE,
  label.size = 7,
  pt.size = 0.4,
  repel = TRUE
) +
  labs(title = "Clústeres (UMAP)") +
  theme(
    axis.title = element_text(size = 20), axis.text = element_text(size = 22),
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 4))) +
  NoLegend()

TSNE_Clusters <- DimPlot(
  datos_integrados,
  reduction = "tsne",
  label = TRUE,
  label.size = 7,
  pt.size = 0.4,
  repel = TRUE
) +
  labs(title = "Clústeres (t-SNE)") +
  theme(
    axis.title = element_text(size = 20), axis.text = element_text(size = 22),
    legend.text = element_text(size = 16),
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 4), keyheight = unit(1, "cm"), keywidth = unit(1, "cm")))

# Otros gráficos.
# - DotPlot del gen más expresado por cada clúster
genes_clasicos <- unique(c(
  "pax-5", "blnk", "LOC110537828", # General linfocitos B
  "LOC110530298", "cd3e", "cd8a", # General de linfocitos T
  "LOC110538697", "LOC110485215", # CL0: IgDs
  "hyal3", # CL1
  "f13a1b", # CL2
  "dock10", "il7r", # CL3
  "LOC110487484", "LOC110496825", "cd79b", # CL4
  "LOC110486747", # # CL5: IgT
  "LOC110534451", "LOC100136274", # CL6
  "LOC118964336", "tcf7", # CL7
  "sh2d1ab", "LOC110531827", # CL8
  "il2rb", "LOC110508453", # CL9
  "hba4", # CL10
  "mmp9", "lect2", # CL11
  "si:dkey-9i23.4", "LOC110531658", # CL12
  "mpeg1.1", "plxdc2", # CL13
  "lyz2", "csf3r", # CL14
  "LOC110535338", "LOC110487813", # CL15
  "LOC110530627", "LOC110499928", "LOC110500099", # CL16
  "LOC110498134", "prf1.3" # CL17
))

dotplot_clasico_trucha <- DotPlot(datos_integrados, features = genes_clasicos, dot.scale = 4) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, face = "bold"),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  labs(title = "Marcadores clásicos de cada población", x = "Genes", y = "Clústers")

# Subconjunto de linfocitos
genes_nuevos_linfocitos <- unique(c(
  # Linfocitos B
  "LOC118938876", # CL0; lncRNA
  "LOC110533327", # CL0, CL5; PC
  "LOC110533868", # CL4; PC
  "LOC110527864", # CL5; PC
  "si:dkey-24p1.1", # CL5; PC
  "LOC118942906", # CL16; lncRNA
  "LOC110526114", # CL16; PC
  # Linfocitos T
  "LOC118937335", # CL7; lncRNA
  "LOC110529458", # CL7; PC
  "LOC118947720", # CL8; PC
  "LOC110499816", # CL8; PC
  "LOC110502101", # CL9; PC
  # Descartados
  "LOC110495722", # CL4; lncRNA; B; NO INTERESA
  "LOC110500016", # CL4; lncRNA; B; NO INTERESA
  "LOC110528322", # CL3; T; NO INTERESA
  "si:ch211-67e16.3", # CL3; T; NO INTERESA
  "LOC110534358", # CL6; PC; T; NO INTERESA
  "LOC110496534" # CL6; lncRNA; T; NO INTERESA
))

b_t_clusters <- c(
  "0. B naïve IgD⁺",
  "3. T CD4⁺ Th17",
  "4. B maduras Igλ⁺/Igκ⁺",
  "5. B naïve IgT⁺",
  "6. T CD4⁺ Th0 naïve",
  "7. T CD4⁺ Th0 memoria",
  "8. T CD8⁺ Tc efector proliferativo",
  "9. T CD8⁺ Tc efector terminal",
  "16. B plasmáticas Igκ⁺/CXCR4⁺"
)
linfocitos_subset <- subset(datos_integrados, idents = b_t_clusters)

dotplot_nuevos_trucha <- DotPlot(linfocitos_subset, features = genes_nuevos_linfocitos, dot.scale = 4, assay = "RNA") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, face = "bold"),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  labs(title = "Gen más expresado por clúster", x = "Genes", y = "Clústers")

# - ClusterTree y matriz de distancias
colnames(datos_integrados@meta.data) <- sub("^resolution_", "res.", colnames(datos_integrados@meta.data))
cluster_tree <- clustree(datos_integrados@meta.data, prefix = "res.", node_text_size = 13, node_size_range = c(12, 30)) +
  theme(text = element_text(size = 36)) +
  scale_edge_color_continuous(low = "blue", high = "red")

# - FeaturePlot de genes más característicos ya conocidos
tema_ejes <- theme(
  axis.title = element_text(size = 20),
  axis.text = element_text(size = 18),
  plot.title = element_text(size = 25, hjust = 0.5, face = "bold")
)

# Linf B CLÁSICOS
cl0bc_feature_plot <- FeaturePlot(datos_integrados, features = "LOC110485215") + ggtitle("Clúster 0 – IgD (LOC110485215)") + tema_ejes
cl5bc_feature_plot <- FeaturePlot(datos_integrados, features = "LOC110486747") + ggtitle("Clúster 5 – IgT (LOC110486747)") + tema_ejes
cl4bc_feature_plot <- FeaturePlot(datos_integrados, features = "LOC110487484") + ggtitle("Clúster 4 – Igλ (LOC110487484)") + tema_ejes
cl16bc_feature_plot <- FeaturePlot(datos_integrados, features = "LOC110530627") + ggtitle("Clúster 16 – CXCR4 (LOC110530627)") + tema_ejes

grid_B <- (cl0bc_feature_plot | cl4bc_feature_plot) / (cl5bc_feature_plot | cl16bc_feature_plot)

# Linf T CLÁSICOS
cl9tc_feature_plot <- FeaturePlot(datos_integrados, features = "LOC110520645") + ggtitle("Clúster 9 – LOC110520645") + tema_ejes
cl7tc_feature_plot <- FeaturePlot(datos_integrados, features = "tcf7") + ggtitle("Clúster 7 – tcf7") + tema_ejes
cl6tc_feature_plot <- FeaturePlot(datos_integrados, features = "LOC110489493") + ggtitle("Clúster 6 – LOC110489493") + tema_ejes
cl8tc_feature_plot <- FeaturePlot(datos_integrados, features = "mki67") + ggtitle("Clúster 8 – MKI67") + tema_ejes

grid_T <- (cl6tc_feature_plot | cl7tc_feature_plot) / (cl8tc_feature_plot | cl9tc_feature_plot)

# Linf B NUEVOS
cl0bn_feature_plot <- FeaturePlot(datos_integrados, features = "LOC118938876") + ggtitle("Clúster 0 – LOC118938876") + tema_ejes
cl4bn_feature_plot <- FeaturePlot(datos_integrados, features = "LOC110533868") + ggtitle("Clúster 4 – LOC110533868") + tema_ejes
cl5bn_feature_plot <- FeaturePlot(datos_integrados, features = "LOC110527864") + ggtitle("Clúster 5 – LOC110527864") + tema_ejes
cl16bn_feature_plot <- FeaturePlot(datos_integrados, features = "LOC110526114") + ggtitle("Clúster 16 – LOC110526114") + tema_ejes

grid_B_nuevos <- (cl0bn_feature_plot | cl4bn_feature_plot) / (cl5bn_feature_plot | cl16bn_feature_plot)

# Linf T NUEVOS
cl7tn_1_feature_plot <- FeaturePlot(datos_integrados, features = "LOC118937335") + ggtitle("Clúster 7 – LOC118937335") + tema_ejes
cl7tn_2_feature_plot <- FeaturePlot(datos_integrados, features = "LOC110490984") + ggtitle("Clúster 7 – FOXO1 (LOC110490984)") + tema_ejes
cl8_tn_feature_plot <- FeaturePlot(datos_integrados, features = "LOC110499816") + ggtitle("Clúster 8 – AURKB (LOC110499816)") + tema_ejes
cl9tn_feature_plot <- FeaturePlot(datos_integrados, features = "LOC110502101") + ggtitle("Clúster 9 – DGKB (LOC110502101)") + tema_ejes

grid_T_nuevos <- (cl7tn_1_feature_plot | cl7tn_2_feature_plot) / (cl8_tn_feature_plot | cl9tn_feature_plot)

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
ggsave("Resultados/Trucha_MSP/UMAP_TSNE_muestras_trucha.png", plot = UMAP_Muestras + TSNE_Muestras, width = 16, height = 8, dpi = 500)
ggsave("Resultados/Trucha_MSP/heatmap_trucha.png", plot = heatmap_marcadores, width = 20, height = 12, dpi = 500)
ggsave("Resultados/Trucha_MSP/dotplot_Clásicos_trucha.jpg", plot = dotplot_clasico_trucha, width = 14, height = 8, dpi = 500)
ggsave("Resultados/Trucha_MSP/dotplot_Nuevos_trucha.jpg", plot = dotplot_nuevos_trucha, width = 14, height = 6, dpi = 500)
ggsave("Resultados/Trucha_MSP/clustertree_trucha.png", plot = cluster_tree, width = 30, height = 11, dpi = 500)
ggsave("Resultados/Trucha_MSP/elbowplot_trucha.jpg", plot = elbow_plot, width = 14, height = 8, dpi = 500)
ggsave("Resultados/Trucha_MSP/ridgeplot_T_trucha.png", plot = ridgeplot_T_trucha, width = 20, height = 8, dpi = 500)
ggsave("Resultados/Trucha_MSP/ridgeplot_B_trucha.png", plot = ridgeplot_B_trucha, width = 20, height = 8, dpi = 500)
ggsave("Resultados/Trucha_MSP/ridgeplot_B_Nuevos_trucha.png", plot = ridgeplot_B_nuevos_trucha, width = 20, height = 8, dpi = 500)
ggsave("Resultados/Trucha_MSP/featureplot_B_Clásicos_trucha.png", plot = grid_B, width = 14, height = 14, dpi = 500)
ggsave("Resultados/Trucha_MSP/featureplot_T_Clásicos_trucha.png", plot = grid_T, width = 14, height = 14, dpi = 500)
ggsave("Resultados/Trucha_MSP/featureplot_B_Nuevos_trucha.png", plot = grid_B_nuevos, width = 14, height = 14, dpi = 500)
ggsave("Resultados/Trucha_MSP/featureplot_T_Nuevos_trucha.png", plot = grid_T_nuevos, width = 14, height = 14, dpi = 500)
ggsave("Resultados/Trucha_MSP/vln_coe1a.png", plot = coe1a_vlnplot, width = 14, height = 8, dpi = 500)
