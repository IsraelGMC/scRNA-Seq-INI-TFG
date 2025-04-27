# 1. Instalar y cargar la librería LoupeR, si no está instalada previamente:
# remotes::install_github("10xGenomics/loupeR")
library(loupeR)

# 2. Preparar el objeto Seurat para exportar el Assay de RNA:
# Se trabaja sobre una copia para no modificar el objeto original.
datos_integrados_loupe <- datos_integrados
DefaultAssay(datos_integrados_loupe) <- "RNA"

# 3. Crear el archivo .loupe a partir del objeto Seurat:
create_loupe_from_seurat(
    datos_integrados_loupe,
    output_name = "Resultados/Trucha_MSP/analisis_assay_RNA",
    force = TRUE # Sobrescribe si el archivo ya existe.
)
