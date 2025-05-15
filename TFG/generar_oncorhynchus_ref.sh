# Generar un archivo GTF filtrado con biotipos de interés (protein_coding, lncRNA y pseudogene):
cellranger mkgtf \
  GCF_013265735.2_USDA_OmykA_1.1_genomic.gtf \
  filtered_genes.gtf \
  --attribute=gene_biotype:protein_coding \
  --attribute=gene_biotype:lncRNA \
  --attribute=gene_biotype:pseudogene

# 2. Filtrar para conservar únicamente anotaciones de tipo 'exon' y 'gene':
awk '$3 == "exon" || $3 == "gene"' filtered_genes.gtf > genes_and_exons.gtf

# 3. Eliminar entradas correspondientes a transcritos de tipo rRNA (ruido):
grep -v 'transcript_biotype "rRNA"' genes_and_exons.gtf > genes_and_exons_no_rRNA.gtf

# 4. Modificar el identificador gene_id para la secuencia mitocondrial para posterior filtrado:
sed '/^NC_001717.1/ s/gene_id "/gene_id "MT-/' genes_and_exons_no_rRNA.gtf > genes_and_exons_no_rRNA_no_MT.gtf

# 5. Eliminar archivos intermedios utilizados en el preprocesamiento:
rm filtered_genes.gtf genes_and_exons.gtf genes_and_exons_no_rRNA.gtf

# 6. Crear el archivo de referencia para Cell Ranger de O. mykiss:
cellranger mkref \
  --genome=oncorhynchus_mykiss_genome_ref \
  --fasta=GCF_013265735.2_USDA_OmykA_1.1_genomic.fna \
  --genes=genes_and_exons_no_rRNA_no_MT.gtf
