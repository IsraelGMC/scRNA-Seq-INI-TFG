# Generar un archivo GTF filtrado con biotipos de interés (protein_coding, lncRNA y pseudogene)
cellranger mkgtf GCA_951406725.2_xcOctVulg1.2_genomic.gtf.gz filtered_genes.gtf \
  --attribute=gene_biotype:protein_coding \
  --attribute=gene_biotype:lncRNA \
  --attribute=gene_biotype:pseudogene

# Filtrar para conservar únicamente anotaciones de tipo exón
awk '$3 == "exon"' filtered_genes.gtf > genes_and_exons.gtf

# Eliminar entradas correspondientes a transcritos de tipo rRNA (ruido)
grep -v 'transcript_biotype "rRNA"' genes_and_exons.gtf > genes_and_exons_no_rRNA.gtf

# Modificar el identificador gene_id para la secuencia mitocondrial para posterior filtrado
sed '/^NC_001717.1/ s/gene_id "/gene_id "MT-/' genes_and_exons_no_rRNA.gtf > genes_and_exons_no_rRNA_no_MT.gtf

# Eliminar archivos intermedios utilizados en el preprocesamiento
rm filtered_genes.gtf genes_and_exons.gtf genes_and_exons_no_rRNA.gtf

# Crear el archivo de referencia para Cell Ranger de O. vulgaris
cellranger mkref \
  --genome=oncorhynchus_mykiss_genome_ref \
  --fasta=GCA_951406725.2_xcOctVulg1.2_genomic.fna \
  --genes=genes_and_exons_no_rRNA_no_MT.gtf
