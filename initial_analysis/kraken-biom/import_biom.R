if (!requireNamespace("BiocManager", quietly = TRUE))

# extracting tax table from kraken biom-file
  
install.packages("BiocManager")

BiocManager::install("phyloseq") # Install phyloseq

install.packages(c("RColorBrewer", "patchwork")) #install patchwork to chart publication-quality plots and readr to read rectangular datasets.


library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("patchwork")

merged_metagenomes <- import_biom("allseq.biom")
class(merged_metagenomes)
View(merged_metagenomes@tax_table@.Data)

merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4)
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

tax_data <- merged_metagenomes@tax_table@.Data
tax_df <- as.data.frame(tax_data)

write.xlsx(tax_df, "tax_id.xlsx")


