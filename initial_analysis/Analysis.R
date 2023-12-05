##########################################################
###############****METAGENOMICS****#######################
##########################################################
# Packages -----------------------------------------------------------------

# a lot of extra package that aren't used included here

pacman::p_load(
  rio,               # for importing data
  here,              # for relative filepaths
  skimr,             # for reviewing data
  janitor,           # for cleaning data
  lubridate,         # for date cleaning
  epikit,            # for creating age categories
  gtsummary,         # for creating tables
  tidyverse,         # for data management and visualization
  openxlsx,          # load excel files
  flextable,         # for making pretty tables
  gtsummary,         # creating tables
  scales,            # percents in tables
  ggExtra,
  gghighlight,
  RColorBrewer,
  viridis,           # for color-blind
  scales,
  apyramid,          # for age/sex pyramids
  tsibble,            # for epiweeks and time series
  ggplot2,
  RColorBrewer,
  patchwork,
  vegan,
  devtools
)

# Imports -----------------------------------------------------------------

# note - some files were manually edited using excel (i.e. deleting tax_id label from tables) to fit phyloseq requirements

# import data required by phyloseq

# sample read counts
otu <- read.table(file = "final_data/collated_FINAL.tsv", sep = "\t", header = T, row.names = 1, comment.char = "")

# taxonomy information
taxonomy <- read.table(file = "final_data/filtered_FINAL.tsv", sep = "\t", header = T, row.names = 1)

# metadata
metadata <- import(here("final_data", "metadata_mNGS.xlsx"))

# transpose lab_id from column to row. if not done, phyloseq errors
row.names(metadata) <- metadata$lab_id
metadata <- metadata %>% select (-lab_id)

# do not call at the start because it masks import from rio
library("phyloseq")

SAMPLE = sample_data(metadata)
OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxonomy))

analysis <- phyloseq(OTU, TAX, SAMPLE)


# Read Counts ----------------------------------------------------------------

sample_sums(analysis)
sort(sample_sums(analysis))

hist(sample_sums(analysis), main="Histogram: Read Counts", xlab="Total Reads", 
     border="blue", col="green", las=1, breaks=50) 

# Filtering ---------------------------------------------------------------

# remove samples with less than 20K reads

filtering_20K = prune_samples(sample_sums(analysis) > 20000, analysis)
sort(sample_sums(filtering_20K))

filtering_20K

# remove taxa with less than 1 count to resolve error for rarefaction curve generation

filtering_20K_prune = prune_taxa(taxa_sums(filtering_20K) > 1, filtering_20K)
filtering_20K_prune

# rarefaction curve

otu.rare = otu_table(filtering_20K_prune)
otu.rare = as.data.frame(t(otu.rare))
sample_names = rownames(otu.rare)

otu.rarecurve = rarecurve(otu.rare, step = 1000, label = F)

# normalization via subsampling to 20k reads

subset <- rarefy_even_depth(filtering_20K, sample.size = 20000, rngseed = 123, replace = FALSE)
sample_sums(subset)

# rarefaction curve

otu.rare = otu_table(subset)
otu.rare = as.data.frame(t(otu.rare))
sample_names = rownames(otu.rare)

otu.rarecurve = rarecurve(otu.rare, step = 100, label = T)
otu.rarecurve

# Diversity Analysis ----------------------------------------------------------

# alpha diversity

# plot all alpha diversity indices
rich = estimate_richness(subset)
rich
plot_richness(subset, x="virstrain_analysis")

# Chao1 & Shannon

(p_alpha <- plot_richness(subset, x = "virstrain_analysis", color = "virstrain_analysis", measures = c("Chao1", "Shannon")))

p_alpha + geom_violin()

p_alpha + geom_boxplot()

pairwise.wilcox.test(rich$Shannon, sample_data(subset)$virstrain_analysis)

# beta diversity

subset_nmds_bray <- ordinate(subset, "NMDS", "bray")
sample_data(subset)['sample_id'] <- row.names(sample_data(subset)) 

plot_ordination(subset, subset_nmds_bray, type="samples", color="virstrain_analysis", label="sample_id") + geom_point(size = 3)

# all methods beta diversity

library("plyr")

dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "NMDS", "MDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, subset, dist){
  ordi = ordinate(subset, method=i, distance=dist)
  plot_ordination(subset, ordi, "samples", color="virstrain_analysis")
}, subset, dist)

names(plist) <- ord_meths

pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"

p_beta = ggplot(pdataframe, aes(Axis_1, Axis_2, color=virstrain_analysis, fill=virstrain_analysis))
p_beta = p_beta + geom_point(size=4)
p_beta = p_beta + facet_wrap(~method, scales="free")
p_beta = p_beta + scale_fill_brewer(type="qual", palette="Set1")
p_beta = p_beta + scale_colour_brewer(type="qual", palette="Set1")
p_beta


# Taxonomy visualization ----------------------------------------------------------------

plot_bar(subset, fill="Kingdom") +
  facet_wrap(~ virstrain_analysis, scales = "free_x")

plot_bar(subset, fill="Phylum") +
  facet_wrap(~ virstrain_analysis, scales = "free_x")

# Agglomerate bacterial reads that do not have deeper classification

subset_glom <- tax_glom(physeq = subset, taxrank = "Phylum")
plot_bar(subset_glom, fill="Phylum") +
  facet_wrap(~ virstrain_analysis, scales = "free_x")

subset_glom <- tax_glom(physeq = subset, taxrank = "Genus")

# filter taxa where the count of values greater than 20 is more than 20% of the length of the data in consideration

subset_fil = filter_taxa(subset_glom, function(x) sum(x > 20) > (0.2*length(x)), TRUE)

# bacteroidota

bacteroi = subset_taxa(subset_fil, Phylum=="Bacteroidota")
tax_table(bacteroi)
title = "Bacteroidota"
plot_bar(bacteroi, "virstrain_analysis", "Abundance", title=title)

plot_bar(bacteroi, "virstrain_analysis", "Abundance", "Genus", title=title)

plot_bar(bacteroi, "Genus", "Abundance", "Genus", 
         title=title, facet_grid="virstrain_analysis~.")

# firmicutes

firmi = subset_taxa(subset_fil, Phylum=="Firmicutes")
title = "Firmicutes"
plot_bar(firmi, "virstrain_analysis", "Abundance", title=title)

plot_bar(firmi, "virstrain_analysis", "Abundance", "Genus", title=title)

plot_bar(firmi, "Genus", "Abundance", "Genus", 
         title=title, facet_grid="virstrain_analysis~.")

# proteobacteria

proteo = subset_taxa(subset_fil, Phylum=="Proteobacteria")

title = "Proteobacteria"
plot_bar(proteo, "virstrain_analysis", "Abundance", title=title)

plot_bar(proteo, "virstrain_analysis", "Abundance", "Genus", title=title)

plot_bar(proteo, "Genus", "Abundance", "Genus", 
         title=title, facet_grid="virstrain_analysis~.")

# actinobacteria

actino = subset_taxa(subset_fil, Phylum=="Actinobacteria")

title = "Actinobacteria"
plot_bar(actino, "virstrain_analysis", "Abundance", title=title)

plot_bar(actino, "virstrain_analysis", "Abundance", "Genus", title=title)

plot_bar(actino, "Genus", "Abundance", "Genus", 
         title=title, facet_grid="virstrain_analysis~.")


# fusobacteria

fusi = subset_taxa(subset_fil, Phylum=="Fusobacteria")

title = "Fusobacteria"
plot_bar(fusi, "virstrain_analysis", "Abundance", title=title)

plot_bar(fusi, "virstrain_analysis", "Abundance", "Genus", title=title)

plot_bar(fusi, "Genus", "Abundance", "Genus", 
         title=title, facet_grid="virstrain_analysis~.")



# Heat Map Trials ----------------------------------------------------------------

subset_bac <- subset_taxa(subset, Kingdom=="Bacteria")
subset_bac <- prune_taxa(names(sort(taxa_sums(subset),TRUE)[1:300]), subset)


plot_heatmap(subset_bac,"NMDS", "bray", "virstrain_analysis", "Phylum") + facet_wrap(~ virstrain_analysis, scales = "free_x")

plot_heatmap(subset_bac, "NMDS", "bray","virstrain_analysis", "Family", low="#66CCFF", high="#000033", na.value="white") + facet_wrap(~ virstrain_analysis, scales = "free_x")


subset_glom <- tax_glom(physeq = subset, taxrank = "Genus")
subset_fil = filter_taxa(subset_glom, function(x) sum(x > 3) > (0.2*length(x)), TRUE)


plot_heatmap(subset_fil,"NMDS", "bray", "virstrain_analysis", "Phylum") + facet_wrap(~ virstrain_analysis, scales = "free_x")

plot_heatmap(subset_fil, "NMDS", "bray","virstrain_analysis", "Genus", low="#FFFFCC", high="#000033", na.value="white") + facet_wrap(~ virstrain_analysis, scales = "free_x")


heatmap(otu_table(subset))
