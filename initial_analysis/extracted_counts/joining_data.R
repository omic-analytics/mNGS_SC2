
# getting taxa information from kraken-biom -------------------------------


# tax_id from kraken-biom - manually filtered via excel
# tax_base from recentrifuge (tax with hits)

tax_table <- import("tax_id.xlsx")
tax_base <- import("tax_base.xlsx")

# change class
tax_base <- tax_base %>%
  mutate(tax_id = as.numeric(tax_id))

# remove taxa from tax table (total taxa), and give taxa info to recentrifuge counts
filtered_taxon <- left_join(tax_base, tax_table, by = "tax_id")
write.xlsx(filtered_taxon, "filtered_taxon.xlsx")


# remove taxa rows with NA values
filtered_taxon <- import("filtered_taxon.xlsx")
filtered_tax <- filtered_taxon %>%
  filter(rowSums(!is.na(.[,-1])) > 0)
write.xlsx(filtered_tax, "filtered_FINAL.xlsx")


# joining data ------------------------------------------------------------


# Joining all recentrifuge counts from sequencing batches

jjr5b <- import("JJR5B.rcf.xlsx")
jk6n9 <- import("JK6N9.rcf.xlsx")
jk59b <- import("JK59B.rcf.xlsx")
k7983 <- import("K7983.rcf.xlsx")
kdf4y <- import("KDF4Y.rcf.xlsx")

colnames(jjr5b)[1] <- "taxa_id"
colnames(jk6n9)[1] <- "taxa_id"
colnames(jk59b)[1] <- "taxa_id"
colnames(k7983)[1] <- "taxa_id"
colnames(kdf4y)[1] <- "taxa_id"

full_join <- full_join(
  jjr5b, jk6n9, by = c("taxa_id")
)

full_join <- full_join(
  full_join, jk59b, by = c("taxa_id")
)

full_join <- full_join(
  full_join, k7983, by = c("taxa_id")
)

full_join <- full_join(
  full_join, kdf4y, by = c("taxa_id")
)

write.xlsx(full_join, "collated.xlsx")


full_join <- full_join %>%
  clean_names()

colnames(full_join)[1] <- "taxa_id"

# remove all rows with 0 value
join_trimmed <- full_join %>%
  filter_if(is.numeric, any_vars(!is.na(.) & . != 0)) %>% 
  mutate_all(~ifelse(is.na(.), 0, .))

write.xlsx(join_trimmed, "collated.xlsx")

# further filtering -------------------------------------------------------

# tax_base2 - extracted taxa from filtered_Final
# filtering taxa from collated counts based on the taxa in filterted_Final
# removing unclassified

tax_base2 <- import("tax_base2.xlsx")
all_samples <- import("collated.xlsx")
colnames(all_samples)[1] <- "tax_id"

all_samples <- all_samples %>%
  mutate(tax_id = as.numeric(tax_id))

collated_FINAL <- left_join(tax_base2, all_samples, by = "tax_id")

write.xlsx(collated_FINAL, "collated_FINAL.xlsx")
