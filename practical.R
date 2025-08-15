#### Protein ####
setwd("C:/Users/Jonas Marx/Desktop/Proteomics CopilotR")
source("functions.R")
source("practical_functions.R")
setwd("C:/Users/Jonas Marx/Desktop/Praktikum")
data = read_data("20250814_150419_20250811_Praktikum_CA_WP_PG.tsv")
data2 = read_data("20250815_125230_20250811_Praktikum_CA_WP_new.tsv")
meta = read_data("metadata-2025-08-15.csv")
org_data = inverseof_log2_transform_data(data, meta)

#3.3.1
proteins = length(unique(data$ProteinNames))
print(paste0("There are ", proteins, " unique proteins!"))

peptides = length(unique(data2$Stripped.Sequence))
print(paste0("There are ", peptides, " unique peptides!"))

precursor = length(unique(data2$EG.PrecursorId))
print(paste0("There are ", precursor, " unique precursors!"))

#3.3.2
boxplot_int_single(data, meta, id=F)
cov_plot(org_data, meta)

#3.3.3
print(meta$sample)
xy_intensity_plot(data, meta, "[1] 20250811_nElute_tims1_YY_Praktikum_GroupA_WP_1_Slot1-01_1_11285.htrms.PG.Log2Quantity", "[2] 20250811_nElute_tims1_YY_Praktikum_GroupA_WP_2_Slot1-02_1_11286.htrms.PG.Log2Quantity")

#3.3.4
corr_plot(data, meta, id=F, method = T, full_range = T)

#3.3.5
venn_plot(data, meta, "[1] 20250811_nElute_tims1_YY_Praktikum_GroupA_WP_1_Slot1-01_1_11285.htrms.PG.Log2Quantity", "[2] 20250811_nElute_tims1_YY_Praktikum_GroupA_WP_2_Slot1-02_1_11286.htrms.PG.Log2Quantity")
venn_plot2(data, meta, "[1] 20250811_nElute_tims1_YY_Praktikum_GroupA_WP_1_Slot1-01_1_11285.htrms.PG.Log2Quantity", "[2] 20250811_nElute_tims1_YY_Praktikum_GroupA_WP_2_Slot1-02_1_11286.htrms.PG.Log2Quantity")

#3.3.6
filtered_data = filter_data(data, meta, 2)
imputed_data = impute_values(filtered_data, meta)
pca_plot(imputed_data, meta)

#3.3.7
volcano_plot(imputed_data, meta, "A", "C")

#3.3.8 -> EMAIL

#3.3.9
diff_genes = different_genes(imputed_data, meta, "A", "C")
diff_genes$Upregulated

enrichment_analysis(diff_genes$Upregulated)



#### Phospho ####
data3 = read_data("20250815_121046_HeLa-EGF_PP_Swiss-EMBL_collapse.tsv")
data4 = read_data("20250815_122110_Mono2_5µg_Report.tsv")

#3.4.1
proteins = length(unique(data4$ProteinNames))
print(paste0("There are ", proteins, " unique proteins!"))

peptides = length(unique(data4$Stripped.Sequence))
print(paste0("There are ", peptides, " unique peptides!"))

precursor = length(unique(data4$EG.PrecursorId))
print(paste0("There are ", precursor, " unique precursors!"))

#3.4.2
phossite_per_peptide_stacked(data4)

#3.4.3
phossite_per_peptide(data4)

#3.4.4 -> Yue

#3.4.5 -> Yue

#3.4.6 -> Copilot

