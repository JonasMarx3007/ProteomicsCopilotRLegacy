library(ggvenn)
library(VennDiagram)

xy_intensity_plot <- function(data, meta, sample1, sample2, workflow = "Protein") {
  if (workflow == "Protein") {
    id_col <- "ProteinNames"
    id_col_name <- "Protein"
  } else if (workflow == "Phosphosite") {
    id_col <- "PTM_Collapse_key"
    id_col_name <- "Phossite"
  } else {
    stop("Invalid workflow specified.")
  }
  
  all_columns <- c(sample1, sample2)
  plot_data <- data[, c(id_col, all_columns), drop = FALSE]
  plot_data <- plot_data[!is.na(plot_data[[sample1]]) & !is.na(plot_data[[sample2]]), ]
  
  cor_val <- cor(plot_data[[sample1]], plot_data[[sample2]], method = "pearson")
  cor_text <- paste0("Pearson r = ", round(cor_val, 3))
  
  p <- ggplot(plot_data, aes(x = !!sym(sample1), y = !!sym(sample2))) +
    geom_point(color = "blue", alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(
      x = paste(sample1, "Intensity"),
      y = paste(sample2, "Intensity"),
      title = paste(id_col_name, "Intensity Scatter Plot"),
      subtitle = cor_text
    ) +
    theme_minimal()
  
  return(p)
}
venn_plot2 <- function(data, meta, sample1, sample2, workflow = "Protein") {
  
  if (workflow == "Protein") {
    id_col <- "ProteinNames"
  } else if (workflow == "Phosphosite") {
    id_col <- "PTM_Collapse_key"
  } else {
    stop("Invalid workflow specified.")
  }
  
  plot_data <- data[, c(id_col, sample1, sample2), drop = FALSE]
  
  present1 <- plot_data[[id_col]][!is.na(plot_data[[sample1]])]
  present2 <- plot_data[[id_col]][!is.na(plot_data[[sample2]])]
  
  all_present <- union(present1, present2)
  present1 <- intersect(present1, all_present)
  present2 <- intersect(present2, all_present)
  
  list_data <- list(
    sample1 = present1,
    sample2 = present2
  )
  
  ggvenn(list_data, fill_color = c("lightblue", "lightgreen"), stroke_size = 0.5, set_name_size = 5)
}
venn_plot <- function(data, meta, sample1, sample2, workflow = "Protein") {
  
  if (workflow == "Protein") {
    id_col <- "ProteinNames"
    id_col_name <- "Protein"
  } else if (workflow == "Phosphosite") {
    id_col <- "PTM_Collapse_key"
    id_col_name <- "Phossite"
  } else {
    stop("Invalid workflow specified.")
  }
  
  plot_data <- data[, c(id_col, sample1, sample2), drop = FALSE]
  
  present1 <- plot_data[[id_col]][!is.na(plot_data[[sample1]])]
  present2 <- plot_data[[id_col]][!is.na(plot_data[[sample2]])]
  
  all_present <- union(present1, present2)
  present1 <- intersect(present1, all_present)
  present2 <- intersect(present2, all_present)
  
  venn.plot <- draw.pairwise.venn(
    area1 = length(present1),
    area2 = length(present2),
    cross.area = length(intersect(present1, present2)),
    category = c(sample1, sample2),
    fill = c("lightblue", "lightgreen"),
    alpha = 0.5,
    cex = 2,
    cat.cex = 1.5,
    cat.pos = c(-20, 20)
  )
  
  return(venn.plot)
}
phossite_per_peptide <- function(data) {
  if (!"Modified.Sequence" %in% colnames(data)) {
    stop("Data must have a 'Modified.Sequence' column")
  }
  
  unique_sequences <- unique(na.omit(data$Modified.Sequence))
  
  phospho_counts <- sapply(unique_sequences, function(seq) {
    str_count(seq, fixed("[Phospho (STY)]"))
  })
  
  df_counts <- data.frame(PhosphoCount = phospho_counts)
  
  ggplot(df_counts, aes(x = PhosphoCount)) +
    geom_histogram(
      binwidth = 1,
      fill = "skyblue",
      color = "black",
      boundary = -0.5
    ) +
    scale_x_continuous(breaks = 0:max(df_counts$PhosphoCount)) +
    labs(
      title = "Distribution of Phosphosites per Unique Sequence",
      x = "Number of Phosphosites",
      y = "Frequency"
    ) +
    theme_minimal(base_size = 14)
}
phossite_per_peptide_stacked <- function(data) {
  if (!"Modified.Sequence" %in% colnames(data)) {
    stop("Data must have a 'Modified.Sequence' column")
  }
  
  unique_sequences <- unique(na.omit(data$Modified.Sequence))
  
  phospho_counts <- sapply(unique_sequences, function(seq) {
    str_count(seq, fixed("[Phospho (STY)]"))
  })
  
  df_summary <- tibble(PhosphoCount = phospho_counts) %>%
    group_by(PhosphoCount) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(Fraction = n / sum(n) * 100) %>%
    mutate(PhosphoCount = factor(PhosphoCount,
                                 levels = sort(unique(PhosphoCount), decreasing = FALSE)))
  
  ggplot(df_summary, aes(x = "", y = Fraction, fill = PhosphoCount)) +
    geom_bar(stat = "identity", width = 0.5) +
    scale_fill_brewer(palette = "Blues", name = "Phospho Sites") +
    labs(
      title = "Fraction of Phosphosites per Unique Sequence",
      x = "",
      y = "Percentage (%)"
    ) +
    theme_minimal(base_size = 14)
}
