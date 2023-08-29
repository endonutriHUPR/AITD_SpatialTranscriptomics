# ============================================================================= #
# ============================================================================= #
############################### QC metric summary ###############################
# ============================================================================= #
# ============================================================================= #

# ============================================================================= #
############################### 0. Initialization ###############################
# ============================================================================= #

########### A) Set working directory

setwd("~/Spatial_gene_expression/")

# New working directory
dir.create("./QC_sequencing")
setwd("./QC_sequencing")


########### B) Libraries

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(readr)
library(webshot)
library(gt)


########### C) Directories where the metrics are (4). Output from Spaceranger

# We load those .csv to a variable

GD2 <- as.data.frame(read_csv("~/Spatial_gene_expression/Mapping/GD2/outs/metrics_summary.csv"))
C1 <- as.data.frame(read_csv("~/Spatial_gene_expression/Mapping/C1/outs/metrics_summary.csv"))
GD3 <- as.data.frame(read_csv("~/Spatial_gene_expression/Mapping/GD3/outs/metrics_summary.csv"))
C2 <- as.data.frame(read_csv("~/Spatial_gene_expression/Mapping/C2/outs/metrics_summary.csv"))
GD1 <- as.data.frame(read_csv("~/Spatial_gene_expression/Mapping/GD1/outs/metrics_summary.csv"))
HT1 <- as.data.frame(read_csv("~/Spatial_gene_expression/Mapping/HT1/outs/metrics_summary.csv"))
HT2 <- as.data.frame(read_csv("~/Spatial_gene_expression/Mapping/HT2/outs/metrics_summary.csv"))
HT3 <- as.data.frame(read_csv("~/Spatial_gene_expression/Mapping/HT3/outs/metrics_summary.csv"))



# ============================================================================= #
############################ 1. Some transformations ############################
# ============================================================================= #


########### A) Join the samples into the same dataset

metrics <- do.call(rbind, list(GD2,C2,GD3,C1,GD1, HT1,HT2,HT3))


########### B) Change column names

colnames(metrics) <- stringr::str_replace_all(string = colnames(metrics),
                                              pattern = " ",
                                              replacement = "_")

########### C) Change percentages

# We transform percentage from x1 by x100
metrics$Fraction_of_Spots_Under_Tissue <- metrics$Fraction_of_Spots_Under_Tissue * 100
metrics$Fraction_Reads_in_Spots_Under_Tissue <- metrics$Fraction_Reads_in_Spots_Under_Tissue * 100
metrics$Reads_Mapped_Confidently_to_Exonic_Regions <- metrics$Reads_Mapped_Confidently_to_Exonic_Regions * 100
metrics$Reads_Mapped_Confidently_to_Intergenic_Regions <- metrics$Reads_Mapped_Confidently_to_Intergenic_Regions * 100
metrics$Reads_Mapped_Confidently_to_Intronic_Regions <- metrics$Reads_Mapped_Confidently_to_Intronic_Regions * 100
metrics$Reads_Mapped_Confidently_to_Genome <- metrics$Reads_Mapped_Confidently_to_Genome *100

# ============================================================================= #
#################################### 2. Tables ##################################
# ============================================================================= #


top_metrics_gex <- metrics[, c("Sample_ID",
                               "Number_of_Reads", 
                               "Number_of_Spots_Under_Tissue", 
                               "Fraction_Reads_in_Spots_Under_Tissue",
                               "Mean_Reads_per_Spot",
                               "Reads_Mapped_Confidently_to_Exonic_Regions",
                               "Median_Genes_per_Spot")]
table_met <- top_metrics_gex %>%
  gt::gt() %>%
  gt::fmt_percent(columns = c("Fraction_Reads_in_Spots_Under_Tissue", "Reads_Mapped_Confidently_to_Exonic_Regions"), 
                  scale_values = FALSE, decimals = 1) %>%
  gt::fmt_number(columns = "Number_of_Reads", scale_by = 1 / 1E6, pattern = "{x}M") %>% 
  gt::tab_header(
    title = gt::md("**GEX QC metrics**"),
    subtitle = ("spaceranger")
  ) %>%
  gt::cols_label(
    Sample_ID = gt::md("**GEM ID**"),
    Number_of_Reads = gt::md("**Number of Reads**"),
    Number_of_Spots_Under_Tissue = gt::md("**Number of_Spots Under Tissue**"),
    Fraction_Reads_in_Spots_Under_Tissue = gt::md("**Fraction Reads in Spots Under Tissue**"),
    Mean_Reads_per_Spot = gt::md("**Mean Reads per Spot**"),
    Reads_Mapped_Confidently_to_Exonic_Regions = gt::md("**Fraction of Reads Mapped to Exonic Reads**"),
    Median_Genes_per_Spot = gt::md("**Median Genes per Spot**")
  ) %>%
  gt::tab_style(style = list(
                  cell_fill(color = "red")),
                locations =  cells_body(vars(Number_of_Spots_Under_Tissue),
                                        rows = Number_of_Spots_Under_Tissue <20)
                
  ) %>%
  gt::tab_style(style = list(
    cell_fill(color = "red")),
    locations =  cells_body(vars(Fraction_Reads_in_Spots_Under_Tissue),
                            rows = Fraction_Reads_in_Spots_Under_Tissue <50)
  )%>%
  gt::tab_style(style = list(
    cell_fill(color = "red")),
    locations =  cells_body(vars(Reads_Mapped_Confidently_to_Exonic_Regions),
                            rows = Reads_Mapped_Confidently_to_Exonic_Regions <30)
  )%>%
  gt::tab_style(style = list(
    cell_fill(color = "red")),
    locations =  cells_body(vars(Reads_Mapped_Confidently_to_Exonic_Regions),
                            rows = Reads_Mapped_Confidently_to_Exonic_Regions <30)
  )%>%
  gt::tab_style(style = list(
    cell_fill(color = "red")),
    locations =  cells_body(vars(Mean_Reads_per_Spot),
                            rows = Mean_Reads_per_Spot <50000)
  )


gtsave(table_met, filename = "./metrics_summmary_table.png")

# ============================================================================= #
#################################### 3. Plots ###################################
# ============================================================================= #


############### Q30
qc_seq_vars <- c("Q30_Bases_in_Barcode",
                 "Q30_Bases_in_RNA_Read",
                 "Q30_Bases_in_UMI")

plot1 <- ggplot2::ggplot(metrics, ggplot2::aes_string(x = "Sample_ID",
                                                      y = "Q30_Bases_in_RNA_Read",
                                                      fill = "Sample_ID")) +
  ggplot2::geom_col() +
  geom_abline(slope=0, intercept=0.65,  col = "red",lty=2) + 
  ggplot2::theme_bw() +
  ggplot2::scale_fill_brewer(palette = "Set2") +
  ggplot2::ylim(0, 1) +
  ggplot2::labs(x = "Libraries (SAMPLE IDs)",
                y = str_c(str_replace_all("Q30_Bases_in_RNA_Read", "_", " "), " (%)")) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 8),
    axis.text.x = ggplot2::element_text(hjust = 1, angle = 45),
    strip.placement = "outside",
    strip.background = ggplot2::element_rect(colour = NA),
    legend.position = "none")

plot2 <- ggplot2::ggplot(metrics, ggplot2::aes_string(x = "Sample_ID",
                                                      y = "Q30_Bases_in_Barcode",
                                                      fill = "Sample_ID")) +
  ggplot2::geom_col() +
  geom_abline(slope=0, intercept=0.8,  col = "red",lty=2) + 
  ggplot2::theme_bw() +
  ggplot2::scale_fill_brewer(palette = "Set2") +
  ggplot2::ylim(0, 1) +
  ggplot2::labs(x = "Libraries (SAMPLE IDs)",
                y = str_c(str_replace_all("Q30_Bases_in_Barcode", "_", " "), " (%)")) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 8),
    axis.text.x = ggplot2::element_text(hjust = 1, angle = 45),
    strip.placement = "outside",
    strip.background = ggplot2::element_rect(colour = NA),
    legend.position = "none")

plot3 <- ggplot2::ggplot(metrics, ggplot2::aes_string(x = "Sample_ID",
                                                      y = "Q30_Bases_in_UMI",
                                                      fill = "Sample_ID")) +
  ggplot2::geom_col() +
  geom_abline(slope=0, intercept=0.65,  col = "red",lty=2) + 
  ggplot2::theme_bw() +
  ggplot2::scale_fill_brewer(palette = "Set2") +
  ggplot2::ylim(0, 1) +
  ggplot2::labs(x = "Libraries (SAMPLE IDs)",
                y = str_c(str_replace_all("Q30_Bases_in_UMI", "_", " "), " (%)")) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 8),
    axis.text.x = ggplot2::element_text(hjust = 1, angle = 45),
    strip.placement = "outside",
    strip.background = ggplot2::element_rect(colour = NA),
    legend.position = "none")


png("./plot_Q30.png", width = 700, height = 500)
cowplot::plot_grid(plotlist = list(plot1,plot2,plot3),
                   nrow = 1,
                   ncol = 3,
                   align = "hv",
                   axis = "trbl")
dev.off()

############ Mapping

plot1 <- ggplot2::ggplot(metrics, ggplot2::aes_string(x = "Sample_ID",
                                                      y = "Reads_Mapped_Confidently_to_Genome",
                                                      fill = "Sample_ID")) +
  ggplot2::geom_col() +
  #geom_abline(slope=0, intercept=0.65,  col = "red",lty=2) + 
  ggplot2::theme_bw() +
  ggplot2::scale_fill_brewer(palette = "Set2") +
  ggplot2::ylim(0, 100) +
  ggplot2::labs(x = "Libraries (Sample IDs)",
                y = str_c(str_replace_all("Reads_Mapped_Confidently_to_Genome", "_", " "), " (%)")) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 10),
    axis.text.x = ggplot2::element_text(hjust = 1, angle = 45),
    axis.text.y = ggplot2::element_text(size = 9),
    strip.placement = "outside",
    strip.background = ggplot2::element_rect(colour = NA),
    legend.position = "none")

plot2 <- ggplot2::ggplot(metrics, ggplot2::aes_string(x = "Sample_ID",
                                                      y = "Reads_Mapped_Confidently_to_Intergenic_Regions",
                                                      fill = "Sample_ID")) +
  ggplot2::geom_col() +
  #geom_abline(slope=0, intercept=0.65,  col = "red",lty=2) + 
  ggplot2::theme_bw() +
  ggplot2::scale_fill_brewer(palette = "Set2") +
  ggplot2::ylim(0, 100) +
  ggplot2::labs(x = "Libraries (Sample IDs)",
                y = str_c(str_replace_all("Reads_Mapped_Confidently_to_Intergenic_Regions", "_", " "), " (%)")) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 10),
    axis.text.x = ggplot2::element_text(hjust = 1, angle = 45),
    axis.text.y = ggplot2::element_text(size = 9),
    strip.placement = "outside",
    strip.background = ggplot2::element_rect(colour = NA),
    legend.position = "none")

plot3 <- ggplot2::ggplot(metrics, ggplot2::aes_string(x = "Sample_ID",
                                                      y = "Reads_Mapped_Confidently_to_Intronic_Regions",
                                                      fill = "Sample_ID")) +
  ggplot2::geom_col() +
  #geom_abline(slope=0, intercept=0.65,  col = "red",lty=2) + 
  ggplot2::theme_bw() +
  ggplot2::scale_fill_brewer(palette = "Set2") +
  ggplot2::ylim(0, 100) +
  ggplot2::labs(x = "Libraries (Sample IDs)",
                y = str_c(str_replace_all("Reads_Mapped_Confidently_to_Intronic_Regions", "_", " "), " (%)")) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 10),
    axis.text.x = ggplot2::element_text(hjust = 1, angle = 45),
    axis.text.y = ggplot2::element_text(size = 9),
    strip.placement = "outside",
    strip.background = ggplot2::element_rect(colour = NA),
    legend.position = "none")

plot4 <- ggplot2::ggplot(metrics, ggplot2::aes_string(x = "Sample_ID",
                                                      y = "Reads_Mapped_Confidently_to_Exonic_Regions",
                                                      fill = "Sample_ID")) +
  ggplot2::geom_col() +
  #geom_abline(slope=0, intercept=0.65,  col = "red",lty=2) + 
  ggplot2::theme_bw() +
  ggplot2::scale_fill_brewer(palette = "Set2") +
  ggplot2::ylim(0, 100) +
  ggplot2::labs(x = "Libraries (Sample IDs)",
                y = str_c(str_replace_all("Reads_Mapped_Confidently_to_Exonic_Regions", "_", " "), " (%)")) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 10),
    axis.text.x = ggplot2::element_text(hjust = 1, angle = 45),
    axis.text.y = ggplot2::element_text(size = 9),
    strip.placement = "outside",
    strip.background = ggplot2::element_rect(colour = NA),
    legend.position = "none")

png("./plot_mapping.png", width = 700, height = 700)

cowplot::plot_grid(plotlist = list(plot1,plot2,plot3,plot4),
                   nrow = 2,
                   ncol = 2,
                   align = "hv",
                   axis = "trbl")
dev.off()

############ Sequence saturation and depth

gg_lib_size <- metrics %>%
  dplyr::mutate(Number_of_Reads_mil = Number_of_Reads / 1000000) %>%
  ggplot2::ggplot(ggplot2::aes(x = Sample_ID, y = Number_of_Reads_mil, fill = Sample_ID)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::theme_bw() +
  ggplot2::scale_fill_brewer(palette = "Set2") +
  ggplot2::labs(x = "Libraries (SAMPLE IDs)", y = "Library size (in millions)") +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 8),
    axis.text.x = ggplot2::element_text(hjust = 1, angle = 45),
    strip.placement = "outside",
    strip.background = ggplot2::element_rect(colour = NA),
    legend.position = "none")

gg_qc_seq_sat <- metrics %>%
  dplyr::mutate(Sequencing_Saturation_perc = Sequencing_Saturation,
                Mean_Reads_per_Spot_tho = Mean_Reads_per_Spot / 1000) %>%
  ggplot2::ggplot(ggplot2::aes(x = Mean_Reads_per_Spot_tho,
                               y = Sequencing_Saturation_perc, color = Sample_ID)) +
  ggplot2::geom_point() +
  ggplot2::theme_bw() +
  ggplot2::scale_color_brewer(palette = "Set2") +
  ggplot2::ylim(0, 1) +
  ggrepel::geom_text_repel(ggplot2::aes(label = Sample_ID), size = 4) +
  ggplot2::labs(x = "Mean Reads per Spot (in thousands)", y = "Sequencing Saturation") +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 10),
    legend.position = "none")

gg_qc_seq_depth_cell <- metrics %>%
  dplyr::mutate(Mean_Reads_per_Spot_tho = Mean_Reads_per_Spot / 1000) %>%
  ggplot2::ggplot(ggplot2::aes(x = Mean_Reads_per_Spot_tho,
                               y = Median_Genes_per_Spot, color = Sample_ID)) +
  ggplot2::geom_point() +
  ggplot2::theme_bw() +
  ggplot2::scale_color_brewer(palette = "Set2") +
  ggrepel::geom_text_repel(ggplot2::aes(label = Sample_ID), size = 4) +
  ggplot2::labs(x = "Mean Reads per spot (in thousands)",
                y = "Mean Detected Genes per spot") +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 10),
    legend.position = "none")

gg_qc_seq_depth <- metrics %>%
  dplyr::mutate(Number_of_Reads_mil = Number_of_Reads / 1000000) %>%
  ggplot2::ggplot(ggplot2::aes(x = Number_of_Reads_mil,
                               y = Total_Genes_Detected, color = Sample_ID)) +
  ggplot2::geom_point() +
  ggplot2::theme_bw() +
  ggplot2::scale_color_brewer(palette = "Set2") +
  ggrepel::geom_text_repel(ggplot2::aes(label = Sample_ID), size = 4) +
  ggplot2::labs(x = "Number of Reads (in millions)", y = "Total Genes Detected") +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 10),
    legend.position = "none")

png("./seq_depth.png", width = 700, height = 700)
cowplot::plot_grid(gg_lib_size, gg_qc_seq_sat, 
                   gg_qc_seq_depth_cell, gg_qc_seq_depth, 
                   nrow = 2,
                   ncol = 2,
                   align = "hv",
                   axis = "trbl")
dev.off()


