pkgs <- c("MicrobiotaProcess", "ggplot2", "ggalluvial", "openxlsx", "gghalves", "ggsignif",
          "fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "vroom", "jhtools", "glue", "openxlsx", "ggsci", "patchwork", "cowplot", "tidyverse", "dplyr")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "kpc"
dataset <- "qzhang"
workdir <- glue::glue("~/projects/{project}/analysis/{dataset}/microbiota/metaphlan3")
docsdir <- glue::glue("~/projects/{project}/docs")
setwd(workdir)

checkdir("figures")

metaphlanFile <- "tables/merged_abundance_table.txt"
metaphlan <- read.table(metaphlanFile, sep = "\t", header = TRUE, check.names = FALSE)
sampleNames <- colnames(metaphlan)[-c(1,2)]
sampleNames_new <- gsub("profiled_", "S", sampleNames)
colnames(metaphlan) <- c("species", "taxid", sampleNames_new)
formatFile <- "tables/merged_abundance_table_format.txt"
write.table(metaphlan, formatFile, sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)

sampleInfoAll <- read.xlsx(glue::glue("{docsdir}/sampleinfo_kpc.xlsx"))
sampleInfo <- sampleInfoAll[match(colnames(metaphlan)[-c(1,2)], sampleInfoAll$sample_id), ]
sampleInfoFile <- glue::glue("{docsdir}/sampleinfo_kpc_format.txt")
write.table(sampleInfo, sampleInfoFile, sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)

mpse1 <- mp_import_metaphlan(profile = formatFile, mapfilename = sampleInfoFile)


cols <- c("#9c9c9c", "#56b8b7", "#6a76b6", "#e39d97", "#e34e68")
names(cols) <- c("wildtype", "Normal", "pre_malignant", "ES_Tumor", "LS_Tumor")


## relative abundance barplot
## Family
tbl <- mpse1 %>% 
  mp_cal_abundance(.abundance = Abundance, force = T, relative = F) %>% 
  mp_extract_abundance(taxa.class = "Family", topn = 10) %>% 
  tidyr::unnest(AbundanceBySample)
idx <- factor(sampleInfo$week, levels = sort(unique(sampleInfo$week)))
samplesOrder <- sampleInfo[order(idx), ]$sample_id
df <- tbl
df$Sample <- factor(df$Sample, levels = samplesOrder)
p_abun <- ggplot(df, aes(x = Sample, y = Abundance, fill = label)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(x = NULL, y = "Relative Abundance (%)") +
  facet_grid(~sample_group, scales = "free") +
  scale_fill_brewer(name = "Class", type = "qual", palette = "Paired") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, size = 6))
ggsave("figures/relative_abundance_Family_samples.pdf", p_abun, height = 4, width = 30, device = "pdf")

## Genus
tbl <- mpse1 %>% 
  mp_cal_abundance(.abundance = Abundance, force = T, relative = F) %>% 
  mp_extract_abundance(taxa.class = "Genus", topn = 10) %>% 
  tidyr::unnest(AbundanceBySample)
idx <- factor(sampleInfo$week, levels = sort(unique(sampleInfo$week)))
samplesOrder <- sampleInfo[order(idx), ]$sample_id
df <- tbl
df$Sample <- factor(df$Sample, levels = samplesOrder)
p_abun <- ggplot(df, aes(x = Sample, y = Abundance, fill = label)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(x = NULL, y = "Relative Abundance (%)") +
  facet_grid(~sample_group, scales = "free") +
  scale_fill_brewer(name = "Class", type = "qual", palette = "Paired") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, size = 6))
ggsave("figures/relative_abundance_Genus_samples.pdf", p_abun, height = 4, width = 30, device = "pdf")

## species
tbl <- mpse1 %>% 
  mp_cal_abundance(.abundance = Abundance, force = T, relative = F) %>% 
  mp_extract_abundance(taxa.class = "OTU", topn = 10) %>% 
  tidyr::unnest(AbundanceBySample)
idx <- factor(sampleInfo$week, levels = sort(unique(sampleInfo$week)))
samplesOrder <- sampleInfo[order(idx), ]$sample_id
df <- tbl
df$Sample <- factor(df$Sample, levels = samplesOrder)
p_abun <- ggplot(df, aes(x = Sample, y = Abundance, fill = label)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(x = NULL, y = "Relative Abundance (%)") +
  facet_grid(~sample_group, scales = "free") +
  scale_fill_brewer(name = "Class", type = "qual", palette = "Paired") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, size = 6))
ggsave("figures/relative_abundance_Species_samples.pdf", p_abun, height = 4, width = 30, device = "pdf")
df$stage <- factor(df$stage, levels = names(cols))
p_abun <- ggplot(df, aes(x = Sample, y = Abundance, fill = label)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(x = NULL, y = "Relative Abundance (%)") +
  facet_grid(~stage, scales = "free") +
  scale_fill_brewer(name = "Class", type = "qual", palette = "Paired") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, size = 6))
ggsave("figures/relative_abundance_Species_stage.pdf", p_abun, height = 4, width = 30, device = "pdf")



## alpha diversity
mpse <- mpse1 %>% mp_cal_alpha(.abundance = Abundance, force = TRUE)
tbl <- mpse %>% mp_extract_sample()
tbl %<>% tidyr::pivot_longer(
  cols = !c("Sample", "sequence_id", "stage", "species", "KPC", "week", "sample_group"), 
  names_to = "measure", 
  values_to = "alpha") %>%
  dplyr::mutate(measure = forcats::fct_relevel(measure, 
                                               c("Observe", "Simpson", "Shannon", "J")))
tbl$stage <- factor(tbl$stage, levels = names(cols))
p <- ggplot(data = tbl, aes(x = stage, y = alpha, fill = stage)) +
  geom_half_violin(color = NA, side = "l", trim = FALSE) +
  geom_boxplot(aes(color = stage), fill = NA, position = position_nudge(x = .22), width = 0.2) +
  geom_half_point(side = "r", shape = 21, size = 1.2) +
  # geom_signif(comparisons = combn(names(cols), 2, simplify = FALSE), 
  #             margin_top = 0.2,
  #             test = "wilcox.test", 
  #             textsize = 2) +
  facet_wrap(facet = vars(measure), scales = "free_y", nrow = 1) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_color_manual(values = cols, guide = "none") +
  labs(x = NULL, y = "alpha diversity index") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        plot.margin = margin(t = 1, r = 2, b = 1, l = 1, unit = "cm"))
ggsave("figures/alpha_diversity.pdf", p, height = 8, width = 21, units = "cm", device = "pdf")


## PCA
tbl <- mpse1 %>%
  mp_decostand(.abundance = Abundance) %>%
  mp_cal_pca(.abundance = hellinger, action = "only")
x <- names(tbl)[grepl("PC1 ", names(tbl))] %>% as.symbol()
y <- names(tbl)[grepl("PC2 ", names(tbl))] %>% as.symbol()
p_pca <- ggplot(tbl) +
  geom_point(aes(x = !!x, y = !!y, color = stage), size = 4) +
  scale_color_manual(values = cols) +
  geom_vline(xintercept = 0, color = "grey20", linetype = 2) +
  geom_hline(yintercept = 0, color = "grey20", linetype = 2) +
  theme_bw()
ggsave("figures/pca.pdf", p_pca, height = 16, width = 21, units = "cm", device = "pdf")

## PCOA
tbl <- mpse1 %>%
  mp_decostand(.abundance = Abundance) %>%
  mp_cal_pcoa(.abundance = hellinger, distmethod = "bray", .dim = 2, action = "only")
x <- names(tbl)[grepl("PCo1 ", names(tbl))] %>% as.symbol()
y <- names(tbl)[grepl("PCo2 ", names(tbl))] %>% as.symbol()
p_pcoa <- ggplot(tbl) +
  geom_point(aes(x = !!x, y = !!y, color = stage), size = 4) +
  scale_color_manual(values = cols) +
  # stat_ellipse(aes(fill = stage), geom = "polygon", alpha = 0.5) +
  geom_vline(xintercept = 0, color = "grey20", linetype = 2) +
  geom_hline(yintercept = 0, color = "grey20", linetype = 2) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave("figures/pcoa.pdf", p_pcoa, height = 16, width = 21, units = "cm", device = "pdf")
