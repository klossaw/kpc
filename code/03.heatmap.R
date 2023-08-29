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


## relative abundance heatmap
## Family
tbl <- mpse1 %>% 
  mp_cal_abundance(.abundance = Abundance, force = T, relative = F) %>% 
  mp_extract_abundance(taxa.class = "Family", topn = 20) %>% 
  tidyr::unnest(AbundanceBySample)
idx <- factor(sampleInfo$week, levels = sort(unique(sampleInfo$week)))
samplesOrder <- sampleInfo[order(idx), ]$sample_id
df <- tbl
df$Sample <- factor(df$Sample, levels = samplesOrder)
p_heat <- ggplot(df, aes(x = Sample, y = label)) + geom_tile(aes(fill = Abundance)) +
  scale_fill_distiller("Abundance", palette = "RdYlBu") + 
  facet_grid(~sample_group, scales = "free") +
  ylab("Family") +
  theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y = element_text(colour = 'black', size = 10, face = 'italic'),
        axis.ticks.x=element_blank(),
        legend.key = element_blank(), 
        strip.background = element_rect(colour="black", fill="white"))
ggsave("figures/heatmap_Family_samples.pdf", height = 4, width = 21)


## Species
tbl <- mpse1 %>% 
  mp_cal_abundance(.abundance = Abundance, force = T, relative = F) %>% 
  mp_extract_abundance(taxa.class = "OTU", topn = 20) %>% 
  tidyr::unnest(AbundanceBySample)
idx <- factor(sampleInfo$week, levels = sort(unique(sampleInfo$week)))
samplesOrder <- sampleInfo[order(idx), ]$sample_id
df <- tbl
df$Sample <- factor(df$Sample, levels = samplesOrder)
p_heat <- ggplot(df, aes(x = Sample, y = label)) + geom_tile(aes(fill = Abundance)) +
  scale_fill_distiller("Abundance", palette = "RdYlBu") + 
  facet_grid(~sample_group, scales = "free") +
  ylab("Family") +
  theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y = element_text(colour = 'black', size = 10, face = 'italic'),
        axis.ticks.x=element_blank(),
        legend.key = element_blank(), 
        strip.background = element_rect(colour="black", fill="white"))
ggsave("figures/heatmap_Species_samples.pdf", height = 4, width = 21)

df$stage <- factor(df$stage, levels = names(cols))
p_heat <- ggplot(df, aes(x = Sample, y = label)) + geom_tile(aes(fill = Abundance)) +
  scale_fill_distiller("Abundance", palette = "RdYlBu") + 
  facet_grid(~stage, scales = "free") +
  ylab("Family") +
  theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y = element_text(colour = 'black', size = 10, face = 'italic'),
        axis.ticks.x=element_blank(),
        legend.key = element_blank(), 
        strip.background = element_rect(colour="black", fill="white"))
ggsave("figures/heatmap_Species_stage.pdf", height = 4, width = 21)
