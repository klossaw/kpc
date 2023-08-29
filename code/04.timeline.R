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

## create outdir
outdir_week <- "figures/timeline_week"
checkdir(outdir_week)
outdir_week_scale <- "figures/timeline_week_scale"
checkdir(outdir_week_scale)

## format input file
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

## generate MPSE object
mpse1 <- mp_import_metaphlan(profile = formatFile, mapfilename = sampleInfoFile)

## get abundance table
tbl <- mpse1 %>% 
  mp_cal_abundance(.abundance = Abundance, force = T, relative = F) %>% 
  mp_extract_abundance(taxa.class = "OTU", topn = NULL) %>% 
  tidyr::unnest(AbundanceBySample)
tbl <- tbl[which(tbl$stage != "wildtype"), ]


## timeline week plot for each species
week <- sampleInfo[match(tbl$Sample, sampleInfo$sample_id), ]$week
df <- cbind(tbl, week)
df$week <- factor(df$week, levels = sort(unique(week)))
for (s in unique(df$label)){
  p <- ggplot(df[df$label == s, ], aes(x = week, y = Abundance, color = sample_group, group = sample_group)) +
    geom_point(size = 1.5) +
    geom_line(size = 1) +
    labs(title = s, x = "Week", y = "Relative Abundance (%)") +
    theme_bw() +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5),
          plot.margin = margin(t = 1, r = 0.5, b = 0.5, l = 0.5, unit = "cm"))
  ggsave(glue::glue("{outdir_week}/{s}.timeline.pdf"), p, height = 8, width = 21, units = "cm", device = "pdf")
}


## timeline: week scale
week_scale <- sampleInfo[match(tbl$Sample, sampleInfo$sample_id), ]$week_scale
df <- cbind(tbl, week_scale)
df <- df[!is.na(df$week_scale), ]
df$week_scale <- factor(df$week_scale, levels = sort(unique(week_scale)))
for (s in unique(df$label)){
  p <- ggplot(df[df$label == s, ], aes(x = week_scale, y = Abundance, color = sample_group, group = sample_group)) +
    geom_point(size = 1.5) +
    geom_line(size = 1) +
    geom_vline(xintercept = 8, linetype = 'dashed', size = 0.3, alpha = 0.7) +
    labs(title = s, x = "Week scale", y = "Relative Abundance (%)") +
    theme_bw() +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5),
          plot.margin = margin(t = 1, r = 0.5, b = 0.5, l = 0.5, unit = "cm"))
  ggsave(glue::glue("{outdir_week_scale}/{s}.timeline.pdf"), p, height = 8, width = 21, units = "cm", device = "pdf")
}
