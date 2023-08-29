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

checkdir("stem")

## prepare STEM inputs
metaphlanFile <- "tables/merged_abundance_table_format.txt"
metaphlan <- read.table(metaphlanFile, sep = "\t", header = TRUE, check.names = FALSE)

sampleInfoAll <- read.xlsx(glue::glue("{docsdir}/sampleinfo_kpc.xlsx"))
sampleInfo <- sampleInfoAll[match(colnames(metaphlan)[-c(1,2)], sampleInfoAll$sample_id), ]
stem <- sampleInfo[!is.na(sampleInfo$STEM), ]
stem <- stem[order(stem$STEM, decreasing = FALSE), ]
for (i in unique(stem$sample_group)){
  smps <- stem[stem$sample_group == i, ]$sample_id
  out <- metaphlan[, c("species", smps)]
  write.table(out, glue::glue("stem/S{i}.txt"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
}

## perform STEM java program
# results are stem/S{i}_stem.txt

# ## find overlaps between each sample's stem profile
# files <- dir(path = "stem",
#              full.names = FALSE,
#              pattern = "_stem.txt")
# data <- map(glue::glue("stem/{files}"), read.table, sep = "\t", header = TRUE)
# names(data) <- files
# op <- reduce(data, inner_join, by = c("species", "Profile"))
# ## no overlap


## merge all the species and its STEM profile info
files <- dir(path = "stem",
             full.name = FALSE,
             pattern = "_stem.txt")
data <- map(glue::glue("stem/{files}"), read.table, sep = "\t", header = TRUE) %>%
  setNames(., files) %>%
  reduce(full_join, by = "species")
df2excel(data, "stem/STEM_profiles_raw.xlsx")





