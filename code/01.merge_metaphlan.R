pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "vroom", "jhtools", "glue", "openxlsx", "ggsci", "patchwork", "cowplot", "tidyverse", "dplyr")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "kpc"
dataset <- "qzhang"
workdir <- glue::glue("~/projects/{project}/analysis/{dataset}/microbiota/metaphlan3")
setwd(workdir)

checkdir("tables")
cmd <- "/cluster/home/jhuang/.conda/envs/metaphlan/bin/merge_metaphlan_tables.py output/profiled*.txt > tables/merged_abundance_table.txt"
system(cmd, wait = TRUE)
