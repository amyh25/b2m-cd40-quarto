# aesthetics.R

# themes --------

theme_set(theme_classic() + theme(
  strip.background = element_blank()
))

# palettes --------

palette_treatment <- c(
  "CD40 agonist" = "red",
  "Isotype" = "grey"
)

palette_genotype <- c(
  "B2m KO" = "burlywood4",
  "Ctrl" = "black"
)

palette_genotype_renamed <- c(
  "B2m KO" = "burlywood4",
  "WT" = "black"
)

# cluster ann palettes ------

cluster_ann_list <- list()

cluster_ann_list$FULL <- read_excel(file.path(results_dir, "cluster_annotations.xlsx")) %>%
  tidyr::fill(Subset, clustering) %>%
  filter(Subset == "so_filterRNA200mito10.rds") %>%
  filter(clustering == "RNA_snn_res.0.5") %>%
  mutate(cluster_ID_name = paste0(cluster, ": ", cluster_ID))
# palette_cluster <- cluster_ann_list$FULL$color %>% set_names(cluster_ann_list$FULL$cluster_ID_name)
palette_cluster <- social_palette %>% set_names(cluster_ann_list$FULL$cluster_ID_name[1:length(social_palette)])
cluster_ID_name_ordered <- cluster_ann_list$FULL %>%
  arrange(cluster_ID_grouped) %>%
  .$cluster_ID_name

cluster_ann_list$Tcell <- readxl::read_excel(file.path(results_dir, "cluster_annotations.xlsx")) %>%
  tidyr::fill(Subset, clustering) %>%
  filter(Subset == "so_filterRNA200mito10_subsetTremovedoublets.rds") %>%
  filter(clustering == "RNA_snn_res.0.4") %>%
  mutate(cluster_ID_name = paste0(cluster, ": ", cluster_ID))
palette_cluster_T <- cluster_ann_list$Tcell$color %>% set_names(cluster_ann_list$Tcell$cluster_ID_name)
cluster_T_name_ordered <- cluster_ann_list$Tcell %>%
  arrange(cluster_ID_grouped) %>%
  .$cluster_ID_name

cluster_ann_list$Myeloid <- readxl::read_excel(file.path(results_dir, "cluster_annotations.xlsx")) %>%
  tidyr::fill(Subset, clustering) %>%
  filter(Subset == "so_filterRNA200mito10_subsetMremovedoublets.rds") %>%
  filter(clustering == "RNA_snn_res.0.5") %>%
  mutate(cluster_ID_name = paste0(cluster, ": ", cluster_ID))
palette_cluster_M <- cluster_ann_list$Myeloid$color %>% set_names(cluster_ann_list$Myeloid$cluster_ID_name)
cluster_M_name_ordered <- cluster_ann_list$Myeloid %>%
  arrange(cluster_ID_grouped) %>%
  .$cluster_ID_name
