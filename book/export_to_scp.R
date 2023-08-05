library(here)
source(here("lib", "setup.R"))
source(here("lib", "helper_Seurat.R"))
source(here("lib", "helper.R"))

list.files(results_dir)
scp_dir <- file.path(results_dir, "scp")

# export counts and data matrices -----

so_raw <- readRDS(file.path(results_dir, "so_raw.rds"))
scp_export_counts(so_raw, scp_dir)
scp_export_data(so_raw, scp_dir)
scp_export_metadata_file(
    so_raw,
    output_dir = scp_dir,
    str_subset_regex_group = "hash.ID|seq|sample|mouse|genotype|treatment|group",
    str_subset_regex_numeric = "nCount_|nFeature|percent"
)

# export clustering -------

so <- readRDS(file.path(results_dir, "so_filterRNA200mito10.rds"))
scp_export_clustering_file(
    so,
    output_dir = scp_dir,
    filename = "scp_clustering_filterRNA200mito10",
    str_subset_regex_group = "clusters|general_cell_types",
    str_subset_regex_numeric = NULL
)

so_T <- readRDS(file.path(results_dir, "so_filterRNA200mito10_subsetT.rds"))
new_metadata <- so_T@meta.data %>%
    select(str_subset(colnames(.), "clusters|general_")) %>%
    mutate(across(.cols = everything(), ~str_replace(., ":", ""))) %>%
    mutate(across(.cols = everything(), ~str_replace(., "\\/", "\\."))) %>%
    mutate(across(.cols = everything(), ~str_replace(., "ï", "i")))
so_T <- AddMetaData(so_T, new_metadata)
scp_export_clustering_file(
    so_T,
    output_dir = scp_dir,
    filename = "scp_clustering_filterRNA200mito10_subsetT",
    str_subset_regex_group = "clusters|general_",
    str_subset_regex_numeric = NULL
)

so_T_rm <- readRDS(file.path(results_dir, "so_filterRNA200mito10_subsetTremovedoublets.rds"))
new_metadata <- so_T_rm@meta.data %>%
    select(str_subset(colnames(.), "clusters|general_")) %>%
    mutate(across(.cols = everything(), ~str_replace(., ":", ""))) %>%
    mutate(across(.cols = everything(), ~str_replace(., "\\/", "\\."))) %>%
    mutate(across(.cols = everything(), ~str_replace(., "ï", "i")))
so_T_rm <- AddMetaData(so_T_rm, new_metadata)
scp_export_clustering_file(
    so_T_rm,
    output_dir = scp_dir,
    filename = "scp_clustering_filterRNA200mito10_subsetTremovedoublets",
    str_subset_regex_group = "clusters|general_",
    str_subset_regex_numeric = NULL
)

so_gdT <- readRDS(file.path(results_dir, "so_filterRNA200mito10_subsetTremovedoublets_subsetgdT.rds"))
new_metadata <- so_gdT@meta.data %>%
    select(str_subset(colnames(.), "clusters|general_")) %>%
    mutate(across(.cols = everything(), ~str_replace(., ":", ""))) %>%
    mutate(across(.cols = everything(), ~str_replace(., "\\/", "\\."))) %>%
    mutate(across(.cols = everything(), ~str_replace(., "ï", "i")))
so_gdT <- AddMetaData(so_gdT, new_metadata)
scp_export_clustering_file(
    so_gdT,
    output_dir = scp_dir,
    filename = "scp_clustering_filterRNA200mito10_subsetTremovedoublets_subsetgdT",
    str_subset_regex_group = "clusters|general_",
    str_subset_regex_numeric = NULL
)

so_M <- readRDS(file.path(results_dir, "so_filterRNA200mito10_subsetMremovedoublets.rds"))
new_metadata <- so_M@meta.data %>%
    select(str_subset(colnames(.), "clusters|general_")) %>%
    mutate(across(.cols = everything(), ~str_replace(., ":", ""))) %>%
    mutate(across(.cols = everything(), ~str_replace(., "\\/", "\\."))) %>%
    mutate(across(.cols = everything(), ~str_replace(., "ï", "i")))
so_M <- AddMetaData(so_M, new_metadata)
scp_export_clustering_file(
    so_M,
    output_dir = scp_dir,
    filename = "scp_clustering_filterRNA200mito10_subsetMremovedoublets",
    str_subset_regex_group = "clusters|general_",
    str_subset_regex_numeric = NULL
)


so_NKdropCD3 <- readRDS(file.path(results_dir, "so_filterRNA200mito10_subsetTremovedoublets_subsetNKdropCD3.rds"))
new_metadata <- so_NKdropCD3@meta.data %>%
    select(str_subset(colnames(.), "clusters|general_")) %>%
    scp_clean_colnames()
so_NKdropCD3 <- AddMetaData(so_NKdropCD3, new_metadata)
scp_export_clustering_file(
    so_NKdropCD3,
    output_dir = scp_dir,
    filename = "scp_clustering_filterRNA200mito10_subsetTremovedoublets_subsetNKdropCD3",
    str_subset_regex_group = "RNA_snn_res",
    str_subset_regex_numeric = NULL
)
