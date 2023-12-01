# helper.R

# process helper functions ------

#' @name run_or_readRDS
#'
#' A kind of hack cache-ing function that I use in my notebook for
#' convenience. Will read in an RDS object at the specified path if it exists,
#' otherwise will run the specified function.
#'
#' @param path a path to look for an RDS object with the saved results
#' @param f a function that will output results that you want saved
#' @param verbose boolean specifying verbosity
#' @return a tibble

run_or_readRDS <- function(path, f, verbose = TRUE, overwrite = FALSE) {
  if (!file.exists(path)) {
    if (verbose) message(paste0("File at ", path, " does not exist. Running function..."))
    out <- f()
    saveRDS(out, path)
    if (verbose) message(paste0("Saving at ", path, "."))
  } else if (overwrite) {
    if (verbose) message(paste0("File found. Overwriting."))
    out <- f()
    saveRDS(out, path)
  } else  {
    if (verbose) message("Found file at ", path, ". Reading...")
    out <- readRDS(path)
    if (verbose) message("Done!")
  }
  return(out)
}

# general helper functions -------

#' @name calculate_proportion
#' @param .tbl a tibble
#' @param var the variable to calculate propotion of
#' @param group_var a grouping variable to calculate proportion within
#' @return a tibble

calculate_proportion <- function(.tbl, var, group_var) {
  var <- enquo(var)
  group_var <- enquo(group_var)
  out <- .tbl %>%
    group_by(!!var,!!group_var) %>%
    summarise(count = n()) %>%
    group_by(!!group_var) %>%
    mutate(total = sum(count), prop = count / total)
  out <- out %>%
    pivot_wider(
      id_cols = !!group_var,
      names_from = !!var,
      values_from = prop,
      values_fill = 0
    ) %>%
    pivot_longer(
      cols = 2:ncol(.),
      names_to = rlang::as_name(var),
      values_to = "prop"
    ) %>%
    left_join(select(out,-prop),
              by = c(rlang::as_name(group_var),
                     rlang::as_name(var)))
  return(out)
}

#' @name mutate_signif_levels
#'
#' Creates a new column named `signif` which contains a character column
#' that denotes the significant level
#'
#' @param .data the data frame to modify
#' @param p_value_sym the symbol denoting the pvalue column
#'
#' @return a modified data frame

mutate_signif_levels <- function(.data, p_value_sym) {
  p_value_sym <- enquo(p_value_sym)
  .data %>%
    mutate(signif = case_when(!!p_value_sym < 0.001 ~ "***",
                              !!p_value_sym < 0.01 ~ "**",
                              !!p_value_sym < 0.05 ~ "*",
                              TRUE ~ ""))
}

# plotting helper functions ------


#' @name plot_proportion
#' @param .tbl a tibble
#' @param x symbol designating column of tibble with group information
#' @param y symbol designating column of tibble with proportion values
#' @param fill
#' @return ggplot object

plot_proportion <- function(.tbl, x = group, y = prop, fill = clusters) {
  x <- enquo(x)
  y <- enquo(y)
  fill <- enquo(fill)
  ggplot(.tbl) +
    aes(!!x, !!y, fill = !!fill) +
    geom_col() +
    remove_x_spine() +
    scale_y_continuous(expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "none")
}

#' @name plot_two_groups
#'
#' @param .tbl a data frame or tibble
#' @param var variable to plot as symbol
#' @param group group to plot as symbol
#'
#' Use with fig.asp = 2, fig.width = 2.
#'
#' @return ggplot object

plot_two_groups <- function(.tbl,
                            group,
                            var = prop * 100,
                            stat.method = "wilcox.test",
                            color_palette = c("black", "lightblue")) {
  group <- enquo(group)
  var <- enquo(var)

  mean_df <- .tbl %>%
    group_by(!!group) %>%
    summarise(mean = mean(!!var),
              se = sd(!!var) / sqrt(n()))
  groups_vec <- unique(as.character(.tbl[[rlang::as_name(group)]]))

  .tbl %>%
    ggplot() +
    aes(!!group, y = !!var) +
    geom_errorbar(
      aes(x = !!group,
          y = mean,
          ymin = ..y..,
          ymax = ..y..),
      data = mean_df,
      inherit.aes = FALSE,
      width = 0.5,
      size = 1
    ) +
    geom_errorbar(
      aes(x = !!group,
          ymin = mean - se,
          ymax = mean + se),
      data = mean_df,
      inherit.aes = FALSE,
      width = 0.25,
      size = 0.5
    ) +
    geom_jitter(
      aes(shape = !!group, fill = !!group),
      color = "black",
      size = 3,
      width = 0.3
    ) +
    ggpubr::stat_compare_means(comparison = list(groups_vec),
                               label = "p.signif",
                               method = stat.method,
                               tip.length = 0) +
    scale_y_continuous(expand = c(0, 0), limits = ~ c(0, max(..1) * 1.1)) +
    scale_shape_manual(values = c(21, 22)) +
    scale_fill_manual(values = color_palette) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())
}

#' @name plot_volcano
#' Generates a pretty volcano plot, by default takes in dataframes
#' formatted like the output of Seurat::FindMarkers
#' Pretty by my standards
#' @param volcano_df data frame with volcano data
#' @param var_fc fold change variable
#' @param var_p pvalue variable
#' @param thresh pval-threshold of where to display gene names up to. Can be
#' a single number, which will set a single cutoff, or a numerical vector of
#' length 2, which will set two cutoffs on the left and the right side of the
#' volcano
#' @return a ggplot2 object displaying a pretty volcano plot

plot_volcano <- function(volcano_df,
                         var_fc = avg_log2FC,
                         var_pval = p_val_adj,
                         var_gene = gene,
                         thresh = Inf) {
  var_fc <- enquo(var_fc)
  var_pval <- enquo(var_pval)
  var_gene <- enquo(var_gene)


  p <- volcano_df %>%
    ggplot() +
    aes(!!var_fc,-log10(!!var_pval)) +
    geom_scattermore(pointsize = 2.1, pixels = c(1024, 1024)) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed")

  # label points based on threshold
  if (!is.numeric(thresh))
    stop("Error: thresh must be numeric")
  if (length(thresh) == 1) {
    if (!(thresh == Inf)) {
      p <- p +
        geom_text_repel(
          aes(label = !!var_gene),
          data = ~ ..1 %>% filter(p_val_adj < thresh),
          max.overlaps = Inf,
          size = 3
        ) +
        annotate(
          "segment",
          x = -Inf,
          xend = Inf,
          y = -log10(thresh),
          yend = -log10(thresh),
          linetype = "dashed"
        )
    }
  } else {
    if (length(thresh) != 2)
      stop("Error: thresh must be have either length 1 or 2")
    p <- p +
      geom_text_repel(
        aes(label = !!var_gene),
        data = ~ ..1 %>%
          filter((p_val_adj < thresh[2] & avg_log2FC > 0) |
                   (p_val_adj < thresh[1] & avg_log2FC < 0)
          ),
        max.overlaps = Inf,
        size = 3
      ) +
      annotate(
        "segment",
        x = 0,
        xend = Inf,
        y = -log10(thresh[2]),
        yend = -log10(thresh[2]),
        linetype = "dashed"
      )  +
      annotate(
        "segment",
        x = -Inf,
        xend = 0,
        y = -log10(thresh[1]),
        yend = -log10(thresh[1]),
        linetype = "dashed"
      )
  }

  p <- p +
    theme(aspect.ratio = 1,
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    labs(caption = "genes above dashed horizontal lines are labelled")

  return(p)
}


plot_gsea_curve <- function(leg_df,
                            group_var,
                            pt_size = 2,
                            pt_alpha = 0.8,
                            line_size = 1) {
  group_var <- enquo(group_var)
  curve_plot <- leg_df %>%
    filter(in_gs) %>%
    ggplot() +
    aes(rank, running_es) +
    geom_point(aes(color = !!group_var),
               size = pt_size, alpha = pt_alpha) +
    geom_hline(yintercept = 0, size = line_size) +
    xlim(min(leg_df$rank), max(leg_df$rank)) +
    ylab("Running ES") +
    scale_x_continuous(expand = c(0, 0)) +
    theme(
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "top"
    )

  return(curve_plot)

}

