# =============================================================================
# Figure 6A-B | Bacteria-metabolite correlation analysis
# Run from the repository root. Input CSV files are stored in ./data/.
# Output figures and session information are saved in ./results/.
# =============================================================================

# 1. Packages -----------------------------------------------------------------
packages <- c("ggplot2", "ggalluvial", "dplyr", "pheatmap")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install required packages first: ",
       paste(missing, collapse = ", "),
       "\ninstall.packages(c(",
       paste(sprintf('"%s"', missing), collapse = ", "), "))")
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggalluvial)
  library(dplyr)
  library(pheatmap)
  library(grid)
})

# 2. Settings -----------------------------------------------------------------
DATA_DIR <- "data"
OUT_DIR <- "results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

FILES <- c(
  bacteria = file.path(DATA_DIR, "属水平肠道菌群.csv"),
  neg = file.path(DATA_DIR, "HNS代谢物-neg.csv"),
  pos = file.path(DATA_DIR, "HNS代谢物-pos.csv")
)

TOP_N <- 30
P_HEATMAP <- 0.05
P_SANKEY <- 0.05
RHO_SANKEY <- 0.40

NEG_COL <- "#2166AC"
MID_COL <- "#F7F7F7"
POS_COL <- "#B2182B"


# 3. Helper functions ---------------------------------------------------------
read_input <- function(file) {
  if (!file.exists(file)) stop("File not found: ", file)
  read.csv(file, header = TRUE, row.names = 1,
           check.names = FALSE, fileEncoding = "GBK")
}

format_metabolite <- function(x) {
  x <- enc2utf8(as.character(x))
  ursane <- grepl("Hydroxyursa-12,20", x, fixed = TRUE)
  flavonoid <- grepl("5,7-dihydroxy-3", x, ignore.case = TRUE) &
    grepl("chromen-4-one", x, ignore.case = TRUE)

  x[ursane] <- "(3β,5ξ,9ξ)-3-Hydroxyursa-\n12,20(30)-dien-28-oic acid"
  x[flavonoid] <- paste0(
    "5,7-dihydroxy-3-(4-hydroxyphenyl)-\n",
    "4H-chromen-4-one"
  )
  x
}

spearman_test <- function(x, y) {
  rho <- cor(x, y, method = "spearman")
  p <- matrix(NA_real_, ncol(x), ncol(y),
              dimnames = list(colnames(x), colnames(y)))

  for (i in seq_len(ncol(x))) {
    for (j in seq_len(ncol(y))) {
      ok <- complete.cases(x[, i], y[, j])
      xx <- x[ok, i]
      yy <- y[ok, j]
      if (length(xx) < 3 || sd(xx) == 0 || sd(yy) == 0) next
      p[i, j] <- suppressWarnings(
        cor.test(xx, yy, method = "spearman")$p.value
      )
    }
  }
  list(rho = rho, p = p)
}

significance_stars <- function(p) {
  out <- matrix("", nrow(p), ncol(p), dimnames = dimnames(p))
  out[!is.na(p) & p < 0.05] <- "*"
  out[!is.na(p) & p < 0.01] <- "**"
  out[!is.na(p) & p < 0.001] <- "***"
  out
}

# 4. Import and select features ----------------------------------------------
otu <- read_input(FILES["bacteria"])
metab <- cbind(read_input(FILES["neg"]), read_input(FILES["pos"]))

samples <- rownames(otu)[rownames(otu) %in% rownames(metab)]
if (length(samples) < 3) stop("Fewer than three shared samples were found.")

otu <- otu[samples, , drop = FALSE]
metab <- metab[samples, , drop = FALSE]

if (!all(vapply(otu, is.numeric, logical(1))) ||
    !all(vapply(metab, is.numeric, logical(1)))) {
  stop("All bacterial and metabolite columns must be numeric.")
}

message("Shared samples: ", length(samples))

otu_var <- vapply(otu, var, numeric(1), na.rm = TRUE)
metab_var <- vapply(metab, var, numeric(1), na.rm = TRUE)
otu_top <- names(sort(otu_var, decreasing = TRUE))[seq_len(min(TOP_N, ncol(otu)))]
metab_top <- names(sort(metab_var, decreasing = TRUE))[
  seq_len(min(TOP_N, ncol(metab)))
]

otu_top_data <- otu[, otu_top, drop = FALSE]
metab_top_data <- metab[, metab_top, drop = FALSE]
cor_result <- spearman_test(otu_top_data, metab_top_data)
rho <- cor_result$rho
p_value <- cor_result$p

# 5. Figure 6A | Heatmap ------------------------------------------------------
keep_rows <- rowSums(p_value < P_HEATMAP, na.rm = TRUE) > 0
keep_cols <- colSums(p_value < P_HEATMAP, na.rm = TRUE) > 0
heat_rho <- rho[keep_rows, keep_cols, drop = FALSE]
heat_p <- p_value[keep_rows, keep_cols, drop = FALSE]

if (nrow(heat_rho) > 1 && ncol(heat_rho) > 1) {
  heat <- pheatmap(
    heat_rho, cluster_rows = TRUE, cluster_cols = TRUE,
    color = colorRampPalette(c(NEG_COL, MID_COL, POS_COL))(100),
    fontsize_row = 8, fontsize_col = 7.4, angle_col = 45,
    show_rownames = TRUE, show_colnames = TRUE,
    labels_col = format_metabolite(colnames(heat_rho)),
    display_numbers = significance_stars(heat_p),
    number_color = "black", fontsize_number = 7,
    border_color = NA, silent = TRUE
  )

  label_row <- which(heat$gtable$layout$name == "col_names")
  if (length(label_row) == 1) {
    heat$gtable$heights[label_row] <- unit(7.6, "cm")
  }

  tiff(file.path(OUT_DIR, "Figure6A_Correlation_Heatmap.tiff"),
       width = 14, height = 12, units = "in", res = 600,
       compression = "lzw")
  grid.newpage()
  grid.draw(heat$gtable)
  dev.off()

  pdf(file.path(OUT_DIR, "Figure6A_Correlation_Heatmap.pdf"),
      width = 14, height = 12, useDingbats = FALSE)
  grid.newpage()
  grid.draw(heat$gtable)
  dev.off()
  message("Figure 6A exported.")
} else {
  warning("Figure 6A skipped: the significant matrix is too small.")
}

# 6. Figure 6B | Sankey plot --------------------------------------------------
rho_df <- as.data.frame(as.table(rho))
p_df <- as.data.frame(as.table(p_value))
names(rho_df) <- c("Bacteria", "Metabolite", "Rho")
names(p_df) <- c("Bacteria", "Metabolite", "P_value")

sankey_df <- left_join(rho_df, p_df,
                       by = c("Bacteria", "Metabolite")) %>%
  filter(!is.na(Rho), !is.na(P_value),
         P_value < P_SANKEY, abs(Rho) > RHO_SANKEY) %>%
  mutate(Direction = ifelse(Rho > 0, "Positive", "Negative"))

if (nrow(sankey_df) > 0) {
  bacteria_order <- sankey_df %>%
    group_by(Bacteria) %>%
    summarise(weight = sum(abs(Rho)), .groups = "drop") %>%
    arrange(weight) %>% pull(Bacteria)

  metabolite_order <- sankey_df %>%
    group_by(Metabolite) %>%
    summarise(weight = sum(abs(Rho)), .groups = "drop") %>%
    arrange(weight) %>% pull(Metabolite)

  sankey_df <- sankey_df %>%
    mutate(
      Bacteria_display = factor(Bacteria, levels = bacteria_order),
      Metabolite_display = factor(
        format_metabolite(Metabolite),
        levels = format_metabolite(metabolite_order)
      )
    )

  message("Figure 6B genera: ", n_distinct(sankey_df$Bacteria))
  message("Figure 6B metabolites: ", n_distinct(sankey_df$Metabolite))
  message("Figure 6B significant links: ", nrow(sankey_df))

  sankey_base <- ggplot(
    sankey_df,
    aes(axis1 = Bacteria_display, axis2 = Metabolite_display, y = abs(Rho))
  ) +
    geom_alluvium(aes(fill = Direction), width = 0.10,
                  alpha = 0.82, color = NA) +
    geom_stratum(width = 0.10, fill = "grey96",
                 color = "black", linewidth = 0.35)

  strata <- ggplot_build(sankey_base)$data[[2]]
  required_cols <- c("x", "y", "ymin", "ymax", "stratum")
  if (!all(required_cols %in% names(strata))) {
    stop("Cannot extract Sankey stratum positions; check ggalluvial version.")
  }

  strata <- mutate(strata, label = as.character(stratum))
  left_labels <- filter(strata, x == 1)
  right_labels <- filter(strata, x == 2) %>%
    arrange(desc(y)) %>%
    mutate(
      label_y = seq(max(strata$ymax, na.rm = TRUE) * 0.985,
                    max(strata$ymax, na.rm = TRUE) * 0.015,
                    length.out = n())
    )

  sankey_plot <- sankey_base +
    geom_segment(
      data = left_labels,
      aes(x = 0.80, xend = 0.95, y = y, yend = y),
      inherit.aes = FALSE, color = "grey55", linewidth = 0.25
    ) +
    geom_text(
      data = left_labels,
      aes(x = 0.77, y = y, label = label),
      inherit.aes = FALSE, hjust = 1, size = 4.6,
      lineheight = 0.95, color = "black"
    ) +
    geom_segment(
      data = right_labels,
      aes(x = 2.05, xend = 2.11, y = y, yend = label_y),
      inherit.aes = FALSE, color = "grey55", linewidth = 0.25
    ) +
    geom_text(
      data = right_labels,
      aes(x = 2.13, y = label_y, label = label),
      inherit.aes = FALSE, hjust = 0, size = 4.6,
      lineheight = 0.9, color = "black"
    ) +
    scale_fill_manual(
      values = c(Negative = NEG_COL, Positive = POS_COL),
      breaks = c("Negative", "Positive")
    ) +
    scale_x_continuous(
      breaks = c(1, 2), labels = c("Bacteria", "Metabolites"),
      limits = c(0.48, 2.86), expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = NULL, labels = NULL,
      expand = expansion(mult = c(0.015, 0.015))
    ) +
    labs(
      title = "Bacteria-Metabolite Correlation Sankey (Top30 metabolites)",
      fill = "Correlation"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(
        size = 12.5, face = "bold", color = "black",
        margin = margin(t = 9)
      ),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(
        size = 15.5, face = "bold", hjust = 0.5,
        margin = margin(b = 10)
      ),
      legend.position = "right",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10.5),
      legend.direction = "vertical",
      legend.box.spacing = unit(0.02, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.key.height = unit(0.58, "cm"),
      legend.key.width = unit(0.58, "cm"),
      plot.margin = margin(10, 4, 10, 6)
    ) +
    coord_cartesian(clip = "off")

  ggsave(
    file.path(OUT_DIR, "Figure6B_Correlation_Sankey.tiff"),
    sankey_plot, width = 24.5, height = 11.5, units = "in",
    dpi = 600, compression = "lzw", limitsize = FALSE
  )
  ggsave(
    file.path(OUT_DIR, "Figure6B_Correlation_Sankey.pdf"),
    sankey_plot, width = 24.5, height = 11.5, units = "in",
    device = "pdf", limitsize = FALSE
  )
  message("Figure 6B exported.")
} else {
  warning("Figure 6B skipped: no links met P < ", P_SANKEY,
          " and |Spearman rho| > ", RHO_SANKEY, ".")
}

# 7. Record software versions -------------------------------------------------
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo.txt"))
message("Analysis completed. Results: ", OUT_DIR)
