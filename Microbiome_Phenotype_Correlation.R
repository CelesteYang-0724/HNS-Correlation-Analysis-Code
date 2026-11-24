# =============================================================================
# Microbiome-phenotype Spearman correlation analysis
# Public-code release | R >= 4.5.2 | Base R only
# Input data are NOT distributed in the public repository.
# Run from the repository root; put the two required CSV files in ./data/.
# =============================================================================

# 1. Configuration ------------------------------------------------------------
DATA_DIR <- "data"
OUT_DIR <- "results"
MICROBIOME_FILE <- "属水平肠道菌群n=6.csv"
PHENOTYPE_FILE <- "炎症因子及行为学指标.csv"
CSV_ENCODING <- ""  # System default. For GBK-encoded CSVs, use "GBK".

P_THRESHOLD <- 0.05          # Unadjusted P values, as in the original script
COR_LIMIT <- 0.60           # Color saturation range: -0.6 to +0.6
FIGURE_WIDTH_IN <- NA_real_ # NA: determine size automatically
FIGURE_HEIGHT_IN <- NA_real_
LABEL_CEX <- 0.90           # Axis label size relative to device point size
STAR_CEX <- 1.15            # Significance annotation size

NEG_COL <- "#2166ac"
MID_COL <- "#f7f7f7"
POS_COL <- "#b2182b"

# 2. Import and validate inputs ----------------------------------------------
read_numeric_table <- function(path) {
  if (!file.exists(path)) {
    stop("Missing input file: ", path,
         "\nCreate ./data/ and place the required CSV files there.")
  }
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                 fileEncoding = CSV_ENCODING)
  if (ncol(df) < 2L || nrow(df) < 3L) {
    stop("Each CSV must contain sample IDs and at least one numeric column: ", path)
  }

  ids <- trimws(as.character(df[[1L]]))
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("The first column must contain nonempty, unique sample IDs: ", path)
  }

  values <- df[-1L]
  for (j in seq_along(values)) {
    original <- values[[j]]
    numeric_values <- suppressWarnings(as.numeric(as.character(original)))
    bad <- !is.na(original) & nzchar(trimws(as.character(original))) &
      is.na(numeric_values)
    if (any(bad)) {
      stop("Nonnumeric value(s) in column '", names(values)[j],
           "' of ", path, ". Check the input table.")
    }
    values[[j]] <- numeric_values
  }
  if (anyDuplicated(names(values))) {
    stop("Duplicate feature/phenotype column names detected: ", path)
  }
  rownames(values) <- ids
  values
}

microbiome <- read_numeric_table(file.path(DATA_DIR, MICROBIOME_FILE))
phenotype <- read_numeric_table(file.path(DATA_DIR, PHENOTYPE_FILE))

# Match by explicit sample ID; never truncate tables or pair rows by position.
shared <- rownames(microbiome)[rownames(microbiome) %in% rownames(phenotype)]
if (length(shared) < 3L) {
  stop("Fewer than three shared sample IDs. Check the first column of both CSVs.")
}
if (length(shared) != nrow(microbiome) ||
    length(shared) != nrow(phenotype)) {
  warning("Only ", length(shared), " shared samples will be analyzed. ",
          "Unmatched sample IDs are excluded.")
}
microbiome <- microbiome[shared, , drop = FALSE]
phenotype <- phenotype[shared, , drop = FALSE]
stopifnot(identical(rownames(microbiome), rownames(phenotype)))

# Display names only: no source CSV is modified.
phenotype_labels <- colnames(phenotype)
phenotype_labels[phenotype_labels %in% c("TNF_alpha", "TNF-α", "TNF-alpha")] <- "TNF-alpha"
phenotype_labels[phenotype_labels %in% c("IL-1β", "IL-1beta")] <- "IL-1beta"
phenotype_labels[phenotype_labels %in% c("Y-Maze", "Y-maze")] <- "Y-maze"
if (anyDuplicated(phenotype_labels)) {
  stop("Phenotype labels are duplicated after standardization. Rename columns explicitly.")
}
colnames(phenotype) <- phenotype_labels
message("Shared samples: ", length(shared), "; bacterial genera: ",
        ncol(microbiome), "; phenotypes: ", ncol(phenotype))

# 3. Spearman correlations and P values --------------------------------------
cor_matrix <- matrix(NA_real_, ncol(microbiome), ncol(phenotype),
                     dimnames = list(colnames(microbiome), colnames(phenotype)))
p_matrix <- cor_matrix
n_matrix <- cor_matrix

for (i in seq_len(ncol(microbiome))) {
  for (j in seq_len(ncol(phenotype))) {
    x <- microbiome[[i]]
    y <- phenotype[[j]]
    valid <- is.finite(x) & is.finite(y)
    n_matrix[i, j] <- sum(valid)
    if (sum(valid) < 3L || length(unique(x[valid])) < 2L ||
        length(unique(y[valid])) < 2L) next

    # Keep cor.test's default exact/tie handling from the source script.
    test <- suppressWarnings(cor.test(x[valid], y[valid], method = "spearman"))
    cor_matrix[i, j] <- unname(test$estimate)
    p_matrix[i, j] <- test$p.value
  }
}
if (all(is.na(cor_matrix))) stop("No valid Spearman correlations could be calculated.")

stars <- matrix("", nrow(p_matrix), ncol(p_matrix), dimnames = dimnames(p_matrix))
stars[!is.na(p_matrix) & p_matrix < 0.05] <- "*"
stars[!is.na(p_matrix) & p_matrix < 0.01] <- "**"
stars[!is.na(p_matrix) & p_matrix < 0.001] <- "***"

# 4. Export complete results and significant associations --------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
all_results <- data.frame(
  Microbiota = rep(rownames(cor_matrix), times = ncol(cor_matrix)),
  Phenotype = rep(colnames(cor_matrix), each = nrow(cor_matrix)),
  Correlation = as.vector(cor_matrix),
  P_value = as.vector(p_matrix),
  N_complete = as.integer(as.vector(n_matrix)),
  Significance = as.vector(stars),
  stringsAsFactors = FALSE
)
all_results$P_adjust_BH <- p.adjust(all_results$P_value, method = "BH")
all_results <- all_results[order(all_results$P_value, na.last = TRUE), ]
significant_results <- subset(all_results, !is.na(P_value) & P_value < P_THRESHOLD)

write.csv(all_results, file.path(OUT_DIR, "all_correlations.csv"),
          row.names = FALSE, na = "")
write.csv(significant_results,
          file.path(OUT_DIR, "significant_correlations.csv"),
          row.names = FALSE, na = "")
write.csv(cor_matrix, file.path(OUT_DIR, "correlation_matrix.csv"), na = "")
write.csv(p_matrix, file.path(OUT_DIR, "p_value_matrix.csv"), na = "")
writeLines(capture.output(sessionInfo()),
           con = file.path(OUT_DIR, "sessionInfo.txt"))
message("P < ", P_THRESHOLD, " (uncorrected) associations: ", nrow(significant_results))

# 5. Draw a heatmap with explicit correlation color scale --------------------
plot_heatmap <- function() {
  nc <- ncol(cor_matrix)
  nr <- nrow(cor_matrix)
  palette <- colorRampPalette(c(NEG_COL, "#92c5de", MID_COL,
                                "#f4a582", POS_COL))(100)
  breaks <- seq(-COR_LIMIT, COR_LIMIT, length.out = length(palette) + 1L)
  z <- pmax(-COR_LIMIT, pmin(COR_LIMIT, cor_matrix))

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  layout(matrix(c(1, 2), nrow = 1), widths = c(5.2, 1))
  par(mar = c(max(9, min(15, 4 + nc * 0.65)),
              max(8, min(14, 4 + nr * 0.35)), 2.0, 1.0))

  # Matrix rows are bacteria; matrix columns are phenotypes.
  image(x = seq_len(nc), y = seq_len(nr), z = t(z),
        col = palette, breaks = breaks, axes = FALSE,
        xlab = "", ylab = "", xlim = c(0.5, nc + 0.5),
        ylim = c(0.5, nr + 0.5), useRaster = FALSE)
  axis(1, at = seq_len(nc), labels = colnames(cor_matrix),
       las = 2, tick = FALSE, cex.axis = LABEL_CEX)
  axis(2, at = seq_len(nr), labels = rownames(cor_matrix),
       las = 2, tick = FALSE, cex.axis = LABEL_CEX)
  abline(h = seq(0.5, nr + 0.5, by = 1), col = "white", lwd = 0.7)
  abline(v = seq(0.5, nc + 0.5, by = 1), col = "white", lwd = 0.7)

  for (i in seq_len(nr)) for (j in seq_len(nc)) {
    if (nzchar(stars[i, j])) {
      text(j, i, stars[i, j], cex = STAR_CEX, font = 2, col = "black")
    }
  }

  par(mar = c(5, 0.5, 3, 3.5))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(-COR_LIMIT, COR_LIMIT))
  for (k in seq_along(palette)) {
    rect(0.12, breaks[k], 0.35, breaks[k + 1L],
         col = palette[k], border = NA)
  }
  axis(4, at = seq(-COR_LIMIT, COR_LIMIT, length.out = 7),
       pos = 0.35, las = 1, cex.axis = 0.9)
  mtext("Spearman rho", side = 3, line = 0.5, cex = 0.95)
}

# 6. Export figure; do not display in RStudio -------------------------------
auto_width <- max(12, ncol(cor_matrix) * 0.80 + 2)
auto_height <- max(9, nrow(cor_matrix) * 0.60)
fig_width <- if (is.na(FIGURE_WIDTH_IN)) auto_width else FIGURE_WIDTH_IN
fig_height <- if (is.na(FIGURE_HEIGHT_IN)) auto_height else FIGURE_HEIGHT_IN

pdf(file.path(OUT_DIR, "correlation_heatmap.pdf"),
    width = fig_width, height = fig_height, useDingbats = FALSE)
plot_heatmap()
dev.off()

png(file.path(OUT_DIR, "correlation_heatmap_600dpi.png"),
    width = fig_width, height = fig_height, units = "in", res = 600,
    type = "cairo-png", pointsize = 12)
plot_heatmap()
dev.off()

message("Analysis complete. All files saved in: ", normalizePath(OUT_DIR))
