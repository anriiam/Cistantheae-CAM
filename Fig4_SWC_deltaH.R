# =============================================================================
# Figure 4 reproduction script
# Supplementary material for Jose et. al: Two putative instances of evolutionary
# loss of CAM photosynthesis within the plant clade Cistantheae
# 
# Reproduces Figure 4 (SWC vs CAM expression in Aeonium and Cistantheae) from
# a single species-level data file. Two within-lineage least-squares
# regressions are fitted and displayed on the figure:
#   (1) Aeonium (from Messerschmid et al. 2026 data, water-withholding subset)
#   (2) Cistantheae (this study, water-withholding treatment)
#
# The console output also reports a sensitivity analysis excluding the two
# Montiopsis species (M. gilliesii, M. umbellatum), which is discussed in the
# main text but is not shown as a regression line on Figure 4.
#
# NOTE: A small number of purely cosmetic elements (r/p summary box, threshold
# label + leftward arrow) are added by hand in Illustrator on the exported
# figure. All statistical content, points, regression lines, 95% CI ribbons,
# threshold shaded region and threshold vertical lines are produced by this
# script.
#
# INPUT:   Fig4_species_level_data.csv  (in the same directory as this script)
# OUTPUT:  Fig4_reproduced.pdf, Fig4_reproduced.png  (in the same directory)
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Session setup
# -----------------------------------------------------------------------------
cat("R version:", R.version.string, "\n\n")

required_pkgs <- c("ggplot2", "ggrepel", "dplyr")
missing_pkgs  <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs)
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

if (interactive()) {
  script_dir <- tryCatch(
    dirname(rstudioapi::getSourceEditorContext()$path),
    error = function(e) getwd()
  )
  if (nzchar(script_dir)) setwd(script_dir)
}

# -----------------------------------------------------------------------------
# 1. Load data + standardize display names
# -----------------------------------------------------------------------------
dat <- read.csv("Fig4_species_level_data.csv", stringsAsFactors = FALSE)

stopifnot(
  nrow(dat) >= 30,
  all(c("taxon","clade","genus","SWC_g_per_g","DH_umol_per_g") %in% names(dat)),
  all(dat$clade %in% c("Aeonium","Cistantheae"))
)

# Expand Cistantheae short codes to standardized "G. species" labels
cist_display <- c(
  "Ccrass"        = "C. crassifolia",
  "Caffcymmosa1"  = "C. aff. cymosa 1",
  "Caffcymosa2"   = "C. aff. cymosa 2",
  "Cafflitoralis" = "C. aff. litoralis",
  "Mgli"          = "M. gilliesii",
  "Mumbellatum"   = "M. umbellatum",
  "Pcalycina"     = "P. calycina"
)
dat$display <- ifelse(dat$taxon %in% names(cist_display),
                      cist_display[dat$taxon],
                      dat$taxon)

aeo      <- subset(dat, clade == "Aeonium")
cist_all <- subset(dat, clade == "Cistantheae")
cist_noM <- subset(cist_all, genus != "Montiopsis")   # sensitivity analysis

cat("Sample sizes:\n")
cat("  Aeonium species:                       ", nrow(aeo),      "\n")
cat("  Cistantheae species (all):             ", nrow(cist_all), "\n")
cat("  Cistantheae species (excl. Montiopsis):", nrow(cist_noM), "\n\n")

# -----------------------------------------------------------------------------
# 2. Correlations (Pearson; equivalent to Gaussian GLM with identity link)
# -----------------------------------------------------------------------------
cor_aeo  <- cor.test(aeo$SWC_g_per_g,      aeo$DH_umol_per_g)
cor_cist <- cor.test(cist_all$SWC_g_per_g, cist_all$DH_umol_per_g)
cor_noM  <- cor.test(cist_noM$SWC_g_per_g, cist_noM$DH_umol_per_g)

cat("=== Aeonium (water-withholding treatments) ===\n"); print(cor_aeo)
cat("\n=== Cistantheae (all species, this study) ===\n"); print(cor_cist)
cat("\n=== Sensitivity: Cistantheae excluding Montiopsis ===\n"); print(cor_noM)

# -----------------------------------------------------------------------------
# 3. Plot elements
# -----------------------------------------------------------------------------
# Colours
col_aeonium      <- "#4A75B8"
col_cistanthe    <- "#D45ABE"
col_montiopsis   <- "#7AB54E"
col_phillipiamra <- "#D45757"
col_cist_fit     <- "#882255"

genus_levels  <- c("Aeonium", "Cistanthe", "Montiopsis", "Phillipiamra",
                   "Calyptridium", "Thingia")
genus_colours <- c(Aeonium      = col_aeonium,
                   Cistanthe    = col_cistanthe,
                   Montiopsis   = col_montiopsis,
                   Phillipiamra = col_phillipiamra,
                   Calyptridium = "#DAA520",
                   Thingia      = "#6E5B96")
present_genera <- intersect(genus_levels, unique(dat$genus))

# ----- Threshold region (homogeneous green shade) ------------------------
# Left edge:  abuts the right side of the M. gilliesii circle but does not
#             cross into it (~0.47 data units to the right of its centre;
#             tune if you change the point size).
# Right edge: left edge of the Aeonium 95% CI ribbon (min Aeonium SWC - 0.5).
mgli_x       <- subset(cist_all, taxon == "Mgli")$SWC_g_per_g       # ~3.69
aeo_ci_start <- min(aeo$SWC_g_per_g) - 0.5                           # ~12.07
shade_left   <- mgli_x + 0.47                                        # ~4.16

# Three C. aff. taxa share SWC = 25.85 and overlap heavily; hand-tuned labels
manual_positions <- data.frame(
  taxon   = c("Caffcymmosa1", "Cafflitoralis", "Caffcymosa2"),
  label_x = c(33,            28,             30),
  label_y = c(55,            28,             5)
)
manual_labels <- merge(cist_all, manual_positions, by = "taxon")
auto_cist    <- subset(cist_all, !(taxon %in% manual_positions$taxon))

# -----------------------------------------------------------------------------
# 4. Build the plot
# -----------------------------------------------------------------------------
p <- ggplot() +
  # Homogeneous green threshold region
  annotate("rect", xmin = shade_left, xmax = aeo_ci_start,
           ymin = -Inf, ymax = Inf,
           fill = "#8CCC8C", alpha = 0.22) +
  # Two bold dashed vertical lines at the threshold boundaries
  geom_vline(xintercept = c(shade_left, aeo_ci_start),
             color = "#555555", linetype = "dashed",
             alpha = 0.70, linewidth = 0.9) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dotted") +
  
  # Regressions with 95% CI
  geom_smooth(data = aeo,
              aes(SWC_g_per_g, DH_umol_per_g),
              method = "lm", formula = y ~ x, se = TRUE, level = 0.95,
              color = col_aeonium, fill = col_aeonium,
              alpha = 0.15, linetype = "dashed", linewidth = 1,
              show.legend = FALSE) +
  geom_smooth(data = cist_all,
              aes(SWC_g_per_g, DH_umol_per_g),
              method = "lm", formula = y ~ x, se = TRUE, level = 0.95,
              color = col_cist_fit, fill = col_cist_fit,
              alpha = 0.13, linetype = "dashed", linewidth = 1,
              show.legend = FALSE) +
  
  # Points, coloured by genus
  geom_point(data = dat,
             aes(SWC_g_per_g, DH_umol_per_g, fill = genus),
             color = "black", shape = 21, size = 3.3, stroke = 0.5) +
  scale_fill_manual(
    values = genus_colours[present_genera],
    breaks = present_genera,
    labels = lapply(present_genera, function(g) bquote(italic(.(g)))),
    name   = NULL
  ) +
  
  # Aeonium labels (short codes, auto-placed)
  geom_text_repel(data = aeo,
                  aes(SWC_g_per_g, DH_umol_per_g, label = display),
                  fontface = "italic", size = 2.8, color = "#2C3E50",
                  max.overlaps = Inf, segment.size = 0.2) +
  # Cistantheae labels (auto-placed subset)
  geom_text_repel(data = auto_cist,
                  aes(SWC_g_per_g, DH_umol_per_g, label = display),
                  fontface = "bold.italic", size = 3.3, color = "#7A2A6B",
                  max.overlaps = Inf, segment.size = 0.2) +
  # Cistantheae labels (manual subset with leader lines)
  geom_segment(data = manual_labels,
               aes(x = SWC_g_per_g, y = DH_umol_per_g,
                   xend = label_x, yend = label_y),
               color = "#7A2A6B", linewidth = 0.35, alpha = 0.6) +
  geom_text(data = manual_labels,
            aes(x = label_x, y = label_y, label = display),
            fontface = "bold.italic", size = 3.3, color = "#7A2A6B",
            hjust = 0.5, vjust = 0.5) +
  
  # r / p annotation (top-right)
  annotate("text", x = Inf, y = Inf,
           label = sprintf("italic(Aeonium)*': r = %+.2f, p = %.3f'",
                           cor_aeo$estimate, cor_aeo$p.value),
           parse = TRUE, hjust = 1.05, vjust = 1.5, size = 3.5) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("'Cistantheae:  r = %+.2f, p = %.3f'",
                           cor_cist$estimate, cor_cist$p.value),
           parse = TRUE, hjust = 1.05, vjust = 3.0, size = 3.5) +
  
  labs(
    title    = expression("SWC vs CAM expression in " *
                          italic("Aeonium") * " vs Cistantheae"),
    subtitle = expression("(water-withholding treatments only for " *
                          italic("Aeonium") * ", with 95% CI)"),
    x = expression(bold(SWC~(g~g^{-1}))),
    y = expression(bold(Delta*H^{"+"}~(mu*mol~g^{-1})))
  ) +
  coord_cartesian(xlim = c(0, 42), ylim = c(-25, 150)) +
  theme_classic(base_size = 12) +
  theme(
    plot.title           = element_text(face = "bold"),
    plot.subtitle        = element_text(face = "bold"),
    panel.grid.major.y   = element_line(color = "#E8E8E8", linewidth = 0.4),
    legend.position      = c(0.98, 0.89),
    legend.justification = c(1, 1),
    legend.background    = element_rect(fill = "white", color = "#CCCCCC"),
    legend.title         = element_blank(),
    legend.text.align    = 0
  )

print(p)
ggsave("Fig4_reproduced.pdf", p, width = 11, height = 8)
ggsave("Fig4_reproduced.png", p, width = 11, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# 5. Session info
# -----------------------------------------------------------------------------
cat("\n=== sessionInfo() ===\n")
print(sessionInfo())
