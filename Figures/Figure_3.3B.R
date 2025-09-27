# =========================================================
# Figure 3.3B — Ravel 2011: direction-specific stacks at α=0.10
# Encoding: pattern = correctness (Correct = plain, Wrong = striped)
# Bars: fixed fill #00A9E0; stripe color = black
# =========================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggpattern)
})

# ---- Methods (order only; no colors now) ----
method_levels <- c("BALOR","BZIOLR","ORM","ANCOM-BC","LinDA","DESeq2","MaAsLin2","corncob","LDM")

# ---- Load saved results for Ravel ----
load("True_Ravel.RData")

# ---- collect per-taxon calls at α=0.10 from objects ----
calls_list <- list(
  if (exists("res_balor_calls"))   transmute(res_balor_calls,  method = "BALOR",    taxon, call, ground_truth),
  if (exists("res_ziolr_calls"))   transmute(res_ziolr_calls,  method = "BZIOLR",   taxon, call, ground_truth),
  if (exists("res2"))              transmute(res2,             method = "ORM",      taxon, call = res, ground_truth),
  if (exists("res_ancombc"))       transmute(res_ancombc,      method = "ANCOM-BC", taxon, call, ground_truth),
  if (exists("res_linda"))         transmute(res_linda,        method = "LinDA",    taxon, call, ground_truth),
  if (exists("res_deseq2"))        transmute(res_deseq2,       method = "DESeq2",   taxon, call, ground_truth),
  if (exists("res_maaslin"))       transmute(res_maaslin,      method = "MaAsLin2", taxon, call, ground_truth),
  if (exists("res_corncob_2"))     transmute(res_corncob_2,    method = "corncob",  taxon, call, ground_truth),
  if (exists("res_ldm"))           transmute(res_ldm,          method = "LDM",      taxon, call, ground_truth)
)
calls_all_dir <- bind_rows(calls_list)

# ---- stack counts by panel and correctness (pattern) ----
df_counts <- calls_all_dir %>%
  filter(call %in% c("healthy","bv"), ground_truth %in% c("healthy","bv")) %>%
  mutate(
    panel = recode(call,
                   healthy = "Calls toward HV-associated",
                   bv      = "Calls toward BV-associated"),
    correctness = if_else(call == ground_truth, "Correct", "Wrong"),
    method = factor(method, levels = method_levels),
    panel  = factor(panel, levels = c("Calls toward HV-associated",
                                      "Calls toward BV-associated")),
    correctness = factor(correctness, levels = c("Correct","Wrong"))   # Wrong on top
  ) %>%
  count(panel, method, correctness, name = "n") %>%
  tidyr::complete(panel, method,
                  correctness = c("Correct","Wrong"),
                  fill = list(n = 0))

# ---- fixed y-limit shared by both rows ----
y_max <- df_counts %>%
  group_by(panel, method) %>%
  summarise(tot = sum(n), .groups = "drop") %>%
  summarise(max_tot = max(tot)) %>%
  pull(max_tot)
y_lim <- c(0, y_max * 1.08)

# Make legend swatches a bit larger (helps readability)
options(ggpattern_key_scale_factor = 1.25)

# ---- plot: fixed fill; pattern = correctness (Wrong stacked on top) ----
p <- ggplot(
  df_counts,
  aes(x = method, y = n,
      pattern = correctness,
      order = as.numeric(correctness))  # position_stack honors this for draw order
) +
  ggpattern::geom_col_pattern(
    width = 0.72,
    color = "black", linewidth = 0.25,
    position = position_stack(reverse = TRUE),
    fill = "#00A9E0",            # <-- single bar color
    pattern_colour = "black",    # <-- black stripes for 'Wrong'
    pattern_fill   = NA,
    pattern_spacing = 0.03,
    pattern_size    = 0.4,
    pattern_angle   = 45,
    key_glyph = ggpattern::draw_key_polygon_pattern
  ) +
  facet_grid(rows = vars(panel), scales = "fixed") +
  scale_pattern_manual(values = c(Correct = "none", Wrong = "stripe")) +
  scale_y_continuous(limits = y_lim, expand = expansion(mult = c(0, .05))) +
  labs(
    x = "DA methods",
    y = "Number of identified DA features at \u03B1 = 0.10",
    pattern = "Call"
  ) +
  # show ONLY the pattern legend
  guides(
    pattern = guide_legend(
      title = 'Call',
      order = 1,
      byrow = TRUE,
      reverse = TRUE,
      keyheight = grid::unit(18, "pt"),
      keywidth  = grid::unit(30, "pt"),
      override.aes = list(fill = "#00A9E0", colour = "black", linewidth = 0.4, group = NA)
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.spacing.y = grid::unit(8, "pt"),
    legend.key.height = grid::unit(18, "pt"),
    legend.key.width  = grid::unit(30, "pt"),
    axis.title.x = element_text(size = 16)
  )

print(p)
ggsave("Figure 3.3B.png",
       p, width = 11, height = 6.2, units = "in", dpi = 300, bg = "white")
