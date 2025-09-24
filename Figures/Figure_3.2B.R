# =========================================================
# Figure 3.2B — HMP Gingival Subset: direction-specific stacks at α=0.10
# Encoding: color = method , pattern = correctness
# Legend: pattern only (Correct = plain, Wrong = striped)
# =========================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

# Install + load ggpattern if needed
if (!requireNamespace("ggpattern", quietly = TRUE)) {
  install.packages("ggpattern")
}
library(ggpattern)

# ---- Methods & colors (same as Fig 4.1) ----
method_levels <- c("BALOR","BZIOLR","ORM","ANCOM-BC","LinDA","DESeq2","MaAsLin2","corncob","LDM")
method_colors <- c(
  BALOR      = "#045275",
  BZIOLR     = "#089099",
  ORM        = "#7CCBA2",
  "ANCOM-BC" = "#FCDE9C",
  LinDA      = "#F0746E",
  DESeq2     = "#7C1D6F",
  MaAsLin2   = "#ED85B0",
  corncob    = "#A3A500FF",
  LDM        = "#EEB849"
)

# ---- Load saved results for HMP ----
load("hmp_v35subset_140925_data.RData")

# ---- collect per-taxon calls at α=0.10 ----
res_bziolr_calls <- get0("res_bziolr_calls", ifnotfound = get0("res_ziolr_calls", ifnotfound = NULL))
pieces <- list(
  if (exists("res_balor_calls"))   transmute(res_balor_calls,  method = "BALOR",    taxon, call, ground_truth),
  if (!is.null(res_bziolr_calls))  transmute(res_bziolr_calls, method = "BZIOLR",   taxon, call, ground_truth),
  if (exists("res2"))              transmute(res2,             method = "ORM",      taxon, call = res, ground_truth),
  if (exists("res_ancombc"))       transmute(res_ancombc,      method = "ANCOM-BC", taxon, call, ground_truth),
  if (exists("res_linda"))         transmute(res_linda,        method = "LinDA",    taxon, call, ground_truth),
  if (exists("res_deseq2"))        transmute(res_deseq2,       method = "DESeq2",   taxon, call, ground_truth),
  if (exists("res_maaslin"))       transmute(res_maaslin,      method = "MaAsLin2", taxon, call, ground_truth),
  if (exists("res_corncob"))       transmute(res_corncob,      method = "corncob",  taxon, call, ground_truth),
  if (exists("res_ldm"))           transmute(res_ldm,          method = "LDM",      taxon, call, ground_truth)
)
calls_all_dir <- bind_rows(pieces)

# ---- counts by panel × method × correctness ----
df_counts <- calls_all_dir %>%
  filter(call %in% c("supragingival","subgingival"),
         ground_truth %in% c("supragingival","subgingival")) %>%
  mutate(
    panel = recode(call,
                   supragingival = "Calls toward Supragingival",
                   subgingival   = "Calls toward Subgingival"),
    correctness = if_else(call == ground_truth, "Correct", "Wrong"),
    method = factor(method, levels = method_levels),
    panel  = factor(panel, levels = c("Calls toward Supragingival",
                                      "Calls toward Subgingival")),
    correctness = factor(correctness, levels = c("Correct","Wrong"))   # Wrong on top
  ) %>%
  count(panel, method, correctness, name = "n") %>%
  tidyr::complete(panel, method, correctness = c("Correct","Wrong"), fill = list(n = 0))

# ---- shared y-limit ----
y_max <- df_counts %>%
  group_by(panel, method) %>% summarise(tot = sum(n), .groups = "drop") %>%
  summarise(max_tot = max(tot)) %>% pull(max_tot)
y_lim <- c(0, y_max * 1.08)

# Make legend swatches a bit larger
options(ggpattern_key_scale_factor = 1.25)

# ---- plot (color = method; pattern = correctness) ----
p_hmp_dir <- ggplot(
  df_counts,
  aes(x = method, y = n,
      fill = method,
      group = correctness,
      pattern = correctness,
      order = as.numeric(correctness))    # ensures drawing order within stacks
) +
  ggpattern::geom_col_pattern(
    position = position_stack(reverse = TRUE),  # with levels Correct,Wrong => Wrong on top
    width = 0.72,
    color = "black", linewidth = 0.25,
    # match Fig 4.1 look
    pattern_colour = "white",
    pattern_fill   = NA,
    pattern_spacing = 0.03,
    pattern_size    = 0.4,
    pattern_angle   = 45,
    key_glyph = ggpattern::draw_key_polygon_pattern
  ) +
  facet_grid(rows = vars(panel), scales = "fixed") +
  scale_fill_manual(values = method_colors, breaks = method_levels, drop = FALSE) +
  scale_pattern_manual(values = c(Correct = "none", Wrong = "stripe")) +
  scale_y_continuous(limits = y_lim, expand = expansion(mult = c(0, .05))) +
  labs(
    x = "DA methods",
    y = "Number of identified DA features at \u03B1 = 0.10",
    pattern = "Result"
  ) +
  # Show ONLY pattern legend; keep boxes separate
  guides(
    fill = "none",
    pattern = guide_legend(
      title = 'Call',
      order = 1, byrow = TRUE, reverse = TRUE,
      keyheight = grid::unit(18, "pt"),
      keywidth  = grid::unit(30, "pt"),
      override.aes = list(fill = "grey38", colour = "black", linewidth = 0.4, group = NA)
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
    legend.key.width  = grid::unit(30, "pt")
  )
p_hmp_dir <- p_hmp_dir +
  theme(
    axis.title.x = element_text(size = 16)
  )
print(p_hmp_dir)
ggsave("Figure_3.2B.png",
       p_hmp_dir, width = 11, height = 6, units = "in", dpi = 300, bg = "white")
