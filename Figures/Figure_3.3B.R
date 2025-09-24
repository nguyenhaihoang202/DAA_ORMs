# =========================================================
# Figure 3.3B — Ravel 2011: direction-specific stacks at α=0.10
# Encoding: color = method, pattern = correctness
# Legend: pattern only (Correct = plain, Wrong = striped)
# =========================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggpattern)
})

# ---- Methods & colors ----
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

# ---- Load saved results for Ravel ----
load("ravel_040925_data.RData")

# ---- collect per-taxon calls at α=0.10 from objects ----
calls_list <- list(
  "BALOR"    = res_balor_calls |> transmute(method = "BALOR",    taxon, call, ground_truth),
  "BZIOLR"   = res_ziolr_calls |> transmute(method = "BZIOLR",   taxon, call, ground_truth),
  "ORM"      = res2            |> transmute(method = "ORM",      taxon, call = res, ground_truth),
  "ANCOM-BC" = res_ancombc     |> transmute(method = "ANCOM-BC", taxon, call, ground_truth),
  "LinDA"    = res_linda       |> transmute(method = "LinDA",    taxon, call, ground_truth),
  "DESeq2"   = res_deseq2      |> transmute(method = "DESeq2",   taxon, call, ground_truth),
  "MaAsLin2" = res_maaslin     |> transmute(method = "MaAsLin2", taxon, call, ground_truth),
  "corncob"  = res_corncob     |> transmute(method = "corncob",  taxon, call, ground_truth),
  "LDM"      = res_ldm         |> transmute(method = "LDM",      taxon, call, ground_truth)
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
    # ensure Correct at bottom, Wrong on top:
    correctness = factor(correctness, levels = c("Correct","Wrong"))
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

# ---- plot: color = method; pattern = correctness (Wrong stacked on top) ----
p <- ggplot(
  df_counts,
  aes(x = method, y = n,
      fill = method,
      pattern = correctness,
      order = as.numeric(correctness))  # position_stack honors this for draw order
) +
  ggpattern::geom_col_pattern(
    width = 0.72,
    color = "black", linewidth = 0.25, position = position_stack(reverse = TRUE),  # last group drawn sits on top
    # clearer pattern settings
    pattern_colour = "white",
    pattern_fill   = NA,      # keep method color visible under stripes
    pattern_spacing = 0.03,   # wider spacing -> more legible
    pattern_size    = 0.4,    # thicker lines
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
  # show ONLY the pattern legend
  guides(
      fill = "none",
      pattern = guide_legend(
        title = 'Call',
        order = 1,
        byrow = TRUE,
        reverse = TRUE,
        keyheight = grid::unit(18, "pt"),
        keywidth  = grid::unit(30, "pt"),
        override.aes = list(fill = "grey38", group = 1),
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
p <- p +
  theme(
    axis.title.x = element_text(size = 16)
  )
print(p)
ggsave("Figure_3.3B.png",
       p, width = 11, height = 6.2, units = "in", dpi = 300, bg = "white")
