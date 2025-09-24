# =========================================================
# Figure 3.3A — Ravel 2011: Correct vs Wrong across α 
# =========================================================
library(tidyverse)
library(patchwork)

# ---- Methods ----
method_levels <- c("BALOR","BZIOLR","ORM","ANCOM-BC","LinDA","DESeq2","MaAsLin2","corncob","LDM")

# ---- Color palette ----
method_colors <- c(
  BALOR     = "#045275",   
  BZIOLR    = "#089099",   
  ORM       = "#7CCBA2",  
  "ANCOM-BC"= "#FCDE9C", 
  LinDA     = "#F0746E", 
  DESeq2    = "#7C1D6F", 
  MaAsLin2  = "#ED85B0", 
  corncob   = "#A3A500FF", 
  LDM       = "#EEB849"  
)

# ---- Load saved results for Ravel ----
load("ravel_040925_data.RData")
grid_all <- side_by_side_grid

# ---- Harmonize method names ----
grid_all <- grid_all %>%
  mutate(model = tolower(model)) %>%
  mutate(model = recode(
    model,
    "balor"       = "BALOR",
    "ziolr"      = "BZIOLR",
    "orm"         = "ORM",
    "ancombc"     = "ANCOM-BC",
    "linda"       = "LinDA",
    "deseq2"      = "DESeq2",
    "maaslin2"    = "MaAsLin2",
    "corncob_lrt" = "corncob",
    "ldm"         = "LDM"
  ))

alpha_keep <- c(0.05, 0.10, 0.20)

df <- grid_all %>%
  filter(alpha %in% alpha_keep) %>%
  mutate(
    method  = factor(model, levels = method_levels),
    alpha_f = factor(alpha, levels = alpha_keep,
                     labels = c("\u03B1 = 0.05", "\u03B1 = 0.10", "\u03B1 = 0.20"))
  ) %>%
  filter(!is.na(method))

# Limits (independent scales per row) 
tp_max <- max(df$TP, na.rm = TRUE); tp_lim <- c(0, tp_max * 1.05)

bar_params <- list(width = 0.72, color = "black", linewidth = 0.25, alpha = 0.90)

facet_opt <- facet_grid(. ~ alpha_f, scales = "fixed")  # or "free_y" in BOTH figures

# ---------- TOP: TP row ----------
p_tp <- ggplot(df, aes(x = method, y = TP, fill = method)) +
  do.call(geom_col, bar_params) +         
  facet_opt +
  scale_fill_manual(values = method_colors, limits = method_levels, drop = FALSE) +
  scale_y_continuous(limits = tp_lim, expand = expansion(mult = c(0, .05))) +
  labs(y = "Correct Calls (TP)", x = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.4),
    panel.border = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(5.5, 5.5, 0, 5.5),
    axis.title.y = element_text(size = 10),
    axis.title.x = element_text(size = 16)
  )

# ---------- BOTTOM: FP row ----------
p_fp <- ggplot(df, aes(x = method, y = FP, fill = method)) +
  do.call(geom_col, bar_params) +         
  facet_opt +
  scale_fill_manual(values = method_colors, limits = method_levels, drop = FALSE, name = NULL) +
  scale_y_reverse(limits = c(2,0), breaks = 0:2,
                  expand = expansion(mult = c(0, .05))) + 
  labs(y = "Wrong Calls (FP)", x = "DA Methods") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.4),
    panel.border = element_blank(),
    strip.text = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(0, 5.5, 5.5, 5.5),
    axis.title.y = element_text(size = 10),
    axis.title.x = element_text(size = 16)
  )

gap_pt <- 2
p_tp <- p_tp + theme(plot.margin = margin(5.5, 5.5, gap_pt, 5.5))
p_fp <- p_fp + theme(plot.margin = margin(gap_pt, 5.5, 5.5, 5.5))

# ---------- Compose ----------
fig_4_1_ravel <- p_tp / p_fp + plot_layout(heights = c(2, 1))
print(fig_4_1_ravel)

# ---- Save ----
ggsave("Figure_3.3A.png", fig_4_1_ravel,
       width = 11, height = 6.2, units = "in", dpi = 300)
