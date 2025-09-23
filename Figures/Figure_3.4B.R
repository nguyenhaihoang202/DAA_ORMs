## ===== 0) Setup =====
library(tidyverse)

files_real <- c(
  Sys.glob(file.path("ravel2011_*.RData")),
  Sys.glob(file.path("ravel2011_*.rds")),
  Sys.glob(file.path("hmp_v35subset_*.RData")),
  Sys.glob(file.path("hmp_v35subset_*.rds"))
) %>% unique()

## ===== 1) Extract elapsed time =====
extract_time_any <- function(obj) {
  out <- list(total_sec = NA_real_, warmup_sec = NA_real_, sample_sec = NA_real_, note = "not found")
  
  if (is.environment(obj) && exists("time", envir = obj, inherits = FALSE)) {
    tt <- get("time", envir = obj)
    if (is.function(tt)) {
      res <- try(tt(), silent = TRUE)
      if (!inherits(res, "try-error") && is.list(res)) {
        out$total_sec  <- as.numeric(res$total  %||% res$total_sec)
        if (!is.null(res$chains) && is.data.frame(res$chains)) {
          out$warmup_sec <- suppressWarnings(sum(res$chains$warmup, na.rm = TRUE))
          out$sample_sec <- suppressWarnings(sum(res$chains$sampling, na.rm = TRUE))
        }
        out$note <- "CmdStanR $time()"
        return(out)
      }
    }
    tt <- try(get("time", envir = obj), silent = TRUE)
    if (is.list(tt)) {
      out$total_sec  <- as.numeric(tt$total  %||% tt$total_sec)
      out$note <- "CmdStanR $time list"
      return(out)
    }
  }
  
  # A list wrapper with $fit inside
  if (is.list(obj) && !is.null(obj$fit)) {
    return(extract_time_any(obj$fit))
  }
  
  # Generic containers 
  if (is.list(obj) && !is.null(obj$time) && is.list(obj$time)) {
    t <- obj$time
    out$total_sec  <- as.numeric(t$total %||% t$total_sec)
    out$warmup_sec <- as.numeric(t$warmup %||% t$warmup_sec)
    out$sample_sec <- as.numeric(t$sampling %||% t$sample_sec)
    out$note <- "found in $time"
    return(out)
  }
  for (nm in c("runtime","timing","times","elapsed")) {
    if (is.list(obj) && !is.null(obj[[nm]])) {
      t <- obj[[nm]]
      if (is.list(t)) {
        out$total_sec  <- as.numeric(t$total %||% t$total_sec %||% t$elapsed)
        out$warmup_sec <- as.numeric(t$warmup %||% t$warmup_sec)
        out$sample_sec <- as.numeric(t$sampling %||% t$sample_sec)
        out$note <- paste0("found in $", nm)
        return(out)
      } else if (is.numeric(t) && length(t) == 1L) {
        out$total_sec <- as.numeric(t); out$note <- paste0("numeric $", nm)
        return(out)
      }
    }
  }
  
  # start/end timestamps
  if (!is.null(obj$start_time) && !is.null(obj$end_time)) {
    out$total_sec <- as.numeric(difftime(obj$end_time, obj$start_time, units = "secs"))
    out$note <- "end_time - start_time"
    return(out)
  }
  
  # attribute
  for (att in c("elapsed_sec","time_total_sec","runtime_sec")) {
    if (!is.null(attr(obj, att))) {
      out$total_sec <- as.numeric(attr(obj, att)); out$note <- paste0("attr ", att)
      return(out)
    }
  }
  
  out
}

## ===== 2) Find BALOR/BZIOLR fits within each file and extract =====
pluck_fit <- function(env, candidates) {
  for (nm in candidates) {
    if (exists(nm, envir = env, inherits = FALSE)) {
      obj <- get(nm, envir = env)
      # direct fit
      if (is.environment(obj) && exists("time", envir = obj, inherits = FALSE)) return(obj)
      # list wrapping a fit
      if (is.list(obj) && !is.null(obj$fit)) {
        f <- obj$fit
        if (is.environment(f) && exists("time", envir = f, inherits = FALSE)) return(f)
      }
      # raw list with $time
      if (is.list(obj) && !is.null(obj$time)) return(obj)
    }
  }
  NULL
}

read_one_file <- function(fp) {
  env <- new.env(parent = emptyenv())
  loaded <- try(load(fp, envir = env), silent = TRUE)
  if (inherits(loaded, "try-error")) {
    return(tibble(file = fp, dataset = "Unknown", run_id = NA_integer_,
                  method = NA_character_, total_sec = NA_real_,
                  warmup_sec = NA_real_, sample_sec = NA_real_,
                  note = "load() failed"))
  }
  
  base <- basename(fp)
  dataset <- case_when(
    str_detect(base, regex("^ravel2011", ignore_case = TRUE)) ~ "Ravel_2011_16S_BV",
    str_detect(base, regex("^hmp_v35subset", ignore_case = TRUE)) ~ "HMP_2012_gingival_V35_subset",
    TRUE ~ "Unknown"
  )
  run_id <- str_extract(base, "(?<=_)\\d+") %>% as.integer()
  
  fit_balor <- pluck_fit(env, c("fit_bayes4","fit_bayes","fit_balor"))
  fit_bziolr <- pluck_fit(env, c("fit_zi","fit_ziolr"))
  
  rows <- list()
  if (!is.null(fit_balor)) {
    t <- extract_time_any(fit_balor)
    rows <- append(rows, list(tibble(
      file = fp, dataset, run_id, method = "BALOR",
      total_sec = t$total_sec, warmup_sec = t$warmup_sec, sample_sec = t$sample_sec, note = t$note
    )))
  }
  if (!is.null(fit_bziolr)) {
    t <- extract_time_any(fit_bziolr)
    rows <- append(rows, list(tibble(
      file = fp, dataset, run_id, method = "BZIOLR",
      total_sec = t$total_sec, warmup_sec = t$warmup_sec, sample_sec = t$sample_sec, note = t$note
    )))
  }
  
  if (!length(rows)) {
    rows <- list(tibble(
      file = fp, dataset, run_id, method = NA_character_,
      total_sec = NA_real_, warmup_sec = NA_real_, sample_sec = NA_real_,
      note = "fit_bayes4/fit_zi not found"
    ))
  }
  bind_rows(rows)
}

rt_real <- files_real %>% map_dfr(read_one_file) %>%
  mutate(
    total_min  = total_sec / 60,
    total_hour = total_sec / 3600,
    method = factor(method, levels = c("BALOR","BZIOLR"))
  )

## save to CSV for plotting pipelines
write_csv(rt_real, file.path("runtime_bayes_real.csv"))

library(readr)
library(dplyr)
library(ggplot2)
library(scales)

# ---- load the CSV ----
rt <- read_csv(file.path("runtime_bayes_real.csv"),
               show_col_types = FALSE) %>%
  mutate(
    total_hr = case_when(
      !is.na(total_hour) ~ total_hour,
      !is.na(total_min)  ~ total_min/60,
      !is.na(total_sec)  ~ total_sec/3600,
      TRUE ~ NA_real_
    ),
    method  = factor(method, levels = c("BZIOLR","BALOR")),
    dataset = recode(dataset,
                     "Ravel_2011_16S_BV"               = "Ravel_2011_16S_BV \n(N=345; N_taxa=160)",
                     "HMP_2012_gingival_V35_subset"    = "HMP_2012_gingival_V35_subset \n(N=76; N_taxa=892)",
                     .default = dataset
    )
  )

# colors for datasets
ds_cols <- c("Ravel_2011_16S_BV \n(N=345; N_taxa=160)" = "#1f77b4",
             "HMP_2012_gingival_V35_subset \n(N=76; N_taxa=892)" = "#f2b01e")

# ---- Bayesian runtime panel (real data, styled) ----
p_bayes <- ggplot(rt, aes(x = total_hr, y = method, color = dataset)) +
  # points for individual runs
  geom_point(size = 2, alpha = 0.6,
             position = position_jitter(height = 0.1, width = 0)) +
  # median tick per method × dataset
  stat_summary(aes(group = interaction(method, dataset), color = dataset),
               fun = median, fun.min = median, fun.max = median,
               geom = "errorbarh", height = 0.34) +
  scale_color_manual(values = ds_cols, breaks = c("Ravel_2011_16S_BV \n(N=345; N_taxa=160)","HMP_2012_gingival_V35_subset \n(N=76; N_taxa=892)"), name = "Dataset") +
  # axis in hours
  scale_x_continuous(
    name   = "Run-time (h)",
    breaks = seq(0, 5, by = 0.5),
    limits = c(0, 5),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  labs(y = NULL) +
  theme_light(base_size = 12) +
  theme(
    axis.title.y   = element_blank(),
    axis.text.x    = element_text(size = 8),
    axis.title.x   = element_text(size = 10),
    legend.text    = element_text(size = 8),
    legend.position= "bottom",
    panel.grid.minor = element_line(color = "grey90", size = 0.25)
  ) +
  coord_cartesian(xlim = c(0, 5), clip = "off")

p_bayes
# Save figure
ggsave(
  filename = paste0("Figure_3.4B", ".png"),
  plot     = p_bayes,
  width    = 170,
  height   = 150,
  dpi      = 300,
  units    = "mm",
  bg       = "white"
)

