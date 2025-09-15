library(tidyverse)
library(stringr)

ALPHA <- 0.10
Z90   <- qnorm(1 - ALPHA/2)

# ---------- helpers ----------
safe_load_env <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  e
}

normalize_cols <- function(df) {
  df <- as_tibble(df)
  nm <- names(df)
  
  # sanity checks
  if (!"taxon" %in% nm)  stop("Data frame lacks 'taxon'. Got: ", paste(nm, collapse=", "))
  if (!"method" %in% nm) stop("Data frame lacks 'method'. Got: ", paste(nm, collapse=", "))
  
  # estimate: prefer 'estimate'; else from 'est' or 'beta'; else NA
  if (!"estimate" %in% nm) {
    if ("est" %in% nm) {
      df <- df %>% mutate(estimate = .data$est) %>% select(-est)
    } else if ("beta" %in% nm) {
      df <- df %>% mutate(estimate = .data$beta)
    } else {
      df$estimate <- NA_real_
    }
  }
  
  # standard columns
  if (!"se" %in% names(df)) df$se <- NA_real_
  if (!"p"  %in% names(df)) df$p  <- NA_real_
  if (!"q"  %in% names(df)) df$q  <- NA_real_
  
  # lwr/upr from estimate ± Z90*se when available
  if (!"lwr" %in% names(df)) {
    df$lwr <- ifelse(!is.na(df$estimate) & !is.na(df$se),
                     df$estimate - Z90 * df$se, NA_real_)
  }
  if (!"upr" %in% names(df)) {
    df$upr <- ifelse(!is.na(df$estimate) & !is.na(df$se),
                     df$estimate + Z90 * df$se, NA_real_)
  }
  
  # significant: prefer q; else CI sign if available; else NA
  if (!"significant" %in% names(df)) {
    df$significant <- if (!all(is.na(df$q))) {
      df$q < ALPHA
    } else if (!all(is.na(df$lwr)) && !all(is.na(df$upr))) {
      (df$lwr > 0) | (df$upr < 0)
    } else {
      NA
    }
  }
  
  # types + order
  df$taxon  <- as.character(df$taxon)
  df$method <- as.character(df$method)
  
  df %>% relocate(taxon, method, estimate, se, lwr, upr, p, q, significant)
}

# ---------- find files (ABS-LOAD-SHIFT tags) ----------
base_files  <- Sys.glob("absLOADshift_N*_M*_seed*.RData") %>% sort()
extra_files <- Sys.glob("absLOADshift_EXTRAONLY_SUP_N*_M*_seed*.RData") %>% sort()

stopifnot(length(base_files)  > 0)
stopifnot(length(extra_files) > 0)

seed_of <- function(x) as.integer(str_match(x, "seed(\\d{3})\\.RData$")[,2])
base_map  <- tibble(file = base_files,  seed = seed_of(base_files))  %>% filter(!is.na(seed))
extra_map <- tibble(file = extra_files, seed = seed_of(extra_files)) %>% filter(!is.na(seed))

seeds_to_merge <- intersect(base_map$seed, extra_map$seed) %>% sort()
message("Found ", length(seeds_to_merge), " seeds in common: ",
        paste(seeds_to_merge, collapse = ", "))

# ---------- merge loop ----------
for (sd in seeds_to_merge) {
  bf <- base_map  %>% filter(seed == sd) %>% pull(file) %>% .[[1]]
  xf <- extra_map %>% filter(seed == sd) %>% pull(file) %>% .[[1]]
  message(">>> Merging seed ", sd, "\n    baseline: ", bf, "\n    extras:   ", xf)
  
  e_base  <- safe_load_env(bf)
  e_extra <- safe_load_env(xf)
  
  if (!exists("res_all2", envir = e_base))
    stop("Baseline file lacks res_all2: ", bf)
  if (!exists("res_extras_only", envir = e_extra))
    stop("Extras file lacks res_extras_only: ", xf)
  
  base_df   <- as_tibble(e_base$res_all2)          %>% mutate(seed = e_base$SEED_MAIN, source = "baseline")
  extras_df <- as_tibble(e_extra$res_extras_only)  %>% mutate(seed = e_extra$SEED_MAIN, source = "extras")
  
  base_std   <- normalize_cols(base_df)
  extras_std <- normalize_cols(extras_df)
  res_merged <- bind_rows(base_std, extras_std)
  
  # carry useful objects from baseline
  SEED_MAIN <- e_base$SEED_MAIN
  n         <- if (exists("n",        envir = e_base)) e_base$n        else NA_integer_
  counts    <- if (exists("counts",   envir = e_base)) e_base$counts   else NULL
  counts2   <- if (exists("counts2",  envir = e_base)) e_base$counts2  else NULL
  groups    <- if (exists("groups",   envir = e_base)) e_base$groups   else NULL
  true_effects <- if (exists("true_effects", envir = e_base)) e_base$true_effects else NULL
  M_base    <- if (!is.null(counts2)) ncol(counts2) else NA_integer_
  
  RUN_TAG <- sprintf("absLOADshift_ADD_ALL_N%d_M%d_seed%03d", n, M_base, SEED_MAIN)
  
  # save both RData and CSV
  save(SEED_MAIN, n, counts, counts2, groups, true_effects,
       res_merged,
       file = paste0(RUN_TAG, ".RData"))
  readr::write_csv(res_merged, paste0(RUN_TAG, ".csv"))
  
  cat(sprintf("Saved: %s.[RData|csv]  rows=%d  methods=%s\n",
              RUN_TAG, nrow(res_merged),
              paste(sort(unique(res_merged$method)), collapse = ", ")))
}

message("Done.")
