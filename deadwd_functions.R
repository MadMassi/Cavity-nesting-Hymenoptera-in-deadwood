###-------------------------------------------------------------------------------------
# Helper functions for data analysis of manuscript "Moisture reduces saproxylic bee, wasp, and parasitoid diversity in lying and standing deadwood"
# Author: Massimo Martini
# Date: 10th April 2026

###-------------------------------------------------------------------------------------
#Sensitivity analysis
compare_models <- function(model1, model2) {
  # Extract summary coefficient tables
  coef1 <- summary(model1)$coefficients$cond
  coef2 <- summary(model2)$coefficients$cond
  
  # Align parameter names (assumes models have the same set/order of fixed effects)
  params <- intersect(rownames(coef1), rownames(coef2)) # to cover cases where parameters differ
  coef1 <- coef1[params, , drop = FALSE]
  coef2 <- coef2[params, , drop = FALSE]
  
  # Build the comparison table
  comp <- data.frame(
    Parameter     = params,
    Estimate_1    = coef1[,"Estimate"],
    SE_1          = coef1[,"Std. Error"],
    P_1           = coef1[,"Pr(>|z|)"],
    Estimate_2    = coef2[,"Estimate"],
    SE_2          = coef2[,"Std. Error"],
    P_2           = coef2[,"Pr(>|z|)"]
  )
  
  # Add significance and direction change columns
  comp$Signif_1   <- comp$P_1 < 0.05
  comp$Signif_2   <- comp$P_2 < 0.05
  comp$SignifChange <- comp$Signif_1 != comp$Signif_2
  comp$DirChange    <- sign(comp$Estimate_1) != sign(comp$Estimate_2)
  
  return(comp)
}



#----------------------------------------------------------------------------------------------------------
# Pretty-print and export summaries for a list of glmmTMB models (v4: pseudo-R2 + fixed-width tables)
# - Likelihood-based pseudo-R2 (McFadden, Cox–Snell, Nagelkerke)
# - Fixed-width ASCII tables inside code fences for perfect column alignment
# - Robust random-effects extraction (broom.mixed preferred)
# - Pretty variable names incl. interactions with " X "
# - Prints family/link, n, logLik, AIC, pseudo-R2s, fixed effects, random effects
#
# Usage:
#   res <- pretty_glmmTMB_summaries(
#     models = mods,
#     file   = "model_summaries.md",
#     digits = 3
#   )


### Helpers for pretty_glmmTMB_summaries #######################################

add_sig_stars <- function(df) {
  pcols <- c("Pr(>|z|)", "Pr(>|t|)")
  pcol  <- intersect(pcols, names(df))
  if (length(pcol) == 0) return(df)
  p <- df[[pcol[1]]]
  df$Sig. <- ifelse(is.na(p), "",
                    ifelse(p < 0.001, "***",
                           ifelse(p < 0.01, "**",
                                  ifelse(p < 0.05, "*",
                                         ifelse(p < 0.1, ".", "")))))
  df
}

# digits-aware formatter using an option the main function will set
fmt_num <- function(x, digits = getOption("pretty_glmmTMB_digits", 3)) {
  formatC(x, format = "f", digits = digits)
}


# Fixed-width table formatters (monospaced, code fences)

format_coef_block_fixed <- function(tab, header, prettify_fun = identity) {
  if (is.null(tab) || nrow(tab) == 0) return(character(0))
  df <- as.data.frame(tab)
  
  # Stars + human-readable term names
  df$Term <- rownames(df)
  df <- add_sig_stars(df)
  df$Term <- vapply(df$Term, prettify_fun, character(1))
  
  # Which stat columns exist?
  pcol    <- intersect(c("Pr(>|z|)", "Pr(>|t|)"), names(df))
  if (length(pcol) > 1) pcol <- pcol[1]
  statcol <- intersect(c("z value", "t value"), names(df))
  
  # Final column order
  cols <- c("Term", "Estimate", "Std. Error", statcol, pcol, "Sig.")
  cols <- cols[cols %in% names(df)]
  df <- df[, cols, drop = FALSE]
  
  # Numeric formatting (p-values handled separately)
  for (nm in setdiff(names(df), pcol)) if (is.numeric(df[[nm]])) df[[nm]] <- fmt_num(df[[nm]])
  if (length(pcol) == 1 && is.numeric(df[[pcol]])) {
    df[[pcol]] <- ifelse(df[[pcol]] < 0.001, "<0.001", fmt_num(df[[pcol]]))
  }
  
  # Compute widths (min width is header width)
  widths <- vapply(names(df), function(nm) max(nchar(nm, type = "width"),
                                               nchar(df[[nm]], type = "width")), numeric(1))
  
  # Row printer with padding
  pad_row <- function(vals) paste(mapply(function(val, w) sprintf(paste0("%-", w, "s"), val), vals, widths), collapse = " | ")
  
  header_line <- pad_row(names(df))
  sep_line    <- paste(mapply(function(w) paste(rep("-", w), collapse = ""), widths), collapse = "-|-")
  rows        <- apply(df, 1, pad_row)
  
  c(paste0("### ", header),
    "",
    "```",
    header_line,
    sep_line,
    rows,
    "```",
    "",
    "Significance codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1",
    "")
}

format_ranef_fixed <- function(mod) {
  ran_df <- try({
    if (requireNamespace("broom.mixed", quietly = TRUE))
      broom.mixed::tidy(mod, effects = "ran_pars")
    else stop("no broom.mixed")
  }, silent = TRUE)
  
  if (!inherits(ran_df, "try-error") && !is.null(ran_df) && nrow(ran_df)) {
    ran_df <- subset(ran_df, effect == "ran_pars")
    ran_cond <- subset(ran_df, component %in% c("cond", NA))
    is_sd <- grepl("^sd__", ran_cond$term)
    ran_sd <- ran_cond[is_sd, , drop = FALSE]
    if (!nrow(ran_sd)) return(c("### Random effects", "", "(No random effects)", ""))
    
    clean_term <- function(x) { x <- sub("^sd__", "", x); ifelse(x == "Intercept", "(Intercept)", x) }
    out <- data.frame(
      Group    = ran_sd$group,
      Term     = vapply(ran_sd$term, clean_term, character(1)),
      `Std.Dev.` = ran_sd$estimate,
      Variance = ran_sd$estimate^2,
      check.names = FALSE
    )
    
    out$`Std.Dev.` <- fmt_num(out$`Std.Dev.`)
    out$Variance   <- fmt_num(out$Variance)
    
    widths <- vapply(names(out), function(nm) max(nchar(nm, type = "width"),
                                                  nchar(out[[nm]], type = "width")), numeric(1))
    
    pad_row <- function(vals) paste(mapply(function(val, w) sprintf(paste0("%-", w, "s"), val), vals, widths), collapse = " | ")
    
    header_line <- pad_row(names(out))
    sep_line    <- paste(mapply(function(w) paste(rep("-", w), collapse = ""), widths), collapse = "-|-")
    rows        <- apply(out, 1, pad_row)
    
    return(c("### Random effects", "", "```", header_line, sep_line, rows, "```", ""))
  }
  c("### Random effects", "", "(Random-effect details unavailable)", "")
}

# Nakagawa R2 helper (marginal / conditional)
.nakagawa_r2 <- function(model) {
  # --- 1) performance::r2_nakagawa() (most explicit + stable) ---
  if (requireNamespace("performance", quietly = TRUE)) {
    r2nk <- try(performance::r2_nakagawa(model), silent = TRUE)
    if (!inherits(r2nk, "try-error") && !is.null(r2nk)) {
      vals <- try(unlist(r2nk), silent = TRUE)
      if (!inherits(vals, "try-error")) {
        m_idx <- grep("marginal",    names(vals), ignore.case = TRUE)[1]
        c_idx <- grep("conditional", names(vals), ignore.case = TRUE)[1]
        if (!is.na(m_idx) && !is.na(c_idx)) {
          m <- suppressWarnings(as.numeric(vals[m_idx]))
          c <- suppressWarnings(as.numeric(vals[c_idx]))
          if (is.finite(m) && is.finite(c)) {
            return(c(marginal = m, conditional = c))
          }
        }
      }
    }
    
    # --- 2) fallback: performance::r2() (what you call manually) ---
    r2_any <- try(performance::r2(model), silent = TRUE)
    if (!inherits(r2_any, "try-error") && !is.null(r2_any)) {
      vals <- try(unlist(r2_any), silent = TRUE)
      if (!inherits(vals, "try-error")) {
        m_idx <- grep("marginal",    names(vals), ignore.case = TRUE)[1]
        c_idx <- grep("conditional", names(vals), ignore.case = TRUE)[1]
        if (!is.na(m_idx) && !is.na(c_idx)) {
          m <- suppressWarnings(as.numeric(vals[m_idx]))
          c <- suppressWarnings(as.numeric(vals[c_idx]))
          if (is.finite(m) && is.finite(c)) {
            return(c(marginal = m, conditional = c))
          }
        }
      }
    }
  }
  
  # --- 3) MuMIn::r.squaredGLMM() as last resort ---
  if (requireNamespace("MuMIn", quietly = TRUE)) {
    r2m <- try(MuMIn::r.squaredGLMM(model), silent = TRUE)
    if (!inherits(r2m, "try-error") && !is.null(r2m) &&
        (is.matrix(r2m) || is.data.frame(r2m)) &&
        all(c("R2m", "R2c") %in% colnames(r2m))) {
      m <- suppressWarnings(as.numeric(r2m[1, "R2m"]))
      c <- suppressWarnings(as.numeric(r2m[1, "R2c"]))
      if (is.finite(m) && is.finite(c)) {
        return(c(marginal = m, conditional = c))
      }
    }
  }
  
  # If everything fails or only non-finite values are found, return NULL
  NULL
}



# Pseudo-R2 calculator (auto-builds null unless provided)

.pseudo_r2 <- function(model, null_model = NULL) {
  if (is.null(null_model)) {
    null_try <- try(suppressWarnings(update(model, . ~ 1)), silent = TRUE)
    if (inherits(null_try, "try-error")) {
      stop("Could not auto-build a null model with update(. ~ 1). Provide null_model explicitly.")
    }
    null_model <- null_try
  }
  ll_full <- as.numeric(logLik(model))
  ll_null <- as.numeric(logLik(null_model))
  n <- tryCatch(nobs(model), error = function(e) NA_integer_)
  if (!is.finite(ll_full) || !is.finite(ll_null)) stop("Non-finite log-likelihood(s): check model convergence.")
  if (!is.finite(n) || n <= 0) stop("Could not determine sample size via nobs().")
  
  R2_McF <- 1 - (ll_full / ll_null)
  R2_CS  <- 1 - exp((2 / n) * (ll_null - ll_full))
  denom  <- 1 - exp((2 / n) * ll_null)
  R2_Nag <- if (abs(denom) < .Machine$double.eps) NA_real_ else (R2_CS / denom)
  
  list(McFadden = R2_McF, Nagelkerke = R2_Nag,
       logLik_full = ll_full, logLik_null = ll_null, n = n,
       AIC_full = tryCatch(AIC(model), error = function(e) NA_real_),
       AIC_null = tryCatch(AIC(null_model), error = function(e) NA_real_))
}


# Robustly get genpois dispersion; returns a single numeric or NA
get_genpois_dispersion <- function(mod, fam_name = NULL) {
  if (is.null(fam_name)) {
    fam_name <- tryCatch(family(mod)$family, error = function(e) "")
  }
  if (!grepl("^genpois", fam_name, ignore.case = TRUE)) return(NA_real_)
  
  # --- 1) Structured attempt from coef(summary(mod))$disp ---
  cs_try <- try(coef(summary(mod)), silent = TRUE)
  if (!inherits(cs_try, "try-error") && !is.null(cs_try$disp)) {
    disp_obj <- cs_try$disp
    # Matrix/data.frame with an "(Intercept)" row
    if (is.matrix(disp_obj) || is.data.frame(disp_obj)) {
      rn <- rownames(disp_obj)
      pick <- if (!is.null(rn) && any(rn %in% c("(Intercept)", "(Intercept).1"))) {
        which(rn %in% c("(Intercept)", "(Intercept).1"))[1]
      } else 1
      est_col <- intersect(c("Estimate", "estimate", "(Intercept)"), colnames(disp_obj))
      if (length(est_col) == 0) {
        num_cols <- colnames(disp_obj)[vapply(as.data.frame(disp_obj), is.numeric, logical(1))]
        if (length(num_cols) > 0) est_col <- num_cols[1]
      }
      if (length(est_col) == 1) {
        est <- suppressWarnings(as.numeric(disp_obj[pick, est_col]))
        if (is.finite(est)) return(exp(est))
      }
    } else if (is.numeric(disp_obj) && length(disp_obj) >= 1) {
      # Named numeric vector case
      est <- suppressWarnings(as.numeric(disp_obj[[1]]))
      if (is.finite(est)) return(exp(est))
    }
  }
  
  # --- 2) Fallback: parse printed summary line ---
  txt  <- capture.output(suppressWarnings(summary(mod)))
  line <- grep("^\\s*Dispersion parameter for genpois", txt, value = TRUE)
  if (length(line)) {
    num <- sub(".*:\\s*", "", line[1])
    val <- suppressWarnings(as.numeric(num))
    if (is.finite(val)) return(val)
  }
  
  NA_real_
}



#get autocorrelation dispersion parameters
format_ar1_disp_fixed <- function(mod) {
  # Parse the 'Dispersion model:' block from summary(glmmTMB)
  txt <- capture.output(suppressWarnings(summary(mod)))
  
  # Locate the 'Dispersion model:' section
  start <- grep("^Dispersion model:", txt)
  if (length(start) == 0) return(character(0))
  start <- start[1] + 1L  # first line after the heading
  
  # Collect lines until a stopper (Number of obs, Conditional model, Zero-inflation, blank)
  stop_idx <- length(txt)
  for (i in start:length(txt)) {
    if (grepl("^Number of obs:", txt[i]) ||
        grepl("^Conditional model:", txt[i]) ||
        grepl("^Zero-inflation model:", txt[i]) ||
        trimws(txt[i]) == "") {
      stop_idx <- i - 1L
      break
    }
  }
  
  block <- txt[start:stop_idx]
  block <- block[nzchar(trimws(block))]  # drop empty lines
  
  if (!length(block)) return(character(0))
  
  # First line is the header ("Groups Name Variance Std.Dev. Corr")
  header_line <- block[1]
  data_lines  <- block[-1]
  if (!length(data_lines)) return(character(0))
  
  # Keep only rows that clearly contain AR(1) info
  data_lines <- data_lines[grepl("\\(ar1\\)", data_lines, ignore.case = TRUE)]
  if (!length(data_lines)) return(character(0))
  
  # Parse each line by whitespace
  parsed <- strsplit(trimws(data_lines), "\\s+")
  
  tab <- do.call(rbind, lapply(parsed, function(v) {
    # Expect at least: Group, Name, Variance, Std.Dev., Corr, "(ar1)"
    # If more columns, Corr part is all remaining after the 4th numeric columns
    if (length(v) < 6) return(NULL)
    grp  <- v[1]
    nm   <- v[2]
    var  <- v[3]
    sd   <- v[4]
    corr <- paste(v[5:length(v)], collapse = " ")
    c(Group = grp, Name = nm, Variance = var, Std.Dev. = sd, Corr = corr)
  }))
  
  if (is.null(tab) || nrow(tab) == 0) return(character(0))
  
  tab <- as.data.frame(tab, stringsAsFactors = FALSE)
  
  # Numeric formatting for Variance / Std.Dev.
  tab$Variance <- fmt_num(as.numeric(tab$Variance))
  tab$Std.Dev. <- fmt_num(as.numeric(tab$Std.Dev.))
  
  # Fixed-width formatting
  widths <- vapply(names(tab), function(nm)
    max(nchar(nm, type = "width"), nchar(tab[[nm]], type = "width")),
    numeric(1)
  )
  
  pad_row <- function(vals) {
    paste(mapply(function(val, w) sprintf(paste0("%-", w, "s"), val),
                 vals, widths),
          collapse = " | ")
  }
  
  header_line <- pad_row(names(tab))
  sep_line    <- paste(mapply(function(w) paste(rep("-", w), collapse = ""),
                              widths),
                       collapse = "-|-")
  rows        <- apply(tab, 1, pad_row)
  
  c("### Dispersion / autocorrelation parameters",
    "",
    "```",
    header_line,
    sep_line,
    rows,
    "```",
    "")
}



### Pretty model printing #########################################################
# Main pretty printer

pretty_glmmTMB_summaries <- function(models,
                                     file   = "model_summaries.md",
                                     digits = 3) {
  stopifnot(is.list(models))
  if (length(models) == 0) stop("`models` is empty.")
  
  # ensure global fmt_num() uses this run's digits
  old_opts <- options(pretty_glmmTMB_digits = digits)
  on.exit(options(old_opts), add = TRUE)
  
  # mapping: original -> abbreviated (we'll reverse it below)
  pretty_map <- c(
    "logTR"         = "sc_logtr",
    "treat."        = "treatment",
    "moisutre"      = "sc_mctr",
    "ant occ."      = "ant_pres",
    "canopy cov."   = "sc_cnpy",
    "CWD"           = "sc_cwd",
    "host abund"    = "sc_loghatot",
    "host rich"     = "sc_hr",
    "enemy abund"   = "sc_pa",
    "enemy rich"    = "sc_pr",
    "geomorph."     = "geo_cat"
    
  )
  # reverse lookup (abbr -> pretty)
  abbr_to_pretty <- stats::setNames(names(pretty_map), unname(pretty_map))
  
  # term prettifier (inside so it sees abbr_to_pretty)
  prettify_term <- function(term) {
    if (term %in% c("(Intercept)", "(Intercept).1")) return(term)
    parts <- strsplit(term, ":")[[1]]
    parts <- vapply(parts, function(p) {
      if (!is.na(abbr_to_pretty[p])) return(abbr_to_pretty[p])
      pp <- p
      for (abbr in names(abbr_to_pretty)) pp <- gsub(paste0("\\b", abbr, "\\b"), abbr_to_pretty[abbr], pp)
      pp
    }, character(1))
    if (length(parts) > 1) paste(parts, collapse = " X ") else parts
  }
  
  con <- file(file, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)

  # --- Header note -----------------------------------------------------------
header_note <- c(
  "# Model results for the manuscript \"Humidity mediates saproxylic bee, wasp, and parasitoid diversity in lying and suspended deadwood\"",
  "",
  "Author: Massimo Martini",
  "Date: 10th March 2026",
  "",
  "---",
  ""
)
writeLines(header_note, con)

# --- R2 caution note --------------------------------------------------------  
  r2_note <- c(
    "## Note on R² metrics",
    "",
    "All R² and pseudo-R² values (McFadden, Nagelkerke, Nakagawa) are reported ",
    "for transparency.",
    ""
  )
  
  writeLines(r2_note, con)
  
  results <- vector("list", length(models))
  names(results) <- if (is.null(names(models))) paste0("Model_", seq_along(models)) else names(models)
  
  for (i in seq_along(models)) {
    mod <- models[[i]]
    if (is.null(mod)) next
    s <- summary(mod)
    
    mdl_name <- if (!is.null(names(models)) && nzchar(names(models)[i])) names(models)[i] else paste0("Model_", i)
    mdl_name <- paste0(LETTERS[i], ". ", mdl_name)
    
    fam <- try(family(mod), silent = TRUE)
    fam_name  <- if (!inherits(fam, "try-error") && !is.null(fam$family)) fam$family else "unknown"
    link_name <- if (!inherits(fam, "try-error") && !is.null(fam$link))   fam$link   else "unknown"
    
    n   <- try(stats::nobs(mod), silent = TRUE); if (inherits(n,   "try-error")) n   <- NA
    ll  <- try(as.numeric(stats::logLik(mod)), silent = TRUE); if (inherits(ll,  "try-error")) ll  <- NA
    aic <- try(stats::AIC(mod), silent = TRUE); if (inherits(aic, "try-error")) aic <- NA
    
    # Pseudo-R2s (build matching null automatically)
    pR2 <- .pseudo_r2(mod)
    
    # --- Genpois dispersion parameter (simple, always printed when available) ---
    extra_info <- NULL
    gp <- get_genpois_dispersion(mod, fam_name)
    if (is.finite(gp)) {
      extra_info <- paste0("**Dispersion (genpois)**: ", fmt_num(gp))
    }
    
    # Nakagawa R2 (marginal / conditional); may be NULL if not computable
    nkR2 <- .nakagawa_r2(mod)
    
    header_lines <- c(
      paste0("# ", mdl_name),
      paste0("**Family**: ", fam_name, " (link = ", link_name, ")"),
      paste0("**Observations**: ", n,
             "   |   **logLik**: ", fmt_num(ll),
             "   |   **AIC**: ", fmt_num(aic))
    )
    
    if (!is.null(extra_info)) {
      header_lines <- c(header_lines, extra_info)
    }
    
    if (!is.null(nkR2)) {
      header_lines <- c(
        header_lines,
        paste0("**Nakagawa R2 (marginal / conditional)**: ",
               fmt_num(nkR2["marginal"]), " / ", fmt_num(nkR2["conditional"]))
      )
    }
    
    header_lines <- c(
      header_lines,
      paste0("**Pseudo R2 (McFadden / Nagelkerke)**: ",
             fmt_num(pR2$McFadden), " / ", fmt_num(pR2$Nagelkerke)),
      ""
    )
    
    writeLines(header_lines, con)
    
    coef_list <- try(coef(summary(mod)), silent = TRUE)
    if (!inherits(coef_list, "try-error")) {
      if (!is.null(coef_list$cond) && nrow(coef_list$cond) > 0)
        writeLines(format_coef_block_fixed(coef_list$cond, "Fixed effects (conditional)", prettify_term), con)
      if (!is.null(coef_list$zi) && nrow(coef_list$zi) > 0)
        writeLines(format_coef_block_fixed(coef_list$zi, "Fixed effects (zero-inflation)", prettify_term), con)
      if (!is.null(coef_list$disp) && nrow(coef_list$disp) > 0)
        writeLines(format_coef_block_fixed(coef_list$disp, "Dispersion model", prettify_term), con)
    } else if (!is.null(s$coefficients)) {
      if (!is.null(s$coefficients$cond)) writeLines(format_coef_block_fixed(s$coefficients$cond, "Fixed effects (conditional)", prettify_term), con)
      if (!is.null(s$coefficients$zi))   writeLines(format_coef_block_fixed(s$coefficients$zi,   "Fixed effects (zero-inflation)", prettify_term), con)
      if (!is.null(s$coefficients$disp)) writeLines(format_coef_block_fixed(s$coefficients$disp, "Dispersion model", prettify_term), con)
    }
    
    # Random effects (conditional)
    writeLines(format_ranef_fixed(mod), con)
    
    # Dispersion / autocorrelation parameters (if present)
    ar1_block <- format_ar1_disp_fixed(mod)
    if (length(ar1_block)) writeLines(ar1_block, con)
    
    # Separator
    writeLines(c(paste(rep("=", 80), collapse = ""), ""), con)
    
    results[[i]] <- list(
      name   = mdl_name,
      family = fam_name,
      link   = link_name,
      nobs   = n,
      logLik = ll,
      AIC    = aic,
      pseudoR2   = c(
        McFadden   = pR2$McFadden,
        Nagelkerke = pR2$Nagelkerke,
        DevExplained = pR2$DevExplained  # keep this if your real .pseudo_r2 defines it
      ),
      nakagawaR2 = nkR2,
      coefs  = if (!inherits(coef_list, "try-error")) coef_list else s$coefficients,
      ranef  = try(
        if (requireNamespace("broom.mixed", quietly = TRUE))
          broom.mixed::tidy(mod, effects = "ran_pars")
        else NULL,
        silent = TRUE
      )
    )
    
  }
  
  message("Wrote summaries to: ", normalizePath(file, winslash = "/"))
  invisible(results)
}



#---------------------------------------------------------------------------------------------
#Remove path analysis sub-models from the global environment once they are already safely inside the various lists 
rm_path_objects <- function(...) {
  objs <- list(...)
  all_names <- unlist(lapply(objs, function(x) {
    nms <- names(x)
    # Keep non-empty names, skip if NA or blank
    nms[!is.na(nms) & trimws(nms) != ""]
  }))
  # Only keep names that are valid identifiers and exist in the global workspace
  valid_names <- all_names[make.names(all_names) == all_names & all_names %in% ls(envir = .GlobalEnv)]
  if (length(valid_names) > 0) rm(list = valid_names, envir = .GlobalEnv)
}




#---------------------------------------------------------------------------------------------
#calculating model dispersion parameters and creating a table for printing
get_disp_row <- function(mod, name, nsim = 500) {
  sim <- simulateResiduals(mod, n = nsim, plot = FALSE, refit = FALSE)
  td  <- suppressWarnings(testDispersion(sim, plot = FALSE))
  tibble(
    Model      = name,
    Family     = tryCatch(family(mod)$family, error = function(e) NA_character_),
    Dispersion = if ("statistic" %in% names(td)) as.numeric(td$statistic) else NA_real_,
    P_value    = if ("p.value"   %in% names(td)) as.numeric(td$p.value)    else NA_real_
  )
}




#---------------------------------------------------------------------------------------------
#setting a pretty names map
pretty_map <- c(
  "LogTR"         = "sc_sr",
  "Stand_vol."  = "sc_sv",
  "Tree_FD"       = "sc_fd",
  "Forest_age"    = "sc_fa",
  "Host_abund."    = "sc_ha",
  "Host_rich."     = "sc_hr",
  "Enemy_abund."   = "sc_pa",
  "Enemy_rich."    = "sc_pr",
  "Enemy_rich."    = "sc_pr10",
  "Slope"         = "sc_slope",
  "Elevation"     = "sc_elev",
  "Eastness"      = "sc_east",
  "Northness"     = "sc_north",
  "Host_abund." = "sc_cells"
)
#printing an anova type I table
make_type1_table <- function(mods, anova_obj, pretty_map,
                             digits = 3,
                             baseline_label = "(baseline: intercept + RE)") {
  stopifnot(is.list(mods), length(mods) >= 1)
  
  # ---- helpers --------------------------------------------------------------
  get_terms <- function(m) attr(terms(stats::formula(m)), "term.labels")
  
  # reverse lookup code -> pretty
  code_to_pretty <- function(code) {
    hits <- names(pretty_map)[match(code, unname(pretty_map))]
    ifelse(is.na(hits), code, hits)
  }
  pretty_one_term <- function(term) {
    parts <- strsplit(term, ":", fixed = TRUE)[[1]]
    # use ASCII 'x' to avoid encoding issues in CSV/Excel
    paste(vapply(parts, code_to_pretty, character(1)), collapse = " x ")
  }
  
  # Normalize/standardize ANOVA column names to a canonical set
  std_names <- function(df) {
    nn_raw <- names(df)
    nn_std <- gsub("[^A-Za-z]", "", nn_raw)  # drop spaces, symbols
    out_names <- character(length(nn_std))
    for (i in seq_along(nn_std)) {
      s <- nn_std[i]
      if (s == "Chisq") out_names[i] <- "Chisq"
      else if (s %in% c("PrChisq","Prchisq","PrChisQ","PrGTChisq")) out_names[i] <- "Pr(>Chisq)"
      else if (s %in% c("ChiDf","ChiDF","ChDf")) out_names[i] <- "Chi Df"
      else if (s %in% c("AIC")) out_names[i] <- "AIC"
      else if (s %in% c("Df","df")) out_names[i] <- "Df"
      else if (s %in% c("logLik","loglik")) out_names[i] <- "logLik"
      else if (s %in% c("deviance","Deviance")) out_names[i] <- "deviance"
      else if (s %in% c("BIC")) out_names[i] <- "BIC"
      else out_names[i] <- nn_raw[i]  # keep original if unknown
    }
    names(df) <- out_names
    df
  }
  
  # ---- align objects & names ------------------------------------------------
  rn <- names(mods)
  if (is.null(rn)) rn <- paste0("m", seq_along(mods) - 1L)
  
  a <- as.data.frame(anova_obj, stringsAsFactors = FALSE)
  a <- std_names(a)
  stopifnot(nrow(a) == length(rn))
  rownames(a) <- rn
  
  # ---- compute Added terms --------------------------------------------------
  term_list <- lapply(mods, get_terms)
  added <- character(length(mods))
  
  # --- baseline row: show random factor structure ---
  re_terms <- lme4::findbars(formula(mods[[1]]))
  if (length(re_terms) > 0) {
    re_names <- vapply(re_terms, function(x) deparse(x[[3]]), character(1))
    re_names <- gsub(":", "/", re_names)
    re_names <- gsub("_id", "", re_names, fixed = TRUE)
    re_names <- unique(re_names)
    for (r in re_names) {
      if (grepl("/", r)) {
        parent <- sub("/.*", "", r)
        re_names <- setdiff(re_names, parent)
      }
    }
    re_names <- c(re_names[grepl("/", re_names)], re_names[!grepl("/", re_names)])
    added[1] <- paste0("Random factors: ", paste(re_names, collapse = " + "))
  } else {
    added[1] <- "(no random factors)"
  }
  
  # --- detect terms added in later models ---
  for (i in 2:length(mods)) {
    new_terms <- setdiff(term_list[[i]], term_list[[i - 1]])
    added[i] <- if (length(new_terms) == 0) "—" else
      paste(vapply(new_terms, pretty_one_term, character(1)), collapse = " + ")
  }
  
  # ---- choose/ensure column order (keep only the essentials) ----------------
  desired <- c("AIC", "Chi Df", "Chisq", "Pr(>Chisq)")
  present <- intersect(desired, names(a))
  extra <- setdiff(names(a), present)  # we won't include extras
  a2 <- a[, present, drop = FALSE]
  
  # --- add significance stars column right after p-value --------------------
  if ("Pr(>Chisq)" %in% names(a2)) {
    sig <- cut(a2[["Pr(>Chisq)"]],
               breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
               labels = c("***", "**", "*", ".", ""))
    idx <- which(names(a2) == "Pr(>Chisq)")
    a2 <- cbind(a2[seq_len(idx)], Signif = sig, a2[-seq_len(idx), drop = FALSE])
  }
  
  # ---- build output ---------------------------------------------------------
  out <- cbind(
    Model = rn,
    Added = added,
    a2,
    deparse.level = 0
  )
  
  # numeric formatting
  num_cols <- setdiff(names(out), c("Model","Added","Signif"))
  
  for (cl in num_cols) {
    if (cl == "Pr(>Chisq)" && is.numeric(out[[cl]])) {
      pvals <- out[[cl]]
      out[[cl]] <- ifelse(
        pvals < 0.001,
        "<0.001",
        sprintf(paste0("%.", digits, "f"), pvals)
      )
    } else if (is.numeric(out[[cl]])) {
      out[[cl]] <- round(out[[cl]], digits)
    }
  }
  
  out
}



#---------------------------------------------------------------------------------------------
# Check linearity of GLMM predictor pearson residuals
check_linearity <- function(model, var = NULL) {

  df <- model.frame(model)

  resids <- residuals(model, type = "pearson")

  response_var <- as.character(formula(model)[[2]])
  cont_vars <- names(df)[sapply(df, is.numeric)]
  cont_vars <- setdiff(cont_vars, response_var)

  if (length(cont_vars) == 0) {
    message("No numeric predictors found.")
    return(invisible(NULL))
  }

  # If a specific variable is requested
  if (!is.null(var)) {

    if (!var %in% cont_vars) {
      stop(paste("Variable", var, "not found or not numeric in model frame."))
    }

    cont_vars <- var
  }

  n_vars <- length(cont_vars)

  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)

  # Dynamic layout
  if (n_vars == 1) {
    par(mfrow = c(1, 1), mar = c(5, 5, 3, 2))
  } else {
    par(mfrow = c(ceiling(n_vars / 2), 2), mar = c(4, 4, 2, 1))
  }

  for (v in cont_vars) {

    x <- df[[v]]

    plot(x, resids,
         xlab = v,
         ylab = "Pearson residuals",
         main = paste("Residuals vs", v),
         pch = 16,
         cex = 0.7)

    abline(h = 0, lty = 2)
    lines(stats::lowess(x, resids), lwd = 2)
  }

  invisible(NULL)
}


# ---------------------------------------------------------------------------------------------------
# Set plotting theme
theme_manuscript <- function(base_size = 28, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.title = element_text(size = base_size + 4),
      axis.text = element_text(size = base_size),
      strip.text = ggtext::element_markdown(size = base_size),
      strip.background = element_rect(fill = "grey90", color = "black", linewidth = 1),
      plot.title = element_text(size = base_size + 6, face = "bold"),
      legend.text = element_text(size = base_size),
      legend.title = element_text(size = base_size + 1, face = "plain"),
      legend.key.width = unit(1.4, "cm"),
      panel.border = element_blank(),
      axis.line    = element_line(color = "black", linewidth = 1),
      panel.spacing = unit(0.6, "lines")
    )
}

theme_set(theme_manuscript())




# ---------------------------------------------------------------------------------------------
# Swap CWD with FWD in models to compare AIC and results

make_fwd <- function(models, old = "sc_cwd", new = "sc_fwd") {
  
  has_term <- function(mod, term) {
    f <- tryCatch(formula(mod), error = function(e) NULL)
    if (is.null(f)) return(FALSE)
    grepl(paste0("\\b", term, "\\b"), paste(deparse(f), collapse = " "))
  }
  
  term_stats <- function(mod, term) {
    s <- tryCatch(summary(mod), error = function(e) NULL)
    if (is.null(s)) return(data.frame(estimate = NA_real_, p_value = NA_real_))
    
    # glmmTMB: conditional table is in $coefficients$cond
    tab <- NULL
    if (!is.null(s$coefficients) && is.list(s$coefficients) && "cond" %in% names(s$coefficients)) {
      tab <- s$coefficients$cond
    } else if (!is.null(s$coefficients)) {
      tab <- s$coefficients
    }
    
    if (is.null(tab) || is.null(rownames(tab)) || !(term %in% rownames(tab))) {
      return(data.frame(estimate = NA_real_, p_value = NA_real_))
    }
    
    est <- as.numeric(tab[term, "Estimate"])
    # p-value column name differs across model types
    pcol <- intersect(colnames(tab), c("Pr(>|z|)", "Pr(>|t|)", "Pr(>Chisq)", "p.value", "p-value"))
    pval <- if (length(pcol) >= 1) as.numeric(tab[term, pcol[1]]) else NA_real_
    
    data.frame(estimate = est, p_value = pval)
  }
  
  # pick only models containing old term
  keep_names <- names(models)[vapply(models, has_term, logical(1), term = old)]
  old_models <- models[keep_names]
  
  # refit with substitution old -> new
  models_fwd <- list()
  
  for (nm in keep_names) {
    mod <- old_models[[nm]]
    
    # update formula: remove old, add new
    mod_new <- tryCatch(
      update(mod, as.formula(paste0(". ~ . - ", old, " + ", new))),
      error = function(e) NULL
    )
    
    if (!is.null(mod_new)) {
      models_fwd[[nm]] <- mod_new
    } else {
      warning(sprintf("Could not refit model '%s' with %s -> %s (skipped).", nm, old, new))
    }
  }
  
  # Build comparison table (only for successfully refit models)
  common_names <- intersect(names(old_models), names(models_fwd))
  
  sens_table <- do.call(rbind, lapply(common_names, function(nm) {
    m_old <- old_models[[nm]]
    m_new <- models_fwd[[nm]]
    
    aic_old <- tryCatch(AIC(m_old), error = function(e) NA_real_)
    aic_new <- tryCatch(AIC(m_new), error = function(e) NA_real_)
    
    st_old <- term_stats(m_old, old)
    st_new <- term_stats(m_new, new)
    
    alpha <- 0.05  # change if you want
    
    data.frame(
      model_name       = nm,
      AIC_old          = aic_old,
      AIC_new          = aic_new,
      delta_AIC        = aic_new - aic_old,
      term_old         = old,
      estimate_old     = st_old$estimate,
      p_old            = st_old$p_value,
      term_new         = new,
      estimate_new     = st_new$estimate,
      p_new            = st_new$p_value,
      
      # did sign flip?
      direction_changed = {
        ok <- is.finite(st_old$estimate) && is.finite(st_new$estimate)
        if (!ok) NA else sign(st_old$estimate) != sign(st_new$estimate)
      },
      
      # did significance status flip?
      signif_changed = {
        ok <- is.finite(st_old$p_value) && is.finite(st_new$p_value)
        if (!ok) NA else (st_old$p_value < alpha) != (st_new$p_value < alpha)
      },
      
      stringsAsFactors = FALSE
    )
  }))
  
  rownames(sens_table) <- NULL
  
  list(models_fwd = models_fwd, sens_table = sens_table)
}









