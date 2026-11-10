##############################################################################
# fonr_erer_analysis_v3_12.R
# Ericameria seedling survival re-analysis (Fort Ord / FONR)
#
# Author: Jon Detka, PhD
# Project: Maritime chaparral revegetation – herbivory × patch context dynamics
#
# ---------------------------------------------------------------------------
# ECOLOGICAL QUESTIONS
# ---------------------------------------------------------------------------
#
# 1) How does herbivory influence seedling survival during revegetation?
#
# 2) Does vegetation patch context (inside vs outside shrub gaps)
#    modify herbivore effects?
#
# 3) Do these relationships shift between planting cohorts
#    (i.e., does herbivory intensity or vulnerability change over time)?
#
# 4) Are observed aboveground declines driven by:
#       (a) temporary top-kill with recovery, or
#       (b) true establishment failure?
#
# 5) Under herbivore exclusion (CAGE), do patch context effects persist?
#    (i.e., are inside vs outside differences present even without herbivory?)
#
# ---------------------------------------------------------------------------
# ANALYTICAL STRATEGY
# ---------------------------------------------------------------------------
#
# This script implements a dual-response survival framework to distinguish
# between transient aboveground suppression and true mortality, and includes
# standalone "CAGE-only" inference to isolate patch effects under exclusion.
#
# Two survival rules are analyzed:
#
# 1) SURVIVAL_RULE = "raw"
#    - Uses observed live/dead status at each survey.
#    - Captures top-kill, recovery, and repeated aboveground suppression.
#    - Primary inference: GEE (plant-clustered logistic regression).
#    - Purpose: Estimate marginal effects on observed survival trajectories.
#
# 2) SURVIVAL_RULE = "last_alive"
#    - Treats plants as alive until their final confirmed alive observation.
#    - Subsequent absence classified as establishment failure.
#    - Produces monotonic time-to-event data.
#    - Primary inference: Cox proportional hazards model (clustered by patch).
#    - Purpose: Estimate hazard of true mortality (viability failure).
#
# This dual structure allows separation of:
#   - Temporary herbivore damage (top-kill)
#   - Permanent establishment failure
#
# ---------------------------------------------------------------------------
# CAGE-ONLY STANDALONE INFERENCE (added v3.3)
# ---------------------------------------------------------------------------
#
# To address: "what happens under herbivore exclusion alone?",
# the script additionally fits models restricted to CAGE plants only:
#
# A) RAW (aboveground status) – CAGE-only GEE
#    - Subset: f.tmt == CAGE
#    - Model: live ~ f.type + time terms (jday_z + sin_z + cos_z)
#    - Optional cohort test: live ~ f.type * f.cohort + time terms
#    - Purpose: isolate patch context effects under exclusion.
#
# B) LAST_ALIVE (viability) – CAGE-only Cox PH
#    - Subset: f.tmt == CAGE (plant-level event table)
#    - Model: Surv(time, event) ~ f.type + cluster(patch)
#    - Optional cohort test: Surv(time, event) ~ f.type * f.cohort + cluster(patch)
#    - Purpose: evaluate establishment hazard differences without herbivory.
#
# ---------------------------------------------------------------------------
# MODEL STRUCTURE
# ---------------------------------------------------------------------------
#
# Fixed effects:
#   - f.type   : Patch context (IN vs OU)
#   - f.tmt    : Herbivory treatment (CAGE vs NONE)
#   - f.cohort : Planting cohort (1 vs 2)
#
# Core interaction:
#   f.type * f.tmt
#   → Tests whether patch context moderates herbivore effects.
#
# Cohort interaction:
#   f.type * f.tmt * f.cohort
#   → Tests whether herbivory and spatial moderation shift between cohorts.
#
# Temporal structure:
#   - jday (days since planting; cohort-relative)
#   - Seasonal sine/cosine terms (calendar seasonality)
#
# ---------------------------------------------------------------------------
# MODEL TYPES AND PURPOSE
# ---------------------------------------------------------------------------
#
# 1) GEE (geepack::geeglm)
#    - Population-averaged marginal effects
#    - Accounts for repeated measures within plant (clustered by plant)
#    - Used for aboveground status ("raw")
#
# 2) Cox proportional hazards (survival::coxph)
#    - Time-to-event model for establishment failure (viability)
#    - Cluster-robust SEs by patch
#    - Used for monotonic survival ("last_alive")
#
# 3) GLMM (lme4::glmer)
#    - Mixed-effects logistic regression
#    - Includes plant and spatial random effects
#    - Used as supporting/sensitivity analysis (may warn under separation)
#
# 4) GAM (mgcv::gam)
#    - Cohort-specific smooth trajectories
#    - Used for visualization of nonlinear survival dynamics
#
# ---------------------------------------------------------------------------
# INFERENCE LOGIC
# ---------------------------------------------------------------------------
#
# Primary inference:
#   - RAW  : GEE results for aboveground survival patterns
#   - LAST : Cox results for establishment hazard
#
# Added standalone inference:
#   - CAGE-only GEE and CAGE-only Cox isolate patch effects under exclusion
#
# Supporting inference:
#   - GLMM (scaled predictors)
#   - Fixed-site GLM
#   - Cohort-specific GAM visualizations
#
# The key biological hypothesis:
#   Herbivory increases mortality risk, and patch context moderates this effect.
#   The strength of herbivory and/or its spatial moderation may shift between cohorts.
#   CAGE-only models test whether patch effects persist even when herbivory is removed.
#
# ---------------------------------------------------------------------------
# CONTRAST TABLES (updated v3.7+: NO emmeans)
# ---------------------------------------------------------------------------
#
# We compute contrasts directly from:
#   - model coefficients, and
#   - the model variance-covariance matrix (robust where applicable).
#
# From full 3-way models we report ratios on the response scale:
#   - Herbivory effect (NONE vs CAGE) within each cohort × patch
#   - Patch effect (OU vs IN) within each cohort × treatment
#   - OVERALL effects (new in v3.10):
#       * Overall patch effect (OU vs IN) within each cohort, averaged over treatment
#       * Overall herbivory effect (NONE vs CAGE) within each cohort, averaged over patch
#
# Output scale:
#   - GEE ("raw"): odds ratios (OR)
#   - Cox ("last_alive"): hazard ratios (HR)
#
# This avoids emmeans issues with coxph link objects and keeps outputs stable.
#
# ---------------------------------------------------------------------------
# EXPECTED ECOLOGICAL INTERPRETATION
# ---------------------------------------------------------------------------
#
# Full models:
# - If f.tmtNONE > 0 in Cox → herbivory increases mortality hazard.
# - If f.typeOU > 0 in Cox → outside patches increase mortality hazard.
# - If f.typeOU:f.tmtNONE < 0 → patch context buffers herbivory.
# - Significant cohort interactions → regime shifts between planting cohorts.
#
# CAGE-only models:
# - If f.typeOU > 0 in CAGE-only Cox → outside patches increase hazard even without herbivory.
# - If f.typeOU < 0 (or OR<1 / HR<1) → inside patches confer establishment advantage under exclusion.
# - If f.typeOU:f.cohort significant (CAGE-only cohort models) → patch effect changes across cohorts
#   even when herbivory is removed.
#
# ---------------------------------------------------------------------------
# OUTPUTS
# ---------------------------------------------------------------------------
#
# For each survival rule (Output/<rule>/):
#   - Model summaries (.txt) and saved fits (.rds)
#   - Standalone CAGE-only model summaries/fits (if enabled)
#   - Contrast tables from full cohort-interaction models (if enabled; *.csv)
#   - Figures: GAM predicted trajectories + observed survival plots
#
##############################################################################

# =========================
#### 0) USER SETTINGS ####
# =========================
INCLUDE_DISTANCE <- FALSE
SURVIVAL_RULES   <- c("last_alive", "raw")  # recommended
GAM_K            <- 5
JDAY_BIN         <- 5

# Modeling toggles
RUN_GLMM_PRIMARY            <- TRUE
RUN_GLMM_NOPLANT_DIAGNOSTIC <- TRUE
RUN_GLMM_COHORT_NOPLANT     <- FALSE   # default OFF (can be slow)
RUN_FIXED_SITE_GLM          <- TRUE
RUN_GEE                     <- TRUE    # gated to SURVIVAL_RULE == "raw"
RUN_COX                     <- TRUE    # gated to SURVIVAL_RULE == "last_alive"
RUN_GAMS                    <- TRUE
RUN_FIGURES                 <- TRUE

# Standalone CAGE-only inference
RUN_CAGE_ONLY_MODELS        <- TRUE

# Contrast tables (computed directly from coef + robust vcov; no emmeans)
RUN_CONTRASTS               <- TRUE

# Optional timeout for slow fits (requires R.utils)
USE_TIMEOUT_FOR_SLOW_FITS <- TRUE
MAX_SECONDS_COHORT_GLMM   <- 10 * 60   # 10 minutes
MAX_SECONDS_GEE_PRIMARY   <- 3 * 60    # 3 minutes
MAX_SECONDS_GEE_COHORT    <- 4 * 60    # 4 minutes
MAX_SECONDS_GEE_CAGE      <- 3 * 60    # 3 minutes
MAX_SECONDS_COX_CAGE      <- 3 * 60    # 3 minutes

# ======================
#### 1) PACKAGES ####
# ======================
pkgs_core <- c(
  "dplyr","readr","lubridate","ggplot2",
  "lme4","mgcv","broom","broom.mixed","car"
)
pkgs_optional <- c("geepack","survival","clubSandwich","DHARMa","R.utils")

need <- unique(c(pkgs_core, pkgs_optional))
missing <- need[!(need %in% installed.packages()[,"Package"])]
if(length(missing)){
  message("Installing missing packages: ", paste(missing, collapse=", "))
  install.packages(missing, dependencies = TRUE)
}

invisible(lapply(pkgs_core, library, character.only = TRUE))

quiet_require <- function(pkg){
  ok <- requireNamespace(pkg, quietly = TRUE)
  if(ok) suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  ok
}
HAS_GEE    <- quiet_require("geepack")
HAS_SURV   <- quiet_require("survival")
HAS_CLUB   <- quiet_require("clubSandwich")
HAS_DHARMA <- quiet_require("DHARMa")
HAS_RUTILS <- quiet_require("R.utils")

options(stringsAsFactors = FALSE)
theme_set(theme_bw())

# ======================
#### 2) PATHS / LOGGING ####
# ======================
# EDIT THIS PATH TO YOUR PROJECT FOLDER
# setwd("C:/Users/jdetk/Documents/fonr_erer/fonr_erer_260223")

proj_root <- "."
data_dir  <- file.path(proj_root, "data")
out_root  <- file.path(proj_root, "Output")
if(!dir.exists(out_root)) dir.create(out_root, recursive = TRUE)

t0 <- Sys.time()
stamp <- function(msg) message(sprintf("[%s | +%.1f min] %s",
                                       format(Sys.time(), "%H:%M:%S"),
                                       as.numeric(difftime(Sys.time(), t0, units="mins")),
                                       msg))

# ===============================
#### 3) RUN README ####
# ===============================
write_run_readme <- function(out_dir, SURVIVAL_RULE){
  lines <- c(
    paste0("SURVIVAL_RULE: ", SURVIVAL_RULE),
    "",
    "Interpretation:",
    if (SURVIVAL_RULE == "raw") {
      "  raw -> observed aboveground status at each visit (captures top-kill + recovery)."
    } else if (SURVIVAL_RULE == "last_alive") {
      "  last_alive -> viability / establishment (alive through last confirmed alive; top-kill treated as temporary)."
    } else {
      paste0("  ", SURVIVAL_RULE)
    },
    "",
    "Primary inference choices (logic):",
    "  - last_alive: Cox time-to-event (plant-level) + clustered SEs (patch),",
    "               because it matches monotonic viability and avoids repeated-binary separation.",
    "  - raw: GEE (plant-clustered) for marginal effects + cohort-interaction GEE,",
    "         because raw includes top-kill + recovery and repeated 0/1 correlation is expected.",
    "",
    "Standalone CAGE-only inference:",
    "  - raw: CAGE-only GEE to isolate patch effects under herbivore exclusion.",
    "  - last_alive: CAGE-only Cox to isolate patch effects on establishment without herbivory.",
    "",
    "Contrasts (no emmeans):",
    "  - Mechanism contrasts from full 3-way models:",
    "      * Herbivory effect (NONE vs CAGE) within cohort × patch",
    "      * Patch effect (OU vs IN) within cohort × treatment",
    "  - Overall (treatment-averaged) ratios within each cohort:",
    "      * Overall patch effect: geometric mean of patch effects under CAGE and NONE",
    "      * Overall herbivory effect: geometric mean of herbivory effects within IN and OU",
    "",
    "Supporting models:",
    "  - Scaled GLMM (may warn under separation): used for comparison/sensitivity.",
    "  - GLMM without plant RE: diagnostic stability check (often fits faster/cleaner).",
    "  - Site-fixed GLM: checks that treatment/context patterns persist controlling for site.",
    "  - Optional cohort-interaction GLMM (NO plant RE): only if explicitly enabled.",
    "",
    paste0("Run generated: ", Sys.time())
  )
  writeLines(lines, con = file.path(out_dir, "README.txt"))
  capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
}

# ==========================================================
#### 4) HELPERS ####
# ==========================================================
coalesce_names <- function(df, base_name){
  cand <- c(paste0(base_name, ".y"), paste0(base_name, ".x"), base_name)
  present <- intersect(cand, names(df))
  if(length(present) == 0) return(df)
  combined <- df[[present[1]]]
  if(length(present) > 1){
    for(col in present[-1]) combined <- dplyr::coalesce(combined, df[[col]])
  }
  df[[base_name]] <- combined
  rm_cols <- intersect(c(paste0(base_name, ".x"), paste0(base_name, ".y")), names(df))
  if(length(rm_cols) > 0) df <- df %>% dplyr::select(-all_of(rm_cols))
  df
}

apply_survival_rule <- function(dat, rule = c("raw","first_death","last_alive")){
  rule <- match.arg(rule)
  if(rule == "raw") return(dat)

  dat <- dat %>% arrange(plant_id, date, jday)

  if(rule == "first_death"){
    dat <- dat %>%
      group_by(plant_id) %>%
      mutate(live = {
        x <- live
        first_dead <- which(x == 0)[1]
        if(!is.na(first_dead)) x[first_dead:length(x)] <- 0
        x
      }) %>% ungroup()
  }

  if(rule == "last_alive"){
    dat <- dat %>%
      group_by(plant_id) %>%
      mutate(
        max_live_jday = ifelse(any(live == 1, na.rm = TRUE),
                               max(jday[live == 1], na.rm = TRUE),
                               NA_real_),
        live = ifelse(!is.na(max_live_jday) & jday <= max_live_jday, 1, live)
      ) %>%
      ungroup() %>% select(-max_live_jday)
  }

  dat
}

fit_glmer_robust <- function(formula, data, family = binomial){
  ctrls <- list(
    lme4::glmerControl(optimizer="bobyqa",    optCtrl=list(maxfun=4e5)),
    lme4::glmerControl(optimizer="nloptwrap", optCtrl=list(maxfun=6e5))
  )
  for(ctrl in ctrls){
    fit <- tryCatch(
      lme4::glmer(formula, data=data, family=family, control=ctrl),
      error = function(e) NULL
    )
    if(!is.null(fit)) return(fit)
  }
  NULL
}

with_timeout_if_available <- function(fun, seconds){
  if(USE_TIMEOUT_FOR_SLOW_FITS && HAS_RUTILS){
    return(R.utils::withTimeout(fun(), timeout = seconds, onTimeout = "error"))
  } else {
    return(fun())
  }
}

# --- Contrast math from coef + vcov ---
# Returns a single-row data.frame with ratio, CI, and z test for L' beta
contrast_from_L <- function(beta, V, L, label, model_type, SURVIVAL_RULE){
  # Ensure names align
  b <- beta
  if(is.null(names(b))) stop("beta must be a named vector.")
  L <- L[names(b)]
  L[is.na(L)] <- 0

  log_est <- sum(L * b)
  se <- sqrt(as.numeric(t(L) %*% V %*% L))
  z  <- log_est / se
  p  <- 2 * pnorm(abs(z), lower.tail = FALSE)
  ratio <- exp(log_est)
  conf.low  <- exp(log_est - 1.96*se)
  conf.high <- exp(log_est + 1.96*se)

  data.frame(
    contrast = label,
    estimate = ratio,
    conf.low = conf.low,
    conf.high = conf.high,
    log_est = log_est,
    SE = se,
    z = z,
    p = p,
    model_type = model_type,
    SURVIVAL_RULE = SURVIVAL_RULE,
    stringsAsFactors = FALSE
  )
}

# Build L vectors for common contrasts under baseline coding:
# levels: f.type IN (ref), OU; f.tmt CAGE (ref), NONE; f.cohort 1 (ref), 2
L_patch_OU_vs_IN <- function(beta_names, cohort, tmt){
  L <- setNames(rep(0, length(beta_names)), beta_names)
  # Base OU vs IN:
  if("f.typeOU" %in% beta_names) L["f.typeOU"] <- 1
  # Add interactions as needed:
  if(tmt == "NONE" && "f.typeOU:f.tmtNONE" %in% beta_names) L["f.typeOU:f.tmtNONE"] <- 1
  if(cohort == "2" && "f.typeOU:f.cohort2" %in% beta_names) L["f.typeOU:f.cohort2"] <- 1
  if(cohort == "2" && tmt == "NONE" && "f.typeOU:f.tmtNONE:f.cohort2" %in% beta_names) L["f.typeOU:f.tmtNONE:f.cohort2"] <- 1
  L
}

L_herb_NONE_vs_CAGE <- function(beta_names, cohort, type){
  L <- setNames(rep(0, length(beta_names)), beta_names)
  # Base NONE vs CAGE:
  if("f.tmtNONE" %in% beta_names) L["f.tmtNONE"] <- 1
  # Add interactions:
  if(type == "OU" && "f.typeOU:f.tmtNONE" %in% beta_names) L["f.typeOU:f.tmtNONE"] <- 1
  if(cohort == "2" && "f.tmtNONE:f.cohort2" %in% beta_names) L["f.tmtNONE:f.cohort2"] <- 1
  if(cohort == "2" && type == "OU" && "f.typeOU:f.tmtNONE:f.cohort2" %in% beta_names) L["f.typeOU:f.tmtNONE:f.cohort2"] <- 1
  L
}

# Overall/treatment-averaged patch effect in a cohort:
# geometric mean of patch effects under CAGE and NONE (equal weights)
L_overall_patch_avg_tmt <- function(beta_names, cohort){
  Lc <- L_patch_OU_vs_IN(beta_names, cohort=cohort, tmt="CAGE")
  Ln <- L_patch_OU_vs_IN(beta_names, cohort=cohort, tmt="NONE")
  0.5*Lc + 0.5*Ln
}

# Overall/patch-averaged herbivory effect in a cohort:
# geometric mean of herbivory effects within IN and OU (equal weights)
L_overall_herb_avg_patch <- function(beta_names, cohort){
  Li <- L_herb_NONE_vs_CAGE(beta_names, cohort=cohort, type="IN")
  Lo <- L_herb_NONE_vs_CAGE(beta_names, cohort=cohort, type="OU")
  0.5*Li + 0.5*Lo
}

# ==========================
#### 5) READ + MERGE ####
# ==========================
c1_path <- file.path(data_dir, "erer_cohort_1.csv")
c2_path <- file.path(data_dir, "erer_cohort_2.csv")
if(!(file.exists(c1_path) && file.exists(c2_path))){
  stop("Cohort files missing. Put 'erer_cohort_1.csv' and 'erer_cohort_2.csv' in data/.")
}

loc  <- readr::locale(encoding = "Latin1")
d.c1 <- readr::read_csv(c1_path, show_col_types = FALSE, locale = loc)
d.c2 <- readr::read_csv(c2_path, show_col_types = FALSE, locale = loc)

for(nm in c("site","type","tmt","patch")){
  d.c1 <- coalesce_names(d.c1, nm)
  d.c2 <- coalesce_names(d.c2, nm)
}

d.base <- bind_rows(d.c1, d.c2)

keep_cols <- c("cohort","plot","site","patch","type","tmt","tag",
               "survey","live","yr","mo","da","date")
d.base <- d.base %>% select(any_of(keep_cols))

# ==========================
#### 6) TIME + SANITY ####
# ==========================
d.base <- d.base %>% mutate(date = lubridate::mdy(trimws(as.character(date))))
if(all(is.na(d.base$date))) stop("All dates NA after parsing. Check date format.")
if(!("live" %in% names(d.base))) stop("No 'live' column found.")

if(!is.numeric(d.base$live) || any(!is.na(d.base$live) & !(d.base$live %in% c(0,1)))){
  d.base <- d.base %>%
    mutate(
      live = as.character(live),
      live = trimws(live),
      live = ifelse(live %in% c("1","TRUE","T","Y","Yes","yes"), 1,
                    ifelse(live %in% c("0","FALSE","F","N","No","no"), 0, NA)),
      live = as.numeric(live)
    )
}
if(any(is.na(d.base$live))) warning("Some 'live' values are NA after coercion.")

d.base <- d.base %>%
  mutate(year = lubridate::year(date),
         doy  = lubridate::yday(date))

min_by_cohort <- d.base %>%
  group_by(cohort) %>%
  summarize(min_date = min(date, na.rm = TRUE), .groups="drop")

d.base <- d.base %>%
  left_join(min_by_cohort, by="cohort") %>%
  mutate(jday = as.numeric(difftime(date, min_date, units="days"))) %>%
  select(-min_date)

d.base <- d.base %>% mutate(live_raw = live)

# =============================================
#### 7) MAIN PIPELINE PER SURVIVAL RULE ####
# =============================================
run_pipeline <- function(d.base, SURVIVAL_RULE, out_dir){

  # initialize always
  glmm_scaled         <- NULL
  glmm_noplant        <- NULL
  glmm_cohort_noplant <- NULL
  fit_fixed_site      <- NULL

  gee_fit             <- NULL
  gee_cohort          <- NULL

  # cage-only models
  gee_cage_only        <- NULL
  gee_cage_only_cohort <- NULL
  cox_cage_only        <- NULL
  cox_cage_only_cohort <- NULL

  cox_fit             <- NULL
  cox_cohort          <- NULL
  plant_tbl           <- NULL

  gam_c1              <- NULL
  gam_c2              <- NULL

  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  write_run_readme(out_dir, SURVIVAL_RULE)

  message("\n=============================")
  message("Running SURVIVAL_RULE = ", SURVIVAL_RULE)
  message("Outputs -> ", out_dir)
  message("=============================\n")

  d.dat <- d.base %>% mutate(live = live_raw)

  d.dat <- d.dat %>%
    mutate(
      tmt  = toupper(trimws(as.character(tmt))),
      type = toupper(trimws(as.character(type))),
      tmt  = dplyr::recode(tmt, "NO C"="NONE", "NO_C"="NONE", "NO-C"="NONE",
                           "NO CAGE"="NONE", .default = tmt),
      type = dplyr::recode(type, "OUT"="OU", "OUTSIDE"="OU", .default = type),

      f.cohort = factor(as.character(cohort)),
      f.site   = factor(as.character(site)),
      f.patch  = factor(as.character(patch)),
      f.plot   = factor(as.character(plot)),
      f.tag    = factor(as.character(tag)),

      f.tmt  = factor(tmt,  levels=c("CAGE","NONE")),
      f.type = factor(type, levels=c("IN","OU")),

      f.patch_u = interaction(f.cohort, f.site, f.patch, drop=TRUE),
      f.plot_u  = interaction(f.patch_u, f.plot, drop=TRUE),
      plant_id  = interaction(f.plot_u, f.tag, drop=TRUE),
      f.plant   = factor(plant_id)
    )

  d.dat <- apply_survival_rule(d.dat, rule = SURVIVAL_RULE)

  d.obs <- d.dat %>%
    filter(!is.na(jday), jday >= 0, !is.na(live)) %>%
    mutate(jday_b = floor(jday / JDAY_BIN) * JDAY_BIN) %>%
    group_by(f.cohort, jday_b, f.type, f.tmt) %>%
    summarise(
      n  = n(),
      p  = mean(live, na.rm=TRUE),
      se = sqrt(p*(1-p)/n),
      lo = pmax(0, p - 1.96*se),
      hi = pmin(1, p + 1.96*se),
      .groups="drop"
    )

  d_in <- d.dat %>%
    filter(!is.na(jday), jday >= 0, !is.na(doy), !is.na(live)) %>%
    droplevels() %>%
    mutate(
      jday_c  = jday - mean(jday, na.rm=TRUE),
      doy_rad = 2*pi*doy/365.25,
      sin_doy = sin(doy_rad),
      cos_doy = cos(doy_rad)
    ) %>%
    mutate(
      jday_z = as.numeric(scale(jday_c)),
      sin_z  = as.numeric(scale(sin_doy)),
      cos_z  = as.numeric(scale(cos_doy))
    )

  # Distance (optional)
  cdist_center <- NULL
  cdist_scale  <- NULL
  if(INCLUDE_DISTANCE){
    if(!("cdist" %in% names(d_in))) stop("INCLUDE_DISTANCE=TRUE but 'cdist' not in data.")
    sc <- scale(d_in$cdist)
    d_in$cdist_z <- as.numeric(sc)
    cdist_center <- as.numeric(attr(sc, "scaled:center"))
    cdist_scale  <- as.numeric(attr(sc, "scaled:scale"))
    d.dat <- d.dat %>% mutate(cdist_z = (cdist - cdist_center) / cdist_scale)
  }

  print(with(d_in, table(f.cohort, f.type, f.tmt)))

  # ==========================
  # A) GLMMs (supporting)
  # ==========================
  if(RUN_GLMM_PRIMARY){
    stamp("Fitting primary GLMM (scaled predictors)...")
    f_primary <- if(INCLUDE_DISTANCE){
      live ~ f.type * f.tmt + jday_z + sin_z + cos_z + cdist_z +
        (1|f.patch_u/f.plot_u) + (1|f.plant)
    } else {
      live ~ f.type * f.tmt + jday_z + sin_z + cos_z +
        (1|f.patch_u/f.plot_u) + (1|f.plant)
    }

    glmm_scaled <- fit_glmer_robust(f_primary, data = d_in)

    if(!is.null(glmm_scaled)){
      sink(file.path(out_dir, "glmm_primary_scaled_summary.txt"))
      cat("MODEL: Primary GLMM (scaled predictors; includes plant RE)\n")
      cat("SURVIVAL_RULE:", SURVIVAL_RULE, "\n\n")
      print(summary(glmm_scaled))
      try(print(car::Anova(glmm_scaled, type=3)), silent=TRUE)
      cat("\nSingular fit? ", lme4::isSingular(glmm_scaled), "\n")
      cat("\nConvergence messages:\n")
      print(glmm_scaled@optinfo$conv$lme4$messages)
      cat("\nVarCorr:\n")
      print(VarCorr(glmm_scaled))
      sink()
      saveRDS(glmm_scaled, file.path(out_dir, "glmm_primary_scaled_fit.rds"))
    } else {
      message("Primary scaled GLMM failed to fit.")
    }
  }

  if(RUN_GLMM_NOPLANT_DIAGNOSTIC){
    stamp("Fitting GLMM without plant RE (diagnostic stability check)...")
    f_noplant <- if(INCLUDE_DISTANCE){
      live ~ f.type * f.tmt + jday_z + sin_z + cos_z + cdist_z +
        (1|f.patch_u/f.plot_u)
    } else {
      live ~ f.type * f.tmt + jday_z + sin_z + cos_z +
        (1|f.patch_u/f.plot_u)
    }
    glmm_noplant <- fit_glmer_robust(f_noplant, data=d_in)

    if(!is.null(glmm_noplant)){
      sink(file.path(out_dir, "glmm_noplant_summary.txt"))
      cat("MODEL: GLMM without (1|f.plant) (diagnostic)\n")
      cat("SURVIVAL_RULE:", SURVIVAL_RULE, "\n\n")
      print(summary(glmm_noplant))
      try(print(car::Anova(glmm_noplant, type=3)), silent=TRUE)
      cat("\nConvergence messages:\n")
      print(glmm_noplant@optinfo$conv$lme4$messages)
      cat("\nVarCorr:\n")
      print(VarCorr(glmm_noplant))
      sink()
      saveRDS(glmm_noplant, file.path(out_dir, "glmm_noplant_fit.rds"))
    } else {
      message("GLMM no-plant failed to fit.")
    }
  }

  # ==========================
  # B) Fixed-site GLM (supporting)
  # ==========================
  if(RUN_FIXED_SITE_GLM){
    stamp("Fitting fixed-site GLM (sensitivity)...")
    f_site <- if(INCLUDE_DISTANCE){
      live ~ f.type * f.tmt + jday_z + sin_z + cos_z + cdist_z + f.site
    } else {
      live ~ f.type * f.tmt + jday_z + sin_z + cos_z + f.site
    }
    fit_fixed_site <- tryCatch(glm(f_site, data=d_in, family=binomial),
                               error = function(e) NULL)
    if(!is.null(fit_fixed_site)){
      sink(file.path(out_dir, "fixed_site_summary.txt"))
      cat("MODEL: Fixed-site GLM (sensitivity)\n")
      cat("SURVIVAL_RULE:", SURVIVAL_RULE, "\n\n")
      print(summary(fit_fixed_site))
      sink()
      saveRDS(fit_fixed_site, file.path(out_dir, "fixed_site_fit.rds"))
    }
  }

  # ==========================
  # C) PRIMARY INFERENCE (GATED)
  # ==========================

  # ---- RAW: GEE ----
  if(RUN_GEE && HAS_GEE && SURVIVAL_RULE == "raw"){

    stamp("Preparing GEE inputs (plant-clustered) ...")
    d_in <- d_in %>%
      mutate(pid = as.integer(f.plant)) %>%
      filter(!is.na(pid)) %>%
      droplevels()

    cl_sizes <- as.integer(table(d_in$pid))
    stamp(sprintf("GEE clusters: %d plants; median obs/plant = %d; max obs/plant = %d",
                  length(cl_sizes), as.integer(stats::median(cl_sizes)), max(cl_sizes)))

    gee_ctrl <- geepack::geese.control(maxit = 50, epsilon = 1e-6, trace = FALSE)

    stamp("Fitting GEE (plant-clustered) ...")
    gee_fit <- tryCatch({
      with_timeout_if_available(
        fun = function(){
          geepack::geeglm(
            live ~ f.type * f.tmt + jday_z + sin_z + cos_z,
            id = pid, data = d_in, family = binomial,
            corstr = "independence",
            control = gee_ctrl
          )
        },
        seconds = MAX_SECONDS_GEE_PRIMARY
      )
    }, error = function(e){
      message("GEE fit failed/timed out: ", e$message)
      NULL
    })

    if(!is.null(gee_fit)){
      sink(file.path(out_dir, "gee_summary.txt"))
      cat("MODEL: GEE (marginal; clustered by plant)\n")
      cat("SURVIVAL_RULE:", SURVIVAL_RULE, "\n\n")
      print(summary(gee_fit))
      sink()
      saveRDS(gee_fit, file.path(out_dir, "gee_fit.rds"))
    }

    if(nlevels(d_in$f.cohort) > 1){
      stamp("Fitting cohort-interaction GEE (RAW; 3-way) ...")
      gee_cohort <- tryCatch({
        with_timeout_if_available(
          fun = function(){
            geepack::geeglm(
              live ~ f.type * f.tmt * f.cohort + jday_z + sin_z + cos_z,
              id = pid, data = d_in, family = binomial,
              corstr = "independence",
              control = gee_ctrl
            )
          },
          seconds = MAX_SECONDS_GEE_COHORT
        )
      }, error = function(e){
        message("Cohort-interaction GEE failed/timed out: ", e$message)
        NULL
      })

      if(!is.null(gee_cohort)){
        sink(file.path(out_dir, "gee_cohort3way_summary.txt"))
        cat("MODEL: GEE cohort interaction (f.type * f.tmt * f.cohort)\n")
        cat("SURVIVAL_RULE:", SURVIVAL_RULE, "\n\n")
        print(summary(gee_cohort))
        sink()
        saveRDS(gee_cohort, file.path(out_dir, "gee_cohort3way_fit.rds"))

        if(RUN_CONTRASTS){
          stamp("Computing contrasts from GEE 3-way model (coef + robust vcov) ...")
          beta <- stats::coef(gee_cohort)
          V    <- stats::vcov(gee_cohort)
          bn   <- names(beta)

          rows <- list()

          for(coh in c("1","2")){
            for(tp in c("IN","OU")){
              L <- L_herb_NONE_vs_CAGE(bn, cohort=coh, type=tp)
              rows[[length(rows)+1]] <- contrast_from_L(
                beta, V, L,
                label = paste0("Herbivory (NONE vs CAGE) | Cohort ", coh, ", ", tp),
                model_type = "GEE_OR",
                SURVIVAL_RULE = SURVIVAL_RULE
              )
            }
            for(tm in c("CAGE","NONE")){
              L <- L_patch_OU_vs_IN(bn, cohort=coh, tmt=tm)
              rows[[length(rows)+1]] <- contrast_from_L(
                beta, V, L,
                label = paste0("Patch (OU vs IN) | Cohort ", coh, ", ", tm),
                model_type = "GEE_OR",
                SURVIVAL_RULE = SURVIVAL_RULE
              )
            }
            # Overall (treatment-averaged) patch effect in cohort
            L <- L_overall_patch_avg_tmt(bn, cohort=coh)
            rows[[length(rows)+1]] <- contrast_from_L(
              beta, V, L,
              label = paste0("Overall patch (OU vs IN), treatment-avg | Cohort ", coh),
              model_type = "GEE_OR",
              SURVIVAL_RULE = SURVIVAL_RULE
            )
            # Overall (patch-averaged) herbivory effect in cohort
            L <- L_overall_herb_avg_patch(bn, cohort=coh)
            rows[[length(rows)+1]] <- contrast_from_L(
              beta, V, L,
              label = paste0("Overall herbivory (NONE vs CAGE), patch-avg | Cohort ", coh),
              model_type = "GEE_OR",
              SURVIVAL_RULE = SURVIVAL_RULE
            )
          }

          tab <- dplyr::bind_rows(rows)
          readr::write_csv(tab, file.path(out_dir, "contrast_table_raw_GEE.csv"))

          # A small "Laurel-style" wide table (just patch effects split into abiotic vs herbiv columns)
          # Here: "Abiotic" = CAGE patch effect; "Herbiv" = NONE patch effect.
          laurel_patch <- tab %>%
            dplyr::filter(grepl("^Patch \\(OU vs IN\\)", contrast)) %>%
            tidyr::extract(contrast, into=c("cohort","tmt"), regex="Cohort ([12]), (CAGE|NONE)") %>%
            dplyr::select(cohort, tmt, estimate, conf.low, conf.high, p)
          if(nrow(laurel_patch) > 0){
            laurel_wide <- laurel_patch %>%
              dplyr::mutate(col = ifelse(tmt=="CAGE","Abiotic (CAGE)","Herbiv (NONE)")) %>%
              dplyr::select(cohort, col, estimate) %>%
              tidyr::pivot_wider(names_from = col, values_from = estimate)
            readr::write_csv(laurel_wide, file.path(out_dir, "table_patch_effects_OR_wide.csv"))
          }
        }
      }
    }

    # ==========================
    # CAGE-only GEE (RAW)
    # ==========================
    if(RUN_CAGE_ONLY_MODELS){
      d_cage <- d_in %>% filter(as.character(f.tmt) == "CAGE") %>% droplevels()
      if(nrow(d_cage) > 0){
        stamp("Fitting CAGE-only GEE (RAW; patch effect under exclusion)...")
        gee_cage_only <- tryCatch({
          with_timeout_if_available(
            fun = function(){
              geepack::geeglm(
                live ~ f.type + jday_z + sin_z + cos_z,
                id = pid, data = d_cage, family = binomial,
                corstr = "independence",
                control = gee_ctrl
              )
            },
            seconds = MAX_SECONDS_GEE_CAGE
          )
        }, error = function(e){
          message("CAGE-only GEE failed/timed out: ", e$message)
          NULL
        })

        if(!is.null(gee_cage_only)){
          sink(file.path(out_dir, "gee_cage_only_summary.txt"))
          cat("MODEL: CAGE-only GEE (RAW)\n")
          cat("SURVIVAL_RULE:", SURVIVAL_RULE, "\n\n")
          print(summary(gee_cage_only))
          sink()
          saveRDS(gee_cage_only, file.path(out_dir, "gee_cage_only_fit.rds"))
        }

        if(nlevels(d_cage$f.cohort) > 1){
          stamp("Fitting CAGE-only cohort-interaction GEE (RAW; patch × cohort)...")
          gee_cage_only_cohort <- tryCatch({
            with_timeout_if_available(
              fun = function(){
                geepack::geeglm(
                  live ~ f.type * f.cohort + jday_z + sin_z + cos_z,
                  id = pid, data = d_cage, family = binomial,
                  corstr = "independence",
                  control = gee_ctrl
                )
              },
              seconds = MAX_SECONDS_GEE_CAGE
            )
          }, error = function(e){
            message("CAGE-only cohort GEE failed/timed out: ", e$message)
            NULL
          })

          if(!is.null(gee_cage_only_cohort)){
            sink(file.path(out_dir, "gee_cage_only_cohort_summary.txt"))
            cat("MODEL: CAGE-only GEE (RAW; f.type * f.cohort)\n")
            cat("SURVIVAL_RULE:", SURVIVAL_RULE, "\n\n")
            print(summary(gee_cage_only_cohort))
            sink()
            saveRDS(gee_cage_only_cohort, file.path(out_dir, "gee_cage_only_cohort_fit.rds"))
          }
        }
      } else {
        message("CAGE-only GEE skipped: no CAGE rows after filtering.")
      }
    }

  } else if(RUN_GEE && SURVIVAL_RULE != "raw"){
    stamp("Skipping GEE: only run for SURVIVAL_RULE == 'raw' (aboveground status).")
  }

  # ---- LAST_ALIVE: Cox ----
  if(RUN_COX && HAS_SURV && SURVIVAL_RULE == "last_alive"){

    plant_tbl <- d_in %>%
      arrange(f.plant, jday) %>%
      group_by(f.plant, f.cohort, f.type, f.tmt, f.patch_u, f.site) %>%
      summarize(
        last_obs_j  = max(jday, na.rm = TRUE),
        last_live_j = if(any(live == 1, na.rm = TRUE)) max(jday[live == 1], na.rm = TRUE) else NA_real_,
        first_dead_after_last_live = {
          if(is.na(last_live_j)) NA_real_
          else {
            dd <- jday[live == 0 & jday > last_live_j]
            if(length(dd) == 0) NA_real_ else min(dd, na.rm = TRUE)
          }
        },
        .groups="drop"
      ) %>%
      mutate(
        time  = ifelse(is.na(first_dead_after_last_live), last_obs_j, first_dead_after_last_live),
        event = ifelse(is.na(first_dead_after_last_live), 0, 1)
      )

    write.csv(plant_tbl, file.path(out_dir, "plant_event_table_last_alive.csv"), row.names=FALSE)

    stamp("Fitting Cox model (plant-level time-to-event; clustered by patch)...")
    cox_fit <- tryCatch(
      survival::coxph(
        survival::Surv(time, event) ~ f.type * f.tmt + cluster(f.patch_u),
        data = plant_tbl
      ),
      error = function(e) NULL
    )

    if(!is.null(cox_fit)){
      sink(file.path(out_dir, "cox_summary.txt"))
      cat("MODEL: Cox PH (viability; last_alive), clustered by patch\n")
      cat("SURVIVAL_RULE:", SURVIVAL_RULE, "\n\n")
      print(summary(cox_fit))
      sink()
      saveRDS(cox_fit, file.path(out_dir, "cox_fit.rds"))
    }

    if(!is.null(cox_fit) && nlevels(plant_tbl$f.cohort) > 1){
      stamp("Fitting cohort-interaction Cox (viability; 3-way) ...")
      cox_cohort <- tryCatch(
        survival::coxph(
          survival::Surv(time, event) ~ f.type * f.tmt * f.cohort + cluster(f.patch_u),
          data = plant_tbl
        ),
        error = function(e) NULL
      )
      if(!is.null(cox_cohort)){
        sink(file.path(out_dir, "cox_cohort3way_summary.txt"))
        cat("MODEL: Cox cohort interaction (f.type * f.tmt * f.cohort), clustered by patch\n")
        cat("SURVIVAL_RULE:", SURVIVAL_RULE, "\n\n")
        print(summary(cox_cohort))
        sink()
        saveRDS(cox_cohort, file.path(out_dir, "cox_cohort3way_fit.rds"))

        if(RUN_CONTRASTS){
          stamp("Computing contrasts from Cox 3-way model (coef + robust vcov) ...")
          beta <- stats::coef(cox_cohort)
          V    <- stats::vcov(cox_cohort)
          bn   <- names(beta)

          rows <- list()

          for(coh in c("1","2")){
            for(tp in c("IN","OU")){
              L <- L_herb_NONE_vs_CAGE(bn, cohort=coh, type=tp)
              rows[[length(rows)+1]] <- contrast_from_L(
                beta, V, L,
                label = paste0("Herbivory (NONE vs CAGE) | Cohort ", coh, ", ", tp),
                model_type = "Cox_HR",
                SURVIVAL_RULE = SURVIVAL_RULE
              )
            }
            for(tm in c("CAGE","NONE")){
              L <- L_patch_OU_vs_IN(bn, cohort=coh, tmt=tm)
              rows[[length(rows)+1]] <- contrast_from_L(
                beta, V, L,
                label = paste0("Patch (OU vs IN) | Cohort ", coh, ", ", tm),
                model_type = "Cox_HR",
                SURVIVAL_RULE = SURVIVAL_RULE
              )
            }
            # Overall (treatment-averaged) patch effect in cohort
            L <- L_overall_patch_avg_tmt(bn, cohort=coh)
            rows[[length(rows)+1]] <- contrast_from_L(
              beta, V, L,
              label = paste0("Overall patch (OU vs IN), treatment-avg | Cohort ", coh),
              model_type = "Cox_HR",
              SURVIVAL_RULE = SURVIVAL_RULE
            )
            # Overall (patch-averaged) herbivory effect in cohort
            L <- L_overall_herb_avg_patch(bn, cohort=coh)
            rows[[length(rows)+1]] <- contrast_from_L(
              beta, V, L,
              label = paste0("Overall herbivory (NONE vs CAGE), patch-avg | Cohort ", coh),
              model_type = "Cox_HR",
              SURVIVAL_RULE = SURVIVAL_RULE
            )
          }

          tab <- dplyr::bind_rows(rows)
          readr::write_csv(tab, file.path(out_dir, "contrast_table_last_alive_Cox.csv"))

          # Laurel-style wide table: HR patch effects under CAGE vs NONE by cohort
          laurel_patch <- tab %>%
            dplyr::filter(grepl("^Patch \\(OU vs IN\\)", contrast)) %>%
            tidyr::extract(contrast, into=c("cohort","tmt"), regex="Cohort ([12]), (CAGE|NONE)") %>%
            dplyr::select(cohort, tmt, estimate, conf.low, conf.high, p)
          if(nrow(laurel_patch) > 0){
            laurel_wide <- laurel_patch %>%
              dplyr::mutate(col = ifelse(tmt=="CAGE","Abiotic (CAGE)","Herbiv (NONE)")) %>%
              dplyr::select(cohort, col, estimate) %>%
              tidyr::pivot_wider(names_from = col, values_from = estimate)
            readr::write_csv(laurel_wide, file.path(out_dir, "table_patch_effects_HR_wide.csv"))
          }
        }
      }
    }

    # ==========================
    # CAGE-only Cox (LAST_ALIVE)
    # ==========================
    if(RUN_CAGE_ONLY_MODELS){
      plant_tbl_cage <- plant_tbl %>% filter(as.character(f.tmt) == "CAGE") %>% droplevels()

      if(nrow(plant_tbl_cage) > 0){
        stamp(sprintf("CAGE-only Cox inputs: %d plants; events=%d",
                      nrow(plant_tbl_cage), sum(plant_tbl_cage$event, na.rm=TRUE)))

        stamp("Fitting CAGE-only Cox (viability; patch effect under exclusion)...")
        cox_cage_only <- tryCatch({
          with_timeout_if_available(
            fun = function(){
              survival::coxph(
                survival::Surv(time, event) ~ f.type + cluster(f.patch_u),
                data = plant_tbl_cage
              )
            },
            seconds = MAX_SECONDS_COX_CAGE
          )
        }, error = function(e){
          message("CAGE-only Cox failed/timed out: ", e$message)
          NULL
        })

        if(!is.null(cox_cage_only)){
          sink(file.path(out_dir, "cox_cage_only_summary.txt"))
          cat("MODEL: CAGE-only Cox PH (viability; clustered by patch)\n")
          cat("SURVIVAL_RULE:", SURVIVAL_RULE, "\n\n")
          print(summary(cox_cage_only))
          sink()
          saveRDS(cox_cage_only, file.path(out_dir, "cox_cage_only_fit.rds"))
        }

        if(nlevels(plant_tbl_cage$f.cohort) > 1){
          stamp("Fitting CAGE-only Cox cohort interaction (viability; patch × cohort)...")
          cox_cage_only_cohort <- tryCatch({
            with_timeout_if_available(
              fun = function(){
                survival::coxph(
                  survival::Surv(time, event) ~ f.type * f.cohort + cluster(f.patch_u),
                  data = plant_tbl_cage
                )
              },
              seconds = MAX_SECONDS_COX_CAGE
            )
          }, error = function(e){
            message("CAGE-only Cox cohort failed/timed out: ", e$message)
            NULL
          })

          if(!is.null(cox_cage_only_cohort)){
            sink(file.path(out_dir, "cox_cage_only_cohort_summary.txt"))
            cat("MODEL: CAGE-only Cox PH (viability; f.type * f.cohort)\n")
            cat("SURVIVAL_RULE:", SURVIVAL_RULE, "\n\n")
            print(summary(cox_cage_only_cohort))
            sink()
            saveRDS(cox_cage_only_cohort, file.path(out_dir, "cox_cage_only_cohort_fit.rds"))
          }
        }
      } else {
        message("CAGE-only Cox skipped: no CAGE plants in plant_tbl.")
      }
    }

  } else if(RUN_COX && SURVIVAL_RULE != "last_alive"){
    stamp("Skipping Cox: only run for SURVIVAL_RULE == 'last_alive' (viability).")
  }

  # ==========================
  # D) GAMs (trajectory visualization)
  # ==========================
  if(RUN_GAMS){
    run_cohort_gam <- function(dat, cohort_id){
      datc <- dat %>%
        filter(as.character(cohort) == as.character(cohort_id),
               jday >= 0, !is.na(live)) %>%
        droplevels()

      if(nrow(datc) == 0) return(NULL)

      f_gam <- if(INCLUDE_DISTANCE){
        as.formula(paste0("live ~ f.type * f.tmt + cdist_z + s(jday, k=", GAM_K, ") + s(f.patch_u, bs='re')"))
      } else {
        as.formula(paste0("live ~ f.type * f.tmt + s(jday, k=", GAM_K, ") + s(f.patch_u, bs='re')"))
      }

      fit <- tryCatch(
        mgcv::gam(f_gam, data = datc, family = binomial(link="logit"), method = "REML"),
        error = function(e) NULL
      )

      if(!is.null(fit)){
        saveRDS(fit, file.path(out_dir, paste0("gam_cohort", cohort_id, ".rds")))
        sink(file.path(out_dir, paste0("gam_cohort", cohort_id, "_summary.txt")))
        cat("MODEL: Cohort-specific GAM\n")
        cat("SURVIVAL_RULE:", SURVIVAL_RULE, "\n\n")
        print(summary(fit))
        print(anova(fit))
        sink()
      }
      fit
    }

    gam_c1 <- run_cohort_gam(d.dat, 1)
    gam_c2 <- run_cohort_gam(d.dat, 2)
  }

  # ==========================
  # E) FIGURES
  # ==========================
  if(RUN_FIGURES){
    # Manuscript-ready labels
    cohort_labs <- c(
      "1" = "Drought-year cohort",
      "2" = "Higher-rainfall cohort"
    )
    
    patch_labs <- c(
      "IN" = "Dieback Patch",
      "OU" = "Legacy Patch"
    )
    
    tmt_labs <- c(
      "CAGE" = "Caged",
      "NONE" = "Uncaged"
    )
    
    theme_manuscript <- theme_bw(base_size = 12) +
      theme(
        plot.title = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )

    predict_gam_marginal <- function(gam_obj, dat_ref, cohort_id, jday_by=5){

      grid <- expand.grid(
        jday   = seq(0, max(dat_ref$jday, na.rm=TRUE), by=jday_by),
        f.type = levels(dat_ref$f.type),
        f.tmt  = levels(dat_ref$f.tmt)
      )

      levs <- levels(gam_obj$model$f.patch_u)
      grid$f.patch_u <- factor(levs[1], levels = levs)

      if(INCLUDE_DISTANCE && "cdist_z" %in% names(dat_ref)){
        grid$cdist_z <- mean(dat_ref$cdist_z, na.rm=TRUE)
      }

      pr <- predict(gam_obj, newdata=grid, type="link", se.fit=TRUE, exclude=c("s(f.patch_u)"))
      grid$eta   <- pr$fit
      grid$se    <- pr$se.fit
      grid$mu    <- plogis(grid$eta)
      grid$lower <- plogis(grid$eta - 1.96*grid$se)
      grid$upper <- plogis(grid$eta + 1.96*grid$se)
      grid$cohort <- factor(cohort_id)
      grid
    }

    if(!is.null(gam_c1) && !is.null(gam_c2)){
      d_c1 <- d.dat %>% filter(as.character(cohort)=="1", jday>=0, !is.na(live))
      d_c2 <- d.dat %>% filter(as.character(cohort)=="2", jday>=0, !is.na(live))

      g1 <- predict_gam_marginal(gam_c1, d_c1, 1, jday_by=5)
      g2 <- predict_gam_marginal(gam_c2, d_c2, 2, jday_by=5)
      g_all <- dplyr::bind_rows(g1, g2)

      p1 <- ggplot(g_all, aes(x = jday, y = mu, color = f.type, linetype = f.tmt)) +
        geom_ribbon(aes(ymin = lower, ymax = upper, fill = f.type),
                    alpha = 0.15, color = NA) +
        geom_line(linewidth = 1) +
        facet_wrap(
          ~cohort,
          ncol = 1,
          labeller = labeller(cohort = cohort_labs)
        ) +
        scale_color_discrete(
          name = "Patch type",
          breaks = names(patch_labs),
          labels = unname(patch_labs)
        ) +
        scale_fill_discrete(
          name = "Patch type",
          breaks = names(patch_labs),
          labels = unname(patch_labs)
        ) +
        scale_linetype_discrete(
          name = "Herbivore treatment",
          breaks = names(tmt_labs),
          labels = unname(tmt_labs)
        ) +
        coord_cartesian(ylim = c(0, 1)) +
        labs(
          x = "Days since planting",
          y = "Predicted probability alive"
        ) +
        theme_manuscript
      
      ggsave(
        file.path(out_dir, "Fig_S1_GAM_smoothed_survival_curves.png"),
        p1, width = 7.5, height = 8, dpi = 300
      )
      
      ggsave(
        file.path(out_dir, "Fig_S1_GAM_smoothed_survival_curves.pdf"),
        p1, width = 7.5, height = 8
      )
    }

    d_fig2 <- d.dat %>%
      filter(!is.na(jday), jday>=0, !is.na(live)) %>%
      mutate(
        cohort = factor(as.character(cohort)),
        jday_b = floor(jday / JDAY_BIN) * JDAY_BIN
      ) %>%
      group_by(cohort, jday_b, f.type, f.tmt) %>%
      summarise(
        n=n(),
        p=mean(live, na.rm=TRUE),
        se=sqrt(p*(1-p)/n),
        lo=pmax(0, p - 1.96*se),
        hi=pmin(1, p + 1.96*se),
        .groups="drop"
      )

    p2 <- ggplot(d_fig2, aes(x = jday_b, y = p, color = f.type, linetype = f.tmt, fill = f.type)) +
      geom_ribbon(aes(ymin = lo, ymax = hi),
                  alpha = 0.15, color = NA) +
      geom_line(linewidth = 1) +
      facet_wrap(
        ~cohort,
        ncol = 1,
        labeller = labeller(cohort = cohort_labs)
      ) +
      scale_color_discrete(
        name = "Patch type",
        breaks = names(patch_labs),
        labels = unname(patch_labs)
      ) +
      scale_fill_discrete(
        name = "Patch type",
        breaks = names(patch_labs),
        labels = unname(patch_labs)
      ) +
      scale_linetype_discrete(
        name = "Herbivore treatment",
        breaks = names(tmt_labs),
        labels = unname(tmt_labs)
      ) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(
        x = "Days since planting",
        y = "Observed proportion alive"
      ) +
      theme_manuscript
    
    ggsave(
      file.path(out_dir, "Fig_1_observed_survival_curves.png"),
      p2, width = 8, height = 9, dpi = 300
    )
    
    ggsave(
      file.path(out_dir, "Fig_1_observed_survival_curves.pdf"),
      p2, width = 8, height = 9
    )

    d_caged <- d.obs %>% filter(as.character(f.tmt)=="CAGE")
    p3 <- ggplot(d_caged, aes(x=jday_b, y=p, color=f.type, fill=f.type)) +
      geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.15, color=NA) +
      geom_line(linewidth=1) +
      facet_wrap(~ factor(f.cohort), ncol=1) +
      coord_cartesian(ylim=c(0,1)) +
      labs(
        title = paste0("Observed survival under herbivore exclusion (CAGE only): IN vs OU (", SURVIVAL_RULE, ")"),
        x="Days since planting (jday)", y="Observed proportion alive"
      )
    ggsave(file.path(out_dir, "Figure3_caged_only_patch_contrast.png"), p3, width=8, height=8, dpi=300)

    d_fig4 <- d.dat %>%
      filter(!is.na(date), !is.na(live)) %>%
      mutate(cohort = factor(as.character(cohort))) %>%
      group_by(cohort, date, f.type, f.tmt) %>%
      summarise(
        n=n(),
        p=mean(live, na.rm=TRUE),
        se=sqrt(p*(1-p)/n),
        lo=pmax(0, p - 1.96*se),
        hi=pmin(1, p + 1.96*se),
        .groups="drop"
      ) %>%
      mutate(cohort_label = ifelse(as.character(cohort)=="1",
                                   "Cohort 1 (Apr 2021 – Jan 2022)",
                                   "Cohort 2 (Jan 2022 – Feb 2023)"))

    cohort_starts <- data.frame(start_date = as.Date(c("2021-04-07","2022-01-04")))

    p4 <- ggplot(d_fig4, aes(x=date, y=p, color=f.type, linetype=f.tmt, fill=f.type)) +
      geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.15, color=NA) +
      geom_line(linewidth=1) +
      geom_vline(data=cohort_starts, aes(xintercept=start_date),
                 linetype="dotted", linewidth=0.6, alpha=0.8) +
      facet_wrap(~cohort_label, ncol=1) +
      scale_x_date(date_breaks="3 months", date_labels="%b %Y") +
      coord_cartesian(ylim=c(0,1)) +
      labs(
        title = paste0("Observed proportions by calendar date (", SURVIVAL_RULE, ")"),
        x="Calendar date", y="Observed proportion alive"
      )
    ggsave(file.path(out_dir, "Figure4_observed_survival_by_date.png"), p4, width=8, height=9, dpi=300)
  }

  # ==========================
  # F) Optional DHARMa for GLMM
  # ==========================
  if(HAS_DHARMA && !is.null(glmm_scaled)){
    simres <- DHARMa::simulateResiduals(glmm_scaled, n=200)
    png(file.path(out_dir, "glmm_DHARMa_residuals.png"), width=900, height=700)
    plot(simres)
    dev.off()
    saveRDS(simres, file.path(out_dir, "glmm_DHARMa_residuals.rds"))
  }

  saveRDS(list(d.dat=d.dat, d_in=d_in),
          file.path(out_dir, "data_processed.rds"))

  message("Done for SURVIVAL_RULE = ", SURVIVAL_RULE)

  invisible(list(
    gee = gee_fit,
    gee_cohort3way = gee_cohort,
    cox = cox_fit,
    cox_cohort3way = cox_cohort,
    plant_tbl = plant_tbl
  ))
}

# =============================================
#### 8) RUN ALL REQUESTED RULES ####
# =============================================
results <- list()
for(rule in SURVIVAL_RULES){
  out_dir <- file.path(out_root, rule)
  results[[rule]] <- run_pipeline(d.base, SURVIVAL_RULE = rule, out_dir = out_dir)
}

saveRDS(results, file.path(out_root, "results_all_rules.rds"))
message("All done. See Output/<rule>/ for summaries, models, figures, and contrast tables.")

# =============================================
#### 9) CONSOLIDATED KEY RATIOS TABLE ####
# =============================================
# Output: Output/KEY_table_ratios_all.csv
key_files <- c(
  file.path(out_root, "raw", "contrast_table_raw_GEE.csv"),
  file.path(out_root, "last_alive", "contrast_table_last_alive_Cox.csv")
)
key_files <- key_files[file.exists(key_files)]

if(length(key_files) > 0){
  key_tabs <- lapply(key_files, function(fp){
    x <- tryCatch(readr::read_csv(fp, show_col_types = FALSE), error = function(e) NULL)
    if(!is.null(x)) x$source_file <- basename(fp)
    x
  })
  key_tabs <- key_tabs[!sapply(key_tabs, is.null)]
  if(length(key_tabs) > 0){
    key_all <- dplyr::bind_rows(key_tabs)
    readr::write_csv(key_all, file.path(out_root, "KEY_table_ratios_all.csv"))
  }
}
# =============================================
#### 10) COHORT 1 POST-CAGE-REMOVAL FOLLOW-UP FIGURE ####
# =============================================

library(tidyverse)
library(lubridate)

# Read files
# Adjust names if needed to match your actual files
c1_main <- read_csv("data/erer_cohort_1.csv")
c1_all  <- read_csv("data/erer_cohort_1_all_surveys.csv")

# Use survey 0 from cleaned cohort 1 as the plot-level truth table
plot_lookup <- c1_main %>%
  filter(survey == 0) %>%
  distinct(plot, site, patch, type, tmt) %>%
  rename(
    type_clean = type,
    tmt_clean = tmt
  )

# Join clean plot/treatment labels onto all-surveys file
c1_followup <- c1_all %>%
  select(-any_of(c("type", "tmt"))) %>%
  left_join(plot_lookup, by = c("plot", "site", "patch")) %>%
  mutate(
    date = mdy(date),
    patch_label = case_when(
      type_clean == "IN" ~ "Dieback gap",
      type_clean == "OU" ~ "Legacy patch",
      TRUE ~ as.character(type_clean)
    ),
    treatment_history = case_when(
      tmt_clean == "CAGE" ~ "Originally caged",
      tmt_clean == "NONE" ~ "Uncaged",
      TRUE ~ as.character(tmt_clean)
    ),
    group = paste(patch_label, treatment_history, sep = " / ")
  )

# Keep only post-cage-removal/follow-up period
c1_post <- c1_followup %>%
  filter(date > as.Date("2022-01-03"))

# Observed survival proportions
post_surv_summary <- c1_post %>%
  group_by(date, survey, patch_label, treatment_history, group) %>%
  summarize(
    n = n(),
    alive = sum(live, na.rm = TRUE),
    prop_alive = alive / n,
    .groups = "drop"
  )

# Plot
fig_c1_post <- ggplot(
  post_surv_summary,
  aes(
    x = date,
    y = prop_alive,
    color = patch_label,
    linetype = treatment_history,
    group = group
  )
) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format()
  ) +
  labs(
    x = "Date",
    y = "Observed survival",
    color = "Patch type",
    linetype = "Original treatment",
    title = "Cohort 1 post-removal follow-up survival",
    subtitle = "Curves begin after the active cage-exclusion treatment ended"
  ) +
  theme_classic()

fig_c1_post

# Save figure
ggsave(
  filename = file.path(out_root, "figures", "FigureS_post_removal_cohort1_survival.png"),
  plot = fig_c1_post,
  width = 7,
  height = 5,
  dpi = 300
)
