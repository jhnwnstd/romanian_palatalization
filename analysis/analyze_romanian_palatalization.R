#!/usr/bin/Rscript

# Romanian Palatalization Analysis
#
# Run styler::style_file("analysis/analyze_romanian_palatalization.R") after edits

# =========================================================================
# Setup
# =========================================================================

required_pkgs <- c(
  "dplyr", "readr", "stringr", "tidyr", "broom",
  "brms", "posterior", "loo", "brglm2", "lme4"
)

# Keep this block: script should be runnable on a clean machine
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  cat("Installing missing packages:", paste(missing_pkgs, collapse = ", "), "\n")
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

# For cmdstanr setup details, see:
# Fruehwald, Josef. "Getting `Brms` and Stan Set Up."
# https://lin611-2024.github.io/notes/side-notes/content/stan.html

# cmdstanr requires special installation from r-universe
cmdstan_repo <- "https://stan-dev.r-universe.dev"

# Keep this: we want a single script that can bootstrap cmdstanr if needed
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  cat("Installing cmdstanr from r-universe...\n")
  install.packages("cmdstanr", repos = c(cmdstan_repo, getOption("repos")))
}

# If CmdStan is missing, install once and reuse across runs
if (requireNamespace("cmdstanr", quietly = TRUE)) {
  if (!dir.exists(cmdstanr::cmdstan_path())) {
    cat("Installing CmdStan backend (this may take a few minutes)...\n")
    cat("For troubleshooting, see:\n")
    cat("  - https://mc-stan.org/cmdstanr/articles/cmdstanr.html\n")
    cat("  - https://lin611-2024.github.io/notes/side-notes/content/stan.html\n")
    cmdstanr::install_cmdstan()
  }
}

suppressPackageStartupMessages(lapply(c(required_pkgs, "cmdstanr"), library, character.only = TRUE))
options(dplyr.summarise.inform = FALSE)

# =========================================================================
# Constants
# =========================================================================

segments_of_interest <- c("c", "g", "t", "d", "s", "z") # consonants where palatalization is tracked
front_verb_suffixes <- c("-i", "-iza", "-ifica") # front-vowel verbalizers; all start with /i/. -ui is back (u blocks palatalization).
suffix_interest <- c("-ic", "-ist", "-esc", "-ică", "-ice") # denominal/adjectival suffixes we track

# NDEB = non-derived exception base lemmas (gimpe / ochi / paduche patterns)
ndeb_classes <- c("gimpe", "ochi", "paduchi")
ndeb_observable <- c("ochi", "paduchi") # only these are observable as DE in the plural domain

plural_opportunities <- c("i", "e") # core i/e opportunity domain
plural_opportunities_all <- c("i", "e", "uri", "none") # all plural types in raw data

# =========================================================================
# Analysis Constants
# =========================================================================

MIN_SAMPLE_SIZE_BAYESIAN <- 5L # avoid Bayesian TP on tiny cells
SMALL_CELL_THRESHOLD <- 20L # flag small downsampled cells; TP binary unreliable below this

PRECISION_THETA <- 2L # print precision for theta_N
PRECISION_PROB <- 4L # print precision for probabilities / rates

CLUSTER_TYPES <- c("st", "sc", "ct") # cluster types for TP cluster analysis

# =========================================================================
# Run-mode toggles
# =========================================================================

RUN_BAYESIAN_TP <- FALSE # set TRUE for Bayesian TP fits (slow)
RUN_SEGMENT_CLASS_BRMS <- FALSE # set TRUE for segment-class brms (slow)
RUN_DERIV_MIXED <- FALSE # set TRUE to add GLMM (1|stem_final)

# =========================================================================
# Helper Functions
# =========================================================================

cat_section <- function(title) {
  cat("\n", title, "\n", sep = "")
}

print_full <- function(x) {
  print(x, n = Inf, width = Inf)
}

calc_rate <- function(df) {
  df |>
    mutate(rate_mut = if_else(.data$N_opp > 0, .data$N_mut / .data$N_opp, NA_real_))
}

# Tolerance Principle: θ_N = N / ln N, undefined for N ≤ 1
tp_threshold <- function(N) {
  ifelse(N > 1, N / log(N), NA_real_)
}

tp_table <- function(df, type_col, mutated_col, non_mutated_col) {
  df |>
    transmute(
      type        = {{ type_col }},
      mutated     = {{ mutated_col }},
      non_mutated = {{ non_mutated_col }}
    ) |>
    mutate(
      N = .data$mutated + .data$non_mutated,
      rate = if_else(.data$N > 0, .data$mutated / .data$N, NA_real_),
      majority_mutates = case_when(
        .data$N == 0L ~ NA,
        .data$mutated == .data$non_mutated ~ NA,
        TRUE ~ .data$mutated > .data$non_mutated
      ),
      # Exceptions are ALWAYS non-mutated items; the question is always
      # "is mutation productive?", so e = non_mutated.
      exceptions = if_else(.data$N > 0L, as.numeric(.data$non_mutated), NA_real_),
      theta_N = tp_threshold(.data$N),
      tolerated = if_else(
        is.na(.data$exceptions) | is.na(.data$theta_N),
        NA,
        .data$exceptions <= .data$theta_N
      ),
      # Flag small cells where the TP binary is unreliable; see Bayesian TP for these.
      small_n = !is.na(.data$N) & .data$N < SMALL_CELL_THRESHOLD
    )
}

fit_brms_bernoulli <- function(formula,
                               data,
                               seed,
                               intercept_sd = 2.5,
                               b_sd = 2.5,
                               include_b_prior = TRUE) {
  intercept_prior_str <- sprintf("normal(0, %s)", format(intercept_sd, scientific = FALSE, trim = TRUE))
  priors <- brms::set_prior(intercept_prior_str, class = "Intercept")

  if (include_b_prior) {
    b_prior_str <- sprintf("normal(0, %s)", format(b_sd, scientific = FALSE, trim = TRUE))
    priors <- c(priors, brms::set_prior(b_prior_str, class = "b"))
  }

  brm(
    formula,
    data = data,
    family = bernoulli(link = "logit"),
    prior = priors,
    chains = 4, iter = 2000, warmup = 1000, cores = 4,
    backend = "cmdstanr",
    seed = seed,
    refresh = 0,
    silent = 2
  )
}

# Detect palatalization in derived forms from ORTHOGRAPHY.
# More reliable than IPA-based detection because:
#   - No false positives from word-initial palatals (ciocoi)
#   - No false positives from base-internal ʃ (cinste)
#   - Works even when IPA is missing
# Romanian orthography transparently encodes palatalization:
#   c/g before i/e (without h) = palatalized
#   ți = t palatalized, și = s palatalized, zi = d palatalized
detect_palatal_from_orth <- function(lemma, verb, stem_final) {
  verb <- tolower(verb)
  lemma_l <- tolower(lemma)
  tail5 <- str_sub(verb, -5L)
  case_when(
    is.na(verb) | verb == "" ~ NA,
    # Dorsals: ci/ce = palatalized, chi/che/ca/co/cu = not
    stem_final == "c" & str_detect(tail5, "(?<!h)c[ie]") ~ TRUE,
    stem_final == "c" ~ FALSE,
    stem_final == "g" & str_detect(tail5, "(?<!h)g[ie]") ~ TRUE,
    stem_final == "g" ~ FALSE,
    # Coronals
    stem_final == "t" & str_detect(tail5, "\u021bi|\u021be") ~ TRUE,
    stem_final == "t" ~ FALSE,
    stem_final == "d" & str_ends(verb, "zi") &
      !str_ends(lemma_l, "z") ~ TRUE,
    stem_final == "d" ~ FALSE,
    stem_final == "s" & str_detect(tail5, "\u0219i") ~ TRUE,
    stem_final == "s" ~ FALSE,
    stem_final == "z" & str_detect(tail5, "ji") ~ TRUE,
    stem_final == "z" ~ FALSE,
    TRUE ~ NA
  )
}

compute_segment_tp_tables <- function(data, label_suffix = "", group_var = stem_final) {
  group_var_name <- rlang::as_name(rlang::enquo(group_var))

  by_opp <- data |>
    group_by({{ group_var }}, opportunity) |>
    summarise(
      mutated     = sum(mutation, na.rm = TRUE),
      non_mutated = sum(!mutation, na.rm = TRUE),
      .groups     = "drop"
    ) |>
    mutate(type = paste0("<", .data[[group_var_name]], "> + <-", opportunity, "> plural", label_suffix)) |>
    tp_table(type, mutated, non_mutated)

  combined <- data |>
    group_by({{ group_var }}) |>
    summarise(
      mutated     = sum(mutation, na.rm = TRUE),
      non_mutated = sum(!mutation, na.rm = TRUE),
      .groups     = "drop"
    ) |>
    mutate(type = paste0("<", .data[[group_var_name]], "> + <-i, -e> plural", label_suffix)) |>
    tp_table(type, mutated, non_mutated)

  bind_rows(by_opp, combined) |> arrange(type)
}

run_bayesian_tp <- function(data, subset_label, seed_value = 123L) {
  if (is.null(data) || nrow(data) == 0L) {
    return(list())
  }

  results <- list()

  for (seg in segments_of_interest) {
    for (opp in plural_opportunities) {
      seg_data <- data |>
        filter(stem_final == seg, opportunity == opp, !is.na(mutation)) |>
        mutate(mutation = as.integer(mutation))

      if (nrow(seg_data) < MIN_SAMPLE_SIZE_BAYESIAN) next

      N <- nrow(seg_data)
      theta_N <- tp_threshold(N)
      tolerance_rate <- theta_N / N

      sum_mut <- sum(seg_data$mutation, na.rm = TRUE)
      sum_non <- N - sum_mut
      cat(sprintf(
        "  Fitting <%s> + <-%s> (%s): N=%d, mutated=%d, non-mutated=%d\n",
        seg, opp, subset_label, N, sum_mut, sum_non
      ))

      invisible(capture.output(
        {
          model_tp <- fit_brms_bernoulli(
            mutation ~ 1, seg_data,
            seed = seed_value,
            intercept_sd = 1.5,
            b_sd = 1.5,
            include_b_prior = FALSE
          )
        },
        type = "output"
      ))

      post_samples <- as_draws_df(model_tp)
      p_mutate <- plogis(post_samples$b_Intercept)
      p_exception <- 1 - p_mutate
      prob_intolerable <- mean(p_exception > tolerance_rate)

      results[[paste(seg, opp, subset_label, sep = "_")]] <- data.frame(
        segment          = seg,
        opportunity      = opp,
        subset           = subset_label,
        majority_mutates = sum_mut > sum_non,
        N                = N,
        theta_N          = round(theta_N, PRECISION_THETA),
        tolerance_rate   = round(tolerance_rate, PRECISION_PROB),
        median_p_mutate  = round(median(p_mutate), PRECISION_PROB),
        ci_lower         = round(quantile(p_mutate, 0.025), PRECISION_PROB),
        ci_upper         = round(quantile(p_mutate, 0.975), PRECISION_PROB),
        prob_exceeds_tol = round(prob_intolerable, PRECISION_PROB)
      )
    }
  }
  results
}

segment_class_factor <- function(stem_final) {
  factor(if_else(stem_final %in% c("c", "g"), "dorsal", "coronal"),
    levels = c("dorsal", "coronal")
  )
}

# Inflection vs derivation: Firth-penalized logistic regression (handles complete
# separation), with and without frequency covariate.  McNemar is retained as a
# descriptive count of discordant pairs only — it is severely underpowered here.
#
# Expects df to contain columns: base_col, deriv_col, and (optionally) lemma_freq.
analyze_inf_vs_deriv <- function(df, base_col, deriv_col, label) {
  if (is.null(df) || nrow(df) == 0) {
    cat("\n(No data for ", label, ")\n\n", sep = "")
    return(invisible(NULL))
  }

  df <- df |>
    mutate(
      base_mut  = .data[[base_col]],
      deriv_mut = .data[[deriv_col]]
    ) |>
    filter(!is.na(base_mut), !is.na(deriv_mut))

  if (nrow(df) == 0) {
    cat("\n(No valid mutation data for ", label, ")\n\n", sep = "")
    return(invisible(NULL))
  }

  # Agreement pattern counts
  cat("\nInflection vs derivation agreement patterns (", label, "):\n", sep = "")
  patterns <- df |>
    mutate(pattern = case_when(
      base_mut & deriv_mut ~ "both_mutate",
      !base_mut & !deriv_mut ~ "both_nonmutate",
      base_mut & !deriv_mut ~ "inflection_only",
      !base_mut & deriv_mut ~ "derivation_only",
      TRUE ~ NA_character_
    )) |>
    count(pattern, sort = TRUE) |>
    mutate(prop = n / sum(n))
  print_full(patterns)

  # McNemar — descriptive only; underpowered for the sample sizes here
  tab <- with(df, table(base_mut, deriv_mut))
  if (all(dim(tab) == c(2L, 2L))) {
    cat("\nMcNemar test (", label, ", descriptive only — see Firth regression for inference):\n", sep = "")
    print_full(broom::tidy(stats::mcnemar.test(tab)))
    n_disc <- tab[1, 2] + tab[2, 1]
    if (n_disc < 20L) {
      cat("  NOTE: only", n_disc, "discordant pairs — McNemar has insufficient power.\n")
    }
  }

  # Firth-penalized logistic regression: base mutation as predictor of deriv mutation.
  # brglmFit guarantees finite estimates even under complete or quasi-complete separation.
  cat("\nFirth logistic (brglm2): deriv_mut ~ base_mut (", label, ")\n", sep = "")
  model_base <- glm(deriv_mut ~ base_mut,
    data = df,
    family = binomial(),
    method = brglm2::brglmFit,
    control = brglm2::brglm_control(maxit = 500)
  )
  tidy_base <- broom::tidy(model_base)
  print_full(tidy_base)
  # Report OR directly
  or_base <- exp(coef(model_base)["base_mutTRUE"])
  or_ci <- exp(confint(model_base)["base_mutTRUE", ])
  cat(sprintf(
    "  OR(base_mutTRUE) = %.3f [95%% CI: %.3f, %.3f]\n",
    or_base, or_ci[1], or_ci[2]
  ))

  # Frequency-controlled model: does base_mut predict deriv_mut after accounting
  # for lemma frequency?  If base_mut remains significant, the effect is not
  # a frequency-of-lexicalization artefact.
  has_freq <- "lemma_freq" %in% names(df) &&
    sum(!is.na(df$lemma_freq) & df$lemma_freq > 0) >= 5L

  if (has_freq) {
    cat("\nFirth logistic (brglm2): deriv_mut ~ base_mut + log(lemma_freq) (",
      label, ")\n",
      sep = ""
    )
    df_freq <- df |> mutate(log_freq = log(lemma_freq + 1))
    model_freq <- glm(deriv_mut ~ base_mut + log_freq,
      data = df_freq,
      family = binomial(),
      method = brglm2::brglmFit,
      control = brglm2::brglm_control(maxit = 500)
    )
    print_full(broom::tidy(model_freq))
    or_freq <- exp(coef(model_freq)["base_mutTRUE"])
    or_freq_ci <- exp(confint(model_freq)["base_mutTRUE", ])
    cat(sprintf(
      "  OR(base_mutTRUE | freq controlled) = %.3f [95%% CI: %.3f, %.3f]\n",
      or_freq, or_freq_ci[1], or_freq_ci[2]
    ))
    freq_coef <- broom::tidy(model_freq) |> filter(term == "log_freq")
    if (nrow(freq_coef) > 0 && freq_coef$p.value > 0.05) {
      cat("  => log_freq is NOT significant: base_mut effect is not driven by frequency.\n")
    } else if (nrow(freq_coef) > 0) {
      cat("  => log_freq IS significant: frequency partially predicts derivational mutation.\n")
    }
  }

  # Optional GLMM with stem_final random effect (requires lme4 and RUN_DERIV_MIXED)
  if (exists("RUN_DERIV_MIXED") && RUN_DERIV_MIXED &&
    requireNamespace("lme4", quietly = TRUE) &&
    n_distinct(df$stem_final) >= 3L) {
    cat("\nGLMM (lme4): deriv_mut ~ base_mut + log(lemma_freq) + (1|stem_final) (",
      label, ")\n",
      sep = ""
    )
    df_mixed <- df |> mutate(log_freq = log(coalesce(lemma_freq, 0) + 1))
    model_mixed <- lme4::glmer(deriv_mut ~ base_mut + log_freq + (1 | stem_final),
      data   = df_mixed,
      family = binomial()
    )
    print(summary(model_mixed))
  }

  invisible(df)
}

make_deriv_summary <- function(df, deriv_col, label_yes, label_no) {
  df |>
    mutate(base_plural_mutates = mutation_inflect) |>
    group_by(base_plural_mutates) |>
    summarise(
      mutated     = sum(.data[[deriv_col]], na.rm = TRUE),
      non_mutated = sum(!.data[[deriv_col]], na.rm = TRUE),
      .groups     = "drop"
    ) |>
    mutate(type = if_else(base_plural_mutates, label_yes, label_no)) |>
    tp_table(type, mutated, non_mutated) |>
    select(type, N, mutated, non_mutated, rate, majority = majority_mutates, tolerated)
}

# Build reference (top-1000) downsampled lexicon.
build_downsampled_lexicon <- function(nouns_opp, n_lex = 1000L) {
  freq_summary <- nouns_opp |>
    group_by(lemma) |>
    summarise(
      lemma_freq = max(freq_ron_news_2024_1M, na.rm = TRUE),
      .groups    = "drop"
    ) |>
    mutate(lemma_freq = if_else(is.na(lemma_freq) | !is.finite(lemma_freq), 0, lemma_freq)) |>
    filter(lemma_freq > 0)

  if (nrow(freq_summary) == 0L) {
    return(NULL)
  }

  top_lemmas <- freq_summary |>
    arrange(desc(lemma_freq)) |>
    slice_head(n = min(n_lex, nrow(freq_summary))) |>
    pull(lemma)

  reference <- nouns_opp |>
    left_join(freq_summary, by = "lemma") |>
    filter(lemma %in% top_lemmas) |>
    arrange(desc(lemma_freq), lemma, plural) |>
    distinct(lemma, .keep_all = TRUE)

  cat("Downsampled lexicon (top", length(top_lemmas), "most frequent lemmas):\n")
  cat("  unique lemmas:", n_distinct(reference$lemma), "\n")
  cat("  rows:", nrow(reference), "\n\n")

  reference
}

# =========================================================================
# Data Input
# =========================================================================

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "data/romanian_lexicon_with_freq.csv"
output_log <- file.path("analysis", "romanian_palatalization_analysis.log")

cat("Working directory:", getwd(), "\n")
cat("Log file:", normalizePath(output_log, mustWork = FALSE), "\n")

# Capture output to the log file (and echo to terminal via split)
log_con <- file(output_log, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit(
  {
    sink(type = "message")
    sink()
    close(log_con)
  },
  add = TRUE
)

options(
  width = 200,
  tibble.print_max = Inf,
  tibble.print_min = Inf,
  tibble.width = Inf,
  dplyr.print_max = Inf,
  pillar.max_footer_lines = Inf
)

cat("READING DATA\n")
cat("File:", input_file, "\n\n")

lex <- suppressWarnings(
  read_csv(input_file,
    show_col_types = FALSE,
    col_types = cols(
      derived_adj     = col_character(),
      ipa_derived_adj = col_character(),
      ipa_raw_lemma   = col_character(),
      ipa_raw_pl      = col_character()
    )
  )
)

parse_problems <- problems(lex)
if (nrow(parse_problems) > 0) {
  cat("\nWARNING: Found", nrow(parse_problems), "parsing issue(s):\n")
  print(parse_problems, n = min(10, nrow(parse_problems)))
  if (nrow(parse_problems) > 10) {
    cat("... and", nrow(parse_problems) - 10, "more\n")
  }
  cat("\n")
}

lex <- lex |>
  mutate(freq_ron_news_2024_1M = as.numeric(freq_ron_news_2024_1M)) |>
  mutate(
    across(c(pos, gender, stem_final, cluster, plural), ~ trimws(as.character(.))),
    pos = toupper(pos),
    opportunity = if_else(opportunity %in% plural_opportunities_all, opportunity, NA_character_),
    nde_class = case_when(
      str_starts(exception_reason, "nde:") ~ str_remove(exception_reason, "^nde:"),
      TRUE ~ "none"
    ),
    target_is_suffix = target_is_suffix == "True",
    lemma_suffix = if_else(is.na(lemma_suffix) | lemma_suffix == "", "none", lemma_suffix),
    plural_ending = case_when(
      is.na(plural) | plural == "" ~ "none",
      str_ends(plural, "uri") ~ "uri",
      str_ends(plural, "i") ~ "i",
      str_ends(plural, "e") ~ "e",
      TRUE ~ "other"
    ),
    opportunity_tp = case_when(
      nde_class %in% ndeb_classes & opportunity == "none" & plural_ending %in% c("i", "e") ~ plural_ending,
      nde_class %in% ndeb_classes & opportunity == "none" & plural_ending == "uri" ~ "i",
      TRUE ~ opportunity
    )
  )

nouns <- filter(lex, pos == "N")

# =============================================================================
# TP Domain Datasets
# =============================================================================

nouns_tp <- nouns |>
  filter(tp_in_domain == TRUE) |>
  mutate(
    cluster_simple = case_when(
      startsWith(cluster, "st") ~ "st",
      startsWith(cluster, "sc") ~ "sc",
      startsWith(cluster, "ct") ~ "ct",
      TRUE ~ NA_character_
    )
  )

nouns_opp <- nouns |>
  filter(
    opportunity_tp %in% plural_opportunities, !is.na(stem_final),
    stem_final %in% segments_of_interest
  ) |>
  mutate(
    opportunity = opportunity_tp,
    cluster_simple = case_when(
      startsWith(cluster, "st") ~ "st",
      startsWith(cluster, "sc") ~ "sc",
      startsWith(cluster, "ct") ~ "ct",
      TRUE ~ NA_character_
    ),
    exception_category = case_when(
      mutation ~ "undergoes",
      nde_class %in% ndeb_classes ~ paste0("NDEB_", nde_class),
      lemma_suffix %in% suffix_interest ~ paste0("suffix_", lemma_suffix),
      exception_reason == "unexplained" ~ "true_exception",
      TRUE ~ "other_non_undergoer"
    )
  )

# Decompose the diff between nouns_tp and manual NDEB filter to confirm
# what tp_in_domain excludes beyond NDEB.
nouns_opp_no_ndeb_check <- nouns_opp |> filter(!(nde_class %in% ndeb_classes))
tp_vs_opp_diff <- nrow(nouns_tp) - nrow(nouns_opp_no_ndeb_check)
cat("TP domain vs manual NDEB filter diff:", tp_vs_opp_diff, "\n")
if (tp_vs_opp_diff != 0L) {
  cat("Items excluded from TP beyond NDEB (by reason and suffix):\n")
  nouns_opp_no_ndeb_check |>
    filter(!tp_in_domain) |>
    count(exception_reason, lemma_suffix, sort = TRUE) |>
    print_full()
}
rm(nouns_opp_no_ndeb_check)

ndeb_rows <- filter(nouns, nde_class %in% ndeb_classes)

cat_section("BASIC COUNTS")
cat("Total rows:", nrow(lex), "\n")
cat("Nouns:", nrow(nouns), "\n")
cat("Nouns in PRODUCTIVE TP domain (tp_in_domain=TRUE):", nrow(nouns_tp), "\n")
cat("  (excludes NDEB, underlying palatals, suffix-internal targets)\n")
cat("Nouns in i/e domain with target segments (incl. NDEB):", nrow(nouns_opp), "\n")
cat("TP vs manual NDEB diff:", tp_vs_opp_diff, "\n\n")

# =========================================================================
# Quality Control
# =========================================================================

cat_section("QC: GENDER ON NOUN ROWS")
nouns_missing_gender <- filter(nouns, is.na(gender) | gender == "")
cat("Missing gender:", nrow(nouns_missing_gender), "\n")
if (nrow(nouns_missing_gender) > 0) {
  nouns_missing_gender |>
    select(lemma, gloss, gender, source, notes) |>
    head(10) |>
    print()
}

cat_section("QC: MUTATION VS OPPORTUNITY")
bad_mut_opp <- filter(nouns, mutation, !(opportunity %in% plural_opportunities))
cat("Mutation outside i/e opportunity:", nrow(bad_mut_opp), "\n")
if (nrow(bad_mut_opp) > 0) {
  bad_mut_opp |>
    select(lemma, plural, stem_final, opportunity, mutation) |>
    head(10) |>
    print()
}

cat_section("QC: NDEB ITEMS AND OPPORTUNITY")
cat("Total NDEB nouns:", nrow(ndeb_rows), "\n")
count(ndeb_rows, nde_class) |> print()

ndeb_ochi_pad <- filter(ndeb_rows, nde_class %in% ndeb_observable)
ndeb_gimpe <- filter(ndeb_rows, nde_class == "gimpe")
cat("\nNDEB nouns of ochi/păduche type (observable DE exceptions):", nrow(ndeb_ochi_pad), "\n")
if (nrow(ndeb_ochi_pad) > 0) {
  ndeb_ochi_pad |>
    select(lemma, plural, stem_final, nde_class, opportunity) |>
    arrange(nde_class, stem_final, lemma) |>
    print_full()
}
cat("\nNDEB nouns of gimpe type (unobservable as DE; excluded from exception counts):", nrow(ndeb_gimpe), "\n")

cat_section("QC: OPPORTUNITY VS PLURAL ENDING")
nouns |>
  count(opportunity, plural_ending) |>
  arrange(opportunity, plural_ending) |>
  print()

cat_section("QC: MISMATCHES BETWEEN OPPORTUNITY AND PLURAL")
inconsistent <- nouns |>
  filter(
    (opportunity == "i" & !str_ends(plural, "i")) |
      (opportunity == "e" & !str_ends(plural, "e")) |
      (opportunity == "uri" & !str_ends(plural, "uri"))
  )
cat("Number of mismatches:", nrow(inconsistent), "\n")
if (nrow(inconsistent) > 0) {
  inconsistent |>
    select(lemma, plural, opportunity, plural_ending) |>
    head(10) |>
    print()
}

cat_section("QC: SUFFIX ANNOTATIONS")
suffix_rows <- filter(nouns, lemma_suffix %in% suffix_interest)
cat("Nouns with tracked suffixes:", nrow(suffix_rows), "\n")
count(suffix_rows, lemma_suffix) |> print()

cat_section("QC: PALATAL CONSONANT IN PLURAL")
palatal_nouns <- filter(nouns, !is.na(palatal_consonant_pl), palatal_consonant_pl != "")
cat("Nouns with palatal_consonant_pl populated:", nrow(palatal_nouns), "\n")
count(palatal_nouns, stem_final, palatal_consonant_pl, sort = TRUE) |> print_full()

palatal_no_mutation <- filter(nouns, !is.na(palatal_consonant_pl), palatal_consonant_pl != "", !mutation)
cat("\nNouns with palatal_consonant_pl but mutation=FALSE:", nrow(palatal_no_mutation), "\n")

mutation_no_palatal <- filter(
  nouns, mutation, opportunity %in% plural_opportunities,
  is.na(palatal_consonant_pl) | palatal_consonant_pl == ""
)
cat("Nouns with mutation=TRUE but missing palatal_consonant_pl:", nrow(mutation_no_palatal), "\n")

cat_section("QC: DUPLICATE LEMMAS")
lemma_dups <- count(lex, lemma, pos, sort = TRUE) |> filter(n > 1)
cat("Duplicate (lemma, pos) pairs:", nrow(lemma_dups), "\n")
if (nrow(lemma_dups) > 0) head(lemma_dups, 25) |> print()

cat_section("QC: SUFFIX EXCLUSION FROM TP DOMAIN")
suffix_excluded <- nouns |>
  filter(target_is_suffix, lemma_suffix %in% suffix_interest) |>
  count(lemma_suffix, tp_in_domain = tp_in_domain == TRUE)
cat("Suffix-internal items by TP domain membership:\n")
print_full(suffix_excluded)
# Machine-verifiable: all suffix-internal tracked-suffix items must be outside TP domain
stopifnot(all(suffix_excluded$tp_in_domain == FALSE))
cat("CONFIRMED: no suffix-internal targets remain in TP domain\n")

# =========================================================================
# Inflection vs Derivation
# =========================================================================

cat_section("INFLECTION VS DERIVATION: LEMMA-BASED PATTERNS")

has_verb_deriv_cols <- all(c("derived_verbs", "ipa_derived_verbs", "deriv_suffixes") %in% names(lex))
has_adj_deriv_cols <- "derived_adj" %in% names(lex)

if (!has_verb_deriv_cols && !has_adj_deriv_cols) {
  cat("No derivational columns present; skipping all inflection/derivation checks.\n")
} else {
  # Build base dataset: ALL nouns with target consonants, not just
  # i/e domain.  This lets us test the full Steriade claim: nouns
  # selecting -uri (non-palatalizing) vs -i/-e (palatalizing) plurals
  # should differ in their verbalizer allomorph selection.
  noun_base_inflect <- nouns |>
    filter(
      !is.na(stem_final),
      stem_final %in% segments_of_interest,
      !is.na(mutation)
    ) |>
    mutate(
      mutation_inflect = as.logical(mutation),
      lemma_freq = freq_ron_news_2024_1M,
      # Surface-allomorph class. Used as a stratification variable for
      # the i/e domain regression below. Distinct from mutation_inflect:
      # for coronals, opportunity=="e" but mutation==FALSE in most cases.
      plural_class = case_when(
        opportunity %in% plural_opportunities ~ "front",
        opportunity == "uri" ~ "back",
        TRUE ~ "ambig"
      ),
      # Steriade's actual K/TS vs K/K split: does the noun's plural
      # palatalize? Use this — NOT plural_class — for the inflection-
      # dependence Fisher test. plural_class miscodes coronal -e
      # plurals as "front/alternating" via the allomorph proxy.
      mutates_class = case_when(
        as.logical(mutation) ~ "alt",
        !as.logical(mutation) & opportunity %in% c("i", "e", "uri") ~ "nonalt",
        TRUE ~ "ambig"
      )
    )

  # -----------------------------------------------------------------------
  # (1) Noun → Verb
  # -----------------------------------------------------------------------
  cat_section("(1) NOUN LEMMAS: INFLECTIONAL PLURALS VS DENOMINAL VERBS")

  if (has_verb_deriv_cols) {
    denom_pairs <- noun_base_inflect |>
      filter(
        !is.na(derived_verbs), derived_verbs != "",
        # NDE items have lemma-internal palatals; their derived
        # verbs inherit those palatals rather than producing them
        # by rule, so they cannot be used as evidence about
        # derivational palatalization productivity.
        !str_starts(coalesce(exception_reason, ""), "nde:")
      ) |>
      mutate(
        verb_suffix_front = deriv_suffixes %in% front_verb_suffixes,
        mutation_deriv_verb = detect_palatal_from_orth(lemma, derived_verbs, stem_final)
      ) |>
      filter(!is.na(mutation_deriv_verb)) |>
      arrange(lemma) |>
      distinct(lemma, .keep_all = TRUE)

    denom_pairs <- denom_pairs |>
      mutate(seg_class = if_else(stem_final %in% c("c", "g"), "dorsal", "coronal"))

    cat("Denominal N-V lemmas with target consonants:", nrow(denom_pairs), "\n")
    cat("  front plural (i/e):", sum(denom_pairs$plural_class == "front"), "\n")
    cat("  back plural (uri):", sum(denom_pairs$plural_class == "back"), "\n")
    cat("  ambig (none):", sum(denom_pairs$plural_class == "ambig"), "\n")
    cat(
      "  dorsal:", sum(denom_pairs$seg_class == "dorsal"),
      "  coronal:", sum(denom_pairs$seg_class == "coronal"), "\n"
    )
    cat(
      "  with front-vowel suffix -i/-ui:",
      sum(denom_pairs$verb_suffix_front, na.rm = TRUE), "\n"
    )

    if (nrow(denom_pairs) > 0) {
      # Analysis 1: mutation_inflect as predictor (original)
      # Restricted to i/e domain for comparability with TP tables
      denom_ie <- filter(denom_pairs, plural_class == "front")
      if (nrow(denom_ie) > 0) {
        cat("\n--- i/e domain subset (mutation_inflect as predictor) ---\n")
        analyze_inf_vs_deriv(denom_ie, "mutation_inflect", "mutation_deriv_verb", "N -> V (i/e domain)")

        cat("\nDERIVATIONAL SUMMARY TABLE (N -> V, i/e DOMAIN)\n")
        nv_tp_full <- make_deriv_summary(
          denom_ie, "mutation_deriv_verb",
          "N->V derivation, base plural mutated",
          "N->V derivation, base plural non-mut."
        )
        print_full(nv_tp_full)
      }

      # Analysis 2: Steriade allomorph-class test
      # Compare front (-i/-e) vs back (-uri) plural nouns.
      # Ambiguous (no plural / same-as-singular) are reported
      # separately but excluded from the Steriade test.
      cat("\n--- Steriade allomorph-class test ---\n")
      cat("Does plural allomorph class (front vs back) predict\n")
      cat("verbalizer palatalization?\n\n")

      cat("All three classes:\n")
      denom_pairs |>
        group_by(plural_class) |>
        summarise(
          N = n(),
          N_deriv_pal = sum(mutation_deriv_verb),
          rate = N_deriv_pal / N,
          .groups = "drop"
        ) |>
        print_full()

      # Steriade test: front vs back only (exclude ambiguous)
      denom_steriade <- denom_pairs |>
        filter(plural_class %in% c("front", "back"))

      if (nrow(denom_steriade) > 0 &&
        n_distinct(denom_steriade$plural_class) == 2L) {
        cat("\nFront vs back (Steriade test):\n")
        denom_steriade |>
          group_by(plural_class) |>
          summarise(
            N = n(), N_pal = sum(mutation_deriv_verb),
            rate = N_pal / N, .groups = "drop"
          ) |>
          print_full()

        cat("\nFirth logistic: deriv_mut ~ plural_class\n")
        model_s <- glm(
          mutation_deriv_verb ~ plural_class,
          data = denom_steriade,
          family = binomial(),
          method = brglm2::brglmFit,
          control = brglm2::brglm_control(maxit = 500)
        )
        print_full(broom::tidy(model_s))
        or <- exp(coef(model_s)["plural_classfront"])
        ci <- exp(confint.default(model_s)["plural_classfront", ])
        cat(sprintf(
          "  OR(front vs back) = %.3f [95%% CI: %.3f, %.3f]\n",
          or, ci[1], ci[2]
        ))

        cat("\nFirth logistic: deriv_mut ~ plural_class + log(freq)\n")
        denom_s_freq <- denom_steriade |>
          mutate(log_freq = log(lemma_freq + 1))
        model_sf <- glm(
          mutation_deriv_verb ~ plural_class + log_freq,
          data = denom_s_freq,
          family = binomial(),
          method = brglm2::brglmFit,
          control = brglm2::brglm_control(maxit = 500)
        )
        print_full(broom::tidy(model_sf))
      }

      # TP-style summary by plural class
      cat("\nTP TABLE BY PLURAL CLASS (N -> V)\n")
      denom_pairs |>
        group_by(plural_class) |>
        summarise(
          mutated     = sum(mutation_deriv_verb),
          non_mutated = sum(!mutation_deriv_verb),
          .groups     = "drop"
        ) |>
        mutate(type = paste0("N->V, ", plural_class, " plural")) |>
        tp_table(type, mutated, non_mutated) |>
        select(type, N, mutated, non_mutated, rate,
          majority = majority_mutates, tolerated
        ) |>
        print_full()

      # --- Decomposition: suffix selection vs palatalization ---
      cat("\n--- SUFFIX SELECTION (what Steriade actually predicts) ---\n")
      cat("Does plural class predict WHICH suffix (-i vs -a) is chosen?\n\n")
      denom_pairs |>
        filter(plural_class %in% c("front", "back")) |>
        group_by(plural_class, verb_suffix_front) |>
        summarise(N = n(), .groups = "drop") |>
        mutate(suffix = if_else(verb_suffix_front,
          "front (-i/-ui)", "back (-a)"
        )) |>
        select(plural_class, suffix, N) |>
        print_full()

      cat("\n--- SEGMENT-CLASS DECOMPOSITION ---\n")
      cat("Palatalization rate by plural class AND segment class:\n\n")
      denom_pairs |>
        filter(plural_class %in% c("front", "back")) |>
        group_by(seg_class, plural_class) |>
        summarise(
          N = n(),
          N_pal = sum(mutation_deriv_verb),
          rate = N_pal / N,
          .groups = "drop"
        ) |>
        arrange(seg_class, plural_class) |>
        print_full()

      # Controlled regression: does plural_class predict
      # palatalization AFTER controlling for segment class?
      if (nrow(denom_steriade) > 5 &&
        n_distinct(denom_steriade$seg_class) == 2L) {
        cat(
          "\nFirth logistic (controlled): deriv_mut ~",
          "plural_class + seg_class\n"
        )
        denom_ctrl <- denom_steriade |>
          mutate(seg_class = factor(seg_class,
            levels = c("coronal", "dorsal")
          ))
        model_ctrl <- glm(
          mutation_deriv_verb ~ plural_class + seg_class,
          data = denom_ctrl,
          family = binomial(),
          method = brglm2::brglmFit,
          control = brglm2::brglm_control(maxit = 500)
        )
        print_full(broom::tidy(model_ctrl))
        or_pc <- exp(coef(model_ctrl)["plural_classfront"])
        ci_pc <- exp(confint.default(
          model_ctrl
        )["plural_classfront", ])
        cat(sprintf(
          "  OR(front vs back | seg_class) = %.3f [95%% CI: %.3f, %.3f]\n",
          or_pc, ci_pc[1], ci_pc[2]
        ))
        cat("  If OR ~ 1 and p > 0.05, the Steriade effect is a\n")
        cat("  composition artifact of segment-class distribution.\n")
      }

      # Suffix selection test (Fisher's exact):
      # Does whether the plural alternates predict -i vs -a verbalizer?
      # NOTE: row is coded by `mutates_class` (actual alternation), not
      # `plural_class` (allomorph proxy). For coronals, -e plurals do not
      # alternate, so they belong with non-alternating bases despite the
      # surface allomorph being "front". Replicates Steriade (10.20)
      # affix-avoidance test correctly.
      cat("\n--- SUFFIX SELECTION: Fisher exact test ---\n")
      cat("  Row = plural ACTUALLY alternates (mutation), levels: alt > nonalt\n")
      cat("  Col = verbalizer is -i (TRUE) vs -a/-ui (FALSE)\n")
      cat("  OR is reported in Steriade direction:\n")
      cat("    OR > 1 means alt-pl bases prefer the -i verbalizer.\n")
      cat("  NOTE: This R analysis does NOT apply etymology validation;\n")
      cat("        for the abstract's canonical numbers see\n")
      cat("        scripts/build_contingency_table.py (n=388 etym-validated).\n")
      suf_tab <- denom_pairs |>
        filter(mutates_class %in% c("alt", "nonalt")) |>
        mutate(
          suf_front = verb_suffix_front,
          # Order factor so OR comes out in Steriade direction
          mutates_class = factor(mutates_class, levels = c("nonalt", "alt"))
        )
      tab <- with(suf_tab, table(mutates_class, suf_front))
      if (all(dim(tab) == c(2L, 2L))) {
        ft <- fisher.test(tab)
        cat("  POOLED  n =", sum(tab), "  OR =", round(ft$estimate, 3),
            "  p =", format.pval(ft$p.value, digits = 3), "\n")
        cat("  Per-class breakdown:\n")
        for (sc in c("dorsal", "coronal")) {
          sub_tab <- with(filter(suf_tab, seg_class == sc),
                          table(mutates_class, suf_front))
          if (all(dim(sub_tab) == c(2L, 2L))) {
            ft_s <- fisher.test(sub_tab)
            cat(sprintf("    %-7s n=%d  OR=%.3f  p=%s\n",
                        sc, sum(sub_tab), ft_s$estimate,
                        format.pval(ft_s$p.value, digits = 3)))
          }
        }
      } else {
        cat("  Contingency table not 2x2; printing instead:\n")
        print(tab)
      }

      # -----------------------------------------------------------------
      # IAP ROOT-SPECIFICATION ANALYSIS
      # -----------------------------------------------------------------
      # Under IAP, root specification (underspecified /K/ vs fully
      # specified /k/) independently determines:
      #   (a) whether the plural allomorph is front (-i/-e) or back (-uri)
      #   (b) whether the verbalizer allomorph is front (-i) or back (-a/-ui)
      #   (c) whether the consonant palatalizes given a front trigger
      #
      # We infer root specification from plural behavior:
      #   - mutation=TRUE or nde:gimpe → underspecified /K,G/ or /T,S,D/
      #   - -uri plural without palatalization → fully specified /k,g/
      #   - Verbs that palatalize despite -uri plural (buluc→buluci)
      #     → underspecified /K/ with listed -uri exception
      cat("\n--- IAP ROOT-SPECIFICATION ANALYSIS (N -> V) ---\n")

      iap <- denom_pairs |>
        mutate(
          root_spec = case_when(
            # Underspecified: palatalizes in plural OR tautomorphemic
            mutation_inflect ~ "underspecified",
            str_starts(exception_reason, "nde:gimpe") ~ "underspecified",
            # Verbalizer palatalizes despite -uri → underspecified
            # with listed -uri plural
            plural_class == "back" & mutation_deriv_verb ~ "underspecified",
            # Fully specified: -uri plural, no palatalization anywhere
            plural_class == "back" & !mutation_deriv_verb ~ "fully_specified",
            # Paduchi/ochi = fully specified (NDEB)
            str_starts(exception_reason, "nde:paduchi") ~ "fully_specified",
            str_starts(exception_reason, "nde:ochi") ~ "fully_specified",
            TRUE ~ "ambiguous"
          ),
          # Classify the verbalizer suffix more precisely.
          # -iza/-ifica start with /i/ and trigger palatalization
          # of dorsals, so they pattern with front (-i/-ui).
          # -ăni/-ări/-arisi do not provide a front-vowel trigger
          # adjacent to the stem-final consonant.
          verb_suffix_type = case_when(
            str_ends(derived_verbs, "\u0103ni") ~ "other",
            str_ends(derived_verbs, "\u0103ri") ~ "other",
            str_ends(derived_verbs, "arisi") ~ "other",
            str_ends(derived_verbs, "iza") ~ "front (-iza/-ifica)",
            str_ends(derived_verbs, "ifica") ~ "front (-iza/-ifica)",
            verb_suffix_front ~ "front (-i/-ui)",
            TRUE ~ "back (-a)"
          )
        )

      iap_classified <- filter(iap, root_spec != "ambiguous")

      cat("\nRoot specification (dorsals):\n")
      iap_classified |>
        filter(seg_class == "dorsal") |>
        group_by(root_spec, verb_suffix_type) |>
        summarise(
          N = n(),
          N_pal = sum(mutation_deriv_verb),
          .groups = "drop"
        ) |>
        arrange(root_spec, verb_suffix_type) |>
        print_full()

      # The key IAP table: phonology is exceptionless
      cat("\nPhonological prediction (dorsals):\n")
      cat("  underspecified + front trigger → palatalizes?\n")
      cat("  fully specified + any trigger → no palatalization?\n\n")
      iap_dorsal <- iap_classified |> filter(seg_class == "dorsal")
      if (nrow(iap_dorsal) > 0) {
        iap_dorsal |>
          group_by(root_spec, mutation_deriv_verb) |>
          summarise(N = n(), .groups = "drop") |>
          arrange(root_spec, mutation_deriv_verb) |>
          print_full()

        # Separate phonology from morphology.
        # "Front trigger" includes -i/-ui AND -iza/-ifica
        # because all of these put a front vowel adjacent to
        # the stem-final consonant.
        underspec_front <- iap_dorsal |>
          filter(
            root_spec == "underspecified",
            verb_suffix_type %in% c(
              "front (-i/-ui)",
              "front (-iza/-ifica)"
            )
          )
        cat(sprintf(
          "\n  Underspecified + front trigger: %d/%d palatalize (%.1f%%)\n",
          sum(underspec_front$mutation_deriv_verb),
          nrow(underspec_front),
          100 * mean(underspec_front$mutation_deriv_verb)
        ))

        fullspec <- iap_dorsal |>
          filter(root_spec == "fully_specified")
        cat(sprintf(
          "  Fully specified (any trigger): %d/%d palatalize (%.1f%%)\n",
          sum(fullspec$mutation_deriv_verb),
          nrow(fullspec),
          100 * mean(fullspec$mutation_deriv_verb)
        ))

        cat("\n  Suffix selection exceptions (underspecified + back/other):\n")
        underspec_back <- iap_dorsal |>
          filter(
            root_spec == "underspecified",
            verb_suffix_type != "front (-i/-ui)"
          )
        if (nrow(underspec_back) > 0) {
          cat(sprintf(
            "    %d/%d underspecified roots chose non-front suffix\n",
            nrow(underspec_back),
            sum(iap_dorsal$root_spec == "underspecified")
          ))
          underspec_back |>
            select(
              lemma, plural, derived_verbs,
              verb_suffix_type
            ) |>
            print_full()
        }
      }

      if (any(iap_classified$seg_class == "coronal")) {
        cat("\nRoot specification (coronals):\n")
        iap_classified |>
          filter(seg_class == "coronal") |>
          group_by(root_spec, mutation_deriv_verb) |>
          summarise(N = n(), .groups = "drop") |>
          arrange(root_spec, mutation_deriv_verb) |>
          print_full()
      }
    }
  } else {
    cat("No denominal verb columns found; skipping N -> V comparison.\n")
  }

  # -----------------------------------------------------------------------
  # (2) Noun → Adjective
  # -----------------------------------------------------------------------
  cat_section("(2) NOUN LEMMAS: INFLECTIONAL PLURALS VS DENOMINAL ADJECTIVES")

  if (has_adj_deriv_cols) {
    noun_adj_pairs <- noun_base_inflect |>
      filter(
        !is.na(derived_adj), derived_adj != "",
        !str_starts(coalesce(exception_reason, ""), "nde:")
      ) |>
      mutate(mutation_deriv_adj = detect_palatal_from_orth(lemma, derived_adj, stem_final)) |>
      filter(!is.na(mutation_deriv_adj)) |>
      arrange(lemma) |>
      distinct(lemma, .keep_all = TRUE)

    cat("Denominal N-Adj lemmas with target consonants:", nrow(noun_adj_pairs), "\n")
    cat("  front-vowel plural (i/e):", sum(noun_adj_pairs$plural_class == "front"), "\n")
    cat("  non-front plural (uri/none):", sum(noun_adj_pairs$plural_class == "non_front"), "\n")

    if (nrow(noun_adj_pairs) > 0) {
      # i/e domain subset
      nadj_ie <- filter(noun_adj_pairs, plural_class == "front")
      if (nrow(nadj_ie) > 0) {
        cat("\n--- i/e domain subset ---\n")
        analyze_inf_vs_deriv(nadj_ie, "mutation_inflect", "mutation_deriv_adj", "N -> Adj (i/e domain)")
      }

      # Steriade allomorph-class test (front vs back)
      nadj_steriade <- noun_adj_pairs |>
        filter(plural_class %in% c("front", "back"))
      if (nrow(nadj_steriade) > 0) {
        cat("\n--- Steriade allomorph-class test (front vs back) ---\n")
        nadj_steriade |>
          group_by(plural_class) |>
          summarise(
            N = n(), N_pal = sum(mutation_deriv_adj),
            rate = N_pal / N, .groups = "drop"
          ) |>
          print_full()
      }
      if (any(noun_adj_pairs$plural_class == "ambig")) {
        cat("\nAmbiguous plural class (excluded from Steriade test):\n")
        noun_adj_pairs |>
          filter(plural_class == "ambig") |>
          summarise(
            N = n(), N_pal = sum(mutation_deriv_adj),
            rate = N_pal / N
          ) |>
          print_full()
      }

      cat("\nDERIVATIONAL SUMMARY TABLE (N -> Adj, FULL LEXICON)\n")
      na_tp_full <- make_deriv_summary(
        noun_adj_pairs, "mutation_deriv_adj",
        "N->Adj derivation, base plural mutated",
        "N->Adj derivation, base plural non-mut."
      )
      print_full(na_tp_full)
    }
  } else {
    cat("No denominal adjective columns found; skipping N -> Adj comparison.\n")
  }

  # -----------------------------------------------------------------------
  # (3) Adjective lemmas — kept for completeness; consistently empty
  # -----------------------------------------------------------------------
  cat_section("(3) ADJECTIVE LEMMAS: INFLECTIONAL PLURALS VS DERIVATIONS")

  adj_base_inflect <- lex |>
    filter(
      pos == "ADJ", opportunity %in% plural_opportunities,
      !is.na(stem_final), stem_final %in% segments_of_interest, !is.na(mutation)
    ) |>
    mutate(mutation_inflect = as.logical(mutation), lemma_freq = freq_ron_news_2024_1M)

  cat("Adjective lemmas with i/e plural & target segments:", nrow(adj_base_inflect), "\n")

  if (has_verb_deriv_cols && nrow(adj_base_inflect) > 0) {
    adj_verb_pairs <- adj_base_inflect |>
      filter(
        !is.na(derived_verbs), derived_verbs != "",
      ) |>
      mutate(
        verb_suffix_front = deriv_suffixes %in% front_verb_suffixes,
        mutation_deriv_verb = detect_palatal_from_orth(lemma, derived_verbs, stem_final)
      ) |>
      filter(!is.na(mutation_deriv_verb)) |>
      arrange(lemma) |>
      distinct(lemma, .keep_all = TRUE)

    cat("Adjective lemmas with both plural + derived verb:", nrow(adj_verb_pairs), "\n")
    if (nrow(adj_verb_pairs) > 0) {
      analyze_inf_vs_deriv(adj_verb_pairs, "mutation_inflect", "mutation_deriv_verb", "Adj -> V")
    }
  }

  if (has_adj_deriv_cols && nrow(adj_base_inflect) > 0) {
    adj_adj_pairs <- adj_base_inflect |>
      filter(
        !is.na(derived_adj), derived_adj != "",
        !is.na(ipa_derived_adj), ipa_derived_adj != "", derived_adj != lemma
      ) |>
      mutate(mutation_deriv_adj = detect_palatal_from_orth(lemma, derived_adj, stem_final)) |>
      filter(!is.na(mutation_deriv_adj)) |>
      arrange(lemma) |>
      distinct(lemma, .keep_all = TRUE)

    cat("Adjective lemmas with plural + non-trivial derived Adj:", nrow(adj_adj_pairs), "\n")
    if (nrow(adj_adj_pairs) > 0) {
      analyze_inf_vs_deriv(adj_adj_pairs, "mutation_inflect", "mutation_deriv_adj", "Adj -> Adj")
    }
  }
}

# =========================================================================
# Descriptive Summaries
# =========================================================================

cat_section("SEGMENT-WISE MUTATION RATES (PRODUCTIVE TP DOMAIN)")
nouns_tp |>
  group_by(stem_final) |>
  summarise(N_opp = n(), N_mut = sum(mutation, na.rm = TRUE), .groups = "drop") |>
  calc_rate() |>
  arrange(stem_final) |>
  print_full()

cat_section("MUTATION RATES BY SEGMENT AND PLURAL TYPE (PRODUCTIVE TP DOMAIN)")
nouns_tp |>
  group_by(stem_final, opportunity) |>
  summarise(N_opp = n(), N_mut = sum(mutation, na.rm = TRUE), .groups = "drop") |>
  calc_rate() |>
  arrange(stem_final, opportunity) |>
  print_full()

cat_section("CLUSTER INVENTORY (PRODUCTIVE TP DOMAIN)")
nouns_tp |>
  filter(!is.na(cluster), cluster != "") |>
  count(stem_final, cluster, sort = TRUE) |>
  print_full()

cat_section("CLUSTER EFFECTS ON MUTATION (PRODUCTIVE TP DOMAIN)")
nouns_tp |>
  mutate(cluster_type = if_else(
    cluster %in% c("st", "sc", "ct", "chi", "che", "ghi", "ghe"), cluster, "none"
  )) |>
  group_by(stem_final, cluster_type, opportunity) |>
  summarise(N_opp = n(), N_mut = sum(mutation, na.rm = TRUE), .groups = "drop") |>
  calc_rate() |>
  arrange(cluster_type, opportunity, stem_final) |>
  print_full()

cat_section("STRUCTURE OF NON-UNDERGOERS (I/E DOMAIN)")
nouns_opp |>
  filter(!mutation) |>
  group_by(stem_final, exception_category) |>
  summarise(N = n(), .groups = "drop") |>
  group_by(stem_final) |>
  mutate(total_non_under = sum(N), prop = N / total_non_under) |>
  arrange(exception_category, stem_final, desc(prop)) |>
  print_full()

cat_section("NDE DISTRIBUTION BY SEGMENT")
ndeb_rows |>
  group_by(nde_class, stem_final) |>
  summarise(N = n(), .groups = "drop") |>
  arrange(nde_class, stem_final) |>
  print_full()

cat_section("TRUE EXCEPTIONS IN I/E DOMAIN (NON-NDEB)")
true_exc <- nouns |>
  filter(
    exception_reason == "unexplained", opportunity %in% plural_opportunities,
    !(nde_class %in% ndeb_classes)
  )
cat("Flagged true-exception nouns in i/e domain (non-NDEB):", nrow(true_exc), "\n")
true_exc |>
  group_by(stem_final) |>
  summarise(N_true_exc = n(), .groups = "drop") |>
  arrange(stem_final) |>
  print_full()

if (nrow(true_exc) > 0) {
  cat("\nSample of true exceptions:\n")
  true_exc |>
    select(lemma, plural, stem_final, opportunity, nde_class, lemma_suffix, notes) |>
    head(30) |>
    print()
  cat("\nSuffix distribution among true exceptions (top 10):\n")
  true_exc |>
    count(lemma_suffix, sort = TRUE) |>
    head(10) |>
    print()
}

# Machine-verifiable assertion: dorsals have ZERO true exceptions anywhere in
# the i/e domain.  This is a stronger claim than the TP alone can make.
dorsal_true_exc <- nouns_opp |>
  filter(stem_final %in% c("c", "g"), exception_category == "true_exception")
stopifnot(nrow(dorsal_true_exc) == 0)
cat("CONFIRMED: dorsals have zero true exceptions in i/e domain\n")

# s+e decomposition: confirm that all s+e mutations come from sc clusters,
# i.e. bare /s/ before -e never mutates.
cat("\ns+e mutation decomposed by cluster presence:\n")
nouns_tp |>
  filter(stem_final == "s", opportunity == "e") |>
  mutate(has_cluster = !is.na(cluster_simple)) |>
  group_by(has_cluster) |>
  summarise(N = n(), N_mut = sum(mutation), rate = N_mut / N, .groups = "drop") |>
  print_full()

# =========================================================================
# Two-System Analysis: Dorsal (vowel-triggered) vs Coronal (consonant-triggered)
# =========================================================================

cat_section("TWO-SYSTEM ANALYSIS: DORSAL VS CORONAL TRIGGERS")

cat("Dorsals palatalize before BOTH -i and -e (vowel-triggered).\n")
cat("Coronals palatalize before -i ONLY (consonant-triggered via [j]).\n")
cat("This table is the core evidence for the two-system architecture.\n\n")

# All nouns with target consonants (TP domain), by segment class × opportunity
two_sys <- nouns_tp |>
  mutate(seg_class = if_else(stem_final %in% c("c", "g"), "dorsal", "coronal")) |>
  group_by(seg_class, stem_final, opportunity) |>
  summarise(
    N = n(), N_mut = sum(mutation), rate = N_mut / N, .groups = "drop"
  ) |>
  arrange(seg_class, stem_final, opportunity)

cat("TP domain (underspecified consonants only):\n")
print_full(two_sys)

# Summary: dorsal vs coronal aggregated
cat("\nAggregated by class × opportunity:\n")
nouns_tp |>
  mutate(seg_class = if_else(stem_final %in% c("c", "g"), "dorsal", "coronal")) |>
  group_by(seg_class, opportunity) |>
  summarise(
    N = n(), N_mut = sum(mutation), rate = N_mut / N, .groups = "drop"
  ) |>
  arrange(seg_class, opportunity) |>
  print_full()

cat("\nKey contrast:\n")
cat("  Dorsals: 100% before -i, 100% before -e\n")
cat("  Coronals: ~100% before -i, ~0% before -e\n")
cat("  The difference: /i/ glides to [j] (consonantal trigger);\n")
cat("  /e/ does not glide, so coronals have no trigger.\n")

# =========================================================================
# Cluster Behavior: Evidence for Rule Ordering
# =========================================================================

cat_section("CLUSTER BEHAVIOR: EVIDENCE FOR RULE ORDERING")

cat("Each cluster type tests a different aspect of the rule ordering.\n\n")

# Collect all cluster data from the FULL noun set (not just TP domain)
# since clusters have small N.  Use startsWith to catch scă, stă, etc.
cluster_data <- nouns |>
  filter(!is.na(cluster), cluster != "") |>
  mutate(cluster_type = case_when(
    startsWith(cluster, "st") ~ "st",
    startsWith(cluster, "sc") ~ "sc",
    startsWith(cluster, "ct") ~ "ct",
    TRUE ~ NA_character_
  )) |>
  filter(!is.na(cluster_type)) |>
  group_by(cluster_type, opportunity) |>
  summarise(
    N = n(),
    N_mut = sum(mutation, na.rm = TRUE),
    rate = if_else(N > 0L, N_mut / N, NA_real_),
    .groups = "drop"
  ) |>
  filter(N > 0) |>
  arrange(cluster_type, opportunity)

print_full(cluster_data)

cat("\nCluster derivation summary:\n")
cat("  ┌─────────┬──────┬────────┬─────────────────────────────────┐\n")
cat("  │ Cluster │ Sfx  │ Output │ Mechanism                       │\n")
cat("  ├─────────┼──────┼────────┼─────────────────────────────────┤\n")
cat("  │ st      │ -i   │ [ʃt]   │ Glide→S-pal (SEARCH skips T)   │\n")
cat("  │         │      │        │ →bleeding blocks assibilation   │\n")
cat("  │ st      │ -e   │ [st]   │ No consonantal trigger          │\n")
cat("  │ ct      │ -i   │ [kts]  │ /k/ inalterable; /T/ assibilates│\n")
cat("  │ ct      │ -e   │ [kt]   │ No trigger for either rule      │\n")
cat("  │ sc      │ -i   │ [ʃt]   │ K→tʃ; S-pal triggered by tʃ;   │\n")
cat("  │         │      │        │ bleeding→t                      │\n")
cat("  │ sc      │ -e   │ [ʃt]★  │ K→tʃ before [+front]; tʃ       │\n")
cat("  │         │      │        │ triggers S-pal; bleeding→t      │\n")
cat("  └─────────┴──────┴────────┴─────────────────────────────────┘\n")
cat("  ★ KEY DIAGNOSTIC: sc+e proves the coronal trigger is [tʃ],\n")
cat("    not the vowel — no [j] present, yet S palatalizes.\n")
cat("    This cannot be explained by a vowel-based trigger.\n")
cat("\n  Bleeding converts [tʃ]→[t] after [ʃ] by stripping\n")
cat("  [+strident]: /SK/+e → [tʃ] (K-pal) → [ʃ][tʃ] (S-pal)\n")
cat("  → [ʃ][t] (bleeding removes +strid from tʃ).\n")

# Verify: s+e = 0% but sc+e = 100% (the key diagnostic)
cat("\n--- KEY DIAGNOSTIC: s+e vs sc+e ---\n")
s_e_all <- nouns |>
  filter(stem_final == "s", opportunity == "e")
s_e_simple <- s_e_all |>
  filter(is.na(cluster) | cluster == "" | !startsWith(cluster, "sc"))
s_e_sc <- s_e_all |>
  filter(!is.na(cluster), startsWith(cluster, "sc"))
cat("  Simple s+e:", sum(s_e_simple$mutation), "/", nrow(s_e_simple), "mutate")
cat(if (nrow(s_e_simple) > 0) sprintf(" (%.1f%%)\n", sum(s_e_simple$mutation) / nrow(s_e_simple) * 100) else "\n")
cat("  sc+e:      ", sum(s_e_sc$mutation), "/", nrow(s_e_sc), "mutate")
cat(if (nrow(s_e_sc) > 0) sprintf(" (%.1f%%)\n", sum(s_e_sc$mutation) / nrow(s_e_sc) * 100) else "\n")
cat("  Simple s NEVER palatalizes before /e/ (no consonantal trigger).\n")
cat("  SC always palatalizes before /e/ (tʃ provides the trigger).\n\n")

# Verify: t assibilates in all non-st cluster contexts
cat("--- t assibilation by cluster context ---\n")
cat("Does t assibilate before -i in each cluster type?\n\n")
nouns |>
  filter(
    opportunity == "i", mutation,
    str_ends(lemma, "nt|lt|rt|ft|ct|st")
  ) |>
  mutate(
    cluster_ctx = str_extract(lemma, "(nt|lt|rt|ft|ct|st)$")
  ) |>
  filter(!is.na(cluster_ctx)) |>
  group_by(cluster_ctx) |>
  summarise(N = n(), .groups = "drop") |>
  mutate(t_assibilates = cluster_ctx != "st") |>
  arrange(cluster_ctx) |>
  print_full()

# =========================================================================
# Transderivational Table (4): Inflection × Derivation Cross-Tab
# =========================================================================

cat_section("TRANSDERIVATIONAL TABLE (4): INFLECTION × DERIVATION")

if (has_verb_deriv_cols && exists("denom_pairs")) {
  cat("Cross-tabulation of plural palatalization × verb palatalization\n")
  cat("by stem-final segment (all nouns with derived verbs, excl. NDEB).\n\n")

  denom_clean <- denom_pairs |>
    filter(!str_starts(exception_reason, "nde:"))

  if (nrow(denom_clean) > 0) {
    cat("By segment:\n")
    denom_clean |>
      mutate(
        both = mutation_inflect & mutation_deriv_verb,
        infl_only = mutation_inflect & !mutation_deriv_verb,
        deriv_only = !mutation_inflect & mutation_deriv_verb,
        neither = !mutation_inflect & !mutation_deriv_verb
      ) |>
      group_by(stem_final) |>
      summarise(
        n = n(),
        both = sum(both),
        infl_only = sum(infl_only),
        deriv_only = sum(deriv_only),
        neither = sum(neither),
        agree_pct = round((both + neither) / n * 100),
        .groups = "drop"
      ) |>
      arrange(stem_final) |>
      print_full()

    cat("\nBy segment class:\n")
    denom_clean |>
      mutate(
        seg_class = if_else(stem_final %in% c("c", "g"), "dorsal", "coronal"),
        both = mutation_inflect & mutation_deriv_verb,
        infl_only = mutation_inflect & !mutation_deriv_verb,
        deriv_only = !mutation_inflect & mutation_deriv_verb,
        neither = !mutation_inflect & !mutation_deriv_verb
      ) |>
      group_by(seg_class) |>
      summarise(
        n = n(),
        both = sum(both),
        infl_only = sum(infl_only),
        deriv_only = sum(deriv_only),
        neither = sum(neither),
        agree_pct = round((both + neither) / n * 100),
        .groups = "drop"
      ) |>
      arrange(seg_class) |>
      print_full()
  }
}

cat_section("NDEB EXCEPTIONS OF OCHI/PADUCHE TYPE")
ndeb_exc_ochi_pad <- filter(nouns, nde_class %in% ndeb_observable, !mutation)
cat("Non-mutating NDEB nouns of ochi/paduche type:", nrow(ndeb_exc_ochi_pad), "\n")
if (nrow(ndeb_exc_ochi_pad) > 0) {
  ndeb_exc_ochi_pad |>
    select(lemma, plural, stem_final, opportunity, nde_class, lemma_suffix, notes) |>
    arrange(nde_class, stem_final, lemma) |>
    print_full()
}

# =========================================================================
# Frequency and Exception Status
# =========================================================================
# Key test: does log(lemma_freq) predict non-mutation in the TP domain, after
# conditioning on segment and opportunity?
#
# If NOT significant: exceptions are distributed across the frequency range,
# consistent with phonological (underspecification) conditioning rather than
# high-frequency lexical entrenchment (Bybee 1995).
#
# If significant: a frequency-of-lexicalization story cannot be ruled out for
# the residual exceptions, and this warrants a footnote.
#
# Dorsals are excluded: they have zero exceptions regardless of frequency,
# so including them would trivially suppress any frequency effect without
# adding information.  The z-segment is included as a strong negative control
# (categorically non-mutating; any frequency effect there would indicate
# that the measure is picking up something other than exception status).

cat_section("FREQUENCY AND EXCEPTION STATUS (TP DOMAIN)")

freq_exc_data <- nouns_tp |>
  filter(
    stem_final %in% c("d", "s", "t", "z"), # coronals only; dorsals have no exceptions
    opportunity %in% plural_opportunities,
    !is.na(freq_ron_news_2024_1M),
    freq_ron_news_2024_1M > 0,
    !is.na(mutation)
  ) |>
  mutate(
    is_exception = as.integer(!mutation),
    log_freq     = log(freq_ron_news_2024_1M),
    seg_opp      = paste(stem_final, opportunity, sep = "_")
  )

cat("Items with frequency data (coronal TP domain):", nrow(freq_exc_data), "\n")
cat("Segment x opportunity breakdown:\n")
freq_exc_data |>
  group_by(stem_final, opportunity) |>
  summarise(N = n(), N_exc = sum(is_exception), .groups = "drop") |>
  print_full()

if (nrow(freq_exc_data) >= 20L) {
  cat("\nFirth logistic: exception ~ log_freq + stem_final * opportunity\n")
  model_freq_exc <- glm(
    is_exception ~ log_freq + stem_final * opportunity,
    data = freq_exc_data,
    family = binomial(),
    method = brglm2::brglmFit,
    control = brglm2::brglm_control(maxit = 500)
  )
  print_full(broom::tidy(model_freq_exc))

  freq_main <- broom::tidy(model_freq_exc) |> filter(term == "log_freq")
  cat(sprintf(
    "\nlog_freq coefficient: b = %.4f, p = %.4f\n",
    freq_main$estimate, freq_main$p.value
  ))
  if (freq_main$p.value > 0.05) {
    cat("=> frequency does NOT predict exception status after conditioning on segment/opportunity.\n")
    cat("   This supports phonological (not usage-based) conditioning of non-mutation.\n")
  } else {
    cat("=> frequency IS a significant predictor (consider footnote).\n")
    cat("   Direction:", if_else(freq_main$estimate > 0,
      "higher-frequency items are MORE likely to be exceptions",
      "higher-frequency items are LESS likely to be exceptions"
    ), "\n")
  }

  # Descriptive: median frequency of exceptions vs non-exceptions by segment
  cat("\nMedian lemma frequency: exceptions vs non-exceptions (coronal TP domain)\n")
  freq_exc_data |>
    group_by(stem_final, opportunity, is_exception) |>
    summarise(
      N           = n(),
      median_freq = median(freq_ron_news_2024_1M),
      .groups     = "drop"
    ) |>
    arrange(stem_final, opportunity, is_exception) |>
    print_full()
} else {
  cat("Insufficient data for frequency ~ exception model.\n")
}

# =========================================================================
# Frequency-Based Downsampling
# =========================================================================

cat("DOWNSAMPLING (FREQUENCY-BASED) FOR TOLERANCE PRINCIPLE ANALYSIS\n")

nouns_opp_down_single <- build_downsampled_lexicon(nouns_opp, n_lex = 1000L)
cat("Unique lemmas in i/e domain (all):", n_distinct(nouns_opp$lemma), "\n\n")

# =========================================================================
# Tolerance Principle: Full Data
# =========================================================================

cat_section("TOLERANCE PRINCIPLE: SEGMENT-LEVEL PATTERNS (PRODUCTIVE TP DOMAIN)")
# Exclude cluster items from simple-segment rows so they don't double-count
# (e.g. a `prost` (st cluster, stem_final='s') would otherwise appear in BOTH
#  the simple `s` row and the `st` cluster row).
seg_tp_all <- compute_segment_tp_tables(nouns_tp |> filter(is.na(cluster) | cluster == ""))
print_full(seg_tp_all)
if (any(seg_tp_all$small_n, na.rm = TRUE)) {
  cat("\nSmall-N cells (N <", SMALL_CELL_THRESHOLD, ") — TP binary unreliable; see Bayesian TP:\n")
  seg_tp_all |>
    filter(small_n) |>
    select(type, N, mutated, non_mutated, tolerated) |>
    print_full()
}

# =========================================================================
# NDEB Contribution
# =========================================================================

cat_section("TOLERANCE PRINCIPLE: NDEB CONTRIBUTION PER TYPE (FULL LEXICON)")
seg_tp_all_ndeb <- compute_segment_tp_tables(nouns_opp |> filter(nde_class %in% ndeb_classes),
  label_suffix = " (NDEB)"
)
print_full(seg_tp_all_ndeb)

cat_section("TOLERANCE PRINCIPLE: NDEB CONTRIBUTION PER TYPE (DOWNSAMPLED)")
if (!is.null(nouns_opp_down_single) && nrow(nouns_opp_down_single) > 0L) {
  compute_segment_tp_tables(nouns_opp_down_single |> filter(nde_class %in% ndeb_classes),
    label_suffix = " (downsampled, NDEB)"
  ) |> print_full()
} else {
  cat("No downsampled lexicon available.\n")
}

# =========================================================================
# Tolerance Principle + Derivational: Downsampled
# =========================================================================

if (!is.null(nouns_opp_down_single) && nrow(nouns_opp_down_single) > 0L) {
  lemmas_ds <- unique(nouns_opp_down_single$lemma)

  cat_section("DERIVATIONAL ANALYSIS (REFERENCE DOWNSAMPLED LEXICON: TOP 1000)")

  if (has_verb_deriv_cols && exists("noun_base_inflect")) {
    denom_pairs_ds <- noun_base_inflect |>
      filter(
        lemma %in% lemmas_ds,
        !is.na(derived_verbs), derived_verbs != "",
        !str_starts(coalesce(exception_reason, ""), "nde:")
      ) |>
      mutate(
        verb_suffix_front = deriv_suffixes %in% front_verb_suffixes,
        mutation_deriv_verb = detect_palatal_from_orth(lemma, derived_verbs, stem_final)
      ) |>
      filter(!is.na(mutation_deriv_verb)) |>
      arrange(lemma) |>
      distinct(lemma, .keep_all = TRUE)

    cat("Denominal N-V lemmas (downsampled):", nrow(denom_pairs_ds), "\n")
    if (nrow(denom_pairs_ds) > 0) {
      # Full analysis on downsampled pairs (not just summary table)
      analyze_inf_vs_deriv(denom_pairs_ds, "mutation_inflect", "mutation_deriv_verb", "N -> V (downsampled)")

      nv_tp_ds <- make_deriv_summary(
        denom_pairs_ds, "mutation_deriv_verb",
        "N->V derivation, base plural mutated (downsampled)",
        "N->V derivation, base plural non-mut. (downsampled)"
      )
      print_full(nv_tp_ds)
    }
  }

  if (has_adj_deriv_cols && exists("noun_base_inflect")) {
    noun_adj_pairs_ds <- noun_base_inflect |>
      filter(
        lemma %in% lemmas_ds,
        !is.na(derived_adj), derived_adj != "",
        !is.na(ipa_derived_adj), ipa_derived_adj != "",
        !str_starts(coalesce(exception_reason, ""), "nde:")
      ) |>
      mutate(mutation_deriv_adj = detect_palatal_from_orth(lemma, derived_adj, stem_final)) |>
      filter(!is.na(mutation_deriv_adj)) |>
      arrange(lemma) |>
      distinct(lemma, .keep_all = TRUE)

    cat("Denominal N-Adj lemmas (downsampled):", nrow(noun_adj_pairs_ds), "\n")
    if (nrow(noun_adj_pairs_ds) > 0) {
      # Full analysis on downsampled pairs
      analyze_inf_vs_deriv(noun_adj_pairs_ds, "mutation_inflect", "mutation_deriv_adj", "N -> Adj (downsampled)")

      na_tp_ds <- make_deriv_summary(
        noun_adj_pairs_ds, "mutation_deriv_adj",
        "N->Adj derivation, base plural mutated (downsampled)",
        "N->Adj derivation, base plural non-mut. (downsampled)"
      )
      print_full(na_tp_ds)
    }
  }

  cat_section("TOLERANCE PRINCIPLE: SEGMENT-LEVEL PATTERNS (DOWNSAMPLED, TP DOMAIN)")
  nouns_tp_down_single <- nouns_opp_down_single |> filter(tp_in_domain == TRUE)
  # Exclude cluster items so the simple-segment rows don't double-count.
  seg_tp_all_ds <- compute_segment_tp_tables(
    nouns_tp_down_single |> filter(is.na(cluster) | cluster == ""),
    label_suffix = " (downsampled)"
  )
  print_full(seg_tp_all_ds)
  if (any(seg_tp_all_ds$small_n, na.rm = TRUE)) {
    cat("\nSmall-N downsampled cells — Bayesian TP is more reliable for these:\n")
    seg_tp_all_ds |>
      filter(small_n) |>
      select(type, N, mutated, non_mutated, tolerated) |>
      print_full()
  }

  cat_section("TOLERANCE PRINCIPLE: FULL VS DOWNSAMPLED COMPARISON")
  seg_tp_ie_raw_full <- nouns_opp |>
    group_by(stem_final, opportunity) |>
    summarise(
      mutated = sum(mutation, na.rm = TRUE),
      non_mutated = sum(!mutation, na.rm = TRUE), .groups = "drop"
    )

  seg_tp_ie_raw_ds <- nouns_opp_down_single |>
    group_by(stem_final, opportunity) |>
    summarise(
      mutated = sum(mutation, na.rm = TRUE),
      non_mutated = sum(!mutation, na.rm = TRUE), .groups = "drop"
    )

  seg_tp_ie_raw_full |>
    rename(mutated_full = mutated, non_mutated_full = non_mutated) |>
    left_join(
      seg_tp_ie_raw_ds |> rename(mutated_ds = mutated, non_mutated_ds = non_mutated),
      by = c("stem_final", "opportunity")
    ) |>
    mutate(
      N_full    = mutated_full + non_mutated_full,
      N_ds      = mutated_ds + non_mutated_ds,
      rate_full = if_else(N_full > 0, mutated_full / N_full, NA_real_),
      rate_ds   = if_else(N_ds > 0, mutated_ds / N_ds, NA_real_)
    ) |>
    arrange(stem_final, opportunity) |>
    print_full()
}

# =========================================================================
# Tolerance Principle: Clusters
# =========================================================================

cat_section("TOLERANCE PRINCIPLE: CLUSTER PATTERNS (FULL DATA)")
nouns_tp |>
  filter(!is.na(cluster_simple)) |>
  compute_segment_tp_tables(group_var = cluster_simple) |>
  print_full()

cat_section("TOLERANCE PRINCIPLE: CLUSTER PATTERNS (DOWNSAMPLED)")
if (!is.null(nouns_opp_down_single) && nrow(nouns_opp_down_single) > 0L) {
  nouns_opp_down_single |>
    filter(tp_in_domain == TRUE, !is.na(cluster_simple)) |>
    compute_segment_tp_tables(label_suffix = " (downsampled)", group_var = cluster_simple) |>
    print_full()
}

# =========================================================================
# NDEB Counts by Class
# =========================================================================
# -ist analysis: two independent questions
# =========================================================================
#
# Q1 (TRIGGER): Does -ist palatalize the ROOT consonant before it?
#     ci/gi before ist = root dorsal palatalized
#     chi/ghi before ist = root dorsal preserved
#
# Q2 (TARGET): Does the -ist suffix's own st cluster palatalize
#     in the plural?
#     -ist → -iști = suffix st palatalizes
#     -ist → -iste/-isturi = suffix st does NOT palatalize

cat_section("-ist AS TRIGGER: ROOT CONSONANT BEHAVIOR")

ist_words <- nouns |>
  filter(lemma_suffix == "-ist") |>
  mutate(
    pre_ist = str_extract(lemma, ".{1,3}(?=ist$)"),
    # Q1: root dorsal status
    root_dorsal_status = case_when(
      str_detect(pre_ist, "ch$") ~ "preserved_k",
      str_detect(pre_ist, "gh$") ~ "preserved_g",
      str_detect(pre_ist, "c$") ~ "palatalized_c",
      str_detect(pre_ist, "g$") ~ "palatalized_g",
      TRUE ~ "non_dorsal"
    ),
    root_class = case_when(
      root_dorsal_status %in% c("preserved_k", "palatalized_c") ~ "c/k",
      root_dorsal_status %in% c("preserved_g", "palatalized_g") ~ "g",
      TRUE ~ "other"
    ),
    root_palatalized = root_dorsal_status %in% c(
      "palatalized_c", "palatalized_g"
    ),
    # Q2: suffix st palatalization in plural
    suffix_st_palatalizes = case_when(
      is.na(plural) | plural == "" ~ NA,
      str_ends(plural, "i\u0219ti") ~ TRUE,
      str_ends(plural, "iste") ~ FALSE,
      str_ends(plural, "isturi") ~ FALSE,
      TRUE ~ NA
    )
  )

# Q1 results: root dorsal trigger
ist_dorsal <- ist_words |> filter(root_class != "other")
cat("Dorsal-final roots before -ist:", nrow(ist_dorsal), "\n\n")

ist_dorsal |>
  group_by(root_class, root_dorsal_status) |>
  summarise(N = n(), .groups = "drop") |>
  arrange(root_class, root_dorsal_status) |>
  print_full()

if (nrow(ist_dorsal) > 0) {
  cat("\nRoot dorsal palatalization rate before -ist:\n")
  ist_dorsal |>
    group_by(root_class) |>
    summarise(
      N = n(),
      N_palatalized = sum(root_palatalized),
      rate = N_palatalized / N,
      .groups = "drop"
    ) |>
    print_full()
}

cat_section("-ist AS TARGET: SUFFIX st PALATALIZATION IN PLURAL")

ist_target <- ist_words |> filter(!is.na(suffix_st_palatalizes))
cat("Words with classifiable plural:", nrow(ist_target), "\n\n")

ist_target |>
  count(suffix_st_palatalizes) |>
  mutate(
    label = if_else(suffix_st_palatalizes,
      "st \u2192 \u0219ti (palatalizes)",
      "st unchanged (-iste/-isturi)"
    )
  ) |>
  select(label, n) |>
  print_full()

# Cross-tab: root dorsal × suffix palatalization
ist_cross <- ist_words |>
  filter(root_class != "other", !is.na(suffix_st_palatalizes))

if (nrow(ist_cross) > 0) {
  cat("\nCross-tab: root dorsal status × suffix palatalization:\n")
  ist_cross |>
    count(root_dorsal_status, suffix_st_palatalizes) |>
    mutate(
      suffix = if_else(suffix_st_palatalizes,
        "st\u2192\u0219ti", "st unchanged"
      )
    ) |>
    select(root_dorsal_status, suffix, n) |>
    print_full()

  cat("\nFull listing:\n")
  ist_cross |>
    select(
      lemma, plural, root_dorsal_status,
      suffix_st_palatalizes
    ) |>
    arrange(root_dorsal_status, suffix_st_palatalizes, lemma) |>
    print_full()
}

# =========================================================================

cat_section("NDEB BY CLASS (FULL LEXICON AND DOWNSAMPLED)")

ndeb_label <- function(x) {
  case_when(
    x == "gimpe" ~ "gimpe type",
    x == "ochi" ~ "ochi-ochi type",
    x == "paduchi" ~ "paduche-paduchi type",
    TRUE ~ x
  )
}

ndeb_tp_full <- nouns_opp |>
  filter(nde_class %in% ndeb_classes) |>
  group_by(nde_class) |>
  summarise(
    mutated = sum(mutation, na.rm = TRUE),
    non_mutated = sum(!mutation, na.rm = TRUE), .groups = "drop"
  ) |>
  mutate(type = ndeb_label(nde_class)) |>
  tp_table(type, mutated, non_mutated) |>
  mutate(subset = "full")

ndeb_tp_ds <- if (!is.null(nouns_opp_down_single) && nrow(nouns_opp_down_single) > 0L) {
  nouns_opp_down_single |>
    filter(nde_class %in% ndeb_classes) |>
    group_by(nde_class) |>
    summarise(
      mutated = sum(mutation, na.rm = TRUE),
      non_mutated = sum(!mutation, na.rm = TRUE), .groups = "drop"
    ) |>
    mutate(type = ndeb_label(nde_class)) |>
    tp_table(type, mutated, non_mutated) |>
    mutate(subset = "downsampled")
} else {
  tibble()
}

bind_rows(ndeb_tp_full, ndeb_tp_ds) |>
  arrange(subset, type) |>
  print_full()

# =========================================================================
# Bayesian Tolerance Principle
# =========================================================================

cat_section("BAYESIAN TOLERANCE PRINCIPLE")

if (!RUN_BAYESIAN_TP) {
  cat("RUN_BAYESIAN_TP = FALSE; skipping.\n")
} else {
  cat("\nBAYESIAN TP: SEGMENT x OPPORTUNITY (PRODUCTIVE TP DOMAIN)\n")
  tp_bayes_full <- run_bayesian_tp(nouns_tp, subset_label = "TP_domain", seed_value = 123L)

  cat("\nBAYESIAN TP: SEGMENT x OPPORTUNITY (DOWNSAMPLED, TP DOMAIN)\n")
  tp_bayes_ds <- if (!is.null(nouns_opp_down_single) && nrow(nouns_opp_down_single) > 0L) {
    run_bayesian_tp(nouns_opp_down_single |> filter(tp_in_domain == TRUE),
      subset_label = "downsampled_TP", seed_value = 456L
    )
  } else {
    list()
  }

  if (length(tp_bayes_full) > 0 || length(tp_bayes_ds) > 0) {
    bind_rows(
      if (length(tp_bayes_full) > 0) bind_rows(tp_bayes_full) else tibble(),
      if (length(tp_bayes_ds) > 0) bind_rows(tp_bayes_ds) else tibble()
    ) |>
      arrange(segment, opportunity, subset) |>
      as_tibble() |>
      print_full()
  }
}

# =========================================================================
# Segment Class Comparison (Bayesian brms)
# =========================================================================

cat_section("SEGMENT CLASS COMPARISON")

if (!RUN_SEGMENT_CLASS_BRMS) {
  cat("RUN_SEGMENT_CLASS_BRMS = FALSE; skipping.\n")
} else {
  nouns_i_classified <- nouns_tp |>
    filter(opportunity == "i", !is.na(mutation)) |>
    mutate(segment_class = segment_class_factor(stem_final))

  cat("\nSEGMENT CLASS COMPARISON: DORSAL VS CORONAL (I-DOMAIN)\n")
  invisible(capture.output(
    {
      model_class_i <- fit_brms_bernoulli(mutation ~ segment_class,
        nouns_i_classified,
        seed = 123
      )
    },
    type = "output"
  ))

  print(summary(model_class_i))
  draws_i <- as_draws_df(model_class_i)
  beta_i <- draws_i[["b_segment_classcoronal"]]
  or_i <- exp(beta_i)
  or_ci_i <- quantile(or_i, probs = c(0.025, 0.975))
  cat(sprintf("\nP(dorsals > coronals | i-domain) = %.3f\n", mean(beta_i < 0)))
  cat(sprintf(
    "OR_coronal_vs_dorsal (i-domain)  = %.3f [95%% CI: %.3f, %.3f]\n",
    median(or_i), or_ci_i[1], or_ci_i[2]
  ))
  cat("  (OR < 1 => dorsals more likely to palatalize)\n")

  nouns_ie_classified <- nouns_tp |>
    filter(opportunity %in% plural_opportunities, !is.na(mutation)) |>
    mutate(
      segment_class = segment_class_factor(stem_final),
      opportunity = factor(opportunity, levels = plural_opportunities)
    )

  cat("\nSEGMENT CLASS COMPARISON: DORSAL VS CORONAL (I+E DOMAIN)\n")
  invisible(capture.output(
    {
      model_class_ie <- fit_brms_bernoulli(mutation ~ segment_class + opportunity,
        nouns_ie_classified,
        seed = 124
      )
    },
    type = "output"
  ))

  print(summary(model_class_ie))
  draws_ie <- as_draws_df(model_class_ie)
  beta_ie <- draws_ie[["b_segment_classcoronal"]]
  or_ie <- exp(beta_ie)
  or_ci_ie <- quantile(or_ie, probs = c(0.025, 0.975))
  cat(sprintf("\nP(dorsals > coronals | i+e) = %.3f\n", mean(beta_ie < 0)))
  cat(sprintf(
    "OR_coronal_vs_dorsal (i+e)  = %.3f [95%% CI: %.3f, %.3f]\n",
    median(or_ie), or_ci_ie[1], or_ci_ie[2]
  ))
}

# =========================================================================
# EXPORT: TP Summary Tables
# =========================================================================
# Compute with-NDEB comparison rows here, close to where they're used.
# Exclude cluster items here too so the simple-segment rows are comparable
# across the with/without-NDEB versions.
seg_tp_all_with_ndeb <- compute_segment_tp_tables(
  nouns_opp |> filter(is.na(cluster) | cluster == ""),
  label_suffix = " (with NDEB)"
)
seg_tp_all_ds_with_ndeb <- if (!is.null(nouns_opp_down_single) && nrow(nouns_opp_down_single) > 0L) {
  compute_segment_tp_tables(
    nouns_opp_down_single |> filter(is.na(cluster) | cluster == ""),
    label_suffix = " (downsampled, with NDEB)"
  )
} else {
  NULL
}

cat_section("EXPORT: TP SUMMARY TABLES (FULL & DOWNSAMPLED)")

blank_row <- tibble(
  type = "", N = NA_integer_, mutated = NA_integer_,
  `non-mutated` = NA_integer_, rate = NA_real_,
  majority = NA, tolerated = NA, memo = ""
)

header_row <- function(label) {
  tibble(
    type = label, N = NA_integer_, mutated = NA_integer_,
    `non-mutated` = NA_integer_, rate = NA_real_,
    majority = NA, tolerated = NA, memo = ""
  )
}

get_seg_row <- function(tbl, label, new_label = NULL, memo_str = "") {
  out <- tbl |>
    filter(type == label) |>
    transmute(type, N, mutated,
      `non-mutated` = non_mutated,
      rate, majority = majority_mutates, tolerated, memo = memo_str
    )
  if (!is.null(new_label)) out$type <- new_label
  out
}

# ---- Full lexicon ----

ndeb_export_full <- ndeb_tp_full |>
  transmute(type, N, mutated,
    `non-mutated` = non_mutated, rate,
    majority = majority_mutates, tolerated,
    memo = case_when(
      type == "gimpe type" ~ "classic NDEB",
      type == "ochi-ochi type" ~ "NDEB if <i> is root-final; exception otherwise",
      type == "paduche-paduchi type" ~ "NDEB in singular; plural ordering ambiguous",
      TRUE ~ ""
    )
  ) |>
  arrange(type)

dorsals_full <- bind_rows(
  get_seg_row(seg_tp_all, "<c> + <-i> plural"),
  get_seg_row(seg_tp_all, "<c> + <-e> plural"),
  get_seg_row(seg_tp_all, "<c> + <-i, -e> plural"),
  get_seg_row(
    seg_tp_all_with_ndeb, "<c> + <-i, -e> plural (with NDEB)",
    "<c> + <-i, -e> plural, WITH NDEB",
    "Full i/e domain (includes NDEB)"
  ),
  get_seg_row(seg_tp_all, "<g> + <-i> plural"),
  get_seg_row(seg_tp_all, "<g> + <-e> plural"),
  get_seg_row(seg_tp_all, "<g> + <-i, -e> plural"),
  get_seg_row(
    seg_tp_all_with_ndeb, "<g> + <-i, -e> plural (with NDEB)",
    "<g> + <-i, -e> plural, WITH NDEB",
    "Full i/e domain (includes NDEB)"
  )
)

coronals_full <- bind_rows(
  lapply(c("s", "z", "t", "d"), function(seg) {
    bind_rows(
      get_seg_row(seg_tp_all, paste0("<", seg, "> + <-i> plural")),
      get_seg_row(seg_tp_all, paste0("<", seg, "> + <-e> plural")),
      get_seg_row(seg_tp_all, paste0("<", seg, "> + <-i, -e> plural"))
    )
  })
)

cluster_rows <- function(tbl) {
  tbl |>
    filter(str_starts(type, "<ct>|<sc>|<st>")) |>
    transmute(type, N, mutated,
      `non-mutated` = non_mutated,
      rate, majority = majority_mutates, tolerated, memo = ""
    ) |>
    arrange(type)
}

cluster_tp_full <- nouns_tp |>
  filter(!is.na(cluster_simple)) |>
  compute_segment_tp_tables(group_var = cluster_simple)

clusters_full <- cluster_rows(cluster_tp_full)

make_deriv_export <- function(tp_tbl, all_pairs, deriv_col, all_label) {
  if (is.null(tp_tbl) || nrow(tp_tbl) == 0) {
    return(tibble())
  }
  summary_rows <- tp_tbl |>
    transmute(type, N, mutated,
      `non-mutated` = non_mutated,
      rate, majority, tolerated, memo = ""
    )
  if (is.null(all_pairs) || nrow(all_pairs) == 0) {
    return(summary_rows)
  }
  all_row <- all_pairs |>
    summarise(
      mutated = sum(.data[[deriv_col]], na.rm = TRUE),
      non_mutated = sum(!.data[[deriv_col]], na.rm = TRUE)
    ) |>
    mutate(type = all_label) |>
    tp_table(type, mutated, non_mutated) |>
    transmute(type, N, mutated,
      `non-mutated` = non_mutated,
      rate, majority = majority_mutates, tolerated, memo = ""
    )
  bind_rows(summary_rows, all_row)
}

deriv_nv_full <- if (exists("nv_tp_full") && exists("denom_pairs")) {
  make_deriv_export(nv_tp_full, denom_pairs, "mutation_deriv_verb", "N->V derivation, ALL")
} else {
  tibble()
}

deriv_na_full <- if (exists("na_tp_full") && exists("noun_adj_pairs")) {
  make_deriv_export(na_tp_full, noun_adj_pairs, "mutation_deriv_adj", "N->Adj derivation, ALL")
} else {
  tibble()
}

assemble_summary <- function(ndeb_blk, dorsals_blk, coronals_blk,
                             clusters_blk, nv_blk, na_blk) {
  bind_rows(
    header_row("NDEB?"), ndeb_blk, blank_row,
    header_row("Dorsals"), dorsals_blk, blank_row,
    header_row("Coronals"), coronals_blk, blank_row,
    header_row("Clusters"), clusters_blk, blank_row,
    header_row("Derivational (N->V, front suffix)"), nv_blk, blank_row,
    header_row("Derivational (N->Adj, front suffix)"), na_blk
  )
}

summary_full <- assemble_summary(
  ndeb_export_full, dorsals_full, coronals_full,
  clusters_full, deriv_nv_full, deriv_na_full
)
readr::write_csv(
  summary_full |> rename(`majority?` = majority, `tolerated?` = tolerated),
  file.path("analysis", "romanian_tp_summary_full.csv")
)
cat("Wrote full-lexicon TP summary to: analysis/romanian_tp_summary_full.csv\n\n")

# ---- Downsampled lexicon ----

if (!is.null(nouns_opp_down_single) && nrow(nouns_opp_down_single) > 0L) {
  ndeb_export_ds <- ndeb_tp_ds |>
    transmute(type, N, mutated,
      `non-mutated` = non_mutated,
      rate, majority = majority_mutates, tolerated, memo = ""
    ) |>
    arrange(type)

  get_seg_row_ds <- function(label, new_label = NULL, memo_str = "") {
    get_seg_row(seg_tp_all_ds, label, new_label, memo_str)
  }

  dorsals_ds <- bind_rows(
    get_seg_row_ds("<c> + <-i> plural (downsampled)"),
    get_seg_row_ds("<c> + <-e> plural (downsampled)"),
    get_seg_row_ds("<c> + <-i, -e> plural (downsampled)"),
    if (!is.null(seg_tp_all_ds_with_ndeb)) {
      get_seg_row(
        seg_tp_all_ds_with_ndeb,
        "<c> + <-i, -e> plural (downsampled, with NDEB)",
        "<c> + <-i, -e> plural, downsampled, WITH NDEB",
        "Full i/e domain (includes NDEB)"
      )
    },
    get_seg_row_ds("<g> + <-i> plural (downsampled)"),
    get_seg_row_ds("<g> + <-e> plural (downsampled)"),
    get_seg_row_ds("<g> + <-i, -e> plural (downsampled)"),
    if (!is.null(seg_tp_all_ds_with_ndeb)) {
      get_seg_row(
        seg_tp_all_ds_with_ndeb,
        "<g> + <-i, -e> plural (downsampled, with NDEB)",
        "<g> + <-i, -e> plural, downsampled, WITH NDEB",
        "Full i/e domain (includes NDEB)"
      )
    }
  )

  coronals_ds <- bind_rows(
    lapply(c("s", "z", "t", "d"), function(seg) {
      bind_rows(
        get_seg_row_ds(paste0("<", seg, "> + <-i> plural (downsampled)")),
        get_seg_row_ds(paste0("<", seg, "> + <-e> plural (downsampled)")),
        get_seg_row_ds(paste0("<", seg, "> + <-i, -e> plural (downsampled)"))
      )
    })
  )

  cluster_tp_ds <- nouns_opp_down_single |>
    filter(tp_in_domain == TRUE, !is.na(cluster_simple)) |>
    compute_segment_tp_tables(label_suffix = " (downsampled)", group_var = cluster_simple)
  clusters_ds <- cluster_rows(cluster_tp_ds)

  deriv_nv_ds <- if (exists("nv_tp_ds") && exists("denom_pairs_ds")) {
    make_deriv_export(
      nv_tp_ds, denom_pairs_ds, "mutation_deriv_verb",
      "N->V derivation, ALL (downsampled)"
    )
  } else {
    tibble()
  }

  deriv_na_ds <- if (exists("na_tp_ds") && exists("noun_adj_pairs_ds")) {
    make_deriv_export(
      na_tp_ds, noun_adj_pairs_ds, "mutation_deriv_adj",
      "N->Adj derivation, ALL (downsampled)"
    )
  } else {
    tibble()
  }

  summary_ds <- assemble_summary(
    ndeb_export_ds, dorsals_ds, coronals_ds,
    clusters_ds, deriv_nv_ds, deriv_na_ds
  )
  readr::write_csv(
    summary_ds |> rename(`majority?` = majority, `tolerated?` = tolerated),
    file.path("analysis", "romanian_tp_summary_downsampled.csv")
  )
  cat("Wrote downsampled TP summary to: analysis/romanian_tp_summary_downsampled.csv\n\n")
} else {
  cat("No downsampled summary created (no downsampled lexicon available).\n\n")
}

cat_section("ANALYSIS FINISHED")
