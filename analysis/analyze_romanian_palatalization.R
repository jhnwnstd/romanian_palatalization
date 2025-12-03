#!/usr/bin/Rscript

# Romanian Palatalization Analysis
#
# Run styler::style_file("analysis/analyze_romanian_palatalization.R") after edits

# =========================================================================
# Setup
# =========================================================================

required_pkgs <- c(
  "dplyr", "readr", "stringr", "tidyr", "broom",
  "brms", "posterior", "loo"
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
front_verb_suffixes <- c("-i", "-ui") # front-vowel verbalizers (for N→V / Adj→V)
suffix_interest <- c("-ic", "-ist", "-esc", "-ică", "-ice") # denominal/adjectival suffixes we explicitly track

# NDEB = non-derived exception base lemmas (gimpe / ochi / păduche patterns)
# These are lexically special noun bases that behave like exceptions but form structured families.
ndeb_classes <- c("gimpe", "ochi", "paduchi")
ndeb_observable <- c("ochi", "paduchi") # only these are observable as DE in the plural domain; gimpe-type is hidden

plural_opportunities <- c("i", "e") # core i/e opportunity domain
plural_opportunities_all <- c("i", "e", "uri", "none") # all plural types that can show up in the raw data

# =========================================================================
# Analysis Constants
# =========================================================================

MIN_SAMPLE_SIZE_BAYESIAN <- 5L  # avoid running Bayesian TP on tiny cells that give unstable posteriors
SMALL_CELL_THRESHOLD <- 20L     # flag small downsampled cells for interpretability

PRECISION_THETA <- 2L # print precision for θ_N
PRECISION_PROB <- 4L  # print precision for probabilities / rates

CLUSTER_TYPES <- c("st", "sc", "ct") # cluster types we care about for TP-style cluster analysis

# =========================================================================
# Run-mode toggles
# =========================================================================

RUN_BAYESIAN_TP         <- TRUE  # set FALSE to skip Bayesian TP fits
RUN_SEGMENT_CLASS_BRMS  <- TRUE  # set FALSE to skip segment-class brms models
RUN_DERIVATION_ANALYSES <- TRUE  # set FALSE to skip inflection vs derivation analyses

# =========================================================================
# Helper Functions
# =========================================================================

cat_section <- function(title) {
  cat("\n", title, "\n", sep = "")
}

cat_subsection <- function(title) {
  cat("\n", title, "\n", sep = "")
}

print_full <- function(x) {
  print(x, n = Inf, width = Inf)
}

calc_rate <- function(df) {
  df |>
    mutate(
      rate_mut = if_else(.data$N_opp > 0, .data$N_mut / .data$N_opp, NA_real_)
    )
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
      # we only treat a pattern as majority-mutating if there is a clear majority
      majority_mutates = case_when(
        .data$N == 0L ~ NA,
        .data$mutated == .data$non_mutated ~ NA,
        TRUE ~ .data$mutated > .data$non_mutated
      ),
      # "exceptions" are whichever side is in the minority, given the majority direction
      exceptions = case_when(
        is.na(.data$majority_mutates) ~ NA_real_,
        .data$majority_mutates ~ as.numeric(.data$non_mutated),
        !.data$majority_mutates ~ as.numeric(.data$mutated)
      ),
      theta_N = tp_threshold(.data$N),
      tolerated = if_else(
        is.na(.data$exceptions) | is.na(.data$theta_N),
        NA,
        .data$exceptions <= .data$theta_N
      )
    )
}

# NOTE: updated helper so priors no longer reference undefined Stan identifiers like `intercept_sd`
fit_brms_bernoulli <- function(formula,
                               data,
                               seed,
                               intercept_sd = 2.5,
                               b_sd = 2.5,
                               include_b_prior = TRUE) {
  # Build priors as strings (brms evaluates them at Stan generation time)
  intercept_prior_str <- sprintf("normal(0, %s)", format(intercept_sd, scientific = FALSE, trim = TRUE))
  priors <- brms::set_prior(
    intercept_prior_str,
    class = "Intercept"
  )

  if (include_b_prior) {
    b_prior_str <- sprintf("normal(0, %s)", format(b_sd, scientific = FALSE, trim = TRUE))
    priors <- c(
      priors,
      brms::set_prior(
        b_prior_str,
        class = "b"
      )
    )
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

# Approximate palatalization from IPA for derived forms
# This is deliberately coarse: it lets us compare derived forms to inflection even when spelling is opaque.
detect_palatal_from_ipa <- function(stem_final, ipa_str) {
  case_when(
    stem_final == "c" ~ str_detect(ipa_str, "t͡ʃ|tʃ"),
    stem_final == "g" ~ str_detect(ipa_str, "d͡ʒ|dʒ"),
    stem_final == "t" ~ str_detect(ipa_str, "t͡s|ts"),
    stem_final == "d" ~ str_detect(ipa_str, "z|dʲ"),
    stem_final == "s" ~ str_detect(ipa_str, "ʃ"),
    stem_final == "z" ~ str_detect(ipa_str, "ʒ"),
    TRUE ~ NA
  )
}

# TP tables for segments (or clusters) and opportunities
# This keeps all the "how to compute TP" logic in one place so that all segment- and cluster-level
# tables are constructed in a consistent way and can be compared directly.
compute_segment_tp_tables <- function(data, label_suffix = "", group_var = stem_final) {
  group_var_name <- rlang::as_name(rlang::enquo(group_var))

  # By opportunity (i vs e)
  by_opp <- data |>
    group_by({{ group_var }}, opportunity) |>
    summarise(
      mutated     = sum(mutation, na.rm = TRUE),
      non_mutated = sum(!mutation, na.rm = TRUE),
      .groups     = "drop"
    ) |>
    mutate(
      type = paste0("<", .data[[group_var_name]], "> + <-", opportunity, "> plural", label_suffix)
    ) |>
    tp_table(type, mutated, non_mutated)

  # Combined i+e
  combined <- data |>
    group_by({{ group_var }}) |>
    summarise(
      mutated     = sum(mutation, na.rm = TRUE),
      non_mutated = sum(!mutation, na.rm = TRUE),
      .groups     = "drop"
    ) |>
    mutate(
      type = paste0("<", .data[[group_var_name]], "> + <-i, -e> plural", label_suffix)
    ) |>
    tp_table(type, mutated, non_mutated)

  bind_rows(by_opp, combined) |>
    arrange(type)
}

# Bayesian TP: quantify how often the posterior "violates" the classical TP, given sampling uncertainty.
run_bayesian_tp <- function(data, subset_label, seed_value = 123L) {
  if (is.null(data) || nrow(data) == 0L) {
    return(list())
  }

  results <- list()

  for (seg in segments_of_interest) {
    for (opp in plural_opportunities) {
      seg_data <- data |>
        filter(
          stem_final == seg,
          opportunity == opp,
          !is.na(mutation)
        )

      if (nrow(seg_data) < MIN_SAMPLE_SIZE_BAYESIAN) next

      N <- nrow(seg_data)
      theta_N <- tp_threshold(N)
      tolerance_rate <- theta_N / N

      # Diagnostic: verify data before fitting
      sum_mut_check <- sum(seg_data$mutation, na.rm = TRUE)
      sum_non_check <- sum(!seg_data$mutation, na.rm = TRUE)

      cat(sprintf(
        "  Fitting <%s> + <-%s> (%s): N=%d, mutated=%d, non-mutated=%d\n",
        seg, opp, subset_label, N, sum_mut_check, sum_non_check
      ))

      invisible(capture.output(
        {
          # Intercept-only model: we suppress b-priors so brms doesn't create priors for non-existent b’s
          model_tp <- fit_brms_bernoulli(
            mutation ~ 1,
            seg_data,
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
      p_nonmutate <- 1 - p_mutate

      sum_mut <- sum(seg_data$mutation, na.rm = TRUE)
      sum_non <- sum(!seg_data$mutation, na.rm = TRUE)
      majority <- sum_mut > sum_non

      p_exception <- if (majority) p_nonmutate else p_mutate
      prob_intolerable <- mean(p_exception > tolerance_rate)

      results[[paste(seg, opp, subset_label, sep = "_")]] <- data.frame(
        segment = seg,
        opportunity = opp,
        subset = subset_label,
        majority_mutates = majority,
        N = N,
        theta_N = round(theta_N, PRECISION_THETA),
        tolerance_rate = round(tolerance_rate, PRECISION_PROB),
        median_p_mutate = round(median(p_mutate), PRECISION_PROB),
        ci_lower = round(quantile(p_mutate, 0.025), PRECISION_PROB),
        ci_upper = round(quantile(p_mutate, 0.975), PRECISION_PROB),
        prob_exceeds_tol = round(prob_intolerable, PRECISION_PROB)
      )
    }
  }

  results
}

# -------------------------------------------------------------------------
# Analytic helpers to keep later code DRY
# -------------------------------------------------------------------------

# Collapse detailed suffix tags into a small set of theoretically relevant groups,
# so that regression models don't chase micro-suffix idiosyncrasies.
suffix_group_factor <- function(lemma_suffix) {
  group <- dplyr::case_when(
    lemma_suffix == "-ic" ~ "ic",
    lemma_suffix == "-ist" ~ "ist",
    lemma_suffix %in% c("-ică", "-ice") ~ "ica_ice",
    lemma_suffix == "none" ~ "none",
    TRUE ~ "other"
  )
  factor(group, levels = c("none", "ic", "ist", "ica_ice", "other"))
}

# Collapse individual consonants into dorsal vs coronal classes,
# so the segment-class models reflect the theoretical contrast of interest.
segment_class_factor <- function(stem_final) {
  cls <- dplyr::if_else(stem_final %in% c("c", "g"), "dorsal", "coronal")
  factor(cls, levels = c("dorsal", "coronal"))
}

# DRY helper for inflection vs derivation comparisons
# Analyzes agreement patterns between inflectional and derivational palatalization,
# runs McNemar test and logistic regression.
analyze_inf_vs_deriv <- function(df, base_col, deriv_col, label) {
  if (is.null(df) || nrow(df) == 0) {
    cat("\n(No data for ", label, ")\n\n", sep = "")
    return(invisible(NULL))
  }

  # Rename columns to standardized names for easier processing
  df <- df |>
    mutate(
      base_mut = .data[[base_col]],
      deriv_mut = .data[[deriv_col]]
    ) |>
    filter(!is.na(base_mut), !is.na(deriv_mut))

  if (nrow(df) == 0) {
    cat("\n(No valid mutation data for ", label, ")\n\n", sep = "")
    return(invisible(NULL))
  }

  # Agreement patterns
  cat("\nInflection vs derivation agreement patterns (", label, "):\n", sep = "")
  patterns <- df |>
    mutate(
      pattern = dplyr::case_when(
        base_mut & deriv_mut ~ "both_mutate",
        !base_mut & !deriv_mut ~ "both_nonmutate",
        base_mut & !deriv_mut ~ "inflection_only",
        !base_mut & deriv_mut ~ "derivation_only",
        TRUE ~ NA_character_
      )
    ) |>
    count(pattern, sort = TRUE) |>
    mutate(prop = n / sum(n))
  print_full(patterns)

  # McNemar test
  tab <- with(df, table(base_mut, deriv_mut))
  if (all(dim(tab) == c(2L, 2L))) {
    cat("\nMcNemar test (", label, "):\n", sep = "")
    print_full(broom::tidy(stats::mcnemar.test(tab)))
  } else {
    cat("\nContingency table (", label, ") not 2×2; skipping McNemar test.\n", sep = "")
    print(tab)
  }

  # Logistic regression
  cat("\nLogistic regression: does derivational palatalization track inflection? (", label, ")\n", sep = "")
  model <- glm(deriv_mut ~ base_mut, data = df, family = binomial())
  print_full(broom::tidy(model))

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
    mutate(
      type = if_else(
        base_plural_mutates,
        label_yes,
        label_no
      )
    ) |>
    tp_table(type, mutated, non_mutated) |>
    select(
      type,
      N,
      mutated,
      non_mutated,
      rate,
      majority = majority_mutates,
      tolerated
    )
}

build_downsampled_lexica <- function(nouns_opp, sample_sizes = c(1000L, 2500L, 5000L, 10000L)) {
  nouns_opp_freq <- nouns_opp |>
    group_by(lemma) |>
    summarise(
      lemma_freq = max(freq_ron_wikipedia_2021_1M, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      lemma_freq = if_else(is.na(lemma_freq) | !is.finite(lemma_freq), 0, lemma_freq)
    )

  nouns_opp_freq_pos <- filter(nouns_opp_freq, lemma_freq > 0)

  if (nrow(nouns_opp_freq_pos) == 0L) {
    return(list(
      reference = NULL,
      tp_table = NULL,
      freq_table = nouns_opp_freq_pos
    ))
  }

  nouns_opp_with_freq <- nouns_opp |>
    left_join(nouns_opp_freq_pos, by = "lemma")

  seg_tp_ie_ds_list <- vector("list", length(sample_sizes))
  reference <- NULL

  for (i in seq_along(sample_sizes)) {
    n_lex <- sample_sizes[[i]]
    target_n_lex <- min(n_lex, nrow(nouns_opp_freq_pos))

    top_lemmas <- nouns_opp_freq_pos |>
      arrange(desc(lemma_freq)) |>
      slice_head(n = target_n_lex) |>
      pull(lemma)

    nouns_opp_down <- nouns_opp_with_freq |>
      filter(lemma %in% top_lemmas) |>
      arrange(desc(lemma_freq), lemma, plural) |>
      distinct(lemma, .keep_all = TRUE)

    if (is.null(reference) && n_lex == 1000L) {
      reference <- nouns_opp_down
      cat("Reference downsampled lexicon (top 1000 most frequent lemmas):\n")
      cat("  target_lexemes:", target_n_lex, "\n")
      cat("  unique lemmas:", n_distinct(reference$lemma), "\n")
      cat("  rows:", nrow(reference), "\n\n")
    }

    seg_tp_ie_ds_list[[i]] <- nouns_opp_down |>
      group_by(stem_final, opportunity) |>
      summarise(
        mutated     = sum(mutation, na.rm = TRUE),
        non_mutated = sum(!mutation, na.rm = TRUE),
        .groups     = "drop"
      ) |>
      mutate(sample_lexemes = n_lex)
  }

  seg_tp_ie_ds_all <- bind_rows(seg_tp_ie_ds_list) |>
    mutate(
      N = mutated + non_mutated,
      rate = if_else(N > 0, mutated / N, NA_real_),
      majority_mutates = case_when(
        N == 0L ~ NA,
        mutated == non_mutated ~ NA,
        TRUE ~ mutated > non_mutated
      ),
      exceptions = case_when(
        is.na(majority_mutates) ~ NA_real_,
        majority_mutates ~ as.numeric(non_mutated),
        !majority_mutates ~ as.numeric(mutated)
      ),
      theta_N = tp_threshold(N),
      tolerated = if_else(is.na(exceptions) | is.na(theta_N), NA, exceptions <= theta_N)
    )

  list(
    reference = reference,
    tp_table = seg_tp_ie_ds_all,
    freq_table = nouns_opp_freq_pos
  )
}

# =========================================================================
# Data Input
# =========================================================================

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "data/romanian_lexicon_with_freq.csv"
output_log <- file.path("analysis", "romanian_palatalization_analysis.log")

cat("Working directory:", getwd(), "\n")
cat("Log file:", normalizePath(output_log, mustWork = FALSE), "\n")

# Log to file, but keep echoing to console for interactive runs
sink(output_log, split = TRUE)
on.exit(sink(NULL), add = TRUE)

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
  read_csv(
    input_file,
    show_col_types = FALSE,
    # Explicit types for columns that sometimes get parsed as non-character
    col_types = cols(
      derived_adj = col_character(),
      ipa_derived_adj = col_character(),
      ipa_raw_lemma = col_character(),
      ipa_raw_pl = col_character()
    )
  )
)

# Surface any CSV parsing problems early, so later patterns aren't silently skewed.
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
  mutate(
    # Some frequency fields come in as character; normalize so they can be ranked / downsampled.
    freq_ron_wikipedia_2021_1M = as.numeric(freq_ron_wikipedia_2021_1M)
  ) |>
  mutate(
    # Normalize basic categorical tags so downstream filters don't have to worry about stray whitespace / cases.
    across(c(pos, gender, stem_final, cluster, plural), ~ trimws(as.character(.))),
    pos = toupper(pos),
    # Ensure "opportunity" is NA outside the explicitly supported set, so the TP domain is well-defined.
    opportunity = if_else(opportunity %in% plural_opportunities_all, opportunity, NA_character_),
    # Make NDEB and suffix tags explicit "none" rather than NA/blank, so grouping is well-behaved.
    nde_class = if_else(is.na(nde_class) | nde_class == "", "none", nde_class),
    lemma_suffix = if_else(is.na(lemma_suffix) | lemma_suffix == "", "none", lemma_suffix),
    # Coarse-grain plural endings, used to determine whether an item sits in the i/e domain even when
    # the original "opportunity" annotation is absent or "none".
    plural_ending = case_when(
      is.na(plural) | plural == "" ~ "none",
      str_ends(plural, "uri") ~ "uri",
      str_ends(plural, "i") ~ "i",
      str_ends(plural, "e") ~ "e",
      TRUE ~ "other"
    ),
    # TP-specific "opportunity" that pulls NDEB lemmas into the i/e domain whenever their surface plural
    # transparently reflects that environment. This lets us count NDEB as exceptions where they are *visible*.
    opportunity_tp = case_when(
      nde_class %in% ndeb_classes &
        opportunity == "none" &
        plural_ending %in% c("i", "e") ~ plural_ending,

      # NDEB with plural in -uri: treat them as i-type environments, since the triggering front vowel is still present.
      nde_class %in% ndeb_classes &
        opportunity == "none" &
        plural_ending == "uri" ~ "i",

      # everyone else: keep the original opportunity value
      TRUE ~ opportunity
    )
  )

nouns <- filter(lex, pos == "N")

nouns_opp <- nouns |>
  filter(
    opportunity_tp %in% plural_opportunities,
    !is.na(stem_final),
    stem_final %in% segments_of_interest
  ) |>
  mutate(
    # Use TP-opportunity for all downstream grouping (includes NDEB items in the domain where they are theoretically relevant).
    opportunity = opportunity_tp,
    # Keep a simple cluster variable that only tracks the clusters we care about in TP-style cluster analyses.
    cluster_simple = if_else(
      cluster %in% CLUSTER_TYPES,
      cluster,
      NA_character_
    ),
    # High-level exception categories: these are used to see *where* the non-undergoers live.
    exception_category = case_when(
      mutation ~ "undergoes",
      nde_class %in% ndeb_classes ~ paste0("NDEB_", nde_class),
      lemma_suffix %in% suffix_interest ~ paste0("suffix_", lemma_suffix),
      is_true_exception ~ "true_exception",
      TRUE ~ "other_non_undergoer"
    )
  )

# Core i/e-domain, excluding NDEB, to approximate the productive grammar without the structured DE families.
nouns_opp_no_ndeb <- nouns_opp |>
  filter(!(nde_class %in% ndeb_classes))

ndeb_rows <- filter(nouns, nde_class %in% ndeb_classes)

cat_section("BASIC COUNTS")
cat("Total rows:", nrow(lex), "\n")
cat("Nouns:", nrow(nouns), "\n")
cat("Nouns in i/e domain with target segments (incl. NDEB):", nrow(nouns_opp), "\n")
cat("Nouns in i/e domain with target segments (NO NDEB):", nrow(nouns_opp_no_ndeb), "\n\n")

# =========================================================================
# Quality Control
# =========================================================================

cat_section("QC: GENDER ON NOUN ROWS")
# Gender is needed mainly for sanity: errors here can hint at mis-labeled POS or mis-parsed rows.
nouns_missing_gender <- filter(nouns, is.na(gender) | gender == "")
cat("Missing gender:", nrow(nouns_missing_gender), "\n")
if (nrow(nouns_missing_gender) > 0) {
  nouns_missing_gender |>
    select(lemma, gloss, gender, source, notes) |>
    head(10) |>
    print()
}

cat_section("QC: MUTATION VS OPPORTUNITY")
# Any mutated item outside the i/e opportunity domain is suspect and usually indicates inconsistent annotation.
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

# Ochi/păduche types are observable as DE exceptions; gimpe-type are hidden and should not inflate exception counts.
cat("\nNDEB nouns of ochi/păduche type (observable DE exceptions):", nrow(ndeb_ochi_pad), "\n")
if (nrow(ndeb_ochi_pad) > 0) {
  ndeb_ochi_pad |>
    select(lemma, plural, stem_final, nde_class, opportunity) |>
    arrange(nde_class, stem_final, lemma) |>
    print_full()
}

cat("\nNDEB nouns of gimpe type (unobservable as DE; excluded from exception counts):", nrow(ndeb_gimpe), "\n")

cat_section("QC: OPPORTUNITY VS PLURAL ENDING")
# This table is a quick way to see whether the i/e/uri tags line up with actual surface plurals.
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
# Make sure the marked "triggering" suffixes and "target_is_suffix" are internally consistent.
suffix_rows <- filter(nouns, lemma_suffix %in% suffix_interest)
cat("Nouns with tracked suffixes:", nrow(suffix_rows), "\n")
count(suffix_rows, lemma_suffix) |> print()

suffix_flag_true_nomut <- filter(suffix_rows, suffix_triggers_plural_mutation, !mutation)
cat("\nSuffix marked as trigger but lemma does not mutate:", nrow(suffix_flag_true_nomut), "\n")

suffix_flag_true_notarget <- filter(suffix_rows, suffix_triggers_plural_mutation, !target_is_suffix)
cat("\nSuffix marked as trigger but not marked as target site:", nrow(suffix_flag_true_notarget), "\n")

cat_section("QC: PALATAL CONSONANT IN PLURAL")
# palatal_consonant_pl is used as a more fine-grained sanity check on the binary mutation flag.
palatal_nouns <- filter(nouns, !is.na(palatal_consonant_pl), palatal_consonant_pl != "")
cat("Nouns with palatal_consonant_pl populated:", nrow(palatal_nouns), "\n")

palatal_summary <- count(palatal_nouns, stem_final, palatal_consonant_pl, sort = TRUE)
cat("\nPalatal consonant distribution by stem_final:\n")
print_full(palatal_summary)

palatal_no_mutation <- filter(nouns, !is.na(palatal_consonant_pl), palatal_consonant_pl != "", !mutation)
cat("\nNouns with palatal_consonant_pl but mutation=FALSE:", nrow(palatal_no_mutation), "\n")
if (nrow(palatal_no_mutation) > 0) {
  palatal_no_mutation |>
    select(lemma, plural, stem_final, palatal_consonant_pl, mutation, opportunity) |>
    head(10) |>
    print_full()
}

mutation_no_palatal <- filter(
  nouns,
  mutation,
  opportunity %in% plural_opportunities,
  is.na(palatal_consonant_pl) | palatal_consonant_pl == ""
)
cat("\nNouns with mutation=TRUE in i/e domain but missing palatal_consonant_pl:", nrow(mutation_no_palatal), "\n")
if (nrow(mutation_no_palatal) > 0) {
  mutation_no_palatal |>
    select(lemma, plural, stem_final, opportunity, mutation, palatal_consonant_pl) |>
    head(10) |>
    print_full()
}

cat_section("QC: DUPLICATE LEMMAS")
# Duplicated (lemma, pos) pairs hint at double entries from different sources; they matter for lemma-based counts.
lemma_dups <- count(lex, lemma, pos, sort = TRUE) |> filter(n > 1)
cat("Duplicate (lemma, pos) pairs:", nrow(lemma_dups), "\n")
if (nrow(lemma_dups) > 0) head(lemma_dups, 25) |> print()

# =========================================================================
# Inflection vs Derivation
# =========================================================================

cat_section("INFLECTION VS DERIVATION: LEMMA-BASED PATTERNS (NOUNS & ADJECTIVES)")

has_verb_deriv_cols <- all(c("derived_verbs", "ipa_derived_verbs", "deriv_suffixes") %in% names(lex))
has_adj_deriv_cols <- all(c("derived_adj", "ipa_derived_adj") %in% names(lex))

if (!RUN_DERIVATION_ANALYSES) {
  cat("RUN_DERIVATION_ANALYSES = FALSE; skipping all inflection/derivation checks.\n")
} else if (!has_verb_deriv_cols && !has_adj_deriv_cols) {
  cat("No derivational columns present; skipping all inflection/derivation checks.\n")
} else {
  # Focus on lemmas that participate in the i/e domain and have a clear mutation flag.
  noun_base_inflect <- nouns |>
    filter(
      opportunity %in% plural_opportunities,
      !is.na(stem_final),
      stem_final %in% segments_of_interest,
      !is.na(mutation)
    ) |>
    mutate(mutation_inflect = as.logical(mutation))

  # -----------------------------------------------------------------------
  # (1) Noun lemmas: inflectional plurals vs denominal verbs
  # -----------------------------------------------------------------------
  # This probes whether N→V derivations respect the same palatalization pattern
  # as the plural, i.e. whether derivation "copies" inflection.
  cat_subsection("(1) NOUN LEMMAS: INFLECTIONAL PLURALS VS DENOMINAL VERBS")

  if (has_verb_deriv_cols) {
    denom_pairs <- noun_base_inflect |>
      filter(
        !is.na(derived_verbs), derived_verbs != "",
        !is.na(ipa_derived_verbs), ipa_derived_verbs != ""
      ) |>
      mutate(
        verb_suffix_front = deriv_suffixes %in% front_verb_suffixes,
        mutation_deriv_verb = detect_palatal_from_ipa(stem_final, ipa_derived_verbs)
      ) |>
      # Discard cases where the IPA is too underspecified to tell what happened.
      filter(!is.na(mutation_deriv_verb)) |>
      arrange(lemma) |>
      distinct(lemma, .keep_all = TRUE)

    cat("Denominal N–V lemmas in i/e domain with usable IPA:", nrow(denom_pairs), "\n")
    cat("  (of which with front-vowel verbal suffix -i/-ui:", sum(denom_pairs$verb_suffix_front, na.rm = TRUE), ")\n")

    if (nrow(denom_pairs) > 0) {
      # Use DRY helper for inflection vs derivation analysis
      analyze_inf_vs_deriv(denom_pairs, "mutation_inflect", "mutation_deriv_verb", "N → V")

      cat("\nDERIVATIONAL SUMMARY TABLE (N → V, FULL LEXICON; GOOGLE SHEET FORMAT)\n")
      nv_tp_full <- make_deriv_summary(
        denom_pairs,
        "mutation_deriv_verb",
        "N→V derivation, base plural mutated",
        "N→V derivation, base plural non-mut."
      )
      print_full(nv_tp_full)
    }
  } else {
    cat("No denominal verb columns found; skipping N → V comparison.\n")
  }

  # -----------------------------------------------------------------------
  # (2) Noun lemmas: inflectional plurals vs denominal adjectives
  # -----------------------------------------------------------------------
  # Same logic for N→Adj: do derived adjectives behave like the plural?
  cat_subsection("(2) NOUN LEMMAS: INFLECTIONAL PLURALS VS DENOMINAL ADJECTIVES")

  if (has_adj_deriv_cols) {
    noun_adj_pairs <- noun_base_inflect |>
      filter(
        !is.na(derived_adj), derived_adj != "",
        !is.na(ipa_derived_adj), ipa_derived_adj != ""
      ) |>
      mutate(mutation_deriv_adj = detect_palatal_from_ipa(stem_final, ipa_derived_adj)) |>
      filter(!is.na(mutation_deriv_adj)) |>
      arrange(lemma) |>
      distinct(lemma, .keep_all = TRUE)

    cat("Denominal N–Adj lemmas in i/e domain with usable IPA:", nrow(noun_adj_pairs), "\n")

    if (nrow(noun_adj_pairs) > 0) {
      # Use DRY helper for inflection vs derivation analysis
      analyze_inf_vs_deriv(noun_adj_pairs, "mutation_inflect", "mutation_deriv_adj", "N → Adj")

      cat("\nDERIVATIONAL SUMMARY TABLE (N → Adj, FULL LEXICON; GOOGLE SHEET FORMAT)\n")
      na_tp_full <- make_deriv_summary(
        noun_adj_pairs,
        "mutation_deriv_adj",
        "N→Adj derivation, base plural mutated",
        "N→Adj derivation, base plural non-mut."
      )
      print_full(na_tp_full)
    }
  } else {
    cat("No denominal adjective columns found; skipping N → Adj comparison.\n")
  }

  # -----------------------------------------------------------------------
  # (3) Adjective lemmas: inflectional plurals vs derivations
  # -----------------------------------------------------------------------
  # This mirrors the noun analysis, but for adjectives as bases.
  cat_subsection("(3) ADJECTIVE LEMMAS: INFLECTIONAL PLURALS VS DERIVATIONS")

  adj_base_inflect <- lex |>
    filter(
      pos == "ADJ",
      opportunity %in% plural_opportunities,
      !is.na(stem_final),
      stem_final %in% segments_of_interest,
      !is.na(mutation)
    ) |>
    mutate(mutation_inflect = as.logical(mutation))

  cat("Adjective lemmas with i/e plural & target segments:", nrow(adj_base_inflect), "\n")

  # Adjectives → Verbs
  if (has_verb_deriv_cols && nrow(adj_base_inflect) > 0) {
    adj_verb_pairs <- adj_base_inflect |>
      filter(
        !is.na(derived_verbs), derived_verbs != "",
        !is.na(ipa_derived_verbs), ipa_derived_verbs != ""
      ) |>
      mutate(
        verb_suffix_front = deriv_suffixes %in% front_verb_suffixes,
        mutation_deriv_verb = detect_palatal_from_ipa(stem_final, ipa_derived_verbs)
      ) |>
      filter(!is.na(mutation_deriv_verb)) |>
      arrange(lemma) |>
      distinct(lemma, .keep_all = TRUE)

    cat("Adjective lemmas with both plural + derived verb:", nrow(adj_verb_pairs), "\n")

    if (nrow(adj_verb_pairs) > 0) {
      # Use DRY helper for inflection vs derivation analysis
      analyze_inf_vs_deriv(adj_verb_pairs, "mutation_inflect", "mutation_deriv_verb", "Adj → V")
    }
  } else {
    cat("No adjective lemmas with both plural and derived verbs; skipping Adj → V comparison.\n")
  }

  # Adjectives → Adjectives
  if (has_adj_deriv_cols && nrow(adj_base_inflect) > 0) {
    adj_adj_pairs <- adj_base_inflect |>
      filter(
        !is.na(derived_adj), derived_adj != "",
        !is.na(ipa_derived_adj), ipa_derived_adj != "",
        derived_adj != lemma
      ) |>
      mutate(mutation_deriv_adj = detect_palatal_from_ipa(stem_final, ipa_derived_adj)) |>
      filter(!is.na(mutation_deriv_adj)) |>
      arrange(lemma) |>
      distinct(lemma, .keep_all = TRUE)

    cat("Adjective lemmas with plural + non-trivial derived Adj:", nrow(adj_adj_pairs), "\n")

    if (nrow(adj_adj_pairs) > 0) {
      # Use DRY helper for inflection vs derivation analysis
      analyze_inf_vs_deriv(adj_adj_pairs, "mutation_inflect", "mutation_deriv_adj", "Adj → Adj")
    }
  } else {
    cat("No adjective lemmas with both plural and non-trivial derived adjectives; skipping Adj → Adj comparison.\n")
  }
}

# =========================================================================
# Descriptive Summaries
# =========================================================================

cat_section("SEGMENT-WISE MUTATION RATES (I+E COMBINED, NO NDEB)")
seg_summary <- nouns_opp_no_ndeb |>
  group_by(stem_final) |>
  summarise(N_opp = n(), N_mut = sum(mutation, na.rm = TRUE), .groups = "drop") |>
  calc_rate() |>
  arrange(stem_final)
print_full(seg_summary)

cat_section("MUTATION RATES BY SEGMENT AND PLURAL TYPE (I VS E, NO NDEB)")
seg_by_opp <- nouns_opp_no_ndeb |>
  group_by(stem_final, opportunity) |>
  summarise(N_opp = n(), N_mut = sum(mutation, na.rm = TRUE), .groups = "drop") |>
  calc_rate() |>
  arrange(stem_final, opportunity)
print_full(seg_by_opp)

cat_section("CLUSTER INVENTORY IN I/E DOMAIN (NO NDEB)")
cluster_inventory <- nouns_opp_no_ndeb |>
  filter(!is.na(cluster), cluster != "") |>
  count(stem_final, cluster, sort = TRUE)
print_full(cluster_inventory)

cat_section("CLUSTER EFFECTS ON MUTATION (NO NDEB)")
# Here we separate out a small set of clusters that may modulate palatalization
# (st/sc/ct and orthographic chi/che/ghi/ghe) to see whether they depress/enhance rates.
cluster_summary <- nouns_opp_no_ndeb |>
  mutate(
    cluster_type = if_else(
      cluster %in% c("st", "sc", "ct", "chi", "che", "ghi", "ghe"),
      cluster,
      "none"
    )
  ) |>
  group_by(stem_final, cluster_type, opportunity) |>
  summarise(N_opp = n(), N_mut = sum(mutation, na.rm = TRUE), .groups = "drop") |>
  calc_rate() |>
  arrange(cluster_type, opportunity, stem_final)
print_full(cluster_summary)

cat_section("STRUCTURE OF NON-UNDERGOERS (I/E DOMAIN)")
# This shows how the non-undergoers are distributed across suffixal classes, NDEB types, and true exceptions.
non_under_summary <- nouns_opp |>
  filter(!mutation) |>
  group_by(stem_final, exception_category) |>
  summarise(N = n(), .groups = "drop") |>
  group_by(stem_final) |>
  mutate(total_non_under = sum(N), prop = N / total_non_under) |>
  arrange(exception_category, stem_final, desc(prop))
print_full(non_under_summary)

cat_section("NDE DISTRIBUTION BY SEGMENT")
# NDEB statistics are useful to see which segments are most affected by lexically idiosyncratic patterns.
nde_by_seg <- ndeb_rows |>
  group_by(nde_class, stem_final) |>
  summarise(N = n(), .groups = "drop") |>
  arrange(nde_class, stem_final)
print_full(nde_by_seg)

cat_section("PALATALIZATION RATES FOR <-că, -gă> NOUNS")
# These are classic textbook examples of "c/g before -ă" patterns; we break them out separately.
cg_cag_summary <- nouns |>
  filter(
    pos == "N",
    str_ends(lemma, "că") | str_ends(lemma, "gă"),
    opportunity %in% plural_opportunities,
    stem_final %in% c("c", "g")
  ) |>
  group_by(stem_final, opportunity) |>
  summarise(
    N_opp = n(),
    N_mut = sum(mutation, na.rm = TRUE),
    rate_mut = if_else(N_opp > 0, N_mut / N_opp, NA_real_),
    .groups = "drop"
  )
print_full(cg_cag_summary)

cat_section("SUFFIX PATTERNS (ALL TRACKED SUFFIXES)")
suffix_summary_all <- nouns |>
  filter(lemma_suffix %in% suffix_interest) |>
  group_by(lemma_suffix, stem_final, opportunity) |>
  summarise(
    N = n(),
    N_mut = sum(mutation, na.rm = TRUE),
    rate_mut = if_else(N > 0, N_mut / N, NA_real_),
    .groups = "drop"
  ) |>
  arrange(lemma_suffix, opportunity, stem_final)
print_full(suffix_summary_all)

cat_section("SUFFIX PATTERNS WHERE SUFFIX IS TARGET")
suffix_target_summary <- nouns |>
  filter(lemma_suffix %in% suffix_interest, target_is_suffix) |>
  group_by(lemma_suffix, stem_final, opportunity) |>
  summarise(
    N = n(),
    N_mut = sum(mutation, na.rm = TRUE),
    rate_mut = if_else(N > 0, N_mut / N, NA_real_),
    .groups = "drop"
  ) |>
  arrange(lemma_suffix, opportunity, stem_final)
print_full(suffix_target_summary)

cat_section("COMPARISON: 'HAS SUFFIX' VS 'SUFFIX IS TARGET'")
suffix_diff <- anti_join(
  suffix_summary_all,
  suffix_target_summary,
  by = c("lemma_suffix", "stem_final", "opportunity")
)
if (nrow(suffix_diff) > 0) {
  cat("Suffix present but not annotated as target in these cells:\n")
  print_full(suffix_diff)
} else {
  cat("All tracked suffix rows are also marked as targets.\n")
}

cat_section("TRUE EXCEPTIONS IN I/E DOMAIN (NON-NDEB)")
true_exc <- nouns |>
  filter(
    opportunity %in% plural_opportunities,
    !mutation,
    !(nde_class %in% ndeb_classes)
  )
cat("Non-mutating, non-NDEB nouns:", nrow(true_exc), "\n")

true_exc_by_seg <- true_exc |>
  group_by(stem_final) |>
  summarise(N_true_exc = n(), .groups = "drop") |>
  arrange(stem_final)
print_full(true_exc_by_seg)

cat_section("NDEB EXCEPTIONS OF OCHI/PĂDUCHE TYPE (OUTSIDE ALIGNMENT-BASED I/E OPPORTUNITY)")
ndeb_exc_ochi_pad <- nouns |>
  filter(
    nde_class %in% ndeb_observable,
    !mutation
  )

cat("Non-mutating NDEB nouns of ochi/păduche type:", nrow(ndeb_exc_ochi_pad), "\n")

if (nrow(ndeb_exc_ochi_pad) > 0) {
  ndeb_exc_ochi_pad |>
    select(lemma, plural, stem_final, opportunity, nde_class, lemma_suffix, notes) |>
    arrange(nde_class, stem_final, lemma) |>
    print_full()
}

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

# =========================================================================
# Frequency-Based Downsampling
# =========================================================================

cat("DOWNSAMPLING (FREQUENCY-BASED) FOR TOLERANCE PRINCIPLE ANALYSIS\n")

# The TP is sensitive to overall N. Here we rebuild the grammar from the top N most
# frequent lemmas to approximate a learner exposed primarily to high-frequency items.
sample_lexeme_sizes <- c(1000L, 2500L, 5000L, 10000L)

downsampled <- build_downsampled_lexica(nouns_opp, sample_lexeme_sizes)
nouns_opp_down_single <- downsampled$reference
seg_tp_ie_ds_all <- downsampled$tp_table
nouns_opp_freq_pos <- downsampled$freq_table

cat("Unique lemmas in i/e domain (all):", n_distinct(nouns_opp$lemma), "\n")
cat("Unique lemmas with freq > 0:", nrow(nouns_opp_freq_pos), "\n\n")

if (!is.null(seg_tp_ie_ds_all)) {
  cat("Segment × opportunity mutation / TP summary across frequency-filtered lexicons (top N most frequent):\n")
  print_full(seg_tp_ie_ds_all)
  cat("\n")
} else {
  cat("No positive-frequency lemmas in i/e domain; skipping frequency downsampling.\n\n")
}

# =========================================================================
# Tolerance Principle: Full Data
# =========================================================================

cat_section("TOLERANCE PRINCIPLE: SEGMENT-LEVEL PATTERNS (FULL DATA)")

seg_tp_all <- compute_segment_tp_tables(nouns_opp)
print_full(seg_tp_all)

cat_section("TOLERANCE PRINCIPLE: SEGMENT-LEVEL PATTERNS (FULL DATA, NO NDEB)")

seg_tp_all_no_ndeb <- compute_segment_tp_tables(nouns_opp_no_ndeb, label_suffix = " (no NDEB)")
print_full(seg_tp_all_no_ndeb)

# =========================================================================
# NDEB CONTRIBUTION PER TYPE (FULL LEXICON & DOWNSAMPLED)
# =========================================================================

cat_section("TOLERANCE PRINCIPLE: NDEB CONTRIBUTION PER TYPE (FULL LEXICON)")

nouns_opp_ndeb <- nouns_opp |>
  filter(nde_class %in% ndeb_classes)

seg_tp_all_ndeb <- compute_segment_tp_tables(nouns_opp_ndeb, label_suffix = " (NDEB)")
print_full(seg_tp_all_ndeb)

cat_section("TOLERANCE PRINCIPLE: NDEB CONTRIBUTION PER TYPE (REFERENCE DOWNSAMPLED LEXICON)")

if (!is.null(nouns_opp_down_single) && nrow(nouns_opp_down_single) > 0L) {
  nouns_opp_down_ndeb <- nouns_opp_down_single |>
    filter(nde_class %in% ndeb_classes)

  seg_tp_all_ds_ndeb <- compute_segment_tp_tables(nouns_opp_down_ndeb, label_suffix = " (downsampled, NDEB)")
  print_full(seg_tp_all_ds_ndeb)
} else {
  cat("No reference downsampled lexicon available; skipping NDEB-by-type counts (downsampled).\n")
}

# =========================================================================
# Tolerance Principle: Downsampled Data
# =========================================================================

if (!is.null(nouns_opp_down_single) && nrow(nouns_opp_down_single) > 0L) {
  cat_section("DERIVATIONAL SUMMARY TABLES (REFERENCE DOWNSAMPLED LEXICON: TOP 1000 MOST FREQUENT LEMMAS)")

  lemmas_ds <- unique(nouns_opp_down_single$lemma)

  if (has_verb_deriv_cols && RUN_DERIVATION_ANALYSES) {
    denom_pairs_ds <- noun_base_inflect |>
      filter(
        lemma %in% lemmas_ds,
        !is.na(derived_verbs), derived_verbs != "",
        !is.na(ipa_derived_verbs), ipa_derived_verbs != ""
      ) |>
      mutate(
        verb_suffix_front = deriv_suffixes %in% front_verb_suffixes,
        mutation_deriv_verb = detect_palatal_from_ipa(stem_final, ipa_derived_verbs)
      ) |>
      filter(!is.na(mutation_deriv_verb)) |>
      arrange(lemma) |>
      distinct(lemma, .keep_all = TRUE)

    cat("  Denominal N–V lemmas in reference downsampled lexicon:", nrow(denom_pairs_ds), "\n")

    if (nrow(denom_pairs_ds) > 0) {
      nv_tp_ds <- make_deriv_summary(
        denom_pairs_ds,
        "mutation_deriv_verb",
        "N→V derivation, base plural mutated (downsampled)",
        "N→V derivation, base plural non-mut. (downsampled)"
      )
      print_full(nv_tp_ds)
    }
  }

  if (has_adj_deriv_cols && RUN_DERIVATION_ANALYSES) {
    noun_adj_pairs_ds <- noun_base_inflect |>
      filter(
        lemma %in% lemmas_ds,
        !is.na(derived_adj), derived_adj != "",
        !is.na(ipa_derived_adj), ipa_derived_adj != ""
      ) |>
      mutate(mutation_deriv_adj = detect_palatal_from_ipa(stem_final, ipa_derived_adj)) |>
      filter(!is.na(mutation_deriv_adj)) |>
      arrange(lemma) |>
      distinct(lemma, .keep_all = TRUE)

    cat("  Denominal N–Adj lemmas in reference downsampled lexicon:", nrow(noun_adj_pairs_ds), "\n")

    if (nrow(noun_adj_pairs_ds) > 0) {
      na_tp_ds <- make_deriv_summary(
        noun_adj_pairs_ds,
        "mutation_deriv_adj",
        "N→Adj derivation, base plural mutated (downsampled)",
        "N→Adj derivation, base plural non-mut. (downsampled)"
      )
      print_full(na_tp_ds)
    }
  }

  cat_section("TOLERANCE PRINCIPLE: SEGMENT-LEVEL PATTERNS (REFERENCE DOWNSAMPLED LEXICON)")

  seg_tp_all_ds <- compute_segment_tp_tables(nouns_opp_down_single, label_suffix = " (downsampled)")
  print_full(seg_tp_all_ds)

  cat_section("TOLERANCE PRINCIPLE: SEGMENT-LEVEL PATTERNS (REFERENCE DOWNSAMPLED, NO NDEB)")

  nouns_opp_down_single_no_ndeb <- nouns_opp_down_single |>
    filter(!(nde_class %in% ndeb_classes))

  seg_tp_all_ds_no_ndeb <- compute_segment_tp_tables(nouns_opp_down_single_no_ndeb, label_suffix = " (downsampled, no NDEB)")
  print_full(seg_tp_all_ds_no_ndeb)

  cat_section("TOLERANCE PRINCIPLE: SEGMENT-LEVEL COMPARISON (I VS E; FULL VS DOWNSAMPLED)")

  # This lets us see whether the downsampled grammar "looks like" the full one.
  seg_tp_ie_raw_full <- nouns_opp |>
    group_by(stem_final, opportunity) |>
    summarise(mutated = sum(mutation, na.rm = TRUE), non_mutated = sum(!mutation, na.rm = TRUE), .groups = "drop")

  seg_tp_ie_raw_ds <- nouns_opp_down_single |>
    group_by(stem_final, opportunity) |>
    summarise(mutated = sum(mutation, na.rm = TRUE), non_mutated = sum(!mutation, na.rm = TRUE), .groups = "drop")

  seg_tp_ie_compare <- seg_tp_ie_raw_full |>
    select(stem_final, opportunity, mutated, non_mutated) |>
    rename(mutated_full = mutated, non_mutated_full = non_mutated) |>
    left_join(
      seg_tp_ie_raw_ds |>
        select(stem_final, opportunity, mutated, non_mutated) |>
        rename(mutated_ds = mutated, non_mutated_ds = non_mutated),
      by = c("stem_final", "opportunity")
    ) |>
    mutate(
      N_full = mutated_full + non_mutated_full,
      N_ds = mutated_ds + non_mutated_ds,
      rate_full = if_else(N_full > 0, mutated_full / N_full, NA_real_),
      rate_ds = if_else(N_ds > 0, mutated_ds / N_ds, NA_real_)
    ) |>
    arrange(stem_final, opportunity)
  print_full(seg_tp_ie_compare)

  small_ds <- filter(seg_tp_ie_compare, !is.na(N_ds), N_ds < SMALL_CELL_THRESHOLD)
  if (nrow(small_ds) > 0) {
    cat(sprintf("\nCells with N_ds < %d in segment × opportunity (downsampled):\n", SMALL_CELL_THRESHOLD))
    select(small_ds, stem_final, opportunity, N_ds, rate_ds) |> print()
  }
} else {
  cat("\nNo reference downsampled lexicon; skipping single-sample TP comparison.\n")
}

# =========================================================================
# Tolerance Principle: Cluster Patterns
# =========================================================================

cat_section("TOLERANCE PRINCIPLE: CLUSTER PATTERNS (FULL DATA)")

nouns_opp_clusters <- nouns_opp |>
  filter(!is.na(cluster_simple))

cluster_tp_all <- compute_segment_tp_tables(nouns_opp_clusters, group_var = cluster_simple)
print_full(cluster_tp_all)

cat_section("TOLERANCE PRINCIPLE: CLUSTER PATTERNS (REFERENCE DOWNSAMPLED LEXICON)")

if (!is.null(nouns_opp_down_single) && nrow(nouns_opp_down_single) > 0L) {
  nouns_opp_down_clusters <- nouns_opp_down_single |>
    filter(!is.na(cluster_simple))

  cluster_tp_all_ds <- compute_segment_tp_tables(nouns_opp_down_clusters, label_suffix = " (downsampled)", group_var = cluster_simple)
  print_full(cluster_tp_all_ds)
} else {
  cat("No reference downsampled lexicon available; skipping cluster TP (downsampled).\n")
}

# =========================================================================
# Tolerance Principle: NDEB Counts
# =========================================================================

cat_section("NDEB BY CLASS (FULL LEXICON AND REFERENCE DOWNSAMPLED)")

# These labels match the terminology used in the write-up and spreadsheets.
ndeb_label <- function(x) {
  dplyr::case_when(
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
    mutated     = sum(mutation, na.rm = TRUE),
    non_mutated = sum(!mutation, na.rm = TRUE),
    .groups     = "drop"
  ) |>
  mutate(type = ndeb_label(nde_class)) |>
  tp_table(type, mutated, non_mutated) |>
  mutate(subset = "full")

ndeb_tp_ds <- if (!is.null(nouns_opp_down_single) &&
  nrow(nouns_opp_down_single) > 0L) {
  nouns_opp_down_single |>
    filter(nde_class %in% ndeb_classes) |>
    group_by(nde_class) |>
    summarise(
      mutated     = sum(mutation, na.rm = TRUE),
      non_mutated = sum(!mutation, na.rm = TRUE),
      .groups     = "drop"
    ) |>
    mutate(type = ndeb_label(nde_class)) |>
    tp_table(type, mutated, non_mutated) |>
    mutate(subset = "downsampled")
} else {
  tibble() # empty, safe for bind_rows()
}

ndeb_tp_all <- bind_rows(ndeb_tp_full, ndeb_tp_ds) |>
  arrange(subset, type)

print_full(ndeb_tp_all)

# =========================================================================
# Bayesian Tolerance Principle
# =========================================================================

cat_section("BAYESIAN TOLERANCE PRINCIPLE")

if (!RUN_BAYESIAN_TP) {
  cat("RUN_BAYESIAN_TP = FALSE; skipping Bayesian TP analysis.\n")
} else {
  cat("\nBAYESIAN TOLERANCE PRINCIPLE: SEGMENT × OPPORTUNITY (FULL DATA, NO NDEB)\n")
  tolerance_bayesian_results <- run_bayesian_tp(nouns_opp_no_ndeb, subset_label = "full_noNDEB", seed_value = 123L)

  cat("\nBAYESIAN TOLERANCE PRINCIPLE: SEGMENT × OPPORTUNITY (REFERENCE DOWNSAMPLED, NO NDEB)\n")
  if (!is.null(nouns_opp_down_single) && nrow(nouns_opp_down_single) > 0L) {
    # Recompute NDEB-filtered downsampled subset locally to keep dependencies explicit
    nouns_opp_down_single_no_ndeb_for_bayes <- nouns_opp_down_single |>
      filter(!(nde_class %in% ndeb_classes))

    tolerance_bayesian_results_ds <- run_bayesian_tp(
      nouns_opp_down_single_no_ndeb_for_bayes,
      subset_label = "downsampled_noNDEB",
      seed_value = 456L
    )
  } else {
    tolerance_bayesian_results_ds <- list()
    cat("No reference downsampled lexicon available; skipping Bayesian TP for downsampled subset.\n")
  }

  if (length(tolerance_bayesian_results) > 0 || length(tolerance_bayesian_results_ds) > 0) {
    tolerance_bayesian_df <- bind_rows(
      if (length(tolerance_bayesian_results) > 0) bind_rows(tolerance_bayesian_results) else tibble(),
      if (length(tolerance_bayesian_results_ds) > 0) bind_rows(tolerance_bayesian_results_ds) else tibble()
    ) |>
      arrange(segment, opportunity, subset)

    cat("\nBayesian TP results (full vs reference downsampled):\n")
    print_full(as_tibble(tolerance_bayesian_df))
  } else {
    cat("\nInsufficient data for Bayesian TP analysis (both subsets)\n\n")
  }
}

# =========================================================================
# Segment Class Comparison
# =========================================================================

cat_section("SEGMENT CLASS COMPARISON")

if (!RUN_SEGMENT_CLASS_BRMS) {
  cat("RUN_SEGMENT_CLASS_BRMS = FALSE; skipping segment-class brms models.\n")
} else {
  cat("\nSEGMENT CLASS COMPARISON: DORSAL VS CORONAL (I-DOMAIN, NO NDEB)\n")

  nouns_i_classified <- nouns_opp |>
    filter(
      opportunity == "i",
      !is.na(mutation),
      !(nde_class %in% ndeb_classes) # Exclude NDEB so the class contrast reflects the productive grammar
    ) |>
    mutate(
      segment_class = segment_class_factor(stem_final),
      suffix_group = suffix_group_factor(lemma_suffix)
    )

  invisible(capture.output(
    {
      model_class_i <- fit_brms_bernoulli(
        mutation ~ segment_class + suffix_group,
        nouns_i_classified,
        seed = 123
      )
    },
    type = "output"
  ))

  cat("\nSUMMARY: SEGMENT CLASS MODEL (I-DOMAIN ONLY)\n")
  print(summary(model_class_i))

  draws_i <- as_draws_df(model_class_i)
  beta_seg_i <- draws_i[["b_segment_classcoronal"]] # log-odds(coronal) - log-odds(dorsal)
  prob_dorsals_gt_coronals_i <- mean(beta_seg_i < 0)
  or_coronal_vs_dorsal_i <- exp(beta_seg_i)
  or_ci_i <- quantile(or_coronal_vs_dorsal_i, probs = c(0.025, 0.975))

  cat("\nSEGMENT CLASS COMPARISON (DORSAL VS CORONAL; I-DOMAIN)\n")
  cat(sprintf("P(dorsals > coronals | i-domain) = %.3f\n", prob_dorsals_gt_coronals_i))
  cat(sprintf(
    "OR_coronal_vs_dorsal (i-domain) = %.3f [95%% CI: %.3f, %.3f]\n",
    median(or_coronal_vs_dorsal_i),
    or_ci_i[1],
    or_ci_i[2]
  ))
  cat("  (OR < 1 ⇒ dorsals more likely to palatalize among /i/ plurals)\n")

  cat("\nSEGMENT CLASS COMPARISON: DORSAL VS CORONAL (I+E DOMAIN, NO NDEB)\n")

  nouns_ie_classified <- nouns_opp |>
    filter(
      opportunity %in% plural_opportunities,
      !is.na(mutation),
      !(nde_class %in% ndeb_classes)
    ) |>
    mutate(
      segment_class = segment_class_factor(stem_final),
      opportunity = factor(
        opportunity,
        levels = plural_opportunities # c("i", "e")
      ),
      suffix_group = suffix_group_factor(lemma_suffix)
    )

  invisible(capture.output(
    {
      model_class_ie <- fit_brms_bernoulli(
        mutation ~ segment_class + opportunity + suffix_group,
        nouns_ie_classified,
        seed = 124
      )
    },
    type = "output"
  ))

  cat("\nSUMMARY: SEGMENT CLASS MODEL (I+E DOMAIN)\n")
  print(summary(model_class_ie))

  draws_ie <- as_draws_df(model_class_ie)
  beta_seg_ie <- draws_ie[["b_segment_classcoronal"]] # still coronal vs dorsal
  prob_dorsals_gt_coronals_ie <- mean(beta_seg_ie < 0)
  or_coronal_vs_dorsal_ie <- exp(beta_seg_ie)
  or_ci_ie <- quantile(or_coronal_vs_dorsal_ie, probs = c(0.025, 0.975))

  cat("\nSEGMENT CLASS COMPARISON (DORSAL VS CORONAL; I+E DOMAIN)\n")
  cat(sprintf("P(dorsals > coronals | i+e) = %.3f\n", prob_dorsals_gt_coronals_ie))
  cat(sprintf(
    "OR_coronal_vs_dorsal (i+e) = %.3f [95%% CI: %.3f, %.3f]\n",
    median(or_coronal_vs_dorsal_ie),
    or_ci_ie[1],
    or_ci_ie[2]
  ))
  cat("  (OR < 1 ⇒ dorsals more likely to palatalize across /i/ and /e/ plurals)\n")
}

cat_section("ANALYSIS FINISHED")
