#!/usr/bin/Rscript

# Romanian Palatalization Analysis

# Run styler::style_file("analysis/analyze_romanian_palatalization.R") after edits

# =========================================================================
# Setup
# =========================================================================

required_pkgs <- c(
  "dplyr", "readr", "stringr", "tidyr", "broom",
  "brms", "posterior", "loo"
)

# Check and install missing packages
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

# Check and install cmdstanr separately (requires r-universe repo)
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  cat("Installing cmdstanr from r-universe...\n")
  install.packages("cmdstanr", repos = c(cmdstan_repo, getOption("repos")))
}

# Check if CmdStan backend is installed
# If installation fails, see: https://mc-stan.org/cmdstanr/articles/cmdstanr.html
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

segments_of_interest <- c("c", "g", "t", "d", "s", "z")
front_verb_suffixes <- c("-i", "-ui")
suffix_interest <- c("-ic", "-ist", "-esc", "-ică", "-ice")

# NDEB = non-derived exception base lemmas (gimpe / ochi / păduche patterns)
ndeb_classes <- c("gimpe", "ochi", "paduchi")
ndeb_observable <- c("ochi", "paduchi")

plural_opportunities <- c("i", "e")
plural_opportunities_all <- c("i", "e", "uri", "none")

# =========================================================================
# Analysis Constants
# =========================================================================

MIN_SAMPLE_SIZE_BAYESIAN <- 5L # minimum N for Bayesian TP
SMALL_CELL_THRESHOLD <- 20L # flag small downsampled cells

PRECISION_THETA <- 2L # e.g. theta_N
PRECISION_PROB <- 4L # rates / probabilities

CLUSTER_TYPES <- c("st", "sc", "ct") # clusters included in cluster TP

# =========================================================================
# Helper Functions
# =========================================================================

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
      majority_mutates = case_when(
        .data$N == 0L ~ NA,
        .data$mutated == .data$non_mutated ~ NA,
        TRUE ~ .data$mutated > .data$non_mutated
      ),
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

      cat(sprintf("  Fitting <%s> + <-%s> (%s)...\n", seg, opp, subset_label))

      invisible(capture.output(
        {
          model_tp <- brm(
            mutation ~ 1,
            data = seg_data,
            family = bernoulli(link = "logit"),
            prior = prior(normal(0, 1.5), class = "Intercept"),
            chains = 4, iter = 2000, warmup = 1000, cores = 4,
            backend = "cmdstanr",
            seed = seed_value,
            refresh = 0,
            silent = 2
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

# =========================================================================
# Data Input
# =========================================================================

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "data/romanian_lexicon_with_freq.csv"
output_log <- file.path("analysis", "romanian_palatalization_analysis.log")

cat("Working directory:", getwd(), "\n")
cat("Log file:", normalizePath(output_log, mustWork = FALSE), "\n")

sink(output_log, split = TRUE)

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
    col_types = cols(
      derived_adj = col_character(),
      ipa_derived_adj = col_character(),
      ipa_raw_lemma = col_character(),
      ipa_raw_pl = col_character()
    )
  )
)

# Check for parsing problems
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
    freq_ron_wikipedia_2021_1M = as.numeric(freq_ron_wikipedia_2021_1M)
  ) |>
  mutate(
    across(c(pos, gender, stem_final, cluster, plural), ~ trimws(as.character(.))),
    pos = toupper(pos),
    opportunity = if_else(opportunity %in% plural_opportunities_all, opportunity, NA_character_),
    nde_class = if_else(is.na(nde_class) | nde_class == "", "none", nde_class),
    lemma_suffix = if_else(is.na(lemma_suffix) | lemma_suffix == "", "none", lemma_suffix),
    plural_ending = case_when(
      is.na(plural) | plural == "" ~ "none",
      str_ends(plural, "uri") ~ "uri",
      str_ends(plural, "i") ~ "i",
      str_ends(plural, "e") ~ "e",
      TRUE ~ "other"
    ),
    # TP-specific opportunity that includes NDEB lemmas in the i/e domain
    opportunity_tp = case_when(
      nde_class %in% ndeb_classes &
        opportunity == "none" &
        plural_ending %in% c("i", "e") ~ plural_ending,

      # NDEB with plural in -uri: treat as i-type for TP purposes
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
    # Use TP-opportunity for all downstream grouping (includes NDEB)
    opportunity = opportunity_tp,
    cluster_simple = if_else(
      cluster %in% CLUSTER_TYPES,
      cluster,
      NA_character_
    ),
    exception_category = case_when(
      mutation ~ "undergoes",
      nde_class %in% ndeb_classes ~ paste0("NDEB_", nde_class),
      lemma_suffix %in% suffix_interest ~ paste0("suffix_", lemma_suffix),
      is_true_exception ~ "true_exception",
      TRUE ~ "other_non_undergoer"
    )
  )

# Core i/e-domain, excluding NDEB (for clean grammar analysis)
nouns_opp_no_ndeb <- nouns_opp |>
  filter(!(nde_class %in% ndeb_classes))

ndeb_rows <- filter(nouns, nde_class %in% ndeb_classes)

cat("BASIC COUNTS\n")
cat("Total rows:", nrow(lex), "\n")
cat("Nouns:", nrow(nouns), "\n")
cat("Nouns in i/e domain with target segments (incl. NDEB):", nrow(nouns_opp), "\n")
cat("Nouns in i/e domain with target segments (NO NDEB):", nrow(nouns_opp_no_ndeb), "\n\n")

# =========================================================================
# Quality Control
# =========================================================================

cat("QC: GENDER ON NOUN ROWS\n")
nouns_missing_gender <- filter(nouns, is.na(gender) | gender == "")
cat("Missing gender:", nrow(nouns_missing_gender), "\n")
if (nrow(nouns_missing_gender) > 0) {
  nouns_missing_gender |>
    select(lemma, gloss, gender, source, notes) |>
    head(10) |>
    print()
}

cat("\nQC: MUTATION VS OPPORTUNITY\n")
bad_mut_opp <- filter(nouns, mutation, !(opportunity %in% plural_opportunities))
cat("Mutation outside i/e opportunity:", nrow(bad_mut_opp), "\n")
if (nrow(bad_mut_opp) > 0) {
  bad_mut_opp |>
    select(lemma, plural, stem_final, opportunity, mutation) |>
    head(10) |>
    print()
}

cat("\nQC: NDEB ITEMS AND OPPORTUNITY\n")
cat("Total NDEB nouns:", nrow(ndeb_rows), "\n")
count(ndeb_rows, nde_class) |> print()

ndeb_ochi_pad <- filter(ndeb_rows, nde_class %in% ndeb_observable)
ndeb_gimpe <- filter(ndeb_rows, nde_class == "gimpe")

cat("\nNDEB nouns of ochi/păduche type (observable DE exceptions):", nrow(ndeb_ochi_pad), "\n")
if (nrow(ndeb_ochi_pad) > 0) {
  ndeb_ochi_pad |>
    select(lemma, plural, stem_final, nde_class, opportunity) |>
    arrange(nde_class, stem_final, lemma) |>
    print(n = Inf, width = Inf)
}

cat("\nNDEB nouns of gimpe type (unobservable as DE; excluded from exception counts):", nrow(ndeb_gimpe), "\n")

cat("\nQC: OPPORTUNITY VS PLURAL ENDING\n")
nouns |>
  count(opportunity, plural_ending) |>
  arrange(opportunity, plural_ending) |>
  print()

cat("\nQC: MISMATCHES BETWEEN OPPORTUNITY AND PLURAL\n")
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

cat("\nQC: SUFFIX ANNOTATIONS\n")
suffix_rows <- filter(nouns, lemma_suffix %in% suffix_interest)
cat("Nouns with tracked suffixes:", nrow(suffix_rows), "\n")
count(suffix_rows, lemma_suffix) |> print()

suffix_flag_true_nomut <- filter(suffix_rows, suffix_triggers_plural_mutation, !mutation)
cat("\nSuffix marked as trigger but lemma does not mutate:", nrow(suffix_flag_true_nomut), "\n")

suffix_flag_true_notarget <- filter(suffix_rows, suffix_triggers_plural_mutation, !target_is_suffix)
cat("\nSuffix marked as trigger but not marked as target site:", nrow(suffix_flag_true_notarget), "\n")

cat("\nQC: PALATAL CONSONANT IN PLURAL\n")
palatal_nouns <- filter(nouns, !is.na(palatal_consonant_pl), palatal_consonant_pl != "")
cat("Nouns with palatal_consonant_pl populated:", nrow(palatal_nouns), "\n")

palatal_summary <- count(palatal_nouns, stem_final, palatal_consonant_pl, sort = TRUE)
cat("\nPalatal consonant distribution by stem_final:\n")
print(palatal_summary, n = Inf, width = Inf)

palatal_no_mutation <- filter(nouns, !is.na(palatal_consonant_pl), palatal_consonant_pl != "", !mutation)
cat("\nNouns with palatal_consonant_pl but mutation=FALSE:", nrow(palatal_no_mutation), "\n")
if (nrow(palatal_no_mutation) > 0) {
  palatal_no_mutation |>
    select(lemma, plural, stem_final, palatal_consonant_pl, mutation, opportunity) |>
    head(10) |>
    print(n = Inf, width = Inf)
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
    print(n = Inf, width = Inf)
}

cat("\nQC: DUPLICATE LEMMAS\n")
lemma_dups <- count(lex, lemma, pos, sort = TRUE) |> filter(n > 1)
cat("Duplicate (lemma, pos) pairs:", nrow(lemma_dups), "\n")
if (nrow(lemma_dups) > 0) head(lemma_dups, 25) |> print()

# =========================================================================
# Inflection vs Derivation
# =========================================================================

cat("\nINFLECTION VS DERIVATION: LEMMA-BASED PATTERNS (NOUNS & ADJECTIVES)\n")

has_verb_deriv_cols <- all(c("derived_verbs", "ipa_derived_verbs", "deriv_suffixes") %in% names(lex))
has_adj_deriv_cols <- all(c("derived_adj", "ipa_derived_adj") %in% names(lex))

if (!has_verb_deriv_cols && !has_adj_deriv_cols) {
  cat("No derivational columns present; skipping all inflection/derivation checks.\n")
} else {
  noun_base_inflect <- nouns |>
    filter(
      opportunity %in% plural_opportunities,
      !is.na(stem_final),
      stem_final %in% segments_of_interest,
      !is.na(mutation)
    ) |>
    mutate(mutation_inflect = as.logical(mutation))

  # Nouns → Verbs
  cat("\n(1) NOUN LEMMAS: INFLECTIONAL PLURALS VS DENOMINAL VERBS\n")

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
      filter(!is.na(mutation_deriv_verb)) |>
      arrange(lemma) |>
      distinct(lemma, .keep_all = TRUE)

    cat("Denominal N–V lemmas in i/e domain with usable IPA:", nrow(denom_pairs), "\n")
    cat("  (of which with front-vowel verbal suffix -i/-ui:", sum(denom_pairs$verb_suffix_front, na.rm = TRUE), ")\n")

    if (nrow(denom_pairs) > 0) {
      nv_patterns <- denom_pairs |>
        mutate(
          pattern = case_when(
            mutation_inflect & mutation_deriv_verb ~ "both_mutate",
            !mutation_inflect & !mutation_deriv_verb ~ "both_nonmutate",
            mutation_inflect & !mutation_deriv_verb ~ "inflection_only",
            !mutation_inflect & mutation_deriv_verb ~ "derivation_only"
          )
        ) |>
        count(pattern, sort = TRUE) |>
        mutate(prop = n / sum(n))

      cat("\nInflection vs derivation agreement patterns (N → V):\n")
      print(nv_patterns, n = Inf, width = Inf)

      tab_nv <- with(denom_pairs, table(
        inflection_mutates = mutation_inflect,
        verb_derivation_mutates = mutation_deriv_verb
      ))

      if (all(dim(tab_nv) == c(2L, 2L))) {
        cat("\nMcNemar test (inflection_mutates vs verb_derivation_mutates):\n")
        print(broom::tidy(stats::mcnemar.test(tab_nv)), n = Inf, width = Inf)
      } else {
        cat("\nContingency table (N → V) not 2×2; skipping McNemar test.\n")
        print(tab_nv)
      }

      cat("\nLogistic regression: does denominal verb palatalization track inflection?\n")
      model_inf_vs_verb <- glm(mutation_deriv_verb ~ mutation_inflect, data = denom_pairs, family = binomial())
      print(broom::tidy(model_inf_vs_verb), n = Inf, width = Inf)

      # >>> ADDED: TP-style summary for N→V (full lexicon) <<<
      cat("\nDERIVATIONAL SUMMARY TABLE (N → V, FULL LEXICON; GOOGLE SHEET FORMAT)\n")
      nv_tp_full <- denom_pairs |>
        mutate(base_plural_mutates = mutation_inflect) |>
        group_by(base_plural_mutates) |>
        summarise(
          mutated = sum(mutation_deriv_verb, na.rm = TRUE),
          non_mutated = sum(!mutation_deriv_verb, na.rm = TRUE),
          .groups = "drop"
        ) |>
        mutate(
          type = if_else(
            base_plural_mutates,
            "N→V derivation, base plural mutated",
            "N→V derivation, base plural non-mut."
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

      print(nv_tp_full, n = Inf, width = Inf)
      # <<< END ADDED >>>
    }
  } else {
    cat("No denominal verb columns found; skipping N → V comparison.\n")
  }

  # Nouns → Adjectives
  cat("\n(2) NOUN LEMMAS: INFLECTIONAL PLURALS VS DENOMINAL ADJECTIVES\n")

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
      na_patterns <- noun_adj_pairs |>
        mutate(
          pattern = case_when(
            mutation_inflect & mutation_deriv_adj ~ "both_mutate",
            !mutation_inflect & !mutation_deriv_adj ~ "both_nonmutate",
            mutation_inflect & !mutation_deriv_adj ~ "inflection_only",
            !mutation_inflect & mutation_deriv_adj ~ "derivation_only"
          )
        ) |>
        count(pattern, sort = TRUE) |>
        mutate(prop = n / sum(n))

      cat("\nInflection vs derivation agreement patterns (N → Adj):\n")
      print(na_patterns, n = Inf, width = Inf)

      tab_na <- with(noun_adj_pairs, table(
        inflection_mutates = mutation_inflect,
        adj_derivation_mutates = mutation_deriv_adj
      ))

      if (all(dim(tab_na) == c(2L, 2L))) {
        cat("\nMcNemar test (inflection_mutates vs adj_derivation_mutates):\n")
        print(broom::tidy(stats::mcnemar.test(tab_na)), n = Inf, width = Inf)
      } else {
        cat("\nContingency table (N → Adj) not 2×2; skipping McNemar test.\n")
        print(tab_na)
      }

      cat("\nLogistic regression: does denominal adjective palatalization track inflection?\n")
      model_inf_vs_adj <- glm(mutation_deriv_adj ~ mutation_inflect, data = noun_adj_pairs, family = binomial())
      print(broom::tidy(model_inf_vs_adj), n = Inf, width = Inf)

      # >>> ADDED: TP-style summary for N→Adj (full lexicon) <<<
      cat("\nDERIVATIONAL SUMMARY TABLE (N → Adj, FULL LEXICON; GOOGLE SHEET FORMAT)\n")
      na_tp_full <- noun_adj_pairs |>
        mutate(base_plural_mutates = mutation_inflect) |>
        group_by(base_plural_mutates) |>
        summarise(
          mutated = sum(mutation_deriv_adj, na.rm = TRUE),
          non_mutated = sum(!mutation_deriv_adj, na.rm = TRUE),
          .groups = "drop"
        ) |>
        mutate(
          type = if_else(
            base_plural_mutates,
            "N→Adj derivation, base plural mutated",
            "N→Adj derivation, base plural non-mut."
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

      print(na_tp_full, n = Inf, width = Inf)
      # <<< END ADDED >>>
    }
  } else {
    cat("No denominal adjective columns found; skipping N → Adj comparison.\n")
  }

  # Adjectives
  cat("\n(3) ADJECTIVE LEMMAS: INFLECTIONAL PLURALS VS DERIVATIONS\n")

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
      av_patterns <- adj_verb_pairs |>
        mutate(
          pattern = case_when(
            mutation_inflect & mutation_deriv_verb ~ "both_mutate",
            !mutation_inflect & !mutation_deriv_verb ~ "both_nonmutate",
            mutation_inflect & !mutation_deriv_verb ~ "inflection_only",
            !mutation_inflect & mutation_deriv_verb ~ "derivation_only"
          )
        ) |>
        count(pattern, sort = TRUE) |>
        mutate(prop = n / sum(n))

      cat("\nInflection vs derivation agreement patterns (Adj → V):\n")
      print(av_patterns, n = Inf, width = Inf)

      tab_av <- with(adj_verb_pairs, table(
        inflection_mutates = mutation_inflect,
        verb_derivation_mutates = mutation_deriv_verb
      ))

      if (all(dim(tab_av) == c(2L, 2L))) {
        cat("\nMcNemar test (Adj plural vs Adj→V derivation):\n")
        print(broom::tidy(stats::mcnemar.test(tab_av)), n = Inf, width = Inf)
      } else {
        cat("\nContingency table (Adj → V) not 2×2; skipping McNemar test.\n")
        print(tab_av)
      }

      cat("\nLogistic regression: does verb derivation from adjectives track adjectival inflection?\n")
      model_adj_inf_vs_verb <- glm(mutation_deriv_verb ~ mutation_inflect, data = adj_verb_pairs, family = binomial())
      print(broom::tidy(model_adj_inf_vs_verb), n = Inf, width = Inf)
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
      aa_patterns <- adj_adj_pairs |>
        mutate(
          pattern = case_when(
            mutation_inflect & mutation_deriv_adj ~ "both_mutate",
            !mutation_inflect & !mutation_deriv_adj ~ "both_nonmutate",
            mutation_inflect & !mutation_deriv_adj ~ "inflection_only",
            !mutation_inflect & mutation_deriv_adj ~ "derivation_only"
          )
        ) |>
        count(pattern, sort = TRUE) |>
        mutate(prop = n / sum(n))

      cat("\nInflection vs derivation agreement patterns (Adj → Adj):\n")
      print(aa_patterns, n = Inf, width = Inf)

      tab_aa <- with(adj_adj_pairs, table(
        inflection_mutates = mutation_inflect,
        adj_derivation_mutates = mutation_deriv_adj
      ))

      if (all(dim(tab_aa) == c(2L, 2L))) {
        cat("\nMcNemar test (Adj plural vs Adj→Adj derivation):\n")
        print(broom::tidy(stats::mcnemar.test(tab_aa)), n = Inf, width = Inf)
      } else {
        cat("\nContingency table (Adj → Adj) not 2×2; skipping McNemar test.\n")
        print(tab_aa)
      }

      cat("\nLogistic regression: does adjectival derivation track adjectival inflection?\n")
      model_adj_inf_vs_adj <- glm(mutation_deriv_adj ~ mutation_inflect, data = adj_adj_pairs, family = binomial())
      print(broom::tidy(model_adj_inf_vs_adj), n = Inf, width = Inf)
    }
  } else {
    cat("No adjective lemmas with both plural and non-trivial derived adjectives; skipping Adj → Adj comparison.\n")
  }
}

# =========================================================================
# Descriptive Summaries
# =========================================================================

cat("\nSEGMENT-WISE MUTATION RATES (I+E COMBINED, NO NDEB)\n")
seg_summary <- nouns_opp_no_ndeb |>
  group_by(stem_final) |>
  summarise(N_opp = n(), N_mut = sum(mutation, na.rm = TRUE), .groups = "drop") |>
  calc_rate() |>
  arrange(stem_final)
print(seg_summary, n = Inf, width = Inf)

cat("\nMUTATION RATES BY SEGMENT AND PLURAL TYPE (I VS E, NO NDEB)\n")
seg_by_opp <- nouns_opp_no_ndeb |>
  group_by(stem_final, opportunity) |>
  summarise(N_opp = n(), N_mut = sum(mutation, na.rm = TRUE), .groups = "drop") |>
  calc_rate() |>
  arrange(stem_final, opportunity)
print(seg_by_opp, n = Inf, width = Inf)

cat("\nCLUSTER INVENTORY IN I/E DOMAIN (NO NDEB)\n")
cluster_inventory <- nouns_opp_no_ndeb |>
  filter(!is.na(cluster), cluster != "") |>
  count(stem_final, cluster, sort = TRUE)
print(cluster_inventory, n = Inf, width = Inf)

cat("\nCLUSTER EFFECTS ON MUTATION (NO NDEB)\n")
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
print(cluster_summary, n = Inf, width = Inf)

cat("\nSTRUCTURE OF NON-UNDERGOERS (I/E DOMAIN)\n")
non_under_summary <- nouns_opp |>
  filter(!mutation) |>
  group_by(stem_final, exception_category) |>
  summarise(N = n(), .groups = "drop") |>
  group_by(stem_final) |>
  mutate(total_non_under = sum(N), prop = N / total_non_under) |>
  arrange(exception_category, stem_final, desc(prop))
print(non_under_summary, n = Inf, width = Inf)

cat("\nNDE DISTRIBUTION BY SEGMENT\n")
nde_by_seg <- ndeb_rows |>
  group_by(nde_class, stem_final) |>
  summarise(N = n(), .groups = "drop") |>
  arrange(nde_class, stem_final)
print(nde_by_seg, n = Inf, width = Inf)

cat("\nPALATALIZATION RATES FOR <-că, -gă> NOUNS\n")
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
print(cg_cag_summary, n = Inf, width = Inf)

cat("\nSUFFIX PATTERNS (ALL TRACKED SUFFIXES)\n")
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
print(suffix_summary_all, n = Inf, width = Inf)

cat("\nSUFFIX PATTERNS WHERE SUFFIX IS TARGET\n")
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
print(suffix_target_summary, n = Inf, width = Inf)

cat("\nCOMPARISON: 'HAS SUFFIX' VS 'SUFFIX IS TARGET'\n")
suffix_diff <- anti_join(
  suffix_summary_all,
  suffix_target_summary,
  by = c("lemma_suffix", "stem_final", "opportunity")
)
if (nrow(suffix_diff) > 0) {
  cat("Suffix present but not annotated as target in these cells:\n")
  print(suffix_diff, n = Inf, width = Inf)
} else {
  cat("All tracked suffix rows are also marked as targets.\n")
}

cat("\nTRUE EXCEPTIONS IN I/E DOMAIN (NON-NDEB)\n")
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
print(true_exc_by_seg, n = Inf, width = Inf)

cat("\nNDEB EXCEPTIONS OF OCHI/PĂDUCHE TYPE (OUTSIDE ALIGNMENT-BASED I/E OPPORTUNITY)\n")
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
    print(n = Inf, width = Inf)
}

if (nrow(true_exc) > 0) {
  cat("\nSample of true exceptions:\n")
  true_exc |>
    select(lemma, plural, stem_final, opportunity, nde_class, lemma_suffix, notes) |>
    head(30) |>
    print(n = Inf, width = Inf)

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

# Use the TOP N most frequent lemmas
sample_lexeme_sizes <- c(1000L, 2500L, 5000L, 10000L)

nouns_opp_freq <- nouns_opp |>
  group_by(lemma) |>
  summarise(lemma_freq = max(freq_ron_wikipedia_2021_1M, na.rm = TRUE), .groups = "drop") |>
  mutate(lemma_freq = if_else(is.na(lemma_freq) | !is.finite(lemma_freq), 0, lemma_freq))

nouns_opp_freq_pos <- filter(nouns_opp_freq, lemma_freq > 0)

nouns_opp_with_freq <- nouns_opp |>
  left_join(nouns_opp_freq_pos, by = "lemma")

cat("Unique lemmas in i/e domain (all):", n_distinct(nouns_opp$lemma), "\n")
cat("Unique lemmas with freq > 0:", nrow(nouns_opp_freq_pos), "\n\n")

seg_tp_ie_ds_all <- NULL
nouns_opp_down_single <- NULL

if (nrow(nouns_opp_freq_pos) == 0L) {
  cat("No positive-frequency lemmas in i/e domain; skipping frequency downsampling.\n\n")
} else {
  seg_tp_ie_ds_list <- vector("list", length(sample_lexeme_sizes))
  idx <- 1L

  for (n_lex in sample_lexeme_sizes) {
    # Take the TOP N most frequent lemmas
    target_n_lex <- min(n_lex, nrow(nouns_opp_freq_pos))

    top_lemmas <- nouns_opp_freq_pos |>
      arrange(desc(lemma_freq)) |>
      slice_head(n = target_n_lex) |>
      pull(lemma)

    nouns_opp_down <- nouns_opp_with_freq |>
      filter(lemma %in% top_lemmas) |>
      arrange(desc(lemma_freq), lemma, plural) |>
      distinct(lemma, .keep_all = TRUE)

    # Use the 1000 most frequent lemmas as the reference downsampled lexicon
    if (is.null(nouns_opp_down_single) && n_lex == 1000L) {
      nouns_opp_down_single <- nouns_opp_down
      cat("Reference downsampled lexicon (top 1000 most frequent lemmas):\n")
      cat("  target_lexemes:", target_n_lex, "\n")
      cat("  unique lemmas:", n_distinct(nouns_opp_down_single$lemma), "\n")
      cat("  rows:", nrow(nouns_opp_down_single), "\n\n")
    }

    seg_tp_ie_ds_list[[idx]] <- nouns_opp_down |>
      group_by(stem_final, opportunity) |>
      summarise(
        mutated = sum(mutation, na.rm = TRUE),
        non_mutated = sum(!mutation, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(sample_lexemes = n_lex)

    idx <- 1L + idx
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

  cat("Segment × opportunity mutation / TP summary across frequency-filtered lexicons (top N most frequent):\n")
  print(seg_tp_ie_ds_all, n = Inf, width = Inf)
  cat("\n")
}

# =========================================================================
# Tolerance Principle: Full Data
# =========================================================================

cat("\nTOLERANCE PRINCIPLE: SEGMENT-LEVEL PATTERNS (FULL DATA)\n")

seg_tp_all <- compute_segment_tp_tables(nouns_opp)
print(seg_tp_all, n = Inf, width = Inf)

cat("\nTOLERANCE PRINCIPLE: SEGMENT-LEVEL PATTERNS (FULL DATA, NO NDEB)\n")

nouns_opp_no_ndeb <- nouns_opp |>
  filter(!(nde_class %in% ndeb_classes))

seg_tp_all_no_ndeb <- compute_segment_tp_tables(nouns_opp_no_ndeb, label_suffix = " (no NDEB)")
print(seg_tp_all_no_ndeb, n = Inf, width = Inf)

# =========================================================================
# NDEB CONTRIBUTION PER TYPE (FULL LEXICON & DOWNSAMPLED)
# =========================================================================

cat("\nTOLERANCE PRINCIPLE: NDEB CONTRIBUTION PER TYPE (FULL LEXICON)\n")

nouns_opp_ndeb <- nouns_opp |>
  filter(nde_class %in% ndeb_classes)

seg_tp_all_ndeb <- compute_segment_tp_tables(nouns_opp_ndeb, label_suffix = " (NDEB)")
print(seg_tp_all_ndeb, n = Inf, width = Inf)

cat("\nTOLERANCE PRINCIPLE: NDEB CONTRIBUTION PER TYPE (REFERENCE DOWNSAMPLED LEXICON)\n")

if (!is.null(nouns_opp_down_single) && nrow(nouns_opp_down_single) > 0L) {
  nouns_opp_down_ndeb <- nouns_opp_down_single |>
    filter(nde_class %in% ndeb_classes)

  seg_tp_all_ds_ndeb <- compute_segment_tp_tables(nouns_opp_down_ndeb, label_suffix = " (downsampled, NDEB)")
  print(seg_tp_all_ds_ndeb, n = Inf, width = Inf)
} else {
  cat("No reference downsampled lexicon available; skipping NDEB-by-type counts (downsampled).\n")
}

# =========================================================================
# Tolerance Principle: Downsampled Data
# =========================================================================

if (!is.null(nouns_opp_down_single) && nrow(nouns_opp_down_single) > 0L) {
  # >>> ADDED: Derivational summaries for reference downsampled lexicon <<<
  cat("\nDERIVATIONAL SUMMARY TABLES (REFERENCE DOWNSAMPLED LEXICON: TOP 1000 MOST FREQUENT LEMMAS)\n")

  lemmas_ds <- unique(nouns_opp_down_single$lemma)

  if (has_verb_deriv_cols) {
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
      nv_tp_ds <- denom_pairs_ds |>
        mutate(base_plural_mutates = mutation_inflect) |>
        group_by(base_plural_mutates) |>
        summarise(
          mutated = sum(mutation_deriv_verb, na.rm = TRUE),
          non_mutated = sum(!mutation_deriv_verb, na.rm = TRUE),
          .groups = "drop"
        ) |>
        mutate(
          type = if_else(
            base_plural_mutates,
            "N→V derivation, base plural mutated (downsampled)",
            "N→V derivation, base plural non-mut. (downsampled)"
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

      print(nv_tp_ds, n = Inf, width = Inf)
    }
  }

  if (has_adj_deriv_cols) {
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
      na_tp_ds <- noun_adj_pairs_ds |>
        mutate(base_plural_mutates = mutation_inflect) |>
        group_by(base_plural_mutates) |>
        summarise(
          mutated = sum(mutation_deriv_adj, na.rm = TRUE),
          non_mutated = sum(!mutation_deriv_adj, na.rm = TRUE),
          .groups = "drop"
        ) |>
        mutate(
          type = if_else(
            base_plural_mutates,
            "N→Adj derivation, base plural mutated (downsampled)",
            "N→Adj derivation, base plural non-mut. (downsampled)"
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

      print(na_tp_ds, n = Inf, width = Inf)
    }
  }
  # <<< END ADDED >>>

  cat("\nTOLERANCE PRINCIPLE: SEGMENT-LEVEL PATTERNS (REFERENCE DOWNSAMPLED LEXICON)\n")

  seg_tp_all_ds <- compute_segment_tp_tables(nouns_opp_down_single, label_suffix = " (downsampled)")
  print(seg_tp_all_ds, n = Inf, width = Inf)

  cat("\nTOLERANCE PRINCIPLE: SEGMENT-LEVEL PATTERNS (REFERENCE DOWNSAMPLED, NO NDEB)\n")

  nouns_opp_down_single_no_ndeb <- nouns_opp_down_single |>
    filter(!(nde_class %in% ndeb_classes))

  seg_tp_all_ds_no_ndeb <- compute_segment_tp_tables(nouns_opp_down_single_no_ndeb, label_suffix = " (downsampled, no NDEB)")
  print(seg_tp_all_ds_no_ndeb, n = Inf, width = Inf)

  cat("\nTOLERANCE PRINCIPLE: SEGMENT-LEVEL COMPARISON (I VS E; FULL VS DOWNSAMPLED)\n")

  # Generate raw counts for comparison
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
  print(seg_tp_ie_compare, n = Inf, width = Inf)

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

cat("\nTOLERANCE PRINCIPLE: CLUSTER PATTERNS (FULL DATA)\n")

nouns_opp_clusters <- nouns_opp |>
  filter(!is.na(cluster_simple))

cluster_tp_all <- compute_segment_tp_tables(nouns_opp_clusters, group_var = cluster_simple)
print(cluster_tp_all, n = Inf, width = Inf)

cat("\nTOLERANCE PRINCIPLE: CLUSTER PATTERNS (REFERENCE DOWNSAMPLED LEXICON)\n")

if (!is.null(nouns_opp_down_single) && nrow(nouns_opp_down_single) > 0L) {
  nouns_opp_down_clusters <- nouns_opp_down_single |>
    filter(!is.na(cluster_simple))

  cluster_tp_all_ds <- compute_segment_tp_tables(nouns_opp_down_clusters, label_suffix = " (downsampled)", group_var = cluster_simple)
  print(cluster_tp_all_ds, n = Inf, width = Inf)
} else {
  cat("No reference downsampled lexicon available; skipping cluster TP (downsampled).\n")
}

# =========================================================================
# Tolerance Principle: NDEB Counts
# =========================================================================

cat("\nNDEB BY CLASS (FULL LEXICON AND REFERENCE DOWNSAMPLED)\n")

# Helper to label rows the way you want in the spreadsheet
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

print(ndeb_tp_all, n = Inf, width = Inf)

# =========================================================================
# Bayesian Tolerance Principle
# =========================================================================

cat("\nBAYESIAN TOLERANCE PRINCIPLE: SEGMENT × OPPORTUNITY (FULL DATA, NO NDEB)\n")
tolerance_bayesian_results <- run_bayesian_tp(nouns_opp_no_ndeb, subset_label = "full_noNDEB", seed_value = 123L)

cat("\nBAYESIAN TOLERANCE PRINCIPLE: SEGMENT × OPPORTUNITY (REFERENCE DOWNSAMPLED, NO NDEB)\n")
if (!is.null(nouns_opp_down_single) && nrow(nouns_opp_down_single) > 0L) {
  tolerance_bayesian_results_ds <- run_bayesian_tp(nouns_opp_down_single_no_ndeb, subset_label = "downsampled_noNDEB", seed_value = 456L)
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
  print(as_tibble(tolerance_bayesian_df), n = Inf, width = Inf)
} else {
  cat("\nInsufficient data for Bayesian TP analysis (both subsets)\n\n")
}

# =========================================================================
# Segment Class Comparison
# =========================================================================

cat("\nSEGMENT CLASS COMPARISON: DORSAL VS CORONAL (I-DOMAIN, NO NDEB)\n")

nouns_i_classified <- nouns_opp |>
  filter(
    opportunity == "i",
    !is.na(mutation),
    !(nde_class %in% ndeb_classes) # Exclude NDEB from core grammar analysis
  ) |>
  mutate(
    segment_class = factor(
      if_else(stem_final %in% c("c", "g"), "dorsal", "coronal"),
      levels = c("dorsal", "coronal") # dorsal = baseline
    ),
    suffix_group = case_when(
      lemma_suffix == "-ic" ~ "ic",
      lemma_suffix == "-ist" ~ "ist",
      lemma_suffix %in% c("-ică", "-ice") ~ "ica_ice",
      lemma_suffix == "none" ~ "none",
      TRUE ~ "other"
    ),
    suffix_group = factor(
      suffix_group,
      levels = c("none", "ic", "ist", "ica_ice", "other")
    )
  )

invisible(capture.output(
  {
    model_class_i <- brm(
      mutation ~ segment_class + suffix_group,
      data = nouns_i_classified,
      family = bernoulli(link = "logit"),
      prior = c(
        prior(normal(0, 2.5), class = "Intercept"),
        prior(normal(0, 2.5), class = "b")
      ),
      chains = 4, iter = 2000, warmup = 1000, cores = 4,
      backend = "cmdstanr",
      seed = 123,
      refresh = 0,
      silent = 2
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
    !(nde_class %in% ndeb_classes) # Exclude NDEB from core grammar analysis
  ) |>
  mutate(
    segment_class = factor(
      if_else(stem_final %in% c("c", "g"), "dorsal", "coronal"),
      levels = c("dorsal", "coronal") # dorsal = baseline
    ),
    opportunity = factor(
      opportunity,
      levels = plural_opportunities # c("i", "e")
    ),
    suffix_group = case_when(
      lemma_suffix == "-ic" ~ "ic",
      lemma_suffix == "-ist" ~ "ist",
      lemma_suffix %in% c("-ică", "-ice") ~ "ica_ice",
      lemma_suffix == "none" ~ "none",
      TRUE ~ "other"
    ),
    suffix_group = factor(
      suffix_group,
      levels = c("none", "ic", "ist", "ica_ice", "other")
    )
  )

invisible(capture.output(
  {
    model_class_ie <- brm(
      # main effects: segment_class + opportunity + suffix_group
      # (no interaction ⇒ segment effect is "overall" across i+e, controlling for opportunity)
      mutation ~ segment_class + opportunity + suffix_group,
      data = nouns_ie_classified,
      family = bernoulli(link = "logit"),
      prior = c(
        prior(normal(0, 2.5), class = "Intercept"),
        prior(normal(0, 2.5), class = "b")
      ),
      chains = 4, iter = 2000, warmup = 1000, cores = 4,
      backend = "cmdstanr",
      seed = 124,
      refresh = 0,
      silent = 2
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

cat("\nANALYSIS FINISHED\n")

sink()
