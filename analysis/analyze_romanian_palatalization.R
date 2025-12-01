#!/usr/bin/Rscript

# Romanian Palatalization Analysis
# Analyzes nominal plural palatalization patterns in Romanian

# Run styler::style_file("analysis/analyze_romanian_palatalization.R") after edits

# =========================================================================
# Setup
# =========================================================================

required_pkgs <- c(
  "dplyr", "readr", "stringr", "tidyr", "broom",
  "brms", "posterior", "loo", "cmdstanr"
)

# For cmdstanr setup details, see:
# Fruehwald, Josef. "Getting `Brms` and Stan Set Up."
# https://lin611-2024.github.io/notes/side-notes/content/stan.html

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", "),
    ". Please install them before running this script."
  )
}

suppressPackageStartupMessages(lapply(required_pkgs, library, character.only = TRUE))
options(dplyr.summarise.inform = FALSE)

# =========================================================================
# Constants
# =========================================================================

segments_of_interest <- c("c", "g", "t", "d", "s", "z")
front_verb_suffixes <- c("-i", "-ui")
suffix_interest <- c("-ic", "-ist", "-esc", "-ică", "-ice")
nde_classes <- c("gimpe", "ochi", "paduchi")
nde_observable <- c("ochi", "paduchi")
plural_opportunities <- c("i", "e")
plural_opportunities_all <- c("i", "e", "uri", "none")

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
    stem_final == "d" ~ str_detect(ipa_str, "z|dʒ|dʲ"),
    stem_final == "s" ~ str_detect(ipa_str, "ʃ"),
    stem_final == "z" ~ str_detect(ipa_str, "ʒ"),
    TRUE ~ NA
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

lex <- read_csv(
  input_file,
  show_col_types = FALSE,
  col_types = cols(
    derived_adj = col_character(),
    ipa_derived_adj = col_character()
  )
) |>
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
    )
  )

nouns <- filter(lex, pos == "N")

nouns_opp <- nouns |>
  filter(
    opportunity %in% plural_opportunities,
    !is.na(stem_final),
    stem_final %in% segments_of_interest
  ) |>
  mutate(
    exception_category = case_when(
      mutation ~ "undergoes",
      nde_class %in% nde_classes ~ paste0("NDE_", nde_class),
      lemma_suffix %in% suffix_interest ~ paste0("suffix_", lemma_suffix),
      is_true_exception ~ "true_exception",
      TRUE ~ "other_non_undergoer"
    )
  )

nde_rows <- filter(nouns, nde_class %in% nde_classes)

cat("BASIC COUNTS\n")
cat("Total rows:", nrow(lex), "\n")
cat("Nouns:", nrow(nouns), "\n")
cat("Nouns in i/e domain with target segments:", nrow(nouns_opp), "\n\n")

# =========================================================================
# Frequency-Based Downsampling
# =========================================================================

cat("DOWNSAMPLING (FREQUENCY-BASED) FOR TOLERANCE PRINCIPLE ANALYSIS\n")

sample_tokens_vec <- c(5000L, 10000L, 20000L)
n_sims_per_size <- 50L
base_seed_downsample <- 123L

nouns_opp_freq <- nouns_opp |>
  group_by(lemma) |>
  summarise(lemma_freq = max(freq_ron_wikipedia_2021_1M, na.rm = TRUE), .groups = "drop") |>
  mutate(lemma_freq = if_else(is.na(lemma_freq) | !is.finite(lemma_freq), 0, lemma_freq))

nouns_opp_freq_pos <- filter(nouns_opp_freq, lemma_freq > 0)

cat("Unique lemmas in i/e domain (all):", n_distinct(nouns_opp$lemma), "\n")
cat("Unique lemmas with freq > 0:", nrow(nouns_opp_freq_pos), "\n\n")

seg_tp_ie_ds_all <- NULL
nouns_opp_down_single <- NULL

if (nrow(nouns_opp_freq_pos) == 0L) {
  cat("No positive-frequency lemmas in i/e domain; skipping frequency downsampling.\n\n")
} else {
  total_freq <- sum(nouns_opp_freq_pos$lemma_freq)
  lemma_probs <- nouns_opp_freq_pos$lemma_freq / total_freq

  seg_tp_ie_ds_list <- vector("list", length(sample_tokens_vec) * n_sims_per_size)
  idx <- 1L

  for (tokens in sample_tokens_vec) {
    for (sim in seq_len(n_sims_per_size)) {
      set.seed(base_seed_downsample + tokens * 1000L + sim)

      sampled_idx <- sample(
        seq_len(nrow(nouns_opp_freq_pos)),
        size = tokens,
        replace = TRUE,
        prob = lemma_probs
      )

      sampled_lemmas <- unique(nouns_opp_freq_pos$lemma[sampled_idx])
      nouns_opp_down <- filter(nouns_opp, lemma %in% sampled_lemmas)

      if (is.null(nouns_opp_down_single)) {
        nouns_opp_down_single <- nouns_opp_down
        cat("Reference downsampled lexicon (first draw):\n")
        cat("  sample_tokens:", tokens, "\n")
        cat("  unique lemmas:", n_distinct(nouns_opp_down_single$lemma), "\n")
        cat("  rows:", nrow(nouns_opp_down_single), "\n\n")
      }

      seg_tp_ie_ds_list[[idx]] <- nouns_opp_down |>
        group_by(stem_final, opportunity) |>
        summarise(
          mutated = sum(mutation == TRUE, na.rm = TRUE),
          non_mutated = sum(mutation == FALSE, na.rm = TRUE),
          .groups = "drop"
        ) |>
        mutate(sample_tokens = tokens, sim = sim)

      idx <- idx + 1L
    }
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

  seg_tp_ie_ds_summary <- seg_tp_ie_ds_all |>
    group_by(sample_tokens, stem_final, opportunity) |>
    summarise(
      N_mean = mean(N),
      rate_mean = mean(rate, na.rm = TRUE),
      rate_sd = sd(rate, na.rm = TRUE),
      prop_tolerated = mean(tolerated, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(sample_tokens, stem_final, opportunity)

  cat("Segment × opportunity mutation / TP summary across downsampled lexicons:\n")
  print(seg_tp_ie_ds_summary, n = Inf, width = Inf)
  cat("\n")
}

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

cat("\nQC: NDE ITEMS AND OPPORTUNITY\n")
cat("Total NDE nouns:", nrow(nde_rows), "\n")
count(nde_rows, nde_class) |> print()

nde_ochi_pad <- filter(nde_rows, nde_class %in% nde_observable)
nde_gimpe <- filter(nde_rows, nde_class == "gimpe")

cat("\nNDE nouns of ochi/păduche type (observable DE exceptions per email):", nrow(nde_ochi_pad), "\n")
if (nrow(nde_ochi_pad) > 0) {
  nde_ochi_pad |>
    select(lemma, plural, stem_final, nde_class, opportunity) |>
    arrange(nde_class, stem_final, lemma) |>
    print(n = Inf, width = Inf)
}

cat("\nNDE nouns of gimpe type (unobservable as DE; excluded from exception counts):", nrow(nde_gimpe), "\n")

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

cat("\nSEGMENT-WISE MUTATION RATES (I+E COMBINED)\n")
seg_summary <- nouns_opp |>
  group_by(stem_final) |>
  summarise(N_opp = n(), N_mut = sum(mutation, na.rm = TRUE), .groups = "drop") |>
  calc_rate() |>
  arrange(stem_final)
print(seg_summary, n = Inf, width = Inf)

cat("\nMUTATION RATES BY SEGMENT AND PLURAL TYPE (I VS E)\n")
seg_by_opp <- nouns_opp |>
  group_by(stem_final, opportunity) |>
  summarise(N_opp = n(), N_mut = sum(mutation, na.rm = TRUE), .groups = "drop") |>
  calc_rate() |>
  arrange(stem_final, opportunity)
print(seg_by_opp, n = Inf, width = Inf)

cat("\nCLUSTER INVENTORY IN I/E DOMAIN\n")
cluster_inventory <- nouns_opp |>
  filter(!is.na(cluster), cluster != "") |>
  count(stem_final, cluster, sort = TRUE)
print(cluster_inventory, n = Inf, width = Inf)

cat("\nCLUSTER EFFECTS ON MUTATION\n")
cluster_summary <- nouns_opp |>
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
nde_by_seg <- nde_rows |>
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

cat("\nTRUE EXCEPTIONS IN I/E DOMAIN (NON-NDE)\n")
true_exc <- nouns |>
  filter(
    opportunity %in% plural_opportunities,
    !mutation,
    !(nde_class %in% nde_classes)
  )
cat("Non-mutating, non-NDE nouns:", nrow(true_exc), "\n")

true_exc_by_seg <- true_exc |>
  group_by(stem_final) |>
  summarise(N_true_exc = n(), .groups = "drop") |>
  arrange(stem_final)
print(true_exc_by_seg, n = Inf, width = Inf)

cat("\nNDE EXCEPTIONS OF OCHI/PĂDUCHE TYPE (OUTSIDE ALIGNMENT-BASED I/E OPPORTUNITY)\n")
nde_observable <- c("ochi", "paduchi")
nde_exc_ochi_pad <- nouns |>
  filter(
    nde_class %in% nde_observable,
    !mutation
  )

cat("Non-mutating NDE nouns of ochi/păduche type:", nrow(nde_exc_ochi_pad), "\n")

if (nrow(nde_exc_ochi_pad) > 0) {
  nde_exc_ochi_pad |>
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
# Tolerance Principle: Full Data
# =========================================================================

cat("\nTOLERANCE PRINCIPLE: SEGMENT-LEVEL PATTERNS (FULL DATA)\n")

seg_tp_ie_raw <- nouns_opp |>
  group_by(stem_final, opportunity) |>
  summarise(
    mutated = sum(mutation, na.rm = TRUE),
    non_mutated = sum(!mutation, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(type = paste0("<", stem_final, "> + <-", opportunity, "> plural"))

seg_tp_ie <- tp_table(seg_tp_ie_raw, type, mutated, non_mutated)

seg_tp_comb <- nouns_opp |>
  group_by(stem_final) |>
  summarise(
    mutated = sum(mutation, na.rm = TRUE),
    non_mutated = sum(!mutation, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(type = paste0("<", stem_final, "> + <-i, -e> plural")) |>
  tp_table(type, mutated, non_mutated)

seg_tp_all <- bind_rows(seg_tp_ie, seg_tp_comb) |> arrange(type)
print(seg_tp_all, n = Inf, width = Inf)

# =========================================================================
# Tolerance Principle: Downsampled Data
# =========================================================================

if (!is.null(nouns_opp_down_single) && nrow(nouns_opp_down_single) > 0L) {
  cat("\nTOLERANCE PRINCIPLE: SEGMENT-LEVEL PATTERNS (REFERENCE DOWNSAMPLED LEXICON)\n")

  seg_tp_ie_ds_raw <- nouns_opp_down_single |>
    group_by(stem_final, opportunity) |>
    summarise(
      mutated = sum(mutation, na.rm = TRUE),
      non_mutated = sum(!mutation, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(type = paste0("<", stem_final, "> + <-", opportunity, "> plural (downsampled)"))

  seg_tp_ie_ds <- tp_table(seg_tp_ie_ds_raw, type, mutated, non_mutated)

  seg_tp_comb_ds <- nouns_opp_down_single |>
    group_by(stem_final) |>
    summarise(
      mutated = sum(mutation, na.rm = TRUE),
      non_mutated = sum(!mutation, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(type = paste0("<", stem_final, "> + <-i, -e> plural (downsampled)")) |>
    tp_table(type, mutated, non_mutated)

  seg_tp_all_ds <- bind_rows(seg_tp_ie_ds, seg_tp_comb_ds) |> arrange(type)
  print(seg_tp_all_ds, n = Inf, width = Inf)

  cat("\nTOLERANCE PRINCIPLE: SEGMENT-LEVEL COMPARISON (I VS E; FULL VS DOWNSAMPLED)\n")
  seg_tp_ie_compare <- seg_tp_ie_raw |>
    select(stem_final, opportunity, mutated, non_mutated) |>
    rename(mutated_full = mutated, non_mutated_full = non_mutated) |>
    left_join(
      seg_tp_ie_ds_raw |>
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

  small_ds <- filter(seg_tp_ie_compare, !is.na(N_ds), N_ds < 20)
  if (nrow(small_ds) > 0) {
    cat("\nCells with N_ds < 20 in segment × opportunity (downsampled):\n")
    select(small_ds, stem_final, opportunity, N_ds, rate_ds) |> print()
  }
} else {
  cat("\nNo reference downsampled lexicon; skipping single-sample TP comparison.\n")
}

# =========================================================================
# Tolerance Principle: Cluster Patterns
# =========================================================================

cat("\nTOLERANCE PRINCIPLE: CLUSTER PATTERNS (FULL DATA)\n")

cluster_tp_ie <- nouns_opp |>
  mutate(cluster_simple = if_else(cluster %in% c("st", "sc", "ct"), cluster, NA_character_)) |>
  filter(!is.na(cluster_simple)) |>
  group_by(cluster_simple, opportunity) |>
  summarise(
    mutated = sum(mutation, na.rm = TRUE),
    non_mutated = sum(!mutation, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(type = paste0("<", cluster_simple, "> + <-", opportunity, "> plural")) |>
  tp_table(type, mutated, non_mutated)

cluster_tp_comb <- nouns_opp |>
  mutate(cluster_simple = if_else(cluster %in% c("st", "sc", "ct"), cluster, NA_character_)) |>
  filter(!is.na(cluster_simple)) |>
  group_by(cluster_simple) |>
  summarise(
    mutated = sum(mutation, na.rm = TRUE),
    non_mutated = sum(!mutation, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(type = paste0("<", cluster_simple, "> + <-i, -e> plural")) |>
  tp_table(type, mutated, non_mutated)

cluster_tp_all <- bind_rows(cluster_tp_ie, cluster_tp_comb) |> arrange(type)
print(cluster_tp_all, n = Inf, width = Inf)

cat("\nTOLERANCE PRINCIPLE: CLUSTER PATTERNS (REFERENCE DOWNSAMPLED LEXICON)\n")

if (!is.null(nouns_opp_down_single) && nrow(nouns_opp_down_single) > 0L) {
  cluster_tp_ie_ds <- nouns_opp_down_single |>
    mutate(cluster_simple = if_else(cluster %in% c("st", "sc", "ct"), cluster, NA_character_)) |>
    filter(!is.na(cluster_simple)) |>
    group_by(cluster_simple, opportunity) |>
    summarise(
      mutated = sum(mutation, na.rm = TRUE),
      non_mutated = sum(!mutation, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(type = paste0("<", cluster_simple, "> + <-", opportunity, "> plural (downsampled)")) |>
    tp_table(type, mutated, non_mutated)

  cluster_tp_comb_ds <- nouns_opp_down_single |>
    mutate(cluster_simple = if_else(cluster %in% c("st", "sc", "ct"), cluster, NA_character_)) |>
    filter(!is.na(cluster_simple)) |>
    group_by(cluster_simple) |>
    summarise(
      mutated = sum(mutation, na.rm = TRUE),
      non_mutated = sum(!mutation, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(type = paste0("<", cluster_simple, "> + <-i, -e> plural (downsampled)")) |>
    tp_table(type, mutated, non_mutated)

  cluster_tp_all_ds <- bind_rows(cluster_tp_ie_ds, cluster_tp_comb_ds) |> arrange(type)
  print(cluster_tp_all_ds, n = Inf, width = Inf)
} else {
  cat("No reference downsampled lexicon available; skipping cluster TP (downsampled).\n")
}

# =========================================================================
# Bayesian Tolerance Principle
# =========================================================================

cat("\nBAYESIAN TOLERANCE PRINCIPLE: SEGMENT × OPPORTUNITY (FULL DATA)\n")

tolerance_bayesian_results <- list()

for (seg in segments_of_interest) {
  for (opp in plural_opportunities) {
    seg_data <- filter(nouns_opp, stem_final == seg, opportunity == opp, !is.na(mutation))
    if (nrow(seg_data) < 5) next

    N <- nrow(seg_data)
    theta_N <- tp_threshold(N)
    tolerance_rate <- theta_N / N

    cat(sprintf("  Fitting <%s> + <-%s>...\n", seg, opp))

    invisible(capture.output(
      {
        model_tp <- brm(
          mutation ~ 1,
          data = seg_data,
          family = bernoulli(link = "logit"),
          prior = prior(normal(0, 1.5), class = "Intercept"),
          chains = 4, iter = 2000, warmup = 1000, cores = 4,
          backend = "cmdstanr",
          seed = 123,
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

    tolerance_bayesian_results[[paste(seg, opp, "full", sep = "_")]] <- data.frame(
      segment = seg,
      opportunity = opp,
      subset = "full",
      majority_mutates = majority,
      N = N,
      theta_N = round(theta_N, 2),
      tolerance_rate = round(tolerance_rate, 4),
      median_p_mutate = round(median(p_mutate), 4),
      ci_lower = round(quantile(p_mutate, 0.025), 4),
      ci_upper = round(quantile(p_mutate, 0.975), 4),
      prob_exceeds_tol = round(prob_intolerable, 4)
    )
  }
}

cat("\nBAYESIAN TOLERANCE PRINCIPLE: SEGMENT × OPPORTUNITY (REFERENCE DOWNSAMPLED)\n")

tolerance_bayesian_results_ds <- list()

if (!is.null(nouns_opp_down_single) && nrow(nouns_opp_down_single) > 0L) {
  for (seg in segments_of_interest) {
    for (opp in plural_opportunities) {
      seg_data_ds <- filter(nouns_opp_down_single, stem_final == seg, opportunity == opp, !is.na(mutation))
      if (nrow(seg_data_ds) < 5) next

      N_ds <- nrow(seg_data_ds)
      theta_N_ds <- tp_threshold(N_ds)
      tolerance_rate_ds <- theta_N_ds / N_ds

      cat(sprintf("  Fitting <%s> + <-%s> (downsampled)...\n", seg, opp))

      invisible(capture.output(
        {
          model_tp_ds <- brm(
            mutation ~ 1,
            data = seg_data_ds,
            family = bernoulli(link = "logit"),
            prior = prior(normal(0, 1.5), class = "Intercept"),
            chains = 4, iter = 2000, warmup = 1000, cores = 4,
            backend = "cmdstanr",
            seed = 456,
            refresh = 0,
            silent = 2
          )
        },
        type = "output"
      ))

      post_samples_ds <- as_draws_df(model_tp_ds)
      p_mutate_ds <- plogis(post_samples_ds$b_Intercept)
      p_nonmutate_ds <- 1 - p_mutate_ds

      sum_mut_ds <- sum(seg_data_ds$mutation, na.rm = TRUE)
      sum_non_ds <- sum(!seg_data_ds$mutation, na.rm = TRUE)
      majority_ds <- sum_mut_ds > sum_non_ds

      p_exception_ds <- if (majority_ds) p_nonmutate_ds else p_mutate_ds
      prob_intolerable_ds <- mean(p_exception_ds > tolerance_rate_ds)

      tolerance_bayesian_results_ds[[paste(seg, opp, "downsampled", sep = "_")]] <- data.frame(
        segment = seg,
        opportunity = opp,
        subset = "downsampled",
        majority_mutates = majority_ds,
        N = N_ds,
        theta_N = round(theta_N_ds, 2),
        tolerance_rate = round(tolerance_rate_ds, 4),
        median_p_mutate = round(median(p_mutate_ds), 4),
        ci_lower = round(quantile(p_mutate_ds, 0.025), 4),
        ci_upper = round(quantile(p_mutate_ds, 0.975), 4),
        prob_exceeds_tol = round(prob_intolerable_ds, 4)
      )
    }
  }
} else {
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

cat("\nSEGMENT CLASS COMPARISON: DORSAL VS CORONAL (I-DOMAIN)\n")

nouns_i_classified <- nouns_opp |>
  filter(opportunity == "i", !is.na(mutation)) |>
  mutate(
    segment_class = factor(
      if_else(stem_final %in% c("c", "g"), "dorsal", "coronal"),
      levels = c("dorsal", "coronal")
    ),
    suffix_group = case_when(
      lemma_suffix == "-ic" ~ "ic",
      lemma_suffix == "-ist" ~ "ist",
      lemma_suffix %in% c("-ică", "-ice") ~ "ica_ice",
      lemma_suffix == "none" ~ "none",
      TRUE ~ "other"
    ),
    suffix_group = factor(suffix_group, levels = c("none", "ic", "ist", "ica_ice", "other"))
  )

invisible(capture.output(
  {
    model_class <- brm(
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

print(summary(model_class))

draws <- as_draws_df(model_class)
beta_seg <- draws[["b_segment_classcoronal"]]
prob_dorsals_gt_coronals <- mean(beta_seg < 0)
or_coronal_vs_dorsal <- exp(beta_seg)
or_ci <- quantile(or_coronal_vs_dorsal, probs = c(0.025, 0.975))

cat("\nSEGMENT CLASS COMPARISON (DORSAL VS CORONAL)\n")
cat(sprintf("P(dorsals > coronals) = %.3f\n", prob_dorsals_gt_coronals))
cat(sprintf(
  "OR_coronal_vs_dorsal = %.3f [95%% CI: %.3f, %.3f]\n",
  median(or_coronal_vs_dorsal),
  or_ci[1],
  or_ci[2]
))
cat("  (OR < 1 ⇒ dorsals more likely to palatalize)\n")

cat("\nANALYSIS FINISHED\n")

sink()
