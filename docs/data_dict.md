This project investigates the palatalization of stem final consonants (dorsals and coronals) in Romanian. Below are the data fields and what entries they take and how to derive or verify them. Assume that your input is a csv (`romanian_lexicon_raw_dex.csv`) that contains the extracted fields listed below. The pipeline generates all the derived fields through a sequence of derivation functions organized in `src/romanian_processor_lib.py` and orchestrated by `scripts/romanian_processor_main.py`.

All text normalization uses `wiktionary_normalizer.py`, which provides:

* `normalize_orthography()` for lemmas/plurals/glosses (lowercasing, NFC, unifying comma- vs. cedilla-variants of <ș>/<ş> and <ț>/<ţ>, folding stressed vowels like <á, é, í> to their plain counterparts, whitespace cleanup).
* `normalize_ipa()` for IPA transcriptions (NFC, optional removal of stress marks, removal of parentheses, elimination of partial variants that start with "-" after a `|`, insertion of tie-bars in affricates tʃ/dʒ/ts/dz → t͡ʃ/d͡ʒ/t͡s/d͡z, and g→ɡ normalization).

These normalizers are applied to the raw Wiktionary harvest before the palatalization-specific derivations.

**Note on boolean storage**: Boolean-valued derived fields (`mutation`, `target_is_suffix`, `suffix_triggers_plural_mutation`, `is_true_exception`) are stored as lowercase strings (`"true"` / `"false"`) in the CSV output for consistency and cross-platform compatibility. R's `readr::read_csv()` automatically converts these to logical types (`TRUE`/`FALSE`) when loading the data.

### Core extracted fields

`lemma` is extracted from Wiktionary and is the orthographic form of the word. Every entry will trivially have this filled (after normalization via `normalize_orthography()`).

`gloss` is extracted from Wiktionary. Some entries have occasional gaps in the source.

`pos` is extracted from Wiktionary and is the part of speech of the lemma. It can only be `N` (Noun) or `ADJ` (Adjective). Every entry should have this filled.

`gender` is extracted from Wiktionary and is the grammatical gender of the lemma. Only `N` can have a gender entry. It can either be `MASC` (masculine), `FEM` (feminine), `NEUT` (Neuter), or `NA` for `ADJ`s.

### Stem structure: `stem_final` and `cluster`

`stem_final` is derived from the lemma orthography (after lowercasing and Unicode normalization, including unifying comma- vs. cedilla-variants of <ș>/<ş> and <ț>/<ţ>). It is the palatalization target consonant at the right edge of the nominal stem, restricted to the set {<c>, <g>, <t>, <d>, <s>, <z>}.

The logic is:

* If the normalized lemma ends in `<chi>`, `<che>`, `<ghi>`, or `<ghe>`, then:
  * `stem_final` is `<c>` for `<chi>/<che>` and `<g>` for `<ghi>/<ghe>`.
  * In these cases we do not strip the final vowel; the velar+front-vowel sequence is treated as part of the stem and is also recorded in the `cluster` field (see below).
* Otherwise, we first form a lemma stem by stripping a single final vowel (a, ă, â, e, i, î, o, u) if present and the lemma has length > 1. We then work in this lemma stem:
  * If the lemma stem ends in one of the clusters `<st>`, `<sc>`, or `<ct>`, we treat the cluster target as `stem_final`: `<s>` for `<st>/<sc>`, and `<t>` for `<ct>`. This matches the fact that `<s>` is the palatalizing consonant in `<st>/<sc> → <ști>/<ște>`, and `<t>` in `<ct> → <cți>/<cțe>`.
  * Otherwise, we search from right to left in the lemma stem until we find the *rightmost* occurrence of one of {<c>, <g>, <t>, <d>, <s>, <z>}. That consonant is `stem_final`.
* If no such consonant is found in the lemma stem (e.g. vowel-final stems without any of these obstruents), `stem_final` is left empty.

Examples:

* `<abulică>` → lemma stem `<abulic>`, `stem_final = c`.
* `<prost>` → lemma stem `<prost>`, final cluster `<st>`, `stem_final = s`.
* `<arhitect>` → lemma stem `<arhitect>`, final cluster `<ct>`, `stem_final = t`.
* `<păduche>` → lemma ends in `<che>`, so `stem_final = c` and the `<che>` sequence is preserved and recorded as a cluster.

`cluster` is derived from the lemma orthography and `stem_final`. It encodes the final orthographic sequence at the right edge of the lemma stem that is anchored at `stem_final` and is relevant for palatalization or non-derived-environment behavior. As above, we work with the lemma stem obtained by normalizing the lemma and stripping at most one final vowel, *except* that lemmas ending in `<chi>`, `<che>`, `<ghi>`, or `<ghe>` keep those endings intact.

* If the (possibly unstripped) lemma ends in one of the velar–front-vowel sequences `<chi>`, `<che>`, `<ghi>`, or `<ghe>`, then:
  * `cluster = "chi"` for lemmas ending in `<chi>`,
  * `cluster = "che"` for lemmas ending in `<che>`,
  * `cluster = "ghi"` for lemmas ending in `<ghi>`,
  * `cluster = "ghe"` for lemmas ending in `<ghe>`,
  * and `stem_final` is `<c>` for `<chi>/<che>` or `<g>` for `<ghi>/<ghe>` as specified above.
* Otherwise, we form the lemma stem by stripping a single final vowel (if present and the lemma has length > 1) and examine only the consonantal right edge:
  * If the lemma stem ends in one of the clusters `<st>`, `<sc>`, or `<ct>`, then:
    * `cluster = "st"` for stems ending in `<st>`,
    * `cluster = "sc"` for stems ending in `<sc>`,
    * `cluster = "ct"` for stems ending in `<ct>`,
    * and `stem_final` is `<s>` for `<st>/<sc>` or `<t>` for `<ct>`, as above.
* If the lemma stem does not end in `<chi>`, `<che>`, `<ghi>`, `<ghe>`, `<st>`, `<sc>`, or `<ct>`, then `cluster` is left empty.

In other words, `cluster` records the final palatalization-relevant orthographic package anchored at `stem_final`:

* multi-consonant clusters `<st>`, `<sc>`, `<ct>`, and
* velar–front-vowel complexes `<chi>`, `<che>`, `<ghi>`, `<ghe>`.

### Plurals and plural quality

`plural` is extracted from Wiktionary. Nouns are the main targets of the analysis, but adjectives can also have plural forms in the raw data; when present for `ADJ`, they are run through the same validation and alignment machinery as noun plurals.

To enforce that a plural is morphologically plausible for the lemma, we first normalize lemma and plural (lowercasing and Unicode normalization, including unifying comma- vs. cedilla-variants of <ș>/<ş> and <ț>/<ţ>). We then compute a multi-metric plausibility score that combines several string similarity measures to robustly detect spurious plural entries while retaining genuine irregular forms.

The plausibility score averages four complementary metrics:

1. sequence similarity (SequenceMatcher ratio),
2. longest common substring ratio,
3. character bigram overlap (Jaccard similarity),
4. common prefix/suffix ratio (longest shared prefix or suffix, normalized by word length).

Additionally, the length ratio between plural and lemma must fall between 0.5 and 2.0 to pass initial sanity checks.

Global thresholds are:

* a reject threshold, starting at 0.35 and later re-estimated as the 5th percentile of observed plausibility scores once we have at least 500 lemma–plural pairs;
* a borderline threshold, starting at 0.45 and later aligned to the 20th percentile, with a small cushion to keep it above the reject threshold.

The logic:

* If the length ratio is too extreme or the plausibility is below the reject threshold (or both the sequence similarity and LCS ratio are close to zero), the plural is rejected outright: `plural` is cleared and any `ipa_normalized_pl` present is also cleared.
* Once thresholds are calibrated, surviving pairs get a `plural_validity` label:
  * `borderline` if plausibility is between the reject and borderline thresholds;
  * `ok` otherwise.
* In downstream derivations (`mutation`, `orth_change`, `opportunity`), only entries with a non-empty `plural` participate in plural-based computations; rejected plurals behave as if missing.

Highly irregular but genuine plurals such as `<cască ~ căști>`, `<iască ~ iești>`, `<pască ~ păști>`, and `<sască ~ săști>` are retained by this algorithm, while misparsed labels like `<-Ă adj>` or unrelated forms like `<ruta ~ virnanți>` are filtered out.

### Mutation and orthographic alternations

`mutation` is derived from comparing lemmas with plural to their plurals at the position of the `stem_final` consonant (possibly as part of a final cluster `<st>`, `<sc>`, or `<ct>`). Mutation is boolean (`"True"` / `"False"`). It is defined for both `N` and `ADJ` with a plural that has passed the plural quality check; for other parts of speech or missing/filtered plurals the field is left empty.

Conceptually, the derivation proceeds in two stages.

(1) Anchor and align lemma and plural.

* Normalize lemma and plural as above.
* Compute `stem_final` and `cluster` from the lemma as specified in those fields. Let the target segment be `stem_final`, and, if `cluster` is non-empty, the target cluster `<st>`, `<sc>`, or `<ct>` at the right edge of the lemma stem.
* Perform a global edit-distance alignment (Needleman–Wunsch) between the normalized lemma and plural, yielding two aligned strings and a mapping from lemma indices to alignment indices. Using this mapping, identify the aligned positions corresponding to the target segment (or cluster) in the lemma.
* Define a change window as a contiguous span in the aligned strings that:
  * (i) covers the target segment (or cluster) in the lemma,
  * (ii) covers all positions where the aligned lemma and plural differ, and
  * (iii) includes one additional alignment column of right-hand context when available (this extra context character preserves phonological environment information).
* Let `lemma_sub` be the subsequence of the aligned lemma within this window with gaps removed, and `plural_sub` the corresponding subsequence of the aligned plural with gaps removed.
* If the normalized lemma and plural are identical, or if `lemma_sub` and `plural_sub` are identical, then no orthographic change is detected at the target site and `mutation = "False"`.

This same change window is also used to build `orth_change` (see below) and, when exported as a separate substring, can be reused in the plural environment computation. **Implementation note**: The plural-side change window substring is computed once in `derive_mutation_and_orth_change()` and cached in the row dictionary as `row["_plural_change_window"]` for efficient reuse by `derive_opportunity()`, then dropped from the final output to avoid clutter.

(2) Inspect the orthographic alternation at the target site.

Mutation is `True` if and only if the lemma–plural mapping at the `stem_final` slot (and any associated final cluster) matches one of the palatalizing patterns in `MUTATION_PATTERNS`, and the resulting plural-side pattern ends in a front vowel `<i>` or `<e>`.

The patterns currently recognized are:

* For `<c>`:
  * `c → ci`, `c → ce` (core c before front vowels),
  * `că → ci`, `că → ce` (suffixal -că / -ică into palatalized -ci/-ce).
* For `<g>`:
  * `g → gi`, `g → ge`,
  * `gă → gi`, `gă → ge`.
* For `<t>`:
  * `t → ți`, `t → țe` (t → ț + front vowel),
  * `ct → cți`, `ct → cțe` (cluster palatalization),
  * `te → ți`, `te → țe`, `tă → ți`, `tă → țe` (various suffixal and feminine patterns).
* For `<d>`:
  * `d → zi`, `d → ze`,
  * `de → zi`, `de → ze`,
  * `de → di` (rare allophonic palatalization),
  * `dă → zi`, `dă → ze`.
* For `<s>`:
  * `s → și`, `s → șe`,
  * `st → ști`, `st → ște`,
  * `sc → ști`, `sc → ște`,
  * `scă → ști`, `scă → ște` (including franciscă-type patterns).
* For `<z>`:
  * `z → ji`, `z → je` (obraz/arbuz-type patterns).

Formally, for each `stem_final` we require that `lemma_sub` ends in the lemma-side pattern and `plural_sub` ends in the corresponding plural-side pattern from `MUTATION_PATTERNS`, and that the plural-side pattern's final segment is `<i>` or `<e>`. If no such pattern matches, `mutation = "False"`.

Note: For `<z>`, only `z → ji/je` is counted as mutation. Plain `z → zi/ze` (i.e. a <z> plus `<i>/<e>` with no change to <j>) is *not* treated as a palatalization and is left as `mutation = "False"`.

`orth_change` is derived by comparing the lemma with the plural in an alignment anchored at the `stem_final` consonant (and, if present, its final cluster `<st>`, `<sc>`, or `<ct>`). It records the minimal orthographic change at the palatalization site as a string of the form `X→Y`, where `X` is the lemma-side substring and `Y` is the plural-side substring in the change window.

The algorithm:

1. Normalize lemma and plural; compute `stem_final` and `cluster` as above.
2. Compute the change window via alignment (as in the `mutation` step).
3. Construct `lemma_sub` and `plural_sub` by stripping gaps in the window.
4. If there is no difference, leave `orth_change` empty.
5. Otherwise, canonicalize the window into a short pattern anchored on the rightmost occurrence of `stem_final`:
   * take the tail of `lemma_sub` from the last occurrence of `stem_final` (or, if that fails, a short final slice), and the corresponding tail from `plural_sub`,
   * if these differ, set `orth_change = "lemma_tail→plural_tail"`.

For genuine palatalizations, the code overwrites this canonical window with a *clean abstract pattern* from `MUTATION_PATTERNS` itself (e.g. `"c→ci"`, `"st→ști"`, `"z→ji"`), so that `orth_change` is stable and matches `ORTH_TO_PALATAL_IPA`. For non-palatal changes at the stem edge (e.g. vowel alternations without consonant change), `orth_change` may still record a small local change window, but `mutation` remains `"False"`.

### Plural environment: `opportunity`

`opportunity` is derived from the lemma and plural orthography (after normalization). It encodes whether the plural morphology creates a front-vowel environment for the `stem_final` consonant at the right edge of the stem, and if so, which front-vowel plural class it belongs to. It is defined for both `N` and `ADJ` with a usable plural; for other parts of speech or missing/filtered plurals, `opportunity = "none"`.

The only possible values are:

* `"i"` – the plural places `stem_final` (or its `cluster`) directly before a plural `<i>` at the stem edge (C+`i`), i.e. an -i plural class (e.g. `<abulic ~ abulici>`, `<prost ~ proști>`, `<arhitect ~ arhitecți>`, `<acid ~ acizi>`).
* `"e"` – the plural places `stem_final` (or its `cluster`) directly before a plural `<e>` at the stem edge (C+`e`), i.e. an -e plural class (e.g. `<abac ~ abace>`, `<Drăgaică ~ Drăgaice>`).
* `"uri"` – the plural is a neuter -uri type where a `<u>` intervenes between `stem_final` and the plural `<i>` (C+`u…i`, e.g. `<afterpic ~ afterpicuri>`), so there is no direct C+front-vowel configuration at the stem edge. This class is structurally non-palatalizing for noun plurals.
* `"none"` – no usable palatalizing plural environment can be identified: either the plural is missing/filtered, or the lemma–plural pair does not realize `stem_final` (or its `cluster`) directly before `<i>` or `<e>`, and is not a standard -uri pattern.

Derivation proceeds in the following order:

1. If a local plural-side change window is available (e.g. `plural_change_window`, the plural-side counterpart to `lemma_sub`/`plural_sub`), the environment is read directly from that substring:
   * if it starts with `<u>` and contains `<i>`, set `opportunity = "uri"`;
   * else, the first occurrence of `<i>` or `<e>` inside that window determines `opportunity = "i"` or `"e"`;
   * if neither appears, keep `"none"`.
2. Otherwise, check `orth_change`:
   * if `orth_change` matches a canonical palatalization pattern in `ORTH_TO_PALATAL_IPA` (e.g. `"d→zi"`, `"t→ți"`, `"s→și"`, `"c→ci"`, `"g→gi"`, `"z→ji"`), read the environment from the ending of the plural side:
     * plural side ends in `<i>` → `"i"`;
     * plural side ends in `<e>` → `"e"`.
3. Fallback: alignment-based heuristics:
   * perform a global alignment between lemma and plural; locate the alignment position corresponding to the end of the target (`stem_final` or cluster);
   * inspect up to a few non-gap characters immediately following this position in the aligned plural string:
     * if the following substring starts with `<u>` and contains `<i>`, set `"uri"`;
     * else, if the first new vowel after the target is `<i>` contributed by the plural (not present in the lemma at that position), set `"i"`;
     * else, if it is `<e>` contributed by the plural, set `"e"`;
     * else leave `"none"`.

`"i"` and `"e"` mark C+front-vowel plural environments; `"uri"` marks the neuter -uri class; `"none"` covers all other cases. Note that `opportunity` does not distinguish derived vs. non-derived environments; that distinction is handled by the separate `nde_class` field.

### Palatalization in IPA: `palatal_consonant_pl`

`palatal_consonant_pl` is derived from the `orth_change` field by mapping canonical orthographic alternations to a single IPA symbol representing the palatal(ized) consonant in the plural. It is a compact way of encoding "what segment did this alternation produce?" in the plural.

The mapping is defined by `ORTH_TO_PALATAL_IPA` and currently includes:

* `c→ce`, `c→ci`, `că→ce`, `că→ci` → `t͡ʃ`
* `g→ge`, `g→gi`, `gă→ge`, `gă→gi` → `d͡ʒ`
* `t→țe`, `t→ți`, `ct→cți`, `ct→cțe`, `te→țe`, `te→ți`, `tă→țe`, `tă→ți` → `t͡s`
* `s→șe`, `s→și`, `st→ște`, `st→ști`, `sc→ște`, `sc→ști`, `scă→ște`, `scă→ști` → `ʃ` or `ʃt` depending on the cluster (simple s vs st/sc/scă)
* `d→ze`, `d→zi`, `de→zi`, `de→di`, `dă→zi` → broadly `z` (or `dʲ` in the special `de→di` mapping)
* `z→je`, `z→ji` → `ʒ`

Formally:

* If `orth_change` is empty, then `palatal_consonant_pl` is empty.
* If `orth_change` is in `ORTH_TO_PALATAL_IPA`, then `palatal_consonant_pl` is set to the corresponding IPA symbol; otherwise it is left empty.

In the core processor, `palatal_consonant_pl` is defined primarily from `orth_change`. **Backfill from IPA**: If `orth_change` is empty but `ipa_normalized_pl` explicitly shows a palatalized consonant (via regex matching t͡ʃ, d͡ʒ, ʃ, ʒ, t͡s, or z), this field may be backfilled from the IPA transcription to recover palatalization information when orthographic alternation detection fails. IPA-based checks are otherwise used for separate QC and phonetic analysis rather than for setting this field.

### IPA fields

These are conceptually as before, but now we can be explicit about normalization:

`ipa_raw_lemma` and `ipa_raw_pl` are scraped directly from Wiktionary (when present) and preserved as-is (including stress marks, syllable boundaries, and morphological segmentation markers `|` and `-`). They may be empty.

`ipa_normalized_lemma` and `ipa_normalized_pl` are derived primarily via `normalize_ipa()` applied to the raw IPA fields:

* stress marks (ˈ, ˌ, ', `) and syllable dots are stripped,
* optional segments in parentheses are removed,
* affricates are given tie-bars (tʃ → t͡ʃ, dʒ → d͡ʒ, ts → t͡s, dz → d͡z),
* Latin `g` is converted to IPA `ɡ`,
* variants separated by `|` are retained only if they are full forms (segments starting with "-" after a `|` are treated as partials and filtered out).

In addition, when no Wiktionary IPA is available for a form, the pipeline can fall back to the Romanian G2P (`to_ipa()` in `romanian_processor_lib.py`) to generate a broad, stressless IPA; the G2P output can then be normalized by the same IPA normalizer and stored in the normalized IPA fields.

The same G2P+normalizer combination is used explicitly for derived verbs and adjectives (see below).

### Denominal verbs: `derived_verbs`, `deriv_suffixes`, `ipa_derived_verbs`

`derived_verbs` is extracted from Wiktionary. It lists denominal or related verbs associated with the lemma, in Romanian orthography, as a `|`-separated list. Each token is a bare verb lemma (no glosses or extra markup), e.g.:

* `bloca|debloca`
* `claca|clăcui`
* `locui|înlocui`

This field is purely orthographic and is the input to both `deriv_suffixes` and `ipa_derived_verbs`.

`deriv_suffixes` is derived from `derived_verbs` by inspecting the right edge of each verb lemma. It encodes the infinitival verbal suffix of each derived verb, keeping the order and arity of `derived_verbs`.

* For each verb in `derived_verbs`, we normalize the orthography and then:
  * if it ends in `<ui>`, we label it `-ui`,
  * else if it ends in `<i>`, we label it `-i`,
  * else if it ends in `<a>`, we label it `-a`,
  * otherwise, the verb is treated as out of scope (e.g. infinitives in `-e`) and dropped.
* Only verbs ending in `-a`, `-i`, or `-ui` are retained; the rest are discarded for the purposes of this project.
* The resulting suffixes are concatenated with `|` in the same order as the surviving `derived_verbs`. If nothing survives, all three fields (`derived_verbs`, `deriv_suffixes`, `ipa_derived_verbs`) are cleared.

Example:

* `derived_verbs = "bloca|debloca"` → `deriv_suffixes = "-a|-a"`.
* `derived_verbs = "claca|clăcui"` → `deriv_suffixes = "-a|-i"`.

The only allowed atomic suffix values are `-a`, `-i`, and `-ui`.

`ipa_derived_verbs` is derived from `derived_verbs` using the Romanian G2P pipeline (`to_ipa()`) plus an injected IPA normalizer:

* Each surviving verb lemma is passed through `to_ipa()`, then normalized (stress removal, affricate tie-bar insertion, g→ɡ, etc.) via the `_ipa_normalizer` that the main script sets from `normalize_ipa()`.
* The resulting IPA strings are joined with `|` in the same order as the verbs in `derived_verbs`.

### Derived adjectives: `derived_adj`, `ipa_derived_adj`

`derived_adj` is extracted from Wiktionary. It lists adjectival derivatives associated with the lemma (e.g. suffixed adjectives in <-esc>, <-ic>, etc.) in Romanian orthography, usually as a single adjective lemma but in principle allowing `|`-separated multiple entries if needed. This field is orthographic only and serves as input for `ipa_derived_adj`.

`ipa_derived_adj` is derived from `derived_adj` using the same G2P + IPA-normalizer pipeline as `ipa_derived_verbs`:

* Each adjectival lemma in `derived_adj` (or each token if multiple adjectives are separated by `|`) is passed through `to_ipa()`, then normalized by `_ipa_normalizer`.
* The resulting forms are joined with `|` in the same order.
* `derived_adj` itself is normalized to a `|`-separated, whitespace-trimmed list of lemmas.

### Lemma-level suffixes: `lemma_suffix`, `target_is_suffix`, `suffix_triggers_plural_mutation`

`lemma_suffix` is derived from the lemma orthography (after lowercasing and Unicode NFC normalization, including unifying comma- vs. cedilla-variants of <ș>/<ş> and <ț>/<ţ>). It identifies whether the lemma ends in one of several productive derivational suffixes that are relevant for palatalization patterns. Romanian diacritics (ă, â, î, ș, ț) are preserved during matching to ensure accurate suffix detection and prevent false matches (e.g., distinguishing -ică from -ica).

The possible non-empty values are stored with a leading hyphen:

* `-ică` – if the normalized lemma ends in `<ică>`,
* `-iști` – if it ends in `<iști>`,
* `-ice` – if it ends in `<ice>`,
* `-ist` – if it ends in `<ist>`,
* `-esc` – if it ends in `<esc>`,
* `-ic` – if it ends in `<ic>` (but not in `<ică>`),
* `-el` – if it ends in `<el>`.

Suffix detection proceeds in this order, checking each suffix sequentially. The first matching suffix is assigned. Because `-ică` is checked before `-ic`, lemmas ending in `<ică>` (e.g. `<fizică>`) correctly match the `-ică` suffix pattern rather than being subsumed by the shorter `-ic` pattern.

`target_is_suffix` is derived from the lemma orthography, `stem_final`, `cluster`, and `lemma_suffix`. It is boolean (stored as `"True"`/`"False"`) and indicates whether the palatalization target consonant (stem_final and any associated final cluster) is wholly contained inside the derivational suffix identified in `lemma_suffix`, rather than in the lexical root.

Derivation:

* If `lemma_suffix` is empty or `stem_final` is empty, then `target_is_suffix = "False"`.
* Otherwise:
  * Let `suffix_text` be the suffix without the leading hyphen (e.g. `-ică` → `ică`).
  * Let `suffix_start` be the index where `suffix_text` begins in the lemma (`len(lemma) - len(suffix_text)`).
  * Define the target span:
    * for simple final consonants (empty `cluster`): the single index of `stem_final` in the lemma stem;
    * for non-empty clusters:
      * if `cluster` is one of `<chi>`, `<che>`, `<ghi>`, `<ghe>`, the target span is the last `len(cluster)` characters of the full lemma;
      * otherwise (clusters `<st>`, `<sc>`, `<ct>`), the span is over the lemma stem (with one final vowel stripped).
  * If the entire target span lies within `[suffix_start, len(lemma)]`, then `target_is_suffix = "True"`; otherwise `"False"`.

For example:

* `<fizic ~ fizici>`: `lemma_suffix = "-ic"`, `stem_final = c`, and `target_is_suffix = "True"`.
* `<lingvistic ~ lingvistice>`: `lemma_suffix = "-ic"`, but the palatalizing target is `<t>` in the root cluster, so `target_is_suffix = "False"`.

`suffix_triggers_plural_mutation` is derived from `lemma_suffix`, `pos`, `opportunity`, `mutation`, and `target_is_suffix`. It is boolean (`"True"` / `"False"`) and indicates whether a noun with a productive suffix actually shows palatalization at the stem edge in the plural paradigm, and whether that palatalization target lies within the suffix itself.

Derivation:

* Default: `suffix_triggers_plural_mutation = "False"`.
* It is set to `"True"` iff all of the following hold:
  * `pos == "N"`,
  * `opportunity ∈ {"i", "e"}` (the plural creates a C+front-vowel environment at `stem_final`/`cluster`),
  * `mutation == "True"` (palatalization occurs),
  * `target_is_suffix == "True"` (the palatalizing consonant is inside the suffix),
  * `lemma_suffix ∈ {"-ic", "-ist", "-esc", "-ică"}` (these suffixes are currently tracked as potentially palatalization-triggering).

In combination:

* `lemma_suffix` tells you which productive suffix (if any) the lemma bears.
* `target_is_suffix` tells you whether the palatalizing consonant belongs to that suffix.
* `suffix_triggers_plural_mutation` tells you whether the suffix itself undergoes palatalization in the plural paradigm.

### NDE patterns and exceptions: `nde_class`, `exception_reason`, `is_true_exception`

These fields carve out non-derived environment blocking (NDEB) patterns and then identify "true exceptions" to the productive palatalization rules.

`nde_class` is derived dynamically from the lemma and plural (for nouns only). It classifies lemmas that belong to non-derived environment blocking paradigms, where the C+front-vowel sequence at the stem edge is not a normal plural-driven environment and therefore is not informative for plural palatalization.

The only possible non-empty values are:

* `gimpe` – stems where the relevant C+<i>/<e> sequence (e.g. `<ci>`, `<ce>`, `<gi>`, `<ge>`) is tautomorphemic and underlyingly present in the stem, so the plural does not newly create a palatalizing environment. In the current code these are identified very conservatively as:
  * `pos == "N"`,
  * `stem_final ∈ {"c", "g"}`,
  * `cluster` is empty,
  * `lemma == plural` (segmentally identical sg~pl),
  * lemma ends in one of `ci`, `ce`, `gi`, or `ge`,
  * `mutation != "True"` (no palatalization).
* `ochi` – stems like `<ochi>`, with segmentally identical singular~plural (e.g. `<ochi ~ ochi>`), where the structure is ambiguous between /ok-i/ and /oki/. These are identified when:
  * `pos == "N"`,
  * `stem_final ∈ {"c", "g"}`,
  * `lemma == plural`,
  * `cluster ∈ {"chi", "ghi"}`,
  * `mutation != "True"`.
* `paduchi` – stems like `<păduche ~ păduchi>`, where plural morphology clearly attaches but palatalization under-applies in a special pattern. Identified when:
  * `pos == "N"`,
  * `stem_final ∈ {"c", "g"}`,
  * `cluster ∈ {"che", "ghe"}`,
  * the plural ends in the corresponding palatal sequence:
    * `<chi>` or `<ghi>`, or
    * `<chiuri>` / `<ghiuri>` (to cover forms like `<mănunche ~ mănunchiuri>`),
  * `mutation != "True"`.

If none of these conditions hold, `nde_class` is left empty.

All three NDE types are non-mutating dorsal nouns with both forms present (`pos == "N"`, `stem_final ∈ {c, g}`, `mutation == "False"` or empty, and a non-empty `plural`). They are designed to carve out paradigms where the presence of a C+front-vowel string at the stem edge is not informative about the status of the plural rule.

`exception_reason` is derived from the lemma, plural, `nde_class`, and mutation status. It gives a short label for non-palatalizing nouns that belong to recognized NDEB classes. It is only defined for `N` with a usable plural and `mutation != "True"`; otherwise it is left empty.

The possible values are:

* `nde:gimpe` – lemma belongs to the gimpe-type NDEB class.
* `nde:ochi` – lemma belongs to the ochi-type NDEB class.
* `nde:paduchi` – lemma belongs to the păduche-type NDEB class.
* empty – no specific exception type (the entry is either an unexplained exception or actually palatalizes).

Derivation:

* If `pos != "N"` or the plural is missing or `mutation == "True"`, then `exception_reason` is left empty.
* Otherwise (non-mutating noun with plural):
  * if `nde_class == "gimpe"`, set `exception_reason = "nde:gimpe"`;
  * else if `nde_class == "ochi"`, set `exception_reason = "nde:ochi"`;
  * else if `nde_class == "paduchi"`, set `exception_reason = "nde:paduchi"`;
  * else leave empty.

Note: Z-final stems with z→ji/je palatalization (e.g. `<obraz ~ obraji>`, `<arbuz ~ arbuji>`) are correctly classified as `mutation == "True"` by the mutation machinery and are *not* marked with `exception_reason`. Plain z→zi/ze patterns (where z does not change to j) are correctly classified as `mutation == "False"` and count as regular exceptions if no NDE pattern explains them.

Note: `exception_reason` is set for all NDE cases regardless of the `opportunity` value, because some NDE patterns (ochi, paduchi) may have `opportunity = "none"` due to alignment issues even though they morphologically involve front vowels.

`is_true_exception` is derived from `pos`, `opportunity`, `mutation`, and `nde_class`. It is boolean (`"True"` / `"False"`) and marks genuine unexplained failures of plural palatalization in derived front-vowel environments, after filtering out all known NDEB patterns.

* `is_true_exception = "True"` iff all of the following hold:
  * `pos == "N"`,
  * `opportunity ∈ {"i", "e"}` (the plural creates a C+`i` or C+`e` environment at `stem_final`/`cluster`),
  * `mutation == "False"` (no palatalization at that site),
  * `nde_class` is empty (not `"gimpe"`, `"ochi"`, or `"paduchi"`).
* In all other cases (no front-vowel opportunity, mutation `True`, or any NDE pattern), `is_true_exception = "False"`.

Entries with `mutation == "True"` are normal undergoers (including z→ji/je) and are never marked as exceptions. All NDE cases (gimpe, ochi, păduche) are excluded from `is_true_exception = "True"` but remain identifiable via `nde_class` and `exception_reason` for separate analysis.

### Other fields

`etym_lang` is extracted from Wiktionary and then normalized via `normalize_orthography()`.

`source` is extracted from Wiktionary. It is the URL the lemma entry came from. It should always be present.

`notes` is derived from Wiktionary and DEX online data. It records which entries required plural or gender confirmation from DEX online due to missing or ambiguous Wiktionary data. If no such confirmation was needed, the field is left empty.