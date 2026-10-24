# Curated Czech decisions

Hand decisions about individual Czech lexemes, read (never written) by
the pipeline (`docs/multilingual_design.md` R6, R9).

## `cell_adjudications.csv`

Starts as a header-only template. One row per decided cell, keyed by
(`lexeme_id`, `cell`) — `lexeme_id` from `data/czech/czech_lexicon_full.csv`
(or the tracked `czech_lexicon_raw.csv`); the column contract is
documented in `src/lexpipe/langs/czech/adjudication.py`. (The first
template was wide, one row per lexeme with `confirmed_<CELL>` and
`guard_forms_<CELL>` columns; it held no decisions when it was converted
to this long form.)

How to fill a row:

1. Take a row of `analysis/czech/czech_manual_worklist.csv` (cells where
   the Wiktionary forms and the automated MorfFlex QC disagree).
2. Open its `ijp_url` **in a browser, by hand**. Do not script it: IJP's
   robots.txt disallows `/?slovo=`, and copying its content needs the
   consent of ÚJČ AV ČR.
3. Record a coded decision only: `cell` (NOM.PL, VOC.SG, LOC.SG, DAT.SG,
   LOC.PL or NOM.PL.M.AN) and in `confirmed` the Wiktionary forms of
   that cell that IJP gives as standard (`|`-joined, in the worklist's
   order), or `-` if it gives none of them. Never paste IJP text.
4. Fill `citation` (source and date) and the guard columns with what you
   saw: `guard_lemma`, `guard_gender`, `guard_animacy` and `guard_forms`
   (the cell's forms) exactly as in the worklist / lexicon.

The features stage checks the guards on every run; a decision whose
lexeme or cell changed in a later dump stops the run until it is
re-checked.

## `suffixoid_classes.csv`

Decisions about -log/-urg/-graf/-fil nouns (R27) where the shared
classifiers of `lexpipe.suffixoids` (the same for Italian and Romanian)
are wrong: whether the ending is the learned suffixoid, and the
person/thing class of the first gloss. Keyed by `lexeme_id`; guard
columns `guard_lemma`, `guard_readings` and `guard_first_gloss`.

| column | content |
|---|---|
| `lexeme_id` | the lexeme |
| `learned_suffixoid` | `True`, `False` (a lookalike) or empty: no decision |
| `semantic_class` | `person`, `animate`, `thing`, `unknown`, or empty: no decision |
| `decided_by`, `note` | who decided and why |
| `guard_lemma`, `guard_readings`, `guard_first_gloss` | the lexicon's `lemma`, `readings`, `first_gloss` |
