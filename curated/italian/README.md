# Curated Italian decisions

Hand decisions about individual Italian lexemes. The pipeline only reads
these files (`RunContext` refuses to write under `curated/`); the
`features` stage loads them with `lexpipe.curated` and stops when a
row's guard columns no longer match the current data (a new kaikki dump
changed or renumbered the lexeme). Re-check such a row by hand, then
update its guard columns (`docs/runbooks.md`, "Moving to a new kaikki
dump", step 6).

Every file is keyed by `lexeme_id` (one row per lexeme). Columns whose
name starts with `guard_` record what the curator saw; the others are
the decision.

## plural_adjudications.csv

Tier-2 plural QC: lexemes flagged in
`analysis/italian/italian_manual_worklist.csv`, looked up **by hand** in
Treccani or De Mauro. Never automate these lookups: Treccani opts out of
text and data mining and plants trap URLs. Store a coded value and a
citation, never copied dictionary text (R9).

| column | content |
|---|---|
| `lexeme_id` | the lexeme |
| `guard_word`, `guard_pos`, `guard_gender` | the lexicon's `word`, `pos`, `gender` |
| `guard_plurals` | the Wiktionary forms of `cell` as listed in `italian_plural_forms.csv`, `\|`-joined in source order |
| `cell` | `pl`, `m.pl` or `f.pl` (empty: every cell) |
| `adjudicated_plurals` | the forms the dictionary accepts, `\|`-joined; other forms of the cell are left out of the cell outcome (`adjudicated_out`) |
| `decision` | coded value: `velar_only`, `palatal_only`, `velar_preferred`, `palatal_preferred`, `both`, `invariable`, `other` |
| `dictionary` | `treccani` or `demauro` |
| `citation` | the entry URL and the date it was consulted |
| `adjudicator`, `date`, `note` | who decided, when, and why (free text, no copied entry text) |

## suffixoid_classes.csv

Decisions about -logo/-fago/-urgo lexemes (R27) where the shared
classifiers of `lexpipe.suffixoids` (the same for Romanian and Czech:
`curated/romanian_compare/` and `curated/czech/` hold their lists) are
wrong: whether the ending is the learned suffixoid, and the
person/thing class of the first gloss. Formerly `person_thing.csv`
(values `person`, `thing`, `not_suffixoid`); its three decisions moved
here unchanged.

| column | content |
|---|---|
| `lexeme_id` | the lexeme |
| `learned_suffixoid` | `False` for a lookalike (the ending is not the suffixoid: *spurgo* < *spurgare*), `True` to force the learned reading, empty: no decision |
| `semantic_class` | `person`, `animate`, `thing` or `unknown` (the shared classes), empty: no decision |
| `decided_by`, `note` | who decided and why |
| `guard_word`, `guard_pos`, `guard_first_gloss` | the lexicon's `word`, `pos`, `first_gloss` |
