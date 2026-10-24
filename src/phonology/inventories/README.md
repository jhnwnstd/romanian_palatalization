# Feature inventories

This directory is package data for `phonology`. It holds one feature
inventory per language, `<name>_features.json` (`italian_features.json`,
`czech_features.json`): the surface segments of the language, each with
a value for every feature. An analysis loads its inventory by name and
declares its own layers on top (`FeaturePatch`, `UnderspecifiedSegment`;
see `phonology/inventory.py`).

```python
from phonology.inventory import inventory_path, load_inventory_file

inventory = load_inventory_file("italian", PATCHES, UNDERSPEC)  # FeatureInventory
json_path = inventory_path("italian")  # for AnalysisProfile.inventory_json
```

## Lookup

`inventory_path(name)` returns the first of these that is a file:

1. `src/phonology/inventories/<name>_features.json`;
2. for `romanian` only, `romanian_features.json` at the repository root.

The Romanian inventory predates this directory. It stays at the
repository root, where the Romanian analysis, the tests and the docs
refer to it, so do not move or copy it here. A copy here would take
precedence over the root file.

`inventory_path(name, directory=d)` and `load_inventory_file(...,
directory=d)` look only in `d`. Tests use this for temporary inventories.
A name must be a lowercase identifier (`[a-z][a-z0-9_]*`). Any other
name raises `ValueError`. A missing file raises `FileNotFoundError`,
which lists every path that was tried.

## File name and git tracking

`.gitignore` ignores `*.json` but re-includes `!*_features.json`. Only
files with the `_features.json` suffix are tracked. Before adding a new
file, run:

```sh
git check-ignore -v src/phonology/inventories/czech_features.json
```

The rule it reports must be the negation `!*_features.json`, not
`*.json`.

## Format

The format follows `romanian_features.json`, the reference inventory.
The file is UTF-8, with IPA written literally (no `\u` escapes) and
2-space indentation. It holds one object with three keys, in this order:

```json
{
  "metadata": {
    "name": "Romanian phoneme inventory",
    "notes": "Core contrasts: seven oral vowels and ..."
  },
  "features": ["Syllabic", "Consonantal", "...", "Long", "Stress"],
  "segments": {
    "i": {"Syllabic": "+", "Consonantal": "-", "...": "0", "Stress": "-"},
    "t͡ʃ": {"Syllabic": "-", "Consonantal": "+", "...": "0", "Stress": "-"}
  }
}
```

`FeatureInventory.load` reads only `features` and `segments`. The
`metadata` key is documentation.

### `metadata`

- `name`: a human-readable name for the inventory.
- `notes`: free text covering the core contrasts and every convention a
  reader could not guess. The Romanian notes, for example, say that
  velars are `0` for the tongue-body features and that `Tense` is `0`
  on vowels.

### `features`

Every inventory lists the same 27 features, in this order. All
languages therefore share one feature geometry, and natural classes
written for one analysis mean the same thing in another.

| # | Feature | # | Feature | # | Feature |
|---|---------|---|---------|---|---------|
| 1 | `Syllabic` | 10 | `Nasal` | 19 | `DORSAL` |
| 2 | `Consonantal` | 11 | `Lateral` | 20 | `High` |
| 3 | `Sonorant` | 12 | `Trill` | 21 | `Low` |
| 4 | `Approximant` | 13 | `Strident` | 22 | `Front` |
| 5 | `Voice` | 14 | `LABIAL` | 23 | `Back` |
| 6 | `SpreadGl` | 15 | `Labiodental` | 24 | `Round` |
| 7 | `ConstrGl` | 16 | `CORONAL` | 25 | `Tense` |
| 8 | `Continuant` | 17 | `Anterior` | 26 | `Long` |
| 9 | `DelRel` | 18 | `Distributed` | 27 | `Stress` |

Copy the list verbatim from `romanian_features.json`. If a language
needs a new feature, add it to every inventory, because the analyses
compare across languages.

### `segments`

Each key is a segment label. Each value is an object with every feature
in `features` order, and each feature value is the string `"+"`, `"-"`
or `"0"`. `"0"` means unspecified: in unification it is a wildcard,
which a rule may fill (see `phonology/segments.py`).

Labels:

- A label is exactly the token that `tokenize_ipa` produces for the
  segment under the analysis's `TokenizerConfig` (`profile.tokenizer`).
  A segment the tokenizer never emits cannot be reached from IPA.
  - Affricates take the tie bar U+0361: `t͡s`, `t͡ʃ`, `d͡ʒ`.
  - Multi-codepoint segments such as Czech `r̝` and `ɪ̯` need
    `keep_combining=True`. Otherwise the diacritic becomes a token of
    its own.
  - Long segments get entries of their own (`aː` with `Long` `"+"`)
    and need `length=LengthHandling.ATTACH`. Romanian drops `ː`, so
    its inventory has no long segments and `Long` is `"-"` throughout.
- Labels are NFC.
- Labels must not collide with an analysis's `UnderspecifiedSegment`
  labels. Romanian declares `K`, `G`, `S`, `T` and `D`.

Contents:

- List surface segments only. Underspecified segments (`/K/`, `/T/`)
  and analysis-specific feature overrides belong in the analysis module,
  as `UnderspecifiedSegment` and `FeaturePatch`. That keeps this file an
  auditable base that analyses do not edit.

Conventions in `romanian_features.json`; follow them unless the notes
say otherwise:

- Every segment is `+` or `-` for each articulator feature (`LABIAL`,
  `CORONAL`, `DORSAL`; capitalised).
- `Anterior` and `Distributed` are `0` unless the segment is
  `CORONAL` `+`.
- `Labiodental` is `0` except on `LABIAL` consonants.
- The tongue-body features (`High`, `Low`, `Front`, `Back`, `Round`)
  are reserved for vocoids (vowels and glides). Velars are `DORSAL`
  `+` and `0` for all five.
- `DelRel` is `-` on oral stops, `+` on affricates and `0` elsewhere
  (nasals included).
- `Tense` is `0` in a language without a tense/lax contrast. `Long`
  and `Stress` are `-` where length and stress are not segmental.

## Checklist for a new inventory

1. Name the file `<name>_features.json` and confirm it is tracked
   (`git check-ignore -v`).
2. Copy `features` from `romanian_features.json` verbatim.
3. Add one entry per surface segment, labelled as the analysis's
   tokenizer emits it, with all 27 features in order and values
   `"+"`, `"-"` or `"0"`.
4. Record the conventions in `metadata.notes`.
5. Run `tests/test_phonology_inventory.py`, which checks every
   `*_features.json` here against this format.
