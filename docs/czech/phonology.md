# Czech palatalization: the Logical Phonology analysis

A provisional analysis of the Czech small (second) palatalization in
nominal inflection, in the framework of Gorman & Winstead's Romanian
analysis (`latex.tex`: underspecified targets, prespecified segments
inalterable by unification failure, feature-filling rules in Logical
Phonology). It is checked against every velar, r and t/d/n cell form of
the Czech data layer (`docs/czech/README.md`).

| what | where |
|---|---|
| feature inventory | `src/phonology/inventories/czech_features.json` |
| analysis (patches, underspecified segments, palatalizers, rules, variants, row loader, UR builder, transcription) | `src/phonology/analyses/czech_palatalization.py` |
| tests | `tests/test_phonology_czech.py` |
| validator reports, one per variant | `analysis/czech/czech_rule_validation[_<variant>].txt` |
| mismatch tables, one per variant | `analysis/czech/czech_rule_mismatches[_<variant>].csv` |
| rates by class × cell, by stem final, by cluster, and the Czech mismatch classes, all variants | `analysis/czech/czech_rule_validation_summary.txt` |

```
python scripts/validate_lang_rules.py --lang cs                      # lexical (default)
python scripts/validate_lang_rules.py --lang cs --variant all_underspecified
python scripts/validate_lang_rules.py --lang cs --variant lexical_loo
python scripts/validate_lang_rules.py --lang cs --variant spellout
PYTHONPATH=src python -m phonology.analyses.czech_palatalization     # the summary
PYTHONPATH=src python -m phonology.analyses.czech_palatalization --trace kluk český
```

The inputs are the generated `data/czech/czech_cells.csv` and
`czech_lexicon_complete.csv` (data stages `lexicon` … `features`). The
analysis reads only their Wiktionary-derived columns and evidence at
τ = dictionary, never a corpus count, so its outputs carry Wiktionary
forms and transcriptions only (CC BY-SA 4.0, R9).

## 1. The inventory

`czech_features.json` has the 27 features of `romanian_features.json`
in the same order and 39 segments, labelled as cs-IPA (the
`{{cs-IPA}}` template of English Wiktionary, the convention of the
lemma transcriptions in the data) writes them. Every bundle is unique.

**Shared with Romanian, copied exactly:** a o u p b t d k ɡ t͡s t͡ʃ d͡ʒ
f v s z ʃ ʒ m n l r j.

**New segments and why:**

| segment | bundle | justification |
|---|---|---|
| ɪ, ɛ | Romanian i, e | cs-IPA writes short i/y as ɪ and e/ě as ɛ; Czech has no tense/lax contrast, so the bundles are Romanian i and e (Tense 0) |
| iː ɛː aː oː uː | the short vowel with Long + | vowel length is contrastive (the tokenizer attaches ː) |
| r̩, l̩ | r, l with Syllabic + | syllabic liquids (vlk, krk, bratr) |
| u̯ | Romanian w | the offglide of ou̯ au̯ ɛu̯ |
| c ɟ ɲ (ť ď ň) | t d n with Anterior −, Distributed + | see decision (a) |
| x (ch) | k with Continuant +, DelRel 0 | the voiceless velar fricative |
| ɦ (h) | x with Voice + | see decision (c) |
| r̝ (ř) | r with Sonorant −, Approximant −, Anterior −, Distributed + | see decision (b) |

**Decisions** (the open questions of the draft inventory,
`data/research_v2/czech_features_draft.json`; the criterion was that
the small palatalization should be feature-filling wherever possible):

- **(a) ť ď ň are plain [−anterior, +distributed] coronals**, not
  coronal + dorsal complexes (the draft's DORSAL +, High +, Front +).
  Then t > ť, d > ď and n > ň fill exactly two features (Anterior,
  Distributed) on /T D N/, and those two features are the floating
  palatalizer's own place (CVW's element I). The complex encoding is
  feature-filling too, but it needs DORSAL and five tongue-body
  features cleared on every coronal target and supplied by the rule,
  and it would put tongue-body values on consonants, which the
  Romanian conventions reserve for vocoids. Italian ɲ is coded the same
  way. The phonetic affrication of the palatal stops (JIPA 2012) is
  left to phonetic interpretation.
- **(b) ř is an obstruent with the palatal place:** r with
  [−sonorant, −approximant] (it devoices and takes part in voicing
  assimilation: JIPA 2012) and [−anterior, +distributed]. The draft's
  alternative kept ř [+anterior] and changed only manner (Strident +),
  which makes r > ř a change of manner that no other target shares;
  with the place coding, r > ř is the same place fill as t > ť and
  n > ň (CVW: r ∪ I = ř, t ∪ I = ť, n ∪ I = ň), and the obstruency of
  ř is a redundancy filled by a default (a [−anterior] trill is an
  obstruent). Trill + and Lateral − keep it apart from ʒ; Strident is
  0, as on r.
- **(c) ɦ is a dorsal: the voiced x.** A laryngeal ɦ (the Romanian h
  with Voice +, the draft's choice) would make h > z a change of
  Consonantal and of the articulator, not a fill. Phonologically ɦ is
  the voicing partner of x (sníh [sɲiːx]; cs-IPA writes the voiced ch
  of dechberoucí as ɦ), and it alternates with z/ʒ exactly like ɡ, so
  /Ɦ/ is placeless like /K G X/ and h > z is the same fill as g > z.
  Its glottal articulation is phonetic. (Poĺomská 2025, fn. 2, leaves
  the question open.)
- **(d) Long vowels are separate segments with Long +**; Stress is −
  throughout (stress is fixed and initial, so it is stripped before
  comparison) and Tense 0, as in Romanian.

**Not listed** (folded by the tokenizer onto their phonemes): the
allophones r̝̊ (ř after a voiceless obstruent or finally), ŋ and ɱ (n,
m before velars and labiodentals), ɣ (= ɦ under decision (c)) and d͡z
(voiced c), and the glottal stop (skipped). Variant symbols of the
cs.wiktionary transcriptions (ɔ ʊ ɔʊ̯ e) and of loans (g h w) fold too.

**Patches** (`PATCHES`, in the analysis, not the JSON): k ɡ x ɦ and
t d c ɟ are [−strident]. As in Romanian (t d), the prespecified value is
what makes a segment inalterable to assibilation; on the velars it
also lets the -ský/-cký /Ḱ/ keep the value while /K/ clears it.

## 2. The analysis

### The trigger is a suffix-borne palatalizer

The vowel cannot be the trigger. -i/-í, -e/-ě and -ích palatalize;
-y/-ý [ɪ iː] (potoky, mladý), -é (občané), INS.SG -em, LOC.PL -ech and
the -e of the coronal vocative (studente) do not, although they share
the vowels. Following Cavirani & Vanden Wyngaerd (2026, CVW) and Ketner
(2005), the palatalizer is a floating feature bundle that belongs to the
P endings. The UR builder inserts it as a placeholder segment between
the stem and every P ending, and only there. The rules locate it with a
broad, next-segment SEARCH, so only the consonant right before it is
affected, and a deletion rule removes it afterwards.

| placeholder | bundle | carried by |
|---|---|---|
| PAL1 | [−syll, −cons, +cor, −ant, +dist, −strid] | CVW's PAL1 'I': DAT/LOC.SG -e/-ě, NOM.PL -i/-í, LOC.PL -ích, adjectival NOM.PL.M.AN -í, the Cr vocative -e |
| PAL3 | [−syll, −cons, +cor, −ant, +dist, +strid] | CVW's PAL3 'H.I': the vocative -e of velar stems (bože, hříšníče) |

No real segment is [−consonantal, +coronal] (j and u̯ are CORONAL −), so
the natural class [−cons, +cor] picks out the placeholders.

### Underspecified targets

| label | base | cleared | surfaces as |
|---|---|---|---|
| K | k | CORONAL DORSAL Strident DelRel | t͡s / t͡ʃ, or k |
| G | ɡ | CORONAL DORSAL Strident DelRel Continuant | z / ʒ, or ɡ |
| X | x | CORONAL DORSAL Strident | ʃ, or x |
| Ɦ | ɦ | CORONAL DORSAL Strident | z / ʒ, or ɦ |
| R | r | Anterior Distributed Sonorant Approximant | r̝, or r |
| T D N | t d n | Anterior Distributed | c ɟ ɲ, or t d n |
| Ḱ | k | CORONAL DORSAL (Strident − kept) | c (ť), or k |
| S, TS | s, t͡s | Anterior Distributed | ʃ, t͡ʃ, or s, t͡s |

The velars are placeless (CVW: velars have no place element; the
palatalizer gives them one), so every velar outcome is a fill. /G/ is
also open for Continuant, the device of Romanian /D/: the palatalized
stop is a fricative (lize), and Czech has no *d͡z output. /Ḱ S TS/ are
the -sk-/-ck- of the adjectival suffix -ský/-cký (and adjectival nouns:
hospodský): /Ḱ/ keeps [−strident], so it cannot be assibilated and
takes the palatalizer's own place instead.

### Rules

All palatalization rules are unification (feature-filling); nothing is
feature-changing. A prespecified segment resists by unification failure,
since the supply contradicts one of its values.

| # | rule | LP statement | effect |
|---|---|---|---|
| 1 | big-pal | [−son, +cons] ⊔ {+cor, −dor, +strid, −ant, +dist} / _ PAL3 | k > t͡ʃ, g ɦ > ʒ, x > ʃ |
| 2 | x-pal | [−son, +cons, +cont, −voice] ⊔ {+cor, −dor, +strid, −ant, +dist} / _ PAL | ch > š (Češi, mouše) |
| 3 | velar-pal | [−son, +cons] ⊔ {+cor, −dor, +strid, +ant, −dist} / _ PAL1 | k > t͡s, g ɦ > z (kluci, lize, knize) |
| 4 | coronal-pal | [+cons] ⊔ {+cor, −dor, −ant, +dist} / _ PAL | t d n r > ť ď ň ř; /Ḱ/ > ť |
| 5 | s-pal | [+cons, +cor, +strid] ⊔ {−ant, +dist} / _ [+cons, +cor, −ant, +dist] | /S TS/ > š č before the derived ť |
| 6 | pal-deletion | [−syll, −cons, +cor] → ∅ | the placeholder |
| 7 | dorsal-default | [−son, +cons] ⊔ {−cor, +dor, −strid} | unpalatalized /K G X Ɦ Ḱ/ |
| 8 | dorsal-stop-default | [−son, +dor] ⊔ {−cont, −delrel} | k, ɡ |
| 9 | sibilant-continuant | [+cor, +strid] ⊔ {+cont} | palatalized /G/ is z, ʒ |
| 10 | affricate-release | [+cor, +strid, −cont] ⊔ {+delrel} | palatalized /K/ is t͡s, t͡ʃ |
| 11 | coronal-place-default | [+cor] ⊔ {+ant, −dist} | unpalatalized /T D N R S TS/ |
| 12 | trill-obstruent | [+trill, −ant] ⊔ {−son, −approx} | ř |
| 13 | trill-sonorant | [+trill] ⊔ {+son, +approx} | r |

The order does the work that Pāṇinian specificity does elsewhere. x-pal
precedes velar-pal, which then fails on the [−anterior] /X/, so the
voiced /Ɦ/, which x-pal does not target, comes out z. velar-pal precedes
coronal-pal, which then fails on the [+anterior] /K G Ɦ/. coronal-pal
feeds s-pal. On the full data every derived segment is exactly an
inventory bundle (no feature left open and no reliance on the
projection's tolerance; checked over all 48,143 rows and four variants).

**What is and is not spreading.** For the coronals (t d n r > ť ď ň ř)
and for /Ḱ/, the supply is the palatalizer's own place: the rule unifies
the trigger's features into the target. For the velars it is not.
k g h > c z z need the alveolar sibilant place [+ant, −dist], which the
palatalizer does not carry. In CVW this is phonetic interpretation (a
placeless stop plus I is pronounced [ʦ]); in a feature system it is a
supply that does not come from the trigger. The h/ch asymmetry (h > z
but ch > š) has to be stated as well: it is x-pal's target class
[+cont, −voice]. The small palatalization is feature-filling throughout
but is a copy of the trigger only for the coronals.

### Derivations

| | kluk+i | kniha+e | Čech+i | doktor+i | český+í | most+ě | bůh VOC | park+u |
|---|---|---|---|---|---|---|---|---|
| UR | kluK·PAL1·ɪ | kɲɪꞪ·PAL1·ɛ | t͡ʃɛX·PAL1·ɪ | doktoR·PAL1·ɪ | t͡ʃɛSḰ·PAL1·iː | mosT·PAL1·ɛ | boꞪ·PAL3·ɛ | parKu |
| big-pal | | | | | | | boʒ | |
| x-pal | | | t͡ʃɛʃ | | | | | |
| velar-pal | klut͡s | kɲɪz | | | | | | |
| coronal-pal | | | | doktor̝ | t͡ʃɛsc | mosc | | |
| s-pal | | | | | t͡ʃɛʃc | | | |
| defaults | affricate-release | | | trill-obstruent | | | | dorsal-default, dorsal-stop-default |
| SR | klut͡sɪ | kɲɪzɛ | t͡ʃɛʃɪ | doktor̝ɪ | t͡ʃɛʃciː | moscɛ | boʒɛ | parku |

(Only the fired rules are shown; the defaults change no symbol in
these columns, they fill the release, the trill's manner and the
velar's place.) The golden set of the tests also has matka+e →
/matt͡sɛ/ (cs-IPA writes [mat͡sɛ], with t + t͡s degeminated; the
comparison degeminates), vrah+i → vrazɪ, biolog+ové → bɪjoloɡovɛː (no
palatalizer, so the /G/ defaults to ɡ), potok+e → potot͡sɛ, moucha+e →
mou̯ʃɛ (cs-IPA writes the diphthong ou̯) and sestra+e → sɛstr̝ɛ (cs-IPA:
[sɛstr̝̊ɛ], with ř devoiced after t).

### Clusters

- **Root sk, st, zd: only the final consonant changes.** deska → desce
  [st͡s], Etrusk → Etruscích, bělásek → bělásci, most → mostě [sc],
  hvězda → hvězdě [zɟ]. The s/z is prespecified.
- **The suffix -ský/-cký: šť čť.** český → čeští [ʃc], německý →
  němečtí [t͡ʃc]. /Ḱ/ becomes ť (coronal-pal), and s-pal then makes the
  suffix's /S TS/ postalveolar before it. This is the Romanian sk
  analysis (latex.tex: S-palatalization fed by the dorsal rule) without
  the bleed, since Czech keeps the palatal stop. Which clusters are
  /SḰ/ is read off the lemma (-ský/-cký), the same in every variant.
  The other possible account, a separate adjectival palatalizer, is
  contradicted by noun NOM.PL -i, which gives sc (bělásci), and by
  vojsko → vojště, where a case ending gives šť.
- **ck in roots** (fracek → fracci [t͡st͡s]) and lexical č + k
  (blizoučký → blizoučcí [t͡ʃt͡s]): the c/č is prespecified and k > t͡s.

## 3. Rows, transcriptions and the comparison

**Rows** (`iter_rows`): every primary (own-edition), aligned,
not-rejected row of `czech_cells.csv` with a velar, r or
coronal_inferred target and either a P ending or a hard-paradigm N
ending (the soft LOC.SG -i of kmen, kmeni, is left out), spelled with
letters only (numeral compounds such as 10letý are left out), one row
per lexeme × cell × form. That gives 48,143 rows of 21,337 lexemes
(25,714 P and 22,429 N): velar 20,352 noun and 2,813 adjective rows,
r 4,590 and 44, coronal_inferred 17,023 and 3,321.

**Attested transcription.** The data has no transcriptions of inflected
forms, so the attested IPA is the spelling transcribed by
`transcribe()`, a re-implementation of the cs-IPA conventions: i y [ɪ],
í ý [iː], ou [ou̯]; ť ď ň before ě i í; bě pě vě fě [Cjɛ], mě [mɲɛ]; a
hiatus j after i (biologie [bɪjoloɡɪjɛ]); final devoicing and regressive
voicing assimilation (v undergoes but does not trigger; ř triggers only
when voiceless); degemination (but not of affricates: hrdliččí); ŋ
before velars; syllabic r l; and ř devoiced after a voiceless obstruent.
Against the English Wiktionary lemma transcriptions it gives exactly the
same string for 25,896 of 28,031 one-word lemmas (92.4%). Most of the
rest are loans whose di ti ni are hard by a template hint (diplomat
[dɪplomat]); with the hard letters of each lemma's own IPA applied, as
the analysis does inside the stem (`hard_coronals`), the agreement is
27,487 of 28,031 (98.1%). The remainder are loan respellings (byte
[bajt]) and cs-IPA quirks (srnka [srn̩ka]).

**UR.** The form itself, with the stem-final change undone (kluci →
kluk-i, čeští → česk-í), transcribed the same way, except that the
ending's effect on the stem edge (ť ď ň before ě i í) and degemination
at the edge (měkký: k + /K/) are suppressed. Then the target is
swapped for its underspecified twin (per variant), and the ending's
palatalizer goes in before the ending. As in Romanian, spell-out is not
modelled: the rules are tested given the attested ending. Everything but
the stem edge is the attested form's own transcription, so the IPA match
tests the derived stem edge against the spelling.

**Comparison** (`CZECH_NORMALISERS`): stress, syllable dots, the glottal
stop, the allophones r̝̊ ŋ ɱ ɣ d͡z and degemination (t + affricate,
identical consonants) are ignored; length is not. **Caveat:** t/d/n rows
(coronal_inferred) compare against a change inferred from the spelling
(the data layer's convention, ti = [cɪ] at the stem edge). Their match
only says that the rules agree with that inference, and for loans
(studenti, diplomati) the inference itself is unverified (no orthoepic
source; design §10).

## 4. Variants: which targets are underspecified

| variant | a target is underspecified iff | what the boolean match tests |
|---|---|---|
| `lexical` (default, `build_profile`) | its lexeme shows alternation evidence in E0 at τ = dictionary (`p_cells_E0_dictionary`): G&W with an Inkelas learner | **circular**: a P form with the change is itself the evidence, so the match tests the rule mechanics (trigger, locality, cluster outcomes, exact bundles), not a lexical prediction |
| `all_underspecified` | always | how far one set of rules gets with no lexical specification: the rule TP (a) as a derivation |
| `lexical_loo` | the evidence comes from a cell outside the tested cell and its copies (NOM.PL = VOC.PL, feminine DAT.SG = LOC.SG) | the non-circular version of `lexical`: can specification learned elsewhere in the nominal paradigm predict this cell (spell-out TP (c) as a derivation)? |
| `spellout` | its lexeme never takes a competing N allomorph (-ové/-é, -u, -ách) in E0 (`n_cells_E0_dictionary`) | G&W's spell-out claim as a policy: non-palatalizing allomorphs select prespecified roots; not circular (reads the N choices, never the P outcome) |

N endings carry no palatalizer, so every variant predicts no change
before them, and all 22,429 N rows match in every variant.

## 5. Results

| variant | boolean match | IPA match | P rows: boolean |
|---|---|---|---|
| lexical | 48,132 / 48,143 (99.98%) | 48,118 (99.95%) | 25,703 / 25,714 |
| all_underspecified | 48,104 (99.92%) | 48,090 (99.89%) | 25,675 / 25,714 |
| lexical_loo | 25,223 (52.39%) | 25,223 (52.39%) | 2,794 / 25,714 |
| spellout | 43,632 (90.63%) | 43,619 (90.60%) | 21,203 / 25,714 |

By class × cell (boolean %; for `lexical` and `all_underspecified`, also
IPA %; N rows are 100% everywhere and left out):

| class | cell | P rows | lexical | all_under. | lexical_loo | spellout |
|---|---|---|---|---|---|---|
| velar | DAT.SG f | 3,795 | 100 / 100 | 100 / 100 | 0.0 | 100 |
| velar | LOC.SG | 3,841 | 100 / 100 | 99.9 / 99.9 | 0.7 | 99.0 |
| velar | LOC.PL | 3,539 | 100 / 100 | 100 / 100 | 36.9 | 24.6 |
| velar | NOM.PL m.an | 1,293 | 100 / 100 | 99.8 / 99.8 | 99.1 | 66.3 |
| velar | VOC.SG (big) | 13 | 53.8 / 53.8 | 38.5 / 38.5 | 53.8 | 84.6 |
| velar (adj) | NOM.PL.M.AN | 2,813 | 100 / 99.5 | 99.6 / 99.1 | 0.4 | 99.6 |
| r | DAT.SG f | 328 | 100 / 100 | 100 / 100 | 0.0 | 100 |
| r | LOC.SG | 362 | 100 / 100 | 99.7 / 99.7 | 0.3 | 91.7 |
| r | NOM.PL m.an | 372 | 100 / 100 | 99.7 / 99.7 | 17.2 | 85.8 |
| r | VOC.SG (Cr) | 69 | 94.2 / 94.2 | 91.3 / 91.3 | 92.8 | 85.5 |
| r (adj) | NOM.PL.M.AN | 44 | 100 / 100 | 90.9 / 90.9 | 9.1 | 90.9 |
| t/d/n | DAT.SG f | 2,113 | 100 / 100 | 100 / 100 | 0.0 | 100 |
| t/d/n | LOC.SG | 2,600 | 100 / 100 | 100 / 100 | 0.1 | 82.2 |
| t/d/n | NOM.PL m.an | 1,191 | 100 / 100 | 100 / 100 | 1.0 | 33.7 |
| t/d/n (adj) | NOM.PL.M.AN | 3,321 | 100 / 100 | 100 / 100 | 0.0 | 100 |

(LOC.PL P rows of r and t/d/n: 2 and 18. Full tables, also by stem final
and by cluster: `analysis/czech/czech_rule_validation_summary.txt`.)

**Mismatch classes** (`mismatch_class`; the shared categories can only
say RULE_FAILURE):

- `lexical`, 25 rows. 11 *unexplained forms*: P endings without the
  change in the data, which the data layer already flags as
  cs.wiktionary table errors without corpus tokens (búchre, lúzre,
  epistemologe, kmoche, kolegáčke, lulane). 14 *cluster outcomes*: 8
  -ský adjectives listed with -scí (badatelscí, dosovscí), 5 -cký
  adjectives with -ccí (chymiccí, epistemiccí) and vojsko → vojště (šť
  before a case ending).
- `all_underspecified`, 53 rows: the same 14 cluster outcomes and 39
  unexplained forms (now also those of lexemes with no alternation
  anywhere: metodickí, bezkorí, nenávistníki, adblockě, captchě).
- `lexical_loo`, 22,920 rows. 19,353 *no-evidence* rows of 13,121
  lexemes: the feminine DAT/LOC.SG and the adjective NOM.PL.M.AN are
  the only P cells of their lexemes, so nothing elsewhere in the
  nominal paradigm tells a learner that they alternate. 3,553
  *mixed-stem* rows: the stem alternates here but takes an N allomorph
  elsewhere, mostly velar LOC.PL -ích of stems with LOC.SG -u (2,004:
  park, parku ~ parcích) or NOM.PL -ové (189), and t/d/n NOM.PL -i
  (783) and LOC.SG -ě (460). Also 10 unexplained forms and 4 rows
  without a UR.
- `spellout`, 4,524 rows: 4,491 *mixed-stem* rows (velar LOC.PL 2,670,
  t/d/n NOM.PL 786 (-i ~ -é doublets: Afričani ~ Afričané), t/d/n LOC.SG
  462, velar NOM.PL 436 (Češi ~ Čechové), r NOM.PL 53), plus 20
  unexplained forms and 13 cluster outcomes.

## 6. What this means for the theory

1. **The trigger is morphological, and the mechanics hold.** A single
   floating palatalizer on the P endings, unified into placeless velars
   and place-open coronals, derives the small palatalization of 99.95%
   of the 48,143 cell forms exactly (100% of the N rows, and every P
   row except 11 table errors and 14 cluster variants), with no rule
   that mentions a front vowel and no feature-changing rule. The
   non-trigger front vowels (-em, -ech, -é, coronal VOC -e) come out
   right because they carry no palatalizer. This is CVW's conclusion,
   and the TP version of it in `docs/czech/README.md` (R17a fails while
   every suffix row is near-exceptionless).
2. **A single rule needs no lexical prespecification in inflection.**
   With every target underspecified, the P rows fail only on the 39
   'unexplained' table forms (e = 39 of N = 25,714, against
   θ = N/ln N ≈ 2,532). Czech inflection gives no evidence for
   inalterable (prespecified) velar or r stems. Once a P ending is
   chosen the change is automatic, loans and names included, as IJP 370
   says. The lexical variation is all in the choice of ending.
3. **Allomorph choice does not follow root specification.** G&W derive
   Steriade's co-variation from one root specification that both blocks
   the mutation and selects the non-palatalizing allomorph. As a policy,
   `spellout` gets 90.6%, worse than having no lexical specification
   (99.9%). Its 4,491 failures are mixed stems: park takes LOC.SG -u
   yet palatalizes in LOC.PL -ích (parcích), biolog takes -ové yet
   palatalizes in biolozích, and Češi stands beside Čechové. A
   prespecified root cannot select -u and still alternate before -ích,
   so in Czech the selection of -u/-ové/-ách must refer to the cell and
   the lexical class (IJP: -u "because -ě requires an alternation", a
   preference, not a root property), not to the root's specification.
   This agrees with the spell-out TP (c) of the README, where the
   prediction holds in one direction per cell.
4. **An Inkelas learner needs the cell itself.** If specification may
   only be learned from alternations elsewhere in nominal inflection,
   `lexical_loo` mispredicts 89% of the P rows: most Czech lexemes have
   a single P cell (feminine DAT/LOC.SG, adjective NOM.PL), so evidence
   from the rest of the paradigm does not exist. A learner that
   prespecifies by default and underspecifies only on evidence would
   fail on unseen cells, while default underspecification (point 2)
   would not. This matches R21 (bases without alternation evidence still
   palatalize in 95–99% of their derivatives) and the simulated learners
   in `docs/czech/README.md`, for whom the palatalization rules are
   tolerated.
5. **Where lexical specification does show, it is partial.** The one
   place where stems differ in how they palatalize, not just in whether
   they take a P ending, is the cluster -sk-/-ck-: the adjectival suffix
   gives šť čť, roots give sc. Here a G&W-style partial prespecification
   ([−strident] on the suffix's /Ḱ/) is the cleanest LP statement, and
   it holds for 2,619 of 2,632 suffix rows.

## 7. Hooks implemented in the analysis module

Nothing in the shared modules was changed. What the shared hooks could
not carry was done locally:

- **Rows.** `LexRow` has no cell, lexeme, ending or evidence fields, so
  `CzechRow` subclasses it, and every Czech hook checks for it
  (`as_czech`). The `row_factory` hook reads one `czech_cells.csv` row;
  the lexeme columns it needs (lemma IPA, E0 evidence) are joined in by
  `iter_rows`, since a row factory sees one CSV row only. For the same
  reason `phonology.validation.trace_lemmas` cannot find Czech rows (and
  its default `pos_filter="N"` would drop the adjectives), so the module
  has its own `trace()` and `--trace`.
- **Categories.** `RefinedCategory` and `Normalisation` are closed
  vocabularies with Romanian labels. `categorise` therefore returns only
  UR_BUILD_FAILED or RULE_FAILURE, `mismatch_class` names the Czech
  causes for the summary, and the normaliser cascade does all its
  folding in `strip`, with no labelled steps. Czech labels in the shared
  enums would let the validator report say what the summary says.
- **Per-cell IPA rates.** `render_report` tabulates the boolean match
  only (per stem final × opportunity and per cluster bucket), and its
  example verdicts ("rules→T, data→F") follow the predicted boolean
  even for IPA-only mismatches. The per-cell IPA rates are in the
  summary (`main()` of the module), since the driver has no hook for a
  cross-variant summary.
- **Distance.** The feature component of `diagnostics.distance`
  compares the word's last consonant. For the endings -ích, -ách, -em
  and -ech that is the ending's consonant, so feature_edit is 0 there.
  A target-locator hook would let it compare the stem edge.
- The module imports `lexpipe` (registry paths, `csvio.split_list`,
  `PROJECT_ROOT`), as the Italian analysis does.

## 8. Limits and open points

- The UR is built from the attested form, so spell-out (which ending)
  is given, and the variants test specification, not allomorph choice.
  `spellout` evaluates G&W's selection claim only through its
  consequence for palatalization.
- The attested IPA comes from the spelling. For velars and r the change
  is visible in spelling; for t/d/n it is inferred (above).
- Only nominal inflection at τ = dictionary. The big palatalization is
  modelled only in the velar vocative, and labial glide insertion
  (rybě [bjɛ]) is not modelled: labial, l, s and z stems are not in the
  rows.
- The weak palatalizer of Poĺomská (PAL0: r > ř only after a consonant,
  bratře vs doktore) is not a separate bundle. The data layer's P/N
  class of VOC -e (P after Cr) decides whether PAL1 is inserted, so the
  Cr condition sits in the allomorph choice.
- The -sk/-ck type comes from the lemma ending (-ský/-cký). The 14
  cluster mismatches are lexemes listed with the other outcome.
- `lexical_loo` writes a 2.7 MB mismatch table (22,920 rows, mostly the
  expected no-evidence rows).

## References

- Cavirani, E. & Vanden Wyngaerd, G. (2026). A representational
  analysis of Czech palatalisation. *Phonology* 43: e10.
  doi:10.1017/S0952675726100347 (CC BY 4.0;
  `data/research_v2/texts/cavirani_vw_2026.txt`).
- Gorman, K. & Winstead, J. Romanian palatalization without tears.
  Manuscript (`latex.tex`).
- IJP: Internetová jazyková příručka, ÚJČ AV ČR, articles 126, 223, 225,
  226, 228, 370.
- Inkelas, S. (1995). The consequences of optimization for
  underspecification. *NELS* 25.
- Ketner, K. (2005). Palatalisation and featural morphemes in Czech.
  *Cambridge Occasional Papers in Linguistics* 2: 129–147.
- Poĺomská, A. (2025). Phonological relations between palatalizers and
  the phonemic system: a case study on Czech. *Linguistics Beyond and
  Within* 11: 171–205.
- Šimáčková, Š., Podlipský, V. J. & Chládková, K. (2012). Czech spoken
  in Bohemia and Moravia. *JIPA* 42(2): 225–232.
- Yang, C. (2016). *The Price of Productivity*. MIT Press.
- Ylonen, T. (2022). Wiktextract: Wiktionary as machine-readable
  structured data. LREC 2022.
