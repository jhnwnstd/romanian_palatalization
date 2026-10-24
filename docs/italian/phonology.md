# Italian phonology: feature inventory and a provisional Logical Phonology analysis

This note documents the Italian feature inventory
(`src/phonology/inventories/italian_features.json`) and the provisional
analysis of Italian velar palatalization in Logical Phonology (LP)
(`src/phonology/analyses/italian_palatalization.py`). The analysis carries
Gorman & Winstead's account of Romanian (`latex.tex`) over to Italian:
alternating dorsals are underspecified, inalterable ones are prespecified,
and one feature-filling unification rule does the palatalizing. The
validator's results for three specification policies are in
`analysis/italian/`. Rule numbers (R…) refer to
[`../multilingual_design.md`](../multilingual_design.md).

```sh
python scripts/validate_lang_rules.py --lang it                   # lexical (default)
python scripts/validate_lang_rules.py --lang it --variant all_underspecified
python scripts/validate_lang_rules.py --lang it --variant stress
PYTHONPATH=src python -m phonology.analyses.italian_palatalization  # cross-variant summary
PYTHONPATH=src python -m phonology.analyses.italian_palatalization --trace amico osco
```

| output | content |
|---|---|
| `analysis/italian/italian_rule_validation.txt`, `italian_rule_mismatches.csv` | validator report and mismatches, `lexical` policy |
| `…_all_underspecified.txt/.csv`, `…_stress.txt/.csv` | the same for the other two policies |
| `analysis/italian/italian_rule_validation_summary.txt` | every policy side by side: boolean and IPA match per cell and per class, the Height-effect check, the Italian mismatch types with examples, coverage, and the coronal and lexically palatal sanity rows |

All of these hold only Wiktionary-derived forms and transcriptions (CC BY-SA)
and aggregates over them. The validation uses the dictionary threshold
(τ = dictionary) only, so no Leipzig count enters a tracked file (R9).

**In brief.** One rule derives every attested plural of every root with
alternation evidence. Its trigger is the coronal vocoid class /i j/, which
gives the Height effect: no underspecified root ever palatalizes before
the feminine -e (2,553 rows). The rule also produces the geminate, cluster
and sc outputs without extra statements. The only failures under the
lexical policy are doublets and three Wiktionary editing errors. That
policy is circular for the -i cell, so it tests mechanics, not
prediction. A single unconditioned rule misses 1,080 of 7,338 rows,
mainly post-tonic velars and NDE stems. A stress-based grammar with no
listing misses 174, mainly the amico type (half of it vowel + -ico
adjectives), profugo-type g-stems, the carico type and -logo things.

## 1. The inventory (`italian_features.json`)

The file follows `src/phonology/inventories/README.md`. It has the 27
features of `romanian_features.json` in the same order and 51 segments,
whose feature bundles are all distinct, both in the file and after the
analysis patches. `tests/test_phonology_italian.py` checks every claim in
this section.

**Shared with Romanian, copied exactly.** i u a, p b t d k ɡ, t͡s t͡ʃ d͡ʒ,
f v s z ʃ ʒ, m n l r, j w. These include the Romanian conventions:

- velars are DORSAL + and 0 for the tongue-body features (High, Low,
  Front, Back, Round), which are reserved for vocoids;
- DelRel is - on oral stops, + on affricates and 0 elsewhere;
- t and d have Strident 0 in the base. The Romanian analysis patches them
  to -; the Italian analysis does not (§3).

**Deviations from Romanian.** The only one is e and o, which have Tense +
(close-mid) where Romanian has 0. Romanian has no tense/lax contrast,
but Italian has seven stressed vowels, with e/ɛ and o/ɔ contrasting
(pesca 'fishing' /e/ vs pesca 'peach' /ɛ/). The contrast holds only under
stress. That is a distributional fact; the inventory does not encode it.
Tense stays 0 on i u a, as in Romanian, so Tense is specified on the mid
vowels only.

**New segments.**

| segment | built from | justification |
|---|---|---|
| ɛ, ɔ | e, o with Tense - | open-mid vowels. Wiktionary writes them (ˈfwɔ.ko, ˈmɛ.di.ko: 6,047 and 4,681 occurrences in the transcriptions); Romanian folds them to e o, Italian cannot |
| d͡z | t͡s with Voice + | lexical voiced affricate (zero /ˈd͡zɛ.ro/, zona, mezzo /ˈmɛd.d͡zo/). It is never derived; Romanian has none, and the Romanian analysis derives that gap |
| ɲ | n with Anterior -, Distributed + | gn (ragno /ˈraɲ.ɲo/). Option A of the research notes: coded like the postalveolars ʃ t͡ʃ d͡ʒ of the base (coronal, DORSAL -). Option B, the Hayes dorso-palatal with High/Front +, would break the convention that tongue-body features belong to vocoids |
| ʎ | l with Anterior -, Distributed + | gli (figlio /ˈfiʎ.ʎo/); option A as for ɲ |
| 20 geminates | the singleton with Long + | pː bː tː dː kː ɡː t͡sː d͡zː t͡ʃː d͡ʒː fː vː sː ʃː mː nː ɲː lː ʎː rː. Consonant length is contrastive (fato/fatto). z j w do not geminate, and marginal ʒ has no geminate |

**Absent.** Romanian ɨ ə h are absent; h is a silent letter in Italian (a
single transcription has an optional "(h)"). [ŋ ɱ] are allophones of n
and Wiktionary writes /n/. ʒ is kept as a marginal segment: 18
transcriptions have it outside d͡ʒ, in loans (garage /ɡaˈraʒ/, abat-jour)
and in a few Tuscan-style forms (frogia /ˈfrɔ.ʒa/).

**Stress and length.** Stress is - on every segment. Stress is not
segmental here: the tokenizer drops ˈ and ˌ, and the stress policy reads
stress from the transcription through the data layer's `velar_post_tonic`.
Length is + only on geminates. Phonetic lengthening of stressed open
vowels is allophonic and not encoded.

**Labels are tokens.** Every label is exactly the token
`ITALIAN_TOKENIZER` produces: ties U+0361, ː attached. Wiktionary writes a
geminate across the syllable dot (ˈsak.ko, ˈfat.t͡ʃa, viʎˈʎak.ko: 3,487 t.t,
869 k.k, 393 ʎ.ʎ…). `merge_geminates` rewrites it as one long segment
(ˈsakːo, ˈfat͡ʃːa, viˈʎːakːo) before tokenizing. The stop before an
affricate is read as the affricate's first half (t.t͡ʃ → t͡ʃː, d.d͡z →
d͡zː).

## 2. Tokenizer and comparison convention

`ITALIAN_TOKENIZER` has these settings:

- it skips stress marks, syllable dots and whitespace;
- it attaches ː to the preceding segment (`LengthHandling.ATTACH`);
- it keeps combining marks on their base;
- it folds only ASCII g → ɡ, and leaves ɛ and ɔ distinct (Romanian folds
  them to e and o).

`ITALIAN_NORMALISERS` only strips: it merges geminates, drops stress marks,
dots and spaces, and keeps ː. It has no smoothing steps, because the
attested plural is rebuilt from the same lemma transcription as the
underlying representation (UR) (§6).

## 3. Patches (`PATCHES`)

| segment | patch | why |
|---|---|---|
| i, j | CORONAL +, Anterior -, Distributed + | as in Romanian: a high front vocoid is coronal (postalveolar place). In Italian this patch *is* the Height effect: the dorsal rule's trigger is the coronal vocoid class, which holds /i j/ and not /e/ (CORONAL -) |
| t͡ʃ, d͡ʒ | DelRel 0, DORSAL + | the Romanian patch. Palatalization leaves Dorsal untouched, so the derived affricate is dorso-coronal; DelRel is left free so the segment derived from /k/ (DelRel -) projects onto t͡ʃ |
| t͡ʃː, d͡ʒː | DelRel 0, DORSAL + | the same for geminates (tombacco/tombacci) |
| ʃː | DORSAL + | the palatalized velar of an sc unit coalesces into [ʃː] and keeps its +Dorsal (§5). Singleton ʃ is never derived and is not patched |

Romanian's other patches have no Italian counterpart:

- t and d Strident - protect Romanian /t d/ from assibilation, and Italian
  has no assibilation;
- t͡s DelRel 0 lets Romanian derive [t͡s], and Italian derives no coronal
  affricate.

## 4. Underspecified segments: dorsals only (`UNDERSPEC`)

| label | base | cleared | used for |
|---|---|---|---|
| K, G | k, ɡ | CORONAL | as in Romanian: amico, medico, asparago, belga |
| Kː, Gː | kː, ɡː | CORONAL | geminates: tombacco/tombacci [t͡ʃː] |
| Kˢ | k | CORONAL, Continuant, Long | the velar of an alternating sc unit: osco/osci, falisco/falisci [ʃː] |

**No underspecified coronal.** No rule writes coronal features into t d s
z, so Italian's lack of a coronal rule (gatto/gatti, naso/nasi) is a fact
about the inventory: the inventory of underspecified segments has no
/T D S Z/. The coronal disjuncts of G&W's spell-out condition P are then
vacuous.

The sc unit therefore puts its extra openness on the dorsal, not on the
s. The /s/ stays prespecified, and that is what keeps the analysis
NDEB-safe. Lexical [s.t͡ʃ] does occur: giuscibernetica
/ˌd͡ʒus.t͡ʃi.berˈnɛ.ti.ka/, at a compound boundary. A feature-changing rule
on prespecified /s/ would wrongly reach it; the rules below cannot.

## 5. Rules (`RULES`, in order)

LP notation: [ ] is a natural class, { } a set of valued features, ⊔
unification, \ subtraction. Every unification rule is feature-filling, so
a prespecified segment resists by unification failure, as in Romanian.

1. **dorsal-pal** (latex.tex e:velar-pal):
   [−Syll, +Cons, −Son, −Approx, −Nas] ⊔ {+Cor, −Ant, +Dist, +Str} /
   \_ [−Cons, +Cor, −Ant, +Dist], with a broad terminator (only the next
   segment counts). Dorsal is untouched, so the output is the dorso-coronal
   t͡ʃ or d͡ʒ.
   - The trigger is the coronal vocoid class, which the patches restrict
     to /i j/. [+High, +Front] picks out the same two segments. The place
     statement follows G&W's revision ("/i/ was arguably a place trigger
     from the start").
   - Romanian's dorsal rule instead takes [+Syll, +Front] and so fires
     before -e. The two languages differ in this one condition.
2. **sc-coalescence**: [+Cor, −Ant, +Dist, +Str] ⊔ {+Cont, +Long} / [s] \_
   (left broad terminator). Only a palatalized /Kˢ/ is open for Continuant
   and Long. Lexical t͡ʃ, ʃ and ʃː resist, which is why giuscibernetica is
   safe.
3. **s-deletion**: [s] → ∅ / \_ [ʃː]. Italian has no [sʃ] sequence (there
   are none in the 40,481 transcriptions), so the rule never meets a
   lexical one. Together with rule 2 this gives Italian /sk/ + i → [ʃː]:
   the two consonants of the unit fuse into one long fricative. Romanian
   instead keeps both and bleeds the affricate (muscă → [muʃte], e:bleed-rev).
4. **dorsal-default**: [obstruent] ⊔ {−Cor, +Dor, −Cont}. This is e:default
   literally (the Romanian code supplies only −Cor, because /K/ keeps the
   other two values). Unpalatalized /K G Kː Gː Kˢ/ surface as velar stops.
5. **length-default**: [+Cons] ⊔ {−Long}. It reaches only an unpalatalized
   /Kˢ/.
6. **glide-absorption**: /j/ → ∅ / \_ /i/. The ending -i absorbs a root
   glide: occhio /ɔkːj+i/ → [ɔkːi], studio → studi. Italian has no [ji].
   This is not a palatalization rule.

`PALATALIZATION_RULE_NAMES` = {dorsal-pal, sc-coalescence, s-deletion}.
The validator's boolean prediction is "one of these fired".

## 6. URs, rows and the attested plural

As in the Romanian validator, spell-out is not modelled. The UR is built
as follows:

1. take the lemma's transcription and merge its geminates;
2. drop the final a/e/o;
3. swap the stem-final dorsal for its underspecified twin when the
   variant's policy says so (§7);
4. add the cell's ending: -i for the i cell, -e for the e cell.

For an sc unit the policy picks /Kˢ/. In the occhio class the target is
the velar before the root vocoid.

**Rows** (`iter_rows`, the driver contract) come from
`data/italian/italian_cells.csv`, which has one row per lexeme × plural
cell × form. The analysis keeps rows that are:

- in-domain, aligned, unfiltered and not curator-rejected;
- before -i or -e;
- dorsal targets or NDE classes (occhio, psiche);
- for lemmas with a transcription in `italian_lexicon_complete.csv`
  (primary IPA, normalized by the Italian `lexpipe` normalizer).

A row without a transcription cannot be derived. Scoring it with the
fallback would only restate the policy, so such rows are counted and left
out:

| class | in scope | validated |
|---|---|---|
| dorsal | 13,886 | 7,128 |
| nde:occhio | 308 | 207 |
| nde:psiche | 4 | 3 |
| total | 14,198 | 7,338 (4,504 lexemes, 4,066 roots) |

**The attested plural** (`plural_ipa`) is a reconstruction, because the
generated data have no plural transcriptions (design §3: "plural IPA
derived from the lemma"). It is built from the lemma transcription:

1. the final vowel is replaced by the ending;
2. the stem-final velar is replaced by the palatal the spelling shows
   (ci gi → t͡ʃ d͡ʒ, cci → t͡ʃː, sci → ʃː) or kept (chi ghi);
3. an occhio-class root glide is dropped before -i.

Stress marks and dots are kept for readability. The IPA match therefore
tests the derived segments at the stem edge against the orthographic
outcome; the rest of each word is the lemma's own transcription. When the
transcription does not fit the spelling, the reconstruction refuses (the
row has no UR, and the fallback scores it). No validated row does on the
current data: ceviche /seˈvi.t͡ʃe/, a loan whose -che is [t͡ʃe], once did,
but the data layer no longer puts it in the psiche class
(`features.velar_che`).

**The oracle** (`expected_firings`) is the attested outcome: dorsal-pal
for a palatalized form, plus the two sc rules for an sc unit. Italian
palatalization is lexical, so no productivity claim can stand in. The
distance metric's rule component therefore measures disagreement with the
data. `fallback_predict` applies the variant's policy and never the data.

## 7. Three specification policies (`PROFILE_VARIANTS`)

The rules are the same in every variant; only the policy differs.

- **`lexical`** (`build_profile()`, Gorman & Winstead). A target is
  underspecified iff its **root** shows alternation evidence in E0 (the
  nominal plural cells) at τ = dictionary.
  - The evidence is computed with the TP stage's own code
    (`lexpipe.langs.italian.roots`).
  - Learner evidence pools all forms of a root (R4), so amica shares
    amico's /K/.
  - **This is circular for the -i cell.** A root is underspecified
    because its -i plural palatalizes, so the -i boolean match tests the
    rule mechanics (trigger, locality, projection, cluster and geminate
    outputs), not a lexical prediction.
  - The lexical policy can still fail in three places, none of which it
    looks at: the -e cell of an underspecified root, doublets (one UR has
    one output), and the surface IPA.
  - 5,477 rows (2,641 roots) are underspecified.
- **`all_underspecified`**. Every dorsal target is underspecified, NDE
  stems included (7,338 rows). This measures how far a single rule /K/ →
  [t͡ʃ] / \_ i gets: the R17b "palatalize" row turned into derivations.
- **`stress`**. A target is underspecified iff the stem-final velar is
  **not** immediately post-tonic (5,619 rows). Six rows have unknown
  stress and count as not post-tonic. This is a learner-style grammar
  with no listing: the stress-nested TP rows (retain post-tonic,
  palatalize elsewhere) turned into derivations.

## 8. Derivations

The golden derivations (lexical policy, `tests/test_phonology_italian.py`)
run on the data's transcriptions. A blank cell means the rule did not
apply.

| UR | amiKi | amiKe | fwɔki | mɛdiKi | laɡi | asparaGi | bɛlGi | ɔkːji | ɔsKˢi | tombaKːi |
|---|---|---|---|---|---|---|---|---|---|---|
| dorsal-pal | amit͡ʃi | | | mɛdit͡ʃi | | asparad͡ʒi | bɛld͡ʒi | | ɔst͡ʃi | tombat͡ʃːi |
| sc-coalescence | | | | | | | | | ɔsʃːi | |
| s-deletion | | | | | | | | | ɔʃːi | |
| dorsal-default | | amike | | | | | | | | |
| glide-absorption | | | | | | | | ɔkːi | | |
| **SR** | amit͡ʃi | amike | fwɔki | mɛdit͡ʃi | laɡi | asparad͡ʒi | bɛld͡ʒi | ɔkːi | ɔʃːi | tombat͡ʃːi |

Each column illustrates one point:

- amica + e: the root is /amiK/, yet dorsal-pal does not apply because e
  is not a coronal vocoid. This is the Height effect with an
  underspecified root.
- fuoco and lago: the roots are prespecified, so nothing applies.
- occhio: no palatalization rule applies. The prespecified velar is
  inalterable, and the root glide is absorbed.
- osco: the sc unit coalesces, giving osci [ˈɔʃːi].

Under the other policies:

- `all_underspecified` derives fuoco → *[fwɔt͡ʃi], occhio → *[ɔt͡ʃːi]
  (before the root glide) and bosco → *[bɔʃːi].
- `stress` derives amico → *[amiki] and belga → *[bɛlɡi], while medico →
  [mɛdit͡ʃi].

## 9. Results (real data, 7,338 rows)

From `analysis/italian/italian_rule_validation_summary.txt`. "bool" is
the boolean mutation match, "IPA" the surface match. The two rates are
equal in every variant: wherever the policy put the right specification,
the derivation reproduced the attested segment exactly. There is no
IPA-only mismatch.

| policy | bool = IPA | boolean mismatches | over / under | c+i | g+i | c+e | g+e |
|---|---|---|---|---|---|---|---|
| lexical | 99.52% | 35 rows, 30 roots | 32 / 3 | 99.46% | 96.29% | 99.91% | 100% |
| stress | 97.63% | 174 rows, 153 roots | 96 / 78 | 96.67% | 84.86% | 99.91% | 100% |
| all_underspecified | 85.28% | 1,080 rows, 969 roots | 1,077 / 3 | 78.01% | 36.57% | 97.67% | 98.39% |

**Height effect.** Underspecified dorsal targets before -e number 2,553
rows (2,501 roots) under `lexical`, 2,706 rows (2,654 roots) under
`stress` and 3,445 rows (3,375 roots) under `all_underspecified`. The
rules palatalize none of them. The only -e "palatalizations" under
`all_underspecified` are 76 NDE rows (anarchia → *anart͡ʃie, macchia →
*mat͡ʃːje), triggered by the root vocoid, not by -e. The three -e mismatches every policy shares (elettrochimica →
elettrochimice, fotomeccanica, idromeccanica) are the data layer's QC
suspects (`e_palatal`): probable Wiktionary editing errors.

**Per class (-i cell, boolean = IPA).**

| class | rows | lexical | all_underspecified | stress |
|---|---|---|---|---|
| R26 suffix -ico | 1,576 | 100% | 100% | 98.60% |
| suffix -iaco | 13 | 100% | 100% | 100% |
| root, velar not post-tonic | 1,202 | 99.17% | 94.18% | 94.18% |
| root, velar post-tonic | 811 | 98.89% | 5.92% | 94.08% |
| -logo person | 71 | 92.96% | 92.96% | 92.96% |
| -logo thing | 9 | 100% | 0% | 0% |
| -urgo person | 14 | 71.43% | 35.71% | 64.29% |
| nde:occhio | 131 | 100% | 0% | 100% |
| nde:psiche | 3 | 100% | 0% | 66.67% |

**Mismatch types** (`mismatch_class`, named from the data side):

- **`lexical`**:
  - doublets: 32 rows, 27 roots, among them chirurgo, manico, stomaco,
    sarcofago, astrologo, bolscevico, osco and falisco. A root with both
    a palatal and a velar -i plural gets one UR, so one variant fails.
  - the 3 editing errors.
- **`all_underspecified`**:
  - post-tonic velar kept (fuoco type): 759 rows, 688 roots;
  - NDE stems palatalized: 207 occhio rows, 180 roots, plus 3 psiche;
  - doublets: 32;
  - g-stems not post-tonic kept (profugo, naufrago, prodigo, obbligo,
    arcipelago…): 23;
  - other proparoxytone velars kept (abaco, anfibraco, basilico,
    indaco…): 21;
  - deverbal from -are (carico, incarico, pizzico, strascico,
    solletico…): 16;
  - -logo things and adjectives (catalogo, dialogo, prologo; analogo):
    12;
  - 3 editing errors; four rows of unknown class or stress
    (neuronofago, riepilogo, pogginco twice).
- **`stress`**:
  - post-tonic palatalizers (amico type): 62 rows, 54 roots;
  - doublets: 32;
  - the same proparoxytone velars as above: 23 g-stems, 21 others, 16
    deverbal, 12 -logo things and adjectives;
  - sineddoche (psiche class, velar not post-tonic);
  - 3 editing errors and the same four rows of unknown class or stress.

  Of the 53 amico-type lemmas, 29 are vowel + -ico adjectives whose
  transcription writes a diphthong: arcaico /arˈkaj.ko/, laico, eroico,
  ebraico, voltaico, paranoico, mosaico… The others are amico, nemico,
  greco and compounds, porco, belga, daco, meteco, taumaturgo, tombacco.

**Sanity rows.** 12,783 coronal rows and 1,507 lexically palatal rows
(luce, valigia, bacio) pass under every policy: no rule touches a
prespecified coronal or a lexical palatal. The coronal part is true by
construction, since the data layer codes every aligned coronal plural as
unchanged. What it does check is that the inventory, tokenizer and UR
builder handle 14,290 more transcriptions without an unknown segment.

## 10. What this means for the theory

1. **The rule mechanics hold.** One feature-filling rule derives every
   attested Italian plural of every root it is given as underspecified.
   That covers single velars, geminates, C + velar clusters (belgi), the
   sc coalescence ([ʃː], not Romanian [ʃt]) and the occhio glide, with
   identical IPA throughout. NDEB follows from prespecification, as in
   Romanian:
   - the occhio class (velar + root vocoid: 207 rows) is never
     palatalized under `lexical`;
   - root-internal velars before front vowels (chilo, ghiro) are never
     targets.
2. **The Height effect is a trigger fact, not a specification fact.**
   Masculine -i triggers and feminine -e does not, even on the same
   underspecified root: amici/amiche, medici/mediche, pubblici/pubbliche,
   and 2,553 rows in all. Romanian's CORONAL+ patch on /i j/ derives this
   with no reference to the cell or the gender. The cross-language
   difference is the dorsal rule's condition: Romanian [+front vowel],
   Italian [coronal vocoid].
3. **The lexical policy's success is not evidence for underspecification
   over inflection dependence.** It is circular for the -i cell: the
   policy reads the outcome. Its non-circular content is points 1 and 2
   and its one real failure, **doublets**. Twenty-seven roots list both
   outcomes (chirurghi ~ chirurgi, manici ~ manichi, osci ~ oschi), and
   a single root specification cannot derive both. G&W need two listed
   URs per doublet (one prespecified, one underspecified) or a variable
   choice between them. The data cannot tell these apart; R19's coding
   (dominant variant, exception, excluded) decides how much it matters.
4. **No single rule.** A rule that palatalizes every target before -i
   (`all_underspecified`) fails on 1,080 rows, 969 roots. Seventy percent
   of these are post-tonic velars (fuoco, lago, bosco, sacco), and 19%
   are NDE stems, which only prespecification protects. This is the
   derivational counterpart of the R17b row: neither palatalize nor
   retain is productive over -i.
5. **Stress gets most of the way without listing.** Underspecifying
   exactly the velars that are not post-tonic (`stress`) leaves 174 rows
   and 153 roots (27 of them doublets). These are the exceptions that
   the it_stress TP rows find tolerable:
   - the amico type: 54 roots, more than half of them vowel + -ico
     adjectives, so the laico syllabification question decides a large
     part of it;
   - Flor's voicing effect in g-stems (profugo, naufrago);
   - the deverbal carico type;
   - -logo things.

   For a learner in Inkelas's sense, who underspecifies only after seeing
   alternation, `stress` describes a default specification for roots
   without evidence. That default is not what the lexical account
   assumes: a pure lexical learner prespecifies a novel stem. The lexicon
   cannot decide between the two defaults. Nonce studies can (Krämer's
   palatalize-all vs palatalize-none groups; Giavazzi's stress-conditioned
   productions).
6. **No allomorph selection.** Italian has one plural exponent, -i. A
   prespecified root is marked only by keeping its velar, never by a
   different allomorph (Romanian -uri). G&W's spell-out condition P
   therefore has no plural correlate in Italian. Its testable remnant is
   that learned /i/-initial derivatives read the same /K/ vs /k/ contrast
   (the co-variation tables in `italian_covariation.txt`).

## 11. Hooks implemented locally (requests for the shared core)

The shared modules were not edited. These gaps were filled inside the
analysis module and are proposed for the core:

1. **Per-row evidence.** `LexRow` has no place for a root, its alternation
   evidence, stress or class, and the `ur_builder` hook receives only the
   row. `ItalianRow` subclasses `LexRow` with these fields. A plain
   `LexRow` gets no UR. *Request:* an `extras` mapping on `LexRow`, or
   hooks generic over the row type.
2. **Category vocabulary.** `RefinedCategory` has no Italian labels, and
   adding them means editing the shared enum. The Italian categoriser
   returns only UR_BUILD_FAILED or RULE_FAILURE; `mismatch_class` and the
   summary file name the Italian causes. *Request:* a free label on
   `Categorisation`, or a profile-supplied vocabulary.
3. **Geminates across the syllable dot.** A `TokenizerConfig` cannot merge
   two tokens. `merge_geminates` therefore runs before tokenizing, in the
   UR builder, in the reconstruction and in the comparison strip. Because
   the distance metric tokenizes the attested IPA with
   `profile.tokenizer` alone, the attested plural is stored already
   merged. *Request:* an optional pre-tokenization normaliser on
   `TokenizerConfig`.
4. **IPA rates per cell and per class.** `render_report` prints only
   boolean rates per cell and per cluster. The summary file
   (`render_summary`, `python -m phonology.analyses.italian_palatalization`)
   prints both rates per cell and per class for every policy.
5. **Edit distance counts marks.** `compare_ipa` computes its edit distance
   on the unstripped strings. The Italian attested plural keeps the
   lemma's ˈ and syllable dots, so the IPA component of the distance is
   inflated by a constant-ish amount; rankings are unaffected.
6. **`trace_lemmas`** reads a CSV through `row_factory` with the Romanian
   filter defaults, but Italian rows need joined data (lemma IPA, root
   evidence). `trace()` and `--trace` replace it.
7. **Data layer.** The data have no plural transcriptions, but kaikki's
   form-of entries carry some (amici /aˈmi.t͡ʃi/). Projecting them in the
   `extract` stage would give an independent attested plural. Until
   then, the IPA match tests only the derived stem edge.
8. **Driver.** `scripts/validate_lang_rules.py` writes the default run
   without a suffix. `--variant lexical` would duplicate it, so it is not
   run.

## 12. Open points

- **The 6,860 rows without a lemma transcription** are not validated
  (stress is missing for 37% of in-domain lexemes; data layer notes). The
  it.wiktionary accented hyphenation gives stress for many of them but no
  segments. A spelling-based G2P restricted to the stem edge could bring
  them in.
- **Vowel + -ico adjectives.** These are the bulk of the `stress`
  policy's amico type. Is Wiktionary's diphthong (/ˈlaj.ko/) the right
  syllabification, or are they proparoxytones with an -ico suffix?
- **Doublets.** Two URs or a variable rule? R19's three codings bound
  the effect on the TP.
- **Verbs** (E1) are not validated. On an underspecified root the rule
  would apply before 2sg -i and 1pl -iamo (its trigger includes /j/) in
  every conjugation, but the first conjugation never palatalizes
  (cerchi, cerchiamo; traffichi beside the noun traffici). That needs
  category-sensitive specification (an underlying theme vowel, or a
  verb-class condition) and is the obvious next test of this rule set.
