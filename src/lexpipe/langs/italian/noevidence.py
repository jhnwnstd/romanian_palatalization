"""Stage ``noevidence``: roots without alternation evidence (R21).

The plural-less prediction: a learner who underspecifies a velar only
after seeing it alternate has no reason to underspecify a root that
never alternates in its inflection, so its derivatives should resist
palatalization. Italian candidates are the feminine-only -ca/-ga stems
(their only plural, -che/-ghe, never alternates) and the uncountable
-co/-go stems (critique: lega → leghista, velar kept; but musica →
musicista, palatalized).

Base classes (R21: from the corpus, not from dictionary paradigms)
------------------------------------------------------------------
At each τ, from the root's E0 realizations before -i (masculine
plurals of its lexemes; :mod:`.roots`):

    alternating      some -i form attested at τ is palatalized
                     (doublets included: an alternating form exists)
    non-alternating  -i forms attested at τ, all with the velar kept
    no evidence      no -i form attested at τ: feminine-only roots
                     (only -e plurals), roots without a plural, and
                     roots whose -i plural is a dictionary form below τ

A root without a noun/adjective -i form whose lemma has a determiner or
pronoun entry with a masculine plural attested at τ (poco N is
invariable, but pochi DET/PRON is the plural of the same word:
``base_det_plural`` of :mod:`.derivations`) takes its class from that
plural instead (pochissimo then has a non-alternating base).

Derivatives: the accepted links of ``italian_derivations.csv``
(candidates generated in both spellings and looked up in the headwords
and the corpus, R21), attested at τ, with a palatal or velar edge. The
unit is (base root, suffix, derivative stratum), coded P (only
palatalized), N (only velar) or both (R23).

Tables: 3 × 2 (base class × derivative P/N, ``both`` counted apart)
per suffix of the learned and the native sets and per derivative
stratum (R22), for the named links and with the generated candidates;
τ = dictionary and ≥ 1 token in full, the other thresholds in short.
The no-evidence class is broken down by the base descriptors of
:mod:`.derivations` (``ic_suffix``: the base has the -ic- suffix;
``masc_counterpart``: a masculine lexeme on the same stem has an -i
plural, musica/musico; feminine-only; no -i cell; -i cell below τ),
and the derivatives whose named base is a feminine lexeme are listed
separately (musicista, politichese). Example words, no per-word counts
(R9).

Output: ``analysis/italian/italian_noevidence.txt`` (tracked).
"""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from lexpipe import attest
from lexpipe.context import RunContext, run_stage
from lexpipe.langs.italian import derivations, etymology, roots, spellout
from lexpipe.langs.italian import suffixes as sx
from lexpipe.langs.italian.layout import Layout

STAGE: Final[str] = "noevidence"
# 2: gated report writer; shared stratum classifier upstream.
# 3: determiner/pronoun plurals of the base lemma are evidence (pochi).
# 4: the candidate tier with the palatal-base rejections.
VERSION: Final[int] = 4
SOURCE: Final[str] = "wiktionary_en"
ALTERNATING: Final[str] = "alternating"
NON_ALTERNATING: Final[str] = "non-alternating"
NO_EVIDENCE: Final[str] = "no evidence"
CLASSES: Final[tuple[str, ...]] = (ALTERNATING, NON_ALTERNATING, NO_EVIDENCE)
FULL_TAUS: Final[tuple[str, ...]] = ("dictionary", "ge1")
TIERS: Final[Mapping[str, frozenset[str]]] = {
    "named": frozenset({derivations.TIER_NAMED}),
    "named+candidates": frozenset(
        {derivations.TIER_NAMED, derivations.TIER_CANDIDATE}
    ),
    # + the candidates only the palatal-base rule rejected (R16: that
    # rule removes palatal candidates only; spellout.OTHER_BASE_TIER)
    "named+candidates+other_base": frozenset(
        {
            derivations.TIER_NAMED,
            derivations.TIER_CANDIDATE,
            spellout.OTHER_BASE_TIER,
        }
    ),
}
EXAMPLES: Final[int] = 15


def base_class(record: roots.RootRecord, tau: str) -> str:
    """The R21 class of a root at τ (module docstring)."""
    attested = [
        r
        for r in record.realizations
        if r.evidence == "E0" and r.vowel == "i" and r.attested.get(tau)
    ]
    if any(r.palatal is True for r in attested):
        return ALTERNATING
    if attested:
        return NON_ALTERNATING
    return NO_EVIDENCE


def det_plural_class(found: Sequence[spellout.Derivative], tau: str) -> str:
    """The class a determiner/pronoun plural of the base gives at τ
    (module docstring), or ""."""
    for item in found:
        if (
            item.base_det_plural_palatal is not None
            and tau in item.base_det_plural_taus
        ):
            return (
                ALTERNATING
                if item.base_det_plural_palatal
                else NON_ALTERNATING
            )
    return ""


@dataclass(frozen=True, slots=True)
class Unit:
    """One base × suffix × derivative-stratum unit with its base class."""

    root_id: str
    lemma: str
    suffix: str
    group: str
    stratum: str
    klass: str
    outcome: str
    words: tuple[str, ...]
    subtype: str
    feminine_named: bool
    masc_counterpart: str
    ic_suffix: bool


def subtype(
    record: roots.RootRecord, found: Sequence[spellout.Derivative], tau: str
) -> str:
    """Why a no-evidence root has no evidence."""
    if not any(
        r.vowel == "i" for r in record.realizations if r.evidence == "E0"
    ):
        if not record.has_noun or all(
            r.vowel != "e" for r in record.realizations if r.evidence == "E0"
        ):
            return "no plural"
        return "feminine-only (-e plurals)"
    return f"-i plural below {tau}"


def build_units(
    records: roots.Records,
    items: Sequence[spellout.Derivative],
    tau: str,
    tiers: frozenset[str],
) -> list[Unit]:
    """Units of every base class at τ (module docstring), sorted."""
    grouped: dict[tuple[str, str, str], list[spellout.Derivative]] = (
        defaultdict(list)
    )
    for item in items:
        if (
            item.tier in tiers
            and item.palatalized is not None
            and item.attested[tau]
            and item.group in (sx.LEARNED, sx.NATIVE)
        ):
            grouped[
                (item.base_root_id, item.suffix, item.deriv_stratum)
            ].append(item)
    units: list[Unit] = []
    for (root_id, suffix, stratum), found in sorted(grouped.items()):
        record = records.roots[root_id]
        klass = base_class(record, tau)
        if klass == NO_EVIDENCE:
            klass = det_plural_class(found, tau) or klass
        kinds = {item.palatalized for item in found}
        outcome = (
            "both" if len(kinds) == 2 else ("P" if True in kinds else "N")
        )
        units.append(
            Unit(
                root_id=root_id,
                lemma=record.lemma,
                suffix=suffix,
                group=found[0].group,
                stratum=stratum,
                klass=klass,
                outcome=outcome,
                words=tuple(sorted({item.word for item in found})),
                subtype=(
                    subtype(record, found, tau) if klass == NO_EVIDENCE else ""
                ),
                feminine_named=any(
                    item.base_pos == "noun" and item.base_gender == "f"
                    for item in found
                ),
                masc_counterpart=next(
                    (i.masc_counterpart for i in found if i.masc_counterpart),
                    "",
                ),
                ic_suffix=any(item.ic_suffix for item in found),
            )
        )
    return units


def _cells(units: Iterable[Unit]) -> Counter[tuple[str, str]]:
    return Counter((u.klass, u.outcome) for u in units)


def table_lines(units: Sequence[Unit], label: str) -> list[str]:
    """One 3 × 2 table (P/N columns, both apart)."""
    counts = _cells(units)
    lines = [
        f"    {label}",
        f"      {'base':<16} {'P':>5} {'N':>5} {'both':>5}",
    ]
    for klass in CLASSES:
        lines.append(
            f"      {klass:<16} {counts[(klass, 'P')]:>5} "
            f"{counts[(klass, 'N')]:>5} {counts[(klass, 'both')]:>5}"
        )
    return lines


def _examples(
    units: Sequence[Unit], test: Callable[[Unit], bool], label: str
) -> str:
    chosen = sorted(
        f"{u.lemma} → {'/'.join(u.words)} ({u.outcome})"
        for u in units
        if test(u)
    )
    shown = "; ".join(chosen[:EXAMPLES])
    more = f"; … (+{len(chosen) - EXAMPLES})" if len(chosen) > EXAMPLES else ""
    return f"      {label} [{len(chosen)}]: {shown}{more}"


def breakdown_lines(units: Sequence[Unit]) -> list[str]:
    """The no-evidence class by base descriptor (module docstring)."""
    none = [u for u in units if u.klass == NO_EVIDENCE]
    lines = ["    no-evidence units by base descriptor (P/N/both):"]
    groups: list[tuple[str, Callable[[Unit], bool]]] = [
        ("-ic- suffix base", lambda u: u.ic_suffix),
        (
            "masculine counterpart (other root)",
            lambda u: u.masc_counterpart == derivations.OTHER_ROOT,
        ),
        (
            "masculine counterpart (same root)",
            lambda u: u.masc_counterpart == derivations.SAME_ROOT,
        ),
        ("neither", lambda u: not u.ic_suffix and not u.masc_counterpart),
    ]
    for label, test in groups:
        chosen = [u for u in none if test(u)]
        counts = Counter(u.outcome for u in chosen)
        lines.append(
            f"      {label:<36} {counts['P']:>4}/{counts['N']:>4}/{counts['both']:>3}"
        )
    for name in sorted({u.subtype for u in none}):
        chosen = [u for u in none if u.subtype == name]
        counts = Counter(u.outcome for u in chosen)
        lines.append(
            f"      {name:<36} {counts['P']:>4}/{counts['N']:>4}/{counts['both']:>3}"
        )
    lines.append("    no-evidence examples (words only):")
    lines.append(_examples(none, lambda u: u.outcome != "N", "palatalized"))
    lines.append(_examples(none, lambda u: u.outcome == "N", "velar"))
    return lines


def feminine_lines(units: Sequence[Unit]) -> list[str]:
    """Derivatives whose named base is a feminine lexeme."""
    chosen = [u for u in units if u.feminine_named]
    lines = [
        "    derivatives of a named feminine base (musicista, politichese):"
    ]
    for klass in CLASSES:
        part = [u for u in chosen if u.klass == klass]
        counts = Counter(u.outcome for u in part)
        lines.append(
            f"      root {klass:<16} P {counts['P']}, N {counts['N']}, "
            f"both {counts['both']}"
        )
    lines.append(_examples(chosen, lambda u: True, "units"))
    return lines


def section(units: Sequence[Unit], title: str, full: bool) -> list[str]:
    """Tables of one τ and tier."""
    lines = [title, "  " + "-" * 70]
    for group, names in (
        ("learned set", sx.LEARNED_SET),
        ("native control", sx.NATIVE_SET),
    ):
        chosen = [u for u in units if u.suffix in names]
        lines.append(f"  {group}:")
        lines += table_lines(chosen, "all suffixes, all strata")
        lines += table_lines(
            [u for u in chosen if u.stratum == etymology.INTERNAL],
            "all suffixes, internal formations",
        )
        if not full:
            continue
        for stratum in etymology.STRATA:
            if stratum == etymology.INTERNAL:
                continue
            part = [u for u in chosen if u.stratum == stratum]
            if part:
                lines += table_lines(part, f"all suffixes, stratum {stratum}")
        for name in names:
            part = [u for u in chosen if u.suffix == name]
            if part:
                lines += table_lines(part, f"{name}, all strata")
        lines += breakdown_lines(chosen)
        lines += feminine_lines(chosen)
    lines.append("")
    return lines


def report(
    records: roots.Records, items: Sequence[spellout.Derivative]
) -> str:
    """The text report (module docstring)."""
    classes = Counter(
        (tau, base_class(record, tau))
        for record in records.roots.values()
        if record.dorsal
        for tau in attest.TAUS
    )
    lines = [
        "Italian roots without alternation evidence x derivative "
        "palatalization (design R21)",
        "=" * 72,
        "",
        "Base classes from the corpus at tau (E0 masculine plurals of the",
        "root): alternating (a palatalized -i form attested), non-",
        "alternating (attested -i forms keep the velar), no evidence (no",
        "-i form attested: feminine-only, no plural, plural below tau).",
        "Units: base root x suffix x derivative stratum; P = only",
        "palatalized derivatives, N = only velar, both. Derivatives:",
        "accepted links attested at tau (dictionary-listed at tau =",
        "dictionary); 'named+candidates' adds the derivatives generated in",
        "both spellings and found in the headwords or the corpus (R21).",
        "",
        "Dorsal roots by class:",
    ]
    for tau in attest.TAUS:
        lines.append(
            f"  {tau:<10} "
            + ", ".join(f"{k} {classes[(tau, k)]}" for k in CLASSES)
        )
    lines.append("")
    for tau in attest.TAUS:
        for tier, tiers in TIERS.items():
            units = build_units(records, items, tau, tiers)
            lines += section(
                units, f"tau = {tau}; links: {tier}", tau in FULL_TAUS
            )
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# The stage
# ---------------------------------------------------------------------------


def _body(ctx: RunContext, layout: Layout) -> int:
    records = spellout.load_records(layout)
    items = spellout.load_derivatives(layout.derivations)
    text = report(records, items)
    ctx.write_tracked_text(layout.noevidence, text, sources=SOURCE)
    print(f"  no-evidence report: {len(text.splitlines()):,} lines")
    return 0


def outputs(layout: Layout) -> list[Path]:
    return [layout.noevidence]


def run(ctx: RunContext) -> int | None:
    """Run the ``noevidence`` stage (skipped when its manifest matches)."""
    layout = Layout(ctx)
    return run_stage(
        ctx,
        STAGE,
        inputs=spellout.inputs(layout),
        params={"version": VERSION},
        outputs=outputs(layout),
        body=lambda: _body(ctx, layout),
    )
