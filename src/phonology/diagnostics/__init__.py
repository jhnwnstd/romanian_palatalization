"""Diagnostic tools for the phonology framework.

Each module here answers a different question about a failing
derivation:

- :mod:`compare` — one canonical IPA comparison stack for
  predicted vs attested; every other module routes through it.
- :mod:`explain` — per-rule firing/blocking narratives (why did
  rule X fire or not fire on segment Y?).
- :mod:`distance` — numeric distance-to-working metric so
  mismatches can be sorted "close" to "far".
- :mod:`ordering` — is there a permutation of the rule tuple that
  would produce the expected SR? Distinguishes ordering bugs from
  content bugs.
- :mod:`perturb` — is there a single-feature tweak (a patch, a
  clear-set entry, a rule target/supply feature) that would flip
  the derivation from failing to succeeding?
- :mod:`categorise` — improved mismatch categorisation that peels
  off data-side transcription noise before naming residual rule
  failures.

The diagnostics package depends on the core engine but not the other
way around. :mod:`compare` and :mod:`categorise` also hold the
Romanian defaults of two :class:`~phonology.profile.AnalysisProfile`
hooks (the normaliser cascade and the categoriser), so
:mod:`phonology.profile`, :mod:`phonology.g2p` and
:mod:`phonology.validation` import from this package.
"""

from __future__ import annotations
