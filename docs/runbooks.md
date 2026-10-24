# Runbooks

Two procedures change pinned inputs on purpose: re-baselining the
Romanian outputs, and moving the Italian and Czech pipelines to a new
kaikki.org dump. Both end with a commit that holds nothing else, so the
change can be reviewed and reverted on its own. A third procedure
rebuilds the cross-language comparison after any of the three languages'
outputs change (§3).

Rule numbers (R1, R2, …) refer to
[`multilingual_design.md`](multilingual_design.md).

## 1. Re-baselining Romanian

`tests/test_romanian_replay.py` pins the Romanian pipeline byte for byte.
It records the SHA-256 of its inputs (`data/romanian_lexicon_raw_dex.csv`,
`data/leipzig/freq/ron_news_2024_1M_freq.csv`) and of the outputs of
stage 4 (`romanian_lexicon_complete.csv`) and stage 5
(`romanian_lexicon_with_freq.csv`). It also checks that the validator
reproduces `analysis/palatalization_rule_validation.txt` and
`analysis/palatalization_mismatches.csv`.

Re-baseline only when an input changes on purpose (a new DEX QC run, a
new Leipzig release) or when a deliberate code change moves an output.
Never re-baseline to make an unexplained difference go away.

1. **Start clean.** `git status` should show no unrelated edits under
   `scripts/`, `src/` or `analysis/`. Run the suite once and note which
   replay assertions fail and why:

   ```bash
   python -m pytest -q tests/test_romanian_replay.py tests/test_repo_guards.py
   ```

2. **Check the Leipzig directory.** `data/leipzig/freq/` must hold
   exactly `ron_news_2024_1M_freq.csv` (R2), because stage 5 attaches
   every `*_freq.csv` there. `tests/test_repo_guards.py` checks this.
   Italian and Czech tables belong in `data/<language>/leipzig/freq/`.

3. **Regenerate** stages 4–6 and the validator outputs:

   ```bash
   python scripts/romanian_processor_main.py      # stage 4
   python scripts/attach_frequencies.py           # stage 5
   python scripts/validate_palatalization_rules.py
   python scripts/build_contingency_table.py
   Rscript analysis/analyze_romanian_palatalization.R   # stage 6, if R is installed
   ```

   Stage 4 calibrates thresholds in module state, so run each stage in
   its own interpreter, as above. Do not run them from one Python
   session.

4. **Review the diffs**: `git diff --stat analysis/`, then the
   validator report, the TP summaries and the contingency table. Every
   changed number needs an explanation in the commit message.

5. **Re-pin.** Recompute the hashes (`sha256sum data/romanian_lexicon_*.csv
   data/leipzig/freq/ron_news_2024_1M_freq.csv`). Update the constants at
   the top of `tests/test_romanian_replay.py` and the date in its
   "Baseline confirmed" comment. Then run the replay test again: it
   must pass on a fresh replay, not only against the files you just
   wrote.

6. **Commit** the regenerated `analysis/` files, the updated test
   constants and nothing else, with a message that starts with
   "Re-baseline Romanian:" and states the reason.

7. **Then rebuild the comparison** (§3) in a separate commit. The
   harmonized recount (`--lang ro`) reads
   `data/romanian_lexicon_with_freq.csv`, so its tables, and the
   cross-language summary built from them, move with the re-baseline.

## 2. Moving to a new kaikki dump

kaikki.org regenerates `raw-wiktextract-data.jsonl.gz` and the edition
extracts in place, at the same URLs. The pipeline therefore pins each
file by size and SHA-256 in `src/lexpipe/sources.py` and
`docs/sources_lock.tsv`. The lock also records the server's ETag and
Last-Modified from when the file was first verified.

1. **Notice.** `verify-sources` sends one HEAD request per kaikki file:

   ```bash
   python scripts/lang_pipeline.py --lang it verify-sources
   python scripts/lang_pipeline.py --lang cs verify-sources
   ```

   The line `upstream changed since <date>: ETag … -> …` means a new
   extraction is being served. Your local files are still the pinned
   ones, and nothing breaks until you choose to move.

2. **Download beside the old dump.** Keep the old file until the move
   is reviewed (`mv data/sources/kaikki/raw-wiktextract-data.jsonl.gz
   …/raw-wiktextract-data.<old-date>.jsonl.gz`). Then update `size` for
   the source in `src/lexpipe/sources.py` to the new `Content-Length`.
   Run `fetch`, which downloads resumably and refuses a file of the
   wrong size:

   ```bash
   python scripts/lang_pipeline.py --lang it fetch
   ```

3. **Re-pin.** Delete the source's row from `docs/sources_lock.tsv` and
   run `verify-sources`. It records the new SHA-256, ETag and
   Last-Modified (status `recorded`). Do not edit hashes by hand.

4. **Re-split with a manifest and check the counts.** Split the dump
   per language (top-level `lang_code`; the byte prefilter tolerates
   any JSON spacing):

   ```bash
   python -c "import sys; sys.path.insert(0, 'src'); \
     from lexpipe import kaikki; \
     print(kaikki.split_raw_dump('data/sources/kaikki/raw-wiktextract-data.jsonl.gz', \
                                 'data/sources/kaikki', ['it', 'cs']))"
   ```

   Compare the printed entry counts with the previous ones: 624,122
   Italian and 72,441 Czech lines in the 2026-10-01 dump. Investigate
   any change above 2% before going on: it usually means a
   serialisation or schema change, not a lexicographic one. The
   plugins' extraction stages record the split files' SHA-256 in their
   manifests (`data/<language>/manifests/*.manifest.tsv`), so every
   derived stage re-runs against the new split.

   The split files are registered as **derived** sources,
   `kaikki_raw_it` and `kaikki_raw_cs` in `src/lexpipe/sources.py`.
   `fetch` never downloads them: it prints the split command instead.
   Each entry pins the file's size and SHA-256, and the SHA-256 of the
   `kaikki_raw` it was split from (`derived_from_sha256`). After
   re-splitting, update `size`, `sha256` and `derived_from_sha256` in
   both entries, using `sha256sum` and the new `kaikki_raw` row of the
   lock. Delete their two lock rows, then run `verify-sources`. A split
   left over from the old dump is reported as `stale`.

5. **Run the it-noun drift test** (R12, once the Italian plugin has
   it). The port of `Module:it-common`'s `make_plural` must reproduce
   kaikki's plural for every noun without an explicit plural argument.
   Each mismatch is either a parser bug or a change of the Wiktionary
   module. Record the module revision you ported against, and fix the
   port or the parser before going on.

6. **Check the curated guards.** `curated/<language>/*.csv` rows are
   keyed by `lexeme_id` and carry guard fields (expected head-template
   arguments or plural set). A new dump can change a lexeme_id (its
   templates changed, or a homograph appeared) or invalidate a guard.
   The loader stops on every row whose guard no longer matches; fix
   those rows by hand. Never let the pipeline rewrite `curated/`.

7. **Re-run and review.** Run the plugins
   (`python scripts/lang_pipeline.py --lang it`, likewise `cs`), then
   the cross-language step (§3), and diff the tracked
   `analysis/<language>/` and `analysis/comparison/` tables.

8. **Write a re-baseline note** in `docs/<language>/` (or the commit
   message). It records the old and new dump dates, ETag,
   Last-Modified, SHA-256, the per-language entry counts, the drift-test
   result, the curated rows you changed and the main table changes.
   Commit `src/lexpipe/sources.py`, `docs/sources_lock.tsv`, the note
   and the regenerated tracked tables together.

## 3. Rebuilding the cross-language comparison

`analysis/comparison/README.md` (the cross-language report) and
`analysis/comparison/crosslang_summary.csv` are built by
`lexpipe.compare` from the tracked outputs of the three languages. Its
tables are generated; its prose (the methods, the check of the
learnability paragraph, the contradictions) is written by hand and
cites the tables' numbers. Rebuild it whenever a language's tracked
outputs change: a re-baseline (§1), a new dump (§2), or a fix to a
plugin.

1. **Bring the languages up to date**, in this order. Each run skips the
   stages whose manifests still match.

   ```bash
   python scripts/lang_pipeline.py --lang it
   python scripts/lang_pipeline.py --lang cs
   python scripts/lang_pipeline.py --lang ro
   ```

   The Logical Phonology validators are not pipeline stages. Re-run
   them when the data or the analyses change:

   ```bash
   python scripts/validate_lang_rules.py --lang it
   python scripts/validate_lang_rules.py --lang it --variant stress
   ```

   Likewise run `--variant all_underspecified` for Italian, and the
   Czech default with its variants `all_underspecified`, `lexical_loo`
   and `spellout`. The comparison reads their headline lines. Then
   rebuild the cross-variant summaries that `docs/<language>/phonology.md`
   quotes, and update the numbers there that move:

   ```bash
   PYTHONPATH=src python -m phonology.analyses.italian_palatalization
   PYTHONPATH=src python -m phonology.analyses.czech_palatalization
   ```

   `tests/test_validation_summaries.py` fails while a summary or a
   phonology note lags the per-variant reports.

2. **Rebuild and check.**

   ```bash
   PYTHONPATH=src python -m lexpipe.compare
   PYTHONPATH=src python -m lexpipe.compare --check   # must print nothing, exit 0
   ```

   The step reads only tracked files and takes about a second. It is
   skipped while its manifest
   (`data/romanian_compare/manifests/crosslang.manifest.tsv`) matches,
   and `--force` rebuilds anyway. It stops with `ReportFormatError`
   when a report's format changed. In that case update the matching
   parser in `src/lexpipe/compare.py` and its excerpt test in
   `tests/lexpipe/test_compare.py`; never edit the generated block by
   hand. `tests/lexpipe/test_compare.py::test_tracked_outputs_are_current`
   fails while the outputs are stale.

3. **Review the prose against the new numbers.**
   - Diff `crosslang_summary.csv` and the generated block.
   - Search the hand-written sections of the README (the Summary, the
     learnability section, the contradictions) for every number whose
     row changed, and update it.
   - A changed verdict (a rule that now passes or fails, a learner share
     that crosses 50%) needs the verdict table at the top revisited too.

4. **Commit** the regenerated `analysis/comparison/` files and the
   edited prose together, with a message that names the language
   change behind them.
