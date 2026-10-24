"""Tests for lexpipe.sources: registry integrity, verify, fetch.

No test touches the network: downloads go through a fake session
that serves bytes from memory and records the requests it gets.
"""

from __future__ import annotations

import dataclasses
import hashlib
import sys
from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path, PurePosixPath

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lexpipe import sources, spec  # noqa: E402
from lexpipe.sources import (  # noqa: E402
    DownloadError,
    Source,
    VerifyStatus,
)

# Sizes verified on 2026-10-23 (docs/multilingual_design.md).
EXPECTED_SIZES = {
    "kaikki_raw": 3_082_192_147,
    "kaikki_raw_it": 665_207_473,  # derived: the 'it' split of kaikki_raw
    "kaikki_raw_cs": 223_701_337,  # derived: the 'cs' split of kaikki_raw
    "kaikki_it_extract": 43_849_715,
    "kaikki_cs_extract": 39_633_204,
    "morphit": 19_068_565,
    "morfflex": 250_486_476,
    "derinet": 439_063_096,
    "czetyl": 191_934,
    "cnk_syn2015_lemma": 1_273_063,
    "cnk_syn2015_word": 3_089_732,
    "leipzig_ita_news_2024_1M": 261_325_052,
    "leipzig_ces_news_2024_1M": 214_278_386,
}

PAYLOAD = bytes(range(256)) * 40  # 10,240 bytes
TOY = Source(
    name="toy",
    url="https://example.org/toy.bin",
    path="data/sources/toy/toy.bin",
    size=len(PAYLOAD),
    licence="CC0",
    citation="none",
)
TOY_REGISTRY = {"toy": TOY}


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------


def test_registry_names_are_unique() -> None:
    names = [s.name for s in sources._REGISTRY]
    assert len(names) == len(set(names)) == len(sources.SOURCES)
    assert all(key == s.name for key, s in sources.SOURCES.items())


def test_registry_entries_are_complete() -> None:
    paths = [s.path for s in sources.SOURCES.values()]
    assert len(paths) == len(set(paths))
    for s in sources.SOURCES.values():
        rel = PurePosixPath(s.path)
        assert rel.parts[:2] == ("data", "sources") and ".." not in rel.parts
        assert s.size > 0
        if s.derived:
            assert not s.url and s.recipe and s.sha256
        else:
            assert s.url.startswith("https://")
        assert s.licence and s.citation
        assert s.languages and set(s.languages) <= set(spec.LANGUAGES)


def test_registry_sizes() -> None:
    assert {n: s.size for n, s in sources.SOURCES.items()} == EXPECTED_SIZES


def test_sources_for_language() -> None:
    italian = {s.name for s in sources.sources_for("it")}
    czech = {s.name for s in sources.sources_for("cs")}
    assert italian == {
        "kaikki_raw",
        "kaikki_raw_it",
        "kaikki_it_extract",
        "morphit",
        "leipzig_ita_news_2024_1M",
    }
    assert {"kaikki_raw", "kaikki_raw_cs", "morfflex"} <= czech
    assert sources.sources_for("ro") == ()  # the comparison reads no source
    assert not italian & (czech - {"kaikki_raw"})
    assert sources.sources_for("xx") == ()


@pytest.mark.parametrize(
    "path", ["/abs/x", "data/other/x", "data/sources/../x", "data/sources"]
)
def test_source_paths_must_stay_under_data_sources(path: str) -> None:
    with pytest.raises(ValueError, match="data/sources"):
        Source(
            name="bad",
            url="https://example.org/x",
            path=path,
            size=1,
            licence="x",
            citation="x",
        )


def test_unknown_source_lists_known_names() -> None:
    with pytest.raises(KeyError, match="kaikki_raw"):
        sources.get_source("nope")


# ---------------------------------------------------------------------------
# verify
# ---------------------------------------------------------------------------


def _verify(tmp_path: Path) -> sources.VerifyResult:
    return sources.verify(
        "toy",
        root=tmp_path,
        registry=TOY_REGISTRY,
        lock_path=tmp_path / "lock.tsv",
        today=date(2026, 10, 23),
    )


def test_verify_lifecycle(tmp_path: Path) -> None:
    target = sources.local_path(TOY, tmp_path)
    lock = tmp_path / "lock.tsv"

    missing = _verify(tmp_path)
    assert missing.status is VerifyStatus.MISSING and not missing.ok
    assert not lock.exists()

    target.parent.mkdir(parents=True)
    target.write_bytes(PAYLOAD[:-1])
    short = _verify(tmp_path)
    assert short.status is VerifyStatus.SIZE_MISMATCH and not short.ok
    assert not lock.exists()  # a wrong-size file is never recorded

    target.write_bytes(PAYLOAD)
    digest = hashlib.sha256(PAYLOAD).hexdigest()
    first = _verify(tmp_path)
    assert first.status is VerifyStatus.RECORDED and first.ok
    assert first.sha256 == digest
    assert lock.read_text(encoding="utf-8") == (
        "name\tpath\tbytes\tsha256\trecorded\tetag\tlast_modified\n"
        f"toy\tdata/sources/toy/toy.bin\t{len(PAYLOAD)}\t{digest}\t"
        "2026-10-23\t\t\n"
    )

    second = _verify(tmp_path)
    assert second.status is VerifyStatus.PINNED and second.ok

    target.write_bytes(PAYLOAD[::-1])  # same size, other bytes
    changed = _verify(tmp_path)
    assert changed.status is VerifyStatus.SHA256_MISMATCH and not changed.ok
    assert sources.read_lock(lock)["toy"].sha256 == digest


def test_lock_keeps_other_entries_sorted(tmp_path: Path) -> None:
    lock = tmp_path / "lock.tsv"
    other = sources.LockEntry("aaa", "data/sources/a", 1, "0" * 64, "x")
    sources.write_lock({"aaa": other}, lock)
    target = sources.local_path(TOY, tmp_path)
    target.parent.mkdir(parents=True)
    target.write_bytes(PAYLOAD)
    assert _verify(tmp_path).status is VerifyStatus.RECORDED
    assert list(sources.read_lock(lock)) == ["aaa", "toy"]


def test_read_lock_checks_header(tmp_path: Path) -> None:
    lock = tmp_path / "lock.tsv"
    lock.write_text("name\tsha\n", encoding="utf-8")
    with pytest.raises(ValueError, match="expected columns"):
        sources.read_lock(lock)


# ---------------------------------------------------------------------------
# fetch / download through a fake session
# ---------------------------------------------------------------------------


@dataclass
class FakeResponse:
    status_code: int
    body: bytes = b""
    headers: Mapping[str, str] = field(default_factory=dict)
    closed: bool = False

    def iter_content(self, chunk_size: int) -> Iterator[bytes]:
        for start in range(0, len(self.body), chunk_size):
            yield self.body[start : start + chunk_size]

    def close(self) -> None:
        self.closed = True


@dataclass
class FakeSession:
    """Serves PAYLOAD like a server honouring Range, unless told not to."""

    status: int = 200
    honour_range: bool = True
    truncate_to: int | None = None
    requests: list[dict[str, str]] = field(default_factory=list)

    def get(
        self,
        url: str,
        *,
        headers: Mapping[str, str],
        stream: bool,
        timeout: tuple[float, float],
    ) -> FakeResponse:
        assert stream
        self.requests.append(dict(headers))
        if self.status != 200:
            return FakeResponse(self.status, b"<html>Bot Check</html>")
        body = PAYLOAD[: self.truncate_to]
        range_header = headers.get("Range")
        if range_header and self.honour_range:
            start = int(range_header.removeprefix("bytes=").rstrip("-"))
            return FakeResponse(
                206,
                body[start:],
                {
                    "Content-Range": f"bytes {start}-{len(body) - 1}/{len(body)}"
                },
            )
        return FakeResponse(200, body)


def _fetch(tmp_path: Path, session: FakeSession) -> Path:
    return sources.fetch(
        "toy",
        root=tmp_path,
        registry=TOY_REGISTRY,
        session=session,
        log=lambda _: None,
    )


def _part(tmp_path: Path) -> Path:
    target = sources.local_path(TOY, tmp_path)
    return target.with_name(target.name + ".part")


def test_fetch_downloads_a_missing_file(tmp_path: Path) -> None:
    session = FakeSession()
    path = _fetch(tmp_path, session)
    assert path.read_bytes() == PAYLOAD
    assert not _part(tmp_path).exists()
    (headers,) = session.requests
    assert "Range" not in headers
    assert headers["Accept-Encoding"] == "identity"
    assert sources.CONTACT_URL in headers["User-Agent"]


def test_fetch_resumes_a_partial_download(tmp_path: Path) -> None:
    part = _part(tmp_path)
    part.parent.mkdir(parents=True)
    part.write_bytes(PAYLOAD[:4000])
    session = FakeSession()
    assert _fetch(tmp_path, session).read_bytes() == PAYLOAD
    assert session.requests[0]["Range"] == "bytes=4000-"


def test_fetch_restarts_when_range_is_ignored(tmp_path: Path) -> None:
    part = _part(tmp_path)
    part.parent.mkdir(parents=True)
    part.write_bytes(b"x" * 4000)
    session = FakeSession(honour_range=False)
    assert _fetch(tmp_path, session).read_bytes() == PAYLOAD


def test_fetch_finishes_a_complete_part_file(tmp_path: Path) -> None:
    part = _part(tmp_path)
    part.parent.mkdir(parents=True)
    part.write_bytes(PAYLOAD)
    session = FakeSession()
    assert _fetch(tmp_path, session).read_bytes() == PAYLOAD
    assert session.requests == []


def test_fetch_leaves_a_present_file_alone(tmp_path: Path) -> None:
    target = sources.local_path(TOY, tmp_path)
    target.parent.mkdir(parents=True)
    target.write_bytes(PAYLOAD)
    session = FakeSession()
    assert _fetch(tmp_path, session) == target
    assert session.requests == []


def test_fetch_never_touches_a_wrong_size_file(tmp_path: Path) -> None:
    target = sources.local_path(TOY, tmp_path)
    target.parent.mkdir(parents=True)
    target.write_bytes(PAYLOAD[:100])  # e.g. another download in progress
    session = FakeSession()
    with pytest.raises(DownloadError, match="left untouched"):
        _fetch(tmp_path, session)
    assert target.read_bytes() == PAYLOAD[:100]
    assert session.requests == []


def test_fetch_reports_http_errors(tmp_path: Path) -> None:
    with pytest.raises(DownloadError, match="HTTP 402"):
        _fetch(tmp_path, FakeSession(status=402))
    assert not sources.local_path(TOY, tmp_path).exists()
    assert not _part(tmp_path).exists()


def test_short_download_is_kept_for_resuming(tmp_path: Path) -> None:
    with pytest.raises(DownloadError, match="kept for resuming"):
        _fetch(tmp_path, FakeSession(truncate_to=5000))
    assert not sources.local_path(TOY, tmp_path).exists()
    assert _part(tmp_path).read_bytes() == PAYLOAD[:5000]
    # The next run resumes where the last one stopped.
    session = FakeSession()
    assert _fetch(tmp_path, session).read_bytes() == PAYLOAD
    assert session.requests[0]["Range"] == "bytes=5000-"


def test_user_agent_contact(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(sources.CONTACT_ENV, raising=False)
    assert sources.CONTACT_URL in sources.user_agent()
    monkeypatch.setenv(sources.CONTACT_ENV, "lab@example.org")
    assert sources.user_agent().endswith("; lab@example.org)")


# ---------------------------------------------------------------------------
# Registry pins, rolling sources and the upstream columns of the lock
# ---------------------------------------------------------------------------


def test_morphit_comes_from_the_pinned_mirror() -> None:
    morphit = sources.get_source("morphit")
    assert morphit.url == (
        "https://raw.githubusercontent.com/giodegas/morphit-lemmatizer/"
        "master/master/morph-it_048_utf8.txt"
    )
    assert morphit.path == "data/sources/morphit/morph-it_048_utf8.txt"
    assert morphit.sha256 == (
        "7375216351ee59cee4326b455503f3a6c4e9542cd63ecb86ff5cd58ca5c8e800"
    )
    assert "CC BY-SA 2.0" in morphit.licence and "LGPL" in morphit.licence
    assert "docs.sslmit.unibo.it" in morphit.notes
    assert "data/sources/morphit/morph-it.tgz" in morphit.notes
    assert "2,617,548" in morphit.notes


def test_only_kaikki_files_are_rolling() -> None:
    rolling = {s.name for s in sources.SOURCES.values() if s.rolling}
    assert rolling == {"kaikki_raw", "kaikki_it_extract", "kaikki_cs_extract"}


def test_leipzig_terms_are_marked_unverified() -> None:
    for name in ("leipzig_ita_news_2024_1M", "leipzig_ces_news_2024_1M"):
        assert "unverified" in sources.get_source(name).licence


def test_source_sha256_must_be_hex() -> None:
    with pytest.raises(ValueError, match="64 hex"):
        Source(
            name="bad",
            url="https://example.org/x",
            path="data/sources/x/x",
            size=1,
            licence="x",
            citation="x",
            sha256="abc",
        )


def _write_toy(tmp_path: Path, payload: bytes = PAYLOAD) -> Path:
    target = sources.local_path(TOY, tmp_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(payload)
    return target


def test_registry_pin_is_checked_before_the_lock(tmp_path: Path) -> None:
    _write_toy(tmp_path)
    digest = hashlib.sha256(PAYLOAD).hexdigest()
    lock = tmp_path / "lock.tsv"
    for pin, expected in (
        ("0" * 64, VerifyStatus.SHA256_MISMATCH),
        (digest, VerifyStatus.RECORDED),
    ):
        pinned = dataclasses.replace(TOY, sha256=pin)
        result = sources.verify(
            "toy", root=tmp_path, registry={"toy": pinned}, lock_path=lock
        )
        assert result.status is expected
    assert sources.read_lock(lock)["toy"].sha256 == digest


def test_legacy_lock_is_read_and_upgraded(tmp_path: Path) -> None:
    _write_toy(tmp_path)
    digest = hashlib.sha256(PAYLOAD).hexdigest()
    lock = tmp_path / "lock.tsv"
    lock.write_text(
        "name\tpath\tbytes\tsha256\trecorded\n"
        f"toy\tdata/sources/toy/toy.bin\t{len(PAYLOAD)}\t{digest}\t"
        "2026-10-01\n",
        encoding="utf-8",
    )
    entry = sources.read_lock(lock)["toy"]
    assert (entry.etag, entry.last_modified) == ("", "")
    remote = sources.RemoteMeta('"e1"', "Mon, 19 Oct 2026", len(PAYLOAD))
    result = sources.verify(
        "toy",
        root=tmp_path,
        registry=TOY_REGISTRY,
        lock_path=lock,
        remote=remote,
    )
    assert result.status is VerifyStatus.PINNED and result.ok
    assert "recorded" in result.detail and not result.upstream_changed
    text = lock.read_text(encoding="utf-8")
    assert text.splitlines()[0].endswith("\tetag\tlast_modified")
    assert text.splitlines()[1].endswith(
        '\t2026-10-01\t"e1"\tMon, 19 Oct 2026'
    )


def test_lock_row_for_another_release(tmp_path: Path) -> None:
    _write_toy(tmp_path)
    lock = tmp_path / "lock.tsv"
    digest = hashlib.sha256(PAYLOAD).hexdigest()
    old = sources.LockEntry("toy", "data/sources/toy/old.tgz", 7, digest, "x")
    sources.write_lock({"toy": old}, lock)
    result = _verify(tmp_path)
    assert result.status is VerifyStatus.LOCK_MISMATCH and not result.ok
    assert "old.tgz" in result.detail


def test_upstream_of_another_size_is_not_recorded(tmp_path: Path) -> None:
    _write_toy(tmp_path)
    lock = tmp_path / "lock.tsv"
    remote = sources.RemoteMeta('"new"', "today", len(PAYLOAD) + 1)
    result = sources.verify(
        "toy",
        root=tmp_path,
        registry=TOY_REGISTRY,
        lock_path=lock,
        remote=remote,
    )
    assert result.status is VerifyStatus.RECORDED
    assert result.upstream_changed and "another release" in result.detail
    assert sources.read_lock(lock)["toy"].etag == ""


def test_upstream_change_keeps_the_lock(tmp_path: Path) -> None:
    _write_toy(tmp_path)
    lock = tmp_path / "lock.tsv"

    def check(etag: str) -> sources.VerifyResult:
        return sources.verify(
            "toy",
            root=tmp_path,
            registry=TOY_REGISTRY,
            lock_path=lock,
            remote=sources.RemoteMeta(etag, "lm", len(PAYLOAD)),
        )

    assert check('"a"').status is VerifyStatus.RECORDED
    same = check('"a"')
    assert same.status is VerifyStatus.PINNED and not same.upstream_changed
    changed = check('"b"')
    assert changed.status is VerifyStatus.PINNED and changed.ok
    assert changed.upstream_changed and '"a" -> "b"' in changed.detail
    assert sources.read_lock(lock)["toy"].etag == '"a"'
    # Without a HEAD answer nothing is compared or rewritten.
    assert _verify(tmp_path).detail == "matches lock file"
    # Headers the server leaves out are not differences.
    silent = sources.verify(
        "toy",
        root=tmp_path,
        registry=TOY_REGISTRY,
        lock_path=lock,
        remote=sources.RemoteMeta(),
    )
    assert not silent.upstream_changed and "sent no ETag" in silent.detail
    partial = sources.verify(
        "toy",
        root=tmp_path,
        registry=TOY_REGISTRY,
        lock_path=lock,
        remote=sources.RemoteMeta(last_modified="lm"),
    )
    assert not partial.upstream_changed


@dataclass
class FakeHead:
    status: int = 200
    headers: Mapping[str, str] = field(default_factory=dict)
    calls: list[dict[str, object]] = field(default_factory=list)

    def head(
        self,
        url: str,
        *,
        headers: Mapping[str, str],
        timeout: tuple[float, float],
        allow_redirects: bool,
    ) -> FakeResponse:
        self.calls.append(
            {"url": url, "allow_redirects": allow_redirects, **headers}
        )
        return FakeResponse(self.status, headers=self.headers)


def test_head_reads_the_upstream_headers() -> None:
    session = FakeHead(
        headers={
            "ETag": '"6ad8a097-b7b68513"',
            "last-modified": "Wed, 21 Oct 2026 07:11:19 GMT",
            "Content-Length": "3082192147",
        }
    )
    meta = sources.head("https://example.org/x", session=session)
    assert meta == sources.RemoteMeta(
        '"6ad8a097-b7b68513"', "Wed, 21 Oct 2026 07:11:19 GMT", 3_082_192_147
    )
    (call,) = session.calls
    assert call["allow_redirects"] is True
    assert sources.CONTACT_URL in str(call["User-Agent"])
    assert sources.head("https://e.org/y", session=FakeHead()).size is None


def test_head_reports_http_errors() -> None:
    with pytest.raises(DownloadError, match="HEAD .* HTTP 403"):
        sources.head("https://example.org/x", session=FakeHead(status=403))


# ---------------------------------------------------------------------------
# Derived sources (the English Wiktionary subsets)
# ---------------------------------------------------------------------------

PARENT_DIGEST = hashlib.sha256(PAYLOAD).hexdigest()
CHILD = PAYLOAD[:1000]
CHILD_DIGEST = hashlib.sha256(CHILD).hexdigest()
DERIVED = Source(
    name="toy_it",
    url="",
    path="data/sources/toy/toy_it.bin",
    size=len(CHILD),
    licence="CC0",
    citation="none",
    sha256=CHILD_DIGEST,
    kind=sources.SourceKind.DERIVED,
    derived_from="toy",
    derived_from_sha256=PARENT_DIGEST,
    recipe="python make_toy.py it",
)
DERIVED_REGISTRY = {"toy": TOY, "toy_it": DERIVED}


def test_registered_subsets() -> None:
    from lexpipe import kaikki

    raw = sources.get_source("kaikki_raw")
    for code in ("it", "cs"):
        subset = sources.get_source(f"kaikki_raw_{code}")
        assert subset.kind is sources.SourceKind.DERIVED and subset.derived
        assert subset.derived_from == "kaikki_raw"
        assert subset.languages == (code,)
        assert not subset.rolling and not subset.url
        assert PurePosixPath(subset.path) == PurePosixPath(
            kaikki.split_path(PurePosixPath(raw.path).parent, code).as_posix()
        )
        assert "split_raw_dump" in subset.recipe
        assert raw.path in subset.recipe and f"['{code}']" in subset.recipe
        assert len(subset.sha256) == len(subset.derived_from_sha256) == 64
    names = list(sources.SOURCES)
    assert names.index("kaikki_raw") < names.index("kaikki_raw_it")
    assert names.index("kaikki_raw") < names.index("kaikki_raw_cs")


@pytest.mark.parametrize(
    "changes",
    [
        {"url": "https://example.org/x"},  # derived files have no URL
        {"rolling": True},
        {"derived_from": ""},
        {"derived_from_sha256": ""},
        {"sha256": ""},
        {"recipe": "  "},
        {"derived_from": "toy_it"},  # itself
        {"derived_from_sha256": "abc"},
    ],
)
def test_derived_source_validation(changes: dict[str, object]) -> None:
    with pytest.raises(ValueError):
        dataclasses.replace(DERIVED, **changes)


def test_download_cannot_have_derived_fields() -> None:
    with pytest.raises(ValueError, match="only derived"):
        dataclasses.replace(TOY, recipe="make it")
    assert dataclasses.replace(TOY, kind="download").kind is (  # type: ignore[arg-type]
        sources.SourceKind.DOWNLOAD
    )


def test_registry_order_and_uniqueness() -> None:
    sources.check_registry([TOY, DERIVED])
    with pytest.raises(ValueError, match="before it"):
        sources.check_registry([DERIVED, TOY])
    with pytest.raises(ValueError, match="twice"):
        sources.check_registry([TOY, TOY])


def _fetch_derived(tmp_path: Path, session: FakeSession) -> Path:
    return sources.fetch(
        "toy_it",
        root=tmp_path,
        registry=DERIVED_REGISTRY,
        session=session,
        log=lambda _: None,
    )


def test_fetch_never_downloads_a_derived_source(tmp_path: Path) -> None:
    session = FakeSession()
    with pytest.raises(sources.DerivedSourceMissing) as excinfo:
        _fetch_derived(tmp_path, session)
    message = str(excinfo.value)
    assert "derived from toy" in message and "python make_toy.py it" in message
    assert "verify-sources" in message
    assert isinstance(excinfo.value, DownloadError)
    assert session.requests == []
    assert not sources.local_path(DERIVED, tmp_path).exists()

    target = sources.local_path(DERIVED, tmp_path)
    target.parent.mkdir(parents=True)
    target.write_bytes(CHILD[:10])
    with pytest.raises(DownloadError, match="make it again"):
        _fetch_derived(tmp_path, session)
    assert target.read_bytes() == CHILD[:10]  # untouched
    target.write_bytes(CHILD)
    assert _fetch_derived(tmp_path, session) == target
    assert session.requests == []


def _verify_derived(
    tmp_path: Path, parent_digest: str | None
) -> sources.VerifyResult:
    lock = tmp_path / "lock.tsv"
    if parent_digest is not None:
        parent = sources.LockEntry(
            "toy", TOY.path, TOY.size, parent_digest, "2026-10-23"
        )
        entries = sources.read_lock(lock)
        entries["toy"] = parent
        sources.write_lock(entries, lock)
    return sources.verify(
        "toy_it",
        root=tmp_path,
        registry=DERIVED_REGISTRY,
        lock_path=lock,
        today=date(2026, 10, 23),
    )


def test_verify_derived_source(tmp_path: Path) -> None:
    missing = _verify_derived(tmp_path, None)
    assert missing.status is VerifyStatus.MISSING
    assert "python make_toy.py it" in missing.detail
    target = sources.local_path(DERIVED, tmp_path)
    target.parent.mkdir(parents=True)
    target.write_bytes(CHILD)

    unchecked = _verify_derived(tmp_path, None)
    assert unchecked.status is VerifyStatus.RECORDED and unchecked.ok
    assert "no lock row yet" in unchecked.detail
    assert not unchecked.upstream_changed
    pinned = _verify_derived(tmp_path, PARENT_DIGEST)
    assert pinned.status is VerifyStatus.PINNED and pinned.ok
    assert "made from the locked toy" in pinned.detail
    assert sources.read_lock(tmp_path / "lock.tsv")["toy_it"].sha256 == (
        CHILD_DIGEST
    )


def test_derived_source_of_another_parent_release_is_stale(
    tmp_path: Path,
) -> None:
    target = sources.local_path(DERIVED, tmp_path)
    target.parent.mkdir(parents=True)
    target.write_bytes(CHILD)
    stale = _verify_derived(tmp_path, "0" * 64)
    assert stale.status is VerifyStatus.STALE and not stale.ok
    assert "make it again" in stale.detail
    assert "toy_it" not in sources.read_lock(tmp_path / "lock.tsv")
    # Its own bytes are still checked first.
    target.write_bytes(CHILD[::-1])
    wrong = _verify_derived(tmp_path, PARENT_DIGEST)
    assert wrong.status is VerifyStatus.SHA256_MISMATCH


def test_derive_instructions() -> None:
    text = sources.derive_instructions(sources.get_source("kaikki_raw_cs"))
    assert "kaikki.split_raw_dump(" in text and "['cs']" in text
    with pytest.raises(ValueError, match="downloaded"):
        sources.derive_instructions(sources.get_source("kaikki_raw"))
