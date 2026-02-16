from __future__ import annotations

import re
import math
from collections import Counter


_AMBIGUOUS_PATTERN = re.compile(r"[^ATGCatgc]+")
_RC_MAP = str.maketrans("ACGTacgt", "TGCAtgca")


def count_invalid_bases(seq: str) -> int:
    return len(_AMBIGUOUS_PATTERN.findall(seq))


def scan_extreme_gc_windows(
    sequence: str,
    window_size: int,
    step: int,
    gc_min: float,
    gc_max: float,
) -> list[tuple[int, int, float]]:
    seq = sequence.upper()
    windows: list[tuple[int, int, float]] = []
    n = len(seq)
    if window_size <= 0 or step <= 0 or n < window_size:
        return windows

    for start in range(0, n - window_size + 1, step):
        win = seq[start : start + window_size]
        gc = (win.count("G") + win.count("C")) / max(len(win), 1) * 100.0
        if gc < gc_min or gc > gc_max:
            windows.append((start, start + window_size, gc))
    return windows


def merge_intervals_with_gap(intervals: list[tuple[int, int, float]], gap: int = 0) -> list[tuple[int, int, float]]:
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda item: item[0])
    merged: list[tuple[int, int, float]] = []
    cur_start, cur_end, cur_score = intervals[0]
    for start, end, score in intervals[1:]:
        if start <= cur_end + gap:
            cur_end = max(cur_end, end)
            cur_score = (cur_score + score) / 2.0
        else:
            merged.append((cur_start, cur_end, cur_score))
            cur_start, cur_end, cur_score = start, end, score
    merged.append((cur_start, cur_end, cur_score))
    return merged


def scan_homopolymers(sequence: str, at_run: int = 5, gc_run: int = 4) -> list[tuple[str, int, int]]:
    seq = sequence.upper()
    at_pattern = re.compile(rf"A{{{at_run},}}|T{{{at_run},}}", re.IGNORECASE)
    gc_pattern = re.compile(rf"G{{{gc_run},}}|C{{{gc_run},}}", re.IGNORECASE)

    hits: list[tuple[str, int, int]] = []
    for match in at_pattern.finditer(seq):
        start, end = match.span()
        hits.append((match.group(0)[0].upper(), start, end))
    for match in gc_pattern.finditer(seq):
        start, end = match.span()
        hits.append((match.group(0)[0].upper(), start, end))
    return sorted(hits, key=lambda x: x[1])


def scan_ambiguous(sequence: str) -> list[tuple[int, int]]:
    raw = [(m.start(), m.end()) for m in _AMBIGUOUS_PATTERN.finditer(sequence)]
    if not raw:
        return []
    merged: list[tuple[int, int]] = []
    start, end = raw[0]
    for s, e in raw[1:]:
        if s <= end:
            end = max(end, e)
        else:
            merged.append((start, end))
            start, end = s, e
    merged.append((start, end))
    return merged


def scan_palindromes(
    sequence: str,
    min_len: int = 8,
    max_len: int = 16,
) -> list[tuple[int, int, int]]:
    seq = sequence.upper()
    if min_len < 2:
        min_len = 2
    if max_len < min_len:
        max_len = min_len

    hits: list[tuple[int, int, int]] = []
    n = len(seq)
    for length in range(min_len, max_len + 1):
        half = length // 2
        if half == 0:
            continue
        for start in range(0, n - length + 1):
            window = seq[start : start + length]
            if _is_perfect_palindrome(window):
                hits.append((start, start + length, length))
    return hits


def scan_inverted_repeats(
    sequence: str,
    min_arm: int = 8,
    max_arm: int = 12,
    max_spacer: int = 20,
) -> list[tuple[int, int, int, int]]:
    seq = sequence.upper()
    if min_arm < 2:
        min_arm = 2
    if max_arm < min_arm:
        max_arm = min_arm
    if max_spacer < 0:
        max_spacer = 0

    n = len(seq)
    hits: list[tuple[int, int, int, int]] = []
    for arm in range(min_arm, max_arm + 1):
        for left_start in range(0, max(0, n - arm - arm) + 1):
            left = seq[left_start : left_start + arm]
            if not left:
                continue
            left_rc = left.translate(_RC_MAP)[::-1]
            right_start_limit = min(left_start + arm + max_spacer + arm, n - arm)
            for right_start in range(left_start + arm, right_start_limit + 1):
                if right_start + arm > n:
                    break
                right = seq[right_start : right_start + arm]
                if right == left_rc:
                    hits.append((left_start, right_start + arm, arm, right_start - (left_start + arm)))
    return hits


def _is_acgt(sequence: str) -> bool:
    return all(ch in {"A", "C", "G", "T"} for ch in sequence)


def scan_tandem_repeats(
    sequence: str,
    min_motif_len: int = 2,
    max_motif_len: int = 6,
    min_copies: int = 3,
    min_total_len: int = 12,
) -> list[tuple[int, int, int, int]]:
    seq = sequence.upper()
    n = len(seq)
    if min_motif_len < 1:
        min_motif_len = 1
    if max_motif_len < min_motif_len:
        max_motif_len = min_motif_len
    if min_copies < 2:
        min_copies = 2

    hits: list[tuple[int, int, int, int]] = []
    for motif_len in range(min_motif_len, max_motif_len + 1):
        i = 0
        min_run_len = max(min_total_len, motif_len * min_copies)
        while i <= n - min_run_len:
            motif = seq[i : i + motif_len]
            if len(motif) < motif_len or not _is_acgt(motif):
                i += 1
                continue

            count = 1
            j = i + motif_len
            while j + motif_len <= n and seq[j : j + motif_len] == motif:
                count += 1
                j += motif_len

            run_len = count * motif_len
            if run_len >= min_run_len and count >= min_copies:
                hits.append((i, j, motif_len, count))
                i = j
                continue

            i += 1

    hits = sorted(hits, key=lambda item: (item[0], item[1]))
    return hits


def scan_low_complexity(
    sequence: str,
    window_size: int = 30,
    step: int = 10,
    max_entropy: float = 1.2,
) -> list[tuple[int, int, float]]:
    seq = sequence.upper()
    n = len(seq)
    if window_size <= 0 or step <= 0 or n < window_size:
        return []

    raw_hits: list[tuple[int, int, float]] = []
    for start in range(0, n - window_size + 1, step):
        window = seq[start : start + window_size]
        entropy = _calc_shannon_entropy(window)
        if entropy <= max_entropy:
            raw_hits.append((start, start + window_size, entropy))

    if not raw_hits:
        return []

    return merge_intervals_with_gap(raw_hits, gap=step)


def _calc_shannon_entropy(sequence: str) -> float:
    seq = [ch for ch in sequence.upper() if ch in {"A", "C", "G", "T"}]
    n = len(seq)
    if n == 0:
        return 0.0
    counts = Counter(seq)
    entropy = 0.0
    for value in counts.values():
        if value == 0:
            continue
        p = value / n
        entropy -= p * math.log2(p)
    return entropy


def _is_perfect_palindrome(sequence: str) -> bool:
    seq = sequence.upper()
    return seq == seq.translate(_RC_MAP)[::-1]
