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


def _motif_distance(fragment: str, motif: str) -> int:
    if len(fragment) != len(motif):
        return math.inf  # type: ignore[return-value]
    if not fragment:
        return 0
    return sum(1 for a, b in zip(fragment, motif) if a != b)


def scan_promoter_like(
    sequence: str,
    min_spacer: int = 13,
    max_spacer: int = 19,
    min_quality: float = 0.75,
) -> list[tuple[int, int, str]]:
    seq = sequence.upper()
    if len(seq) < 30:
        return []

    # Sigma70-like bacterial promoter motifs (consensus-flexible search).
    p35_candidates = ("TTGACA", "TTGACG", "TTGACT", "TTGAAA", "TTGAGT", "CTGACA")
    p10_candidates = ("TATAAT", "TATAAC", "TATTAT", "TATGAT", "TATACT", "TATCAT", "TATGAA")

    hits: list[tuple[int, int, str]] = []
    seen: set[tuple[int, int, str]] = set()
    for p35_start in range(0, len(seq) - 20):
        p35 = seq[p35_start : p35_start + 6]
        if len(p35) < 6:
            continue
        best_35 = min(_motif_distance(p35, c) for c in p35_candidates)
        if best_35 > 2:
            continue

        for spacer in range(min_spacer, max_spacer + 1):
            p10_start = p35_start + 6 + spacer
            p10_end = p10_start + 6
            if p10_end > len(seq):
                break

            p10 = seq[p10_start:p10_end]
            best_10 = min(_motif_distance(p10, c) for c in p10_candidates)
            if best_10 > 2:
                continue

            # prefer canonical spacing and stronger -10 motif
            quality = (6 - best_35 + 6 - best_10) / 12.0
            local_window = seq[max(0, p10_start - 6) : min(len(seq), p10_start + 12)]
            if not local_window:
                continue
            at_frac = (local_window.count("A") + local_window.count("T")) / len(local_window)
            if quality < min_quality or at_frac < 0.55:
                continue

            left = p35_start
            right = p10_end
            motif = f"{seq[left:right]} (spacer={spacer}, -35 mism={best_35}, -10 mism={best_10}, AT={at_frac:.2f})"
            key = (left, right, motif)
            if key not in seen:
                seen.add(key)
                hits.append((left, right, motif))

    hits.sort(key=lambda item: item[0])
    return hits


def scan_rbs_like(sequence: str, min_downstream: int = 3, max_downstream: int = 18) -> list[tuple[int, int, str]]:
    seq = sequence.upper()
    if len(seq) < 4:
        return []
    motifs = [
        "AGGAGG",
        "GAGGAG",
        "GGAGG",
        "AGGA",
    ]
    start_codons = ("ATG", "GTG", "TTG", "CTG")
    hits: list[tuple[int, int, str]] = []
    for motif in motifs:
        idx = 0
        while True:
            idx = seq.find(motif, idx)
            if idx == -1:
                break
            region_end = min(len(seq), idx + len(motif) + max_downstream)
            downstream = seq[idx + len(motif) : region_end]
            has_start = False
            for codon in start_codons:
                pos = downstream.find(codon)
                if pos != -1 and min_downstream <= pos <= max_downstream:
                    has_start = True
                    break
            if has_start:
                hits.append((idx, idx + len(motif), motif))
            idx += 1
    return sorted(hits, key=lambda item: item[0])


def scan_terminator_like(
    sequence: str,
    min_arm: int = 6,
    max_arm: int = 15,
    max_spacer: int = 12,
    poly_t_min: int = 6,
) -> list[tuple[int, int, str]]:
    seq = sequence.upper()
    if len(seq) < 20:
        return []

    hits: list[tuple[int, int, str]] = []
    seen: set[tuple[int, int, str]] = set()
    repeats = scan_inverted_repeats(seq, min_arm=min_arm, max_arm=max_arm, max_spacer=max_spacer)

    # 1) Canonical rho-independent style: stem-loop (inverted repeat) + downstream poly-T run
    for start, end, arm, spacer in repeats:
        scan_end = min(len(seq), end + 45)
        tail = seq[end : scan_end]
        match = re.search(rf"T{{{poly_t_min},}}", tail)
        if match is None:
            continue
        tail_start = end + match.start()
        tail_end = end + match.end()
        # poly-T should start shortly after loop; this catches truncated contexts in windowed regions
        if tail_start - end > 30:
            continue

        # avoid weak stems: very short repeats need stronger T stretch
        if arm < 8 and (tail_end - tail_start) < 7:
            continue

        motif = f"hairpin_arm={arm}_spacer={spacer}_polyT={tail_start}-{tail_end}"
        key = (start, tail_end, motif)
        if key not in seen:
            seen.add(key)
            hits.append((start, tail_end, motif))

    # 2) Fallback for terminator-like poly-T with nearby local inverted context
    repeats_for_fallback = sorted(repeats, key=lambda item: item[1], reverse=True)
    for match in re.finditer(rf"T{{{poly_t_min + 1},}}", seq):
        tail_start = match.start()
        tail_end = match.end()
        if tail_end - tail_start < 7:
            continue

        poly_len = tail_end - tail_start
        pre_window = seq[max(0, tail_start - 45) : tail_start]
        at_frac = (pre_window.count("A") + pre_window.count("T")) / max(len(pre_window), 1)
        if at_frac < 0.68:
            continue

        has_stem = any(
            r_start < tail_start <= r_end + 20 and (tail_start - r_end) <= 25
            for r_start, r_end, r_arm, _ in repeats_for_fallback
            if r_arm >= 6 and r_end <= tail_start
        )
        if not has_stem:
            continue

        # keep endpoint a little broader
        hit_end = min(len(seq), tail_end + 5)
        motif = f"polyT={poly_len}b (tail={tail_start}-{tail_end},ATwin={at_frac:.2f})"
        key = (tail_start, hit_end, motif)
        if key not in seen:
            seen.add(key)
            hits.append((tail_start, hit_end, motif))

    hits.sort(key=lambda item: item[0])
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
