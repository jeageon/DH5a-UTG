from __future__ import annotations

from typing import Iterable, Optional

from ..models.data_schemas import NegativeFeature


def dedupe_features(features: Iterable[NegativeFeature]) -> list[NegativeFeature]:
    unique = {}
    for feature in features:
        key = (
            feature.feature_type,
            feature.start,
            feature.end,
            feature.source,
            feature.description,
            feature.strand,
        )
        if key not in unique:
            unique[key] = feature
    return list(unique.values())


def merge_by_type(
    features: Iterable[NegativeFeature],
    merge_gaps: Optional[dict[str, int]] = None,
) -> list[NegativeFeature]:
    merge_gaps = merge_gaps or {}
    feature_map: dict[str, list[NegativeFeature]] = {}
    for feature in features:
        feature_map.setdefault(feature.feature_type, []).append(feature)

    merged: list[NegativeFeature] = []
    for feature_type, bucket in feature_map.items():
        bucket_sorted = sorted(bucket, key=lambda x: (x.start, x.end))
        gap = merge_gaps.get(feature_type, 0)
        if not bucket_sorted:
            continue
        cur = bucket_sorted[0]
        cur_start = cur.start
        cur_end = cur.end
        cur_score = cur.score
        cur_attrs = dict(cur.attributes)
        cur_desc = cur.description
        for item in bucket_sorted[1:]:
            if item.start <= cur_end + gap:
                cur_end = max(cur_end, item.end)
                cur_score = max(s for s in (cur_score, item.score) if s is not None) if (
                    cur_score is not None or item.score is not None
                ) else None
                if item.description != cur_desc:
                    cur_desc = f"{cur_desc}; {item.description}"
                cur_attrs.update({k: v for k, v in item.attributes.items() if k not in cur_attrs})
            else:
                merged.append(
                    NegativeFeature(
                        feature_type=feature_type,
                        start=cur_start,
                        end=cur_end,
                        description=cur_desc,
                        source=cur.source,
                        score=cur_score,
                        strand=cur.strand,
                        attributes=cur_attrs,
                    )
                )
                cur_start = item.start
                cur_end = item.end
                cur_score = item.score
                cur_desc = item.description
                cur_attrs = dict(item.attributes)
        merged.append(
            NegativeFeature(
                feature_type=feature_type,
                start=cur_start,
                end=cur_end,
                description=cur_desc,
                source=cur.source,
                score=cur_score,
                strand=cur.strand,
                attributes=cur_attrs,
            )
        )
    return dedupe_features(merged)
