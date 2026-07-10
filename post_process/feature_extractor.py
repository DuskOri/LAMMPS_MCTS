"""从片段序列中提取简单结构特征。"""

try:
    from mcts.fragment_registry import get_fragment_registry
except ImportError:
    get_fragment_registry = None


def extract_sequence_features(sequence):
    """把一条候选片段序列整理成可写表的特征字典。"""
    fragments = [fragment for fragment in sequence if fragment not in ("Start", "End")]
    length = len(fragments)
    count_a = fragments.count("A")
    count_ap = fragments.count("A'")
    count_b = fragments.count("B")
    count_c = fragments.count("C")

    family_counts = {}
    role_counts = {}
    if get_fragment_registry is not None:
        registry = get_fragment_registry()
        for fragment in fragments:
            if not registry.has_fragment(fragment):
                continue
            spec = registry.metadata(fragment)
            family_counts[spec.family] = family_counts.get(spec.family, 0) + 1
            for role in spec.roles:
                role_counts[role] = role_counts.get(role, 0) + 1

    imide_count = role_counts.get("imide", count_a + count_ap)
    aromatic_count = role_counts.get("aromatic", count_b)
    flexible_count = role_counts.get("flexible", count_c)
    urethane_count = role_counts.get("urethane", 0)
    phenolic_count = role_counts.get("phenolic", 0)
    rigid_count = role_counts.get("rigid", 0)
    family_text = ",".join(
        f"{family}:{count}" for family, count in sorted(family_counts.items())
    )

    return {
        "length": length,
        "count_A": count_a,
        "count_A_prime": count_ap,
        "count_B": count_b,
        "count_C": count_c,
        "count_urethane": urethane_count,
        "count_phenolic": phenolic_count,
        "count_rigid": rigid_count,
        "count_flexible": flexible_count,
        "imide_ratio": _safe_divide(imide_count, length),
        "urethane_ratio": _safe_divide(urethane_count, length),
        "phenolic_ratio": _safe_divide(phenolic_count, length),
        "aromatic_ratio": _safe_divide(aromatic_count, length),
        "flexible_ratio": _safe_divide(flexible_count, length),
        "family_text": family_text,
        "sequence_text": "-".join(fragments),
    }


def batch_extract_features(sequences):
    """批量提取特征，用于候选列表统一写出。"""
    return [extract_sequence_features(sequence) for sequence in sequences]


def _safe_divide(numerator, denominator):
    """避免空序列造成除零错误。"""
    if denominator == 0:
        return 0.0
    return numerator / denominator
