"""基于候选序列特征的简单解释性分析。"""

from .feature_extractor import extract_sequence_features


def explain_candidate(sequence, score=None):
    """给单条序列生成一段可读的结构解释。"""
    features = extract_sequence_features(sequence)
    notes = []

    if features["aromatic_ratio"] >= 0.4:
        notes.append("芳环片段比例较高，链段刚性较强。")
    if features["imide_ratio"] >= 0.3:
        notes.append("酰亚胺片段占比较高，可能有利于提高主链极性和规整性。")
    if features["flexible_ratio"] >= 0.3:
        notes.append("柔性醚键比例较高，可能降低链段刚性。")
    if not notes:
        notes.append("序列组成较均衡，需要结合后续 MD 结果判断。")

    return {
        "sequence": " ".join(sequence),
        "score": score,
        "features": features,
        "notes": notes,
    }


def summarize_candidates(candidates):
    """对候选列表逐条生成解释结果。"""
    return [
        explain_candidate(candidate["sequence"], candidate.get("score"))
        for candidate in candidates
    ]
