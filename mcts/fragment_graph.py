"""片段连接规则的读取接口。"""

try:
    from .fragment_registry import get_fragment_registry
except ImportError:
    from fragment_registry import get_fragment_registry


def get_allowed_fragments(fragment):
    """根据当前片段，返回下一步可以接上的候选片段。"""
    return get_fragment_registry().allowed_fragments(fragment)


def is_valid_transition(source, target):
    """判断 source -> target 这一步是否满足预设的化学连接图。"""
    return target in get_allowed_fragments(source)


def validate_sequence(sequence):
    """检查一整条片段序列是否能在 PolyGRAPH 中走通。"""
    if not sequence:
        return False

    current = "Start"
    for fragment in sequence:
        if not is_valid_transition(current, fragment):
            return False
        current = fragment
        if current == "End":
            break

    return current == "End"


def current_graph():
    """返回当前运行时使用的片段连接图，主要供调试和导出使用。"""
    return get_fragment_registry().graph()
