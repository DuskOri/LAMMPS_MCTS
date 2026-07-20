"""片段连接规则的读取接口。"""

try:
    from .fragment_registry import get_fragment_registry
except ImportError:
    from fragment_registry import get_fragment_registry


def get_allowed_fragments(fragment):
    """根据当前片段，返回下一步可以接上的候选片段。"""
    return get_fragment_registry().allowed_fragments(fragment)
