"""MCTS 中用于描述聚合物片段序列的状态。"""

from dataclasses import dataclass, field

try:
    from .fragment_graph import get_allowed_fragments
except ImportError:
    from fragment_graph import get_allowed_fragments


@dataclass(frozen=True)
class PolymerState:
    """保存一条正在组装的片段路径。"""

    sequence: tuple[str, ...] = field(default_factory=lambda: ("Start",))
    max_steps: int = 8

    @property
    def current_fragment(self):
        """当前路径最后一个片段。"""
        return self.sequence[-1]

    @property
    def fragments(self):
        """去掉 Start 后给外部模块使用的片段列表。"""
        return [fragment for fragment in self.sequence if fragment != "Start"]

    def is_terminal(self):
        """遇到 End 或达到最大步数时，本轮组装结束。"""
        return self.current_fragment == "End" or len(self.fragments) >= self.max_steps

    def legal_actions(self):
        """给出当前状态还能选择的下一步片段。"""
        if self.is_terminal():
            return []

        actions = get_allowed_fragments(self.current_fragment)
        if len(self.fragments) >= self.max_steps - 1 and "End" in actions:
            return ["End"]
        return actions

    def take_action(self, action):
        """执行一次片段选择，返回新的状态对象。"""
        if action not in self.legal_actions():
            raise ValueError(f"Illegal action {action!r} from {self.current_fragment!r}")
        return PolymerState(sequence=self.sequence + (action,), max_steps=self.max_steps)

    def to_mcts_output(self):
        """转换为聚合物生成模块需要的输出格式。"""
        fragments = self.fragments
        if not fragments or fragments[-1] != "End":
            fragments = fragments + ["End"]
        return fragments
