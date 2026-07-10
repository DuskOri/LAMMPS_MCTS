"""聚合物片段序列的蒙特卡洛树搜索主流程。"""

import random

try:
    from .node import MCTSNode
    from .state import PolymerState
    from .fragment_registry import get_fragment_registry
except ImportError:
    from node import MCTSNode
    from state import PolymerState
    from fragment_registry import get_fragment_registry


def default_reward(fragment_sequence):
    """低热导目标的启发式打分，分数越高代表预估热导越低。"""
    fragments = [fragment for fragment in fragment_sequence if fragment != "End"]
    if not fragments:
        return 0.0

    registry = get_fragment_registry()
    rigid_count = 0
    flexible_count = 0
    hetero_count = 0
    for fragment in fragments:
        if not registry.has_fragment(fragment):
            continue
        roles = set(registry.metadata(fragment).roles)
        rigid_count += int("rigid" in roles or "aromatic" in roles)
        flexible_count += int("flexible" in roles or "ether" in roles)
        hetero_count += int(bool(roles.intersection({"urethane", "imide", "phenolic"})))
    length_count = len(fragments)

    estimated_kappa = (
        1.0
        + 0.22 * rigid_count
        + 0.06 * length_count
        + 0.04 * hetero_count
        - 0.28 * flexible_count
    )
    estimated_kappa = max(estimated_kappa, 0.05)
    return 1.0 / estimated_kappa


class MCTSEngine:
    """负责执行选择、扩展、模拟、回传四个步骤。"""

    def __init__(
        self,
        max_steps=8,
        start_fragment="Start",
        exploration_weight=1.41421356237,
        rollout_limit=32,
        reward_fn=None,
        random_seed=None,
    ):
        self.max_steps = max_steps
        self.start_fragment = start_fragment
        self.exploration_weight = exploration_weight
        self.rollout_limit = rollout_limit
        self.reward_fn = reward_fn or default_reward
        self.reward_cache = {}
        if random_seed is not None:
            random.seed(random_seed)

    def search(self, iterations=200):
        """运行 MCTS，并返回一条访问次数较稳定的片段序列。"""
        root = MCTSNode(
            PolymerState(sequence=(self.start_fragment,), max_steps=self.max_steps)
        )

        for _ in range(iterations):
            node = self._select(root)
            reward = self._simulate(node.state)
            node.backpropagate(reward)

        best_child = root.most_visited_child()
        if best_child is None:
            return root.state.to_mcts_output()
        return self._best_terminal_sequence(best_child)

    def ranked_candidates(self, iterations=200, top_k=5):
        """返回若干条候选序列，便于后续批量建模和筛选。"""
        root = MCTSNode(
            PolymerState(sequence=(self.start_fragment,), max_steps=self.max_steps)
        )

        for _ in range(iterations):
            node = self._select(root)
            reward = self._simulate(node.state)
            node.backpropagate(reward)

        candidates = []
        self._collect_candidates(root, candidates)
        candidates.sort(key=lambda item: item["score"], reverse=True)
        return candidates[:top_k]

    def _select(self, node):
        """从根节点向下走到可扩展节点或终止节点。"""
        while not node.is_terminal():
            if not node.is_fully_expanded():
                return node.expand()
            node = node.best_child(self.exploration_weight)
        return node

    def _simulate(self, state):
        """从当前状态随机 rollout 到终止状态，并计算分数。"""
        rollout_state = state
        for _ in range(self.rollout_limit):
            if rollout_state.is_terminal():
                break
            actions = rollout_state.legal_actions()
            if not actions:
                break
            rollout_state = rollout_state.take_action(random.choice(actions))
        return self._evaluate_sequence(rollout_state.to_mcts_output())

    def _best_terminal_sequence(self, node):
        """从给定节点继续选择访问次数最多的路径。"""
        current = node
        while current.children:
            current = current.most_visited_child()
        return current.state.to_mcts_output()

    def _collect_candidates(self, node, candidates):
        """遍历树，把已经访问过的片段路径整理成候选结果。"""
        if node.visits > 0 and node.state.current_fragment != "Start":
            sequence = node.state.to_mcts_output()
            candidates.append(
                {
                    "sequence": sequence,
                    "score": node.average_reward,
                    "visits": node.visits,
                    "average_reward": node.average_reward,
                }
            )

        for child in node.children:
            self._collect_candidates(child, candidates)

    def _evaluate_sequence(self, sequence):
        """计算序列分数，并缓存已经评估过的序列。"""
        key = tuple(sequence)
        if key not in self.reward_cache:
            self.reward_cache[key] = self.reward_fn(sequence)
        return self.reward_cache[key]


if __name__ == "__main__":
    engine = MCTSEngine(max_steps=6, random_seed=7)
    print(engine.search(iterations=100))
    for candidate in engine.ranked_candidates(iterations=100, top_k=3):
        print(candidate)
