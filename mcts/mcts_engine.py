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
        max_steps=5,
        start_fragment="Start",
        exploration_weight=1.41421356237,
        rollout_limit=32,
        reward_fn=None,
        random_seed=None,
        progress_callback=None,
    ):
        self.max_steps = max_steps
        self.start_fragment = start_fragment
        self.exploration_weight = exploration_weight
        self.rollout_limit = rollout_limit
        self.reward_fn = reward_fn or default_reward
        self.reward_cache = {}
        self.evaluated_candidates = {}
        self.progress_callback = progress_callback
        self.current_iteration = 0
        self.total_iterations = 0
        if random_seed is not None:
            random.seed(random_seed)

    def ranked_candidates(self, iterations=200, top_k=5):
        """返回若干条候选序列，便于后续批量建模和筛选。"""
        self.evaluated_candidates = {}
        root = MCTSNode(
            PolymerState(sequence=(self.start_fragment,), max_steps=self.max_steps)
        )

        self.total_iterations = int(iterations)
        for iteration in range(1, self.total_iterations + 1):
            self.current_iteration = iteration
            self._notify_iteration()
            node = self._select(root)
            reward = self._simulate(node.state)
            node.backpropagate(reward)

        candidates = self._rank_evaluated_candidates()
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

    def _evaluate_sequence(self, sequence):
        """计算序列分数，并缓存已经评估过的序列。"""
        key = tuple(sequence)
        if key not in self.reward_cache:
            self.reward_cache[key] = self.reward_fn(sequence)
        self._record_evaluated_sequence(key, self.reward_cache[key])
        return self.reward_cache[key]

    def _record_evaluated_sequence(self, key, reward):
        """记录 rollout 中真正计算过 reward 的完整序列。"""
        item = self.evaluated_candidates.get(key)
        if item is None:
            item = {
                "sequence": list(key),
                "score": reward,
                "visits": 1,
                "average_reward": reward,
                "total_reward": reward,
            }
            self.evaluated_candidates[key] = item
        else:
            item["visits"] += 1
            item["total_reward"] += reward
            item["average_reward"] = item["total_reward"] / item["visits"]
            item["score"] = item["average_reward"]

        self._notify(
            "candidate",
            f"candidate evaluated at iteration {self.current_iteration}",
            {
                "iteration": self.current_iteration,
                "total_iterations": self.total_iterations,
                "candidate": {
                    "sequence": list(item["sequence"]),
                    "score": item["score"],
                    "visits": item["visits"],
                    "average_reward": item["average_reward"],
                },
            },
        )

    def _notify_iteration(self):
        """向界面报告当前 MCTS 循环次数。"""
        self._notify(
            "search",
            f"MCTS 第 {self.current_iteration}/{self.total_iterations} 轮：正在选择并扩展节点",
            {
                "iteration": self.current_iteration,
                "total_iterations": self.total_iterations,
            },
        )

    def _notify(self, event, message, details=None):
        """发送搜索进度，同时兼容旧的双参数回调。"""
        if self.progress_callback is None:
            return
        try:
            self.progress_callback(event, message, details or {})
        except TypeError:
            self.progress_callback(event, message)

    def _rank_evaluated_candidates(self):
        """把搜索阶段实际评价过的序列整理成候选列表。"""
        return [
            {
                "sequence": list(item["sequence"]),
                "score": item["score"],
                "visits": item["visits"],
                "average_reward": item["average_reward"],
            }
            for item in self.evaluated_candidates.values()
        ]


if __name__ == "__main__":
    engine = MCTSEngine(max_steps=5, random_seed=7)
    for candidate in engine.ranked_candidates(iterations=100, top_k=3):
        print(candidate)
