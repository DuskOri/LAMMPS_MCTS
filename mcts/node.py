"""MCTS 树节点及 UCB 选择公式。"""

import math
import random


class MCTSNode:
    """树搜索中的单个节点。"""

    def __init__(self, state, parent=None, action=None):
        self.state = state
        self.parent = parent
        self.action = action
        self.children = []
        self.untried_actions = state.legal_actions()
        self.visits = 0
        self.total_reward = 0.0

    @property
    def average_reward(self):
        """节点平均收益，未访问时按 0 处理。"""
        if self.visits == 0:
            return 0.0
        return self.total_reward / self.visits

    def is_fully_expanded(self):
        """所有候选动作都扩展过，则该节点已经完全展开。"""
        return len(self.untried_actions) == 0

    def is_terminal(self):
        """判断节点对应状态是否已经结束。"""
        return self.state.is_terminal()

    def expand(self):
        """随机取一个未尝试动作，并生成子节点。"""
        action = self.untried_actions.pop(random.randrange(len(self.untried_actions)))
        child = MCTSNode(
            state=self.state.take_action(action),
            parent=self,
            action=action,
        )
        self.children.append(child)
        return child

    def best_child(self, exploration_weight=1.41421356237):
        """按 UCB 分数选择子节点。"""
        if not self.children:
            raise ValueError("Cannot select from a node without children.")

        parent_visits = max(self.visits, 1)
        return max(
            self.children,
            key=lambda child: child.ucb_score(parent_visits, exploration_weight),
        )

    def ucb_score(self, parent_visits, exploration_weight=1.41421356237):
        """UCB = 平均奖励项 + 探索项。"""
        if self.visits == 0:
            return float("inf")

        reward_term = self.average_reward
        exploration_term = exploration_weight * math.sqrt(
            math.log(parent_visits) / self.visits
        )
        return reward_term + exploration_term

    def backpropagate(self, reward):
        """把一次模拟得到的分数沿父节点方向回传。"""
        self.visits += 1
        self.total_reward += reward
        if self.parent is not None:
            self.parent.backpropagate(reward)
