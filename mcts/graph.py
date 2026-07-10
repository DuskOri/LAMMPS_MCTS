# mcts/graph.py
# 规定链段的连接规则

try:
    from .fragment_registry import get_fragment_registry
except ImportError:
    from fragment_registry import get_fragment_registry


PolyGRAPH = get_fragment_registry().graph()

# 化学逻辑：A代表酸酐，B代表芳环二胺，C代表柔性基团
# Start为开始节点，End为结束

# 不妨令开头为为酸酐，
# 在组装过程中，酸酐A只能与芳环B连接，芳环B可以继续与自身或柔性基团C互相链接
# C 可以连回 B，可以自连，也可以直接结束
# B 可以自连，可以连到 C，也可以直接结束
# 最终选择End代表直接结束组装过程。

# 虚拟环境为/conda activate polymer_mcts
