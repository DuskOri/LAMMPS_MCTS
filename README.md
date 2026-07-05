# PolymerMCTS-MD：基于蒙特卡洛树搜索的高分子热导率搜索框架

[![Python](https://img.shields.io/badge/python-3.11-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## 项目简介
本项目实现了一个**自主探索高热导率线性非晶聚合物**的计算框架。核心思想是将**蒙特卡洛树搜索（MCTS）** 与**分子动力学（MD）** 模拟结合。

代码思路基于 *Macromolecules* 期刊 2022 年的前沿工作 (Nagoya et al.)。该项目利用 MCTS 在一个预设的有向图（化学合成规则）中搜索最优的化学片段组合。

本项目使用 MIT License 协议。

## 核心功能
- **MCTS 自主搜索算法**：自定义实现了 UCB 公式，在自定义的化学空间有向图中实现搜索。
- **化学片段翻译器**：将 MCTS 输出的抽象标签（如 `['A', 'B', 'C', 'End']`）自动翻译为真实的 RDKit SMILES 化学结构。
- **自动化聚合物生成**：
  - 自动将单个单体骨架**聚合 10 次（DP=10）**，生成一条满足 MD 计算要求的高分子长链。
  - 自动生成 3D 空间构象，并导出为 **`.pdb`** 格式文件。
- **LAMMPS 自动化热导率计算**：
  - **NPT 压缩**、**NVT 平衡**、**NVE 系综下 GK 模拟**的自动提交流程。
  - 集成 **“双热流估算方法”**，单条轨迹即可输出低统计误差的热导率\(\kappa\)。

## 项目目录结构
```text
LAMMPS_MCTS/
├── main.py                      # 主程序入口
├── config.yaml                  # 存放 MD 模拟参数（温度、步数、压力等）
├── mcts_core/                   # MCTS 算法核心
│   ├── fragment_graph.py        # 有向图规则定义，[待完善]前端自定义有向图
│   ├── state.py                 # 聚合物拼接状态
│   ├── node.py                  # MCTS 树节点及 UCB 公式
│   └── mcts_engine.py           # MCTS 四大步骤主循环
├── generator/                   # 聚合物结构生成器
│   ├── polymer_builder.py       # 基于 SMILES 聚合为长链
│   └── packmol_runner.py        # 调用 Packmol 生成多链无定形盒子
├── md_engine/                   # MD 计算引擎
│   ├── lammps_runner.py         # 自动生成 input 脚本并提交 LAMMPS 作业
│   └── gk_analyzer.py           # 解析 LAMMPS 日志，计算快速 GK 热导率
├── post_process/                # [待完善] 后处理模块
│   ├── feature_extractor.py     # [待完善] 提取 SMILES 指纹
│   └── explain_model.py         # [待完善] 训练 GBDT 模型与 SHAP 分析
├── templates/                   # LAMMPS 输入脚本模板
│   ├── minimize.in
│   ├── compress.in
│   └── gk_eval.in
└── requirements.txt             # Python 依赖文件