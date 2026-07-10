# Changelog

本文档记录 LAMMPS_MCTS 项目的主要功能变更。

## 0.6.3 - 2026-07-08

### Added

- 有向图编辑器新增 SVG 图形预览。
- 新增“载入示例”按钮。
- 新增连接边标签列表。
- 新增 JSON 配置折叠查看区域。

### Changed

- 有向图页面从纯 JSON 输出调整为图形优先显示。
- 优化自定义片段图的节点、箭头和画布样式。

## 0.6.2 - 2026-07-08

### Changed

- 参考用户提供的页面风格，二次调整前端视觉。
- 将前端主题改为粉紫色工具界面。
- 优化导航为宽屏侧边栏、窄屏顶部导航。
- 调整首页标题、说明文字、面板、按钮、表格和状态区样式。
- 修复窄屏下页面级横向滚动问题。

## 0.6.1 - 2026-07-08

### Changed

- 扩展 README 的详细使用说明。
- 增加前端控制台、命令行运行、参数配置、自定义有向图、AI 片段生成和热导率计算的操作说明。
- 增加常见问题排查内容。
- 更新开发留痕文档，记录 README 使用说明补充工作。
- 修正 README 目录结构展示。

## 0.6.0 - 2026-07-08

### Added

- 新增轻量静态前端 `frontend/`。
- 新增本地前端接口服务 `web_server.py`。
- 新增新手向运行页面。
- 新增 DP、迭代次数、候选数量等参数修改入口。
- 新增一键运行入口和阶段状态显示。
- 新增 MCTS 候选序列表格。
- 新增热导率结果表格。
- 新增配置摘要和片段预览。
- 新增常用命令复制入口。
- 新增有向图配置编辑器。
- 新增前端 API 生成化合物片段配置功能。
- 支持复制或下载 `custom_fragments.json` 和 `ai_fragments.json`。
- 新增 `frontend/assets/` 原创 SVG 素材。
- 新增项目 Logo、聚合物流程图和功能图标。

### Changed

- 调整前端视觉风格为更克制的科研控制台界面。
- 优化窄屏布局，导航、指标卡和操作按钮更紧凑。
- 修正最低热导率显示逻辑，只统计有效正值结果。
- 优化首页运行区、阶段区和数据库入口的视觉层次。

## 0.5.0 - 2026-07-08

### Added

- 新增开发周期记录 `DEVELOPMENT_LOG.md`。
- 整理项目阶段性工作内容和验证记录。
- 记录项目从 MCTS 基础搜索到 AI 片段接口的开发过程。

## 0.4.0 - 2026-07-07

### Added

- 新增 AI 片段生成接口 `mcts/ai_fragment_api.py`。
- 支持 OpenAI-compatible Chat Completions API。
- 支持从环境变量读取 API key。
- 支持 AI 生成 `fragments` 和 `transitions`。
- 支持校验 AI 返回的双端连接 SMILES。
- 新增 `AI_FRAGMENT_API.md`。
- 新增 `ai_fragments.example.json`。

### Changed

- `main.py` 的 chemistry 配置流程接入 AI 生成片段。
- `config.yaml` 增加 AI 接口配置项。

## 0.3.0 - 2026-07-07

### Added

- 新增统一片段库 `mcts/fragment_registry.py`。
- 支持聚酰亚胺、聚氨酯、酚醛树脂和自定义片段。
- 新增 `custom_fragments.json`。
- 支持自定义片段和自定义连接图。

### Changed

- MCTS 连接图从片段库动态读取。
- RDKit 建链模块从片段库读取 SMILES。
- 低热导启发式评分扩展为基于片段角色计算。
- CSV 特征输出增加多聚合物体系相关字段。

## 0.2.0 - 2026-07-07

### Added

- 新增 PDB 到 LAMMPS data 转换模块。
- 新增 Packmol 初始体系生成模块。
- 新增快速热导率 LAMMPS 输入脚本生成模块。
- 新增 LAMMPS 运行封装模块。
- 新增热导率后处理模块。
- 新增低热导数据库输出。

### Changed

- 主流程扩展为 MCTS、RDKit、Packmol、LAMMPS 输入生成和结果后处理的完整流程。
- `config.yaml` 增加 system、thermal_conductivity、lammps 和 output 配置。

## 0.1.0 - 2026-07-07

### Added

- 在 2026-07-03 已有 `Poly_Build.py` 和 `graph.py` 的基础上，整理 MCTS 搜索流程。
- 实现 UCB 选择、扩展、rollout 和回传。
- 实现候选序列排序输出。
- 实现低热导启发式奖励。
- 接入并修复 RDKit 聚合物链构建和 PDB 输出。

### Fixed

- 修复 PDB 保存流程中的变量问题。
- 增加 RDKit 构象生成失败后的重试逻辑。
- 补充 `requirements.txt`。

## 0.0.1 - 2026-07-03

### Added

- 建立初始 `mcts/Poly_Build.py`。
- 建立初始 `mcts/graph.py`。
- 定义聚酰亚胺片段 `A`、`A'`、`B`、`C`。
- 定义 `Start`、`End` 和片段连接规则。
- 建立基于 RDKit 的片段拼接原型。
