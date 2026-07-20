# Changelog

## 2026-07-20

- 将 MCTS 路径定义为 Start 到 End 的重复单元，默认最多 5 个真实片段。
- 用快速 Green-Kubo 分子内热流自相关流程替换零力 NEMD 输入。
- 增加 LAMMPS data 力场系数检查，未参数化体系不再产生伪热导率 reward。
- 增加 RDKit UFF 筛选级参数化，可将 Packmol 坐标直接写成带非键和成键参数的 LAMMPS data。
- 将低热导 reward 归一化到 0 到 1，并默认关闭 thermal 失败后的启发式回落。
- 暂不增加家族官能团规则，候选重复单元由人工筛选。
- 重构前端为简洁的科研控制台，统一导航、表单、运行状态、片段编辑和数据库布局，并补充移动端适配。
- 取消手动 DP 配置，改为按重复单元轮廓长度自动选择聚合度，使聚合物链接近默认目标 `120 Å`。
- 前端增加 reward 来源、探索权重、Green-Kubo 档位和 LAMMPS 超时配置。
- 增加 MCTS 循环轮次、实时候选和实时 Top K 展示。
- 修正内置聚合物家族与自定义片段的加载边界。
- 增加快速筛选与标准复核两套 Green-Kubo 参数档位。
- 清理未使用的旧 MCTS 搜索入口、兼容属性和过期前端结果列。

## 2026-07-19

- 将 LAMMPS 热导率评价接入每一轮 MCTS reward 回传。
- `ranked_candidates()` 改为只返回搜索阶段实际评价过的完整候选。
- 删除搜索结束后的重复 LAMMPS 批处理和热导后处理分支。
- 将候选按低热导 reward 排序并写入反馈记录和低热导数据库。
- 联调 MCTS 轮次、热导评价次数、缓存和最终 Top K 导出流程。

## 2026-07-18

- 复查 `main.py` 完整调用链，统一 MCTS、结构生成、热导评价和结果导出的职责。
- 明确低热导目标的 reward 方向，并核对 UCB 奖励项与探索项的使用方式。
- 梳理 `feedback_max_evaluations`、`iterations` 和 `top_k` 的参数语义。
- 确定 thermal 模式采用逐轮评价与回传，不再采用搜索结束后统一评价的旧流程。

本文档记录 LAMMPS_MCTS 项目的主要功能变更。

## 0.6.6 - 2026-07-10

### Added

- 新增 `/api/clear-data` 本地接口。
- 前端数据库页新增候选序列复选框。
- 新增“清除选中”和“清空候选”按钮。
- 清除前自动备份 CSV 到 `outputs/archive/`。

### Changed

- 数据库页增加候选选择数量提示。
- 清除操作不删除结构文件和 LAMMPS 输入文件，只处理 CSV 表格。

## 0.6.5 - 2026-07-10

### Changed

- 将前端说明页中的运行命令改为跨设备通用写法。
- 移除前端命令中的开发机器绝对路径。
- 更新前端脚本版本号，避免缓存旧命令。
- 在 LAMMPS 示例命令中增加 `PATH` 或配置文件路径提示。

## 0.6.4 - 2026-07-10

### Changed

- 将 README 中的运行命令改为跨机器通用写法。
- 移除 README 中开发机器的固定绝对路径。
- 增加 Packmol 和 LAMMPS 可执行文件路径的通用配置说明。
- 将主程序、前端服务和 LAMMPS 手动运行命令统一为更适合 GitHub 的示例。

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
