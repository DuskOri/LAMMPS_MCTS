# LAMMPS_MCTS 开发记录

本文档用于记录 LAMMPS_MCTS 项目的阶段性开发过程，便于后续整理软著材料、项目说明和答辩文档。记录内容按开发时间和功能推进顺序编写，重点描述每个阶段完成的工作、形成的文件和验证结果。

## 2026-07-03 初始聚合物片段构建原型

### 目标

建立项目最早的聚合物片段表示方式，使 MCTS 输出的片段序列能够被转换为 RDKit 可处理的聚合物分子。

### 已完成内容

- 建立 `mcts/Poly_Build.py` 初始版本。
- 定义聚酰亚胺相关片段 SMILES。
- 使用 `[1*]` 和 `[2*]` 表示片段左右连接位点。
- 实现基础片段拼接逻辑。
- 实现聚合物链端基封端逻辑。
- 实现 PDB 文件输出的初步流程。

### 阶段结果

- 项目具备从片段序列生成聚合物链的基础能力。
- 后续 MCTS 搜索结果可以对接该建链模块。

## 2026-07-03 初始片段连接图

### 目标

定义最早的聚合物片段搜索空间，为 MCTS 搜索提供片段连接规则。

### 已完成内容

- 建立 `mcts/graph.py` 初始版本。
- 定义 `Start` 和 `End` 节点。
- 定义聚酰亚胺片段 `A`、`A'`、`B`、`C`。
- 定义片段之间的基础连接关系。

### 阶段结果

- 项目具备最早的聚酰亚胺片段连接图。
- MCTS 后续可以基于该图进行片段序列搜索。

## 2026-07-07 项目接手、梳理与基础修复

### 目标

对已有 LAMMPS_MCTS 项目进行代码检查，在 2026-07-03 已有 `Poly_Build.py` 和 `graph.py` 的基础上，修复影响运行的基础问题，并明确项目后续开发方向。

### 工作内容

- 检查项目目录结构和核心文件。
- 阅读已有的 `mcts/Poly_Build.py` 和 `mcts/graph.py`。
- 修复 `mcts/Poly_Build.py` 中 PDB 保存相关问题。
- 修正 RDKit 三维构象生成流程，增加构象嵌入失败后的随机坐标重试。
- 整理 `requirements.txt`，写入项目运行所需 Python 依赖。

### 阶段结果

- RDKit 建链和 PDB 保存流程可以独立运行。
- 项目具备继续扩展 MCTS 和分子模拟模块的基础。

## 2026-07-07 MCTS 搜索模块完善

### 目标

补全 MCTS 搜索流程，使项目可以从片段连接图出发自动搜索候选聚合物序列。

### 工作内容

- 在已有 `graph.py` 连接图基础上，完善 MCTS 节点、状态、搜索引擎和片段图调用。
- 实现 UCB 选择逻辑。
- 明确 UCB 公式为奖励项与探索项之和。
- 设置搜索根节点为 `Start`。
- 增加 `max_steps`、`iterations`、`top_k`、`exploration_weight` 等配置项。
- 将候选序列按得分排序输出。

### 阶段结果

- MCTS 可以从 `Start` 节点开始搜索。
- 搜索结果可写入 `outputs/candidates.csv`。
- 候选序列可继续传入 RDKit 建链模块。

## 2026-07-07 低热导目标与反馈逻辑

### 目标

将热导率结果作为 MCTS 的反馈，使搜索方向面向低热导聚合物。

### 工作内容

- 明确项目目标为寻找低热导聚合物。
- 将奖励函数方向调整为低热导优先。
- 设计热导反馈奖励：

```text
reward = 1 / (kappa + offset)
```

- 增加启发式低热导评分函数。
- 增加热导反馈评估模块 `mcts/thermal_feedback.py`。
- 将候选数据库按热导率从低到高排序。

### 阶段结果

- MCTS 奖励越高代表热导率越低。
- `low_k_database.csv` 用于保存低热导候选数据库。
- 热导率结果可以反馈到 MCTS 节点。

## 2026-07-07 RDKit 聚合物构建模块扩展

### 目标

在 2026-07-03 已有 `Poly_Build.py` 的基础上，把 MCTS 输出的片段序列稳定转化为可用于后续模拟的聚合物结构。

### 工作内容

- 保留并修复已有 `mcts/Poly_Build.py` 的片段拼接逻辑。
- 完善 `generator/polymer_builder.py`。
- 将 MCTS 片段序列标准化，自动补齐 `End`。
- 将聚合度参数接入构建流程。
- 设置 `degree_of_polymerization`，后续调整为 DP=10。
- 输出候选聚合物的 SMILES 和 PDB 文件。

### 阶段结果

- 候选序列可以转化为聚合物链。
- PDB 文件输出到 `outputs/pdb/`。
- DP=10 时，序列片段会重复 10 次构建长链。

## 2026-07-07 PDB 到 LAMMPS Data 转换

### 目标

实现从 PDB 文件到 LAMMPS data 文件的自动转换，为 LAMMPS 模拟准备输入结构。

### 工作内容

- 新增 `md_engine/lammps_data_writer.py`。
- 解析 PDB 中的原子坐标和连接信息。
- 根据元素类型生成 atom types。
- 根据键连接生成 bond types。
- 输出 `atom_style full` 风格的 LAMMPS data 文件。

### 阶段结果

- 单链 PDB 可以转换为 `candidate_xxx_single.data`。
- Packmol 输出的体系 PDB 可以转换为 `candidate_xxx_system.data`。

## 2026-07-07 Packmol 初始体系生成

### 目标

自动生成多条聚合物链组成的初始模拟体系。

### 工作内容

- 新增 `generator/packmol_runner.py`。
- 编写 Packmol 输入文件。
- 配置 Packmol 可执行文件路径。
- 调用 Packmol 进行多链装盒。
- 将 Packmol 输出体系转换为 LAMMPS data。

### 阶段结果

- `candidate_001` 和 `candidate_002` 成功生成多链初始体系。
- 输出文件包括 Packmol 输入、体系 PDB 和体系 data。
- 对较大候选结构，观察到 Packmol 收敛时间明显增加，后续可通过增大盒子或降低分子数优化。

## 2026-07-07 Packmol 安装与环境配置

### 目标

解决 Windows 环境中 Packmol 不可用的问题。

### 工作内容

- 尝试通过 conda 安装 Packmol。
- 由于 win-64 依赖源不可用，改为源码编译。
- 获取 Packmol 源码。
- 使用 MinGW gfortran 编译生成 `packmol.exe`。
- 将 `packmol.exe` 放入 `polymer_mcts` 环境的 `Scripts` 目录。
- 在 `config.yaml` 中写入 Packmol 路径。

### 阶段结果

- Packmol 可以被主程序调用。
- 项目可以在当前 Windows 机器上自动生成初始体系。

## 2026-07-07 快速热导率 LAMMPS 模块

### 目标

在给定结构基础上生成快速热导率计算输入脚本。

### 工作内容

- 参考用户提供的热导率计算示例。
- 新增 `md_engine/thermal_conductivity_writer.py`。
- 生成 NEMD 风格快速热导率 LAMMPS 输入。
- 设置 hot/cold 区域。
- 输出热流、温度剖面和轨迹文件。
- 默认忽略分子间作用力，使用 `pair_style zero`。
- 增加 `lammps_runner.py` 用于自动运行 LAMMPS。

### 阶段结果

- 成功生成 `candidate_001_fast_tc.in` 和 `candidate_002_fast_tc.in`。
- 手动运行 `candidate_001_fast_tc.in` 成功完成 LAMMPS 热导计算。
- 生成 `flux.profile`、`temp.profile` 和 `dump.lammpstrj`。

## 2026-07-07 热导率后处理与数据库输出

### 目标

解析 LAMMPS 输出结果，形成可排序的低热导候选数据库。

### 工作内容

- 新增 `post_process/thermal_analyzer.py`。
- 读取热流输出文件。
- 读取温度剖面文件。
- 计算热流平均值和温度梯度。
- 估算热导率。
- 新增 `post_process/result_writer.py`，统一写出 CSV。

### 阶段结果

- 输出 `thermal_results.csv`。
- 输出 `low_k_database.csv`。
- 候选按低热导优先排序。

## 2026-07-07 多聚合物片段库扩展

### 目标

使项目不局限于聚酰亚胺，支持更多聚合物体系。

### 工作内容

- 新增 `mcts/fragment_registry.py`。
- 将片段 SMILES、片段家族、片段角色和连接规则集中管理。
- 保留聚酰亚胺片段。
- 增加聚氨酯片段。
- 增加酚醛树脂片段。
- 增加自定义片段入口。
- 修改 MCTS 图读取方式，使搜索空间来自片段库。
- 修改 RDKit 建链模块，使 SMILES 也来自片段库。

### 阶段结果

- 支持 `polyimide`、`polyurethane`、`phenolic` 和 `custom`。
- MCTS 可以在多个聚合物家族之间搜索。
- `Start` 节点可连接多个聚合物体系的起始片段。

## 2026-07-07 自定义片段输入

### 目标

允许用户输入自己的化合物片段，扩展搜索空间。

### 工作内容

- 新增 `custom_fragments.json`。
- 定义自定义片段 JSON 格式。
- 要求片段 SMILES 包含 `[1*]` 和 `[2*]`。
- 支持自定义 transitions 连接图。
- 在主程序中读取自定义片段并合并进片段库。

### 阶段结果

- 用户可通过 JSON 文件添加自己的片段。
- 自定义片段可参与 MCTS 搜索和 RDKit 建链。

## 2026-07-07 AI API 片段生成接口

### 目标

增加人工智能接口，使 AI 可以自动生成片段和连接图。

### 工作内容

- 新增 `mcts/ai_fragment_api.py`。
- 实现 OpenAI-compatible Chat Completions 接口。
- 支持通过环境变量读取 API key。
- 设计 AI 返回 JSON 格式。
- 校验 AI 返回的片段是否包含 `[1*]` 和 `[2*]`。
- 如果 RDKit 可用，进一步校验 SMILES 能否解析。
- 将 AI 输出保存为 `ai_fragments.json`。
- 增加 `ai_fragments.example.json` 示例文件。
- 增加 `AI_FRAGMENT_API.md` 说明文档。
- 在 `config.yaml` 中加入 AI 接口配置。

### 阶段结果

- AI 片段接口默认关闭，不影响离线运行。
- 打开 `ai_enabled` 后，AI 生成片段可以自动进入 MCTS 搜索空间。
- 本地模拟 AI 文件合并测试通过。

## 2026-07-07 README 与说明文档整理

### 目标

完善项目文档，形成可阅读、可交付的说明材料。

### 工作内容

- 重写 `README.md`。
- 增加项目功能说明。
- 增加安装和启动方式。
- 增加 MCTS 低热导奖励说明。
- 增加自定义片段说明。
- 增加 AI 接口说明。
- 增加 Packmol、LAMMPS 和输出文件说明。

### 阶段结果

- README 可作为项目使用说明。
- AI 接口有独立文档。
- 项目当前代码和文档总行数超过 3000 行。

## 2026-07-08 开发留痕整理

### 目标

补充项目开发周期记录，使项目过程更清晰。

### 工作内容

- 整理从基础修复到 AI 接口的阶段性开发过程。
- 按日期记录每个阶段的目标、工作内容和阶段结果。
- 形成 `DEVELOPMENT_LOG.md`。
- 准备同步整理 `CHANGELOG.md`。

### 阶段结果

- 项目具备按时间展开的开发记录。
- 后续可用于软著材料、项目说明和答辩准备。

## 验证记录

开发过程中执行过以下验证：

```text
syntax ok 21
syntax ok 22
syntax ok 23
```

阶段性运行结果包括：

- RDKit 聚合物构建通过。
- Packmol 对部分候选体系装盒通过。
- LAMMPS 快速热导率输入生成通过。
- `candidate_001_fast_tc.in` LAMMPS 运行通过。
- AI 示例片段合并入 MCTS 片段库通过。

## 2026-07-08 轻量前端控制台

### 目标

为项目增加一个简单的可视化入口，使候选序列、热导率结果、片段配置和常用命令可以在浏览器中查看。

### 工作内容

- 新建 `frontend/` 目录。
- 新增 `frontend/index.html`。
- 新增 `frontend/styles.css`。
- 新增 `frontend/app.js`。
- 新增 `web_server.py`，提供本地前端 API。
- 前端读取 `config.yaml`、`outputs/candidates.csv`、`outputs/thermal_results.csv` 等文件。
- 前端支持修改 DP、迭代次数、候选数量和片段空间。
- 前端支持一键启动主流程。
- 前端支持轮询显示当前运行阶段。
- 增加候选序列筛选功能。
- 增加配置摘要展示。
- 增加自定义片段和 AI 片段示例展示。
- 增加主流程、前端服务器和 LAMMPS 命令复制功能。
- 增加有向图配置编辑器。
- 增加 API 生成化合物配置功能。
- 支持将前端配置导出为 `custom_fragments.json` 或 `ai_fragments.json`。
- 增加 `frontend/assets/` 原创 SVG 素材。
- 设计项目标识、聚合物流程图和功能图标。
- 优化首页运行区、阶段区和数据库入口的视觉层次。

### 阶段结果

- 项目具备新手向前端控制台。
- 通过本地服务可以访问 `http://localhost:8010/frontend/`。
- 前端可以通过本地 API 修改常用配置并启动主流程。
- 前端可以人工配置片段连接有向图。
- 前端可以调用 OpenAI-compatible API 生成化合物片段配置。
- 前端视觉风格调整为更偏科研软件的浅色控制台，窄屏下指标卡和导航更紧凑。
- 前端未引入第三方素材，Logo 和图标均为项目内自绘 SVG。

## 2026-07-08 README 使用说明补充

### 目标

补充更完整的使用说明，使初次使用者可以根据 README 独立完成前端启动、参数配置、MCTS 搜索、片段图编辑、AI 片段生成和热导率结果查看。

### 工作内容

- 扩展 `README.md` 中的详细使用说明。
- 增加前端控制台的分步骤使用流程。
- 增加命令行运行流程。
- 增加 DP、MCTS 迭代次数、最大片段步数和反馈模式的调参说明。
- 增加自定义有向图片段输入说明。
- 增加 AI API 片段生成操作说明。
- 增加热导率计算和结果文件查看说明。
- 增加常见问题排查说明。
- 修正 README 目录结构中 `frontend/` 与 `web_server.py` 的展示位置。

### 阶段结果

- README 从功能简介扩展为较完整的用户手册。
- 项目留痕文档记录了本次文档增强工作。
- 后续整理软著材料时，可以直接引用 README 的使用流程和开发记录。

## 2026-07-08 前端视觉风格二次调整

### 目标

参考用户提供的粉紫色侧边栏页面，将前端控制台从普通表单页面调整为更像小型科研工具的页面，同时保留新手向操作入口。

### 工作内容

- 调整 `frontend/index.html` 中的导航分组和首页标题文案。
- 重写 `frontend/styles.css` 的页面主题。
- 将顶部导航在宽屏下调整为左侧导航，在窄屏下自动变为顶部导航。
- 使用粉紫色背景、蓝色主按钮和分组面板增强页面识别度。
- 保留原有 Logo、聚合物流程图和功能图标。
- 收紧页面横向宽度，修复窄屏下的横向滚动问题。
- 验证页面 SVG 资源均可正常加载。
- 验证浏览器控制台无报错。

### 阶段结果

- 前端更接近参考页面的粉紫色工具风格。
- 页面仍保留 DP、迭代次数、有向图、API 生成和结果数据库等原有功能。
- 移动/窄屏显示不再出现页面级横向滚动。

## 2026-07-08 有向图编辑器可视化

### 目标

将前端有向图配置从纯表单和 JSON 输出改为更直观的图形预览，使使用者能直接看到 `Start`、片段节点、`End` 和片段之间的可选连接方向。

### 工作内容

- 在 `frontend/index.html` 的有向图编辑区增加图形预览面板。
- 增加“载入示例”按钮，方便新用户快速理解有向图配置方式。
- 在 `frontend/app.js` 中增加 SVG 绘图逻辑。
- 自动绘制 `Start -> 片段` 的入口边。
- 自动绘制用户配置的片段跳转边和结束边。
- 增加连接边标签列表，便于对照图形和规则。
- 将 JSON 配置折叠到高级查看区域。
- 调整 `frontend/styles.css`，增加有向图画布、节点、箭头和边列表样式。

### 阶段结果

- 有向图页面可以直接显示片段连接图。
- 示例图可生成 5 个节点和多条连接边。
- 页面仍保持原有保存自定义图和 API 片段生成能力。

## 当前行数统计

最近一次统计结果：

```text
总行数：6585
Python：3426
Markdown：1348
HTML：264
CSS：747
JavaScript：629
SVG：51
JSON：44
YAML：73
TXT：3
```
