# LAMMPS_MCTS

LAMMPS_MCTS 是一个面向低热导聚合物筛选的自动化建模与模拟程序。项目将蒙特卡洛树搜索、RDKit 分子构建、Packmol 初始体系生成、LAMMPS 快速热导率计算和结果后处理串联起来，用于在较大的化学片段空间中搜索潜在低热导聚合物候选。

当前目标是寻找热导率较低的聚合物结构。MCTS 中的奖励值越高，表示对应候选的估计热导率或反馈热导率越低。

## 功能概述

项目当前支持以下功能：

- 从 `Start` 节点开始执行 MCTS 搜索。
- 使用 UCB 公式进行奖励项和探索项平衡。
- 将低热导目标写成 MCTS 奖励反馈。
- 支持聚酰亚胺、聚氨酯、酚醛树脂和自定义片段。
- 支持 AI API 自动生成片段和连接图。
- 使用 RDKit 将片段序列拼接为聚合物链。
- 支持设置聚合度 `degree_of_polymerization`。
- 输出聚合物 PDB 文件。
- 将 PDB 转换为 LAMMPS data 文件。
- 调用 Packmol 生成多链初始体系。
- 生成快速热导率 LAMMPS 输入脚本。
- 可选择自动运行 LAMMPS。
- 解析热流和温度剖面，估算热导率。
- 建立按低热导优先排序的候选数据库。

## 目录结构

```text
LAMMPS_MCTS/
├── main.py
├── config.yaml
├── requirements.txt
├── custom_fragments.json
├── ai_fragments.example.json
├── AI_FRAGMENT_API.md
├── generator/
│   ├── polymer_builder.py
│   └── packmol_runner.py
├── mcts/
│   ├── Poly_Build.py
│   ├── ai_fragment_api.py
│   ├── fragment_registry.py
│   ├── fragment_graph.py
│   ├── graph.py
│   ├── mcts_engine.py
│   ├── node.py
│   ├── state.py
│   └── thermal_feedback.py
├── md_engine/
│   ├── lammps_data_writer.py
│   ├── lammps_runner.py
│   └── thermal_conductivity_writer.py
├── post_process/
│   ├── explain_model.py
│   ├── feature_extractor.py
│   ├── result_writer.py
│   └── thermal_analyzer.py
├── frontend/
│   ├── index.html
│   ├── styles.css
│   ├── app.js
│   └── assets/
└── web_server.py
```

## 环境准备

推荐使用已经配置好的 conda 环境：

```powershell
conda activate polymer_mcts
```

如果当前终端不能正常激活 conda，也可以直接使用环境里的 Python：

```powershell
C:\DuskORI\Application\miniconda\envs\polymer_mcts\python.exe main.py
```

Python 依赖写在 `requirements.txt`：

```text
rdkit
pyyaml
```

外部程序：

- Packmol：用于生成多链初始体系。
- LAMMPS：用于运行快速热导率计算。

当前配置中的路径为：

```yaml
system:
  packmol_executable: C:/DuskORI/Application/miniconda/envs/polymer_mcts/Scripts/packmol.exe

lammps:
  executable: C:/DuskORI/Application/LAMMPS/bin/lmp.exe
```

## 快速启动

进入项目目录：

```powershell
cd C:\DuskORI\Files\Code\LAMMPS_MCTS\LAMMPS_MCTS
```

运行主程序：

```powershell
C:\DuskORI\Application\miniconda\envs\polymer_mcts\python.exe main.py
```

默认流程会执行：

1. 读取 `config.yaml`。
2. 加载内置片段和自定义片段。
3. 从 `Start` 节点开始 MCTS 搜索。
4. 生成候选聚合物片段序列。
5. 按 DP 构建聚合物链。
6. 输出 PDB。
7. 调用 Packmol 生成初始体系。
8. 写出快速热导率 LAMMPS 输入文件。
9. 如果配置允许，则运行 LAMMPS。
10. 解析热导率结果并写入 CSV。

## 前端控制台

项目提供一个轻量静态前端，目录为：

```text
frontend/
```

启动本地前端：

```powershell
cd C:\DuskORI\Files\Code\LAMMPS_MCTS\LAMMPS_MCTS
C:\DuskORI\Application\miniconda\envs\polymer_mcts\python.exe web_server.py
```

浏览器打开：

```text
http://localhost:8010/frontend/
```

前端当前提供：

- 新手向运行页。
- DP、迭代次数、候选数量、最大片段步数修改。
- 已有有向图/片段空间选择。
- 自定义有向图创建。
- 一键运行完整流程。
- 当前运行阶段显示。
- MCTS 候选序列和低热导数据库查看。
- 热导率结果表。
- API 生成化合物片段配置。
- 主程序、前端服务器和 LAMMPS 命令复制。
- 项目内原创 Logo、聚合物流程图和功能图标。

新版前端通过 `web_server.py` 提供本地 API，可以修改 `config.yaml`、保存自定义有向图，并在后台启动主流程。该服务只绑定本机使用，API key 不会写入项目文件。

有向图编辑器可以人工输入片段 key、SMILES、family、roles，并添加有向边。生成结果可复制或下载为 `custom_fragments.json`。API 生成器可以调用 OpenAI-compatible 接口生成片段和连接图，输出可下载为 `ai_fragments.json`。API key 只保存在浏览器当前页面，不写入项目文件。

`frontend/assets/` 中的 SVG 素材为项目内自绘素材，没有引入第三方图片文件，后续整理软著或商用说明时更容易说明来源。

## 详细使用说明

本项目有两种使用方式：一种是通过浏览器前端操作，适合初次使用和演示；另一种是直接运行 `main.py`，适合批量计算和命令行调试。

### 方式一：使用前端控制台

1. 打开 PowerShell，进入项目目录：

```powershell
cd C:\DuskORI\Files\Code\LAMMPS_MCTS\LAMMPS_MCTS
```

2. 启动本地服务：

```powershell
C:\DuskORI\Application\miniconda\envs\polymer_mcts\python.exe web_server.py
```

3. 在浏览器访问：

```text
http://localhost:8010/frontend/
```

4. 在“运行”页面设置参数：

- `聚合度 DP`：控制候选聚合物片段序列重复次数，当前默认 DP=10。
- `MCTS 迭代次数`：控制搜索次数，数值越大搜索越充分，耗时也越长。
- `输出候选数`：控制最终保留多少个候选结构。
- `最大片段步数`：控制一次候选序列中片段节点的最大数量。
- `选择已有有向图/片段空间`：可以选择全部内置、聚酰亚胺、聚氨酯、酚醛树脂或自定义片段。
- `自动运行 LAMMPS`：勾选后会在生成输入文件后尝试直接调用 LAMMPS。
- `热导结果反馈 MCTS`：勾选后使用热导率结果作为搜索反馈，计算成本更高。

5. 点击“保存参数”，前端会把常用参数写入 `config.yaml`。

6. 点击“一键运行”，程序会按阶段执行：

```text
MCTS 搜索 -> RDKit 建链 -> Packmol 初始体系 -> 写出 LAMMPS 输入 -> 运行 LAMMPS -> 热导率后处理 -> 保存数据库
```

7. 运行结束后，在“数据库”页面查看候选序列和热导率结果。低热导数据库会优先展示热导率较低的候选。

### 方式二：命令行运行

命令行方式适合不打开前端时使用：

```powershell
cd C:\DuskORI\Files\Code\LAMMPS_MCTS\LAMMPS_MCTS
C:\DuskORI\Application\miniconda\envs\polymer_mcts\python.exe main.py
```

命令行运行会读取 `config.yaml`，然后自动完成搜索、建链、装盒、LAMMPS 输入生成和结果整理。是否自动运行 LAMMPS 由配置项控制：

```yaml
lammps:
  run_enabled: false
```

如果只想快速生成候选结构和 LAMMPS 输入文件，可以保持 `run_enabled: false`。如果本机 LAMMPS 路径已经配置正确，可以改为 `true`。

### 修改搜索参数

常用参数集中在 `config.yaml`：

```yaml
mcts:
  iterations: 300
  max_steps: 6
  top_k: 5
  feedback_mode: heuristic

polymer:
  degree_of_polymerization: 10
```

建议调参顺序：

1. 初次测试时使用较小的 `iterations`，例如 20 到 50。
2. 确认 RDKit、Packmol 和 LAMMPS 输入生成正常后，再提高到 300 或更高。
3. 如果候选序列过短，可以提高 `max_steps`。
4. 如果需要更长聚合物链，可以提高 `degree_of_polymerization`。
5. 如果只是快速筛选，优先使用 `feedback_mode: heuristic`。
6. 如果要把热导率计算结果反馈给 MCTS，再使用 `feedback_mode: thermal`。

### 选择和创建有向图

内置片段空间包括：

- `polyimide`：聚酰亚胺片段。
- `polyurethane`：聚氨酯片段。
- `phenolic`：酚醛树脂片段。
- `custom`：用户自定义片段。

前端“有向图”页面可以创建新的片段和连接边。每个片段至少需要：

- `key`：片段编号，例如 `CUSTOM_AR`。
- `SMILES`：带两个连接点的片段 SMILES，例如 `[1*]c1ccc([2*])cc1`。
- `family`：片段所属类别，例如 `custom`。
- `roles`：片段特征标签，例如 `rigid,aromatic`。

有向边表示某个片段后面可以接哪些片段。例如：

```json
{
  "CUSTOM_AR": ["CUSTOM_FLEX", "End"],
  "CUSTOM_FLEX": ["CUSTOM_AR", "End"]
}
```

保存后会写入 `custom_fragments.json`，下次运行时自动合并进搜索空间。

### 使用 AI API 生成片段

前端“有向图”页面提供 API 生成入口。使用时需要填写：

- API URL：OpenAI-compatible Chat Completions 地址。
- Model：模型名称。
- API Key：当前页面临时使用，不会写入项目文件。
- 片段数量：期望生成的片段数。
- 生成要求：对片段类型、结构倾向和连接规则的描述。

生成结果应包含 `fragments` 和 `transitions` 两部分。建议生成后先人工检查 SMILES 和连接关系，再下载为 `ai_fragments.json` 或复制到自定义片段文件中。

### 运行热导率计算

项目会先生成快速热导率 LAMMPS 输入文件，例如：

```text
outputs/lammps/candidate_001_fast_tc.in
```

手动运行某个候选：

```powershell
C:\DuskORI\Application\LAMMPS\bin\lmp.exe -in outputs\lammps\candidate_001_fast_tc.in
```

当前快速热导率脚本用于初筛和流程验证。脚本默认忽略分子间非键作用，适合快速比较候选，不适合直接作为最终物性结论。若要获得更可靠的热导率，需要补充真实力场、结构弛豫、平衡过程和多次重复采样。

### 查看结果

常用结果文件如下：

- `outputs/candidates.csv`：MCTS 搜索到的候选片段序列。
- `outputs/build_results.csv`：RDKit 建链结果和 PDB 路径。
- `outputs/system_results.csv`：Packmol 初始体系和 LAMMPS data 结果。
- `outputs/thermal_inputs.csv`：快速热导率 LAMMPS 输入脚本路径。
- `outputs/thermal_results.csv`：热流、温度梯度和热导率估算。
- `outputs/low_k_database.csv`：按低热导优先排序的候选数据库。

如果只关心最终低热导候选，优先查看 `outputs/low_k_database.csv`。如果要排查流程问题，则按上述文件顺序逐步检查。

### 常见问题

- 前端打不开：确认 `web_server.py` 是否正在运行，端口是否为 `8010`。
- RDKit 报错：确认当前环境为 `polymer_mcts`，并且已安装 RDKit。
- Packmol 失败：可以适当增大 `system.box_size` 或降低 `system.molecule_count`。
- LAMMPS 没有运行：确认 `lammps.executable` 路径正确，并将 `lammps.run_enabled` 改为 `true`。
- 没有热导率结果：先检查 `outputs/thermal_inputs.csv` 是否生成，再检查 LAMMPS 是否成功输出 `flux.profile` 和 `temp.profile`。
- AI 片段不可用：确认 API key、API URL、模型名称和返回 JSON 格式是否正确。

## 主要配置

核心配置在 `config.yaml`。

```yaml
mcts:
  iterations: 300
  max_steps: 6
  start_fragment: Start
  top_k: 5
  exploration_weight: 1.41421356237
  rollout_limit: 32
  random_seed: 7
  feedback_mode: heuristic

polymer:
  degree_of_polymerization: 10
  build_pdb: true
```

其中：

- `iterations`：MCTS 搜索迭代次数。
- `max_steps`：片段序列最大长度，不是聚合度。
- `start_fragment`：搜索起点，当前为 `Start`。
- `top_k`：输出候选数量。
- `degree_of_polymerization`：聚合度。比如序列 `A-B` 且 DP=10，表示 `A-B` 重复 10 次。

## MCTS 与低热导目标

MCTS 使用 UCB 公式：

```text
UCB = 平均奖励项 + 探索项
```

程序中的奖励值被设计为低热导优先。对于真实或估计热导率 `kappa`，反馈奖励可写成：

```text
reward = 1 / (kappa + offset)
```

因此：

- 热导率越低，reward 越高。
- MCTS 越倾向于保留低热导候选。
- 候选数据库按低热导优先排序。

当前支持两种反馈模式：

```yaml
mcts:
  feedback_mode: heuristic
```

`heuristic` 使用片段结构特征进行快速估计，适合大规模初筛。

```yaml
mcts:
  feedback_mode: thermal
  feedback_run_lammps: true
```

`thermal` 会在 MCTS 过程中把快速热导率计算结果反馈到节点。该模式计算成本较高，建议先用较小的 `feedback_max_evaluations` 测试。

## 聚合物片段库

当前内置片段家族包括：

- `polyimide`：聚酰亚胺。
- `polyurethane`：聚氨酯。
- `phenolic`：酚醛树脂。
- `custom`：用户自定义片段。

配置入口：

```yaml
chemistry:
  active_families: polyimide,polyurethane,phenolic
  custom_fragment_file: custom_fragments.json
  allow_custom_self_transitions: true
```

自定义片段写在 `custom_fragments.json`：

```json
{
  "fragments": [
    {
      "key": "CUSTOM_AR",
      "smiles": "[1*]c1ccc([2*])cc1",
      "family": "custom",
      "roles": ["rigid", "aromatic"]
    }
  ],
  "transitions": {
    "CUSTOM_AR": ["CUSTOM_AR", "End"]
  }
}
```

片段 SMILES 必须包含两个虚拟连接点：

- `[1*]`：左端连接位点。
- `[2*]`：右端连接位点。

程序拼接时会删除相邻片段的虚拟端点，并在真实原子之间建立单键。

## AI 片段接口

项目提供 OpenAI-compatible Chat Completions 接口，可让 AI 自动生成片段和连接图。

默认关闭：

```yaml
chemistry:
  ai_enabled: false
```

启用方式：

```powershell
$env:OPENAI_API_KEY="your_api_key"
```

然后修改配置：

```yaml
chemistry:
  ai_enabled: true
  ai_api_url: https://api.openai.com/v1/chat/completions
  ai_api_key_env: OPENAI_API_KEY
  ai_model: gpt-4.1-mini
  ai_output_file: ai_fragments.json
  ai_refresh: true
  ai_temperature: 0.2
  ai_max_fragments: 8
  ai_timeout_seconds: 60
  ai_validate_rdkit: true
  ai_prompt: Generate low thermal conductivity polymer fragments for MCTS screening.
```

AI 应返回 JSON：

```json
{
  "fragments": [
    {
      "key": "AI_FLEX_ETHER",
      "smiles": "[1*]CCOCC[2*]",
      "family": "ai_custom",
      "roles": ["flexible", "ether"]
    }
  ],
  "transitions": {
    "AI_FLEX_ETHER": ["AI_FLEX_ETHER", "End"]
  }
}
```

程序会校验：

- 是否为 JSON。
- 是否包含 `fragments` 和 `transitions`。
- 每个 SMILES 是否包含 `[1*]` 和 `[2*]`。
- 如果 RDKit 可用，是否能解析 SMILES。

通过校验后，AI 片段会和内置片段、自定义片段一起进入 MCTS 搜索空间。

更多说明见 `AI_FRAGMENT_API.md`。

## Packmol 初始体系

当配置为：

```yaml
system:
  prepare_initial_system: true
  molecule_count: 10
  box_size: 60.0
  tolerance: 2.0
```

程序会：

1. 将单链 PDB 转换为单分子 LAMMPS data。
2. 调用 Packmol 将多条链装入盒子。
3. 将 Packmol 输出 PDB 转换为体系 LAMMPS data。

输出文件示例：

```text
outputs/systems/candidate_001_single.data
outputs/systems/candidate_001_system.pdb
outputs/systems/candidate_001_system.data
```

## 快速热导率计算

热导率输入脚本由 `md_engine/thermal_conductivity_writer.py` 生成。

默认写出：

```text
outputs/lammps/candidate_001_fast_tc.in
```

手动运行：

```powershell
C:\DuskORI\Application\LAMMPS\bin\lmp.exe -in outputs\lammps\candidate_001_fast_tc.in
```

自动运行：

```yaml
lammps:
  run_enabled: true
```

快速热导率脚本使用 NEMD 思路，设置 hot/cold 区域，输出热流和温度剖面。

当前快速脚本主要用于流程验证和候选初筛。由于自动生成的 data 文件尚未包含完整真实力场，脚本默认使用 `pair_style zero` 并忽略分子间非键相互作用。若要做严格物理结论，需要后续接入真实力场参数、充分弛豫和重复采样。

## 输出文件

默认输出目录：

```text
outputs/
```

主要文件：

```text
outputs/candidates.csv
outputs/build_results.csv
outputs/system_results.csv
outputs/thermal_inputs.csv
outputs/lammps_runs.csv
outputs/thermal_results.csv
outputs/low_k_database.csv
outputs/feedback_records.csv
outputs/candidate_notes.txt
```

说明：

- `candidates.csv`：MCTS 候选序列、访问次数、平均奖励和结构特征。
- `build_results.csv`：RDKit 建链结果、SMILES 和 PDB 路径。
- `system_results.csv`：Packmol 初始体系和 LAMMPS data 结果。
- `thermal_inputs.csv`：快速热导率 LAMMPS 输入脚本路径。
- `lammps_runs.csv`：LAMMPS 运行记录。
- `thermal_results.csv`：热流、温度梯度和热导率估算。
- `low_k_database.csv`：按低热导优先排序的候选数据库。
- `feedback_records.csv`：MCTS 热导反馈记录。

## 典型运行命令

只生成结构和 LAMMPS 输入：

```powershell
cd C:\DuskORI\Files\Code\LAMMPS_MCTS\LAMMPS_MCTS
C:\DuskORI\Application\miniconda\envs\polymer_mcts\python.exe main.py
```

手动跑某个候选热导率：

```powershell
C:\DuskORI\Application\LAMMPS\bin\lmp.exe -in outputs\lammps\candidate_001_fast_tc.in
```

## 当前注意事项

- 当前酚醛树脂按线型片段模型处理，不是完整交联网络模型。
- 当前聚氨酯片段用于快速搜索，不显式模拟真实逐步加成反应。
- AI 生成片段需要人工复核，不建议直接作为最终化学结论。
- 快速热导率结果适合候选筛选，不适合作为生产级模拟结果。
- 生产计算需要补充真实力场、充分弛豫、统计误差和重复采样。
