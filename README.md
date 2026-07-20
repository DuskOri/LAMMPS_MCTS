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
- 按重复单元轮廓长度自动确定聚合度，使链长接近目标值，默认约 `120 Å`。
- 输出保留键级的聚合物 MOL 文件和用于装箱的 PDB 文件。
- 将 Packmol 坐标与 RDKit UFF 参数组合为筛选级 LAMMPS data 文件。
- 调用 Packmol 生成多链初始体系。
- 生成冷热区直接 NEMD 热导率 LAMMPS 输入脚本。
- 可选择自动运行 LAMMPS。
- 解析中央区热流和分层温度梯度，估算热导率。
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

推荐使用 conda 环境运行。环境名称可以自定义，下面以 `polymer_mcts` 为例：

```powershell
conda create -n polymer_mcts python=3.11
conda activate polymer_mcts
pip install -r requirements.txt
```

如果当前终端不能正常激活 conda，也可以使用 `conda run`：

```powershell
conda run -n polymer_mcts python main.py
```

Python 依赖写在 `requirements.txt`：

```text
rdkit
pyyaml
```

外部程序：

- Packmol：用于生成多链初始体系。
- LAMMPS：用于运行快速热导率计算。

项目默认通过系统 `PATH` 查找 Packmol 和 LAMMPS，因此仓库配置不包含任何开发者电脑的绝对路径：

```yaml
system:
  packmol_executable: packmol

lammps:
  executable: lmp
```

请先在终端运行 `packmol` 和 `lmp -help`，确认两个命令能够被当前环境找到。Windows、Linux 和 macOS 均使用上述配置；如果安装包提供了不同的命令名，只需把对应字段改成实际命令名。确实无法配置系统 `PATH` 时，也可以在自己的 `config.yaml` 中填写本机可执行文件地址，但不要把个人绝对路径提交到公共仓库。

## 快速启动

进入项目目录：

```powershell
cd LAMMPS_MCTS
```

运行主程序：

```powershell
conda activate polymer_mcts
python main.py
```

默认流程会执行：

1. 读取 `config.yaml`。
2. 加载内置片段和自定义片段。
3. 从 `Start` 节点开始 MCTS 搜索。
4. 生成候选聚合物片段序列。
5. 根据目标链长自动计算聚合度并构建聚合物链。
6. 输出 PDB。
7. 调用 Packmol 生成初始体系。
8. 写出包含密度平衡和直接 NEMD 的 LAMMPS 输入文件。
9. 如果配置允许，则依次执行高温熔融、预压缩、NPT 密度平台判断、冷却和平衡。
10. 在平衡后的密度下建立冷热区温差并运行 NEMD 采样。
11. 解析热导率结果并写入 CSV。

## 前端控制台

项目提供一个轻量静态前端，位于仓库的 `frontend/` 目录。

启动本地前端：

```powershell
cd LAMMPS_MCTS
conda activate polymer_mcts
python web_server.py
```

浏览器打开：

```text
http://localhost:8010/frontend/
```

前端当前提供：

- 新手向运行页。
- 目标链长、迭代次数、候选数量、最大片段步数修改。
- 已有有向图/片段空间选择。
- 内置家族与自定义片段空间互斥；选择聚酰亚胺时不会加载 `CUSTOM_*` 片段。
- 自定义有向图创建。
- 一键运行完整流程。
- 当前运行阶段显示。
- 按“第 N / 总轮数”显示循环进度，并在每次评价后即时更新实时 Top K。
- NEMD 提供“快速筛选”和“标准复核”两档；快速档用于 MCTS，标准档用于最终候选复算。
- 可开关密度平台平衡，并设置预压缩密度、平台判断阈值和最大判断窗口数。
- 前端可单独设置 LAMMPS 超时秒数；`0` 表示不限制运行时间。
- MCTS 候选序列和低热导数据库查看。
- 热导率结果表。
- 候选序列选择性清除。
- API 生成化合物片段配置。
- 主程序、前端服务器和 LAMMPS 命令复制。
- 项目内原创 Logo、聚合物流程图和功能图标。

新版前端通过 `web_server.py` 提供本地 API，可以修改 `config.yaml`、保存自定义有向图，并在后台启动主流程。服务默认只绑定当前设备，API key 不会写入项目文件。`localhost` 是操作系统统一提供的本机回环名称，不是某台开发电脑的文件路径。

有向图编辑器可以人工输入片段 key、SMILES、family、roles，并添加有向边。生成结果可复制或下载为 `custom_fragments.json`。API 生成器可以调用 OpenAI-compatible 接口生成片段和连接图，输出可下载为 `ai_fragments.json`。API key 只保存在浏览器当前页面，不写入项目文件。

`frontend/assets/` 中的 SVG 素材为项目内自绘素材，没有引入第三方图片文件，后续整理软著或商用说明时更容易说明来源。

## 详细使用说明

本项目有两种使用方式：一种是通过浏览器前端操作，适合初次使用和演示；另一种是直接运行 `main.py`，适合批量计算和命令行调试。

### 方式一：使用前端控制台

1. 打开 PowerShell，进入项目目录：

```powershell
cd LAMMPS_MCTS
```

2. 启动本地服务：

```powershell
conda activate polymer_mcts
python web_server.py
```

3. 在浏览器访问：

```text
http://localhost:8010/frontend/
```

4. 在“运行”页面设置参数：

- `目标链长`：默认 `120 Å`。程序根据每个候选重复单元的化学键轮廓长度自动计算聚合度。
- `MCTS 迭代次数`：控制搜索次数，数值越大搜索越充分，耗时也越长。
- `输出候选数`：控制最终保留多少个候选结构。
- `重复单元最大片段数`：控制一个候选重复单元中的真实片段数量，`Start` 和 `End` 不计数。
- `选择已有有向图/片段空间`：每次选择聚酰亚胺、聚氨酯、酚醛树脂中的一个，或者使用自定义片段。
- `自动运行 LAMMPS`：勾选后会在生成输入文件后尝试直接调用 LAMMPS。
- `热导结果反馈 MCTS`：勾选后使用热导率结果作为搜索反馈，计算成本更高。

5. 点击“保存参数”，前端会把常用参数写入 `config.yaml`。

6. 点击“一键运行”，程序会按阶段执行：

```text
MCTS 搜索 -> RDKit 建链 -> Packmol 初始体系 -> 高温熔融与预压缩 -> NPT 密度平台 -> 冷却平衡 -> 直接 NEMD -> 保存数据库
```

7. 运行结束后，在“数据库”页面查看候选序列和热导率结果。低热导数据库会优先展示热导率较低的候选。

8. 如果需要清理候选序列，可以进入“数据库”页面，勾选候选行后点击“清除选中”。也可以点击“清空候选”清除全部候选序列。清除前程序会自动把原始 CSV 备份到 `outputs/archive/`。

### 方式二：命令行运行

命令行方式适合不打开前端时使用：

```powershell
cd LAMMPS_MCTS
conda activate polymer_mcts
python main.py
```

命令行运行会读取 `config.yaml`，然后自动完成搜索、建链、装盒、LAMMPS 输入生成和结果整理。是否自动运行 LAMMPS 由配置项控制：

```yaml
mcts:
  feedback_run_lammps: true
```

如果只想快速检查片段搜索和结构导出，可以临时改为 `feedback_run_lammps: false`。如果本机 LAMMPS 路径已经配置正确，保持 `true` 即可在 MCTS 搜索阶段用热导率反馈 reward。

### 修改搜索参数

常用参数集中在 `config.yaml`：

```yaml
mcts:
  iterations: 300
  max_steps: 5
  top_k: 5
  feedback_mode: thermal
  feedback_every_iteration: true
  feedback_max_evaluations: 50
  feedback_run_lammps: true

polymer:
  target_chain_length_angstrom: 120.0
  min_degree_of_polymerization: 1
  max_degree_of_polymerization: 100
```

建议调参顺序：

1. 初次测试时使用较小的 `iterations`，例如 20 到 50。
2. 确认 RDKit、Packmol 和 LAMMPS 输入生成正常后，再提高到 300 或更高。
3. 如果候选重复单元过短，可以提高 `max_steps`，达到上限后程序会自动补上 `End`。
4. 如果需要更长聚合物链，可以提高 `target_chain_length_angstrom`；实际 DP 会随候选重复单元自动变化。
5. 当前默认使用 `feedback_mode: thermal`，热导率结果会在 MCTS 每轮评价后回传到节点。
6. 如果只想快速调试流程，可以临时改成 `feedback_mode: heuristic`。

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

### 热导率反馈

当前主流程把热导率计算放在 MCTS 搜索阶段。配置为：

```yaml
mcts:
  feedback_mode: thermal
  feedback_run_lammps: true
```

时，每一轮 MCTS rollout 会生成候选序列，调用 LAMMPS 快速热导脚本，读取热导率并转换为 reward 回传到节点。

搜索结束后，程序只会对最终 top_k 候选导出结构和 LAMMPS 输入文件，例如：

```text
outputs/lammps/candidate_001_direct_nemd.in
```

需要人工复查时，可以手动运行某个导出的候选：

```powershell
lmp -in outputs/lammps/candidate_001_direct_nemd.in
```

默认直接 NEMD 脚本保留分子间相互作用。它沿 x 方向设置固定端、热源区、中央传热区、冷源区和另一固定端；冷热 Langevin fix 挂在 `all` 组上，并通过 `fix_modify ... temp/region` 使用对应区域温度。输入 data 必须包含真实的 Pair、Bond、Angle、Dihedral 等力场系数；缺少所需系数时程序会停止该次热导评价并记录原因。

NEMD 采样前会先建立密度平台。程序在高温下熔融体系，用 `fix deform` 将明显疏松的 Packmol 盒子保守预压缩到筛选起点，再通过常压 NPT 让盒子尺寸和密度自行变化。预压缩密度不是聚合物的指定最终密度；最终密度由当前结构、力场、温度和压力共同决定。

密度平台由连续三个采样窗口的平均密度判断。快速档默认相对波动小于 `5%` 时停止，最多检查 `10` 个窗口；标准档默认阈值为 `2%`，最多检查 `20` 个窗口。达到最大窗口数仍未收敛时，日志会输出警告并继续后续冷却与筛选，正式计算前应延长平衡时间重新检查。

### 查看结果

常用结果文件如下：

- `outputs/candidates.csv`：MCTS 搜索到的候选片段序列。
- `outputs/build_results.csv`：RDKit 建链结果和 PDB 路径。
- `outputs/system_results.csv`：Packmol 初始体系和 LAMMPS data 结果。
- `outputs/thermal_inputs.csv`：最终 top_k 候选的 LAMMPS 输入脚本路径。
- `outputs/feedback_records.csv`：MCTS 搜索阶段的热导率反馈记录。
- `outputs/low_k_database.csv`：按搜索阶段 reward 排序的低热导候选数据库。
- `outputs/mcts_feedback/eval_NNNN/lammps/direct_nemd_density.dat`：本轮高温 NPT 的密度窗口均值。
- `outputs/mcts_feedback/eval_NNNN/lammps/direct_nemd_temperature_profile.dat`：x 方向分层温度。
- `outputs/mcts_feedback/eval_NNNN/lammps/direct_nemd_heat_flux.dat`：中央测量区的 x 方向热流密度。

如果只关心最终低热导候选，优先查看 `outputs/low_k_database.csv`。如果要排查流程问题，则按上述文件顺序逐步检查。

LAMMPS 运行时会把阶段命令、thermo 数据、警告和错误实时输出到启动程序的终端。每次 MCTS 热导评价的完整日志保存在 `outputs/mcts_feedback/eval_NNNN/lammps/direct_nemd.log`；LAMMPS 自己生成的通用 `log.lammps` 位于启动命令时的工作目录，通常是项目根目录。排错时优先查看候选目录中的 `direct_nemd.log`，它不会被下一候选覆盖。

前端支持清除候选序列。清除操作只重写对应 CSV 表格，不删除 `outputs/pdb/`、`outputs/systems/` 或 `outputs/lammps/` 中已经生成的结构和输入文件。需要彻底清理结构文件时，建议先手动备份再处理。

### 常见问题

- 前端打不开：确认 `web_server.py` 是否正在运行，端口是否为 `8010`。
- 更新项目代码后前端仍表现异常：先停止旧的 `web_server.py` 进程，再重新运行 `python web_server.py`。Python 服务不会自动重新导入磁盘上刚修改的后端模块。
- RDKit 报错：确认当前环境为 `polymer_mcts`，并且已安装 RDKit。
- Packmol 失败：先检查反馈信息是否停在 `preparing Packmol system`。若提示 `Packmol must be run with: packmol < inputfile.inp`，请重启前端服务以加载当前版本；若提示装箱空间不足，再适当增大 `system.box_size` 或降低 `system.molecule_count`。
- LAMMPS 没有运行：先确认终端可以执行 `lmp -help`，再检查 `lammps.executable` 命令名，并将 `mcts.feedback_run_lammps` 改为 `true`。
- 没有热导率结果：先检查 `outputs/feedback_records.csv` 中的失败原因，再检查 LAMMPS 是否输出 `*_temperature_profile.dat` 和 `*_heat_flux.dat`。温度梯度线性度或温差不足时结果会主动判为失败。
- AI 片段不可用：确认 API key、API URL、模型名称和返回 JSON 格式是否正确。

## 主要配置

核心配置在 `config.yaml`。

```yaml
mcts:
  iterations: 300
  max_steps: 5
  start_fragment: Start
  top_k: 5
  exploration_weight: 1.41421356237
  rollout_limit: 32
  random_seed: 7
  feedback_mode: thermal
  feedback_every_iteration: true
  feedback_max_evaluations: 50
  feedback_run_lammps: true

polymer:
  target_chain_length_angstrom: 120.0
  min_degree_of_polymerization: 1
  max_degree_of_polymerization: 100
  build_pdb: true
```

其中：

- `iterations`：MCTS 搜索迭代次数。
- `max_steps`：一个重复单元允许的最大真实片段数，不包含 `Start` 和 `End`，也不是聚合度。
- `start_fragment`：搜索起点，当前为 `Start`。
- `top_k`：输出候选数量。
- `target_chain_length_angstrom`：目标链轮廓长度，默认 `120 Å`。
- `min_degree_of_polymerization`、`max_degree_of_polymerization`：自动计算 DP 的上下限，用于防止异常短或异常长的结构。

程序先沿 `[1*]` 到 `[2*]` 之间的最短化学键路径计算重复单元轮廓长度，键长优先采用 RDKit UFF 平衡键长。随后选择使整条链最接近目标长度的整数 DP。`build_results.csv` 会保存实际 DP、重复单元长度和估算链长，便于人工复核。这里的 `120 Å` 指化学键路径的轮廓长度，不是折叠构象两端的空间直线距离。

## MCTS 与低热导目标

MCTS 使用 UCB 公式：

```text
UCB = 平均奖励项 + 探索项
```

程序中的奖励值被设计为低热导优先。对于真实或估计热导率 `kappa`，反馈奖励可写成：

```text
reward = 1 / (1 + kappa / reward_scale)
```

因此：

- 热导率越低，reward 越高。
- MCTS 越倾向于保留低热导候选。
- 候选数据库按低热导优先排序。

当前支持两种反馈模式：

```yaml
mcts:
  feedback_mode: thermal
  feedback_every_iteration: true
  feedback_max_evaluations: 50
  feedback_run_lammps: true
```

`thermal` 会在 MCTS 过程中把快速热导率计算结果反馈到节点。`feedback_every_iteration: true` 时，程序会让热导反馈次数至少等于 `iterations`，即每一轮 MCTS 都有热导 reward 参与回传。若同一条序列已经算过，程序会直接使用缓存结果，不会重复运行同一条序列。

如果只想快速检查流程，可以临时使用：

```yaml
mcts:
  feedback_mode: heuristic
```

`heuristic` 使用片段结构特征进行快速估计，不调用 LAMMPS，适合调试和大规模预筛。

## 聚合物片段库

当前内置片段家族包括：

- `polyimide`：聚酰亚胺。
- `polyurethane`：聚氨酯。
- `phenolic`：酚醛树脂。
- `custom`：用户自定义片段。

配置入口：

```yaml
chemistry:
  active_families: polyimide
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
outputs/lammps/candidate_001_direct_nemd.in
```

手动运行：

```powershell
lmp -in outputs/lammps/candidate_001_direct_nemd.in
```

自动反馈：

```yaml
mcts:
  feedback_mode: thermal
  feedback_run_lammps: true
```

快速热导率脚本先执行密度平衡，再进入直接 NEMD。默认过程为高温升温与熔融、疏松盒子预压缩、高温常压 NPT 密度平台检测、降温 NPT、目标温度 NPT/NVT，然后按照参考 `in.pur.lmp` 在 x 方向建立固定端、`450 K` 热源区、中央传热区和 `150 K` 冷源区。主体使用 NVE 积分，先运行稳态阶段，再同时输出中央区 `Jx` 和 `0.5 A` 分层温度。后处理拟合中间 60% 区域的温度梯度，并按 `kappa = |Jx| / |dT/dx|` 换算为 `W/(m K)`。

快速档默认使用 `100000` 步稳态和 `100000` 步生产采样，要求温度梯度拟合 `R2 >= 0.70`；标准档对应参考脚本使用 `500000 + 500000` 步，并要求 `R2 >= 0.85`。两档都要求拟合区温差至少 `5 K`，否则本轮标记失败，不向 MCTS 回传伪热导率。旧 `rapid_green_kubo` 方法仍保留用于历史配置兼容，但不再作为默认 reward。

热导流程参数位于 `thermal_conductivity` 段：

```yaml
thermal_conductivity:
  method: direct_nemd
  density_equilibration: true
  melt_temperature: 550.0
  precompression_density: 0.7
  density_block_steps: 2000
  density_sample_every: 100
  density_max_blocks: 10
  density_plateau_tolerance: 0.05
  nemd_hot_temperature: 450.0
  nemd_cold_temperature: 150.0
  nemd_steady_steps: 100000
  production_steps: 100000
  nemd_min_gradient_r2: 0.7
```

`precompression_density` 只控制预压缩起点，而且不会让本来更致密的体系反向膨胀。`density_plateau_tolerance` 是三个连续窗口的相对波动阈值，不是目标密度。最终密度会写入 LAMMPS 日志，完整的窗口历史写入 `direct_nemd_density.dat`。

默认的 `system.force_field: uff_screening` 会从 RDKit MOL 文件读取键级，计算 Gasteiger 电荷，并写入 UFF 的非键、键、键角和二面角参数。Packmol 只负责多链装箱，随后程序把这些参数和装箱坐标组合成可直接读取的 LAMMPS data。若选择 `topology_only`，或 data 缺少必要的 `Pair Coeffs`、`Bond Coeffs` 等段，热导评价会明确停止，不会生成零作用力的伪结果。

UFF 路径用于快速筛选和程序联调，不等同于论文采用的 GAFF 参数化，也不应直接作为定量热导率结论。正式计算建议接入 GAFF/AM1-BCC 或经过验证的聚合物力场，并增加平衡时间、生产步数和独立重复轨迹。

当前暂不执行聚合物家族官能团检查。内置和自定义候选只受片段有向图、最大片段数及 RDKit 结构解析约束，正式研究前需要人工筛除不符合反应与合成规则的重复单元。

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
outputs/low_k_database.csv
outputs/feedback_records.csv
outputs/candidate_notes.txt
```

说明：

- `candidates.csv`：MCTS 候选序列、访问次数、平均奖励和结构特征。
- `build_results.csv`：RDKit 建链结果、SMILES 和 PDB 路径。
- `system_results.csv`：Packmol 初始体系和 LAMMPS data 结果。
- `thermal_inputs.csv`：最终 top_k 候选的 LAMMPS 输入脚本路径。
- `low_k_database.csv`：按搜索阶段 reward 排序的低热导候选数据库。
- `feedback_records.csv`：MCTS 搜索阶段的热导反馈记录。

## 典型运行命令

只生成结构和 LAMMPS 输入：

```powershell
cd LAMMPS_MCTS
conda activate polymer_mcts
python main.py
```

手动跑某个候选热导率：

```powershell
lmp -in outputs/lammps/candidate_001_direct_nemd.in
```

## 当前注意事项

- 当前酚醛树脂按线型片段模型处理，不是完整交联网络模型。
- 当前聚氨酯片段用于快速搜索，不显式模拟真实逐步加成反应。
- AI 生成片段需要人工复核，不建议直接作为最终化学结论。
- 快速热导率结果适合候选筛选，不适合作为生产级模拟结果。
- 生产计算需要补充真实力场、充分弛豫、统计误差和重复采样。

## 当前行数统计

源码、配置和说明文件合计 `10768` 行，其中 Python `5349` 行、Markdown `1955` 行、HTML `435` 行、CSS `1684` 行、JavaScript `1146` 行，其余为 SVG、JSON、YAML 和 TXT。统计不包含 `.git/`、`outputs/` 与 `__pycache__/` 中的运行产物。
