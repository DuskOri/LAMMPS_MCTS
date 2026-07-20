# AI Fragment API

本项目提供一个 OpenAI-compatible Chat Completions 接口，用于自动生成聚合物片段和 MCTS 连接图。

接口默认关闭。关闭时程序不联网，也不会消耗 API 调用。

## 配置

在 `config.yaml` 中使用：

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

`ai_api_url` 兼容 OpenAI 风格接口。如果使用本地大模型或其他服务，只要它支持 `/v1/chat/completions` 返回格式，也可以替换这个地址。

## API Key

PowerShell 中设置：

```powershell
$env:OPENAI_API_KEY="your_api_key"
```

也可以换一个环境变量名：

```yaml
chemistry:
  ai_api_key_env: MY_POLYMER_AI_KEY
```

程序只读取环境变量，不把 key 写入输出文件。

## AI 输出格式

AI 必须返回 JSON：

```json
{
  "fragments": [
    {
      "key": "AI_FLEX_ETHER",
      "smiles": "[1*]CCOCC[2*]",
      "family": "ai_custom",
      "label": "flexible ether segment",
      "group": "soft_segment",
      "roles": ["flexible", "ether"],
      "notes": "AI generated soft segment"
    }
  ],
  "transitions": {
    "AI_FLEX_ETHER": ["AI_FLEX_ETHER", "End"]
  }
}
```

片段要求：

- 每个 `smiles` 必须正好包含一个 `[1*]` 和一个 `[2*]`。
- `key` 会作为 MCTS 节点名称，建议使用 `AI_` 前缀。
- `roles` 会参与低热导启发式评分。
- `transitions` 表示连接图，目标只能是另一个片段 key 或 `End`。

## 主流程

运行：

```powershell
cd LAMMPS_MCTS
conda activate polymer_mcts
python main.py
```

流程如下：

1. 读取 `config.yaml`。
2. 如果 `ai_enabled: true`，调用 AI 接口生成片段和连接图。
3. 保存到 `ai_fragments.json`。
4. 校验片段是否含有 `[1*]` 和 `[2*]`。
5. 如果安装了 RDKit，会进一步检查 SMILES 能否解析。
6. 将 AI 片段、`custom_fragments.json` 和内置片段合并。
7. MCTS 从合并后的 `Start` 连接图开始搜索。

如果 `ai_refresh: false` 且 `ai_fragments.json` 已存在，程序会复用已有文件，不重复请求 API。

## 建议

AI 生成的片段适合扩大搜索空间，但不应不经检查直接用于正式结论。

推荐工作流：

1. 先用 AI 生成一批片段。
2. 人工检查 `ai_fragments.json` 中的 SMILES 和 transitions。
3. 用 `DP=1` 或 `DP=2` 做小规模 RDKit 拼接测试。
4. 确认可建模后，再提高 DP、top_k 和 Packmol 分子数。
5. 最终热导率结论仍以 LAMMPS 计算和后处理结果为准。
