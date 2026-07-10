"""AI 片段生成接口。

该模块使用 OpenAI-compatible Chat Completions 接口生成聚合物片段库。
接口默认不参与主流程，只有在 config.yaml 中打开 chemistry.ai_enabled 后才会调用。

AI 返回的数据必须是 JSON，结构与 custom_fragments.json 保持一致：

{
  "fragments": [
    {
      "key": "CUSTOM_NAME",
      "smiles": "[1*]... [2*]",
      "family": "custom",
      "label": "readable name",
      "group": "soft_segment",
      "roles": ["flexible"]
    }
  ],
  "transitions": {
    "CUSTOM_NAME": ["CUSTOM_NAME", "End"]
  }
}

这里不把 AI 结果直接写入代码，而是保存成 JSON 文件。这样便于人工复核，也方便
把同一批片段反复用于 MCTS 搜索、Packmol 建模和 LAMMPS 快速筛选。
"""

import json
import os
import re
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

try:
    from rdkit import Chem
except ImportError:
    Chem = None


DEFAULT_SYSTEM_PROMPT = (
    "You are a polymer informatics assistant. Generate chemically plausible linear "
    "polymer fragments for RDKit chain construction. Every fragment SMILES must "
    "contain exactly two dummy connection atoms, [1*] and [2*]. Return JSON only."
)


DEFAULT_USER_PROMPT = (
    "Generate low-thermal-conductivity polymer fragments and a transition graph. "
    "Include aromatic, flexible, heteroatom-rich, polyurethane, phenolic, or custom "
    "segments when useful. Avoid charged fragments unless requested."
)


@dataclass
class AIFragmentResult:
    """记录一次 AI 片段生成的结果。"""

    success: bool
    data: dict
    output_path: str
    message: str = ""


def generate_ai_fragment_file(config, project_root="."):
    """调用 AI API 并把片段 JSON 写入文件。

    如果 ai_refresh 为 false 且输出文件已存在，则直接复用已有文件。
    这样可以避免每次运行 main.py 都消耗 API 调用。
    """
    project_root = Path(project_root)
    output_path = _resolve_path(project_root, config.get("ai_output_file", "ai_fragments.json"))
    refresh = _as_bool(config.get("ai_refresh", False))
    if output_path.exists() and not refresh:
        return AIFragmentResult(True, _read_json(output_path), str(output_path), "reuse existing ai fragment file")

    api_key = os.environ.get(str(config.get("ai_api_key_env", "OPENAI_API_KEY")))
    if not api_key:
        return AIFragmentResult(False, {}, str(output_path), "AI API key environment variable is not set")

    request_body = _build_request_body(config)
    raw_text = _call_chat_completions(
        api_url=str(config.get("ai_api_url", "https://api.openai.com/v1/chat/completions")),
        api_key=api_key,
        request_body=request_body,
        timeout_seconds=float(config.get("ai_timeout_seconds", 60)),
    )
    data = parse_ai_fragment_response(raw_text)
    data = validate_ai_fragment_data(data, require_rdkit=_as_bool(config.get("ai_validate_rdkit", True)))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_payload = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "model": config.get("ai_model", ""),
        "source": "ai_fragment_api",
        **data,
    }
    output_path.write_text(
        json.dumps(output_payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return AIFragmentResult(True, output_payload, str(output_path), "generated")


def parse_ai_fragment_response(raw_text):
    """从 AI 回复里提取 JSON。

    有些模型会把 JSON 放在 ```json 代码块里，这里会自动剥离代码块。
    如果回复里混有解释文字，则截取第一段大括号包围的 JSON。
    """
    text = raw_text.strip()
    fence_match = re.search(r"```(?:json)?\s*(.*?)\s*```", text, flags=re.DOTALL | re.IGNORECASE)
    if fence_match:
        text = fence_match.group(1).strip()

    if not text.startswith("{"):
        start = text.find("{")
        end = text.rfind("}")
        if start >= 0 and end > start:
            text = text[start : end + 1]

    data = json.loads(text)
    if "custom_fragments" in data and "fragments" not in data:
        data["fragments"] = data["custom_fragments"]
    if "custom_transitions" in data and "transitions" not in data:
        data["transitions"] = data["custom_transitions"]
    return data


def validate_ai_fragment_data(data, require_rdkit=True):
    """校验 AI 返回的片段和连接图。

    这里做的是输入边界检查，不判断片段是否一定适合真实合成。
    后续 MCTS 和 RDKit 建链仍会继续筛掉失败结构。
    """
    fragments = data.get("fragments", [])
    transitions = data.get("transitions", {})
    if not isinstance(fragments, list):
        raise ValueError("AI fragment data must contain a list field named fragments.")
    if not isinstance(transitions, Mapping):
        raise ValueError("AI fragment data transitions must be an object.")

    cleaned_fragments = []
    seen_keys = set()
    for item in fragments:
        if not isinstance(item, Mapping):
            continue
        key = _clean_key(item.get("key", ""))
        smiles = str(item.get("smiles", "")).strip()
        if not key or not smiles:
            continue
        if key in seen_keys:
            continue
        if smiles.count("[1*]") != 1 or smiles.count("[2*]") != 1:
            continue
        if require_rdkit and Chem is not None and Chem.MolFromSmiles(smiles) is None:
            continue

        roles = item.get("roles", [])
        if isinstance(roles, str):
            roles = [role.strip() for role in roles.split(",") if role.strip()]
        elif not isinstance(roles, Sequence):
            roles = []

        cleaned_fragments.append(
            {
                "key": key,
                "smiles": smiles,
                "family": str(item.get("family", "ai_custom")).strip() or "ai_custom",
                "label": str(item.get("label", key)).strip() or key,
                "group": str(item.get("group", "ai")).strip() or "ai",
                "roles": [str(role).strip() for role in roles if str(role).strip()],
                "notes": str(item.get("notes", "AI generated fragment")).strip(),
            }
        )
        seen_keys.add(key)

    if not cleaned_fragments:
        raise ValueError("AI did not return any valid two-ended fragment.")

    valid_keys = {item["key"] for item in cleaned_fragments}
    cleaned_transitions = {}
    for source, targets in transitions.items():
        source_key = _clean_key(source)
        if source_key not in valid_keys or not isinstance(targets, Sequence) or isinstance(targets, str):
            continue
        valid_targets = []
        for target in targets:
            target_key = _clean_key(target)
            if target_key == "End" or target_key in valid_keys:
                valid_targets.append(target_key)
        if valid_targets:
            cleaned_transitions[source_key] = valid_targets

    for key in valid_keys:
        cleaned_transitions.setdefault(key, [other for other in valid_keys if other != key] + ["End"])

    return {"fragments": cleaned_fragments, "transitions": cleaned_transitions}


def _build_request_body(config):
    """组装 Chat Completions 请求体。"""
    max_fragments = int(config.get("ai_max_fragments", 8))
    user_prompt = str(config.get("ai_prompt", DEFAULT_USER_PROMPT)).strip() or DEFAULT_USER_PROMPT
    user_prompt = (
        f"{user_prompt}\n\n"
        f"Return at most {max_fragments} fragments. "
        "Use keys with prefix AI_. "
        "The JSON must contain fragments and transitions. "
        "Every transition target must be another key or End."
    )

    return {
        "model": config.get("ai_model", "gpt-4.1-mini"),
        "temperature": float(config.get("ai_temperature", 0.2)),
        "messages": [
            {"role": "system", "content": str(config.get("ai_system_prompt", DEFAULT_SYSTEM_PROMPT))},
            {"role": "user", "content": user_prompt},
        ],
    }


def _call_chat_completions(api_url, api_key, request_body, timeout_seconds):
    """发送 OpenAI-compatible Chat Completions 请求并返回文本内容。"""
    payload = json.dumps(request_body).encode("utf-8")
    request = urllib.request.Request(
        api_url,
        data=payload,
        method="POST",
        headers={
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
            response_data = json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="ignore")
        raise RuntimeError(f"AI API request failed: HTTP {exc.code} {detail}") from exc

    return response_data["choices"][0]["message"]["content"]


def _resolve_path(project_root, filename):
    """把配置中的相对路径解析到项目根目录。"""
    path = Path(str(filename))
    if path.is_absolute():
        return path
    return project_root / path


def _read_json(path):
    """读取 JSON 文件。"""
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _clean_key(value):
    """清理片段 key，只保留便于写 CSV 和连接图的字符。"""
    key = str(value).strip()
    key = re.sub(r"[^A-Za-z0-9_'\-]", "_", key)
    return key[:64]


def _as_bool(value):
    """把配置值转换成 bool。"""
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in ("1", "true", "yes", "on")
    return bool(value)

