const API_BASE = location.port === "8010" ? "" : "http://localhost:8010";

const state = {
  config: {},
  results: {
    candidates: [],
    thermal_results: [],
    low_k_database: []
  },
  graphDraft: {
    fragments: [],
    transitions: {}
  },
  selectedCandidates: new Set(),
  apiGenerated: null,
  pollTimer: null,
  latestStatus: null
};

const ROLE_OPTIONS = [
  "rigid",
  "aromatic",
  "flexible",
  "ether",
  "imide",
  "urethane",
  "phenolic",
  "aliphatic",
  "diol",
  "bridge",
  "polar"
];

const ROLE_TEXT = {
  rigid: "刚性片段",
  aromatic: "芳香结构",
  flexible: "柔性链段",
  ether: "醚键",
  imide: "酰亚胺",
  urethane: "氨基甲酸酯",
  phenolic: "酚醛/酚环",
  aliphatic: "脂肪链",
  diol: "二醇软段",
  bridge: "桥接单元",
  polar: "极性基团"
};

const stageOrder = ["search", "build", "packmol", "write_lammps", "compress", "run_lammps", "mcts_lammps", "export", "save"];
const stageText = {
  idle: "等待运行",
  search: "MCTS 搜索",
  build: "RDKit 建链",
  packmol: "Packmol 初始体系",
  write_lammps: "写出 LAMMPS 输入",
  compress: "高温熔融与体系压缩",
  run_lammps: "Green-Kubo 计算",
  mcts_lammps: "热导 reward 回传",
  export: "导出最终 Top K",
  save: "保存数据库",
  done: "完成",
  failed: "失败"
};

const thermalProfileDefaults = {
  quick: {
    timeout: 1800,
    precompressionDensity: 0.7,
    densityPlateauPercent: 5,
    densityMaxBlocks: 10
  },
  standard: {
    timeout: 7200,
    precompressionDensity: 0.7,
    densityPlateauPercent: 2,
    densityMaxBlocks: 20
  }
};

document.addEventListener("DOMContentLoaded", () => {
  bindTabs();
  bindActions();
  loadInitialData();
});

function bindTabs() {
  document.querySelectorAll(".nav-button").forEach((button) => {
    button.setAttribute("aria-selected", String(button.classList.contains("active")));
    button.addEventListener("click", () => {
      const tab = button.dataset.tab;
      document.querySelectorAll(".nav-button").forEach((item) => {
        item.classList.remove("active");
        item.setAttribute("aria-selected", "false");
      });
      document.querySelectorAll(".tab-panel").forEach((item) => item.classList.remove("active"));
      button.classList.add("active");
      button.setAttribute("aria-selected", "true");
      document.getElementById(tab).classList.add("active");
      if (tab === "database") loadResults();
    });
  });
}

function bindActions() {
  document.getElementById("saveConfigBtn").addEventListener("click", saveConfig);
  document.getElementById("reloadConfigBtn").addEventListener("click", loadConfig);
  document.getElementById("runNowBtn").addEventListener("click", runPipeline);
  document.getElementById("refreshResultsBtn").addEventListener("click", loadResults);
  document.getElementById("candidateFilter").addEventListener("input", renderCandidateTable);
  document.getElementById("iterationsInput").addEventListener("input", updateIdleIterationTarget);
  document.getElementById("gkProfileSelect").addEventListener("change", updateThermalProfileFields);
  document.getElementById("clearSelectedCandidatesBtn").addEventListener("click", clearSelectedCandidates);
  document.getElementById("clearAllCandidatesBtn").addEventListener("click", clearAllCandidates);

  document.getElementById("addGraphNodeBtn").addEventListener("click", addGraphNode);
  document.getElementById("loadGraphExampleBtn").addEventListener("click", loadGraphExample);
  document.getElementById("clearGraphBtn").addEventListener("click", clearGraphDraft);
  document.getElementById("addEdgeBtn").addEventListener("click", addGraphEdge);
  document.getElementById("addEndEdgeBtn").addEventListener("click", addEndEdge);
  document.getElementById("saveGraphBtn").addEventListener("click", saveGraphDraft);
  document.getElementById("generateByApiBtn").addEventListener("click", generateByApi);
  document.getElementById("applyApiGraphBtn").addEventListener("click", applyApiGeneratedToGraph);
  document.getElementById("downloadApiJsonBtn").addEventListener("click", () => downloadText("ai_fragments.json", getApiJson()));
  document.getElementById("graphNodeRoles").addEventListener("input", syncRoleChipsFromField);
  document.getElementById("graphNodeSmiles").addEventListener("input", refreshRoleHint);
  document.getElementById("inferRolesBtn").addEventListener("click", inferRolesFromCurrentSmiles);
  document.getElementById("inferRolesByApiBtn").addEventListener("click", inferRolesByApi);
  document.querySelectorAll("input[name='rewardMode']").forEach((input) => {
    input.addEventListener("change", updateRoleHelperText);
  });
}

async function loadInitialData() {
  await loadConfig();
  await loadResults();
  await refreshStatus();
  renderCommands();
  renderRoleChips();
  updateRoleHelperText();
  renderGraphEditor();
}

async function apiGet(path) {
  const response = await fetch(`${API_BASE}${path}`, { cache: "no-store" });
  if (!response.ok) throw new Error(`GET ${path} failed`);
  return response.json();
}

async function apiPost(path, payload) {
  const response = await fetch(`${API_BASE}${path}`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload || {})
  });
  if (!response.ok) throw new Error(`POST ${path} failed`);
  return response.json();
}

async function loadConfig() {
  try {
    state.config = await apiGet("/api/config");
    fillConfigForm(state.config);
  } catch (error) {
    showStatus("failed", "没有连接到本地接口，请先运行 python web_server.py");
  }
}

function fillConfigForm(config) {
  document.getElementById("targetLengthInput").value = config.target_length_angstrom || 120;
  document.getElementById("iterationsInput").value = config.iterations || 300;
  document.getElementById("explorationWeightInput").value = config.exploration_weight || 1.41421356237;
  document.getElementById("topKInput").value = config.top_k || 5;
  document.getElementById("maxStepsInput").value = config.max_steps || 5;
  const builtinFamilies = ["polyimide", "polyurethane", "phenolic", "custom"];
  const activeFamily = builtinFamilies.includes(config.active_families)
    ? config.active_families
    : "polyimide";
  document.getElementById("graphPresetSelect").value = activeFamily;
  document.getElementById("gkProfileSelect").value = config.gk_profile === "standard" ? "standard" : "quick";
  document.getElementById("lammpsTimeoutInput").value = Number(config.timeout_seconds) || 0;
  document.getElementById("densityEquilibrationInput").checked = config.density_equilibration !== false;
  document.getElementById("precompressionDensityInput").value = Number(config.precompression_density) || 0.7;
  document.getElementById("densityPlateauPercentInput").value = Number(config.density_plateau_percent) || 5;
  document.getElementById("densityMaxBlocksInput").value = Number(config.density_max_blocks) || 10;
  setRewardMode(config.feedback_mode === "thermal" && config.feedback_run_lammps !== false ? "thermal" : "heuristic");
}

function readConfigForm() {
  const rewardMode = getRewardMode();
  const useThermalReward = rewardMode === "thermal";
  return {
    target_length_angstrom: Number(document.getElementById("targetLengthInput").value || 120),
    iterations: Number(document.getElementById("iterationsInput").value || 300),
    exploration_weight: Number(document.getElementById("explorationWeightInput").value || 1.41421356237),
    top_k: Number(document.getElementById("topKInput").value || 5),
    max_steps: Number(document.getElementById("maxStepsInput").value || 5),
    active_families: document.getElementById("graphPresetSelect").value,
    gk_profile: document.getElementById("gkProfileSelect").value,
    timeout_seconds: Number(document.getElementById("lammpsTimeoutInput").value || 0),
    density_equilibration: document.getElementById("densityEquilibrationInput").checked,
    precompression_density: Number(document.getElementById("precompressionDensityInput").value || 0.7),
    density_plateau_percent: Number(document.getElementById("densityPlateauPercentInput").value || 5),
    density_max_blocks: Number(document.getElementById("densityMaxBlocksInput").value || 10),
    run_lammps: useThermalReward,
    feedback_mode: useThermalReward ? "thermal" : "heuristic",
    feedback_run_lammps: useThermalReward,
    feedback_every_iteration: useThermalReward,
    feedback_max_evaluations: Number(document.getElementById("iterationsInput").value || 300)
  };
}

function updateThermalProfileFields() {
  const profile = document.getElementById("gkProfileSelect").value;
  const defaults = thermalProfileDefaults[profile] || thermalProfileDefaults.quick;
  document.getElementById("lammpsTimeoutInput").value = defaults.timeout;
  document.getElementById("precompressionDensityInput").value = defaults.precompressionDensity;
  document.getElementById("densityPlateauPercentInput").value = defaults.densityPlateauPercent;
  document.getElementById("densityMaxBlocksInput").value = defaults.densityMaxBlocks;
}

function getRewardMode() {
  const selected = document.querySelector("input[name='rewardMode']:checked");
  return selected ? selected.value : "thermal";
}

function setRewardMode(mode) {
  const value = mode === "heuristic" ? "heuristic" : "thermal";
  const input = document.querySelector(`input[name='rewardMode'][value='${value}']`);
  if (input) input.checked = true;
  updateRoleHelperText();
}

async function saveConfig() {
  const payload = readConfigForm();
  try {
    const result = await apiPost("/api/config", payload);
    state.config = result.config;
    showStatus("idle", "参数已保存，可以运行。");
  } catch (error) {
    showStatus("failed", `保存失败：${error.message}`);
  }
}

async function runPipeline() {
  const payload = readConfigForm();
  try {
    await apiPost("/api/run", payload);
    showStatus("search", "已启动，正在搜索。");
    startPolling();
  } catch (error) {
    showStatus("failed", `启动失败：${error.message}`);
  }
}

function startPolling() {
  if (state.pollTimer) clearInterval(state.pollTimer);
  state.pollTimer = setInterval(refreshStatus, 1800);
  refreshStatus();
}

async function refreshStatus() {
  try {
    const status = await apiGet("/api/status");
    updateStatus(status);
    if (!status.running && state.pollTimer) {
      clearInterval(state.pollTimer);
      state.pollTimer = null;
      await loadResults();
    }
  } catch (error) {
    showStatus("failed", "本地接口未连接。");
  }
}

function updateStatus(status) {
  state.latestStatus = status;
  showStatus(status.stage, status.message || "");
  setRunButtonState(Boolean(status.running));
  renderLoopProgress(status);
  renderLiveTopK(status);
  if (status.counts) {
    const liveCount = Array.isArray(status.live_candidates) ? status.live_candidates.length : 0;
    document.getElementById("candidateCount").textContent = status.running
      ? liveCount
      : (status.counts.candidates || liveCount || 0);
    document.getElementById("buildCount").textContent = status.counts.build_results || 0;
    document.getElementById("inputCount").textContent = status.counts.thermal_inputs || 0;
  }
}

function updateIdleIterationTarget() {
  if (state.latestStatus?.running) return;
  const total = Math.max(0, Number(document.getElementById("iterationsInput").value) || 0);
  renderLoopProgress({
    ...(state.latestStatus || {}),
    iteration: 0,
    total_iterations: total,
    evaluation_count: 0
  });
}

function renderLoopProgress(status) {
  const iteration = Math.max(0, Number(status.iteration) || 0);
  const inputTotal = Number(document.getElementById("iterationsInput")?.value) || 0;
  const total = Math.max(0, Number(status.total_iterations) || inputTotal || Number(state.config.iterations) || 0);
  const evaluations = Math.max(0, Number(status.evaluation_count) || 0);
  const percent = total > 0 ? Math.min(100, (iteration / total) * 100) : 0;
  document.getElementById("iterationLabel").textContent = `第 ${iteration} / ${total} 轮`;
  document.getElementById("evaluationLabel").textContent = `已评价 ${evaluations} 个候选`;
  document.getElementById("iterationProgressBar").style.width = `${percent}%`;
}

function renderLiveTopK(status) {
  const topK = Math.max(1, Number(status.top_k) || Number(state.config.top_k) || 5);
  const rows = (status.live_candidates || []).slice(0, topK).map((row, index) => ({
    rank: index + 1,
    sequence_text: row.sequence_text || (row.sequence || []).join(" -> "),
    score: formatLiveNumber(row.score),
    conductivity_w_mk: row.conductivity_w_mk === undefined ? "--" : formatLiveNumber(row.conductivity_w_mk),
    visits: row.visits || 1,
    iteration: row.iteration || "--",
    success: row.success === undefined ? "--" : (row.success ? "成功" : "失败")
  }));
  renderTable("liveTopKTable", rows, [
    "rank",
    "sequence_text",
    "score",
    "conductivity_w_mk",
    "visits",
    "iteration",
    "success"
  ]);

  const meta = document.getElementById("liveTopKMeta");
  meta.textContent = rows.length
    ? `当前显示 ${rows.length} / ${status.live_candidates.length} 个已评价序列`
    : "等待产生候选";

  const validKappa = rows
    .filter((row) => row.success === "成功")
    .map((row) => Number(row.conductivity_w_mk))
    .filter((value) => Number.isFinite(value) && value >= 0);
  if (validKappa.length) {
    document.getElementById("bestKappa").textContent = `${Math.min(...validKappa).toFixed(4)} W/mK`;
  }
}

function formatLiveNumber(value) {
  const number = Number(value);
  if (!Number.isFinite(number)) return "--";
  return number.toFixed(6).replace(/0+$/, "").replace(/\.$/, "");
}

function showStatus(stage, message) {
  const displayStage = resolveDisplayStage(stage, message);
  document.getElementById("statusTitle").textContent = stageText[displayStage] || displayStage;
  document.getElementById("statusMessage").textContent = message || "";
  const badge = document.getElementById("statusBadge");
  if (badge) {
    badge.classList.remove("is-running", "is-done", "is-failed");
    const label = badge.querySelector("b");
    if (stage === "failed") {
      badge.classList.add("is-failed");
      label.textContent = "异常";
    } else if (stage === "done") {
      badge.classList.add("is-done");
      label.textContent = "完成";
    } else if (stage && stage !== "idle") {
      badge.classList.add("is-running");
      label.textContent = "运行中";
    } else {
      label.textContent = "就绪";
    }
  }
  renderStageList(displayStage);
}

function resolveDisplayStage(stage, message) {
  if (String(stage).startsWith("export_")) return "export";
  if (stage !== "mcts_lammps") return stage;

  const text = String(message || "").toLowerCase();
  if (text.includes("building polymer")) return "build";
  if (text.includes("packmol")) return "packmol";
  if (text.includes("writing lammps")) return "write_lammps";
  if (text.includes("compressing") || text.includes("density plateau")) return "compress";
  if (text.includes("running lammps")) return "run_lammps";
  return "mcts_lammps";
}

function setRunButtonState(running) {
  const button = document.getElementById("runNowBtn");
  const label = button.querySelector(".run-label");
  button.disabled = running;
  button.classList.toggle("is-running", running);
  label.textContent = running ? "正在运行" : "运行搜索";
}

function renderStageList(currentStage) {
  const currentIndex = stageOrder.indexOf(currentStage);
  document.querySelectorAll("#stageList li").forEach((item) => {
    const index = stageOrder.indexOf(item.dataset.stage);
    item.classList.remove("active", "done", "failed");
    if (currentStage === "failed") {
      item.classList.add("failed");
    } else if (index >= 0 && currentIndex >= 0 && index < currentIndex) {
      item.classList.add("done");
    } else if (item.dataset.stage === currentStage) {
      item.classList.add("active");
    } else if (currentStage === "done") {
      item.classList.add("done");
    }
  });
}

async function loadResults() {
  try {
    state.results = await apiGet("/api/results");
    renderAllTables();
  } catch (error) {
    state.results = { candidates: [], thermal_results: [], low_k_database: [] };
    renderAllTables();
  }
}

function renderAllTables() {
  renderTable("lowKTable", state.results.low_k_database || state.results.thermal_results || [], [
    "low_k_rank",
    "input_path",
    "conductivity_w_mk",
    "success",
    "message"
  ]);
  renderCandidateTable();
  renderTable("thermalTable", state.results.thermal_results || [], [
    "low_k_rank",
    "input_path",
    "conductivity_w_mk",
    "success",
    "message"
  ]);
  renderMetricsFromResults();
}

function renderCandidateTable() {
  const filter = document.getElementById("candidateFilter")?.value?.trim().toLowerCase() || "";
  const rows = (state.results.candidates || []).filter((row) => {
    if (!filter) return true;
    return Object.values(row).join(" ").toLowerCase().includes(filter);
  });
  renderCandidateSelectionTable(rows, [
    "rank",
    "sequence_text",
    "score",
    "visits",
    "length",
    "family_text",
    "flexible_ratio",
    "aromatic_ratio"
  ]);
  updateCandidateSelectionCount();
}

function renderCandidateSelectionTable(rows, preferredColumns) {
  const table = document.getElementById("candidateTable");
  const thead = table.querySelector("thead");
  const tbody = table.querySelector("tbody");
  thead.innerHTML = "";
  tbody.innerHTML = "";
  if (!rows || rows.length === 0) {
    tbody.innerHTML = "<tr><td>暂无数据</td></tr>";
    return;
  }

  const allColumns = Object.keys(rows[0]);
  const columns = preferredColumns.filter((key) => allColumns.includes(key));
  allColumns.forEach((key) => {
    if (!columns.includes(key) && columns.length < 9) columns.push(key);
  });

  const headerRow = document.createElement("tr");
  const selectTh = document.createElement("th");
  const selectAll = document.createElement("input");
  selectAll.type = "checkbox";
  selectAll.title = "选择当前显示的候选";
  selectAll.checked = rows.every((row) => state.selectedCandidates.has(candidateRowKey(row)));
  selectAll.addEventListener("change", () => {
    rows.forEach((row) => {
      const key = candidateRowKey(row);
      if (selectAll.checked) state.selectedCandidates.add(key);
      else state.selectedCandidates.delete(key);
    });
    renderCandidateTable();
  });
  selectTh.appendChild(selectAll);
  headerRow.appendChild(selectTh);
  columns.forEach((column) => {
    const th = document.createElement("th");
    th.textContent = column;
    headerRow.appendChild(th);
  });
  thead.appendChild(headerRow);

  rows.forEach((row) => {
    const key = candidateRowKey(row);
    const tr = document.createElement("tr");
    const selectTd = document.createElement("td");
    const checkbox = document.createElement("input");
    checkbox.type = "checkbox";
    checkbox.checked = state.selectedCandidates.has(key);
    checkbox.addEventListener("change", () => {
      if (checkbox.checked) state.selectedCandidates.add(key);
      else state.selectedCandidates.delete(key);
      updateCandidateSelectionCount();
    });
    selectTd.appendChild(checkbox);
    tr.appendChild(selectTd);
    columns.forEach((column) => {
      const td = document.createElement("td");
      td.textContent = row[column] || "";
      tr.appendChild(td);
    });
    tbody.appendChild(tr);
  });
}

function candidateRowKey(row) {
  return `${String(row.rank || "").trim()}|${String(row.sequence_text || "").trim()}`;
}

function updateCandidateSelectionCount() {
  const countNode = document.getElementById("candidateSelectionCount");
  if (countNode) countNode.textContent = `已选 ${state.selectedCandidates.size} 个候选`;
}

async function clearSelectedCandidates() {
  const selected = Array.from(state.selectedCandidates);
  if (selected.length === 0) {
    alert("请先勾选要清除的候选序列。");
    return;
  }
  if (!confirm(`确认清除 ${selected.length} 个候选序列？清除前会自动备份 candidates.csv。`)) return;
  await clearData({ mode: "selected_candidates", candidate_keys: selected });
}

async function clearAllCandidates() {
  if (!confirm("确认清空全部候选序列？清除前会自动备份 candidates.csv。")) return;
  await clearData({ mode: "all_candidates" });
}

async function clearData(payload) {
  try {
    const result = await apiPost("/api/clear-data", payload);
    if (result.ok === false) throw new Error(result.error || "clear data failed");
    state.selectedCandidates.clear();
    await loadResults();
    await refreshStatus();
    showStatus("idle", result.message || "数据已清除。");
  } catch (error) {
    showStatus("failed", `清除失败：${error.message}`);
  }
}

function renderMetricsFromResults() {
  const thermalRows = state.results.thermal_results || [];
  const values = thermalRows
    .filter((row) => String(row.success || "").toLowerCase() !== "false")
    .map((row) => Number(row.conductivity_w_mk))
    .filter((value) => Number.isFinite(value) && value > 0);
  document.getElementById("bestKappa").textContent = values.length ? `${Math.min(...values).toFixed(4)} W/mK` : "--";
}

function renderTable(tableId, rows, preferredColumns) {
  const table = document.getElementById(tableId);
  const thead = table.querySelector("thead");
  const tbody = table.querySelector("tbody");
  thead.innerHTML = "";
  tbody.innerHTML = "";
  if (!rows || rows.length === 0) {
    tbody.innerHTML = "<tr><td>暂无数据</td></tr>";
    return;
  }

  const allColumns = Object.keys(rows[0]);
  const columns = preferredColumns.filter((key) => allColumns.includes(key));
  allColumns.forEach((key) => {
    if (!columns.includes(key) && columns.length < 9) columns.push(key);
  });

  const headerRow = document.createElement("tr");
  columns.forEach((column) => {
    const th = document.createElement("th");
    th.textContent = column;
    headerRow.appendChild(th);
  });
  thead.appendChild(headerRow);

  rows.forEach((row) => {
    const tr = document.createElement("tr");
    columns.forEach((column) => {
      const td = document.createElement("td");
      td.textContent = row[column] || "";
      tr.appendChild(td);
    });
    tbody.appendChild(tr);
  });
}

function renderRoleChips() {
  const container = document.getElementById("roleChips");
  if (!container) return;
  container.innerHTML = "";
  ROLE_OPTIONS.forEach((role) => {
    const label = document.createElement("label");
    label.className = "role-chip";
    label.title = ROLE_TEXT[role] || role;

    const checkbox = document.createElement("input");
    checkbox.type = "checkbox";
    checkbox.value = role;
    checkbox.addEventListener("change", syncRoleFieldFromChips);

    const span = document.createElement("span");
    span.textContent = role;

    label.appendChild(checkbox);
    label.appendChild(span);
    container.appendChild(label);
  });
  syncRoleChipsFromField();
}

function readRolesFromField() {
  return normalizeRoles(document.getElementById("graphNodeRoles").value);
}

function setRolesToField(roles, message) {
  const cleaned = normalizeRoles(roles);
  document.getElementById("graphNodeRoles").value = cleaned.join(",");
  syncRoleChipsFromField();
  refreshRoleHint(message);
}

function syncRoleFieldFromChips() {
  const selected = Array.from(document.querySelectorAll("#roleChips input:checked")).map((item) => item.value);
  document.getElementById("graphNodeRoles").value = selected.join(",");
  refreshRoleHint();
}

function syncRoleChipsFromField() {
  const roles = new Set(readRolesFromField());
  document.querySelectorAll("#roleChips input").forEach((input) => {
    input.checked = roles.has(input.value);
  });
  refreshRoleHint();
}

function refreshRoleHint(message) {
  const node = document.getElementById("roleHint");
  if (!node) return;
  if (message) {
    node.textContent = message;
    return;
  }
  const roles = readRolesFromField();
  const mode = getRewardMode();
  if (roles.length) {
    node.textContent = `当前 roles：${roles.join(", ")}。`;
  } else if (mode === "heuristic") {
    node.textContent = "启发式 reward 依赖 roles；可以点选标签，或按 SMILES/API 自动判断。";
  } else {
    node.textContent = "LAMMPS reward 主要使用热导率结果，roles 可作为备用和结果说明。";
  }
}

function updateRoleHelperText() {
  refreshRoleHint();
}

function inferRolesFromCurrentSmiles() {
  const smiles = document.getElementById("graphNodeSmiles").value.trim();
  if (!smiles) {
    alert("请先填写片段 SMILES。");
    return;
  }
  const roles = inferRolesFromSmiles(smiles);
  setRolesToField(roles, `已按 SMILES 粗判：${roles.join(", ") || "未识别到明显标签"}。`);
}

async function inferRolesByApi() {
  const smiles = document.getElementById("graphNodeSmiles").value.trim();
  const apiKey = document.getElementById("apiKeyInput").value.trim();
  const apiUrl = document.getElementById("apiUrlInput").value.trim();
  const model = document.getElementById("apiModelInput").value.trim();
  if (!smiles) {
    alert("请先填写片段 SMILES。");
    return;
  }
  if (!apiKey || !apiUrl || !model) {
    alert("请先在右侧 API 区域填写 API URL、Model 和 API Key。");
    return;
  }

  refreshRoleHint("正在调用 API 判断 roles...");
  try {
    const response = await fetch(apiUrl, {
      method: "POST",
      headers: { "Content-Type": "application/json", "Authorization": `Bearer ${apiKey}` },
      body: JSON.stringify({
        model,
        temperature: 0,
        messages: [
          {
            role: "system",
            content: "Classify polymer fragment roles. Return JSON only. Use only the provided role vocabulary."
          },
          {
            role: "user",
            content: [
              `SMILES: ${smiles}`,
              `Role vocabulary: ${ROLE_OPTIONS.join(", ")}`,
              "Return format: {\"roles\":[\"flexible\",\"ether\"],\"reason\":\"short reason\"}.",
              "The fragment is used for heuristic low thermal conductivity screening."
            ].join("\n")
          }
        ]
      })
    });
    const data = await response.json();
    const content = data.choices?.[0]?.message?.content || JSON.stringify(data);
    const parsed = parseJsonFromText(content);
    const fallback = inferRolesFromSmiles(smiles);
    const roles = normalizeRoles(parsed.roles || parsed.role || [], smiles);
    setRolesToField(roles.length ? roles : fallback, parsed.reason ? `API 判断完成：${parsed.reason}` : "API 判断完成。");
  } catch (error) {
    const fallback = inferRolesFromSmiles(smiles);
    setRolesToField(fallback, `API 判断失败，已使用本地粗判：${error.message}`);
  }
}

function normalizeRoles(rawRoles, smiles = "") {
  let values = [];
  if (Array.isArray(rawRoles)) {
    values = rawRoles;
  } else if (typeof rawRoles === "string") {
    values = rawRoles.split(",");
  }
  const cleaned = values
    .map((item) => String(item).trim().toLowerCase())
    .filter(Boolean)
    .map((item) => item.replace(/\s+/g, "_"));
  const known = new Set(ROLE_OPTIONS);
  const roles = cleaned.filter((item) => known.has(item));
  if (roles.length === 0 && smiles) roles.push(...inferRolesFromSmiles(smiles));
  return Array.from(new Set(roles));
}

function inferRolesFromSmiles(smiles) {
  const text = String(smiles || "");
  const lower = text.toLowerCase();
  const roles = [];
  if (/c1|n1|o1|s1/.test(lower) || lower.includes("ccc")) roles.push("rigid", "aromatic");
  if (/cccc|cc\[2\*\]|\[1\*\]cc/.test(lower) && !roles.includes("aromatic")) roles.push("flexible", "aliphatic");
  if (lower.includes("o") && !lower.includes("n-c(=o)o")) roles.push("ether");
  if (/ccocc|coc|o\[2\*\]|\[1\*\]o/.test(lower)) roles.push("flexible", "ether");
  if (/n.*c\(=o\)|c\(=o\).*n/.test(lower)) roles.push("polar");
  if (/nc\(=o\)o|oc\(=o\)n/.test(lower)) roles.push("urethane", "polar");
  if (/c\(=o\)n.*c\(=o\)|n1.*c\(=o\)/.test(lower)) roles.push("imide", "rigid");
  if (lower.includes("(o)") || lower.includes("c(o)") || lower.includes("co)")) roles.push("phenolic");
  if (/c\[2\*\]|\[1\*\]c\[2\*\]/.test(lower)) roles.push("bridge");
  if (/cco|occ|cccc/.test(lower)) roles.push("flexible");
  return normalizeRoles(roles);
}

function addGraphNode() {
  const key = document.getElementById("graphNodeKey").value.trim();
  const smiles = document.getElementById("graphNodeSmiles").value.trim();
  const family = document.getElementById("graphNodeFamily").value.trim() || "custom";
  const roles = normalizeRoles(document.getElementById("graphNodeRoles").value, smiles);
  if (!key || !smiles) {
    alert("片段 key 和 SMILES 不能为空。");
    return;
  }
  if (!smiles.includes("[1*]") || !smiles.includes("[2*]")) {
    alert("SMILES 必须包含 [1*] 和 [2*]。");
    return;
  }
  if (roles.length === 0) {
    alert("请至少选择一个 roles，或点击 SMILES/API 自动判断。");
    return;
  }
  state.graphDraft.fragments = state.graphDraft.fragments.filter((item) => item.key !== key);
  state.graphDraft.fragments.push({ key, smiles, family, label: key, group: family, roles, notes: "Configured in frontend." });
  setRolesToField(roles);
  state.graphDraft.transitions[key] = state.graphDraft.transitions[key] || [];
  document.getElementById("edgeSource").value = key;
  renderGraphEditor();
}

function addGraphEdge() {
  const source = document.getElementById("edgeSource").value.trim();
  const target = document.getElementById("edgeTarget").value.trim();
  if (!source || !target) {
    alert("起点和终点不能为空。");
    return;
  }
  state.graphDraft.transitions[source] = state.graphDraft.transitions[source] || [];
  if (!state.graphDraft.transitions[source].includes(target)) state.graphDraft.transitions[source].push(target);
  renderGraphEditor();
}

function addEndEdge() {
  const source = document.getElementById("edgeSource").value.trim() || document.getElementById("graphNodeKey").value.trim();
  if (!source) {
    alert("请先填写片段 key 或有向边起点。");
    return;
  }
  document.getElementById("edgeSource").value = source;
  document.getElementById("edgeTarget").value = "End";
  addGraphEdge();
}

function clearGraphDraft() {
  state.graphDraft = { fragments: [], transitions: {} };
  renderGraphEditor();
}

function renderGraphEditor() {
  document.getElementById("graphEditorOutput").textContent = JSON.stringify(state.graphDraft, null, 2);
  renderGraphPreview();
  renderGraphEdgeList();
}

function loadGraphExample() {
  state.graphDraft = {
    fragments: [
      {
        key: "CUSTOM_AR",
        smiles: "[1*]c1ccc([2*])cc1",
        family: "custom",
        label: "aromatic",
        group: "custom",
        roles: ["rigid", "aromatic"],
        notes: "Example rigid aromatic fragment."
      },
      {
        key: "CUSTOM_ETHER",
        smiles: "[1*]CCOCC[2*]",
        family: "custom",
        label: "ether",
        group: "custom",
        roles: ["flexible", "ether"],
        notes: "Example flexible ether fragment."
      },
      {
        key: "CUSTOM_LINK",
        smiles: "[1*]CNC(=O)O[2*]",
        family: "custom",
        label: "linker",
        group: "custom",
        roles: ["polar", "urethane"],
        notes: "Example polar linker."
      }
    ],
    transitions: {
      CUSTOM_AR: ["CUSTOM_ETHER", "CUSTOM_LINK", "End"],
      CUSTOM_ETHER: ["CUSTOM_AR", "End"],
      CUSTOM_LINK: ["CUSTOM_ETHER", "End"]
    }
  };
  document.getElementById("graphNodeKey").value = "CUSTOM_AR";
  document.getElementById("edgeSource").value = "CUSTOM_AR";
  document.getElementById("edgeTarget").value = "CUSTOM_ETHER";
  renderGraphEditor();
}

function renderGraphPreview() {
  const svg = document.getElementById("graphPreviewSvg");
  if (!svg) return;
  while (svg.firstChild) svg.removeChild(svg.firstChild);

  const fragments = state.graphDraft.fragments || [];
  const transitions = state.graphDraft.transitions || {};
  const nodes = buildGraphPreviewNodes(fragments);
  const positions = placeGraphPreviewNodes(nodes);

  addSvgDefs(svg);
  drawGraphEdges(svg, transitions, positions);
  drawGraphNodes(svg, nodes, positions);

  if (fragments.length === 0) {
    drawEmptyGraphHint(svg);
  }
}

function buildGraphPreviewNodes(fragments) {
  return [
    { key: "Start", label: "Start", type: "terminal", roles: ["root"] },
    ...fragments.map((fragment) => ({
      key: fragment.key,
      label: fragment.label || fragment.key,
      type: "fragment",
      roles: fragment.roles || [],
      family: fragment.family || "custom"
    })),
    { key: "End", label: "End", type: "terminal", roles: ["stop"] }
  ];
}

function placeGraphPreviewNodes(nodes) {
  const positions = {};
  const fragments = nodes.filter((node) => node.type === "fragment");
  positions.Start = { x: 54, y: 150 };
  positions.End = { x: 366, y: 150 };

  if (fragments.length === 0) return positions;
  const xStep = fragments.length === 1 ? 0 : 150 / Math.max(1, fragments.length - 1);
  fragments.forEach((node, index) => {
    const rowOffset = index % 2 === 0 ? -46 : 46;
    positions[node.key] = {
      x: 135 + index * xStep,
      y: 150 + rowOffset
    };
  });
  return positions;
}

function addSvgDefs(svg) {
  const defs = createSvg("defs");
  const marker = createSvg("marker", {
    id: "arrowHead",
    markerWidth: 10,
    markerHeight: 10,
    refX: 8,
    refY: 3,
    orient: "auto",
    markerUnits: "strokeWidth"
  });
  const path = createSvg("path", { d: "M0,0 L0,6 L8,3 z", class: "graph-arrow-head" });
  marker.appendChild(path);
  defs.appendChild(marker);
  svg.appendChild(defs);
}

function drawGraphEdges(svg, transitions, positions) {
  const fragmentKeys = Object.keys(positions).filter((key) => key !== "Start" && key !== "End");
  fragmentKeys.forEach((key) => {
    drawEdge(svg, positions.Start, positions[key], "graph-edge start-edge");
  });

  Object.entries(transitions).forEach(([source, targets]) => {
    if (!positions[source] || !Array.isArray(targets)) return;
    targets.forEach((target) => {
      if (!positions[target]) return;
      drawEdge(svg, positions[source], positions[target], target === "End" ? "graph-edge end-edge" : "graph-edge");
    });
  });
}

function drawEdge(svg, from, to, className) {
  const midX = (from.x + to.x) / 2;
  const bend = from.y === to.y ? 0 : from.y < to.y ? 24 : -24;
  const d = `M ${from.x + 39} ${from.y} C ${midX} ${from.y + bend}, ${midX} ${to.y - bend}, ${to.x - 39} ${to.y}`;
  svg.appendChild(createSvg("path", { d, class: className, "marker-end": "url(#arrowHead)" }));
}

function drawGraphNodes(svg, nodes, positions) {
  nodes.forEach((node) => {
    const point = positions[node.key];
    if (!point) return;
    const group = createSvg("g", { class: `graph-node ${node.type === "terminal" ? "terminal-node" : "fragment-node"}` });
    group.appendChild(createSvg("rect", { x: point.x - 43, y: point.y - 28, width: 86, height: 56, rx: 12 }));
    const title = createSvg("text", { x: point.x, y: point.y - 3, "text-anchor": "middle", class: "node-title" });
    title.textContent = shortenGraphLabel(node.key, 13);
    group.appendChild(title);
    const sub = createSvg("text", { x: point.x, y: point.y + 17, "text-anchor": "middle", class: "node-subtitle" });
    sub.textContent = node.type === "terminal" ? node.label : shortenGraphLabel((node.roles || []).join(",") || node.family || "fragment", 15);
    group.appendChild(sub);
    svg.appendChild(group);
  });
}

function drawEmptyGraphHint(svg) {
  const text = createSvg("text", { x: 210, y: 238, "text-anchor": "middle", class: "graph-empty-text" });
  text.textContent = "先加入片段，或点击“载入示例”查看图形效果";
  svg.appendChild(text);
}

function renderGraphEdgeList() {
  const container = document.getElementById("graphEdgeList");
  if (!container) return;
  const entries = Object.entries(state.graphDraft.transitions || {}).filter(([, targets]) => Array.isArray(targets) && targets.length);
  if (entries.length === 0) {
    container.innerHTML = "<span>暂无连接边</span>";
    return;
  }
  container.innerHTML = "";
  entries.forEach(([source, targets]) => {
    targets.forEach((target) => {
      const item = document.createElement("span");
      item.textContent = `${source} -> ${target}`;
      container.appendChild(item);
    });
  });
}

function createSvg(tag, attrs = {}) {
  const item = document.createElementNS("http://www.w3.org/2000/svg", tag);
  Object.entries(attrs).forEach(([key, value]) => item.setAttribute(key, value));
  return item;
}

function shortenGraphLabel(value, maxLength) {
  const text = String(value || "");
  return text.length > maxLength ? `${text.slice(0, maxLength - 1)}...` : text;
}

async function saveGraphDraft() {
  try {
    await apiPost("/api/custom-graph", state.graphDraft);
    document.getElementById("graphPresetSelect").value = "custom";
    await saveConfig();
    alert("自定义有向图已保存。");
  } catch (error) {
    alert(`保存失败：${error.message}`);
  }
}

async function generateByApi() {
  const output = document.getElementById("apiGeneratorOutput");
  const apiKey = document.getElementById("apiKeyInput").value.trim();
  const apiUrl = document.getElementById("apiUrlInput").value.trim();
  const model = document.getElementById("apiModelInput").value.trim();
  const maxFragments = Number(document.getElementById("apiMaxFragmentsInput").value || 8);
  const prompt = document.getElementById("apiPromptInput").value.trim();
  if (!apiKey || !apiUrl || !model) {
    alert("API URL、Model 和 API Key 都需要填写。");
    return;
  }

  output.textContent = "正在请求 API...";
  try {
    const response = await fetch(apiUrl, {
      method: "POST",
      headers: { "Content-Type": "application/json", "Authorization": `Bearer ${apiKey}` },
      body: JSON.stringify({
        model,
        temperature: 0.2,
        messages: [
          { role: "system", content: "Generate polymer fragments for MCTS. Return JSON only. Each fragment SMILES must contain [1*] and [2*]. Every fragment must include a roles array chosen from the provided vocabulary." },
          { role: "user", content: `${prompt}\nReturn at most ${maxFragments} fragments. JSON fields: fragments, transitions. Role vocabulary: ${ROLE_OPTIONS.join(", ")}.` }
        ]
      })
    });
    const data = await response.json();
    const content = data.choices?.[0]?.message?.content || JSON.stringify(data);
    state.apiGenerated = validateGeneratedConfig(parseJsonFromText(content));
    output.textContent = JSON.stringify(state.apiGenerated, null, 2);
  } catch (error) {
    output.textContent = `API 调用失败：${error.message}`;
  }
}

function parseJsonFromText(text) {
  const fenced = text.match(/```(?:json)?\s*([\s\S]*?)\s*```/i);
  const raw = fenced ? fenced[1] : text.slice(text.indexOf("{"), text.lastIndexOf("}") + 1);
  return JSON.parse(raw);
}

function validateGeneratedConfig(data) {
  const fragments = Array.isArray(data.fragments) ? data.fragments : data.custom_fragments || [];
  const transitions = data.transitions || data.custom_transitions || {};
  const cleanedFragments = fragments
    .filter((item) => item && item.key && item.smiles)
    .filter((item) => item.smiles.includes("[1*]") && item.smiles.includes("[2*]"))
    .map((item) => ({
      key: String(item.key).trim(),
      smiles: String(item.smiles).trim(),
      family: String(item.family || "ai_custom").trim(),
      label: String(item.label || item.key).trim(),
      group: String(item.group || "ai").trim(),
      roles: normalizeRoles(item.roles, item.smiles),
      notes: String(item.notes || "Generated by frontend API call.").trim()
    }));
  const validKeys = new Set(cleanedFragments.map((item) => item.key));
  const cleanedTransitions = {};
  Object.entries(transitions).forEach(([source, targets]) => {
    if (!validKeys.has(source) || !Array.isArray(targets)) return;
    cleanedTransitions[source] = targets.filter((target) => target === "End" || validKeys.has(target));
  });
  cleanedFragments.forEach((item) => {
    cleanedTransitions[item.key] = cleanedTransitions[item.key] || ["End"];
  });
  return { fragments: cleanedFragments, transitions: cleanedTransitions };
}

function applyApiGeneratedToGraph() {
  if (!state.apiGenerated || !Array.isArray(state.apiGenerated.fragments)) {
    alert("请先调用 API 生成片段配置。");
    return;
  }
  state.graphDraft = JSON.parse(JSON.stringify(state.apiGenerated));
  document.getElementById("graphPresetSelect").value = "custom";
  const first = state.graphDraft.fragments[0];
  if (first) {
    document.getElementById("graphNodeKey").value = first.key || "";
    document.getElementById("graphNodeSmiles").value = first.smiles || "";
    document.getElementById("graphNodeFamily").value = first.family || "custom";
    setRolesToField(first.roles || []);
    document.getElementById("edgeSource").value = first.key || "";
    document.getElementById("edgeTarget").value = (state.graphDraft.transitions[first.key] || [])[0] || "End";
  }
  renderGraphEditor();
  refreshRoleHint("已把 API 生成的片段和 roles 应用到有向图编辑器。");
}

function getApiJson() {
  return state.apiGenerated ? JSON.stringify(state.apiGenerated, null, 2) : "{}";
}

function downloadText(filename, text) {
  const blob = new Blob([text], { type: "application/json;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  link.click();
  URL.revokeObjectURL(url);
}

function renderCommands() {
  document.getElementById("serverCommand").textContent = [
    "cd LAMMPS_MCTS",
    "conda activate polymer_mcts",
    "python web_server.py",
    "# 打开 http://localhost:8010/frontend/"
  ].join("\n");
  document.getElementById("mainCommand").textContent = [
    "cd LAMMPS_MCTS",
    "conda activate polymer_mcts",
    "python main.py"
  ].join("\n");
  document.getElementById("lammpsCommand").textContent = [
    "# 需要先把 LAMMPS 加入 PATH，或在 config.yaml 中配置 lammps.executable",
    "lmp -in outputs/lammps/candidate_001_rapid_gk.in"
  ].join("\n");
}
