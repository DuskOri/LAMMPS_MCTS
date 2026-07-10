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
  apiGenerated: null,
  pollTimer: null
};

const stageOrder = ["search", "build", "packmol", "write_lammps", "run_lammps", "analyze", "save"];
const stageText = {
  idle: "等待运行",
  search: "MCTS 搜索",
  build: "RDKit 建链",
  packmol: "Packmol 初始体系",
  write_lammps: "写出 LAMMPS 输入",
  run_lammps: "运行 LAMMPS",
  analyze: "热导率后处理",
  save: "保存数据库",
  done: "完成",
  failed: "失败"
};

document.addEventListener("DOMContentLoaded", () => {
  bindTabs();
  bindActions();
  loadInitialData();
});

function bindTabs() {
  document.querySelectorAll(".nav-button").forEach((button) => {
    button.addEventListener("click", () => {
      const tab = button.dataset.tab;
      document.querySelectorAll(".nav-button").forEach((item) => item.classList.remove("active"));
      document.querySelectorAll(".tab-panel").forEach((item) => item.classList.remove("active"));
      button.classList.add("active");
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

  document.getElementById("addGraphNodeBtn").addEventListener("click", addGraphNode);
  document.getElementById("loadGraphExampleBtn").addEventListener("click", loadGraphExample);
  document.getElementById("clearGraphBtn").addEventListener("click", clearGraphDraft);
  document.getElementById("addEdgeBtn").addEventListener("click", addGraphEdge);
  document.getElementById("addEndEdgeBtn").addEventListener("click", addEndEdge);
  document.getElementById("saveGraphBtn").addEventListener("click", saveGraphDraft);
  document.getElementById("generateByApiBtn").addEventListener("click", generateByApi);
  document.getElementById("downloadApiJsonBtn").addEventListener("click", () => downloadText("ai_fragments.json", getApiJson()));
}

async function loadInitialData() {
  await loadConfig();
  await loadResults();
  await refreshStatus();
  renderCommands();
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
  document.getElementById("dpInput").value = config.dp || 10;
  document.getElementById("iterationsInput").value = config.iterations || 300;
  document.getElementById("topKInput").value = config.top_k || 5;
  document.getElementById("maxStepsInput").value = config.max_steps || 6;
  document.getElementById("graphPresetSelect").value = config.active_families || "polyimide,polyurethane,phenolic";
  document.getElementById("runLammpsInput").checked = Boolean(config.run_lammps);
  document.getElementById("thermalFeedbackInput").checked = config.feedback_mode === "thermal";
}

function readConfigForm() {
  return {
    dp: Number(document.getElementById("dpInput").value || 10),
    iterations: Number(document.getElementById("iterationsInput").value || 300),
    top_k: Number(document.getElementById("topKInput").value || 5),
    max_steps: Number(document.getElementById("maxStepsInput").value || 6),
    active_families: document.getElementById("graphPresetSelect").value,
    run_lammps: document.getElementById("runLammpsInput").checked,
    feedback_mode: document.getElementById("thermalFeedbackInput").checked ? "thermal" : "heuristic"
  };
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
  showStatus(status.stage, status.message || "");
  if (status.counts) {
    document.getElementById("candidateCount").textContent = status.counts.candidates || 0;
    document.getElementById("buildCount").textContent = status.counts.build_results || 0;
    document.getElementById("inputCount").textContent = status.counts.thermal_inputs || 0;
  }
}

function showStatus(stage, message) {
  document.getElementById("statusTitle").textContent = stageText[stage] || stage;
  document.getElementById("statusMessage").textContent = message || "";
  renderStageList(stage);
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
    "heat_flux",
    "temperature_gradient",
    "success",
    "message"
  ]);
  renderCandidateTable();
  renderTable("thermalTable", state.results.thermal_results || [], [
    "low_k_rank",
    "input_path",
    "conductivity_w_mk",
    "heat_flux",
    "temperature_gradient",
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
  renderTable("candidateTable", rows, [
    "rank",
    "sequence_text",
    "score",
    "visits",
    "length",
    "family_text",
    "flexible_ratio",
    "aromatic_ratio"
  ]);
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

function addGraphNode() {
  const key = document.getElementById("graphNodeKey").value.trim();
  const smiles = document.getElementById("graphNodeSmiles").value.trim();
  const family = document.getElementById("graphNodeFamily").value.trim() || "custom";
  const roles = document.getElementById("graphNodeRoles").value.split(",").map((item) => item.trim()).filter(Boolean);
  if (!key || !smiles) {
    alert("片段 key 和 SMILES 不能为空。");
    return;
  }
  if (!smiles.includes("[1*]") || !smiles.includes("[2*]")) {
    alert("SMILES 必须包含 [1*] 和 [2*]。");
    return;
  }
  state.graphDraft.fragments = state.graphDraft.fragments.filter((item) => item.key !== key);
  state.graphDraft.fragments.push({ key, smiles, family, label: key, group: family, roles, notes: "Configured in frontend." });
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
          { role: "system", content: "Generate polymer fragments for MCTS. Return JSON only. Each fragment SMILES must contain [1*] and [2*]." },
          { role: "user", content: `${prompt}\nReturn at most ${maxFragments} fragments. JSON fields: fragments, transitions.` }
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
      roles: Array.isArray(item.roles) ? item.roles : String(item.roles || "").split(",").map((role) => role.trim()).filter(Boolean),
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
  document.getElementById("serverCommand").textContent = "cd C:\\DuskORI\\Files\\Code\\LAMMPS_MCTS\\LAMMPS_MCTS\nC:\\DuskORI\\Application\\miniconda\\envs\\polymer_mcts\\python.exe web_server.py\n# 打开 http://localhost:8010/frontend/";
  document.getElementById("mainCommand").textContent = "cd C:\\DuskORI\\Files\\Code\\LAMMPS_MCTS\\LAMMPS_MCTS\nC:\\DuskORI\\Application\\miniconda\\envs\\polymer_mcts\\python.exe main.py";
  document.getElementById("lammpsCommand").textContent = "C:\\DuskORI\\Application\\LAMMPS\\bin\\lmp.exe -in outputs\\lammps\\candidate_001_fast_tc.in";
}
