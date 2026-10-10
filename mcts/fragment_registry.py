"""聚合物片段库和连接规则。

本模块集中保存 MCTS 搜索可以使用的片段。这样做有两个目的：
一是避免片段 SMILES、片段分类、片段连接图分散在多个文件里；
二是给自定义化合物输入留出统一入口，后面新增片段时不用改搜索算法。

片段 SMILES 里保留两个虚拟连接点：
    [1*] 表示片段左端连接位点
    [2*] 表示片段右端连接位点

Poly_Build 里的拼接反应会把前一个片段的 [2*] 和后一个片段的 [1*]
去掉，然后在两个相邻真实原子之间建立单键。自定义片段也遵守同样规则。
"""

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence


@dataclass(frozen=True)
class FragmentSpec:
    """描述一个可被 MCTS 选择、可被 RDKit 拼接的聚合物片段。"""

    key: str
    smiles: str
    family: str
    label: str
    group: str = "middle"
    notes: str = ""
    roles: tuple[str, ...] = field(default_factory=tuple)


BUILTIN_FRAGMENTS: tuple[FragmentSpec, ...] = (
    # 聚酰亚胺片段。A/A' 提供酰亚胺结构，B 提供芳香刚性链段，C 提供醚键柔性链段。
    FragmentSpec(
        key="A",
        smiles="O=C1c2ccc([1*])cc2C(=O)N1[2*]",
        family="polyimide",
        label="phthalimide",
        group="imide",
        notes="单酰亚胺型刚性片段，用于保持最早版本 PI 搜索空间。",
        roles=("rigid", "imide"),
    ),
    FragmentSpec(
        key="A'",
        smiles="O=C1c2cc3C(=O)N([2*])C(=O)c3cc2C(=O)N1[1*]",
        family="polyimide",
        label="pyromellitic diimide",
        group="imide",
        notes="均苯四甲酸二酰亚胺结构，提供更高芳香环和酰亚胺比例。",
        roles=("rigid", "imide"),
    ),
    FragmentSpec(
        key="B",
        smiles="[1*]c1ccc([2*])cc1",
        family="polyimide",
        label="p-phenylene",
        group="aromatic",
        notes="对苯撑片段，作为 PI 主链里的芳香桥接单元。",
        roles=("rigid", "aromatic"),
    ),
    FragmentSpec(
        key="C",
        smiles="[1*]O[2*]",
        family="polyimide",
        label="ether",
        group="flexible",
        notes="醚键柔性片段，通常用于降低链段刚性和热传导连续性。",
        roles=("flexible", "ether"),
    ),
    # 聚氨酯片段。这里把二异氰酸酯和二醇反应后的氨基甲酸酯结构写入片段，
    # 建链时仍使用统一的虚拟端点拼接，不在 MCTS 阶段显式模拟加成反应。
    FragmentSpec(
        key="PU_MDI",
        smiles="[1*]NC(=O)Oc1ccc(Cc2ccc(OC(=O)N[2*])cc2)cc1",
        family="polyurethane",
        label="MDI urethane hard segment",
        group="urethane_hard",
        notes="含二苯甲烷骨架的聚氨酯硬段，芳香性较强。",
        roles=("rigid", "urethane", "aromatic"),
    ),
    FragmentSpec(
        key="PU_HDI",
        smiles="[1*]NC(=O)OCCCCCCOC(=O)N[2*]",
        family="polyurethane",
        label="HDI urethane segment",
        group="urethane_hard",
        notes="脂肪族聚氨酯片段，比 MDI 片段更柔顺。",
        roles=("urethane", "aliphatic"),
    ),
    FragmentSpec(
        key="PU_BDO",
        smiles="[1*]CCCC[2*]",
        family="polyurethane",
        label="butanediol soft segment",
        group="urethane_soft",
        notes="丁二醇碳链软段，端点用碳连接以避免 N-O 非典型连接。",
        roles=("flexible", "diol"),
    ),
    FragmentSpec(
        key="PU_PEG",
        smiles="[1*]CCOCC[2*]",
        family="polyurethane",
        label="PEG soft segment",
        group="urethane_soft",
        notes="聚醚型软段，保留内部醚键，端点用碳连接。",
        roles=("flexible", "ether", "diol"),
    ),
    # 酚醛树脂片段。酚环和亚甲基桥是最小可解释单元，适合做快速搜索。
    FragmentSpec(
        key="PF_PHENOL",
        smiles="[1*]c1cc(O)cc([2*])c1",
        family="phenolic",
        label="phenol ring",
        group="phenol",
        notes="酚环片段，可与亚甲基桥交替连接。",
        roles=("rigid", "phenolic", "aromatic"),
    ),
    FragmentSpec(
        key="PF_METHYLENE",
        smiles="[1*]C[2*]",
        family="phenolic",
        label="methylene bridge",
        group="bridge",
        notes="酚醛树脂中的亚甲基桥接单元。",
        roles=("flexible", "bridge"),
    ),
    FragmentSpec(
        key="PF_RESOL",
        smiles="[1*]c1c(O)cc(CO)cc1[2*]",
        family="phenolic",
        label="methylol phenol",
        group="phenol",
        notes="带羟甲基侧基的酚环，用于表示 resol 型酚醛片段。",
        roles=("phenolic", "aromatic"),
    ),
)


BUILTIN_TRANSITIONS: dict[str, list[str]] = {
    # 聚酰亚胺：保持原始 PI 搜索空间，同时允许柔性醚键回到芳香桥。
    "A": ["B"],
    "A'": ["B"],
    "B": ["A", "A'", "B", "C", "End"],
    "C": ["B"],
    # 聚氨酯：硬段和软段交替，硬段之间也允许少量连接以扩大搜索空间。
    "PU_MDI": ["PU_BDO", "PU_PEG", "End"],
    "PU_HDI": ["PU_BDO", "PU_PEG", "End"],
    "PU_BDO": ["PU_MDI", "PU_HDI", "End"],
    "PU_PEG": ["PU_MDI", "PU_HDI", "End"],
    # 酚醛树脂：酚环与桥接单元交替，保留 resol 片段作为酚环变体。
    "PF_PHENOL": ["PF_METHYLENE", "PF_RESOL", "End"],
    "PF_RESOL": ["PF_METHYLENE", "PF_PHENOL", "End"],
    "PF_METHYLENE": ["PF_PHENOL", "PF_RESOL"],
}


BUILTIN_STARTS: dict[str, list[str]] = {
    "polyimide": ["A", "A'"],
    "polyurethane": ["PU_MDI", "PU_HDI"],
    "phenolic": ["PF_PHENOL", "PF_RESOL"],
}


class FragmentRegistry:
    """保存当前运行使用的片段集合和搜索连接图。"""

    def __init__(self):
        self._fragments: MutableMapping[str, FragmentSpec] = {}
        self._transitions: MutableMapping[str, list[str]] = {}
        self._starts_by_family: MutableMapping[str, list[str]] = {}
        self.configure()

    def configure(
        self,
        active_families: Sequence[str] | None = None,
        custom_fragments: Sequence[Mapping[str, object]] | None = None,
        custom_transitions: Mapping[str, Sequence[str]] | None = None,
        allow_custom_self_transitions: bool = True,
    ) -> None:
        """根据配置重建片段库。

        active_families 控制启用哪些内置聚合物家族。为空时默认全部启用。
        custom_fragments 可以从 config.yaml 传入，至少需要 key 和 smiles。
        custom_transitions 可以显式写自定义片段之间的连接关系。
        """
        if isinstance(active_families, str):
            families = tuple(
                item.strip() for item in active_families.split(",") if item.strip()
            )
        else:
            families = tuple(active_families or ("polyimide", "polyurethane", "phenolic"))
        family_set = set(families)

        self._fragments = {}
        self._transitions = {"Start": [], "End": []}
        self._starts_by_family = {}

        for spec in BUILTIN_FRAGMENTS:
            if spec.family not in family_set:
                continue
            self._fragments[spec.key] = spec

        for family in families:
            starts = [key for key in BUILTIN_STARTS.get(family, []) if key in self._fragments]
            if starts:
                self._starts_by_family[family] = starts
                self._transitions["Start"].extend(starts)

        for source, targets in BUILTIN_TRANSITIONS.items():
            if source not in self._fragments:
                continue
            self._transitions[source] = [
                target for target in targets if target == "End" or target in self._fragments
            ]

        if isinstance(custom_fragments, Mapping):
            custom_rows = [custom_fragments]
        elif isinstance(custom_fragments, str):
            if custom_fragments.strip():
                raise ValueError("custom_fragments must be a list of objects")
            custom_rows = []
        else:
            custom_rows = custom_fragments or ()

        # 前端把内置家族与自定义片段作为互斥搜索空间。只有明确选择
        # custom 时才加载自定义文件，避免 CUSTOM_* 混入聚酰亚胺等内置家族。
        if "custom" not in family_set:
            custom_rows = ()
            custom_transitions = {}

        if not isinstance(custom_transitions, Mapping):
            if custom_transitions is not None:
                raise ValueError("custom_transitions must be an object")
            custom_transition_rows = {}
        else:
            custom_transition_rows = custom_transitions

        custom_keys = self._load_custom_fragments(custom_rows)
        self._load_custom_transitions(
            custom_keys=custom_keys,
            custom_transitions=custom_transition_rows,
            allow_self=allow_custom_self_transitions,
        )

    def _load_custom_fragments(self, rows: Iterable[Mapping[str, object]]) -> list[str]:
        """把用户在 config.yaml 里写的片段加入库中。"""
        custom_keys: list[str] = []
        for row in rows:
            if not isinstance(row, Mapping):
                raise ValueError("custom fragments must be objects")
            key = str(row.get("key", "")).strip()
            smiles = str(row.get("smiles", "")).strip()
            if not key or not smiles or key in ("Start", "End") or key in self._fragments:
                raise ValueError(f"invalid or duplicate fragment key: {key!r}")
            from rdkit import Chem
            molecule = Chem.MolFromSmiles(smiles)
            if molecule is None or len(Chem.GetMolFrags(molecule)) != 1:
                raise ValueError(f"fragment {key}: invalid or disconnected SMILES")
            dummies = [atom for atom in molecule.GetAtoms() if atom.GetAtomicNum() == 0]
            if sorted(atom.GetIsotope() for atom in dummies) != [1, 2] or any(atom.GetDegree() != 1 for atom in dummies):
                raise ValueError(f"fragment {key}: exactly one terminal [1*] and [2*] are required")

            family = str(row.get("family", "custom")).strip() or "custom"
            label = str(row.get("label", key)).strip() or key
            group = str(row.get("group", "custom")).strip() or "custom"
            notes = str(row.get("notes", "")).strip()

            roles_raw = row.get("roles", ())
            if isinstance(roles_raw, str):
                roles = tuple(item.strip() for item in roles_raw.split(",") if item.strip())
            else:
                roles = tuple(str(item).strip() for item in roles_raw or () if str(item).strip())

            self._fragments[key] = FragmentSpec(
                key=key,
                smiles=smiles,
                family=family,
                label=label,
                group=group,
                notes=notes,
                roles=roles,
            )
            self._starts_by_family.setdefault(family, []).append(key)
            self._transitions.setdefault("Start", []).append(key)
            custom_keys.append(key)
        return custom_keys

    def _load_custom_transitions(
        self,
        custom_keys: Sequence[str],
        custom_transitions: Mapping[str, Sequence[str]],
        allow_self: bool,
    ) -> None:
        """整理自定义片段的连接图。

        如果用户没有给出 custom_transitions，默认允许自定义片段之间互相连接并可结束。
        这个默认规则比较宽，适合先验证片段能不能被 RDKit 正常拼接。
        做正式课题筛选时，可以在配置文件里写更严格的连接图。
        """
        for source, raw_targets in custom_transitions.items():
            source_key = str(source).strip()
            if source_key != "Start" and source_key not in self._fragments:
                raise ValueError(f"unknown transition source: {source_key}")
            if not isinstance(raw_targets, (list, tuple)):
                raise ValueError(f"transition targets for {source_key} must be a list")

            targets = []
            for target in raw_targets:
                target_key = str(target).strip()
                if target_key == "End" or target_key in self._fragments:
                    if target_key not in targets:
                        targets.append(target_key)
                else:
                    raise ValueError(f"unknown transition target: {target_key}")
            self._transitions[source_key] = targets

        if not custom_keys:
            return

        for source in custom_keys:
            if source in self._transitions:
                continue
            targets = list(custom_keys)
            if not allow_self:
                targets = [target for target in targets if target != source]
            targets.append("End")
            self._transitions[source] = targets

    def smiles_map(self) -> dict[str, str]:
        """返回 key -> SMILES 的普通字典，兼容旧代码的使用方式。"""
        return {key: spec.smiles for key, spec in self._fragments.items()}

    def graph(self) -> dict[str, list[str]]:
        """返回当前 MCTS 使用的连接图。"""
        graph = {key: list(value) for key, value in self._transitions.items()}
        graph.setdefault("Start", [])
        graph.setdefault("End", [])
        return graph

    def allowed_fragments(self, fragment: str) -> list[str]:
        """给出某个片段后面允许接入的候选片段。"""
        return list(self._transitions.get(fragment, []))

    def has_fragment(self, key: str) -> bool:
        """判断片段是否已经注册。"""
        return key in self._fragments

    def metadata(self, key: str) -> FragmentSpec:
        """取得片段的完整说明信息。"""
        return self._fragments[key]

    def keys(self) -> list[str]:
        """返回所有真实片段，不包含 Start/End。"""
        return list(self._fragments.keys())

    def families(self) -> list[str]:
        """返回当前启用的聚合物家族名称。"""
        return sorted({spec.family for spec in self._fragments.values()})


DEFAULT_REGISTRY = FragmentRegistry()


def configure_fragment_registry(config: Mapping[str, object] | None = None) -> FragmentRegistry:
    """根据配置文件刷新全局片段库。"""
    config = config or {}
    DEFAULT_REGISTRY.configure(
        active_families=config.get("active_families"),
        custom_fragments=config.get("custom_fragments"),
        custom_transitions=config.get("custom_transitions"),
        allow_custom_self_transitions=bool(config.get("allow_custom_self_transitions", True)),
    )
    return DEFAULT_REGISTRY


def get_fragment_registry() -> FragmentRegistry:
    """返回当前进程使用的全局片段库。"""
    return DEFAULT_REGISTRY


# 片段按线型重复单元处理，每个片段必须保留 [1*]、[2*] 两个端点。
# 新片段先用 RDKit 单独检查 SMILES 和一次拼接，再放进搜索图。
# PU 和 PF 目前也是线型近似；交联结构留给后续单独处理。
# roles 只影响启发式奖励，使用 LAMMPS 反馈时不读取这些标签。
# 历史结果会保存片段 key，已有 key 不轻易改名。
# 片段较多时应明确写 transitions，避免搜索到明显无意义的组合。
