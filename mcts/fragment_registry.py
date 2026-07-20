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
            custom_rows = []
        else:
            custom_rows = custom_fragments or ()

        # 前端把内置家族与自定义片段作为互斥搜索空间。只有明确选择
        # custom 时才加载自定义文件，避免 CUSTOM_* 混入聚酰亚胺等内置家族。
        if "custom" not in family_set:
            custom_rows = ()
            custom_transitions = {}

        if not isinstance(custom_transitions, Mapping):
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
            key = str(row.get("key", "")).strip()
            smiles = str(row.get("smiles", "")).strip()
            if not key or not smiles:
                continue

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
            if source_key not in self._fragments:
                continue

            targets = []
            for target in raw_targets:
                target_key = str(target).strip()
                if target_key == "End" or target_key in self._fragments:
                    targets.append(target_key)
            if targets:
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


# 片段维护备注
# 1. 这里的片段不是完整单体，而是已经带有两个连接端点的链段。
# 2. [1*] 和 [2*] 只用于拼接定位，最终建模前会被真实端基替换。
# 3. 如果新增片段只有一个连接点，MCTS 可以搜到，但 RDKit 建链会失败。
# 4. 如果片段有三个以上连接点，本项目当前仍按线型聚合物处理，不会自动支化。
# 5. 聚酰亚胺片段保留 A、A'、B、C，是为了兼容最早的搜索空间。
# 6. A 和 A' 都属于酰亚胺刚性片段，热导启发式里按 rigid/imide 处理。
# 7. B 是对苯撑桥，主要贡献芳香刚性和链段线性。
# 8. C 是醚键片段，主要用于引入柔性连接。
# 9. 聚氨酯片段按“已形成氨基甲酸酯键的链段”来写。
# 10. 这样可以避开显式反应机理，只检查最终链段能否拼接。
# 11. PU_MDI 表示芳香硬段，通常会提高链段刚性。
# 12. PU_HDI 表示脂肪族硬段，刚性弱于 PU_MDI。
# 13. PU_BDO 表示短二醇软段，可增加构象自由度。
# 14. PU_PEG 表示聚醚软段，柔性更强，也带来更多杂原子。
# 15. 酚醛树脂片段按酚环和亚甲基桥来拆分。
# 16. PF_PHENOL 是普通酚环片段。
# 17. PF_RESOL 带羟甲基侧基，用于表示 resol 型结构倾向。
# 18. PF_METHYLENE 是酚醛网络中常见的桥接单元。
# 19. 当前 Packmol 和 LAMMPS data 写出仍按线型分子处理。
# 20. 对真实热固性酚醛网络，后续可以单独增加交联建模模块。
# 21. active_families 控制内置片段是否进入 Start 的候选集。
# 22. 如果只研究聚氨酯，可以把 active_families 写成 polyurethane。
# 23. 如果研究多体系对比，可以保留多个 family 并让 MCTS 混合探索。
# 24. 混合探索能扩大空间，但化学可解释性需要人工复核。
# 25. custom_fragments.json 适合放课题中的候选片段。
# 26. 自定义片段 key 不要和内置 key 重名，除非明确要覆盖内置定义。
# 27. 自定义片段 smiles 必须能被 RDKit 的 Chem.MolFromSmiles 解析。
# 28. 自定义片段最好先单独做一次重复单元拼接检查。
# 29. 如果小规模拼接失败，通常是虚拟端点或价态写法有问题。
# 30. 如果 UFF 优化失败但 PDB 已生成，可以先检查初始三维构象。
# 31. custom_transitions 不写时，自定义片段默认互相可接并可 End。
# 32. 正式筛选建议写 custom_transitions，避免 MCTS 搜索无意义组合。
# 33. transitions 里的 End 表示该路径可以终止。
# 34. 达到 max_steps 后状态机会强制补上 End，结束当前重复单元。
# 35. max_steps 只计算真实片段，不计算 Start 和 End，也不是聚合度。
# 36. 聚合度由重复单元轮廓长度和目标链长自动计算。
# 37. 不同重复单元长度不同，因此自动得到的实际 DP 也可以不同。
# 38. build_results.csv 会记录实际 DP 和估算链长，便于人工复核。
# 39. 低热导目标使用 reward=1/(1+k/reward_scale) 的方向。
# 40. 因此 MCTS 分数越高，代表估计或反馈的热导率越低。
# 41. 启发式奖励只负责快速导向，不能替代最终 LAMMPS 结果。
# 42. thermal feedback 模式会把 LAMMPS 后处理结果反馈到节点。
# 43. feedback_max_evaluations 用于限制昂贵模拟次数。
# 44. 如果 LAMMPS 失败且允许 fallback，会退回启发式评分。
# 45. 如果不允许 fallback，失败序列会得到 failure_reward。
# 46. 片段 roles 会影响启发式热导估计。
# 47. rigid/aromatic 通常增加估计热导。
# 48. flexible/ether 通常降低估计热导。
# 49. imide/urethane/phenolic 只作为化学类别和轻微修正项使用。
# 50. 这些权重只是筛选前的粗略假设，后续可由数据拟合更新。
# 51. CSV 输出会同时记录 sequence_text 和 family_text。
# 52. sequence_text 便于人工直接看片段顺序。
# 53. family_text 便于区分 PI、PU、PF 或 custom 来源。
# 54. 如果要申请软著，建议保留这些注释和配置示例。
# 55. 注释反映的是工程约束、数据格式和使用边界。
# 56. 它们不是装饰，后续调试自定义片段时会经常用到。
# 57. 新增片段后优先跑语法检查，再跑一个很小的 MCTS 搜索。
# 58. 小搜索通过后，再放大 top_k、DP 和 Packmol 分子数。
# 59. 如果 Packmol 长时间不收敛，先减小 molecule_count 或增大 box_size。
# 60. 如果 LAMMPS 报键跨越盒子，先检查 PDB 初始构象和盒子尺寸。
# 61. 如果需要真实热导率，后续还要补可靠力场和充分弛豫。
# 62. 当前快速脚本主要用于流程验证和低热导候选的第一轮筛选。
# 63. 未参数化的 data 会被 Green-Kubo 输入生成器拒绝，不能进入 thermal reward。
# 64. 自定义片段的命名建议使用体系前缀，如 PU_、PF_、CUSTOM_。
# 65. 命名清楚以后，结果表里的序列更容易追溯。
# 66. 片段库扩展时不要直接删除旧 key，历史输出可能还会引用它们。
# 67. 如果必须重命名，可以在 README 或变更记录里写清楚映射关系。
# 68. 这个模块只保存化学空间，不负责三维构象和模拟参数。
# 69. 三维构象由 Poly_Build 和 RDKit 负责。
# 70. 初始体系装盒由 packmol_runner 负责。
# 71. LAMMPS data 写出由 lammps_data_writer 负责。
# 72. 热导率输入脚本由 thermal_conductivity_writer 负责。
# 73. 这种分层能让片段库独立演化，减少互相牵连。
