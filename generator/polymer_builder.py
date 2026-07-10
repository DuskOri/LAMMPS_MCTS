"""把 MCTS 输出的片段序列转换成可保存的聚合物结构。"""

from dataclasses import dataclass
from pathlib import Path

try:
    from rdkit import Chem
except ImportError:
    Chem = None

@dataclass
class PolymerBuildResult:
    """记录一次聚合物构建的结果，后续写表和分析都用它。"""

    sequence: list
    dp: int
    smiles: str
    pdb_path: str
    success: bool
    message: str = ""


def normalize_sequence(sequence):
    """整理片段序列，保证末尾带 End。"""
    fragments = [item for item in sequence if item != "Start"]
    if not fragments or fragments[-1] != "End":
        fragments.append("End")
    return fragments


def build_polymer_from_sequence(sequence, dp=2, output_dir="outputs", name="candidate"):
    """根据片段序列生成封端后的聚合物，并保存 PDB 文件。"""
    build_tools = _load_poly_build_tools()
    if Chem is None:
        return PolymerBuildResult(
            sequence=list(sequence),
            dp=dp,
            smiles="",
            pdb_path="",
            success=False,
            message="RDKit is not installed.",
        )

    if build_tools is None:
        return PolymerBuildResult(
            sequence=list(sequence),
            dp=dp,
            smiles="",
            pdb_path="",
            success=False,
            message="Poly_Build module is not available.",
        )

    build_poly_chain, cap_poly_ends, save_pdb = build_tools
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    normalized = normalize_sequence(sequence)
    polymer = build_poly_chain(normalized, dp=dp)
    if polymer is None:
        return PolymerBuildResult(normalized, dp, "", "", False, "polymer build failed")

    capped_polymer = cap_poly_ends(polymer)
    if capped_polymer is None:
        return PolymerBuildResult(normalized, dp, "", "", False, "end capping failed")

    smiles = Chem.MolToSmiles(capped_polymer)
    pdb_file = output_path / f"{name}.pdb"
    saved = save_pdb(capped_polymer, str(pdb_file))

    return PolymerBuildResult(
        sequence=normalized,
        dp=dp,
        smiles=smiles,
        pdb_path=str(pdb_file) if saved else "",
        success=bool(saved),
        message="ok" if saved else "pdb save failed",
    )


def _load_poly_build_tools():
    """延迟导入 Poly_Build，避免 mcts 包初始化时产生循环导入。"""
    try:
        from mcts.Poly_Build import build_poly_chain, cap_poly_ends, save_pdb
        return build_poly_chain, cap_poly_ends, save_pdb
    except ImportError:
        return None
