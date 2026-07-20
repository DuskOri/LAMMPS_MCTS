"""把 MCTS 输出的片段序列转换成可保存的聚合物结构。"""

from dataclasses import dataclass
import math
from pathlib import Path

try:
    from rdkit import Chem
    from rdkit.Chem import rdForceFieldHelpers
except ImportError:
    Chem = None
    rdForceFieldHelpers = None


@dataclass(frozen=True)
class RepeatUnitLengthEstimate:
    """保存重复单元沿化学键路径的长度估算。"""

    internal_length_angstrom: float
    connection_length_angstrom: float

    @property
    def repeat_increment_angstrom(self):
        """返回每增加一个重复单元带来的轮廓长度。"""
        return self.internal_length_angstrom + self.connection_length_angstrom

    def chain_length(self, dp):
        """估算指定聚合度链的轮廓长度。"""
        count = max(1, int(dp))
        return (
            count * self.internal_length_angstrom
            + (count - 1) * self.connection_length_angstrom
        )

@dataclass
class PolymerBuildResult:
    """记录一次聚合物构建的结果，后续写表和分析都用它。"""

    sequence: list
    dp: int
    smiles: str
    pdb_path: str
    success: bool
    message: str = ""
    mol_path: str = ""
    target_length_angstrom: float = 0.0
    repeat_unit_length_angstrom: float = 0.0
    estimated_chain_length_angstrom: float = 0.0
    dp_mode: str = "auto"


def normalize_sequence(sequence):
    """整理片段序列，保证末尾带 End。"""
    fragments = [item for item in sequence if item != "Start"]
    if not fragments or fragments[-1] != "End":
        fragments.append("End")
    return fragments


def build_polymer_from_sequence(
    sequence,
    dp=None,
    target_length_angstrom=120.0,
    min_dp=1,
    max_dp=100,
    output_dir="outputs",
    name="candidate",
):
    """构建接近目标轮廓长度的聚合物，并保存 MOL 和 PDB。"""
    build_tools = _load_poly_build_tools()
    if Chem is None:
        return PolymerBuildResult(
            sequence=list(sequence),
            dp=int(dp or 0),
            smiles="",
            pdb_path="",
            success=False,
            message="RDKit is not installed.",
            target_length_angstrom=float(target_length_angstrom),
        )

    if build_tools is None:
        return PolymerBuildResult(
            sequence=list(sequence),
            dp=int(dp or 0),
            smiles="",
            pdb_path="",
            success=False,
            message="Poly_Build module is not available.",
            target_length_angstrom=float(target_length_angstrom),
        )

    build_poly_chain, cap_poly_ends, save_pdb = build_tools
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    normalized = normalize_sequence(sequence)
    repeat_unit = build_poly_chain(normalized, dp=1)
    if repeat_unit is None:
        return PolymerBuildResult(
            sequence=normalized,
            dp=int(dp or 0),
            smiles="",
            pdb_path="",
            success=False,
            message="repeat unit build failed",
            target_length_angstrom=float(target_length_angstrom),
        )

    try:
        length_estimate = estimate_repeat_unit_contour_length(repeat_unit)
    except ValueError as error:
        return PolymerBuildResult(
            sequence=normalized,
            dp=int(dp or 0),
            smiles="",
            pdb_path="",
            success=False,
            message=str(error),
            target_length_angstrom=float(target_length_angstrom),
        )

    if dp is None:
        resolved_dp = resolve_degree_of_polymerization(
            length_estimate,
            target_length_angstrom=target_length_angstrom,
            min_dp=min_dp,
            max_dp=max_dp,
        )
        dp_mode = "auto"
    else:
        resolved_dp = max(int(min_dp), min(int(max_dp), int(dp)))
        dp_mode = "manual"

    polymer = repeat_unit if resolved_dp == 1 else build_poly_chain(normalized, dp=resolved_dp)
    if polymer is None:
        return PolymerBuildResult(
            sequence=normalized,
            dp=resolved_dp,
            smiles="",
            pdb_path="",
            success=False,
            message="polymer build failed",
            target_length_angstrom=float(target_length_angstrom),
            repeat_unit_length_angstrom=length_estimate.repeat_increment_angstrom,
            estimated_chain_length_angstrom=length_estimate.chain_length(resolved_dp),
            dp_mode=dp_mode,
        )

    capped_polymer = cap_poly_ends(polymer)
    if capped_polymer is None:
        return PolymerBuildResult(
            sequence=normalized,
            dp=resolved_dp,
            smiles="",
            pdb_path="",
            success=False,
            message="end capping failed",
            target_length_angstrom=float(target_length_angstrom),
            repeat_unit_length_angstrom=length_estimate.repeat_increment_angstrom,
            estimated_chain_length_angstrom=length_estimate.chain_length(resolved_dp),
            dp_mode=dp_mode,
        )

    smiles = Chem.MolToSmiles(capped_polymer)
    mol_file = output_path / f"{name}.mol"
    pdb_file = output_path / f"{name}.pdb"
    Chem.MolToMolFile(capped_polymer, str(mol_file))
    saved = save_pdb(capped_polymer, str(pdb_file))

    return PolymerBuildResult(
        sequence=normalized,
        dp=resolved_dp,
        smiles=smiles,
        pdb_path=str(pdb_file) if saved else "",
        success=bool(saved),
        message="ok" if saved else "pdb save failed",
        mol_path=str(mol_file) if saved else "",
        target_length_angstrom=float(target_length_angstrom),
        repeat_unit_length_angstrom=length_estimate.repeat_increment_angstrom,
        estimated_chain_length_angstrom=length_estimate.chain_length(resolved_dp),
        dp_mode=dp_mode,
    )


def estimate_repeat_unit_contour_length(molecule):
    """沿两个连接位点之间的最短化学键路径估算重复单元长度。"""
    if molecule is None:
        raise ValueError("repeat unit molecule is empty")

    terminals = {}
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() != 0 or atom.GetIsotope() not in (1, 2):
            continue
        neighbors = list(atom.GetNeighbors())
        if len(neighbors) != 1:
            raise ValueError("repeat unit connection site must have one neighbor")
        terminals[atom.GetIsotope()] = neighbors[0].GetIdx()

    if set(terminals) != {1, 2}:
        raise ValueError("repeat unit must contain [1*] and [2*] connection sites")

    left_index = terminals[1]
    right_index = terminals[2]
    path = Chem.rdmolops.GetShortestPath(molecule, left_index, right_index)
    if len(path) < 2:
        raise ValueError("repeat unit connection sites do not form a valid path")

    internal_length = sum(
        _equilibrium_bond_length(molecule, atom_i, atom_j)
        for atom_i, atom_j in zip(path, path[1:])
    )
    connection_length = _covalent_bond_length(
        molecule.GetAtomWithIdx(right_index),
        molecule.GetAtomWithIdx(left_index),
    )
    return RepeatUnitLengthEstimate(
        internal_length_angstrom=internal_length,
        connection_length_angstrom=connection_length,
    )


def resolve_degree_of_polymerization(
    length_estimate,
    target_length_angstrom=120.0,
    min_dp=1,
    max_dp=100,
):
    """选择使链轮廓长度最接近目标值的整数聚合度。"""
    target = float(target_length_angstrom)
    minimum = max(1, int(min_dp))
    maximum = max(minimum, int(max_dp))
    increment = float(length_estimate.repeat_increment_angstrom)
    if target <= 0.0:
        raise ValueError("target chain length must be positive")
    if not math.isfinite(increment) or increment <= 0.0:
        raise ValueError("repeat unit contour length must be positive")

    approximate = (target + length_estimate.connection_length_angstrom) / increment
    center = int(math.floor(approximate + 0.5))
    candidates = {
        max(minimum, min(maximum, center - 1)),
        max(minimum, min(maximum, center)),
        max(minimum, min(maximum, center + 1)),
    }
    return min(
        candidates,
        key=lambda value: (abs(length_estimate.chain_length(value) - target), value),
    )


def _equilibrium_bond_length(molecule, atom_i, atom_j):
    if rdForceFieldHelpers is not None:
        values = rdForceFieldHelpers.GetUFFBondStretchParams(molecule, atom_i, atom_j)
        if values is not None:
            return float(values[1])
    return _covalent_bond_length(
        molecule.GetAtomWithIdx(atom_i),
        molecule.GetAtomWithIdx(atom_j),
    )


def _covalent_bond_length(atom_i, atom_j):
    periodic_table = Chem.GetPeriodicTable()
    return float(
        periodic_table.GetRcovalent(atom_i.GetAtomicNum())
        + periodic_table.GetRcovalent(atom_j.GetAtomicNum())
    )


def _load_poly_build_tools():
    """延迟导入 Poly_Build，避免 mcts 包初始化时产生循环导入。"""
    try:
        from mcts.Poly_Build import build_poly_chain, cap_poly_ends, save_pdb
        return build_poly_chain, cap_poly_ends, save_pdb
    except ImportError:
        return None
