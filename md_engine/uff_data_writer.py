"""用 RDKit UFF 参数生成筛选级 LAMMPS data。"""

from itertools import combinations
from pathlib import Path
import math

from rdkit import Chem
from rdkit.Chem import AllChem, rdForceFieldHelpers

from .lammps_data_writer import ATOMIC_MASSES, read_pdb_topology


def write_uff_lammps_data_from_template(
    template_mol,
    coordinate_pdb,
    filename,
    molecule_count=1,
    box_size=60.0,
):
    """把单链 MOL 的 UFF 参数复制到 Packmol 多链坐标。"""
    molecule_count = int(molecule_count)
    molecule = Chem.MolFromMolFile(str(template_mol), removeHs=False, sanitize=True)
    if molecule is None:
        raise ValueError("template MOL cannot be parsed by RDKit")
    if not rdForceFieldHelpers.UFFHasAllMoleculeParams(molecule):
        raise ValueError("RDKit UFF does not cover every atom in the polymer")

    AllChem.ComputeGasteigerCharges(molecule)
    packed = read_pdb_topology(coordinate_pdb)
    atoms_per_molecule = molecule.GetNumAtoms()
    expected_atoms = atoms_per_molecule * molecule_count
    if len(packed.atoms) != expected_atoms:
        raise ValueError(
            f"packed pdb atom count mismatch: expected {expected_atoms}, got {len(packed.atoms)}"
        )

    atom_rows, atom_type_rows = _atom_parameters(molecule)
    bond_rows, bond_type_rows = _bond_parameters(molecule)
    angle_rows, angle_type_rows = _angle_parameters(molecule)
    dihedral_rows, dihedral_type_rows = _dihedral_parameters(molecule)

    output_path = Path(filename)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "LAMMPS data file: RDKit UFF screening force field",
        "",
        f"{expected_atoms} atoms",
        f"{len(bond_rows) * molecule_count} bonds",
        f"{len(angle_rows) * molecule_count} angles",
        f"{len(dihedral_rows) * molecule_count} dihedrals",
        "0 impropers",
        "",
        f"{len(atom_type_rows)} atom types",
        f"{len(bond_type_rows)} bond types",
        f"{len(angle_type_rows)} angle types",
        f"{len(dihedral_type_rows)} dihedral types",
        "0 improper types",
        "",
        f"0.000000 {float(box_size):.6f} xlo xhi",
        f"0.000000 {float(box_size):.6f} ylo yhi",
        f"0.000000 {float(box_size):.6f} zlo zhi",
        "",
        "Masses",
        "",
    ]

    for type_id, row in enumerate(atom_type_rows, start=1):
        lines.append(f"{type_id} {row['mass']:.8f} # {row['label']}")
    lines.extend(["", "Pair Coeffs", ""])
    for type_id, row in enumerate(atom_type_rows, start=1):
        lines.append(f"{type_id} {row['epsilon']:.10f} {row['sigma']:.10f}")

    _append_coefficients(lines, "Bond Coeffs", bond_type_rows)
    _append_coefficients(lines, "Angle Coeffs", angle_type_rows)
    _append_coefficients(lines, "Dihedral Coeffs", dihedral_type_rows)

    lines.extend(["", "Atoms # full", ""])
    for molecule_index in range(molecule_count):
        atom_offset = molecule_index * atoms_per_molecule
        for local_index, atom_row in enumerate(atom_rows):
            atom_id = atom_offset + local_index + 1
            coordinate = packed.atoms[atom_id - 1]
            lines.append(
                f"{atom_id} {molecule_index + 1} {atom_row['type_id']} "
                f"{atom_row['charge']:.10f} {coordinate.x:.8f} "
                f"{coordinate.y:.8f} {coordinate.z:.8f}"
            )

    _append_topology(lines, "Bonds", bond_rows, molecule_count, atoms_per_molecule)
    _append_topology(lines, "Angles", angle_rows, molecule_count, atoms_per_molecule)
    _append_topology(lines, "Dihedrals", dihedral_rows, molecule_count, atoms_per_molecule)

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(output_path)


def _atom_parameters(molecule):
    type_map = {}
    type_rows = []
    atom_rows = []
    for atom in molecule.GetAtoms():
        index = atom.GetIdx()
        vdw = rdForceFieldHelpers.GetUFFVdWParams(molecule, index, index)
        if vdw is None:
            raise ValueError(f"UFF van der Waals parameters missing for atom {index}")
        x_min, epsilon = vdw
        sigma = x_min / (2.0 ** (1.0 / 6.0))
        key = (atom.GetAtomicNum(), round(epsilon, 10), round(sigma, 10))
        if key not in type_map:
            type_map[key] = len(type_rows) + 1
            symbol = atom.GetSymbol()
            type_rows.append(
                {
                    "mass": ATOMIC_MASSES.get(symbol, atom.GetMass()),
                    "epsilon": epsilon,
                    "sigma": sigma,
                    "label": f"UFF {symbol}",
                }
            )
        atom_rows.append(
            {
                "type_id": type_map[key],
                "charge": _finite_charge(atom),
            }
        )
    return atom_rows, type_rows


def _bond_parameters(molecule):
    rows = []
    type_map = {}
    type_rows = []
    for bond in molecule.GetBonds():
        atom_i = bond.GetBeginAtomIdx()
        atom_j = bond.GetEndAtomIdx()
        values = rdForceFieldHelpers.GetUFFBondStretchParams(molecule, atom_i, atom_j)
        if values is None:
            raise ValueError(f"UFF bond parameters missing for atoms {atom_i}-{atom_j}")
        force_constant, equilibrium = values
        coefficients = (force_constant / 2.0, equilibrium)
        type_id = _coefficient_type(coefficients, type_map, type_rows)
        rows.append((type_id, atom_i + 1, atom_j + 1))
    return rows, type_rows


def _angle_parameters(molecule):
    rows = []
    type_map = {}
    type_rows = []
    for center in molecule.GetAtoms():
        center_index = center.GetIdx()
        neighbors = [atom.GetIdx() for atom in center.GetNeighbors()]
        for atom_i, atom_k in combinations(neighbors, 2):
            values = rdForceFieldHelpers.GetUFFAngleBendParams(
                molecule, atom_i, center_index, atom_k
            )
            if values is None:
                raise ValueError(
                    f"UFF angle parameters missing for atoms {atom_i}-{center_index}-{atom_k}"
                )
            force_constant, equilibrium_degrees = values
            coefficients = (force_constant / 2.0, equilibrium_degrees)
            type_id = _coefficient_type(coefficients, type_map, type_rows)
            rows.append((type_id, atom_i + 1, center_index + 1, atom_k + 1))
    return rows, type_rows


def _dihedral_parameters(molecule):
    rows = []
    type_map = {}
    type_rows = []
    seen = set()
    for central_bond in molecule.GetBonds():
        atom_j = central_bond.GetBeginAtomIdx()
        atom_k = central_bond.GetEndAtomIdx()
        left = [atom.GetIdx() for atom in molecule.GetAtomWithIdx(atom_j).GetNeighbors() if atom.GetIdx() != atom_k]
        right = [atom.GetIdx() for atom in molecule.GetAtomWithIdx(atom_k).GetNeighbors() if atom.GetIdx() != atom_j]
        for atom_i in left:
            for atom_l in right:
                if atom_i == atom_l:
                    continue
                canonical = min(
                    (atom_i, atom_j, atom_k, atom_l),
                    (atom_l, atom_k, atom_j, atom_i),
                )
                if canonical in seen:
                    continue
                seen.add(canonical)
                barrier = rdForceFieldHelpers.GetUFFTorsionParams(
                    molecule, atom_i, atom_j, atom_k, atom_l
                )
                if barrier is None:
                    continue
                sign, multiplicity = _torsion_shape(
                    molecule.GetAtomWithIdx(atom_j),
                    molecule.GetAtomWithIdx(atom_k),
                )
                coefficients = (barrier / 2.0, sign, multiplicity)
                type_id = _coefficient_type(coefficients, type_map, type_rows)
                rows.append(
                    (type_id, atom_i + 1, atom_j + 1, atom_k + 1, atom_l + 1)
                )
    return rows, type_rows


def _torsion_shape(atom_j, atom_k):
    sp2 = Chem.HybridizationType.SP2
    sp3 = Chem.HybridizationType.SP3
    if atom_j.GetHybridization() == sp2 and atom_k.GetHybridization() == sp2:
        return -1, 2
    if atom_j.GetHybridization() == sp3 and atom_k.GetHybridization() == sp3:
        return 1, 3
    return 1, 6


def _coefficient_type(coefficients, type_map, type_rows):
    key = tuple(round(float(value), 10) if isinstance(value, float) else value for value in coefficients)
    if key not in type_map:
        type_map[key] = len(type_rows) + 1
        type_rows.append(coefficients)
    return type_map[key]


def _append_coefficients(lines, title, rows):
    if not rows:
        return
    lines.extend(["", title, ""])
    for type_id, coefficients in enumerate(rows, start=1):
        values = " ".join(_format_coefficient(value) for value in coefficients)
        lines.append(f"{type_id} {values}")


def _append_topology(lines, title, rows, molecule_count, atoms_per_molecule):
    if not rows:
        return
    lines.extend(["", title, ""])
    item_id = 1
    for molecule_index in range(int(molecule_count)):
        atom_offset = molecule_index * atoms_per_molecule
        for row in rows:
            type_id, *atom_ids = row
            shifted = " ".join(str(atom_id + atom_offset) for atom_id in atom_ids)
            lines.append(f"{item_id} {type_id} {shifted}")
            item_id += 1


def _format_coefficient(value):
    if isinstance(value, int):
        return str(value)
    return f"{float(value):.10f}"


def _finite_charge(atom):
    try:
        charge = float(atom.GetProp("_GasteigerCharge"))
    except (KeyError, ValueError):
        return 0.0
    return charge if math.isfinite(charge) else 0.0
