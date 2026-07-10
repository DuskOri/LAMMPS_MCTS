"""把 PDB 坐标和键连接关系整理成 LAMMPS data 文件。"""

from dataclasses import dataclass
from pathlib import Path


ATOMIC_MASSES = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "S": 32.06,
    "P": 30.974,
    "Cl": 35.45,
    "Br": 79.904,
    "I": 126.904,
}


@dataclass
class PDBAtom:
    """保存 PDB 中一个原子的基本信息。"""

    serial: int
    name: str
    element: str
    x: float
    y: float
    z: float


@dataclass
class PDBTopology:
    """PDB 文件解析后的坐标和键表。"""

    atoms: list
    bonds: list


def read_pdb_topology(filename):
    """读取 PDB 文件中的 HETATM/ATOM 坐标和 CONECT 键连接。"""
    atoms = []
    serial_to_index = {}
    bonds = set()

    for line in Path(filename).read_text(encoding="utf-8", errors="ignore").splitlines():
        record = line[:6].strip()
        if record in ("ATOM", "HETATM"):
            atom = _parse_atom_line(line)
            serial_to_index[atom.serial] = len(atoms) + 1
            atoms.append(atom)
        elif record == "CONECT":
            _read_conect_line(line, serial_to_index, bonds)

    return PDBTopology(atoms=atoms, bonds=sorted(bonds))


def write_lammps_data(topology, filename, box=None, comment="generated from pdb"):
    """根据拓扑对象写出 LAMMPS data 文件，使用 atom_style full。"""
    output_file = Path(filename)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    atom_types = _build_atom_types(topology.atoms)
    bond_types = _build_bond_types(topology.atoms, topology.bonds)
    bounds = box or _auto_box(topology.atoms)

    lines = []
    lines.append(f"LAMMPS data file: {comment}")
    lines.append("")
    lines.append(f"{len(topology.atoms)} atoms")
    lines.append(f"{len(topology.bonds)} bonds")
    lines.append("")
    lines.append(f"{len(atom_types)} atom types")
    lines.append(f"{len(bond_types)} bond types")
    lines.append("")
    lines.append(f"{bounds[0][0]:.6f} {bounds[0][1]:.6f} xlo xhi")
    lines.append(f"{bounds[1][0]:.6f} {bounds[1][1]:.6f} ylo yhi")
    lines.append(f"{bounds[2][0]:.6f} {bounds[2][1]:.6f} zlo zhi")
    lines.append("")
    lines.append("Masses")
    lines.append("")

    for element, type_id in atom_types.items():
        mass = ATOMIC_MASSES.get(element, 12.0)
        lines.append(f"{type_id} {mass:.6f} # {element}")

    lines.append("")
    lines.append("Atoms # full")
    lines.append("")

    for atom_id, atom in enumerate(topology.atoms, start=1):
        mol_id = getattr(atom, "molecule_id", 1)
        atom_type = atom_types[atom.element]
        charge = 0.0
        lines.append(
            f"{atom_id} {mol_id} {atom_type} {charge:.6f} "
            f"{atom.x:.6f} {atom.y:.6f} {atom.z:.6f} # {atom.element}"
        )

    lines.append("")
    lines.append("Bonds")
    lines.append("")

    for bond_id, (atom_i, atom_j) in enumerate(topology.bonds, start=1):
        bond_type = bond_types[_bond_key(topology.atoms, atom_i, atom_j)]
        lines.append(f"{bond_id} {bond_type} {atom_i} {atom_j}")

    output_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(output_file)


def write_lammps_data_from_template(template_pdb, coordinate_pdb, filename, molecule_count=1, box_size=60.0):
    """用单链模板提供键表，用 Packmol 输出的 PDB 提供多链坐标。"""
    template = read_pdb_topology(template_pdb)
    coordinates = read_pdb_topology(coordinate_pdb)

    atoms_per_molecule = len(template.atoms)
    if atoms_per_molecule == 0:
        raise ValueError("template pdb contains no atoms")

    expected_atoms = atoms_per_molecule * molecule_count
    if len(coordinates.atoms) != expected_atoms:
        raise ValueError(
            f"packed pdb atom count mismatch: expected {expected_atoms}, got {len(coordinates.atoms)}"
        )

    bonds = []
    for molecule_id in range(molecule_count):
        atom_offset = molecule_id * atoms_per_molecule
        for atom_i, atom_j in template.bonds:
            bonds.append((atom_i + atom_offset, atom_j + atom_offset))

    for index, atom in enumerate(coordinates.atoms):
        atom.molecule_id = index // atoms_per_molecule + 1

    topology = PDBTopology(atoms=coordinates.atoms, bonds=bonds)
    box = ((0.0, box_size), (0.0, box_size), (0.0, box_size))
    return write_lammps_data(topology, filename, box=box, comment="packed polymer system")


def _parse_atom_line(line):
    """解析 PDB 原子行，元素优先使用 77-78 列。"""
    serial = int(line[6:11])
    name = line[12:16].strip()
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    element = line[76:78].strip()
    if not element:
        element = _guess_element(name)
    return PDBAtom(serial=serial, name=name, element=element, x=x, y=y, z=z)


def _read_conect_line(line, serial_to_index, bonds):
    """读取 CONECT 行，并把 PDB serial 转成连续原子编号。"""
    numbers = [int(part) for part in line.split()[1:] if part.isdigit()]
    if len(numbers) < 2:
        return

    source = serial_to_index.get(numbers[0])
    if source is None:
        return

    for target_serial in numbers[1:]:
        target = serial_to_index.get(target_serial)
        if target is None or target == source:
            continue
        bonds.add(tuple(sorted((source, target))))


def _guess_element(atom_name):
    """从原子名推断元素，适配 RDKit 输出的 C1、O2、N3 等命名。"""
    letters = "".join(ch for ch in atom_name if ch.isalpha())
    if not letters:
        return "C"

    if len(letters) >= 2 and letters[:2].capitalize() in ATOMIC_MASSES:
        return letters[:2].capitalize()
    return letters[0].upper()


def _build_atom_types(atoms):
    """按照元素名称建立 LAMMPS atom type 编号。"""
    elements = sorted({atom.element for atom in atoms})
    return {element: index for index, element in enumerate(elements, start=1)}


def _build_bond_types(atoms, bonds):
    """按照两端元素组合建立 bond type 编号。"""
    keys = sorted({_bond_key(atoms, atom_i, atom_j) for atom_i, atom_j in bonds})
    return {key: index for index, key in enumerate(keys, start=1)}


def _bond_key(atoms, atom_i, atom_j):
    """把一根键转成元素对，例如 C-O、C-H。"""
    element_i = atoms[atom_i - 1].element
    element_j = atoms[atom_j - 1].element
    return tuple(sorted((element_i, element_j)))


def _auto_box(atoms, padding=5.0):
    """根据坐标范围自动给单链 data 文件设置盒子边界。"""
    xs = [atom.x for atom in atoms]
    ys = [atom.y for atom in atoms]
    zs = [atom.z for atom in atoms]
    return (
        (min(xs) - padding, max(xs) + padding),
        (min(ys) - padding, max(ys) + padding),
        (min(zs) - padding, max(zs) + padding),
    )
