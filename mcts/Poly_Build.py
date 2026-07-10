# generator/Poly_Build.py
# 由MCTS算法输出的聚合物片段序列，构建成真实的聚合物分子，并保存为PDB文件

import os
from rdkit import Chem
from rdkit.Chem import AllChem

try:
    from .fragment_registry import get_fragment_registry
except ImportError:
    from fragment_registry import get_fragment_registry

# 将 MCTS 标签映射为真实的 SMILES 片段。
# `[1*]` 代表 片段 的左连接位点，`[2*]` 代表 片段 的右连接位点
SMILES_MAP = get_fragment_registry().smiles_map()

# 拼接主函数
def build_poly_chain(fragment_list, dp=10):
    """
    接受一个 MCTS 序列（如 ['A','B','C','End']），拼接成一条长度为 dp 的聚合物链
    """
    # 过滤掉 End
    smiles_map = get_fragment_registry().smiles_map()
    real_fragments = [f for f in fragment_list if f != 'End']
    if not real_fragments:
        return None
    missing = [fragment for fragment in real_fragments if fragment not in smiles_map]
    if missing:
        print(f"ERROR：片段库中没有这些片段: {missing}")
        return None

    # 2. 将传入的链段乘以聚合度 DP，即可得到完整的片段列表
    full_chain = real_fragments * dp
    
    # 3. 定义化学拼接方法deicide
    # 在第一段聚合物链找到它的右位点 [2*]，并把与它相连的那个真实的末端原子标记为 1号 (*:1)。
    # 下一个准备拼进来的单体，找到它的左位点 [1*]，并把与它相连的那个真实的头部原子标记为 2号 (*:2)。
    # >>[*:1]-[*:2]：两个虚拟接头（[1*]和[2*]）删除后，然后在1号原子 和 2号原子 之间打上一个单键（-）
    deicide = AllChem.ReactionFromSmarts('[*:1][2*].[1*][*:2]>>[*:1]-[*:2]')

    # 4. 初始化聚合物链，在RDKIT中解析为第一个片段
    chain_mol = Chem.MolFromSmiles(smiles_map[full_chain[0]])
    if chain_mol is None:
        print(f"ERROR：无法解析起始片段 {full_chain[0]}")
        return None
    
    print(f"{full_chain[0]} 片段解析结果:", Chem.MolToSmiles(chain_mol))
    if 'B' in smiles_map:
        print("B 片段解析结果:", Chem.MolToSmiles(Chem.MolFromSmiles(smiles_map['B'])))

    # 5. 循环拼接剩下的片段
    for frag_name in full_chain[1:]:
        next_mol = Chem.MolFromSmiles(smiles_map[frag_name])
        
        # 运行反应，这会返回一个产物列表
        products = deicide.RunReactants((chain_mol, next_mol))
        
        if products and len(products) > 0:
            # 取第一个反应物组合的第一个产物（RDKIT特性）
            chain_mol = products[0][0]
            # 反应后必须清理分子结构（计算化合价、芳香性等）
            Chem.SanitizeMol(chain_mol)
        else:
            print(f"ERROR：片段 {frag_name} 拼接失败！")
            return None

    return chain_mol

# 封端函数
def cap_poly_ends(chain_mol):
    # 1. 遍历分子中的所有原子
    for atom in chain_mol.GetAtoms():
        # 虚拟原子 [*] 的原子序数是 0
        if atom.GetAtomicNum() == 0:
            atom.SetAtomicNum(1)  # 把原子序数改成 1，也就是氢原子(H)
            atom.SetIsotope(0)    # 之前 [1*] 和 [2*] 里的标记 1 和 2 删除
            
    # 2. 清理和刷新分子属性（重新计算化合价等）
    Chem.SanitizeMol(chain_mol)
    
    # 3. 为 LAMMPS 准备“显式氢”
    # RDKit 默认为了节省空间，普通氢原子默认隐藏，
    # 但做 MD 模拟时，每一个氢原子都必须是节点，必须为显示
    final_mol = Chem.AddHs(chain_mol)
    
    return final_mol


#将分子对象保存为 PDB 文件（包含三维坐标）
def save_pdb(final_mol, filename="output.pdb"):
    if final_mol is None:
        return False

    # 生成三维构象
    params = AllChem.ETKDGv3()
    params.randomSeed = 2026
    conf_id = AllChem.EmbedMolecule(final_mol, params)
    if conf_id < 0:
        params.useRandomCoords = True
        conf_id = AllChem.EmbedMolecule(final_mol, params)
    if conf_id < 0:
        print(f"ERROR 无法生成三维构象: {filename}")
        return False

    # 简单的力场优化，防止原子过度重叠
    try:
        AllChem.UFFOptimizeMolecule(final_mol, confId=conf_id, maxIters=200)
    except ValueError:
        print(f"WARNING UFF优化失败，直接保存未优化构象: {filename}")
    
    # 保存为 PDB 文件
    Chem.MolToPDBFile(final_mol, filename)
    print(f"CORRECT 聚合物链已成功生成并保存为: {filename}")
    return True



# ================= 测试运行 =================
if __name__ == "__main__":
    # 假设 MCTS 生成的重复单元序列
    mcts_sequence = ['A', 'B', "A'",'B', 'C', 'End']
    
    # 生成 DP=2 的聚合物
    polymer = build_poly_chain(mcts_sequence, dp=2)
    cap_poly = cap_poly_ends(polymer)
    if cap_poly:
        # 输出最终生成的聚合物的 SMILES
        final_smiles = Chem.MolToSmiles(cap_poly)
        print("CORRECT 聚合成功！最终长链的 SMILES 为：\n")
        print(final_smiles)
