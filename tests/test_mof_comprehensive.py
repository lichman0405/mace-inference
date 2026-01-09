"""
真实 MOF 材料全面测试 - MACE Inference

测试材料：HKUST-1 (Cu-BTC) 的简化模型
测试任务：
1. 单点能量计算
2. 结构优化
3. NVT 分子动力学
4. NPT 分子动力学
5. 声子计算与热力学性质
6. 体积模量计算
7. CO2 吸附能计算
8. 配位环境分析
"""

import numpy as np
from pathlib import Path
from ase import Atoms
from ase.build import molecule
from ase.io import write, read
import tempfile
import os

# 创建输出目录
OUTPUT_DIR = Path("test_outputs")
OUTPUT_DIR.mkdir(exist_ok=True)


def create_simple_mof():
    """
    创建一个简化的 Cu-paddlewheel MOF 结构
    基于 HKUST-1 的 Cu2 金属节点
    """
    # Cu-paddlewheel 单元 + 有机连接
    # Cu2(O2C-)4 结构
    positions = [
        # Cu paddlewheel 中心
        [5.0, 5.0, 5.0],      # Cu1
        [5.0, 5.0, 7.5],      # Cu2
        # 羧酸氧原子 (连接 Cu)
        [6.5, 5.0, 5.8],      # O1
        [6.5, 5.0, 6.7],      # O2
        [3.5, 5.0, 5.8],      # O3
        [3.5, 5.0, 6.7],      # O4
        [5.0, 6.5, 5.8],      # O5
        [5.0, 6.5, 6.7],      # O6
        [5.0, 3.5, 5.8],      # O7
        [5.0, 3.5, 6.7],      # O8
        # 羧酸碳原子
        [7.2, 5.0, 6.25],     # C1
        [2.8, 5.0, 6.25],     # C2
        [5.0, 7.2, 6.25],     # C3
        [5.0, 2.8, 6.25],     # C4
        # 苯环碳原子 (简化)
        [8.5, 5.0, 6.25],     # C5
        [1.5, 5.0, 6.25],     # C6
        [5.0, 8.5, 6.25],     # C7
        [5.0, 1.5, 6.25],     # C8
        # 氢原子
        [9.2, 5.0, 6.25],     # H1
        [0.8, 5.0, 6.25],     # H2
        [5.0, 9.2, 6.25],     # H3
        [5.0, 0.8, 6.25],     # H4
    ]
    
    symbols = ['Cu', 'Cu', 
               'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O',
               'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
               'H', 'H', 'H', 'H']
    
    mof = Atoms(
        symbols=symbols,
        positions=positions,
        cell=[10.0, 10.0, 12.5],
        pbc=True
    )
    
    return mof


def main():
    from mace_inference import MACEInference
    
    print("=" * 70)
    print("🧪 真实 MOF 材料全面测试 - MACE Inference")
    print("=" * 70)
    
    # 初始化
    print("\n📦 初始化 MACEInference...")
    calc = MACEInference(model='medium', device='cpu')
    print(f"   设备: {calc.device}")
    print(f"   模型: {calc.model_name}")
    
    # 创建 MOF 结构
    print("\n📐 创建 Cu-paddlewheel MOF 结构...")
    mof = create_simple_mof()
    print(f"   原子数: {len(mof)}")
    print(f"   元素: {set(mof.get_chemical_symbols())}")
    print(f"   晶胞: {mof.cell.lengths()}")
    
    # 保存初始结构
    write(OUTPUT_DIR / "mof_initial.xyz", mof)
    print(f"   已保存: {OUTPUT_DIR / 'mof_initial.xyz'}")
    
    results = {}
    
    # =========================================================================
    # 1. 单点能量计算
    # =========================================================================
    print("\n" + "=" * 70)
    print("1️⃣  单点能量计算")
    print("=" * 70)
    
    sp_result = calc.single_point(mof)
    results['single_point'] = sp_result
    
    print(f"   总能量:       {sp_result['energy']:.4f} eV")
    print(f"   每原子能量:   {sp_result['energy_per_atom']:.4f} eV/atom")
    print(f"   最大力:       {sp_result['max_force']:.4f} eV/Å")
    print(f"   RMS 力:       {sp_result['rms_force']:.4f} eV/Å")
    if sp_result.get('pressure_GPa') is not None:
        print(f"   压力:         {sp_result['pressure_GPa']:.4f} GPa")
    print("   ✅ 单点能量计算完成!")
    
    # =========================================================================
    # 2. 结构优化
    # =========================================================================
    print("\n" + "=" * 70)
    print("2️⃣  结构优化")
    print("=" * 70)
    
    print("   正在优化... (fmax=0.1, steps=50)")
    mof_opt = calc.optimize(
        mof, 
        fmax=0.1, 
        steps=50,
        optimizer="LBFGS"
    )
    
    opt_result = calc.single_point(mof_opt)
    results['optimization'] = opt_result
    
    print(f"   优化后能量:   {opt_result['energy']:.4f} eV")
    print(f"   能量变化:     {opt_result['energy'] - sp_result['energy']:.4f} eV")
    print(f"   优化后最大力: {opt_result['max_force']:.4f} eV/Å")
    
    write(OUTPUT_DIR / "mof_optimized.xyz", mof_opt)
    print(f"   已保存: {OUTPUT_DIR / 'mof_optimized.xyz'}")
    print("   ✅ 结构优化完成!")
    
    # =========================================================================
    # 3. NVT 分子动力学
    # =========================================================================
    print("\n" + "=" * 70)
    print("3️⃣  NVT 分子动力学 (300K, 50步)")
    print("=" * 70)
    
    traj_file = str(OUTPUT_DIR / "mof_nvt.traj")
    
    print("   正在运行 NVT MD...")
    final_nvt = calc.run_md(
        mof_opt,
        ensemble="nvt",
        temperature_K=300,
        steps=50,
        timestep=1.0,
        trajectory=traj_file,
        log_interval=10
    )
    
    nvt_result = calc.single_point(final_nvt)
    results['nvt_md'] = nvt_result
    
    print(f"   最终能量:     {nvt_result['energy']:.4f} eV")
    print(f"   轨迹已保存:   {traj_file}")
    print("   ✅ NVT MD 完成!")
    
    # =========================================================================
    # 4. NPT 分子动力学
    # =========================================================================
    print("\n" + "=" * 70)
    print("4️⃣  NPT 分子动力学 (300K, 0.0 GPa, 50步)")
    print("=" * 70)
    
    traj_file_npt = str(OUTPUT_DIR / "mof_npt.traj")
    
    print("   正在运行 NPT MD...")
    final_npt = calc.run_md(
        mof_opt,
        ensemble="npt",
        temperature_K=300,
        pressure_GPa=0.0,
        steps=50,
        timestep=1.0,
        trajectory=traj_file_npt,
        log_interval=10
    )
    
    npt_result = calc.single_point(final_npt)
    results['npt_md'] = npt_result
    
    print(f"   最终能量:     {npt_result['energy']:.4f} eV")
    print(f"   最终体积:     {final_npt.get_volume():.2f} ų")
    print(f"   轨迹已保存:   {traj_file_npt}")
    print("   ✅ NPT MD 完成!")
    
    # =========================================================================
    # 5. 声子计算
    # =========================================================================
    print("\n" + "=" * 70)
    print("5️⃣  声子计算与热力学性质")
    print("=" * 70)
    
    print("   正在计算声子... (supercell=[1,1,1], 因结构较大)")
    phonon_result = calc.phonon(
        mof_opt,
        supercell_matrix=1,  # 使用 1x1x1 因为结构已经够大
        displacement=0.01,
        mesh=[5, 5, 5],
        temperature_range=(0, 500, 50)
    )
    results['phonon'] = phonon_result
    
    print(f"   位移数量:     {len(phonon_result['phonon'].supercells_with_displacements)}")
    
    if 'thermal' in phonon_result:
        thermal = phonon_result['thermal']
        # 找到 300K 对应的索引
        temps = thermal['temperatures']
        idx_300 = np.argmin(np.abs(temps - 300))
        
        print(f"\n   300K 热力学性质:")
        print(f"   自由能:       {thermal['free_energy'][idx_300]:.3f} kJ/mol")
        print(f"   熵:           {thermal['entropy'][idx_300]:.3f} J/(mol·K)")
        print(f"   热容:         {thermal['heat_capacity'][idx_300]:.3f} J/(mol·K)")
    
    print("   ✅ 声子计算完成!")
    
    # =========================================================================
    # 6. 体积模量计算
    # =========================================================================
    print("\n" + "=" * 70)
    print("6️⃣  体积模量计算")
    print("=" * 70)
    
    print("   正在计算体积模量... (5个体积点)")
    bulk_result = calc.bulk_modulus(
        mof_opt,
        n_points=5,
        scale_range=(0.98, 1.02),
        optimize_first=False
    )
    results['bulk_modulus'] = bulk_result
    
    print(f"   平衡体积:     {bulk_result['v0']:.2f} ų")
    print(f"   平衡能量:     {bulk_result['e0']:.4f} eV")
    print(f"   体积模量:     {bulk_result['B_GPa']:.2f} GPa")
    print("   ✅ 体积模量计算完成!")
    
    # =========================================================================
    # 7. CO2 吸附能计算
    # =========================================================================
    print("\n" + "=" * 70)
    print("7️⃣  CO2 吸附能计算")
    print("=" * 70)
    
    # 在 Cu 附近放置 CO2
    adsorption_site = [5.0, 5.0, 3.0]  # Cu paddlewheel 附近
    
    print(f"   吸附位点:     {adsorption_site}")
    print("   正在计算吸附能...")
    
    ads_result = calc.adsorption_energy(
        framework=mof_opt,
        adsorbate="CO2",
        site_position=adsorption_site,
        optimize=True,
        fmax=0.1,
        fix_framework=True
    )
    results['adsorption'] = ads_result
    
    print(f"   MOF 能量:     {ads_result['E_mof']:.4f} eV")
    print(f"   CO2 能量:     {ads_result['E_gas']:.4f} eV")
    print(f"   复合物能量:   {ads_result['E_complex']:.4f} eV")
    print(f"   吸附能:       {ads_result['E_ads']:.4f} eV")
    print(f"   吸附能:       {ads_result['E_ads'] * 96.485:.2f} kJ/mol")
    
    write(OUTPUT_DIR / "mof_co2_complex.xyz", ads_result['complex_structure'])
    print(f"   已保存: {OUTPUT_DIR / 'mof_co2_complex.xyz'}")
    print("   ✅ 吸附能计算完成!")
    
    # =========================================================================
    # 8. 配位环境分析
    # =========================================================================
    print("\n" + "=" * 70)
    print("8️⃣  配位环境分析")
    print("=" * 70)
    
    coord_result = calc.coordination(mof_opt)
    results['coordination'] = coord_result
    
    print(f"   金属中心数量: {coord_result['n_metal_centers']}")
    print(f"   金属原子索引: {coord_result['metal_indices']}")
    
    for metal_idx, info in coord_result['coordination'].items():
        print(f"\n   {info['metal_symbol']} (原子 {metal_idx}):")
        print(f"     配位数:     {info['coordination_number']}")
        print(f"     平均距离:   {info['average_distance']:.3f} Å")
        print(f"     配位原子:   ", end="")
        neighbor_symbols = [n['symbol'] for n in info['neighbors']]
        from collections import Counter
        counts = Counter(neighbor_symbols)
        print(", ".join([f"{sym}×{cnt}" for sym, cnt in counts.items()]))
    
    print("\n   ✅ 配位分析完成!")
    
    # =========================================================================
    # 总结
    # =========================================================================
    print("\n" + "=" * 70)
    print("📊 测试总结")
    print("=" * 70)
    
    print("""
    ┌─────────────────────────────────────────────────────────────────┐
    │  任务                    │  状态   │  关键结果                  │
    ├─────────────────────────────────────────────────────────────────┤""")
    print(f"    │  1. 单点能量计算         │   ✅    │  E = {results['single_point']['energy']:.2f} eV")
    print(f"    │  2. 结构优化             │   ✅    │  ΔE = {results['optimization']['energy'] - results['single_point']['energy']:.2f} eV")
    print(f"    │  3. NVT 分子动力学       │   ✅    │  T = 300 K, 50 步")
    print(f"    │  4. NPT 分子动力学       │   ✅    │  T = 300 K, P = 0 GPa")
    print(f"    │  5. 声子计算             │   ✅    │  热力学性质已计算")
    print(f"    │  6. 体积模量             │   ✅    │  B = {results['bulk_modulus']['B_GPa']:.1f} GPa")
    print(f"    │  7. CO2 吸附能           │   ✅    │  E_ads = {results['adsorption']['E_ads']:.3f} eV")
    print(f"    │  8. 配位分析             │   ✅    │  {results['coordination']['n_metal_centers']} 个金属中心")
    print("""    └─────────────────────────────────────────────────────────────────┘
    """)
    
    print("🎉 所有 8 项推理任务全部完成!")
    print(f"📁 输出文件保存在: {OUTPUT_DIR.absolute()}")
    print("=" * 70)
    
    return results


if __name__ == "__main__":
    results = main()
