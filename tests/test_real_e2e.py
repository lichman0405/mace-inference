"""
真实端到端测试 - 使用真实的 MACE 模型和结构

这个测试验证库的真正功能，不使用 mock。
需要安装 mace-torch 并且有网络连接（首次运行需要下载模型）。
"""

import pytest
import numpy as np
from ase.build import bulk, molecule


class TestRealE2E:
    """端到端真实测试"""

    @pytest.fixture
    def calc(self):
        """创建真实的 MACEInference 实例"""
        from mace_inference import MACEInference
        return MACEInference(model='medium', device='cpu')

    def test_single_point_cu_bulk(self, calc):
        """测试单点能量计算 - Cu 晶体"""
        cu = bulk('Cu', 'fcc', a=3.6)
        result = calc.single_point(cu)
        
        assert 'energy' in result
        assert 'energy_per_atom' in result
        assert 'forces' in result
        assert 'max_force' in result
        
        # Cu 的能量应该是负数
        assert result['energy'] < 0
        print(f"\nCu 能量: {result['energy']:.4f} eV")
        print(f"每原子能量: {result['energy_per_atom']:.4f} eV/atom")

    def test_optimize_structure(self, calc):
        """测试结构优化 - 压缩的 Cu 晶体"""
        cu = bulk('Cu', 'fcc', a=3.5)  # 稍微压缩
        
        # 优化前
        before = calc.single_point(cu)
        
        # 优化
        optimized = calc.optimize(cu, fmax=0.1, steps=20)
        
        # 优化后
        after = calc.single_point(optimized)
        
        # 优化后能量应该更低
        assert after['energy'] <= before['energy']
        print(f"\n优化前能量: {before['energy']:.4f} eV")
        print(f"优化后能量: {after['energy']:.4f} eV")

    def test_bulk_modulus(self, calc):
        """测试体积模量计算 - Cu 晶体"""
        cu = bulk('Cu', 'fcc', a=3.6)
        
        result = calc.bulk_modulus(cu, n_points=5, optimize_first=False)
        
        assert 'v0' in result
        assert 'B_GPa' in result
        
        # Cu 的体积模量约 140 GPa，允许较大误差
        assert 50 < result['B_GPa'] < 300
        print(f"\n体积模量: {result['B_GPa']:.1f} GPa")

    def test_coordination_analysis(self, calc):
        """测试配位分析"""
        from mace_inference.tasks.adsorption import analyze_coordination
        from ase import Atoms
        
        # 创建简单的 Cu-O 结构
        atoms = Atoms(
            'Cu2O4',
            positions=[
                [0, 0, 0],      # Cu
                [5, 5, 5],      # Cu
                [1.5, 0, 0],    # O
                [0, 1.5, 0],    # O
                [6.5, 5, 5],    # O
                [5, 6.5, 5],    # O
            ],
            cell=[10, 10, 10],
            pbc=True
        )
        
        result = analyze_coordination(atoms)
        
        assert 'coordination' in result
        assert 'n_metal_centers' in result
        assert result['n_metal_centers'] == 2
        print(f"\n金属中心数量: {result['n_metal_centers']}")


if __name__ == '__main__':
    # 可以直接运行此脚本
    from mace_inference import MACEInference
    
    print('='*60)
    print('🧪 真实端到端测试 - MACE Inference')
    print('='*60)
    
    # 1. 初始化计算器
    print('\n1️⃣ 初始化 MACEInference...')
    calc = MACEInference(model='medium', device='cpu')
    print(f'   ✅ 成功! 使用设备: {calc.device}')
    
    # 2. 单点能量计算
    print('\n2️⃣ 单点能量计算 (Cu bulk)...')
    cu = bulk('Cu', 'fcc', a=3.6)
    result = calc.single_point(cu)
    print(f'   能量: {result["energy"]:.4f} eV')
    print(f'   每原子能量: {result["energy_per_atom"]:.4f} eV/atom')
    print(f'   最大力: {result["max_force"]:.6f} eV/Å')
    print('   ✅ 单点计算成功!')
    
    # 3. 结构优化
    print('\n3️⃣ 结构优化 (Cu bulk)...')
    cu_opt = bulk('Cu', 'fcc', a=3.5)
    optimized = calc.optimize(cu_opt, fmax=0.1, steps=10)
    opt_result = calc.single_point(optimized)
    print(f'   优化后能量: {opt_result["energy"]:.4f} eV')
    print(f'   优化后最大力: {opt_result["max_force"]:.6f} eV/Å')
    print('   ✅ 结构优化成功!')
    
    # 4. 体积模量计算
    print('\n4️⃣ 体积模量计算 (Cu bulk)...')
    cu_for_bulk = bulk('Cu', 'fcc', a=3.6)
    bulk_result = calc.bulk_modulus(cu_for_bulk, n_points=5, optimize_first=False)
    print(f'   平衡体积: {bulk_result["v0"]:.3f} Å³')
    print(f'   体积模量: {bulk_result["B_GPa"]:.1f} GPa')
    print('   ✅ 体积模量计算成功!')
    
    print('\n' + '='*60)
    print('🎉 所有真实测试通过! 库可以正常工作!')
    print('='*60)
