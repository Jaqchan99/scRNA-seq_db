#!/usr/bin/env python3
"""
测试 Demo 脚本
用于验证 scRNA-seq 数据处理平台的完整功能
"""

import requests
import time
import json
import os
import sys
from pathlib import Path

# 测试配置
BASE_URL = "http://localhost:8100"
TEST_FILE = "test_data.h5ad"

def create_test_data():
    """创建测试用的 h5ad 文件"""
    try:
        import scanpy as sc
        import pandas as pd
        import numpy as np
        
        print("🧪 创建测试数据...")
        
        # 创建模拟的单细胞数据
        n_cells = 100
        n_genes = 500
        
        # 创建表达矩阵
        X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(float)
        
        # 创建细胞元数据
        obs = pd.DataFrame({
            'cell_id': [f'cell_{i}' for i in range(n_cells)],
            'raw_cell_type_label': np.random.choice([
                'T cell', 'B cell', 'NK cell', 'Monocyte', 'Dendritic cell',
                'Neutrophil', 'Eosinophil', 'Basophil', 'Macrophage', 'DC'
            ], n_cells),
            'n_genes': np.random.randint(200, 1000, n_cells),
            'n_counts': np.random.randint(1000, 5000, n_cells),
            'pct_mt': np.random.uniform(5, 25, n_cells),
            'batch': np.random.choice(['batch1', 'batch2', 'batch3'], n_cells)
        })
        
        # 创建基因元数据
        var = pd.DataFrame({
            'gene_symbol': [f'Gene_{i}' for i in range(n_genes)],
            'ensembl_gene_id': [f'ENSMUSG{i:011d}' for i in range(n_genes)]
        })
        
        # 创建 AnnData 对象
        adata = sc.AnnData(X=X, obs=obs, var=var)
        
        # 添加全局元数据
        adata.uns['organism'] = 'mouse'
        adata.uns['genome_build'] = 'GRCm39'
        adata.uns['platform'] = '10X'
        adata.uns['chemistry_version'] = 'v3'
        
        # 保存测试文件
        adata.write_h5ad(TEST_FILE)
        print(f"✅ 测试数据已创建: {TEST_FILE}")
        
        # 确保文件句柄关闭
        del adata
        
        return True
        
    except ImportError as e:
        print(f"❌ 缺少依赖包: {e}")
        print("请安装: pip install scanpy pandas numpy")
        return False
    except Exception as e:
        print(f"❌ 创建测试数据失败: {e}")
        return False

def test_api_endpoints():
    """测试 API 端点"""
    print("\n🚀 开始测试 API 端点...")
    
    try:
        # 1. 测试根端点
        print("1. 测试根端点...")
        response = requests.get(f"{BASE_URL}/")
        if response.status_code == 200:
            print("✅ 根端点正常")
            print(f"   响应: {response.json()}")
        else:
            print(f"❌ 根端点失败: {response.status_code}")
            return False
        
        # 2. 创建提交任务
        print("\n2. 创建提交任务...")
        response = requests.post(f"{BASE_URL}/submissions")
        if response.status_code == 200:
            submission_data = response.json()
            submission_id = submission_data["submission_id"]
            print(f"✅ 提交任务已创建: {submission_id}")
        else:
            print(f"❌ 创建提交任务失败: {response.status_code}")
            return False
        
        # 3. 上传文件
        print("\n3. 上传测试文件...")
        if not os.path.exists(TEST_FILE):
            print(f"❌ 测试文件不存在: {TEST_FILE}")
            return False
        
        with open(TEST_FILE, 'rb') as f:
            files = {'file': (TEST_FILE, f, 'application/octet-stream')}
            response = requests.post(
                f"{BASE_URL}/submissions/{submission_id}/upload",
                files=files
            )
        
        if response.status_code == 200:
            print("✅ 文件上传成功")
            print(f"   响应: {response.json()}")
        else:
            print(f"❌ 文件上传失败: {response.status_code}")
            print(f"   错误: {response.text}")
            return False
        
        # 4. 等待处理完成
        print("\n4. 等待处理完成...")
        max_wait_time = 300  # 5分钟
        wait_time = 0
        
        while wait_time < max_wait_time:
            response = requests.get(f"{BASE_URL}/submissions/{submission_id}/status")
            if response.status_code == 200:
                status_data = response.json()
                status = status_data["status"]
                print(f"   当前状态: {status}")
                
                if status == "completed":
                    print("✅ 处理完成!")
                    break
                elif status in ["failed", "validation_failed"]:
                    print(f"❌ 处理失败: {status}")
                    return False
                
                time.sleep(10)
                wait_time += 10
            else:
                print(f"❌ 查询状态失败: {response.status_code}")
                return False
        
        if wait_time >= max_wait_time:
            print("❌ 处理超时")
            return False
        
        # 5. 下载报告
        print("\n5. 下载处理报告...")
        response = requests.get(f"{BASE_URL}/submissions/{submission_id}/report")
        if response.status_code == 200:
            report_data = response.json()
            print("✅ 报告下载成功")
            print(f"   摘要: {report_data.get('summary', {})}")
            print(f"   基因映射: {report_data.get('gene_mapping', {})}")
            print(f"   细胞类型映射: {report_data.get('cell_type_mapping', {})}")
            
            # 保存报告到文件
            with open(f"test_report_{submission_id}.json", 'w') as f:
                json.dump(report_data, f, indent=2, ensure_ascii=False)
            print(f"   报告已保存: test_report_{submission_id}.json")
        else:
            print(f"❌ 报告下载失败: {response.status_code}")
            return False
        
        # 6. 下载处理结果
        print("\n6. 下载处理结果...")
        response = requests.get(f"{BASE_URL}/submissions/{submission_id}/export")
        if response.status_code == 200:
            with open(f"test_result_{submission_id}.h5ad", 'wb') as f:
                f.write(response.content)
            print("✅ 处理结果下载成功")
            print(f"   文件已保存: test_result_{submission_id}.h5ad")
        else:
            print(f"❌ 处理结果下载失败: {response.status_code}")
            return False
        
        print("\n🎉 所有测试通过!")
        return True
        
    except requests.exceptions.ConnectionError:
        print("❌ 无法连接到服务器，请确保服务已启动")
        return False
    except Exception as e:
        print(f"❌ 测试过程中发生错误: {e}")
        return False

def cleanup_test_files():
    """清理测试文件"""
    print("\n🧹 清理测试文件...")
    
    test_files = [
        TEST_FILE,
        "test_report_*.json",
        "test_result_*.h5ad"
    ]
    
    for pattern in test_files:
        if '*' in pattern:
            import glob
            files = glob.glob(pattern)
            for file in files:
                if os.path.exists(file):
                    try:
                        os.remove(file)
                        print(f"   删除: {file}")
                    except PermissionError:
                        print(f"   跳过（文件被占用）: {file}")
                    except Exception as e:
                        print(f"   删除失败: {file} - {e}")
        else:
            if os.path.exists(pattern):
                try:
                    os.remove(pattern)
                    print(f"   删除: {pattern}")
                except PermissionError:
                    print(f"   跳过（文件被占用）: {pattern}")
                except Exception as e:
                    print(f"   删除失败: {pattern} - {e}")

def main():
    """主函数"""
    print("🧪 scRNA-seq Data Processing Platform 测试脚本")
    print("=" * 50)
    
    # 检查服务是否运行
    try:
        response = requests.get(f"{BASE_URL}/", timeout=5)
        if response.status_code != 200:
            print("❌ 服务未正常运行")
            return False
    except requests.exceptions.ConnectionError:
        print("❌ 无法连接到服务器")
        print("请先启动服务: python main.py")
        return False
    
    # 创建测试数据
    if not create_test_data():
        return False
    
    # 运行测试
    success = test_api_endpoints()
    
    # 清理测试文件
    cleanup_test_files()
    
    if success:
        print("\n🎉 测试完成! 平台运行正常")
        return True
    else:
        print("\n❌ 测试失败! 请检查错误信息")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)

