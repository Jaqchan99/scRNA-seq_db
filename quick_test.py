#!/usr/bin/env python3
"""
快速测试脚本 - 测试修复后的处理逻辑
"""

import requests
import time
import json
import os

BASE_URL = "http://localhost:8000"

def create_minimal_test_data():
    """创建最小的测试数据"""
    print("🧪 创建最小测试数据...")
    
    try:
        import scanpy as sc
        import pandas as pd
        import numpy as np
        
        # 创建最小的测试数据
        n_cells = 3
        n_genes = 5
        
        # 创建表达矩阵
        X = np.array([
            [1, 2, 3, 4, 5],
            [2, 3, 4, 5, 6],
            [3, 4, 5, 6, 7]
        ], dtype=float)
        
        # 创建细胞元数据
        obs = pd.DataFrame({
            'cell_id': ['cell_1', 'cell_2', 'cell_3'],
            'raw_cell_type_label': ['T cell', 'B cell', 'NK cell'],
            'n_genes': [200, 300, 250],
            'n_counts': [1000, 1500, 1200],
            'pct_mt': [10, 15, 12],
            'batch': ['batch1', 'batch1', 'batch2']
        })
        
        # 创建基因元数据
        var = pd.DataFrame({
            'gene_symbol': ['Gene1', 'Gene2', 'Gene3', 'Gene4', 'Gene5'],
            'ensembl_gene_id': ['ENSMUSG00000000001', 'ENSMUSG00000000002', 'ENSMUSG00000000003', 'ENSMUSG00000000004', 'ENSMUSG00000000005']
        })
        
        # 创建 AnnData 对象
        adata = sc.AnnData(X=X, obs=obs, var=var)
        adata.uns['organism'] = 'mouse'
        adata.uns['genome_build'] = 'GRCm39'
        
        # 保存测试文件
        adata.write_h5ad("minimal_test.h5ad")
        print("✅ 最小测试数据创建成功: minimal_test.h5ad")
        
        # 确保文件句柄关闭
        del adata
        
        return True
        
    except Exception as e:
        print(f"❌ 创建测试数据失败: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_complete_flow():
    """测试完整流程"""
    print("🚀 测试完整处理流程")
    print("=" * 40)
    
    # 1. 创建提交任务
    print("1. 创建提交任务...")
    response = requests.post(f"{BASE_URL}/submissions")
    if response.status_code != 200:
        print(f"❌ 创建提交失败: {response.text}")
        return False
    
    submission_id = response.json()["submission_id"]
    print(f"✅ 提交ID: {submission_id}")
    
    # 2. 上传文件
    print("\n2. 上传测试文件...")
    if not os.path.exists("minimal_test.h5ad"):
        print("❌ 测试文件不存在")
        return False
    
    with open("minimal_test.h5ad", "rb") as f:
        files = {"file": ("minimal_test.h5ad", f, "application/octet-stream")}
        response = requests.post(
            f"{BASE_URL}/submissions/{submission_id}/upload",
            files=files
        )
    
    print(f"上传响应状态码: {response.status_code}")
    if response.status_code != 200:
        print(f"❌ 上传失败: {response.text}")
        return False
    
    print("✅ 文件上传成功")
    
    # 3. 监控处理过程
    print("\n3. 监控处理过程...")
    max_wait = 120  # 2分钟
    wait_time = 0
    check_interval = 5
    
    while wait_time < max_wait:
        response = requests.get(f"{BASE_URL}/submissions/{submission_id}/status")
        if response.status_code == 200:
            data = response.json()
            status = data["status"]
            print(f"   状态: {status} (等待 {wait_time}s)")
            
            if status == "completed":
                print("✅ 处理完成!")
                
                # 4. 下载报告
                print("\n4. 下载处理报告...")
                response = requests.get(f"{BASE_URL}/submissions/{submission_id}/report")
                if response.status_code == 200:
                    report_data = response.json()
                    print("✅ 报告下载成功")
                    print(f"   摘要: {json.dumps(report_data.get('summary', {}), indent=2)}")
                else:
                    print(f"❌ 报告下载失败: {response.status_code}")
                
                return True
                
            elif status in ["failed", "validation_failed"]:
                print(f"❌ 处理失败: {status}")
                return False
        
        time.sleep(check_interval)
        wait_time += check_interval
    
    print("❌ 处理超时")
    return False

def main():
    """主函数"""
    print("🧪 快速测试修复后的处理逻辑")
    print("=" * 50)
    
    # 创建测试数据
    if not create_minimal_test_data():
        return False
    
    # 测试完整流程
    success = test_complete_flow()
    
    # 清理测试文件
    try:
        if os.path.exists("minimal_test.h5ad"):
            os.remove("minimal_test.h5ad")
            print("\n🧹 清理测试文件")
    except Exception as e:
        print(f"清理文件时出错: {e}")
    
    if success:
        print("\n🎉 测试成功! 处理逻辑已修复")
        return True
    else:
        print("\n❌ 测试失败! 需要进一步调试")
        return False

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
