#!/usr/bin/env python3
"""
调试测试脚本 - 获取详细错误信息
"""

import requests
import time
import json
import os
import sys
import traceback

# 测试配置
BASE_URL = "http://localhost:8000"
TEST_FILE = "test_data.h5ad"

def create_test_data():
    """创建测试用的 h5ad 文件"""
    try:
        import scanpy as sc
        import pandas as pd
        import numpy as np
        
        print("🧪 创建测试数据...")
        
        # 创建模拟的单细胞数据
        n_cells = 10  # 减少数据量，加快测试
        n_genes = 50
        
        # 创建表达矩阵
        X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(float)
        
        # 创建细胞元数据
        obs = pd.DataFrame({
            'cell_id': [f'cell_{i}' for i in range(n_cells)],
            'raw_cell_type_label': np.random.choice([
                'T cell', 'B cell', 'NK cell', 'Monocyte', 'Dendritic cell'
            ], n_cells),
            'n_genes': np.random.randint(200, 1000, n_cells),
            'n_counts': np.random.randint(1000, 5000, n_cells),
            'pct_mt': np.random.uniform(5, 25, n_cells),
            'batch': np.random.choice(['batch1', 'batch2'], n_cells)
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
        print(f"   文件大小: {os.path.getsize(TEST_FILE)} bytes")
        
        return True
        
    except Exception as e:
        print(f"❌ 创建测试数据失败: {e}")
        traceback.print_exc()
        return False

def test_basic_api():
    """测试基本 API 功能"""
    print("\n🚀 测试基本 API...")
    
    try:
        # 1. 测试根端点
        print("1. 测试根端点...")
        response = requests.get(f"{BASE_URL}/", timeout=10)
        print(f"   状态码: {response.status_code}")
        if response.status_code == 200:
            print("✅ 根端点正常")
            print(f"   响应: {response.json()}")
        else:
            print(f"❌ 根端点失败: {response.text}")
            return False, None
        
        # 2. 创建提交任务
        print("\n2. 创建提交任务...")
        response = requests.post(f"{BASE_URL}/submissions", timeout=10)
        print(f"   状态码: {response.status_code}")
        if response.status_code == 200:
            submission_data = response.json()
            submission_id = submission_data["submission_id"]
            print(f"✅ 提交任务已创建: {submission_id}")
            return True, submission_id
        else:
            print(f"❌ 创建提交任务失败: {response.text}")
            return False, None
        
    except Exception as e:
        print(f"❌ API 测试失败: {e}")
        traceback.print_exc()
        return False, None

def test_file_upload(submission_id):
    """测试文件上传"""
    print(f"\n📤 测试文件上传 (ID: {submission_id})...")
    
    try:
        if not os.path.exists(TEST_FILE):
            print(f"❌ 测试文件不存在: {TEST_FILE}")
            return False
        
        print(f"   上传文件: {TEST_FILE}")
        print(f"   文件大小: {os.path.getsize(TEST_FILE)} bytes")
        
        with open(TEST_FILE, 'rb') as f:
            files = {'file': (TEST_FILE, f, 'application/octet-stream')}
            response = requests.post(
                f"{BASE_URL}/submissions/{submission_id}/upload",
                files=files,
                timeout=30
            )
        
        print(f"   状态码: {response.status_code}")
        if response.status_code == 200:
            print("✅ 文件上传成功")
            print(f"   响应: {response.json()}")
            return True
        else:
            print(f"❌ 文件上传失败: {response.text}")
            return False
        
    except Exception as e:
        print(f"❌ 文件上传异常: {e}")
        traceback.print_exc()
        return False

def monitor_processing(submission_id):
    """监控处理过程"""
    print(f"\n👀 监控处理过程 (ID: {submission_id})...")
    
    max_wait_time = 120  # 2分钟
    wait_time = 0
    check_interval = 5   # 5秒检查一次
    
    while wait_time < max_wait_time:
        try:
            response = requests.get(f"{BASE_URL}/submissions/{submission_id}/status", timeout=10)
            if response.status_code == 200:
                status_data = response.json()
                status = status_data["status"]
                print(f"   当前状态: {status} (等待 {wait_time}s)")
                
                if status == "completed":
                    print("✅ 处理完成!")
                    return True, status_data
                elif status in ["failed", "validation_failed"]:
                    print(f"❌ 处理失败: {status}")
                    print(f"   状态数据: {status_data}")
                    return False, status_data
                
                time.sleep(check_interval)
                wait_time += check_interval
            else:
                print(f"❌ 查询状态失败: {response.status_code} - {response.text}")
                time.sleep(check_interval)
                wait_time += check_interval
                
        except Exception as e:
            print(f"❌ 状态查询异常: {e}")
            time.sleep(check_interval)
            wait_time += check_interval
    
    print("❌ 处理超时")
    return False, None

def main():
    """主函数"""
    print("🔍 scRNA-seq Data Processing Platform 调试测试")
    print("=" * 60)
    
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
    
    # 测试基本 API
    success, submission_id = test_basic_api()
    if not success:
        return False
    
    # 测试文件上传
    if not test_file_upload(submission_id):
        return False
    
    # 监控处理过程
    success, status_data = monitor_processing(submission_id)
    
    if success:
        print("\n🎉 测试成功完成!")
        
        # 尝试下载报告
        try:
            print("\n📋 下载处理报告...")
            response = requests.get(f"{BASE_URL}/submissions/{submission_id}/report", timeout=10)
            if response.status_code == 200:
                report_data = response.json()
                print("✅ 报告下载成功")
                print(f"   摘要: {json.dumps(report_data.get('summary', {}), indent=2)}")
                
                # 保存报告
                with open(f"debug_report_{submission_id}.json", 'w') as f:
                    json.dump(report_data, f, indent=2, ensure_ascii=False)
                print(f"   报告已保存: debug_report_{submission_id}.json")
            else:
                print(f"❌ 报告下载失败: {response.status_code} - {response.text}")
        except Exception as e:
            print(f"❌ 报告下载异常: {e}")
    else:
        print(f"\n❌ 测试失败! 状态: {status_data}")
    
    # 清理测试文件
    try:
        if os.path.exists(TEST_FILE):
            os.remove(TEST_FILE)
            print(f"\n🧹 清理测试文件: {TEST_FILE}")
    except:
        pass
    
    return success

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
