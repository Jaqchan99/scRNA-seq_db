#!/usr/bin/env python3
"""
超简单调试脚本
"""

import requests
import json

BASE_URL = "http://localhost:8000"

def test_connection():
    """测试连接"""
    print("🔍 测试服务器连接...")
    
    try:
        response = requests.get(f"{BASE_URL}/", timeout=5)
        print(f"状态码: {response.status_code}")
        print(f"响应: {response.json()}")
        return True
    except Exception as e:
        print(f"连接失败: {e}")
        return False

def test_create_submission():
    """测试创建提交"""
    print("\n📝 测试创建提交...")
    
    try:
        response = requests.post(f"{BASE_URL}/submissions", timeout=10)
        print(f"状态码: {response.status_code}")
        
        if response.status_code == 200:
            data = response.json()
            submission_id = data["submission_id"]
            print(f"✅ 提交ID: {submission_id}")
            return submission_id
        else:
            print(f"❌ 失败: {response.text}")
            return None
    except Exception as e:
        print(f"异常: {e}")
        return None

def test_status(submission_id):
    """测试状态查询"""
    print(f"\n📊 测试状态查询 (ID: {submission_id})...")
    
    try:
        response = requests.get(f"{BASE_URL}/submissions/{submission_id}/status", timeout=10)
        print(f"状态码: {response.status_code}")
        
        if response.status_code == 200:
            data = response.json()
            print(f"✅ 状态: {data}")
            return True
        else:
            print(f"❌ 失败: {response.text}")
            return False
    except Exception as e:
        print(f"异常: {e}")
        return False

if __name__ == "__main__":
    print("🧪 简单调试测试")
    print("=" * 40)
    
    # 测试连接
    if not test_connection():
        exit(1)
    
    # 测试创建提交
    submission_id = test_create_submission()
    if not submission_id:
        exit(1)
    
    # 测试状态查询
    test_status(submission_id)
    
    print("\n✅ 基本测试完成")
