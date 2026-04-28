#!/usr/bin/env python3
"""
检查数据库中的失败记录
"""

from database import get_db
from models import Submission
from datetime import datetime

def check_failed_submissions():
    """检查失败的提交记录"""
    print("🔍 检查数据库中的提交记录...")
    
    db = next(get_db())
    
    try:
        # 获取最近的提交记录
        submissions = db.query(Submission).order_by(Submission.created_at.desc()).limit(5).all()
        
        print(f"📊 找到 {len(submissions)} 条最近的提交记录:")
        
        for sub in submissions:
            print(f"\n📝 提交ID: {sub.id}")
            print(f"   状态: {sub.status}")
            print(f"   创建时间: {sub.created_at}")
            print(f"   更新时间: {sub.updated_at}")
            print(f"   文件路径: {sub.file_path}")
            
            if sub.error_message:
                print(f"   ❌ 错误信息: {sub.error_message}")
            
            if sub.report_path:
                print(f"   📋 报告路径: {sub.report_path}")
            
            if sub.export_path:
                print(f"   📤 导出路径: {sub.export_path}")
            
            print("-" * 50)
        
        # 检查失败状态的记录
        failed_subs = db.query(Submission).filter(
            Submission.status.in_(['failed', 'validation_failed'])
        ).order_by(Submission.created_at.desc()).limit(3).all()
        
        if failed_subs:
            print(f"\n❌ 找到 {len(failed_subs)} 条失败记录:")
            for sub in failed_subs:
                print(f"\n失败提交: {sub.id}")
                print(f"状态: {sub.status}")
                print(f"错误: {sub.error_message}")
                print("-" * 30)
        
    finally:
        db.close()

if __name__ == "__main__":
    check_failed_submissions()
