import os
import hashlib
from typing import Dict, List, Any
import shutil

class UploadService:
    """文件上传服务"""
    
    def __init__(self):
        self.upload_dir = "uploads"
        self.chunk_dir = "chunks"
        
        # 确保目录存在
        os.makedirs(self.upload_dir, exist_ok=True)
        os.makedirs(self.chunk_dir, exist_ok=True)
    
    def save_chunk(self, submission_id: str, chunk_number: int, chunk_data: bytes) -> Dict[str, Any]:
        """保存文件分片"""
        chunk_path = os.path.join(self.chunk_dir, f"{submission_id}_chunk_{chunk_number}")
        
        with open(chunk_path, "wb") as f:
            f.write(chunk_data)
        
        return {
            "chunk_path": chunk_path,
            "chunk_size": len(chunk_data),
            "status": "saved"
        }
    
    def merge_chunks(self, submission_id: str, total_chunks: int, filename: str) -> Dict[str, Any]:
        """合并文件分片"""
        final_path = os.path.join(self.upload_dir, filename)
        
        with open(final_path, "wb") as final_file:
            for chunk_number in range(total_chunks):
                chunk_path = os.path.join(self.chunk_dir, f"{submission_id}_chunk_{chunk_number}")
                
                if os.path.exists(chunk_path):
                    with open(chunk_path, "rb") as chunk_file:
                        final_file.write(chunk_file.read())
                    
                    # 删除分片文件
                    os.remove(chunk_path)
                else:
                    raise FileNotFoundError(f"分片文件不存在: {chunk_path}")
        
        # 计算文件哈希
        file_hash = self._calculate_file_hash(final_path)
        file_size = os.path.getsize(final_path)
        
        return {
            "file_path": final_path,
            "file_size": file_size,
            "file_hash": file_hash,
            "status": "merged"
        }
    
    def _calculate_file_hash(self, file_path: str) -> str:
        """计算文件 MD5 哈希"""
        hash_md5 = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()
    
    def validate_file(self, file_path: str) -> Dict[str, Any]:
        """验证上传的文件"""
        if not os.path.exists(file_path):
            return {
                "valid": False,
                "error": "文件不存在"
            }
        
        # 检查文件大小
        file_size = os.path.getsize(file_path)
        max_size = 1024 * 1024 * 1024  # 1GB
        
        if file_size > max_size:
            return {
                "valid": False,
                "error": f"文件过大: {file_size / (1024*1024):.1f}MB > {max_size / (1024*1024):.1f}MB"
            }
        
        # 检查文件扩展名
        if not file_path.endswith('.h5ad'):
            return {
                "valid": False,
                "error": "仅支持 .h5ad 文件"
            }
        
        return {
            "valid": True,
            "file_size": file_size,
            "file_hash": self._calculate_file_hash(file_path)
        }
    
    def cleanup_chunks(self, submission_id: str):
        """清理分片文件"""
        chunk_files = [f for f in os.listdir(self.chunk_dir) if f.startswith(f"{submission_id}_chunk_")]
        
        for chunk_file in chunk_files:
            chunk_path = os.path.join(self.chunk_dir, chunk_file)
            if os.path.exists(chunk_path):
                os.remove(chunk_path)
    
    def get_upload_progress(self, submission_id: str, total_chunks: int) -> Dict[str, Any]:
        """获取上传进度"""
        uploaded_chunks = 0
        
        for chunk_number in range(total_chunks):
            chunk_path = os.path.join(self.chunk_dir, f"{submission_id}_chunk_{chunk_number}")
            if os.path.exists(chunk_path):
                uploaded_chunks += 1
        
        progress = uploaded_chunks / total_chunks if total_chunks > 0 else 0
        
        return {
            "uploaded_chunks": uploaded_chunks,
            "total_chunks": total_chunks,
            "progress": progress,
            "percentage": round(progress * 100, 2)
        }

