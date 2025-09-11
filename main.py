from fastapi import FastAPI, UploadFile, File, HTTPException, BackgroundTasks
from fastapi.responses import FileResponse, JSONResponse
from fastapi.middleware.cors import CORSMiddleware
import uvicorn
import os
import uuid
from datetime import datetime
from typing import Optional
import json

from database import get_db, init_db
from models import Submission, MappingLog
from upload_service import UploadService
from validation_service import ValidationService
from mapping_service import MappingService
from export_service import ExportService

app = FastAPI(title="scRNA-seq Data Processing Platform", version="1.0.0")

# CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# 创建必要的目录
os.makedirs("uploads", exist_ok=True)
os.makedirs("exports", exist_ok=True)
os.makedirs("reports", exist_ok=True)

# 初始化数据库
init_db()

# 服务实例
upload_service = UploadService()
validation_service = ValidationService()
mapping_service = MappingService()
export_service = ExportService()

@app.on_event("startup")
async def startup_event():
    """应用启动时初始化"""
    print("🚀 scRNA-seq Data Processing Platform started!")

@app.post("/submissions")
async def create_submission():
    """创建新的提交任务"""
    submission_id = str(uuid.uuid4())
    
    # 创建数据库记录
    db = next(get_db())
    submission = Submission(
        id=submission_id,
        status="pending",
        created_at=datetime.now(),
        updated_at=datetime.now()
    )
    db.add(submission)
    db.commit()
    db.close()
    
    return {
        "submission_id": submission_id,
        "status": "pending",
        "message": "提交任务已创建，请上传 h5ad 文件"
    }

@app.post("/submissions/{submission_id}/upload")
async def upload_file(
    submission_id: str,
    file: UploadFile = File(...),
    background_tasks: BackgroundTasks = BackgroundTasks()
):
    """上传 h5ad 文件并开始处理"""
    
    # 验证文件类型
    if not file.filename.endswith('.h5ad'):
        raise HTTPException(status_code=400, detail="仅支持 .h5ad 文件")
    
    # 更新数据库状态
    db = next(get_db())
    submission = db.query(Submission).filter(Submission.id == submission_id).first()
    if not submission:
        raise HTTPException(status_code=404, detail="提交任务不存在")
    
    submission.status = "uploading"
    submission.updated_at = datetime.now()
    db.commit()
    
    # 保存文件
    file_path = f"uploads/{submission_id}.h5ad"
    with open(file_path, "wb") as buffer:
        content = await file.read()
        buffer.write(content)
    
    submission.file_path = file_path
    submission.status = "processing"
    db.commit()
    db.close()
    
    # 后台处理
    background_tasks.add_task(process_submission, submission_id)
    
    return {
        "submission_id": submission_id,
        "status": "processing",
        "message": "文件上传成功，正在处理中..."
    }

@app.get("/submissions/{submission_id}/status")
async def get_status(submission_id: str):
    """查询处理状态"""
    db = next(get_db())
    submission = db.query(Submission).filter(Submission.id == submission_id).first()
    if not submission:
        raise HTTPException(status_code=404, detail="提交任务不存在")
    
    status_info = {
        "submission_id": submission_id,
        "status": submission.status,
        "created_at": submission.created_at.isoformat(),
        "updated_at": submission.updated_at.isoformat()
    }
    
    db.close()
    return status_info

@app.get("/submissions/{submission_id}/report")
async def get_report(submission_id: str):
    """下载处理报告"""
    db = next(get_db())
    submission = db.query(Submission).filter(Submission.id == submission_id).first()
    if not submission:
        raise HTTPException(status_code=404, detail="提交任务不存在")
    
    if not submission.report_path or not os.path.exists(submission.report_path):
        raise HTTPException(status_code=404, detail="报告尚未生成")
    
    db.close()
    return FileResponse(
        submission.report_path,
        media_type="application/json",
        filename=f"report_{submission_id}.json"
    )

@app.get("/submissions/{submission_id}/export")
async def get_export(submission_id: str):
    """下载处理结果"""
    db = next(get_db())
    submission = db.query(Submission).filter(Submission.id == submission_id).first()
    if not submission:
        raise HTTPException(status_code=404, detail="提交任务不存在")
    
    if not submission.export_path or not os.path.exists(submission.export_path):
        raise HTTPException(status_code=404, detail="导出文件尚未生成")
    
    db.close()
    return FileResponse(
        submission.export_path,
        media_type="application/octet-stream",
        filename=f"processed_{submission_id}.h5ad"
    )

async def process_submission(submission_id: str):
    """后台处理提交的数据"""
    try:
        db = next(get_db())
        submission = db.query(Submission).filter(Submission.id == submission_id).first()
        
        if not submission:
            return
        
        # 1. 结构校验
        print(f"🔍 开始校验提交 {submission_id}")
        validation_result = validation_service.validate_h5ad(submission.file_path)
        
        if not validation_result["valid"]:
            submission.status = "validation_failed"
            submission.error_message = validation_result["errors"]
            db.commit()
            db.close()
            return
        
        # 2. 基因映射
        print(f"🧬 开始基因映射 {submission_id}")
        gene_mapping_result = mapping_service.map_genes(submission.file_path)
        
        # 3. 细胞类型标准化
        print(f"🔬 开始细胞类型标准化 {submission_id}")
        cell_type_result = mapping_service.map_cell_types(submission.file_path)
        
        # 4. 导出结果
        print(f"📤 开始导出结果 {submission_id}")
        export_result = export_service.export_processed_data(
            submission_id, 
            submission.file_path,
            gene_mapping_result,
            cell_type_result
        )
        
        # 更新状态
        submission.status = "completed"
        submission.report_path = export_result["report_path"]
        submission.export_path = export_result["export_path"]
        submission.updated_at = datetime.now()
        db.commit()
        
        print(f"✅ 提交 {submission_id} 处理完成")
        
    except Exception as e:
        print(f"❌ 处理提交 {submission_id} 时出错: {str(e)}")
        db = next(get_db())
        submission = db.query(Submission).filter(Submission.id == submission_id).first()
        if submission:
            submission.status = "failed"
            submission.error_message = str(e)
            db.commit()
        db.close()

@app.get("/")
async def root():
    """根路径，返回 API 信息"""
    return {
        "message": "scRNA-seq Data Processing Platform",
        "version": "1.0.0",
        "status": "running",
        "endpoints": {
            "create_submission": "POST /submissions",
            "upload_file": "POST /submissions/{id}/upload",
            "get_status": "GET /submissions/{id}/status",
            "get_report": "GET /submissions/{id}/report",
            "get_export": "GET /submissions/{id}/export"
        }
    }

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)

