from fastapi import FastAPI, UploadFile, File, HTTPException, BackgroundTasks
from fastapi.responses import FileResponse, JSONResponse
from fastapi.middleware.cors import CORSMiddleware
import uvicorn
import os
import uuid
from datetime import datetime
from typing import Optional
import json

from database import get_db, init_db, Submission, MappingLog
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
    
    # 调试日志: 写入字节数与磁盘大小
    try:
        on_disk_size = os.path.getsize(file_path)
        print(f"📝 上传调试: wrote_bytes={len(content)}, on_disk_size={on_disk_size}, path={file_path}")
    except Exception as e:
        print(f"📝 上传调试: 获取文件大小失败: {e}")
    
    submission.file_path = file_path
    submission.status = "processing"
    db.commit()
    db.close()
    
    # 不再自动处理，等待用户触发
    submission.status = "uploaded"
    db.commit()
    db.close()
    
    return {
        "submission_id": submission_id,
        "status": "uploaded",
        "message": "文件上传成功，请开始校验"
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

@app.post("/submissions/{submission_id}/validate")
async def validate_submission(submission_id: str, background_tasks: BackgroundTasks):
    """执行数据校验"""
    db = next(get_db())
    submission = db.query(Submission).filter(Submission.id == submission_id).first()
    if not submission:
        raise HTTPException(status_code=404, detail="提交任务不存在")
    
    submission.status = "validating"
    db.commit()
    db.close()
    
    background_tasks.add_task(process_validation, submission_id)
    
    return {"status": "validating", "message": "开始数据校验"}

@app.post("/submissions/{submission_id}/continue")
async def continue_processing(submission_id: str, background_tasks: BackgroundTasks):
    """继续处理（执行映射和导出）"""
    db = next(get_db())
    submission = db.query(Submission).filter(Submission.id == submission_id).first()
    if not submission:
        raise HTTPException(status_code=404, detail="提交任务不存在")
    
    if submission.status != "validated":
        raise HTTPException(status_code=400, detail="请先完成数据校验")
    
    submission.status = "mapping"
    db.commit()
    db.close()
    
    background_tasks.add_task(process_mapping_and_export, submission_id)
    
    return {"status": "mapping", "message": "开始基因映射和细胞类型标准化"}

def process_validation(submission_id: str):
    """只执行数据校验"""
    print(f"🔍 开始校验提交 {submission_id}")
    
    try:
        db = next(get_db())
        submission = db.query(Submission).filter(Submission.id == submission_id).first()
        
        if not submission:
            return
        
        if not os.path.exists(submission.file_path):
            submission.status = "failed"
            submission.error_message = f"文件不存在: {submission.file_path}"
            db.commit()
            db.close()
            return
        
        # 执行校验
        validation_result = validation_service.validate_h5ad(submission.file_path)
        
        if not validation_result["valid"]:
            submission.status = "validation_failed"
            submission.error_message = str(validation_result["errors"])
            db.commit()
            db.close()
            return
        
        # 保存校验结果到临时文件
        import json
        validation_report_path = f"reports/validation_{submission_id}.json"
        os.makedirs("reports", exist_ok=True)
        with open(validation_report_path, 'w') as f:
            json.dump({
                "validation_result": validation_result,
                "summary": validation_result.get("metadata", {})
            }, f, indent=2, ensure_ascii=False)
        
        submission.status = "validated"
        submission.report_path = validation_report_path
        db.commit()
        db.close()
        
        print(f"✅ 校验完成: {submission_id}")
        
    except Exception as e:
        print(f"❌ 校验出错: {str(e)}")
        import traceback
        traceback.print_exc()
        try:
            db = next(get_db())
            submission = db.query(Submission).filter(Submission.id == submission_id).first()
            if submission:
                submission.status = "failed"
                submission.error_message = str(e)
                db.commit()
            db.close()
        except:
            pass

def process_mapping_and_export(submission_id: str):
    """执行映射和导出"""
    print(f"🧬 开始映射和导出 {submission_id}")
    
    try:
        db = next(get_db())
        submission = db.query(Submission).filter(Submission.id == submission_id).first()
        
        if not submission:
            return
        
        # 基因映射
        print(f"🧬 开始基因映射 {submission_id}")
        gene_mapping_result = mapping_service.map_genes(submission.file_path)
        
        # 细胞类型映射
        print(f"🔬 开始细胞类型标准化 {submission_id}")
        cell_type_result = mapping_service.map_cell_types(submission.file_path)
        
        # 导出结果
        print(f"📤 开始导出结果 {submission_id}")
        export_result = export_service.export_processed_data(
            submission_id, 
            submission.file_path,
            gene_mapping_result,
            cell_type_result
        )
        
        submission.status = "completed"
        submission.report_path = export_result["report_path"]
        submission.export_path = export_result["export_path"]
        submission.updated_at = datetime.now()
        db.commit()
        db.close()
        
        print(f"✅ 处理完成: {submission_id}")
        
    except Exception as e:
        print(f"❌ 处理出错: {str(e)}")
        import traceback
        traceback.print_exc()
        try:
            db = next(get_db())
            submission = db.query(Submission).filter(Submission.id == submission_id).first()
            if submission:
                submission.status = "failed"
                submission.error_message = str(e)
                db.commit()
            db.close()
        except:
            pass

def process_submission(submission_id: str):
    """后台处理提交的数据"""
    print(f"🚀 开始处理提交 {submission_id}")
    
    try:
        db = next(get_db())
        submission = db.query(Submission).filter(Submission.id == submission_id).first()
        
        if not submission:
            print(f"❌ 找不到提交记录: {submission_id}")
            return
        
        print(f"📁 处理文件: {submission.file_path}")
        
        # 调试日志: 处理前检查文件存在性与大小
        try:
            exists = os.path.exists(submission.file_path)
            size = os.path.getsize(submission.file_path) if exists else -1
            print(f"📝 处理调试: exists={exists}, size={size}, path={submission.file_path}")
        except Exception as e:
            print(f"📝 处理调试: 检查文件失败: {e}")
        
        # 检查文件是否存在
        if not os.path.exists(submission.file_path):
            print(f"❌ 文件不存在: {submission.file_path}")
            submission.status = "failed"
            submission.error_message = f"文件不存在: {submission.file_path}"
            db.commit()
            db.close()
            return
        
        # 1. 结构校验
        print(f"🔍 开始校验提交 {submission_id}")
        try:
            validation_result = validation_service.validate_h5ad(submission.file_path)
            print(f"✅ 校验完成: {validation_result}")
        except Exception as e:
            print(f"❌ 校验过程出错: {e}")
            submission.status = "failed"
            submission.error_message = f"校验过程出错: {str(e)}"
            db.commit()
            db.close()
            return
        
        if not validation_result["valid"]:
            print(f"❌ 校验失败: {validation_result['errors']}")
            submission.status = "validation_failed"
            submission.error_message = str(validation_result["errors"])
            db.commit()
            db.close()
            return
        
        # 2. 基因映射
        print(f"🧬 开始基因映射 {submission_id}")
        try:
            gene_mapping_result = mapping_service.map_genes(submission.file_path)
            print(f"✅ 基因映射完成: {gene_mapping_result}")
        except Exception as e:
            print(f"❌ 基因映射出错: {e}")
            gene_mapping_result = {"error": str(e)}
        
        # 3. 细胞类型标准化
        print(f"🔬 开始细胞类型标准化 {submission_id}")
        try:
            cell_type_result = mapping_service.map_cell_types(submission.file_path)
            print(f"✅ 细胞类型标准化完成: {cell_type_result}")
        except Exception as e:
            print(f"❌ 细胞类型标准化出错: {e}")
            cell_type_result = {"error": str(e)}
        
        # 4. 导出结果
        print(f"📤 开始导出结果 {submission_id}")
        try:
            export_result = export_service.export_processed_data(
                submission_id, 
                submission.file_path,
                gene_mapping_result,
                cell_type_result
            )
            print(f"✅ 导出完成: {export_result}")
        except Exception as e:
            print(f"❌ 导出出错: {e}")
            submission.status = "failed"
            submission.error_message = f"导出过程出错: {str(e)}"
            db.commit()
            db.close()
            return
        
        # 更新状态
        submission.status = "completed"
        submission.report_path = export_result["report_path"]
        submission.export_path = export_result["export_path"]
        submission.updated_at = datetime.now()
        db.commit()
        db.close()
        
        print(f"✅ 提交 {submission_id} 处理完成")
        
    except Exception as e:
        print(f"❌ 处理提交 {submission_id} 时出错: {str(e)}")
        import traceback
        traceback.print_exc()
        
        try:
            db = next(get_db())
            submission = db.query(Submission).filter(Submission.id == submission_id).first()
            if submission:
                submission.status = "failed"
                submission.error_message = str(e)
                db.commit()
            db.close()
        except Exception as db_error:
            print(f"❌ 更新数据库状态失败: {db_error}")

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
    uvicorn.run(app, host="0.0.0.0", port=8100)

