from pydantic import BaseModel
from typing import Optional, List, Dict, Any
from datetime import datetime

class SubmissionCreate(BaseModel):
    """创建提交任务的请求模型"""
    pass

class SubmissionResponse(BaseModel):
    """提交任务响应模型"""
    submission_id: str
    status: str
    message: str

class UploadResponse(BaseModel):
    """上传响应模型"""
    submission_id: str
    status: str
    message: str

class StatusResponse(BaseModel):
    """状态查询响应模型"""
    submission_id: str
    status: str
    created_at: str
    updated_at: str

class ValidationResult(BaseModel):
    """校验结果模型"""
    valid: bool
    errors: List[str]
    warnings: List[str]
    metadata: Dict[str, Any]

class GeneMappingResult(BaseModel):
    """基因映射结果模型"""
    mapped_count: int
    unmapped_count: int
    ambiguous_count: int
    mapping_details: List[Dict[str, Any]]

class CellTypeMappingResult(BaseModel):
    """细胞类型映射结果模型"""
    mapped_count: int
    unmapped_count: int
    ambiguous_count: int
    mapping_details: List[Dict[str, Any]]

class ProcessingReport(BaseModel):
    """处理报告模型"""
    submission_id: str
    status: str
    validation_result: ValidationResult
    gene_mapping_result: GeneMappingResult
    cell_type_mapping_result: CellTypeMappingResult
    processing_time: float
    created_at: str

