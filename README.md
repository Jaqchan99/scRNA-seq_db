<<<<<<< HEAD
# scRNA-seq Data Processing Platform

一个用于单细胞RNA测序数据标准化与前置处理的平台，支持数据入库前的准备和审核。

## 功能特性

- 🧬 **基因名称对齐**: 使用 mygene.info API 将基因符号统一转换为 Ensembl ID
- 🔬 **细胞类型标准化**: 使用 Zooma API 将细胞类型标签映射到 CL Ontology
- 📊 **数据校验**: 结构化和语义校验，确保数据质量
- 📤 **结果导出**: 生成标准化的 metadata 和表达矩阵
- 📋 **详细报告**: JSON + TSV 格式的处理报告

## 技术栈

- **后端**: FastAPI + SQLite
- **数据处理**: scanpy + pandas + numpy
- **基因映射**: mygene.info API
- **细胞类型映射**: Zooma API
- **文件格式**: h5ad (AnnData)

## 快速开始

### 1. 安装依赖

```bash
pip install -r requirements.txt
```

### 2. 启动平台（推荐）

使用一键启动脚本，同时启动后端API和前端界面：

```bash
python start_platform.py
```

这将启动：
- 后端API服务：`http://localhost:8000`
- 前端界面：`http://localhost:6666`

### 3. 手动启动

如果需要分别启动：

```bash
# 启动后端API服务
python main.py

# 在另一个终端启动前端（可选）
cd frontend
python -m http.server 6666
```

### 4. 使用平台

#### 通过Web界面（推荐）
1. 打开浏览器访问 `http://localhost:6666`
2. 按照界面提示上传 .h5ad 文件
3. 等待处理完成并下载结果

#### 使用 API

#### 创建提交任务
```bash
curl -X POST "http://localhost:8000/submissions"
```

#### 上传 h5ad 文件
```bash
curl -X POST "http://localhost:8000/submissions/{submission_id}/upload" \
     -F "file=@your_data.h5ad"
```

#### 查询处理状态
```bash
curl "http://localhost:8000/submissions/{submission_id}/status"
```

#### 下载处理报告
```bash
curl "http://localhost:8000/submissions/{submission_id}/report" -o report.json
```

#### 下载处理结果
```bash
curl "http://localhost:8000/submissions/{submission_id}/export" -o processed_data.h5ad
```

## API 接口

### POST /submissions
创建新的提交任务

**响应**:
```json
{
  "submission_id": "uuid",
  "status": "pending",
  "message": "提交任务已创建，请上传 h5ad 文件"
}
```

### POST /submissions/{submission_id}/upload
上传 h5ad 文件

**参数**:
- `file`: h5ad 文件

**响应**:
```json
{
  "submission_id": "uuid",
  "status": "processing",
  "message": "文件上传成功，正在处理中..."
}
```

### GET /submissions/{submission_id}/status
查询处理状态

**响应**:
```json
{
  "submission_id": "uuid",
  "status": "completed",
  "created_at": "2024-01-01T00:00:00",
  "updated_at": "2024-01-01T00:05:00"
}
```

### GET /submissions/{submission_id}/report
下载处理报告 (JSON 格式)

### GET /submissions/{submission_id}/export
下载处理结果 (h5ad 格式)

## 处理流程

1. **上传**: 用户上传 h5ad 文件
2. **校验**: 检查文件结构和基本语义
3. **基因映射**: 将基因符号映射到 Ensembl ID
4. **细胞类型标准化**: 将细胞类型标签映射到 CL Ontology
5. **导出**: 生成标准化的数据和报告

## 状态说明

- `pending`: 等待上传
- `uploading`: 正在上传
- `processing`: 正在处理
- `completed`: 处理完成
- `failed`: 处理失败
- `validation_failed`: 校验失败

## 输出文件

### 处理后的 h5ad 文件
- 包含原始数据 + 映射结果
- 新增列: `ensembl_gene_id_mapped`, `cl_id`, `cl_label`
- 处理信息保存在 `uns['processing_info']`

### 处理报告 (JSON)
```json
{
  "submission_id": "uuid",
  "summary": {
    "n_cells": 1000,
    "n_genes": 2000,
    "file_size_mb": 15.5
  },
  "gene_mapping": {
    "mapped_count": 1800,
    "unmapped_count": 200,
    "mapping_rate": 90.0
  },
  "cell_type_mapping": {
    "mapped_count": 8,
    "unmapped_count": 2,
    "mapping_rate": 80.0
  },
  "quality_metrics": {
    "sparsity": 0.95,
    "mean_genes_per_cell": 1200
  },
  "recommendations": [
    "建议检查基因符号格式",
    "建议人工审核细胞类型标签"
  ]
}
```

### 映射详情 (TSV)
- `gene_mapping_{submission_id}.tsv`: 基因映射详情
- `cell_type_mapping_{submission_id}.tsv`: 细胞类型映射详情

## 注意事项

1. **网络依赖**: 需要访问 mygene.info 和 Zooma API
2. **文件大小**: 建议单个文件不超过 1GB
3. **处理时间**: 取决于文件大小和网络状况
4. **数据格式**: 仅支持 h5ad 格式
5. **物种支持**: 当前主要支持小鼠数据

## 开发说明

### 项目结构
```
├── main.py                 # FastAPI 主应用
├── database.py            # 数据库配置和模型
├── models.py              # Pydantic 模型
├── upload_service.py      # 文件上传服务
├── validation_service.py  # 数据校验服务
├── mapping_service.py     # 基因和细胞类型映射服务
├── export_service.py      # 数据导出服务
├── start_platform.py     # 一键启动脚本
├── requirements.txt       # 依赖包
├── frontend/              # 前端界面
│   ├── index.html        # 主页面
│   ├── style.css         # 样式文件
│   ├── script.js         # JavaScript逻辑
│   └── README.md         # 前端说明
└── README.md             # 说明文档
```

### 数据库表
- `submissions`: 提交任务记录
- `mapping_logs`: 映射处理日志

### 目录结构
- `uploads/`: 原始上传文件
- `exports/`: 处理后的导出文件
- `reports/`: 处理报告
- `chunks/`: 分片上传临时文件

## 后续计划

- [ ] 支持更多文件格式 (MTX, Seurat RDS)
- [ ] 添加前端界面
- [ ] 支持批量处理
- [ ] 添加人工审核工作台
- [ ] 优化大文件处理性能
- [ ] 添加更多质量检查规则

## 许可证

MIT License

=======
# scRNA-seq_db
>>>>>>> bea7ad2d31ed4d25f52c0107b01f30ca3dab40c7
