import scanpy as sc
import pandas as pd
import numpy as np
from typing import Dict, List, Any
import os

class ValidationService:
    """数据校验服务"""
    
    def __init__(self):
        self.required_obs_columns = [
            'cell_id', 'raw_cell_type_label', 'n_genes', 'n_counts', 'pct_mt'
        ]
        self.required_var_columns = [
            'gene_symbol', 'ensembl_gene_id'
        ]
        self.valid_organisms = ['mouse', 'Mus musculus', 'mm10', 'mm39']
        self.valid_genome_builds = ['GRCm38', 'GRCm39', 'mm10', 'mm39']
    
    def validate_h5ad(self, file_path: str) -> Dict[str, Any]:
        """校验 h5ad 文件"""
        errors = []
        warnings = []
        metadata = {}
        
        try:
            # 检查文件是否存在
            if not os.path.exists(file_path):
                errors.append(f"文件不存在: {file_path}")
                return {"valid": False, "errors": errors, "warnings": warnings, "metadata": metadata}
            
            # 尝试读取 h5ad 文件
            try:
                adata = sc.read_h5ad(file_path)
            except Exception as e:
                errors.append(f"无法读取 h5ad 文件: {str(e)}")
                return {"valid": False, "errors": errors, "warnings": warnings, "metadata": metadata}
            
            # 基本结构检查
            if adata.n_obs == 0:
                errors.append("表达矩阵中没有细胞数据")
            if adata.n_vars == 0:
                errors.append("表达矩阵中没有基因数据")
            
            # 检查 obs 列
            missing_obs_cols = [col for col in self.required_obs_columns if col not in adata.obs.columns]
            if missing_obs_cols:
                warnings.append(f"缺少建议的 obs 列: {missing_obs_cols}")
            
            # 检查 var 列
            missing_var_cols = [col for col in self.required_var_columns if col not in adata.var.columns]
            if missing_var_cols:
                warnings.append(f"缺少建议的 var 列: {missing_var_cols}")
            
            # 检查基因符号列
            if 'gene_symbol' in adata.var.columns:
                gene_symbols = adata.var['gene_symbol'].dropna()
                if len(gene_symbols) == 0:
                    warnings.append("基因符号列为空")
                else:
                    # 检查是否有重复的基因符号
                    duplicates = gene_symbols[gene_symbols.duplicated()]
                    if len(duplicates) > 0:
                        warnings.append(f"发现重复的基因符号: {len(duplicates)} 个")
            
            # 检查细胞类型标签
            if 'raw_cell_type_label' in adata.obs.columns:
                cell_types = adata.obs['raw_cell_type_label'].dropna()
                if len(cell_types) == 0:
                    warnings.append("细胞类型标签为空")
                else:
                    unique_types = cell_types.unique()
                    metadata['unique_cell_types'] = len(unique_types)
                    metadata['cell_type_list'] = unique_types.tolist()[:10]  # 只显示前10个
            
            # 检查表达矩阵
            if hasattr(adata, 'X') and adata.X is not None:
                metadata['n_cells'] = adata.n_obs
                metadata['n_genes'] = adata.n_vars
                
                # 检查稀疏矩阵
                if hasattr(adata.X, 'nnz'):
                    metadata['sparsity'] = 1 - (adata.X.nnz / (adata.n_obs * adata.n_vars))
                
                # 检查是否有负值
                if hasattr(adata.X, 'data'):
                    if np.any(adata.X.data < 0):
                        warnings.append("表达矩阵包含负值")
            
            # 检查全局元数据
            if hasattr(adata, 'uns') and adata.uns:
                metadata['uns_keys'] = list(adata.uns.keys())
                
                # 检查 organism 信息
                organism_info = None
                for key in ['organism', 'species', 'genome_build']:
                    if key in adata.uns:
                        organism_info = adata.uns[key]
                        break
                
                if organism_info:
                    metadata['organism_info'] = organism_info
                    if isinstance(organism_info, str):
                        organism_lower = organism_info.lower()
                        if not any(valid in organism_lower for valid in ['mouse', 'mus']):
                            warnings.append(f"检测到的物种可能不是小鼠: {organism_info}")
                else:
                    warnings.append("未找到物种信息，假设为小鼠")
            
            # 检查批次信息
            batch_cols = [col for col in adata.obs.columns if 'batch' in col.lower()]
            if batch_cols:
                metadata['batch_columns'] = batch_cols
                for col in batch_cols:
                    unique_batches = adata.obs[col].nunique()
                    metadata[f'{col}_unique_count'] = unique_batches
            
            # 检查质量指标
            quality_cols = ['n_genes', 'n_counts', 'pct_mt']
            for col in quality_cols:
                if col in adata.obs.columns:
                    values = adata.obs[col].dropna()
                    if len(values) > 0:
                        metadata[f'{col}_stats'] = {
                            'mean': float(values.mean()),
                            'median': float(values.median()),
                            'min': float(values.min()),
                            'max': float(values.max())
                        }
            
            metadata['file_size_mb'] = os.path.getsize(file_path) / (1024 * 1024)
            
        except Exception as e:
            errors.append(f"校验过程中发生错误: {str(e)}")
        
        return {
            "valid": len(errors) == 0,
            "errors": errors,
            "warnings": warnings,
            "metadata": metadata
        }
    
    def validate_metadata_consistency(self, adata) -> List[str]:
        """检查元数据一致性"""
        warnings = []
        
        # 检查细胞ID唯一性
        if 'cell_id' in adata.obs.columns:
            if adata.obs['cell_id'].duplicated().any():
                warnings.append("细胞ID不唯一")
        
        # 检查基因ID唯一性
        if 'ensembl_gene_id' in adata.var.columns:
            if adata.var['ensembl_gene_id'].duplicated().any():
                warnings.append("Ensembl基因ID不唯一")
        
        return warnings

