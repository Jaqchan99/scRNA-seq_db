import scanpy as sc
import pandas as pd
import numpy as np
import json
import os
from datetime import datetime
from typing import Dict, List, Any

class ExportService:
    """数据导出服务"""
    
    def __init__(self):
        self.export_dir = "exports"
        self.report_dir = "reports"
        
        # 确保目录存在
        os.makedirs(self.export_dir, exist_ok=True)
        os.makedirs(self.report_dir, exist_ok=True)
    
    def export_processed_data(
        self, 
        submission_id: str, 
        file_path: str,
        gene_mapping_result: Dict[str, Any],
        cell_type_mapping_result: Dict[str, Any]
    ) -> Dict[str, str]:
        """导出处理后的数据"""
        print(f"📤 开始导出数据 {submission_id}")
        
        try:
            # 读取原始数据
            adata = sc.read_h5ad(file_path)
            
            # 应用基因映射结果
            adata = self._apply_gene_mapping(adata, gene_mapping_result)
            
            # 应用细胞类型映射结果
            adata = self._apply_cell_type_mapping(adata, cell_type_mapping_result)
            
            # 添加处理元数据
            adata.uns['processing_info'] = {
                'submission_id': submission_id,
                'processed_at': datetime.now().isoformat(),
                # 仅保留可序列化的统计摘要，避免将完整详情放入 uns
                'gene_mapping_mapped': int(gene_mapping_result.get('mapped_count', 0)),
                'gene_mapping_unmapped': int(gene_mapping_result.get('unmapped_count', 0)),
                'gene_mapping_ambiguous': int(gene_mapping_result.get('ambiguous_count', 0)),
                'cell_type_mapped': int(cell_type_mapping_result.get('mapped_count', 0)),
                'cell_type_unmapped': int(cell_type_mapping_result.get('unmapped_count', 0)),
                'cell_type_ambiguous': int(cell_type_mapping_result.get('ambiguous_count', 0)),
                'platform_version': '1.0.0'
            }
            
            # 导出处理后的 h5ad 文件
            export_path = os.path.join(self.export_dir, f"processed_{submission_id}.h5ad")
            adata.write_h5ad(export_path)
            
            # 生成处理报告
            report_path = self._generate_report(
                submission_id, 
                adata, 
                gene_mapping_result, 
                cell_type_mapping_result
            )
            
            print(f"✅ 数据导出完成: {export_path}")
            
            return {
                "export_path": export_path,
                "report_path": report_path
            }
            
        except Exception as e:
            print(f"❌ 数据导出失败: {str(e)}")
            raise e
    
    def _apply_gene_mapping(self, adata, gene_mapping_result: Dict[str, Any]) -> sc.AnnData:
        """应用基因映射结果"""
        print("🧬 应用基因映射结果...")
        
        # 新算法已经在 mapping_service 中直接更新了 adata.var_names
        # 这里只需记录统计信息
        
        mapped_count = gene_mapping_result.get("mapped_count", 0)
        total_genes = gene_mapping_result.get("total_genes", len(adata.var))
        success_rate = gene_mapping_result.get("success_rate", 0)
        
        # 检查是否有 original_gene_name 列（新算法会添加）
        if 'original_gene_name' in adata.var.columns:
            print(f"📊 基因映射统计: {mapped_count}/{total_genes} ({success_rate}%)")
        else:
            # 兼容旧数据：如果没有 original_gene_name，尝试使用 mapping_details
            mapping_details = gene_mapping_result.get("mapping_details", [])
            if mapping_details:
                # 兼容新旧两种格式
                symbol_to_mapped = {}
                for detail in mapping_details:
                    symbol = detail.get("gene_symbol", "")
                    # 新格式用 mapped_symbol，旧格式用 ensembl_id
                    mapped = detail.get("mapped_symbol") or detail.get("ensembl_id")
                    if mapped and detail.get("status") == "mapped":
                        symbol_to_mapped[symbol] = mapped
                
                print(f"📊 基因映射详情: {len(symbol_to_mapped)} 个基因有映射记录")
        
        return adata
    
    def _apply_cell_type_mapping(self, adata, cell_type_mapping_result: Dict[str, Any]) -> sc.AnnData:
        """应用细胞类型映射结果"""
        print("🔬 应用细胞类型映射结果...")
        
        # 创建细胞类型映射详情
        mapping_details = cell_type_mapping_result.get("mapping_details", [])
        
        # 创建映射字典
        raw_to_cl = {}
        raw_to_cl_label = {}
        for detail in mapping_details:
            raw_type = detail["raw_cell_type"]
            cl_id = detail["cl_id"]
            cl_label = detail["cl_label"]
            if cl_id:
                raw_to_cl[raw_type] = cl_id
                raw_to_cl_label[raw_type] = cl_label
        
        # 添加 CL ID 列到 obs
        if 'raw_cell_type_label' in adata.obs.columns:
            # 不覆盖已有 cl_id/cl_label
            if 'cl_id' not in adata.obs.columns:
                adata.obs['cl_id'] = None
            if 'cl_label' not in adata.obs.columns:
                adata.obs['cl_label'] = None
            
            new_cl = adata.obs['raw_cell_type_label'].map(raw_to_cl)
            new_cl_label = adata.obs['raw_cell_type_label'].map(raw_to_cl_label)
            
            # 仅在原值为空时写入
            adata.obs.loc[adata.obs['cl_id'].isna(), 'cl_id'] = new_cl
            adata.obs.loc[adata.obs['cl_label'].isna(), 'cl_label'] = new_cl_label
            
            # 统计映射情况
            mapped_cells = adata.obs['cl_id'].notna().sum()
            total_cells = len(adata.obs)
            
            print(f"📊 细胞类型映射统计: {mapped_cells}/{total_cells} ({mapped_cells/total_cells*100:.1f}%)")
        
        return adata
    
    def _generate_report(
        self, 
        submission_id: str, 
        adata: sc.AnnData,
        gene_mapping_result: Dict[str, Any],
        cell_type_mapping_result: Dict[str, Any]
    ) -> str:
        """生成处理报告"""
        print("📋 生成处理报告...")
        
        report = {
            "submission_id": submission_id,
            "generated_at": datetime.now().isoformat(),
            "summary": {
                "n_cells": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
                "file_size_mb": round(os.path.getsize(f"uploads/{submission_id}.h5ad") / (1024 * 1024), 2)
            },
            "validation": {
                "status": "passed",
                "warnings": [],
                "metadata": {}
            },
            "gene_mapping": {
                "total_genes": gene_mapping_result.get("total_genes")
                    or len(gene_mapping_result.get("mapping_details", [])),
                "mapped_count": gene_mapping_result.get("mapped_count", 0),
                "unmapped_count": gene_mapping_result.get("unmapped_count", 0),
                "ambiguous_count": gene_mapping_result.get("ambiguous_count", 0),
                # 与评估侧一致：成功映射基因数 / 总基因数（而非 / 详情条数）
                "mapping_rate": round(
                    gene_mapping_result.get("mapped_count", 0)
                    / max(
                        gene_mapping_result.get("total_genes")
                        or len(gene_mapping_result.get("mapping_details", [])),
                        1,
                    )
                    * 100,
                    2,
                ),
            },
            "cell_type_mapping": {
                "total_types": len(cell_type_mapping_result.get("mapping_details", [])),
                "mapped_count": cell_type_mapping_result.get("mapped_count", 0),
                "unmapped_count": cell_type_mapping_result.get("unmapped_count", 0),
                "ambiguous_count": cell_type_mapping_result.get("ambiguous_count", 0),
                # 与 evaluate_mapping 一致：映射到 CL 的细胞数 / 总细胞数
                "mapping_rate": (
                    round(
                        cell_type_mapping_result.get("mapped_cells", 0)
                        / max(cell_type_mapping_result.get("total_cells", 0), 1)
                        * 100,
                        2,
                    )
                    if cell_type_mapping_result.get("total_cells", 0) > 0
                    else round(cell_type_mapping_result.get("cell_mapping_rate", 0), 2)
                ),
                "mapped_cells": cell_type_mapping_result.get("mapped_cells", 0),
                "total_cells": cell_type_mapping_result.get("total_cells", 0),
                "cell_mapping_rate": cell_type_mapping_result.get("cell_mapping_rate", 0),
            },
            "quality_metrics": self._calculate_quality_metrics(adata),
            "mapping_details": {
                "gene_mapping": gene_mapping_result.get("mapping_details", []),
                "cell_type_mapping": cell_type_mapping_result.get("mapping_details", [])
            },
            "recommendations": self._generate_recommendations(adata, gene_mapping_result, cell_type_mapping_result)
        }
        
        # 保存 JSON 报告
        report_path = os.path.join(self.report_dir, f"report_{submission_id}.json")
        with open(report_path, 'w', encoding='utf-8') as f:
            json.dump(report, f, indent=2, ensure_ascii=False)
        
        # 生成 TSV 格式的映射详情
        self._generate_tsv_report(submission_id, gene_mapping_result, cell_type_mapping_result)
        
        print(f"📋 处理报告已生成: {report_path}")
        
        return report_path
    
    def _calculate_quality_metrics(self, adata: sc.AnnData) -> Dict[str, Any]:
        """计算质量指标"""
        metrics = {}
        
        # 基本统计
        metrics["n_cells"] = int(adata.n_obs)
        metrics["n_genes"] = int(adata.n_vars)
        
        # 表达矩阵统计
        if hasattr(adata.X, 'data'):
            metrics["sparsity"] = round(1 - (adata.X.nnz / (adata.n_obs * adata.n_vars)), 4)
            metrics["mean_expression"] = round(float(adata.X.data.mean()), 4)
            metrics["median_expression"] = round(float(np.median(adata.X.data)), 4)
        
        # 细胞质量指标
        if 'n_genes' in adata.obs.columns:
            metrics["mean_genes_per_cell"] = round(float(adata.obs['n_genes'].mean()), 2)
            metrics["median_genes_per_cell"] = round(float(adata.obs['n_genes'].median()), 2)
        
        if 'n_counts' in adata.obs.columns:
            metrics["mean_counts_per_cell"] = round(float(adata.obs['n_counts'].mean()), 2)
            metrics["median_counts_per_cell"] = round(float(adata.obs['n_counts'].median()), 2)
        
        if 'pct_mt' in adata.obs.columns:
            metrics["mean_mt_percentage"] = round(float(adata.obs['pct_mt'].mean()), 2)
            metrics["median_mt_percentage"] = round(float(adata.obs['pct_mt'].median()), 2)
        
        return metrics
    
    def _generate_recommendations(
        self, 
        adata: sc.AnnData, 
        gene_mapping_result: Dict[str, Any], 
        cell_type_mapping_result: Dict[str, Any]
    ) -> List[str]:
        """生成处理建议"""
        recommendations = []
        
        # 基因映射建议（mapped / total_genes）
        tg = gene_mapping_result.get("total_genes") or len(gene_mapping_result.get("mapping_details", []))
        gene_mapping_rate = gene_mapping_result.get("mapped_count", 0) / max(tg, 1)
        if gene_mapping_rate < 0.8:
            recommendations.append(f"基因映射率较低 ({gene_mapping_rate*100:.1f}%)，建议检查基因符号格式")

        # 细胞类型：细胞级映射率
        tc = cell_type_mapping_result.get("total_cells", 0)
        if tc > 0:
            cell_type_mapping_rate = cell_type_mapping_result.get("mapped_cells", 0) / tc
        else:
            cell_type_mapping_rate = cell_type_mapping_result.get("mapped_count", 0) / max(
                len(cell_type_mapping_result.get("mapping_details", [])), 1
            )
        if cell_type_mapping_rate < 0.5:
            recommendations.append(f"细胞类型映射率较低 ({cell_type_mapping_rate*100:.1f}%)，建议人工审核细胞类型标签")
        
        # 质量指标建议
        if 'pct_mt' in adata.obs.columns:
            high_mt_cells = (adata.obs['pct_mt'] > 20).sum()
            if high_mt_cells > adata.n_obs * 0.1:
                recommendations.append(f"发现 {high_mt_cells} 个高线粒体基因比例的细胞 (>20%)，建议检查数据质量")
        
        if 'n_genes' in adata.obs.columns:
            low_gene_cells = (adata.obs['n_genes'] < 200).sum()
            if low_gene_cells > adata.n_obs * 0.1:
                recommendations.append(f"发现 {low_gene_cells} 个低基因数的细胞 (<200)，建议检查数据质量")
        
        return recommendations
    
    def _generate_tsv_report(
        self, 
        submission_id: str, 
        gene_mapping_result: Dict[str, Any], 
        cell_type_mapping_result: Dict[str, Any]
    ):
        """生成 TSV 格式的映射详情报告"""
        # 基因映射详情
        gene_details = gene_mapping_result.get("mapping_details", [])
        if gene_details:
            gene_df = pd.DataFrame(gene_details)
            gene_tsv_path = os.path.join(self.report_dir, f"gene_mapping_{submission_id}.tsv")
            gene_df.to_csv(gene_tsv_path, sep='\t', index=False)
        
        # 细胞类型映射详情
        cell_type_details = cell_type_mapping_result.get("mapping_details", [])
        if cell_type_details:
            cell_type_df = pd.DataFrame(cell_type_details)
            cell_type_tsv_path = os.path.join(self.report_dir, f"cell_type_mapping_{submission_id}.tsv")
            cell_type_df.to_csv(cell_type_tsv_path, sep='\t', index=False)

