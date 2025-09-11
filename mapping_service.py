import requests
import pandas as pd
import numpy as np
import scanpy as sc
from typing import Dict, List, Any, Optional
import time
import json

class MappingService:
    """基因和细胞类型映射服务"""
    
    def __init__(self):
        self.mygene_base_url = "https://api.mygene.info/v3"
        self.ols_base_url = "https://www.ebi.ac.uk/ols/api"
        self.zooma_base_url = "https://www.ebi.ac.uk/spot/zooma/v2/api"
        
        # 请求超时设置
        self.timeout = 30
        
        # 缓存
        self.gene_cache = {}
        self.cell_type_cache = {}
    
    def map_genes(self, file_path: str) -> Dict[str, Any]:
        """基因映射"""
        print("🧬 开始基因映射...")
        
        try:
            adata = sc.read_h5ad(file_path)
            mapping_details = []
            mapped_count = 0
            unmapped_count = 0
            ambiguous_count = 0
            
            # 检查是否有基因符号列
            if 'gene_symbol' not in adata.var.columns:
                return {
                    "mapped_count": 0,
                    "unmapped_count": adata.n_vars,
                    "ambiguous_count": 0,
                    "mapping_details": [],
                    "error": "未找到基因符号列"
                }
            
            gene_symbols = adata.var['gene_symbol'].dropna().unique()
            print(f"📊 需要映射 {len(gene_symbols)} 个基因符号")
            
            for i, symbol in enumerate(gene_symbols):
                if i % 100 == 0:
                    print(f"🔄 处理进度: {i}/{len(gene_symbols)}")
                
                # 检查缓存
                if symbol in self.gene_cache:
                    result = self.gene_cache[symbol]
                else:
                    result = self._query_gene_mapping(symbol)
                    self.gene_cache[symbol] = result
                
                mapping_details.append({
                    "gene_symbol": symbol,
                    "ensembl_id": result.get("ensembl_id"),
                    "status": result.get("status"),
                    "confidence": result.get("confidence"),
                    "source": result.get("source")
                })
                
                if result["status"] == "mapped":
                    mapped_count += 1
                elif result["status"] == "ambiguous":
                    ambiguous_count += 1
                else:
                    unmapped_count += 1
                
                # 避免请求过于频繁
                time.sleep(0.1)
            
            print(f"✅ 基因映射完成: 成功 {mapped_count}, 模糊 {ambiguous_count}, 失败 {unmapped_count}")
            
            return {
                "mapped_count": mapped_count,
                "unmapped_count": unmapped_count,
                "ambiguous_count": ambiguous_count,
                "mapping_details": mapping_details
            }
            
        except Exception as e:
            print(f"❌ 基因映射失败: {str(e)}")
            return {
                "mapped_count": 0,
                "unmapped_count": 0,
                "ambiguous_count": 0,
                "mapping_details": [],
                "error": str(e)
            }
    
    def _query_gene_mapping(self, gene_symbol: str) -> Dict[str, Any]:
        """查询单个基因的映射"""
        try:
            # 查询 mygene.info
            url = f"{self.mygene_base_url}/query"
            params = {
                "q": gene_symbol,
                "species": "mouse",
                "fields": "ensembl.gene,symbol,name",
                "size": 1
            }
            
            response = requests.get(url, params=params, timeout=self.timeout)
            
            if response.status_code == 200:
                data = response.json()
                
                if data.get("hits") and len(data["hits"]) > 0:
                    hit = data["hits"][0]
                    ensembl_id = hit.get("ensembl", {}).get("gene")
                    
                    if ensembl_id:
                        return {
                            "ensembl_id": ensembl_id,
                            "status": "mapped",
                            "confidence": "high",
                            "source": "mygene.info"
                        }
                    else:
                        return {
                            "ensembl_id": None,
                            "status": "unmapped",
                            "confidence": "low",
                            "source": "mygene.info"
                        }
                else:
                    return {
                        "ensembl_id": None,
                        "status": "unmapped",
                        "confidence": "low",
                        "source": "mygene.info"
                    }
            else:
                return {
                    "ensembl_id": None,
                    "status": "error",
                    "confidence": "low",
                    "source": "mygene.info",
                    "error": f"HTTP {response.status_code}"
                }
                
        except requests.exceptions.Timeout:
            return {
                "ensembl_id": None,
                "status": "error",
                "confidence": "low",
                "source": "mygene.info",
                "error": "请求超时"
            }
        except Exception as e:
            return {
                "ensembl_id": None,
                "status": "error",
                "confidence": "low",
                "source": "mygene.info",
                "error": str(e)
            }
    
    def map_cell_types(self, file_path: str) -> Dict[str, Any]:
        """细胞类型标准化"""
        print("🔬 开始细胞类型标准化...")
        
        try:
            adata = sc.read_h5ad(file_path)
            mapping_details = []
            mapped_count = 0
            unmapped_count = 0
            ambiguous_count = 0
            
            # 检查是否有细胞类型列
            if 'raw_cell_type_label' not in adata.obs.columns:
                return {
                    "mapped_count": 0,
                    "unmapped_count": adata.n_obs,
                    "ambiguous_count": 0,
                    "mapping_details": [],
                    "error": "未找到细胞类型标签列"
                }
            
            cell_types = adata.obs['raw_cell_type_label'].dropna().unique()
            print(f"📊 需要映射 {len(cell_types)} 种细胞类型")
            
            for i, cell_type in enumerate(cell_types):
                if i % 50 == 0:
                    print(f"🔄 处理进度: {i}/{len(cell_types)}")
                
                # 检查缓存
                if cell_type in self.cell_type_cache:
                    result = self.cell_type_cache[cell_type]
                else:
                    result = self._query_cell_type_mapping(cell_type)
                    self.cell_type_cache[cell_type] = result
                
                mapping_details.append({
                    "raw_cell_type": cell_type,
                    "cl_id": result.get("cl_id"),
                    "cl_label": result.get("cl_label"),
                    "status": result.get("status"),
                    "confidence": result.get("confidence"),
                    "source": result.get("source")
                })
                
                if result["status"] == "mapped":
                    mapped_count += 1
                elif result["status"] == "ambiguous":
                    ambiguous_count += 1
                else:
                    unmapped_count += 1
                
                # 避免请求过于频繁
                time.sleep(0.2)
            
            print(f"✅ 细胞类型映射完成: 成功 {mapped_count}, 模糊 {ambiguous_count}, 失败 {unmapped_count}")
            
            return {
                "mapped_count": mapped_count,
                "unmapped_count": unmapped_count,
                "ambiguous_count": ambiguous_count,
                "mapping_details": mapping_details
            }
            
        except Exception as e:
            print(f"❌ 细胞类型映射失败: {str(e)}")
            return {
                "mapped_count": 0,
                "unmapped_count": 0,
                "ambiguous_count": 0,
                "mapping_details": [],
                "error": str(e)
            }
    
    def _query_cell_type_mapping(self, cell_type: str) -> Dict[str, Any]:
        """查询单个细胞类型的映射"""
        try:
            # 使用 Zooma 进行细胞类型映射
            url = f"{self.zooma_base_url}/services/annotate"
            params = {
                "propertyValue": cell_type,
                "filter": "required:[none],preferred:[cl]",
                "confidence": "medium"
            }
            
            response = requests.get(url, params=params, timeout=self.timeout)
            
            if response.status_code == 200:
                data = response.json()
                
                if data and len(data) > 0:
                    # 取第一个结果
                    result = data[0]
                    annotations = result.get("annotatedProperty", {}).get("annotations", [])
                    
                    if annotations:
                        annotation = annotations[0]
                        cl_id = annotation.get("semanticTag")
                        cl_label = annotation.get("propertyValue")
                        confidence = annotation.get("confidence", "medium")
                        
                        if cl_id and "CL:" in cl_id:
                            return {
                                "cl_id": cl_id,
                                "cl_label": cl_label,
                                "status": "mapped",
                                "confidence": confidence,
                                "source": "zooma"
                            }
                        else:
                            return {
                                "cl_id": None,
                                "cl_label": None,
                                "status": "unmapped",
                                "confidence": "low",
                                "source": "zooma"
                            }
                    else:
                        return {
                            "cl_id": None,
                            "cl_label": None,
                            "status": "unmapped",
                            "confidence": "low",
                            "source": "zooma"
                        }
                else:
                    return {
                        "cl_id": None,
                        "cl_label": None,
                        "status": "unmapped",
                        "confidence": "low",
                        "source": "zooma"
                    }
            else:
                return {
                    "cl_id": None,
                    "cl_label": None,
                    "status": "error",
                    "confidence": "low",
                    "source": "zooma",
                    "error": f"HTTP {response.status_code}"
                }
                
        except requests.exceptions.Timeout:
            return {
                "cl_id": None,
                "cl_label": None,
                "status": "error",
                "confidence": "low",
                "source": "zooma",
                "error": "请求超时"
            }
        except Exception as e:
            return {
                "cl_id": None,
                "cl_label": None,
                "status": "error",
                "confidence": "low",
                "source": "zooma",
                "error": str(e)
            }
    
    def get_mapping_statistics(self, mapping_result: Dict[str, Any]) -> Dict[str, Any]:
        """获取映射统计信息"""
        total = mapping_result["mapped_count"] + mapping_result["unmapped_count"] + mapping_result["ambiguous_count"]
        
        if total == 0:
            return {
                "total": 0,
                "mapped_percentage": 0,
                "unmapped_percentage": 0,
                "ambiguous_percentage": 0
            }
        
        return {
            "total": total,
            "mapped_percentage": round(mapping_result["mapped_count"] / total * 100, 2),
            "unmapped_percentage": round(mapping_result["unmapped_count"] / total * 100, 2),
            "ambiguous_percentage": round(mapping_result["ambiguous_count"] / total * 100, 2)
        }
