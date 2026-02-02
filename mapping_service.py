import requests
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from typing import Dict, List, Any, Optional
import time
import json
import os

class MappingService:
    """基因和细胞类型映射服务"""
    
    def __init__(self, 
                 id_convert_fp: str = './normalized_geneName_dict_gaogeScp.csv',
                 gtf_fp: str = './gtfM1-M35.csv'):
        self.mygene_base_url = "https://api.mygene.info/v3"
        self.ols_base_url = "https://www.ebi.ac.uk/ols/api"
        self.zooma_base_url = "https://www.ebi.ac.uk/spot/zooma/v2/api"
        
        # 请求超时设置
        self.timeout = 30
        
        # 缓存
        self.gene_cache = {}
        self.cell_type_cache = {}
        
        # 基因映射数据文件路径
        self.id_convert_fp = id_convert_fp
        self.gtf_fp = gtf_fp
        
        # 加载全局映射数据（延迟加载，首次使用时加载）
        self._id_convert = None
        self._gtfM1_M35 = None
        self._gtfM32 = None
    
    def _load_mapping_data(self):
        """延迟加载映射数据"""
        if self._id_convert is None:
            if not os.path.exists(self.id_convert_fp):
                raise FileNotFoundError(f"基因映射文件不存在: {self.id_convert_fp}")
            if not os.path.exists(self.gtf_fp):
                raise FileNotFoundError(f"GTF映射文件不存在: {self.gtf_fp}")
            
            print(f"📂 加载基因映射数据: {self.id_convert_fp}, {self.gtf_fp}")
            self._id_convert = pd.read_csv(self.id_convert_fp, index_col=0)
            self._gtfM1_M35 = pd.read_csv(self.gtf_fp, index_col=0)
            
            # 检查gtfM32是否为空
            self._gtfM32 = self._gtfM1_M35[self._gtfM1_M35['dataset'] == 'gencode.vM32']
            if self._gtfM32.empty:
                print("⚠️ 警告: gtfM32为空，没有dataset='gencode.vM32'的行")
    
    def map_genes(self, file_path: str) -> Dict[str, Any]:
        """基因映射 - 使用新的本地映射算法"""
        print("🧬 开始基因映射...")
        
        try:
            # 加载映射数据
            self._load_mapping_data()
            
            adata = ad.read_h5ad(file_path)
            
            # 保存原始基因数量
            total_genes = len(adata.var_names)
            print(f"📊 总基因数: {total_genes}")
            
            # 保存原始基因名到var列
            adata.var['original_gene_name'] = adata.var_names
            var_names = adata.var_names.values
            
            # 新建 id_convert_, 把 @ 项展开
            id_convert_ = self._id_convert.reset_index()
            id_convert_.columns = ['gene_name', 'alignment']
            id_convert_['alignment'] = id_convert_['alignment'].str.split('@')
            id_convert_ = id_convert_.explode('alignment').reset_index(drop=True)
            
            # 把symbol name 转化为 gene ID
            var_name_ = pd.DataFrame(var_names, columns=['gene_name'])
            var_name_ = pd.merge(var_name_, id_convert_, how='left', on='gene_name')
            var_name_['alignment'] = var_name_['alignment'].fillna(var_name_['gene_name'])
            
            # 与 gtfM32 合并
            gene_id_M32 = pd.merge(var_name_, self._gtfM32, how='left', on='alignment').fillna('UNKNOWN')
            
            # 按 gene_name_x 折叠
            gene_id_M32_ = gene_id_M32.groupby('gene_name_x').agg({
                'alignment': lambda x: list(x),
                'gene_name_y': lambda x: list(x)
            }).reset_index()
            
            # 移除gene_name_y中的UNKNOWN
            gene_id_M32_['gene_name_y'] = gene_id_M32_['gene_name_y'].apply(
                lambda x: [item for item in x if item != 'UNKNOWN']
            )
            
            # 第一步：确定初始的gene_name（处理列表类型）
            def get_initial_gene_name(row):
                """根据条件返回初始的基因名称列表/字符串"""
                if row['gene_name_x'] in row['gene_name_y']:
                    return [row['gene_name_x']]
                else:
                    return row['gene_name_y'] if row['gene_name_y'] else [row['gene_name_x']]
            
            gene_id_M32_['gene_name'] = gene_id_M32_.apply(get_initial_gene_name, axis=1)
            
            # 处理多symbol或无symbol的情况
            def process_gene_name(row):
                """统一gene_name的类型为字符串"""
                gene_name_list = row['gene_name']
                if len(gene_name_list) == 1:
                    return gene_name_list[0]
                else:
                    return row['gene_name_x']  # 用原名称替代
            
            gene_id_M32_['gene_name'] = gene_id_M32_.apply(process_gene_name, axis=1)
            
            # 重新索引并获取新的var名称
            results = gene_id_M32_.set_index('gene_name_x').reindex(var_names)['gene_name']
            new_var_names = results.fillna('UNKNOWN').values
            
            # 统计映射成功率
            gene_mapping = pd.DataFrame({
                'original': var_names,
                'mapped': new_var_names
            })
            
            # 成功映射：映射结果不是UNKNOWN
            condition1 = gene_mapping['mapped'] != 'UNKNOWN'
            success_genes = gene_mapping[condition1].shape[0]
            success_rate = (success_genes / total_genes) * 100 if total_genes > 0 else 0.0
            
            # 统计重命名数量
            renamed_count = gene_mapping[
                (gene_mapping['mapped'] != 'UNKNOWN') & 
                (gene_mapping['mapped'] != gene_mapping['original'])
            ].shape[0]
            rename_rate = (renamed_count / total_genes) * 100 if total_genes > 0 else 0.0
            
            print(f"✅ 有效映射基因数: {success_genes} ({success_rate:.2f}%)")
            print(f"✅ 重命名基因数: {renamed_count} ({rename_rate:.2f}%)")
            
            # 更新adata的var名称
            adata.var_names = new_var_names
            adata.var_names_make_unique()
            
            # 保存更新后的数据
            adata.write_h5ad(file_path)
            
            # 构建映射详情（用于返回）
            mapping_details = []
            unmapped_count = total_genes - success_genes
            ambiguous_count = 0  # 新算法中不区分ambiguous
            
            # 生成部分映射详情（前100个，避免返回数据过大）
            sample_size = min(100, len(gene_mapping))
            for idx in range(sample_size):
                row = gene_mapping.iloc[idx]
                mapping_details.append({
                    "gene_symbol": row['original'],
                    "mapped_symbol": row['mapped'],
                    "status": "mapped" if row['mapped'] != 'UNKNOWN' else "unmapped",
                    "confidence": "high" if row['mapped'] != 'UNKNOWN' else "low",
                    "source": "local_mapping"
                })
            
            print(f"✅ 基因映射完成: 成功 {success_genes}, 失败 {unmapped_count}")
            
            return {
                "mapped_count": success_genes,
                "unmapped_count": unmapped_count,
                "ambiguous_count": ambiguous_count,
                "mapping_details": mapping_details,
                "total_genes": total_genes,
                "success_rate": round(success_rate, 2),
                "rename_rate": round(rename_rate, 2)
            }
            
        except Exception as e:
            import traceback
            print(f"❌ 基因映射失败: {str(e)}")
            traceback.print_exc()
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
            
            # 最小改动：若已有 CL 列且存在非空值，则跳过外部映射
            if 'cl_id' in adata.obs.columns:
                existing = int(adata.obs['cl_id'].notna().sum())
                total = int(adata.n_obs)
                if existing > 0:
                    print(f"🔬 跳过外部细胞类型映射（已有 CL 覆盖 {existing}/{total}）")
                    return {
                        "mapped_count": existing,
                        "unmapped_count": total - existing,
                        "ambiguous_count": 0,
                        "mapping_details": [],
                        "note": "skipped_external_mapping_existing_cl"
                    }
            
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
