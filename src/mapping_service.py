import requests
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from typing import Dict, List, Any, Optional, Tuple
import time
import json
import os
import re
from difflib import SequenceMatcher


def _column_name_suggests_cell_annotation(name: str) -> bool:
    """列名是否像细胞类型 / 注释列（用于判断是否需要用户点选）。"""
    n = str(name).lower()
    keys = (
        'cell', 'type', 'annot', 'cluster', 'leiden', 'louvain',
        'subtype', 'label', 'ontology', 'lineage', 'ident', 'seurat',
    )
    return any(k in n for k in keys)


class MappingService:
    """基因和细胞类型映射服务"""
    
    def __init__(self,
                 id_convert_fp: str = './data/raw/normalized_geneName_dict_gaogeScp.csv',
                 gtf_fp: str = './data/raw/gtfM1-M35.csv',
                 cell_type_mapping_fp: str = './cell_type_mapping/cell_type_mapping.csv'):
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
        
        # 细胞类型映射数据文件路径
        self.cell_type_mapping_fp = cell_type_mapping_fp
        
        # 加载全局映射数据（延迟加载，首次使用时加载）
        self._id_convert = None
        self._gtfM1_M35 = None
        self._gtfM32 = None
        
        # 细胞类型映射数据（延迟加载）
        self._cell_type_mapping = None  # {normalized_name: (cl_id, cl_label)}
        self._cell_type_labels = None   # 所有标准化后的名称列表（用于模糊匹配）
    
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
    
    def _load_cell_type_mapping_data(self):
        """延迟加载细胞类型映射数据"""
        if self._cell_type_mapping is not None:
            return
        
        if not os.path.exists(self.cell_type_mapping_fp):
            raise FileNotFoundError(f"细胞类型映射文件不存在: {self.cell_type_mapping_fp}")
        
        print(f"📂 加载细胞类型映射数据: {self.cell_type_mapping_fp}")
        df = pd.read_csv(self.cell_type_mapping_fp, encoding='utf-8-sig')
        
        # 构建映射字典: {normalized_name: (cl_id, cl_label)}
        # 冲突处理策略：同一 synonym 出现在多个 CL 节点时，
        # 优先保留 cl_label 与 synonym 完全一致的那条（即"正式名称"优先于"别名"）
        self._cell_type_mapping = {}
        # 记录冲突：{normalized_syn: [(cl_id, cl_label), ...]}
        conflict_map: Dict[str, List[Tuple[str, str]]] = {}
        self._cell_type_labels = []

        for _, row in df.iterrows():
            synonyms = str(row['synonyms']).split('@')
            cl_id = str(row['cl_id']).strip()
            cl_label = str(row['cl_label']).strip()

            for syn in synonyms:
                syn = syn.strip()
                if not syn or syn == 'nan':
                    continue
                normalized = self._normalize_cell_type_name(syn)
                if not normalized:
                    continue
                if normalized not in conflict_map:
                    conflict_map[normalized] = []
                conflict_map[normalized].append((cl_id, cl_label))

        duplicate_count = 0
        for normalized, candidates in conflict_map.items():
            if len(candidates) == 1:
                self._cell_type_mapping[normalized] = candidates[0]
            else:
                # 多个候选：优先选 cl_label normalize 后与 synonym 完全一致的
                exact = [c for c in candidates
                         if self._normalize_cell_type_name(c[1]) == normalized]
                if exact:
                    self._cell_type_mapping[normalized] = exact[0]
                else:
                    # 次优：选 cl_id 数值最小的（CL ontology 中编号越小越基础/通用）
                    def cl_sort_key(c):
                        try:
                            return int(c[0].replace('CL:', ''))
                        except Exception:
                            return 9999999
                    self._cell_type_mapping[normalized] = min(candidates, key=cl_sort_key)
                duplicate_count += 1
            self._cell_type_labels.append(normalized)

        print(f"✅ 加载完成，共 {len(self._cell_type_mapping)} 个映射条目"
              f"（其中 {duplicate_count} 个 synonym 存在多 CL 冲突，已自动消歧）")
    
    def _find_cell_type_column(self, adata) -> Optional[str]:
        """自动检测细胞类型列名"""
        # 按优先级排列的候选列名
        candidates = [
            'cell_type', 'celltype', 'CellType', 'cell.type',
            'raw_cell_type_label', 'cell_type_label',
            'annotation', 'cell_annotation', 'Annotation',
            'cluster', 'Cluster', 'clusters',
            'leiden', 'louvain',  # 常见聚类结果列名
            'cell_ontology_class', 'cell_ontology_id',
            'free_annotation', 'author_cell_type'
        ]
        
        for col in candidates:
            if col in adata.obs.columns:
                print(f"📍 检测到细胞类型列: {col}")
                return col
        
        # 如果没有找到精确匹配，尝试模糊匹配
        for col in adata.obs.columns:
            col_lower = col.lower()
            if 'cell' in col_lower and 'type' in col_lower:
                print(f"📍 模糊匹配到细胞类型列: {col}")
                return col
            if 'annot' in col_lower:
                print(f"📍 模糊匹配到注释列: {col}")
                return col
        
        return None

    def discover_cell_type_column_info(self, adata) -> Dict[str, Any]:
        """扫描 obs，列出可能用作细胞类型注释的列，供校验阶段与前端选择。

        requires_user_selection:
        - 自动检测失败但存在候选列；或
        - 存在且仅存在多个「列名像注释」且 n_unique 均 ≥ 3 的列（易混淆）。
        """
        recommended = self._find_cell_type_column(adata)
        candidates: List[Dict[str, Any]] = []
        max_n = min(5000, max(50, int(adata.n_obs * 0.5)))
        min_cells = max(1, int(adata.n_obs * 0.01))

        for col in adata.obs.columns:
            ser = adata.obs[col].dropna()
            if len(ser) < min_cells:
                continue
            nu = int(ser.nunique())
            if nu < 2 or nu > max_n:
                continue
            # 跳过高基数浮点连续特征
            if pd.api.types.is_float_dtype(ser.dtype) and nu > max(30, int(adata.n_obs * 0.25)):
                continue
            try:
                top = ser.astype(str).value_counts().head(3).index.tolist()
                examples = [str(x)[:48] for x in top]
            except Exception:
                examples = []
            candidates.append({
                "column": col,
                "n_unique": nu,
                "examples": examples,
            })

        candidates.sort(key=lambda x: (-x["n_unique"], x["column"]))
        candidates = candidates[:48]

        annot_like = [c for c in candidates if _column_name_suggests_cell_annotation(c["column"])]
        requires_user_selection = False
        if recommended is None and len(candidates) >= 1:
            requires_user_selection = True
        elif len(annot_like) >= 2:
            # 至少两个「看起来像」注释的列，容易选错
            nu_ok = [c for c in annot_like if c["n_unique"] >= 3]
            if len(nu_ok) >= 2:
                requires_user_selection = True

        return {
            "recommended": recommended,
            "candidates": candidates,
            "requires_user_selection": requires_user_selection,
        }
    
    def _normalize_cell_type_name(self, name: str) -> str:
        """标准化细胞类型名称"""
        if not name or pd.isna(name):
            return ""
        
        name = str(name).strip()
        
        # 转小写
        name = name.lower()
        
        # 统一分隔符
        name = name.replace('-', ' ').replace('_', ' ').replace('.', ' ')
        
        # 移除多余空格
        name = re.sub(r'\s+', ' ', name).strip()
        
        # 处理常见缩写和变体
        replacements = [
            (r'\bpositive\b', '+'),
            (r'\bnegative\b', '-'),
            (r'\bcells?\b', 'cell'),  # cells -> cell
            (r'\blymphocytes?\b', 'lymphocyte'),
            (r'\bhi\b', 'high'),
            (r'\blo\b', 'low'),
        ]
        for pattern, replacement in replacements:
            name = re.sub(pattern, replacement, name)
        
        return name.strip()
    
    def _fuzzy_match_cell_type(self, query: str, threshold: float = 0.8) -> Optional[Tuple[str, str, float]]:
        """模糊匹配细胞类型
        
        Returns:
            Tuple of (cl_id, cl_label, confidence) or None
        """
        if not self._cell_type_labels:
            return None
        
        query_normalized = self._normalize_cell_type_name(query)
        if not query_normalized:
            return None
        
        best_match = None
        best_score = 0.0
        
        for label in self._cell_type_labels:
            # 使用 SequenceMatcher 计算相似度
            score = SequenceMatcher(None, query_normalized, label).ratio()
            
            if score > best_score:
                best_score = score
                best_match = label
        
        if best_score >= threshold and best_match:
            cl_id, cl_label = self._cell_type_mapping[best_match]
            return (cl_id, cl_label, best_score)
        
        return None
    
    def _match_cell_type(self, cell_type: str) -> Dict[str, Any]:
        """匹配单个细胞类型（本地映射 + 模糊匹配）"""
        # 1. 精确匹配（标准化后）
        normalized = self._normalize_cell_type_name(cell_type)
        if normalized in self._cell_type_mapping:
            cl_id, cl_label = self._cell_type_mapping[normalized]
            return {
                "cl_id": cl_id,
                "cl_label": cl_label,
                "status": "mapped",
                "confidence": "high",
                "source": "local_exact"
            }

        # 2. 短缩写（≤4字符）不做 fuzzy：字符太短时 SequenceMatcher 极易误匹配
        #    直接返回 unmapped，避免引入错误
        if len(normalized.replace(' ', '')) <= 4:
            return {
                "cl_id": None,
                "cl_label": None,
                "status": "unmapped",
                "confidence": "low",
                "source": "local_short_abbrev"
            }

        # 3. 模糊匹配
        fuzzy_result = self._fuzzy_match_cell_type(cell_type, threshold=0.75)
        if fuzzy_result:
            cl_id, cl_label, score = fuzzy_result
            confidence = "high" if score >= 0.9 else "medium" if score >= 0.8 else "low"
            return {
                "cl_id": cl_id,
                "cl_label": cl_label,
                "status": "mapped",
                "confidence": confidence,
                "source": f"local_fuzzy (score={score:.2f})"
            }

        # 4. 未找到
        return {
            "cl_id": None,
            "cl_label": None,
            "status": "unmapped",
            "confidence": "low",
            "source": "local"
        }
    
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
    
    def map_cell_types(self, file_path: str, cell_type_column: Optional[str] = None) -> Dict[str, Any]:
        """细胞类型标准化 - 使用本地映射表（两阶段：精确匹配 + 模糊匹配）

        Args:
            file_path: h5ad 路径
            cell_type_column: 可选，用户在前端指定的 obs 列名；无效时回退自动检测。
        """
        print("🔬 开始细胞类型标准化...")

        try:
            self._load_cell_type_mapping_data()

            adata = ad.read_h5ad(file_path)

            # 自动检测细胞类型列，或尊重用户指定列
            # 注意：即使数据集已有 cl_id 列（如 Tabula Muris），仍需跑算法，
            # 因为前端需要 mapping_details（含每种 cell type 的细胞数量）来渲染图表。
            cell_type_col = None
            if cell_type_column and str(cell_type_column).strip():
                c = str(cell_type_column).strip()
                if c in adata.obs.columns:
                    cell_type_col = c
                    print(f"📍 使用用户指定细胞类型列: {cell_type_col}")
                else:
                    print(f"⚠️ 用户指定列 {c!r} 不在 obs 中，回退自动检测")
            if cell_type_col is None:
                cell_type_col = self._find_cell_type_column(adata)
            if cell_type_col is None:
                print("⚠️ 未找到细胞类型列，尝试使用的列名均不存在")
                return {
                    "mapped_count": 0,
                    "unmapped_count": int(adata.n_obs),
                    "ambiguous_count": 0,
                    "mapping_details": [],
                    "error": "未找到细胞类型标签列"
                }

            # 统计每种 cell type 对应的细胞数量（用于前端图表 Y 轴）
            cell_type_counts = adata.obs[cell_type_col].dropna().value_counts().to_dict()

            # 获取唯一的细胞类型
            cell_types = list(cell_type_counts.keys())
            total_unique_types = len(cell_types)
            print(f"📊 需要映射 {total_unique_types} 种细胞类型（共 {adata.n_obs} 个细胞）")

            mapping_details = []
            mapped_count = 0
            unmapped_count = 0
            ambiguous_count = 0
            type_to_cl = {}

            for i, cell_type in enumerate(cell_types):
                if i % 50 == 0 and i > 0:
                    print(f"🔄 处理进度: {i}/{total_unique_types}")

                cell_type_str = str(cell_type)

                if cell_type_str in self.cell_type_cache:
                    result = self.cell_type_cache[cell_type_str]
                else:
                    result = self._match_cell_type(cell_type_str)
                    self.cell_type_cache[cell_type_str] = result

                type_to_cl[cell_type_str] = result

                # cell_count：该 cell type 在数据集中的细胞数量，供前端图表使用
                mapping_details.append({
                    "raw_cell_type": cell_type_str,
                    "cl_id": result.get("cl_id"),
                    "cl_label": result.get("cl_label"),
                    "status": result.get("status"),
                    "confidence": result.get("confidence"),
                    "source": result.get("source"),
                    "cell_count": int(cell_type_counts.get(cell_type_str, 0)),
                })

                if result["status"] == "mapped":
                    mapped_count += 1
                elif result["status"] == "ambiguous":
                    ambiguous_count += 1
                else:
                    unmapped_count += 1

            # 将映射结果写入 adata.obs（覆盖已有 cl_id/cl_label，保证与算法结果一致）
            adata.obs['original_cell_type'] = adata.obs[cell_type_col].astype(str)
            adata.obs['cl_id'] = adata.obs[cell_type_col].astype(str).map(
                lambda x: type_to_cl.get(x, {}).get('cl_id')
            )
            adata.obs['cl_label'] = adata.obs[cell_type_col].astype(str).map(
                lambda x: type_to_cl.get(x, {}).get('cl_label')
            )

            adata.write_h5ad(file_path)

            success_rate = (mapped_count / total_unique_types * 100) if total_unique_types > 0 else 0.0
            
            # 计算细胞级别的映射率（mapped cells / total cells）
            mapped_cells = sum(cell_type_counts[ct] for ct in cell_types 
                             if type_to_cl.get(str(ct), {}).get('status') == 'mapped')
            total_cells = sum(cell_type_counts.values())
            cell_mapping_rate = (mapped_cells / total_cells * 100) if total_cells > 0 else 0.0

            print(f"✅ 细胞类型映射完成:")
            print(f"   - 成功: {mapped_count} ({success_rate:.1f}%)")
            print(f"   - 失败: {unmapped_count}")
            print(f"   - 模糊: {ambiguous_count}")
            print(f"   - 细胞级映射率: {mapped_cells}/{total_cells} ({cell_mapping_rate:.1f}%)")

            return {
                "mapped_count": mapped_count,
                "unmapped_count": unmapped_count,
                "ambiguous_count": ambiguous_count,
                "mapping_details": mapping_details,
                "cell_type_column": cell_type_col,
                "total_unique_types": total_unique_types,
                "success_rate": round(success_rate, 2),
                "mapped_cells": mapped_cells,
                "total_cells": total_cells,
                "cell_mapping_rate": round(cell_mapping_rate, 2)
            }

        except Exception as e:
            import traceback
            print(f"❌ 细胞类型映射失败: {str(e)}")
            traceback.print_exc()
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
