"""
h5ad 文件检查脚本 —— 快速判断数据集是否已标注、包含哪些元数据
用法: python inspect_h5ad.py your_file.h5ad
"""
import sys
import anndata as ad
import pandas as pd


def inspect(file_path: str):
    print(f"\n{'='*60}")
    print(f"  文件: {file_path}")
    print(f"{'='*60}")

    adata = ad.read_h5ad(file_path)

    # ---- 基本信息 ----
    print(f"\n[基本信息]")
    print(f"  细胞数 (n_obs):  {adata.n_obs}")
    print(f"  基因数 (n_vars): {adata.n_vars}")
    print(f"  稀疏矩阵: {hasattr(adata.X, 'nnz')}")

    # ---- obs 列（细胞级元数据）—— 判断标注的关键 ----
    print(f"\n[obs 列] (共 {len(adata.obs.columns)} 列)")
    for col in adata.obs.columns:
        dtype = adata.obs[col].dtype
        nunique = adata.obs[col].nunique()
        sample = adata.obs[col].dropna().unique()[:5].tolist()
        print(f"  - {col:30s}  dtype={str(dtype):12s}  unique={nunique:6d}  示例: {sample}")

    # ---- 自动检测可能的细胞类型列 ----
    annotation_keywords = [
        'cell_type', 'celltype', 'cell.type', 'annotation', 'cell_annotation',
        'cell_ontology', 'cl_id', 'cl_label', 'cluster', 'leiden', 'louvain',
        'cell_identity', 'author_cell_type', 'free_annotation', 'label',
    ]
    print(f"\n[标注检测] 自动扫描疑似细胞类型列:")
    found_any = False
    for col in adata.obs.columns:
        col_lower = col.lower()
        for kw in annotation_keywords:
            if kw in col_lower:
                values = adata.obs[col].dropna().unique()
                is_numeric = pd.api.types.is_numeric_dtype(adata.obs[col])
                label = "⚠️ 数字型(可能是聚类编号)" if is_numeric else "✅ 文本型(可能是标注)"
                print(f"  ➜ {col:30s}  {label}  ({len(values)} 种)")
                if not is_numeric and len(values) <= 30:
                    print(f"    值: {sorted([str(v) for v in values])}")
                elif not is_numeric:
                    print(f"    前10个: {sorted([str(v) for v in values[:10]])}")
                found_any = True
                break
    if not found_any:
        print("  ❌ 未检测到疑似细胞类型列 —— 该数据集可能未标注")

    # ---- var 列（基因级元数据）----
    print(f"\n[var 列] (共 {len(adata.var.columns)} 列)")
    for col in adata.var.columns:
        dtype = adata.var[col].dtype
        print(f"  - {col:30s}  dtype={str(dtype)}")
    print(f"  基因名前5个: {adata.var_names[:5].tolist()}")

    # ---- uns（全局元数据）----
    if adata.uns:
        print(f"\n[uns 全局元数据] (共 {len(adata.uns)} 个键)")
        for key in list(adata.uns.keys())[:15]:
            val = adata.uns[key]
            val_type = type(val).__name__
            print(f"  - {key:30s}  type={val_type}")
    else:
        print(f"\n[uns] 无全局元数据")

    # ---- obsm（降维结果）----
    if adata.obsm:
        print(f"\n[obsm 降维结果]")
        for key in adata.obsm.keys():
            shape = adata.obsm[key].shape
            print(f"  - {key:30s}  shape={shape}")

    # ---- 总结 ----
    print(f"\n{'='*60}")
    print("  总结:")
    if found_any:
        print("  ✅ 该文件可能包含细胞类型标注信息，请查看上方 [标注检测] 部分")
    else:
        print("  ❌ 未检测到明确的细胞类型标注，可能需要人工确认 obs 列含义")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python inspect_h5ad.py <h5ad文件路径>")
        print("示例: python inspect_h5ad.py data/my_dataset.h5ad")
        sys.exit(1)
    inspect(sys.argv[1])
