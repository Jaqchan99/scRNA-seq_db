import scanpy as sc
adata = sc.read_h5ad("Zilionis_mouse_lung.h5ad")
# print(adata)
print(f"细胞数: {adata.n_obs}, 基因数: {adata.n_vars}")
print("基因名称:",adata.var.index)  # 打印所有基因名称
print("细胞类型（CL）:",adata.obs["Most_likely_Immgen_cell_type"],adata.obs["cell_ontology_id"])

# print("细胞类型:", adata.obs['Most_likely_Immgen_cell_type'].unique())
# print("基因符号:", adata.var['gene_symbol'].head())