# Benchmark 标签提取说明

## 1. 你要提取的到底是什么？

- **不要**用 `cell` 列：那是**细胞 barcode**（如 `10X_P7_11_AAACCTGAGACAGGCT`），用来区分每一个细胞，不是细胞类型。
- **要**提取的是：
  - **文本标签**：人类给的细胞类型名称，例如 `immature T cell`、`CD4+ T cell`，也就是你将来要喂给「cell type mapper」的输入；
  - **CL ID / CL label**（若有）：即标准答案，例如 `CL:0002420`、`immature T cell`，用来评估 mapper 是否映射正确。

所以：**做 benchmark 时，就是把「文本标签」列和「CL ID」列提取出来；若某文件只有文本没有 CL，就只提取文本，后续人工填 CL ID。**

## 2. 你读到的两列分别是什么？

你 inspect 看到的：

- `cell_ontology_class`：**CL 标准名称**（如 `immature T cell`），可以同时当作「mapper 的输入」和「标准名称」；
- `cell_ontology_id`：**CL ID**（如 `CL:0002420`），标准答案。

对这种**已经有 CL 的数据集**，benchmark 用法就是：

- **输入**：`cell_ontology_class`（或该文件里别的自由文本 cell type 列）；
- **标准答案**：`cell_ontology_id`（+ 可选 `cell_ontology_class` 作为标准名称）。

脚本会从每个 h5ad 的 `obs` 里自动识别这些列并抽出唯一 (标签, CL ID) 对。

## 3. 对所有 h5ad 怎么跑？

- 把 90 个 h5ad 的路径放到一个目录里，或写进一个文本文件（每行一个路径）。
- 在项目根目录执行：

```bash
# 方式一：指定目录，会递归找所有 .h5ad
python extract_benchmark_labels.py /path/to/h5ad_folder --output benchmark_candidates.csv

# 方式二：指定一个 txt，里面每行一个 h5ad 路径
python extract_benchmark_labels.py h5ad_list.txt --output benchmark_candidates.csv

# 只跑前 10 个（试跑）
python extract_benchmark_labels.py /path/to/h5ad_folder --output bench_sample.csv --limit 10
```

- 输出 CSV 列：`source_file`, `label_text`, `cl_id`, `cl_label`, `cell_count`, `has_ground_truth`。
  - 若某行 `cl_id` 非空，说明该条**已有 ground truth**，可直接用于评估；
  - 若 `cl_id` 为空，说明只抽到了文本标签，需要**人工查 Cell Ontology 填 CL ID** 后再做 benchmark。

## 4. 小结

1. **做 benchmark 时提取的是**：细胞类型**文本标签**（和可选的 CL 标准名） + **CL ID**（标准答案）；`cell` 列是 barcode，不参与。
2. **已有 `cell_ontology_class` + `cell_ontology_id` 的文件**：脚本会直接抽出 (标签 → CL ID)，当作现成 ground truth。
3. **对所有 h5ad 批量跑**：用 `extract_benchmark_labels.py` 指定目录或路径列表，得到一张表，再对缺 CL 的行做人工标注即可。
