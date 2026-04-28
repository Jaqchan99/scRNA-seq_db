# scRNA-seq Data Preprocessing and Standardisation Platform

A modular platform for standardising single-cell RNA sequencing (scRNA-seq) datasets prior to downstream quality control and analysis. The system performs automated structural validation, gene identifier normalisation, and cell type ontology alignment — all entirely offline using curated local reference tables.

## Key Results

| Metric | Local (ours) | ZOOMA |
|---|---|---|
| **Coverage (Label)** | **100.0%** | 33.64% |
| **Accuracy (Label)** | **97.40%** | 25.28% |
| **Precision (Label)** | 97.40% | 75.14% |
| Accuracy (Cell-weighted) | **98.36%** | 20.62% |
| Benchmark | 538 labels, 56 datasets, 1,072,736 cells | same |

> External validation on held-out GSE166504 (12 cell types): 83.33% accuracy (ours) vs. 41.7% (ZOOMA)

## Features

- **Gene Name Alignment**: Offline normalisation to GENCODE M32 using local conversion tables (95%+ mapping rate)
- **Cell Type Standardisation**: Two-stage local algorithm (exact match + fuzzy match) mapping free-text labels to Cell Ontology (CL) terms
- **Structural Validation**: Schema checking, integrity verification, biological plausibility checks
- **Incremental Update SOP**: Semi-automated pipeline to expand the lookup table with new datasets, with before/after evaluation
- **Reproducible Reporting**: JSON + TSV audit trails for all transformations

## Tech Stack

- **Backend**: FastAPI + SQLite
- **Data Processing**: scanpy + pandas + numpy + anndata
- **Gene Mapping**: Local reference tables (id_convert + GTF M32) — no network dependency
- **Cell Type Mapping**: Local Cell Ontology synonym table (3,240 CL entries, 6,120+ synonyms) — no network dependency
- **Frontend**: Vanilla HTML/CSS/JS (lightweight, no framework dependency)
- **File Format**: h5ad (AnnData)

## Project Structure

```
data_platform/
├── src/                          # Core backend modules
│   ├── main.py                   # FastAPI application entry point
│   ├── database.py               # SQLite database configuration
│   ├── models.py                 # Pydantic data models
│   ├── upload_service.py         # File ingestion and integrity checks
│   ├── validation_service.py     # Structural and semantic validation
│   ├── mapping_service.py        # Gene + cell type mapping algorithms
│   ├── export_service.py         # Standardised output generation
│   ├── read_h5ad.py              # h5ad file inspection utilities
│   └── check_database.py         # Database inspection utilities
├── scripts/                      # Evaluation and utility scripts
│   ├── incremental_update_sop.py # Incremental update orchestration
│   ├── evaluate_mapping.py       # Benchmark evaluation (local + ZOOMA)
│   ├── extract_benchmark_labels.py
│   ├── batch_inspect.py          # Dataset inspection and config generation
│   ├── merge_lookup_table.py     # Lookup table merge utility
│   ├── fix_lookup_errors.py      # Error correction utility
│   └── ...
├── cell_type_mapping/            # Cell Ontology reference data
│   ├── cell_type_mapping.csv     # Parsed CL synonym table (3,240 entries)
│   └── parse_cl_obo.py           # OBO → CSV parser
├── data/                         # (local, gitignored) raw references, h5ad, benchmarks
├── configs/                      # (local, gitignored) batch_inspect / SOP configs
├── results/                      # (local, gitignored) eval outputs, logs, coverage
├── frontend/                     # Web frontend
│   ├── index.html
│   ├── style.css
│   └── script.js
├── 使用说明.md                   # Chinese usage guide
├── requirements.txt
└── README.md
```

## Quick Start

### 1. Install Dependencies

```bash
pip install -r requirements.txt
```

### 2. Launch Platform (Recommended)

```bash
python scripts/start_platform.py
```

This starts:
- Backend API: `http://localhost:8100`
- Frontend: `http://localhost:8080` (default; if the port is busy, the launcher may pick another — check the terminal output)

### 3. Manual Start

```bash
# Start backend
python src/main.py

# In another terminal, start frontend
cd frontend && python -m http.server 6666
```

### 4. Usage

**Web Interface**: Open `http://localhost:8080` (or the URL printed by `scripts/start_platform.py`) and follow the 4-step guided workflow:
1. Upload .h5ad file
2. Review validation results → confirm
3. Review mapping results → confirm
4. Download standardised output

**API**:
```bash
# Create submission
curl -X POST "http://localhost:8100/submissions"

# Upload file
curl -X POST "http://localhost:8100/submissions/{id}/upload" -F "file=@data.h5ad"

# Check status
curl "http://localhost:8100/submissions/{id}/status"

# Download report / processed file
curl "http://localhost:8100/submissions/{id}/report" -o report.json
curl "http://localhost:8100/submissions/{id}/export" -o processed.h5ad
```

## Incremental Update SOP

To expand the cell type lookup table with new datasets (see also **使用说明.md** in Chinese):

```bash
python scripts/incremental_update_sop.py --h5ad_dir /path/to/h5ad_folder
# After filling selected_column in the generated configs/dataset_config_<run_id>.csv:
python scripts/incremental_update_sop.py --resume results/logs/sop_state_<run_id>.json
# Or jump to a phase, e.g. merge:
python scripts/incremental_update_sop.py --h5ad_dir /path/to/h5ad_folder --start_phase merge
```

## Evaluation

Run benchmark evaluation:

```bash
# Local method only
python scripts/evaluate_mapping.py --method local

# Compare local vs. ZOOMA
python scripts/evaluate_mapping.py --method both
```

## Output Files

- **Processed .h5ad**: Original data + `cl_id`, `cl_label`, `original_cell_type` in obs
- **JSON Report**: Summary statistics, gene/cell type mapping rates, quality metrics
- **TSV Mapping Details**: Per-gene and per-cell-type mapping results

## Notes

- Fully offline — no network dependency for mapping operations
- Currently supports murine (Mus musculus) datasets
- .h5ad (AnnData) format only
- Recommended: files under 1GB for optimal performance

## License

MIT License

## Documentation (Chinese)

For a full Chinese usage guide (setup, ports, CLI, data layout), see **[使用说明.md](使用说明.md)**.
