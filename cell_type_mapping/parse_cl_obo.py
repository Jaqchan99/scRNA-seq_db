"""
解析 Cell Ontology OBO 文件，生成细胞类型映射表 CSV
"""
import re
import csv
import os

def parse_obo_file(obo_path: str) -> list:
    """
    解析 OBO 文件，提取 CL 条目的 ID、名称和同义词
    
    Returns:
        list of dict: [{'cl_id': 'CL:0000084', 'cl_label': 'T cell', 'synonyms': ['T-cell', 'T lymphocyte', ...]}]
    """
    entries = []
    current_entry = None
    
    with open(obo_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            
            # 新的 Term 开始
            if line == '[Term]':
                if current_entry and current_entry.get('cl_id', '').startswith('CL:'):
                    entries.append(current_entry)
                current_entry = {'cl_id': '', 'cl_label': '', 'synonyms': []}
                continue
            
            # 其他块类型，保存当前条目
            if line.startswith('[') and line.endswith(']'):
                if current_entry and current_entry.get('cl_id', '').startswith('CL:'):
                    entries.append(current_entry)
                current_entry = None
                continue
            
            if current_entry is None:
                continue
            
            # 解析 ID
            if line.startswith('id: '):
                current_entry['cl_id'] = line[4:].strip()
            
            # 解析 name
            elif line.startswith('name: '):
                name = line[6:].strip()
                # 跳过 obsolete 条目
                if name.startswith('obsolete '):
                    current_entry['cl_id'] = ''  # 标记为无效
                else:
                    current_entry['cl_label'] = name
            
            # 解析 synonym
            elif line.startswith('synonym: '):
                # 格式: synonym: "xxx" TYPE [...]
                match = re.match(r'synonym:\s*"([^"]+)"\s*(\w+)', line)
                if match:
                    syn_text = match.group(1).strip()
                    syn_type = match.group(2)  # EXACT, RELATED, NARROW, BROAD
                    # 只保留 EXACT 和 RELATED 同义词
                    if syn_type in ['EXACT', 'RELATED']:
                        if syn_text and syn_text != current_entry['cl_label']:
                            current_entry['synonyms'].append(syn_text)
    
    # 保存最后一个条目
    if current_entry and current_entry.get('cl_id', '').startswith('CL:'):
        entries.append(current_entry)
    
    return entries


def generate_csv(entries: list, output_path: str):
    """
    生成 CSV 文件
    格式: synonyms (@ 分隔), cl_id, cl_label
    """
    with open(output_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['synonyms', 'cl_id', 'cl_label'])
        
        for entry in entries:
            if not entry['cl_label']:  # 跳过无名称的条目
                continue
            
            # 构建同义词列表：标准名称 + 所有同义词
            all_names = [entry['cl_label']] + entry['synonyms']
            # 去重，保持顺序
            seen = set()
            unique_names = []
            for name in all_names:
                name_lower = name.lower()
                if name_lower not in seen:
                    seen.add(name_lower)
                    unique_names.append(name)
            
            synonyms_str = '@'.join(unique_names)
            writer.writerow([synonyms_str, entry['cl_id'], entry['cl_label']])


def add_common_abbreviations(entries: list) -> list:
    """
    添加常见的缩写和变体
    """
    abbreviations = {
        'T cell': ['T-cell', 'Tcell', 'T cells', 'T-cells'],
        'B cell': ['B-cell', 'Bcell', 'B cells', 'B-cells'],
        'natural killer cell': ['NK cell', 'NK-cell', 'NK', 'NKcell'],
        'dendritic cell': ['DC', 'DCs'],
        'macrophage': ['Mφ', 'MΦ'],
        'regulatory T cell': ['Treg', 'Tregs', 'T reg', 'T-reg'],
        'helper T cell': ['Th cell', 'Th'],
        'cytotoxic T cell': ['CTL', 'Tc cell'],
        'memory T cell': ['Tm cell', 'memory T'],
        'naive T cell': ['naive T', 'naïve T cell'],
        'effector T cell': ['Teff', 'effector T'],
        'CD4-positive, alpha-beta T cell': ['CD4+ T cell', 'CD4 T cell', 'CD4+', 'CD4 positive T cell'],
        'CD8-positive, alpha-beta T cell': ['CD8+ T cell', 'CD8 T cell', 'CD8+', 'CD8 positive T cell'],
        'monocyte': ['mono', 'monocytes'],
        'neutrophil': ['neutrophils', 'PMN'],
        'eosinophil': ['eosinophils', 'eos'],
        'basophil': ['basophils', 'baso'],
        'mast cell': ['mast cells', 'mastocyte'],
        'fibroblast': ['fibroblasts', 'FB'],
        'endothelial cell': ['EC', 'endothelial', 'endothelial cells'],
        'epithelial cell': ['epithelial', 'epithelial cells'],
        'erythrocyte': ['RBC', 'red blood cell', 'red blood cells'],
        'platelet': ['PLT', 'platelets', 'thrombocyte'],
        'hepatocyte': ['hepatocytes'],
        'neuron': ['neurons', 'nerve cell'],
        'astrocyte': ['astrocytes'],
        'oligodendrocyte': ['oligodendrocytes', 'oligo'],
        'microglia': ['microglial cell'],
        'stem cell': ['SC', 'stem cells'],
        'hematopoietic stem cell': ['HSC', 'HSCs'],
        'mesenchymal stem cell': ['MSC', 'MSCs'],
        'progenitor cell': ['progenitor', 'progenitors'],
        'plasma cell': ['plasma cells', 'plasmacyte'],
    }
    
    for entry in entries:
        label_lower = entry['cl_label'].lower()
        for standard_name, abbrevs in abbreviations.items():
            if label_lower == standard_name.lower():
                for abbr in abbrevs:
                    if abbr.lower() not in [s.lower() for s in entry['synonyms']]:
                        entry['synonyms'].append(abbr)
                break
    
    return entries


def main():
    # 路径配置
    script_dir = os.path.dirname(os.path.abspath(__file__))
    obo_path = os.path.join(script_dir, 'cl.obo')
    output_path = os.path.join(script_dir, 'cell_type_mapping.csv')
    
    print(f"📂 解析 OBO 文件: {obo_path}")
    
    # 解析 OBO 文件
    entries = parse_obo_file(obo_path)
    print(f"✅ 解析完成，共 {len(entries)} 个有效 CL 条目")
    
    # 添加常见缩写
    entries = add_common_abbreviations(entries)
    print(f"✅ 添加常见缩写完成")
    
    # 生成 CSV
    generate_csv(entries, output_path)
    print(f"✅ CSV 文件已生成: {output_path}")
    
    # 统计信息
    total_synonyms = sum(len(e['synonyms']) for e in entries)
    print(f"\n📊 统计信息:")
    print(f"   - CL 条目数: {len(entries)}")
    print(f"   - 同义词总数: {total_synonyms}")
    print(f"   - 平均每条目同义词数: {total_synonyms / len(entries):.2f}")


if __name__ == '__main__':
    main()
