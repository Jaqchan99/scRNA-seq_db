// ═══════════════════════════════════════════════════════════════
//  scRNA-seq Platform  —  script.js  v2.0
// ═══════════════════════════════════════════════════════════════
// UI build id — 若控制台看不到此行，说明浏览器仍在使用缓存的旧 script.js
const PLATFORM_UI_BUILD = '20260210-phase2c';

const API_BASE_URL = 'http://localhost:8100';
/** GitHub 新建 Issue（可改为贵仓库的 annotation 模板链接） */
const FEEDBACK_ISSUE_URL = 'https://github.com/Jaqchan99/scRNA-seq_db/issues/new';

let currentSubmissionId = null;
let currentStep = 1;
let cellTypeChart = null;

// Chart.js colour palette (blue-toned, academic)
const CHART_COLORS = [
    '#1e3a5f','#2563eb','#0284c7','#0891b2','#0d9488',
    '#16a34a','#65a30d','#ca8a04','#d97706','#dc2626',
    '#9333ea','#7c3aed','#db2777','#e11d48','#475569',
    '#0f766e','#1d4ed8','#4338ca','#6d28d9','#be185d',
];

// ── Page navigation ──────────────────────────────────────────
function showPage(name) {
    document.querySelectorAll('.page').forEach(p => p.classList.remove('active'));
    document.querySelectorAll('.nav-link').forEach(l => l.classList.remove('active'));
    document.getElementById('page-' + name).classList.add('active');
    document.getElementById('nav-' + name).classList.add('active');
}

// ── Step navigation ──────────────────────────────────────────
function goToStep(n) {
    document.querySelectorAll('.step-content').forEach(c => c.classList.remove('active'));
    document.getElementById('content-' + n).classList.add('active');

    document.querySelectorAll('.step').forEach((s, i) => {
        s.classList.remove('active', 'completed');
        if (i + 1 < n)  s.classList.add('completed');
        if (i + 1 === n) s.classList.add('active');
    });

    currentStep = n;
}

// ── Status bar ───────────────────────────────────────────────
function updateStatus(msg, busy = true) {
    document.getElementById('currentStatus').textContent = msg;
    const dot = document.getElementById('statusDot');
    dot.className = 'status-dot' + (busy ? ' busy' : '');
}

function setStatusId(id) {
    const el = document.getElementById('statusId');
    el.textContent = id ? 'ID: ' + id.slice(0, 8) + '…' : '';
}

// ── Error modal ──────────────────────────────────────────────
function showError(msg) {
    document.getElementById('errorMessage').textContent = msg;
    document.getElementById('errorModal').style.display = 'flex';
    updateStatus('Error — see modal', false);
}

function closeErrorModal() {
    document.getElementById('errorModal').style.display = 'none';
}

// ── Metric card builder ──────────────────────────────────────
function metricCard(label, value, sub, variant) {
    const cls = variant ? ' ' + variant : '';
    return `<div class="metric-card${cls}">
        <div class="metric-label">${label}</div>
        <div class="metric-value">${value}</div>
        ${sub ? `<div class="metric-sub">${sub}</div>` : ''}
    </div>`;
}

// ── Utility ──────────────────────────────────────────────────
function escapeHtml(s) {
    if (s == null) return '';
    return String(s)
        .replace(/&/g, '&amp;')
        .replace(/</g, '&lt;')
        .replace(/"/g, '&quot;');
}

function formatFileSize(bytes) {
    if (bytes < 1024)       return bytes + ' B';
    if (bytes < 1048576)    return (bytes / 1024).toFixed(1) + ' KB';
    if (bytes < 1073741824) return (bytes / 1048576).toFixed(1) + ' MB';
    return (bytes / 1073741824).toFixed(2) + ' GB';
}

/** Final export JSON nests full arrays under report.mapping_details — not under *.mapping_details */
function getGeneMappingDetails(report) {
    const gm = report.gene_mapping || {};
    if (gm.mapping_details && gm.mapping_details.length) return gm.mapping_details;
    return (report.mapping_details && report.mapping_details.gene_mapping) || [];
}

function getCellTypeMappingDetails(report) {
    const cm = report.cell_type_mapping || {};
    if (cm.mapping_details && cm.mapping_details.length) return cm.mapping_details;
    return (report.mapping_details && report.mapping_details.cell_type_mapping) || [];
}

/** Step 2：细胞类型 obs 列选择（Phase2） */
function renderCellTypeColumnPanel(cti) {
    const host = document.getElementById('cellTypeColumnHost');
    if (!host) return;
    if (!cti) {
        host.innerHTML = '';
        return;
    }
    const req = cti.requires_user_selection;
    const rec = cti.recommended || '';
    const cand = cti.candidates || [];
    let opts = `<option value="">${rec ? `Auto (recommended: ${escapeHtml(rec)})` : 'Auto-detect (recommended column if any)'}</option>`;
    cand.forEach((c) => {
        const ex = (c.examples || []).join(', ').slice(0, 60);
        opts += `<option value="${escapeHtml(c.column)}">${escapeHtml(c.column)} — ${c.n_unique} labels${ex ? ' · e.g. ' + escapeHtml(ex) : ''}</option>`;
    });
    const warn = req
        ? '<p class="msg-warning">⚠ Please confirm which <code>obs</code> column holds cell type text labels before mapping.</p>'
        : '<p class="msg-success" style="margin-bottom:10px;">✓ A default column was chosen; change below if needed.</p>';
    host.innerHTML = `
        <div class="detail-panel cell-type-col-panel">
            <p style="font-weight:600;color:#1e293b;margin-bottom:6px;">Cell type column</p>
            ${warn}
            <label class="sr-only" for="cellTypeColumnSelect">Pick column</label>
            <select id="cellTypeColumnSelect" class="select-input cell-type-col-select">${opts}</select>
            <p style="margin-top:10px;font-size:0.85em;color:#64748b;">Or type exact <code>obs</code> column name (overrides list):</p>
            <input type="text" id="cellTypeColumnCustom" class="select-input cell-type-col-input" placeholder="e.g. cell_type" autocomplete="off" />
        </div>`;
}

function openFeedbackIssue() {
    const title = encodeURIComponent('[Platform] Cell type annotation / unmapped labels');
    const body = encodeURIComponent(
        `Short submission ID: ${currentSubmissionId ? currentSubmissionId.slice(0, 8) : 'N/A'}\n\n` +
            `Describe your issue or attach the reviewed CSV (reviewed=yes).\n`
    );
    window.open(`${FEEDBACK_ISSUE_URL}?title=${title}&body=${body}`, '_blank', 'noopener');
}

/** 映射率展示：兼容 mapping_rate 异常大于 100 的旧报告 */
function formatMappingRatePct(block) {
    const r = block.success_rate ?? block.mapping_rate;
    if (r != null && r <= 100 && r >= 0) return Number(r).toFixed(1) + '%';
    const map = Number(block.mapped_count) || 0;
    const un = Number(block.unmapped_count) || 0;
    const tot = map + un;
    if (tot > 0) return ((map / tot) * 100).toFixed(1) + '%';
    return '—';
}

/** Step 4 顶部单行：Status / Cells / Genes / Sparsity + 一行映射率说明 */
function buildFinalSummaryHtml(report) {
    const sum = report.summary || {};
    const qm = report.quality_metrics || {};
    const gm = report.gene_mapping || {};
    const cm = report.cell_type_mapping || {};

    const cells = sum.n_cells ?? qm.n_cells;
    const genes = sum.n_genes ?? qm.n_genes;
    const cellsStr = cells != null ? Number(cells).toLocaleString() : '—';
    const genesStr = genes != null ? Number(genes).toLocaleString() : '—';

    let sparsityStr = '—';
    if (qm.sparsity != null) {
        const s = Number(qm.sparsity);
        sparsityStr = (s <= 1 ? (s * 100) : s).toFixed(1) + '%';
    }

    const gmPct = formatMappingRatePct(gm);
    const cmPct =
        cm.cell_mapping_rate != null
            ? Number(cm.cell_mapping_rate).toFixed(1) + '% (cells)'
            : formatMappingRatePct(cm);

    return `
        <div class="final-summary-grid">
            <div class="metric-card metric-card--compact metric-card--status success">
                <div class="metric-label">Status</div>
                <div class="metric-value metric-value--icon" title="Completed">✓</div>
                <div class="metric-sub">Completed</div>
            </div>
            <div class="metric-card metric-card--compact success">
                <div class="metric-label">Cells</div>
                <div class="metric-value">${cellsStr}</div>
            </div>
            <div class="metric-card metric-card--compact success">
                <div class="metric-label">Genes</div>
                <div class="metric-value">${genesStr}</div>
            </div>
            <div class="metric-card metric-card--compact success">
                <div class="metric-label">Sparsity</div>
                <div class="metric-value">${sparsityStr}</div>
            </div>
        </div>
        <p class="mapping-inline">Gene mapping <strong>${gmPct}</strong> · Cell-level CL coverage <strong>${cmPct}</strong></p>
    `;
}

// ── Init ─────────────────────────────────────────────────────
document.addEventListener('DOMContentLoaded', () => {
    document.getElementById('fileInput').addEventListener('change', handleFileSelect);

    const area = document.getElementById('uploadArea');
    // 勿与内部 Browse 按钮叠加触发，否则会连续打开两次文件选择器
    area.addEventListener('click', (e) => {
        if (e.target.closest('#browseFileBtn')) return;
        document.getElementById('fileInput').click();
    });
    document.getElementById('browseFileBtn').addEventListener('click', (e) => {
        e.stopPropagation();
        document.getElementById('fileInput').click();
    });
    area.addEventListener('dragover',  e => { e.preventDefault(); area.classList.add('dragover'); });
    area.addEventListener('dragleave', e => { e.preventDefault(); area.classList.remove('dragover'); });
    area.addEventListener('drop', e => {
        e.preventDefault();
        area.classList.remove('dragover');
        if (e.dataTransfer.files.length) handleFile(e.dataTransfer.files[0]);
    });

    document.getElementById('uploadBtn').addEventListener('click', handleUpload);
    document.getElementById('continueToMapping').addEventListener('click', startMapping);
    document.getElementById('cancelProcess').addEventListener('click', resetProcess);
    document.getElementById('continueToExport').addEventListener('click', startExport);
    document.getElementById('cancelProcess2').addEventListener('click', resetProcess);
    document.getElementById('resetBtn').addEventListener('click', resetProcess);

    updateStatus('Ready', false);
    console.info('[scRNA platform] UI build', PLATFORM_UI_BUILD, '— API', API_BASE_URL);
});

// ── File selection ────────────────────────────────────────────
function handleFileSelect(e) {
    if (e.target.files[0]) handleFile(e.target.files[0]);
}

function handleFile(file) {
    if (!file.name.endsWith('.h5ad')) {
        showError('Please select a file in .h5ad (AnnData) format.');
        return;
    }
    document.getElementById('fileName').textContent = file.name;
    document.getElementById('fileSize').textContent = formatFileSize(file.size);
    document.getElementById('fileInfo').style.display = 'block';
    const btn = document.getElementById('uploadBtn');
    btn.disabled = false;
    btn._file = file;
}

// ── Step 1: Upload ────────────────────────────────────────────
async function handleUpload() {
    const file = document.getElementById('uploadBtn')._file;
    if (!file) return;

    try {
        updateStatus('Creating submission…');

        const subRes = await fetch(`${API_BASE_URL}/submissions`, { method: 'POST' });
        if (!subRes.ok) throw new Error('Failed to create submission task');
        const subData = await subRes.json();
        currentSubmissionId = subData.submission_id;
        setStatusId(currentSubmissionId);

        updateStatus('Uploading file…');

        const formData = new FormData();
        formData.append('file', file);

        // Attach disease type as a custom header (informational; backend stores it if supported)
        const diseaseType = document.getElementById('diseaseType').value;
        const uploadRes = await fetch(
            `${API_BASE_URL}/submissions/${currentSubmissionId}/upload`,
            { method: 'POST', body: formData }
        );
        if (!uploadRes.ok) {
            const txt = await uploadRes.text();
            throw new Error('Upload failed: ' + txt);
        }

        goToStep(2);
        startValidation();

    } catch (err) {
        showError(err.message);
    }
}

// ── Step 2: Validation ────────────────────────────────────────
async function startValidation() {
    updateStatus('Validating dataset…');

    try {
        const res = await fetch(
            `${API_BASE_URL}/submissions/${currentSubmissionId}/validate`,
            { method: 'POST' }
        );
        if (!res.ok) throw new Error('Failed to start validation');

        const timer = setInterval(async () => {
            try {
                const s = await (await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/status`)).json();
                if (s.status === 'validated') {
                    clearInterval(timer);
                    await showValidationResults();
                    document.getElementById('validationConfirmation').style.display = 'flex';
                    updateStatus('Validation complete — awaiting confirmation', false);
                } else if (s.status === 'validation_failed' || s.status === 'failed') {
                    clearInterval(timer);
                    showError('Validation failed: ' + (s.error_message || 'Unknown error'));
                }
            } catch (e) { clearInterval(timer); showError(e.message); }
        }, 2000);

    } catch (err) { showError(err.message); }
}

async function showValidationResults() {
    try {
        const report = await (await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/report`)).json();
        const summary = report.summary || {};
        const val     = report.validation_result || {};

        // Metrics
        const fileSizeMB = summary.file_size_mb ? Number(summary.file_size_mb).toFixed(2) : '0.00';
        document.getElementById('resultSummary').innerHTML =
            metricCard('Cells',      (summary.n_cells  || 0).toLocaleString(), null, 'success') +
            metricCard('Genes',      (summary.n_genes  || 0).toLocaleString(), null, 'success') +
            metricCard('File Size',  fileSizeMB + ' MB',      null, '');

        // Details
        const errors   = val.errors   || [];
        const warnings = val.warnings || [];
        let html = '';
        if (errors.length === 0 && warnings.length === 0) {
            html = '<p class="msg-success">✓ Validation passed — no errors or warnings found.</p>';
        }
        errors.forEach(e   => { html += `<p class="msg-error">✗ ${e}</p>`; });
        warnings.forEach(w => { html += `<p class="msg-warning">⚠ ${w}</p>`; });
        document.getElementById('resultDetails').innerHTML = html;

        renderCellTypeColumnPanel(report.cell_type_column);

        document.getElementById('validationStatus').style.display = 'none';
        document.getElementById('validationResults').style.display = 'block';

    } catch (err) { showError('Failed to load validation report: ' + err.message); }
}

// ── Step 3: Mapping ───────────────────────────────────────────
async function startMapping() {
    goToStep(3);
    updateStatus('Running gene & cell type mapping…');

    try {
        const custom = (document.getElementById('cellTypeColumnCustom') || {}).value?.trim() || '';
        const sel = (document.getElementById('cellTypeColumnSelect') || {}).value || '';
        const col = custom || sel || null;
        const payload = col ? { cell_type_column: col } : {};

        const res = await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/continue`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(payload),
        });
        if (!res.ok) throw new Error('Failed to start mapping');

        const timer = setInterval(async () => {
            try {
                const s = await (await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/status`)).json();
                if (s.status === 'completed') {
                    clearInterval(timer);
                    await showMappingResults();
                    document.getElementById('mappingConfirmation').style.display = 'flex';
                    updateStatus('Mapping complete — awaiting confirmation', false);
                } else if (s.status === 'failed') {
                    clearInterval(timer);
                    showError('Mapping failed: ' + (s.error_message || 'Unknown error'));
                }
            } catch (e) { clearInterval(timer); showError(e.message); }
        }, 2000);

    } catch (err) { showError(err.message); }
}

async function showMappingResults() {
    try {
        const report   = await (await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/report`)).json();
        const gm       = report.gene_mapping      || {};
        const cm       = report.cell_type_mapping || {};

        const gmRate   = gm.success_rate  ?? gm.mapping_rate  ?? 0;
        const gmRateStr = Number(gmRate).toFixed(1) + '%';
        
        // Cell type mapping: 优先使用后端返回的 cell_mapping_rate，否则前端计算
        let cmCellRate, mappedCells, totalCells;
        if (cm.cell_mapping_rate != null && cm.mapped_cells != null && cm.total_cells != null) {
            cmCellRate = cm.cell_mapping_rate;
            mappedCells = cm.mapped_cells;
            totalCells = cm.total_cells;
        } else {
            // 兼容旧报告：前端计算
            const cmDetails = getCellTypeMappingDetails(report);
            mappedCells = 0;
            totalCells = 0;
            cmDetails.forEach(d => {
                const count = Number(d.cell_count) || 1;
                totalCells += count;
                if (d.status === 'mapped') mappedCells += count;
            });
            cmCellRate = totalCells > 0 ? ((mappedCells / totalCells) * 100) : 0;
        }
        const cmCellRateStr = cmCellRate.toFixed(1) + '%';

        // Metrics
        document.getElementById('mappingSummary').innerHTML =
            metricCard('Genes Mapped',      (gm.mapped_count   || 0).toLocaleString(), gmRateStr + ' success', 'success') +
            metricCard('Cells Mapped',      mappedCells.toLocaleString(), cmCellRateStr + ' success', 'success') +
            metricCard('Unmapped Genes',    (gm.unmapped_count || 0).toLocaleString(), null, gm.unmapped_count > 0 ? 'warning' : '') +
            metricCard('Cell Types Mapped', (cm.mapped_count   || 0) + ' / ' + (cm.total_unique_types || cm.total_types || 0), null, '');

        // Details table (top 10 gene mapping)
        const gd = getGeneMappingDetails(report).slice(0, 10);
        let html = '';
        if (gd.length) {
            html += '<p style="font-weight:600;color:#1e293b;margin-bottom:10px;">Gene Mapping — Top 10 entries</p>';
            html += '<div class="table-wrapper"><table class="data-table"><thead><tr><th>Original Symbol</th><th>Mapped Symbol</th><th>Status</th></tr></thead><tbody>';
            gd.forEach(r => {
                const ok = r.status === 'mapped';
                html += `<tr><td>${r.gene_symbol||''}</td><td>${r.mapped_symbol||'—'}</td>
                    <td><span class="badge ${ok?'badge-success':'badge-warning'}">${r.status||''}</span></td></tr>`;
            });
            html += '</tbody></table></div>';
        }
        document.getElementById('mappingDetails').innerHTML = html;

        document.getElementById('mappingStatus').style.display = 'none';
        document.getElementById('mappingResults').style.display = 'block';

    } catch (err) { showError('Failed to load mapping report: ' + err.message); }
}

// ── Step 4: Export & Results ──────────────────────────────────
async function startExport() {
    goToStep(4);
    updateStatus('Generating output files…');

    try {
        // Results are already generated; just load the final report
        const report = await (await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/report`)).json();
        await showFinalResults(report);
        updateStatus('Processing complete', false);
    } catch (err) { showError(err.message); }
}

async function showFinalResults(report) {
    const details = getCellTypeMappingDetails(report);

    document.getElementById('finalSummary').innerHTML = buildFinalSummaryHtml(report);

    document.getElementById('downloadReport').onclick = () => downloadFile(
        `${API_BASE_URL}/submissions/${currentSubmissionId}/report`,
        `report_${currentSubmissionId.slice(0,8)}.json`
    );
    document.getElementById('downloadData').onclick = () => downloadFile(
        `${API_BASE_URL}/submissions/${currentSubmissionId}/export`,
        `processed_${currentSubmissionId.slice(0,8)}.h5ad`
    );

    const unmappedN = details.filter((d) => d.status !== 'mapped').length;
    const unmappedRow = document.getElementById('unmappedDownloadRow');
    if (unmappedRow) {
        unmappedRow.style.display = unmappedN > 0 ? '' : 'none';
    }
    const du = document.getElementById('downloadUnmapped');
    if (du) {
        du.onclick = () =>
            downloadFile(
                `${API_BASE_URL}/submissions/${currentSubmissionId}/unmapped-labels`,
                `unmapped_labels_${currentSubmissionId.slice(0, 8)}.csv`
            );
    }
    const fb = document.getElementById('openFeedback');
    if (fb) fb.onclick = () => openFeedbackIssue();

    // 先显示面板，再画图：否则父级 display:none 时 Canvas 宽高为 0，Chart.js 不可见
    document.getElementById('resultsContainer').style.display = 'none';
    document.getElementById('finalResults').style.display = 'block';

    if (details.length > 0) {
        document.getElementById('chartSection').style.display = 'block';
        document.getElementById('cellTypeTableSection').style.display = 'block';
        requestAnimationFrame(() => {
            requestAnimationFrame(() => {
                buildCellTypeChart(details);
                buildCellTypeTable(details);
                if (cellTypeChart && typeof cellTypeChart.resize === 'function') {
                    cellTypeChart.resize();
                }
            });
        });
    } else {
        document.getElementById('chartSection').style.display = 'none';
        document.getElementById('cellTypeTableSection').style.display = 'none';
    }
}

// ── Cell Type Distribution Chart ──────────────────────────────
function buildCellTypeChart(details) {
    // Aggregate: sum cell_count per CL label (backend sends cell_count per unique type)
    const counts = {};
    details.forEach(d => {
        const label = (d.cl_label && d.cl_label !== 'None' && d.cl_label !== 'null')
            ? d.cl_label : d.raw_cell_type;
        if (!label) return;
        const n = Number(d.cell_count);
        counts[label] = (counts[label] || 0) + (Number.isFinite(n) && n > 0 ? n : 1);
    });

    // Sort by count descending, cap at 20 for readability
    const sorted = Object.entries(counts)
        .sort((a, b) => b[1] - a[1])
        .slice(0, 20);

    if (sorted.length === 0) {
        if (cellTypeChart) {
            cellTypeChart.destroy();
            cellTypeChart = null;
        }
        document.getElementById('chartLegend').innerHTML =
            '<p class="chart-empty" style="color:#94a3b8;font-size:0.88em;">No cell types to plot.</p>';
        return;
    }

    const labels = sorted.map(([l]) => l.length > 28 ? l.slice(0, 26) + '…' : l);
    const data   = sorted.map(([, v]) => v);
    const colors = sorted.map((_, i) => CHART_COLORS[i % CHART_COLORS.length]);

    const ctx = document.getElementById('cellTypeChart').getContext('2d');
    if (cellTypeChart) cellTypeChart.destroy();

    cellTypeChart = new Chart(ctx, {
        type: 'bar',
        data: {
            labels,
            datasets: [{
                label: 'Cell types',
                data,
                backgroundColor: colors.map(c => c + 'cc'),
                borderColor: colors,
                borderWidth: 1.5,
                borderRadius: 4,
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                legend: { display: false },
                tooltip: {
                    callbacks: {
                        title: (items) => sorted[items[0].dataIndex][0],
                        label: (item) => ' Cells: ' + item.raw,
                    }
                }
            },
            scales: {
                x: {
                    ticks: {
                        font: { size: 11 },
                        maxRotation: 40,
                        color: '#64748b',
                    },
                    grid: { display: false },
                },
                y: {
                    beginAtZero: true,
                    ticks: { font: { size: 11 }, color: '#64748b', precision: 0 },
                    grid: { color: '#f1f5f9' },
                    title: { display: true, text: 'Cells', color: '#94a3b8', font: { size: 11 } }
                }
            }
        }
    });

    // Legend
    const legendEl = document.getElementById('chartLegend');
    legendEl.innerHTML = sorted.slice(0, 10).map(([lbl], i) =>
        `<div class="legend-item">
            <div class="legend-dot" style="background:${colors[i]}"></div>
            <span>${lbl.length > 22 ? lbl.slice(0,20)+'…' : lbl}</span>
        </div>`
    ).join('');
}

// ── Cell Type Table ───────────────────────────────────────────
function buildCellTypeTable(details) {
    const tbody = document.getElementById('cellTypeTableBody');
    tbody.innerHTML = details.map(d => {
        const mapped = d.status === 'mapped';
        const conf   = d.confidence || '';
        const confBadge = conf === 'high'   ? 'badge-success'
                        : conf === 'medium' ? 'badge-info'
                        : conf === 'low'    ? 'badge-warning'
                        : 'badge-info';
        return `<tr>
            <td>${d.raw_cell_type || '—'}</td>
            <td>${d.cl_id ? `<span class="cl-id">${d.cl_id}</span>` : '—'}</td>
            <td>${d.cl_label || '—'}</td>
            <td>${conf ? `<span class="badge ${confBadge}">${conf}</span>` : '—'}</td>
            <td><span class="badge ${mapped ? 'badge-success' : 'badge-warning'}">${d.status||'—'}</span></td>
        </tr>`;
    }).join('');
}

// ── Download helper ───────────────────────────────────────────
function downloadFile(url, filename) {
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
}

// ── Reset ─────────────────────────────────────────────────────
function resetProcess() {
    currentSubmissionId = null;
    currentStep = 1;

    if (cellTypeChart) { cellTypeChart.destroy(); cellTypeChart = null; }

    // Reset step UI
    document.querySelectorAll('.step').forEach((s, i) => {
        s.classList.remove('active', 'completed');
        if (i === 0) s.classList.add('active');
    });
    document.querySelectorAll('.step-content').forEach(c => c.classList.remove('active'));
    document.getElementById('content-1').classList.add('active');

    // Reset Step 1
    document.getElementById('fileInfo').style.display = 'none';
    document.getElementById('uploadBtn').disabled = true;
    document.getElementById('fileInput').value = '';

    // Reset Step 2
    document.getElementById('validationStatus').style.display = 'block';
    document.getElementById('validationResults').style.display = 'none';
    document.getElementById('validationConfirmation').style.display = 'none';

    const cth = document.getElementById('cellTypeColumnHost');
    if (cth) cth.innerHTML = '';

    // Reset Step 3
    document.getElementById('mappingStatus').style.display = 'block';
    document.getElementById('mappingResults').style.display = 'none';
    document.getElementById('mappingConfirmation').style.display = 'none';

    // Reset Step 4
    document.getElementById('resultsContainer').style.display = 'block';
    document.getElementById('finalResults').style.display = 'none';
    document.getElementById('chartSection').style.display = 'none';
    document.getElementById('cellTypeTableSection').style.display = 'none';

    updateStatus('Ready', false);
    setStatusId(null);
}
