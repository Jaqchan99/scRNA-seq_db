// 全局变量
let currentSubmissionId = null;
let currentStep = 1;
const API_BASE_URL = 'http://localhost:8000';

// DOM 元素
const stepElements = document.querySelectorAll('.step');
const contentElements = document.querySelectorAll('.step-content');
const fileInput = document.getElementById('fileInput');
const uploadArea = document.getElementById('uploadArea');
const fileInfo = document.getElementById('fileInfo');
const uploadBtn = document.getElementById('uploadBtn');
const resetBtn = document.getElementById('resetBtn');
const errorModal = document.getElementById('errorModal');
const errorMessage = document.getElementById('errorMessage');

// 初始化
document.addEventListener('DOMContentLoaded', function() {
    initializeEventListeners();
    updateStatus('准备就绪');
});

// 事件监听器
function initializeEventListeners() {
    // 文件上传
    fileInput.addEventListener('change', handleFileSelect);
    uploadArea.addEventListener('click', () => fileInput.click());
    uploadArea.addEventListener('dragover', handleDragOver);
    uploadArea.addEventListener('dragleave', handleDragLeave);
    uploadArea.addEventListener('drop', handleDrop);
    uploadBtn.addEventListener('click', handleUpload);
    resetBtn.addEventListener('click', resetProcess);
}

// 文件选择处理
function handleFileSelect(event) {
    const file = event.target.files[0];
    if (file) {
        processFile(file);
    }
}

// 拖拽处理
function handleDragOver(event) {
    event.preventDefault();
    uploadArea.classList.add('dragover');
}

function handleDragLeave(event) {
    event.preventDefault();
    uploadArea.classList.remove('dragover');
}

function handleDrop(event) {
    event.preventDefault();
    uploadArea.classList.remove('dragover');
    const files = event.dataTransfer.files;
    if (files.length > 0) {
        processFile(files[0]);
    }
}

// 处理文件
function processFile(file) {
    if (!file.name.endsWith('.h5ad')) {
        showError('请选择 .h5ad 格式的文件');
        return;
    }
    
    // 显示文件信息
    document.getElementById('fileName').textContent = file.name;
    document.getElementById('fileSize').textContent = formatFileSize(file.size);
    fileInfo.style.display = 'block';
    uploadBtn.disabled = false;
    
    // 存储文件引用
    uploadBtn.file = file;
}

// 文件上传
async function handleUpload() {
    const file = uploadBtn.file;
    if (!file) return;
    
    try {
        updateStatus('正在创建提交任务...');
        
        // 1. 创建提交任务
        const submissionResponse = await fetch(`${API_BASE_URL}/submissions`, {
            method: 'POST'
        });
        
        if (!submissionResponse.ok) {
            throw new Error('创建提交任务失败');
        }
        
        const submissionData = await submissionResponse.json();
        currentSubmissionId = submissionData.submission_id;
        
        updateStatus('正在上传文件...');
        
        // 2. 上传文件
        const formData = new FormData();
        formData.append('file', file);
        
        const uploadResponse = await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/upload`, {
            method: 'POST',
            body: formData
        });
        
        if (!uploadResponse.ok) {
            const errorText = await uploadResponse.text();
            throw new Error(`文件上传失败: ${errorText}`);
        }
        
        const uploadData = await uploadResponse.json();
        console.log('上传成功:', uploadData);
        
        // 3. 进入校验步骤
        goToStep(2);
        startValidation();
        
    } catch (error) {
        console.error('上传失败:', error);
        showError(`上传失败: ${error.message}`);
    }
}

// 开始校验
async function startValidation() {
    try {
        updateStatus('正在校验数据...');
        
        // 轮询检查状态
        const checkInterval = setInterval(async () => {
            try {
                const statusResponse = await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/status`);
                if (!statusResponse.ok) {
                    throw new Error('状态查询失败');
                }
                
                const statusData = await statusResponse.json();
                console.log('当前状态:', statusData);
                
                if (statusData.status === 'completed') {
                    clearInterval(checkInterval);
                    showValidationResults();
                    goToStep(3);
                    startMapping();
                } else if (statusData.status === 'failed' || statusData.status === 'validation_failed') {
                    clearInterval(checkInterval);
                    showError(`处理失败: ${statusData.error_message || '未知错误'}`);
                }
                // 继续等待 processing 状态
                
            } catch (error) {
                clearInterval(checkInterval);
                console.error('状态检查失败:', error);
                showError(`状态检查失败: ${error.message}`);
            }
        }, 2000); // 每2秒检查一次
        
    } catch (error) {
        console.error('校验失败:', error);
        showError(`校验失败: ${error.message}`);
    }
}

// 显示校验结果
async function showValidationResults() {
    try {
        const reportResponse = await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/report`);
        if (!reportResponse.ok) {
            throw new Error('获取报告失败');
        }
        
        const reportData = await reportResponse.json();
        console.log('校验报告:', reportData);
        
        // 显示校验结果
        const validationResults = document.getElementById('validationResults');
        const resultSummary = document.getElementById('resultSummary');
        const resultDetails = document.getElementById('resultDetails');
        
        // 创建摘要卡片
        const summary = reportData.summary || {};
        const validation = reportData.validation || {};
        
        resultSummary.innerHTML = `
            <div class="summary-card success">
                <h4>细胞数</h4>
                <div class="number">${summary.n_cells || 0}</div>
            </div>
            <div class="summary-card success">
                <h4>基因数</h4>
                <div class="number">${summary.n_genes || 0}</div>
            </div>
            <div class="summary-card success">
                <h4>文件大小</h4>
                <div class="number">${summary.file_size_mb || 0} MB</div>
            </div>
        `;
        
        // 显示校验详情
        const errors = validation.errors || [];
        const warnings = validation.warnings || [];
        
        let detailsHtml = '';
        if (errors.length > 0) {
            detailsHtml += '<h4 style="color: #e74c3c;">错误:</h4><ul>';
            errors.forEach(error => {
                detailsHtml += `<li style="color: #e74c3c;">${error}</li>`;
            });
            detailsHtml += '</ul>';
        }
        
        if (warnings.length > 0) {
            detailsHtml += '<h4 style="color: #f39c12;">警告:</h4><ul>';
            warnings.forEach(warning => {
                detailsHtml += `<li style="color: #f39c12;">${warning}</li>`;
            });
            detailsHtml += '</ul>';
        }
        
        if (errors.length === 0 && warnings.length === 0) {
            detailsHtml = '<p style="color: #27ae60;">✓ 数据校验通过，未发现错误或警告</p>';
        }
        
        resultDetails.innerHTML = detailsHtml;
        validationResults.style.display = 'block';
        
    } catch (error) {
        console.error('显示校验结果失败:', error);
        showError(`显示校验结果失败: ${error.message}`);
    }
}

// 开始映射
async function startMapping() {
    try {
        updateStatus('正在进行基因映射和细胞类型标准化...');
        
        // 轮询检查状态
        const checkInterval = setInterval(async () => {
            try {
                const statusResponse = await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/status`);
                if (!statusResponse.ok) {
                    throw new Error('状态查询失败');
                }
                
                const statusData = await statusResponse.json();
                console.log('当前状态:', statusData);
                
                if (statusData.status === 'completed') {
                    clearInterval(checkInterval);
                    showMappingResults();
                    goToStep(4);
                    showFinalResults();
                } else if (statusData.status === 'failed') {
                    clearInterval(checkInterval);
                    showError(`处理失败: ${statusData.error_message || '未知错误'}`);
                }
                
            } catch (error) {
                clearInterval(checkInterval);
                console.error('状态检查失败:', error);
                showError(`状态检查失败: ${error.message}`);
            }
        }, 2000);
        
    } catch (error) {
        console.error('映射失败:', error);
        showError(`映射失败: ${error.message}`);
    }
}

// 显示映射结果
async function showMappingResults() {
    try {
        const reportResponse = await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/report`);
        if (!reportResponse.ok) {
            throw new Error('获取报告失败');
        }
        
        const reportData = await reportResponse.json();
        console.log('映射报告:', reportData);
        
        // 显示映射结果
        const mappingResults = document.getElementById('mappingResults');
        const mappingSummary = document.getElementById('mappingSummary');
        const mappingDetails = document.getElementById('mappingDetails');
        
        // 创建映射摘要
        const geneMapping = reportData.gene_mapping || {};
        const cellTypeMapping = reportData.cell_type_mapping || {};
        
        mappingSummary.innerHTML = `
            <div class="summary-card success">
                <h4>基因映射成功</h4>
                <div class="number">${geneMapping.mapped_count || 0}</div>
                <p>成功率: ${geneMapping.mapping_rate || 0}%</p>
            </div>
            <div class="summary-card success">
                <h4>细胞类型映射成功</h4>
                <div class="number">${cellTypeMapping.mapped_count || 0}</div>
                <p>成功率: ${cellTypeMapping.mapping_rate || 0}%</p>
            </div>
        `;
        
        // 显示映射详情
        const geneDetails = geneMapping.mapping_details || [];
        const cellTypeDetails = cellTypeMapping.mapping_details || [];
        
        let detailsHtml = '<h4>基因映射详情 (前10条):</h4>';
        if (geneDetails.length > 0) {
            detailsHtml += '<table class="details-table"><thead><tr><th>基因符号</th><th>Ensembl ID</th><th>状态</th><th>置信度</th></tr></thead><tbody>';
            geneDetails.slice(0, 10).forEach(detail => {
                const statusClass = detail.status === 'mapped' ? 'success' : 'warning';
                detailsHtml += `<tr>
                    <td>${detail.gene_symbol || ''}</td>
                    <td>${detail.ensembl_id || 'N/A'}</td>
                    <td class="${statusClass}">${detail.status || 'unknown'}</td>
                    <td>${detail.confidence || 'N/A'}</td>
                </tr>`;
            });
            detailsHtml += '</tbody></table>';
        } else {
            detailsHtml += '<p>暂无基因映射详情</p>';
        }
        
        detailsHtml += '<h4>细胞类型映射详情 (前10条):</h4>';
        if (cellTypeDetails.length > 0) {
            detailsHtml += '<table class="details-table"><thead><tr><th>原始类型</th><th>CL ID</th><th>CL标签</th><th>状态</th></tr></thead><tbody>';
            cellTypeDetails.slice(0, 10).forEach(detail => {
                const statusClass = detail.status === 'mapped' ? 'success' : 'warning';
                detailsHtml += `<tr>
                    <td>${detail.raw_cell_type || ''}</td>
                    <td>${detail.cl_id || 'N/A'}</td>
                    <td>${detail.cl_label || 'N/A'}</td>
                    <td class="${statusClass}">${detail.status || 'unknown'}</td>
                </tr>`;
            });
            detailsHtml += '</tbody></table>';
        } else {
            detailsHtml += '<p>暂无细胞类型映射详情</p>';
        }
        
        mappingDetails.innerHTML = detailsHtml;
        mappingResults.style.display = 'block';
        
    } catch (error) {
        console.error('显示映射结果失败:', error);
        showError(`显示映射结果失败: ${error.message}`);
    }
}

// 显示最终结果
async function showFinalResults() {
    try {
        updateStatus('正在生成最终结果...');
        
        const reportResponse = await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/report`);
        if (!reportResponse.ok) {
            throw new Error('获取最终报告失败');
        }
        
        const reportData = await reportResponse.json();
        console.log('最终报告:', reportData);
        
        // 显示最终结果
        const finalResults = document.getElementById('finalResults');
        const finalSummary = document.getElementById('finalSummary');
        
        const summary = reportData.summary || {};
        const qualityMetrics = reportData.quality_metrics || {};
        
        finalSummary.innerHTML = `
            <div class="summary-card success">
                <h4>处理完成</h4>
                <div class="number">✓</div>
                <p>数据已成功处理</p>
            </div>
            <div class="summary-card success">
                <h4>细胞数</h4>
                <div class="number">${summary.n_cells || 0}</div>
            </div>
            <div class="summary-card success">
                <h4>基因数</h4>
                <div class="number">${summary.n_genes || 0}</div>
            </div>
            <div class="summary-card success">
                <h4>稀疏度</h4>
                <div class="number">${(qualityMetrics.sparsity * 100 || 0).toFixed(1)}%</div>
            </div>
        `;
        
        finalResults.style.display = 'block';
        updateStatus('处理完成');
        resetBtn.style.display = 'inline-block';
        
        // 设置下载按钮
        document.getElementById('downloadReport').onclick = () => downloadReport();
        document.getElementById('downloadData').onclick = () => downloadData();
        
    } catch (error) {
        console.error('显示最终结果失败:', error);
        showError(`显示最终结果失败: ${error.message}`);
    }
}

// 下载报告
async function downloadReport() {
    try {
        const response = await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/report`);
        if (!response.ok) {
            throw new Error('下载报告失败');
        }
        
        const blob = await response.blob();
        const url = window.URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `report_${currentSubmissionId}.json`;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        window.URL.revokeObjectURL(url);
        
    } catch (error) {
        console.error('下载报告失败:', error);
        showError(`下载报告失败: ${error.message}`);
    }
}

// 下载数据
async function downloadData() {
    try {
        const response = await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/export`);
        if (!response.ok) {
            throw new Error('下载数据失败');
        }
        
        const blob = await response.blob();
        const url = window.URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `processed_${currentSubmissionId}.h5ad`;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        window.URL.revokeObjectURL(url);
        
    } catch (error) {
        console.error('下载数据失败:', error);
        showError(`下载数据失败: ${error.message}`);
    }
}

// 步骤切换
function goToStep(step) {
    // 更新步骤指示器
    stepElements.forEach((element, index) => {
        element.classList.remove('active', 'completed');
        if (index + 1 < step) {
            element.classList.add('completed');
        } else if (index + 1 === step) {
            element.classList.add('active');
        }
    });
    
    // 更新内容区域
    contentElements.forEach((element, index) => {
        element.classList.remove('active');
        if (index + 1 === step) {
            element.classList.add('active');
        }
    });
    
    currentStep = step;
}

// 更新状态
function updateStatus(message) {
    document.getElementById('currentStatus').textContent = message;
}

// 显示错误
function showError(message) {
    errorMessage.textContent = message;
    errorModal.style.display = 'flex';
}

// 关闭错误模态框
function closeErrorModal() {
    errorModal.style.display = 'none';
}

// 重置流程
function resetProcess() {
    currentSubmissionId = null;
    currentStep = 1;
    
    // 重置UI
    goToStep(1);
    fileInfo.style.display = 'none';
    uploadBtn.disabled = true;
    uploadBtn.file = null;
    fileInput.value = '';
    
    // 隐藏结果区域
    document.getElementById('validationResults').style.display = 'none';
    document.getElementById('mappingResults').style.display = 'none';
    document.getElementById('finalResults').style.display = 'none';
    resetBtn.style.display = 'none';
    
    updateStatus('准备就绪');
}

// 工具函数
function formatFileSize(bytes) {
    if (bytes === 0) return '0 Bytes';
    const k = 1024;
    const sizes = ['Bytes', 'KB', 'MB', 'GB'];
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i];
}
