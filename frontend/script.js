// Global variables
let currentSubmissionId = null;
let currentStep = 1;
const API_BASE_URL = 'http://localhost:8000';

// DOM elements
const stepElements = document.querySelectorAll('.step');
const contentElements = document.querySelectorAll('.step-content');
const fileInput = document.getElementById('fileInput');
const uploadArea = document.getElementById('uploadArea');
const fileInfo = document.getElementById('fileInfo');
const uploadBtn = document.getElementById('uploadBtn');
const resetBtn = document.getElementById('resetBtn');
const errorModal = document.getElementById('errorModal');
const errorMessage = document.getElementById('errorMessage');

// Initialization
document.addEventListener('DOMContentLoaded', function() {
    initializeEventListeners();
    updateStatus('Ready');
});

// Event listeners
function initializeEventListeners() {
    // File upload
    fileInput.addEventListener('change', handleFileSelect);
    uploadArea.addEventListener('click', () => fileInput.click());
    uploadArea.addEventListener('dragover', handleDragOver);
    uploadArea.addEventListener('dragleave', handleDragLeave);
    uploadArea.addEventListener('drop', handleDrop);
    uploadBtn.addEventListener('click', handleUpload);
    resetBtn.addEventListener('click', resetProcess);
    
    // Confirmation buttons
    document.getElementById('continueToMapping').addEventListener('click', startMapping);
    document.getElementById('cancelProcess').addEventListener('click', cancelProcess);
    document.getElementById('continueToExport').addEventListener('click', startExport);
    document.getElementById('cancelProcess2').addEventListener('click', cancelProcess);
}

// File selection handler
function handleFileSelect(event) {
    const file = event.target.files[0];
    if (file) {
        processFile(file);
    }
}

// Drag and drop handlers
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

// Process file
function processFile(file) {
    if (!file.name.endsWith('.h5ad')) {
        showError('Please select a .h5ad format file');
        return;
    }
    
    // Display file information
    document.getElementById('fileName').textContent = file.name;
    document.getElementById('fileSize').textContent = formatFileSize(file.size);
    fileInfo.style.display = 'block';
    uploadBtn.disabled = false;
    
    // Store file reference
    uploadBtn.file = file;
}

// File upload
async function handleUpload() {
    const file = uploadBtn.file;
    if (!file) return;
    
    try {
        updateStatus('Creating submission task...');
        
        // 1. Create submission task
        const submissionResponse = await fetch(`${API_BASE_URL}/submissions`, {
            method: 'POST'
        });
        
        if (!submissionResponse.ok) {
            throw new Error('Failed to create submission task');
        }
        
        const submissionData = await submissionResponse.json();
        currentSubmissionId = submissionData.submission_id;
        
        updateStatus('Uploading file...');
        
        // 2. Upload file
        const formData = new FormData();
        formData.append('file', file);
        
        const uploadResponse = await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/upload`, {
            method: 'POST',
            body: formData
        });
        
        if (!uploadResponse.ok) {
            const errorText = await uploadResponse.text();
            throw new Error(`File upload failed: ${errorText}`);
        }
        
        const uploadData = await uploadResponse.json();
        console.log('Upload successful:', uploadData);
        
        // 3. Proceed to validation step
        goToStep(2);
        startValidation();
        
    } catch (error) {
        console.error('Upload failed:', error);
        showError(`Upload failed: ${error.message}`);
    }
}

// Start validation
async function startValidation() {
    try {
        updateStatus('Validating data...');
        
        // Call validation API
        const validateResponse = await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/validate`, {
            method: 'POST'
        });
        
        if (!validateResponse.ok) {
            throw new Error('Failed to start validation');
        }
        
        // Poll for status
        const checkInterval = setInterval(async () => {
            try {
                const statusResponse = await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/status`);
                if (!statusResponse.ok) {
                    throw new Error('Failed to query status');
                }
                
                const statusData = await statusResponse.json();
                console.log('Current status:', statusData);
                
                if (statusData.status === 'validated') {
                    clearInterval(checkInterval);
                    await showValidationResults();
                    document.getElementById('validationConfirmation').style.display = 'block';
                    updateStatus('Validation completed, waiting for confirmation');
                } else if (statusData.status === 'validation_failed' || statusData.status === 'failed') {
                    clearInterval(checkInterval);
                    showError(`Validation failed: ${statusData.error_message || 'Unknown error'}`);
                }
                // Continue waiting for validating status
                
            } catch (error) {
                clearInterval(checkInterval);
                console.error('Status check failed:', error);
                showError(`Status check failed: ${error.message}`);
            }
        }, 2000); // Check every 2 seconds
        
    } catch (error) {
        console.error('Validation failed:', error);
        showError(`Validation failed: ${error.message}`);
    }
}

// Show validation results
async function showValidationResults() {
    try {
        const reportResponse = await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/report`);
        if (!reportResponse.ok) {
            throw new Error('Failed to fetch report');
        }
        
        const reportData = await reportResponse.json();
        console.log('Validation report:', reportData);
        
        // Display validation results
        const validationResults = document.getElementById('validationResults');
        const resultSummary = document.getElementById('resultSummary');
        const resultDetails = document.getElementById('resultDetails');
        
        // Hide loading status
        document.getElementById('validationStatus').style.display = 'none';
        
        // Create summary cards
        const summary = reportData.summary || {};
        const validation = reportData.validation_result || {};
        
        resultSummary.innerHTML = `
            <div class="summary-card success">
                <h4>Number of Cells</h4>
                <div class="number">${summary.n_cells || 0}</div>
            </div>
            <div class="summary-card success">
                <h4>Number of Genes</h4>
                <div class="number">${summary.n_genes || 0}</div>
            </div>
            <div class="summary-card success">
                <h4>File Size</h4>
                <div class="number">${summary.file_size_mb || 0} MB</div>
            </div>
        `;
        
        // Display validation details
        const errors = validation.errors || [];
        const warnings = validation.warnings || [];
        
        let detailsHtml = '';
        if (errors.length > 0) {
            detailsHtml += '<h4 style="color: #dc2626; margin-top: 20px;">Errors:</h4><ul style="color: #dc2626; margin-left: 20px;">';
            errors.forEach(error => {
                detailsHtml += `<li style="margin: 8px 0;">${error}</li>`;
            });
            detailsHtml += '</ul>';
        }
        
        if (warnings.length > 0) {
            detailsHtml += '<h4 style="color: #d97706; margin-top: 20px;">Warnings:</h4><ul style="color: #d97706; margin-left: 20px;">';
            warnings.forEach(warning => {
                detailsHtml += `<li style="margin: 8px 0;">${warning}</li>`;
            });
            detailsHtml += '</ul>';
        }
        
        if (errors.length === 0 && warnings.length === 0) {
            detailsHtml = '<p style="color: #0d4f8c; font-size: 1.1em; margin-top: 20px;">✓ Validation passed, no errors or warnings found</p>';
        }
        
        resultDetails.innerHTML = detailsHtml;
        validationResults.style.display = 'block';
        
    } catch (error) {
        console.error('Failed to show validation results:', error);
        showError(`Failed to show validation results: ${error.message}`);
    }
}

// Start mapping
async function startMapping() {
    try {
        goToStep(3);
        updateStatus('Performing gene mapping and cell type standardization...');
        
        // Call continue processing API
        const continueResponse = await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/continue`, {
            method: 'POST'
        });
        
        if (!continueResponse.ok) {
            throw new Error('Failed to start mapping');
        }
        
        // Poll for status
        const checkInterval = setInterval(async () => {
            try {
                const statusResponse = await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/status`);
                if (!statusResponse.ok) {
                    throw new Error('Failed to query status');
                }
                
                const statusData = await statusResponse.json();
                console.log('Current status:', statusData);
                
                if (statusData.status === 'completed') {
                    clearInterval(checkInterval);
                    await showMappingResults();
                    document.getElementById('mappingConfirmation').style.display = 'block';
                    updateStatus('Mapping completed, waiting for confirmation');
                } else if (statusData.status === 'failed') {
                    clearInterval(checkInterval);
                    showError(`Processing failed: ${statusData.error_message || 'Unknown error'}`);
                }
                
            } catch (error) {
                clearInterval(checkInterval);
                console.error('Status check failed:', error);
                showError(`Status check failed: ${error.message}`);
            }
        }, 2000);
        
    } catch (error) {
        console.error('Mapping failed:', error);
        showError(`Mapping failed: ${error.message}`);
    }
}

// Show mapping results
async function showMappingResults() {
    try {
        const reportResponse = await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/report`);
        if (!reportResponse.ok) {
            throw new Error('Failed to fetch report');
        }
        
        const reportData = await reportResponse.json();
        console.log('Mapping report:', reportData);
        
        // Display mapping results
        const mappingResults = document.getElementById('mappingResults');
        const mappingSummary = document.getElementById('mappingSummary');
        const mappingDetails = document.getElementById('mappingDetails');
        
        // Hide loading status
        document.getElementById('mappingStatus').style.display = 'none';
        
        // Create mapping summary
        const geneMapping = reportData.gene_mapping || {};
        const cellTypeMapping = reportData.cell_type_mapping || {};
        
        mappingSummary.innerHTML = `
            <div class="summary-card success">
                <h4>Gene Mapping Success</h4>
                <div class="number">${geneMapping.mapped_count || 0}</div>
                <p>Success Rate: ${geneMapping.mapping_rate || 0}%</p>
            </div>
            <div class="summary-card success">
                <h4>Cell Type Mapping Success</h4>
                <div class="number">${cellTypeMapping.mapped_count || 0}</div>
                <p>Success Rate: ${cellTypeMapping.mapping_rate || 0}%</p>
            </div>
        `;
        
        // Display mapping details
        const geneDetails = geneMapping.mapping_details || [];
        const cellTypeDetails = cellTypeMapping.mapping_details || [];
        
        let detailsHtml = '<h4>Gene Mapping Details (Top 10):</h4>';
        if (geneDetails.length > 0) {
            detailsHtml += '<table class="details-table"><thead><tr><th>Gene Symbol</th><th>Ensembl ID</th><th>Status</th><th>Confidence</th></tr></thead><tbody>';
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
            detailsHtml += '<p>No gene mapping details available</p>';
        }
        
        detailsHtml += '<h4 style="margin-top: 30px;">Cell Type Mapping Details (Top 10):</h4>';
        if (cellTypeDetails.length > 0) {
            detailsHtml += '<table class="details-table"><thead><tr><th>Raw Type</th><th>CL ID</th><th>CL Label</th><th>Status</th></tr></thead><tbody>';
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
            detailsHtml += '<p>No cell type mapping details available</p>';
        }
        
        mappingDetails.innerHTML = detailsHtml;
        mappingResults.style.display = 'block';
        
    } catch (error) {
        console.error('Failed to show mapping results:', error);
        showError(`Failed to show mapping results: ${error.message}`);
    }
}

// Start export
async function startExport() {
    goToStep(4);
    showFinalResults();
}

// Show final results
async function showFinalResults() {
    try {
        updateStatus('Generating final results...');
        
        const reportResponse = await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/report`);
        if (!reportResponse.ok) {
            throw new Error('Failed to fetch final report');
        }
        
        const reportData = await reportResponse.json();
        console.log('Final report:', reportData);
        
        // Display final results
        const finalResults = document.getElementById('finalResults');
        const finalSummary = document.getElementById('finalSummary');
        const resultsContainer = document.getElementById('resultsContainer');
        
        resultsContainer.style.display = 'none';
        
        const summary = reportData.summary || {};
        const qualityMetrics = reportData.quality_metrics || {};
        
        finalSummary.innerHTML = `
            <div class="summary-card success">
                <h4>Processing Completed</h4>
                <div class="number">✓</div>
                <p>Data successfully processed</p>
            </div>
            <div class="summary-card success">
                <h4>Number of Cells</h4>
                <div class="number">${summary.n_cells || 0}</div>
            </div>
            <div class="summary-card success">
                <h4>Number of Genes</h4>
                <div class="number">${summary.n_genes || 0}</div>
            </div>
            <div class="summary-card success">
                <h4>Sparsity</h4>
                <div class="number">${(qualityMetrics.sparsity * 100 || 0).toFixed(1)}%</div>
            </div>
        `;
        
        finalResults.style.display = 'block';
        updateStatus('Processing completed');
        resetBtn.style.display = 'inline-block';
        
        // Set download buttons
        document.getElementById('downloadReport').onclick = () => downloadReport();
        document.getElementById('downloadData').onclick = () => downloadData();
        
    } catch (error) {
        console.error('Failed to show final results:', error);
        showError(`Failed to show final results: ${error.message}`);
    }
}

// Cancel processing
function cancelProcess() {
    if (confirm('Are you sure you want to cancel the processing?')) {
        resetProcess();
    }
}

// Download report
async function downloadReport() {
    try {
        const response = await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/report`);
        if (!response.ok) {
            throw new Error('Failed to download report');
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
        console.error('Failed to download report:', error);
        showError(`Failed to download report: ${error.message}`);
    }
}

// Download data
async function downloadData() {
    try {
        const response = await fetch(`${API_BASE_URL}/submissions/${currentSubmissionId}/export`);
        if (!response.ok) {
            throw new Error('Failed to download data');
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
        console.error('Failed to download data:', error);
        showError(`Failed to download data: ${error.message}`);
    }
}

// Step navigation
function goToStep(step) {
    // Update step indicator
    stepElements.forEach((element, index) => {
        element.classList.remove('active', 'completed');
        if (index + 1 < step) {
            element.classList.add('completed');
        } else if (index + 1 === step) {
            element.classList.add('active');
        }
    });
    
    // Update content area
    contentElements.forEach((element, index) => {
        element.classList.remove('active');
        if (index + 1 === step) {
            element.classList.add('active');
        }
    });
    
    currentStep = step;
}

// Update status
function updateStatus(message) {
    document.getElementById('currentStatus').textContent = message;
}

// Show error
function showError(message) {
    errorMessage.textContent = message;
    errorModal.style.display = 'flex';
}

// Close error modal
function closeErrorModal() {
    errorModal.style.display = 'none';
}

// Reset process
function resetProcess() {
    currentSubmissionId = null;
    currentStep = 1;
    
    // Reset UI
    goToStep(1);
    fileInfo.style.display = 'none';
    uploadBtn.disabled = true;
    uploadBtn.file = null;
    fileInput.value = '';
    
    // Hide result areas
    document.getElementById('validationResults').style.display = 'none';
    document.getElementById('validationConfirmation').style.display = 'none';
    document.getElementById('validationStatus').style.display = 'block';
    document.getElementById('mappingResults').style.display = 'none';
    document.getElementById('mappingConfirmation').style.display = 'none';
    document.getElementById('mappingStatus').style.display = 'block';
    document.getElementById('finalResults').style.display = 'none';
    document.getElementById('resultsContainer').style.display = 'block';
    resetBtn.style.display = 'none';
    
    updateStatus('Ready');
}

// Utility function
function formatFileSize(bytes) {
    if (bytes === 0) return '0 Bytes';
    const k = 1024;
    const sizes = ['Bytes', 'KB', 'MB', 'GB'];
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i];
}

// Global function for HTML calls
window.closeErrorModal = closeErrorModal;