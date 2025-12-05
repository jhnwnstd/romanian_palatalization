#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Full pipeline for Romanian palatalization analysis
#
# Pipeline stages:
#   1. [Optional] Harvest raw data from Wiktionary (hours)
#   2. [Optional] Download Leipzig corpus frequency data
#   3. [Optional] Run DEX QC on harvested data
#   4. Process data and derive palatalization features
#   5. Attach frequency counts to final lexicon
#   6. Run R statistical analysis
#
# Usage:
# ./run_pipeline.sh             # skips data collection stages if data exists
# ./run_pipeline.sh --force     # re-run all stages 
#

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Project paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data"
SCRIPTS_DIR="$SCRIPT_DIR/scripts"
ANALYSIS_DIR="$SCRIPT_DIR/analysis"

# Data files
RAW_CSV="$DATA_DIR/romanian_lexicon_raw.csv"
DEX_QC_CSV="$DATA_DIR/romanian_lexicon_raw_dex.csv"
COMPLETE_CSV="$DATA_DIR/romanian_lexicon_complete.csv"
FINAL_CSV="$DATA_DIR/romanian_lexicon_with_freq.csv"

# Audit files
AUDIT_CSV="$DATA_DIR/dex_qc_audit.csv"
DISAGREEMENTS_CSV="$DATA_DIR/dex_qc_disagreements.csv"

# Analysis files
R_SCRIPT="$ANALYSIS_DIR/analyze_romanian_palatalization.R"
ANALYSIS_LOG="$ANALYSIS_DIR/romanian_palatalization_analysis.txt"

# Parse arguments
FORCE_RUN=0
if [[ "${1:-}" == "--force" ]]; then
    FORCE_RUN=1
fi

# Helper functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[✓]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[!]${NC} $1"
}

log_error() {
    echo -e "${RED}[✗]${NC} $1"
}

log_stage() {
    echo ""
    echo -e "${GREEN}═══════════════════════════════════════════════════${NC}"
    echo -e "${GREEN}  $1${NC}"
    echo -e "${GREEN}═══════════════════════════════════════════════════${NC}"
}

# Check if file exists and show size
check_file() {
    local file="$1"
    if [[ -f "$file" ]]; then
        local size=$(du -h "$file" | cut -f1)
        log_success "Found: $(basename "$file") ($size)"
        return 0
    else
        log_warning "Not found: $(basename "$file")"
        return 1
    fi
}

# Main pipeline
main() {
    log_info "Romanian Palatalization Analysis Pipeline"
    log_info "Working directory: $SCRIPT_DIR"
    echo ""

    # Create data directory if needed
    mkdir -p "$DATA_DIR"

    # ========================================
    # STAGE 1: Harvest raw data (optional)
    # ========================================
    log_stage "STAGE 1: Harvest Raw Data from Wiktionary"

    if [[ -f "$RAW_CSV" ]]; then
        if [[ $FORCE_RUN -eq 1 ]]; then
            log_warning "Force mode: re-running harvest (existing data will be overwritten)"
            log_info "Starting harvest (this may take several hours)..."
            python3 "$SCRIPTS_DIR/romanian_harvester.py"

            if check_file "$RAW_CSV"; then
                log_success "Harvest complete!"
            else
                log_error "Harvest failed - raw CSV not created"
                exit 1
            fi
        else
            log_info "Raw data already exists"
            check_file "$RAW_CSV"
            log_info "Skipping harvest (use --force to re-run)"
        fi
    else
        log_info "Raw data not found - starting harvest (this may take several hours)..."
        python3 "$SCRIPTS_DIR/romanian_harvester.py"

        if check_file "$RAW_CSV"; then
            log_success "Harvest complete!"
        else
            log_error "Harvest failed - raw CSV not created"
            exit 1
        fi
    fi

    # ========================================
    # STAGE 2: Download Leipzig corpora
    # ========================================
    log_stage "STAGE 2: Download Leipzig Frequency Data"

    FREQ_DIR="$DATA_DIR/leipzig/freq"
    if [[ -d "$FREQ_DIR" ]] && [[ -n "$(ls -A "$FREQ_DIR"/*.csv 2>/dev/null)" ]]; then
        if [[ $FORCE_RUN -eq 1 ]]; then
            log_warning "Force mode: re-downloading Leipzig corpora"
            python3 "$SCRIPTS_DIR/download_leipzig.py"

            if [[ -d "$FREQ_DIR" ]] && [[ -n "$(ls -A "$FREQ_DIR"/*.csv 2>/dev/null)" ]]; then
                log_success "Leipzig frequency data ready"
            else
                log_error "Leipzig download failed - no frequency files found"
                exit 1
            fi
        else
            log_success "Leipzig frequency data already exists"
            log_info "Found frequency files:"
            ls -lh "$FREQ_DIR"/*.csv | awk '{print "  - " $9 " (" $5 ")"}'
            log_info "Skipping download (use --force to re-run)"
        fi
    else
        log_info "Downloading and processing Leipzig corpora..."
        python3 "$SCRIPTS_DIR/download_leipzig.py"

        if [[ -d "$FREQ_DIR" ]] && [[ -n "$(ls -A "$FREQ_DIR"/*.csv 2>/dev/null)" ]]; then
            log_success "Leipzig frequency data ready"
            log_info "Found frequency files:"
            ls -lh "$FREQ_DIR"/*.csv | awk '{print "  - " $9 " (" $5 ")"}'
        else
            log_error "Leipzig download failed - no frequency files found"
            exit 1
        fi
    fi

    # ========================================
    # STAGE 3: DEX QC
    # ========================================
    log_stage "STAGE 3: DEX Quality Control"

    if [[ -f "$DEX_QC_CSV" ]]; then
        if [[ $FORCE_RUN -eq 1 ]]; then
            # Check if we have the input (RAW_CSV)
            if [[ ! -f "$RAW_CSV" ]]; then
                log_error "Cannot run DEX QC: raw data not found ($RAW_CSV)"
                log_error "Please run harvest first"
                exit 1
            fi

            log_warning "Force mode: re-running DEX QC"
            python3 "$SCRIPTS_DIR/dex_qc_main_csv.py"

            if check_file "$DEX_QC_CSV"; then
                log_success "DEX QC complete!"
                if check_file "$AUDIT_CSV"; then
                    log_info "Audit report: $AUDIT_CSV"
                fi
                if check_file "$DISAGREEMENTS_CSV"; then
                    log_info "Disagreements: $DISAGREEMENTS_CSV"
                fi
            else
                log_error "DEX QC failed - output not created"
                exit 1
            fi
        else
            log_info "DEX QC data already exists"
            check_file "$DEX_QC_CSV"
            if check_file "$AUDIT_CSV"; then
                log_info "Audit report: $AUDIT_CSV"
            fi
            if check_file "$DISAGREEMENTS_CSV"; then
                log_info "Disagreements: $DISAGREEMENTS_CSV"
            fi
            log_info "Skipping DEX QC (use --force to re-run)"
        fi
    else
        # DEX QC needs to run - check if we have the input (RAW_CSV)
        if [[ ! -f "$RAW_CSV" ]]; then
            log_error "Cannot run DEX QC: raw data not found ($RAW_CSV)"
            log_error "Please run harvest first or provide raw data file"
            exit 1
        fi

        log_info "Running DEX QC to fill missing data..."
        python3 "$SCRIPTS_DIR/dex_qc_main_csv.py"

        if check_file "$DEX_QC_CSV"; then
            log_success "DEX QC complete!"
            if check_file "$AUDIT_CSV"; then
                log_info "Audit report: $AUDIT_CSV"
            fi
            if check_file "$DISAGREEMENTS_CSV"; then
                log_info "Disagreements: $DISAGREEMENTS_CSV"
            fi
        else
            log_error "DEX QC failed - output not created"
            exit 1
        fi
    fi

    # ========================================
    # STAGE 4: Derive palatalization features
    # ========================================
    log_stage "STAGE 4: Derive Palatalization Features"

    log_info "Processing data and deriving features..."
    python3 "$SCRIPTS_DIR/romanian_processor_main.py"

    if check_file "$COMPLETE_CSV"; then
        log_success "Feature derivation complete!"
    else
        log_error "Feature derivation failed - output not created"
        exit 1
    fi

    # ========================================
    # STAGE 5: Attach frequencies
    # ========================================
    log_stage "STAGE 5: Attach Frequency Counts"

    log_info "Attaching Leipzig frequency data to lexicon..."
    python3 "$SCRIPTS_DIR/attach_frequencies.py"

    if check_file "$FINAL_CSV"; then
        log_success "Frequency attachment complete!"
    else
        log_error "Frequency attachment failed - output not created"
        exit 1
    fi

    # ========================================
    # STAGE 6: R Statistical Analysis
    # ========================================
    log_stage "STAGE 6: Run R Statistical Analysis"

    if [[ ! -f "$R_SCRIPT" ]]; then
        log_warning "R analysis script not found: $R_SCRIPT"
        log_info "Skipping analysis stage"
    else
        log_info "Running R statistical analysis..."
        log_info "This may take a few minutes..."
        log_info "Output will be written to: $ANALYSIS_LOG"
        echo ""

        if Rscript "$R_SCRIPT" > /dev/null 2>&1; then
            log_success "R analysis complete!"
            log_info "Results saved to: $ANALYSIS_LOG"
        else
            log_error "R analysis failed"
            log_error "Re-running with output for debugging..."
            Rscript "$R_SCRIPT"
            exit 1
        fi
    fi

    # ========================================
    # PIPELINE COMPLETE
    # ========================================
    log_stage "PIPELINE COMPLETE!"

    echo ""
    log_info "Output files:"
    echo ""
    check_file "$RAW_CSV"
    check_file "$DEX_QC_CSV"
    check_file "$COMPLETE_CSV"
    check_file "$FINAL_CSV"
    echo ""
    log_success "All stages completed successfully!"
    log_info "Final lexicon: $FINAL_CSV"
    echo ""
}

# Run main pipeline
main "$@"
