#!/bin/bash

#################################################################################
# DG-Fragment RT-TDDFT: Proper Two-Step Workflow
#
# 正しい実行フロー:
#   1. Ground State (GS) 計算  → Fragment basis 生成
#   2. Real-Time (RT) 計算    → GS output を読込、時間発展
#
# このスクリプトはそれを自動化します
#
#################################################################################

set -e

SALMON_DIR="/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2"
MYBUILD_DIR="${SALMON_DIR}/mybuild"
SALMON_EXE="${MYBUILD_DIR}/salmon"
TESTDIR="$(cd "$(dirname "$0")" && pwd)/H2"
LOGDIR="${TESTDIR}/logs"
PSEUDO_DIR="${SALMON_DIR}/pseudo"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Check if salmon executable exists
if [ ! -f "$SALMON_EXE" ]; then
    echo -e "${RED}ERROR: SALMON executable not found at $SALMON_EXE${NC}"
    exit 1
fi

# Function: Run GS calculation
run_gs() {
    echo -e "${BLUE}┌─────────────────────────────────────────────────────┐${NC}"
    echo -e "${BLUE}│ STEP 1: Ground State (GS) - LDA Basis Generation    │${NC}"
    echo -e "${BLUE}└─────────────────────────────────────────────────────┘${NC}"
    
    cd "${TESTDIR}"
    mkdir -p "$LOGDIR"
    
    echo -e "${YELLOW}Running: $SALMON_EXE < inputfile_lda_gs${NC}"
    local start=$(date +%s)
    
    $SALMON_EXE < inputfile_lda_gs > stdout_lda_gs.log 2>&1
    local exit_code=$?
    
    local end=$(date +%s)
    local elapsed=$((end - start))
    
    if [ $exit_code -eq 0 ]; then
        echo -e "${GREEN}✓ GS calculation completed${NC}"
        echo "  Elapsed time: ${elapsed} seconds"
        
        # Check if basis files exist
        if ls nutilus_frag*.dat 1>/dev/null 2>&1; then
            echo -e "${GREEN}✓ Fragment basis files generated${NC}"
            ls -lh nutilus_frag*.dat | head -3
        else
            echo -e "${RED}✗ Fragment basis files NOT found${NC}"
            return 1
        fi
    else
        echo -e "${RED}✗ GS calculation failed${NC}"
        tail -20 stdout_lda_gs.log
        return 1
    fi
    
    cp stdout_lda_gs.log "$LOGDIR/"
    echo ""
}

# Function: Run RT calculation
run_rt() {
    local model=$1
    local inputfile=$2
    
    echo -e "${BLUE}┌─────────────────────────────────────────────────────┐${NC}"
    echo -e "${BLUE}│ STEP 2: Real-Time (RT) - ${model}                    │${NC}"
    echo -e "${BLUE}└─────────────────────────────────────────────────────┘${NC}"
    
    cd "${TESTDIR}"
    
    # Check if GS basis exists
    if ! ls nutilus_frag*.dat 1>/dev/null 2>&1; then
        echo -e "${RED}ERROR: Fragment basis files not found${NC}"
        echo "Please run GS calculation first: ./run_h2_test_proper.sh gs"
        return 1
    fi
    
    echo -e "${YELLOW}Running: $SALMON_EXE < ${inputfile}${NC}"
    local start=$(date +%s)
    
    $SALMON_EXE < "$inputfile" > "stdout_${model}_rt.log" 2>&1
    local exit_code=$?
    
    local end=$(date +%s)
    local elapsed=$((end - start))
    
    if [ $exit_code -eq 0 ]; then
        echo -e "${GREEN}✓ RT calculation completed${NC}"
        echo "  Elapsed time: ${elapsed} seconds"
        
        # Extract energy
        if grep -q "Total energy" "stdout_${model}_rt.log"; then
            echo -e "${GREEN}✓ Energy output found:${NC}"
            grep "Total energy" "stdout_${model}_rt.log" | tail -1
        fi
    else
        echo -e "${RED}✗ RT calculation failed${NC}"
        tail -20 "stdout_${model}_rt.log"
        return 1
    fi
    
    cp "stdout_${model}_rt.log" "$LOGDIR/"
    echo ""
}

# Function: Clean
clean() {
    cd "${TESTDIR}"
    echo -e "${YELLOW}Cleaning...${NC}"
    rm -rf work_* logs/ nutilus_frag*.dat Info_frag*.txt Eig_frag*.txt ene.out
    echo -e "${GREEN}Clean complete${NC}"
}

# Main
main() {
    case "${1:-full}" in
        gs)
            run_gs
            echo -e "${GREEN}═════════════════════════════════════════════════════${NC}"
            echo -e "${GREEN}GS calculation complete${NC}"
            echo -e "${YELLOW}Next: Run RT calculations with ./run_h2_test_proper.sh rt${NC}"
            ;;
        rt)
            echo -e "${YELLOW}Running all RT calculations...${NC}"
            mkdir -p "$LOGDIR"
            run_rt "lda"    "inputfile_lda_rt"
            run_rt "plan_a" "inputfile_plan_a_rt"
            run_rt "plan_c" "inputfile_plan_c"
            run_rt "cdri"   "inputfile_cdri"
            
            echo -e "${GREEN}═════════════════════════════════════════════════════${NC}"
            echo -e "${GREEN}All RT calculations complete${NC}"
            echo -e "${YELLOW}Compare results:${NC}"
            grep "Total energy" "$LOGDIR/stdout_*_rt.log"
            ;;
        full)
            echo -e "${BLUE}═════════════════════════════════════════════════════${NC}"
            echo -e "${BLUE}  H2 DG-Fragment RT-TDDFT Complete Test${NC}"
            echo -e "${BLUE}  (GS → LDA RT → HSE RT)${NC}"
            echo -e "${BLUE}═════════════════════════════════════════════════════${NC}"
            echo ""
            
            run_gs
            run_rt "lda"    "inputfile_lda_rt"
            run_rt "plan_a" "inputfile_plan_a_rt"
            run_rt "plan_c" "inputfile_plan_c_rt"
            run_rt "cdri"   "inputfile_cdri_rt"
            
            echo -e "${BLUE}═════════════════════════════════════════════════════${NC}"
            echo -e "${BLUE}  Test Complete!${NC}"
            echo -e "${BLUE}═════════════════════════════════════════════════════${NC}"
            echo -e "${YELLOW}Energy summary:${NC}"
            grep "Total energy" "$LOGDIR/stdout_*_rt.log"
            ;;
        clean)
            clean
            ;;
        *)
            echo "Usage: $0 [gs|rt|full|clean]"
            echo "  gs    - Run GS calculation only"
            echo "  rt    - Run all RT calculations (requires GS output)"
            echo "  full  - Run everything (GS + all RTs)"
            echo "  clean - Remove all output files"
            ;;
    esac
}

main "$@"
