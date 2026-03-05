#!/bin/bash

#################################################################################
# C6H6 Test Suite for DG RT-TDDFT HSE Implementation
# 
# Purpose: Test Plan C (RI/DF) and CD-RI on medium-scale system
# System: C6H6 benzene molecule (12 atoms, ~42 basis functions)
# 
# Usage:
#   ./run_c6h6_test.sh
#   ./run_c6h6_test.sh quick    # Quick test with reduced parameters
#   ./run_c6h6_test.sh clean    # Remove all output files
#
#################################################################################

set -e

SALMON_DIR="/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2"
MYBUILD_DIR="${SALMON_DIR}/mybuild"
SALMON_EXE="${MYBUILD_DIR}/src/salmon"
TESTDIR="$(cd "$(dirname "$0")" && pwd)/C6H6"
LOGDIR="${TESTDIR}/logs"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Check if salmon executable exists
if [ ! -f "$SALMON_EXE" ]; then
    echo -e "${RED}ERROR: SALMON executable not found at $SALMON_EXE${NC}"
    echo "Please build SALMON first:"
    echo "  cd $MYBUILD_DIR"
    echo "  make -j8"
    exit 1
fi

# Function: Run a single test
run_test() {
    local plan=$1
    local inputfile=$2
    local name="${plan}"
    
    echo -e "${BLUE}┌─────────────────────────────────────────────────────┐${NC}"
    echo -e "${BLUE}│ Testing Plan ${name} (C6H6)${NC}"
    echo -e "${BLUE}└─────────────────────────────────────────────────────┘${NC}"
    
    # Create working directory
    local workdir="${TESTDIR}/work_${plan}"
    mkdir -p "$workdir"
    cd "$workdir"
    
    # Copy input file
    cp "${TESTDIR}/${inputfile}" input.inp
    
    # Prepare symbolic links for pseudo potentials
    if [ ! -f "C_rps.dat" ]; then
        ln -sf "${SALMON_DIR}/samples/C_rps.dat" "C_rps.dat"
    fi
    if [ ! -f "H_rps.dat" ]; then
        ln -sf "${SALMON_DIR}/samples/H_rps.dat" "H_rps.dat"
    fi
    
    # Run SALMON
    echo -e "${YELLOW}Running: $SALMON_EXE < input.inp${NC}"
    local start_time=$(date +%s)
    
    $SALMON_EXE < input.inp > stdout.log 2>&1
    local exit_code=$?
    
    local end_time=$(date +%s)
    local elapsed=$((end_time - start_time))
    
    # Check result
    if [ $exit_code -eq 0 ]; then
        echo -e "${GREEN}✓ Plan ${name} completed successfully${NC}"
        echo "  Elapsed time: ${elapsed} seconds ($(printf '%02d:%02d:%02d' $((elapsed/3600)) $((elapsed%3600/60)) $((elapsed%60))))"
    else
        echo -e "${RED}✗ Plan ${name} failed${NC}"
        echo "  Exit code: $exit_code"
        echo "  Last output:"
        tail -20 stdout.log
        return 1
    fi
    
    # Copy results
    mkdir -p "$LOGDIR"
    cp stdout.log "$LOGDIR/stdout_${plan}.log"
    
    # Extract energy and timing info
    if [ -f "stdout.log" ]; then
        echo -e "${YELLOW}Summary:${NC}"
        grep -E "(Total energy|Time step|Memory|converged)" stdout.log | tail -10
    fi
    
    echo ""
}

# Function: Run quick test
run_quick_test() {
    echo -e "${YELLOW}Running QUICK test for C6H6 (reduced parameters)${NC}"
    echo -e "${YELLOW}Testing: LDA baseline + Plan C + CD-RI${NC}"
    
    for plan in lda c cdri; do
        local inputfile="inputfile_${plan}"
        if [ "$plan" = "c" ]; then
            inputfile="inputfile_plan_c"
        elif [ "$plan" = "cdri" ]; then
            inputfile="inputfile_cdri"
        fi
        
        if [ ! -f "${TESTDIR}/${inputfile}" ]; then
            continue
        fi
        
        local workdir="${TESTDIR}/work_${plan}_quick"
        mkdir -p "$workdir"
        cd "$workdir"
        
        # Copy and modify input
        cp "${TESTDIR}/${inputfile}" input.inp
        
        # Reduce parameters for quick test
        sed -i '' 's/nt = .*/nt = 10/g' input.inp
        sed -i '' 's/dl = .*/dl = 0.8, 0.8, 0.8/g' input.inp
        sed -i '' 's/nstate_frag = .*/nstate_frag = 10/g' input.inp
        
        # Link pseudo potentials
        if [ ! -f "C_rps.dat" ]; then
            ln -sf "${SALMON_DIR}/samples/C_rps.dat" "C_rps.dat"
        fi
        if [ ! -f "H_rps.dat" ]; then
            ln -sf "${SALMON_DIR}/samples/H_rps.dat" "H_rps.dat"
        fi
        
        # Run test
        echo -e "${YELLOW}Quick test - Plan ${plan}:${NC}"
        $SALMON_EXE < input.inp > stdout.log 2>&1
        
        if [ $? -eq 0 ]; then
            echo -e "${GREEN}✓ Quick test Plan ${plan} OK${NC}"
        else
            echo -e "${RED}✗ Quick test Plan ${plan} FAILED${NC}"
        fi
    done
}

# Function: Clean up
clean_tests() {
    echo -e "${YELLOW}Cleaning up C6H6 test directories...${NC}"
    rm -rf "${TESTDIR}/work_"*
    rm -rf "${TESTDIR}/logs"
    echo -e "${GREEN}Clean complete${NC}"
}

# Main execution
main() {
    echo -e "${BLUE}═════════════════════════════════════════════════════${NC}"
    echo -e "${BLUE}  C6H6 Test Suite - Plan C vs CD-RI Comparison${NC}"
    echo -e "${BLUE}═════════════════════════════════════════════════════${NC}"
    echo ""
    
    case "${1:-normal}" in
        quick)
            run_quick_test
            ;;
        clean)
            clean_tests
            ;;
        *)
            # Run tests: LDA (baseline) + HSE methods
            mkdir -p "$LOGDIR"
            
            echo -e "${YELLOW}═══════════════════════════════════════════════════${NC}"
            echo -e "${YELLOW}STEP 1: LDA Baseline Test (C6H6)${NC}"
            echo -e "${YELLOW}═══════════════════════════════════════════════════${NC}"
            run_test "lda"  "inputfile_lda"
            
            echo ""
            echo -e "${YELLOW}═══════════════════════════════════════════════════${NC}"
            echo -e "${YELLOW}STEP 2: HSE Methods (Plan C / CD-RI)${NC}"
            echo -e "${YELLOW}═══════════════════════════════════════════════════${NC}"
            
            echo -e "${YELLOW}NOTE: C6H6 test may take 20-40 minutes depending on system${NC}"
            echo ""
            
            run_test "c"    "inputfile_plan_c"
            run_test "cdri" "inputfile_cdri"
            
            # Summary
            echo -e "${BLUE}═════════════════════════════════════════════════════${NC}"
            echo -e "${BLUE}  Test Summary${NC}"
            echo -e "${BLUE}═════════════════════════════════════════════════════${NC}"
            echo -e "${YELLOW}Log files saved to: $LOGDIR${NC}"
            echo ""
            echo "To analyze results:"
            echo "  cd $TESTDIR"
            echo "  python3 analyze_results.py"
            ;;
    esac
}

# Run main function
main "$@"
