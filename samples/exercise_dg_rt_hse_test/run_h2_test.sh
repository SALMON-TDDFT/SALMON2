#!/bin/bash

#################################################################################
# H2 Test Suite for DG RT-TDDFT HSE Implementation
# 
# Purpose: Compare three HSE approaches (Plan A / Plan C / CD-RI)
# System: H2 molecule (simplest test case)
# 
# Usage:
#   ./run_h2_test.sh
#   ./run_h2_test.sh quick    # Quick test with reduced parameters
#   ./run_h2_test.sh clean    # Remove all output files
#
#################################################################################

set -e

SALMON_DIR="/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2"
MYBUILD_DIR="${SALMON_DIR}/mybuild"
SALMON_EXE="${MYBUILD_DIR}/src/salmon"
TESTDIR="$(cd "$(dirname "$0")" && pwd)/H2"
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
    echo -e "${BLUE}│ Testing Plan ${name} ($(echo $inputfile | rev | cut -d'_' -f1- | rev))${NC}"
    echo -e "${BLUE}└─────────────────────────────────────────────────────┘${NC}"
    
    # Create working directory
    local workdir="${TESTDIR}/work_${plan}"
    mkdir -p "$workdir"
    cd "$workdir"
    
    # Copy input file
    cp "${TESTDIR}/${inputfile}" input.inp
    
    # Prepare symbolic links for pseudo potentials
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
        echo "  Elapsed time: ${elapsed} seconds"
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
    echo -e "${YELLOW}Running QUICK test (reduced parameters)${NC}"
    echo -e "${YELLOW}Testing: LDA baseline + Plan C + CD-RI${NC}"
    
    # For quick test, we modify parameters temporarily
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
        sed -i '' 's/nt = .*/nt = 5/g' input.inp
        sed -i '' 's/dl = .*/dl = 0.8, 0.8, 0.8/g' input.inp
        sed -i '' 's/nstate_frag = .*/nstate_frag = 2/g' input.inp
        
        # Link pseudo potentials
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
    echo -e "${YELLOW}Cleaning up test directories...${NC}"
    rm -rf "${TESTDIR}/work_"*
    rm -rf "${TESTDIR}/logs"
    echo -e "${GREEN}Clean complete${NC}"
}

# Main execution
main() {
    echo -e "${BLUE}═════════════════════════════════════════════════════${NC}"
    echo -e "${BLUE}  H2 Test Suite - Plan A / Plan C / CD-RI Comparison${NC}"
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
            # Run all tests: LDA (baseline) + HSE methods
            mkdir -p "$LOGDIR"
            
            echo -e "${YELLOW}═══════════════════════════════════════════════════${NC}"
            echo -e "${YELLOW}STEP 1: LDA Baseline Test (sanity check)${NC}"
            echo -e "${YELLOW}═══════════════════════════════════════════════════${NC}"
            run_test "lda"   "inputfile_lda"
            
            echo ""
            echo -e "${YELLOW}═══════════════════════════════════════════════════${NC}"
            echo -e "${YELLOW}STEP 2: HSE Method Comparison (Plan A/C/CD-RI)${NC}"
            echo -e "${YELLOW}═══════════════════════════════════════════════════${NC}"
            run_test "a"    "inputfile_plan_a"
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
