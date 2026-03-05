#!/bin/bash
#
# DG-Fragment RT Physics Validation Script
# Runs both DG-Fragment and Conventional RT on the same system and compares outputs
#

set -e

SALMON_BIN="/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/mybuild/salmon"
TEST_DIR="/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/samples/exercise_dg_rt_hse_test/H2"

echo "========================================================================"
echo "  DG-Fragment RT Physics Validation"
echo "========================================================================"
echo ""
echo "Test Location: $TEST_DIR"
echo "SALMON Binary: $SALMON_BIN"
echo ""

# Check that necessary files exist
if [ ! -f "$SALMON_BIN" ]; then
    echo "ERROR: SALMON binary not found at $SALMON_BIN"
    exit 1
fi

if [ ! -d "$TEST_DIR" ]; then
    echo "ERROR: Test directory not found at $TEST_DIR"
    exit 1
fi

cd "$TEST_DIR"

# Clean up old test outputs
echo "[Step 1/4] Cleaning previous test outputs..."
rm -f H2_periodic_20_dg_new_param_rt*.data
rm -f dg_physics_test_*.log
echo "  ✓ Cleaned"

# Run DG-Fragment RT test
echo ""
echo "[Step 2/4] Running DG-Fragment RT..."
echo "  Input: inputfile_h2_periodic_20_dg_new_param"
echo "  Duration: ~5 minutes on workstation"
echo ""

START_TIME=$(date +%s)
$SALMON_BIN < inputfile_h2_periodic_20_dg_new_param > dg_physics_test_run.log 2>&1
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

if [ -f "H2_periodic_20_dg_new_param_rt.data" ]; then
    LINES=$(grep -c "^[[:space:]]*[0-9]" H2_periodic_20_dg_new_param_rt.data)
    echo "  ✓ Test completed in ${ELAPSED}s"
    echo "  ✓ Output file generated: $LINES timesteps"
else
    echo "  ✗ ERROR: Output file not generated"
    tail -30 dg_physics_test_run.log
    exit 1
fi

# Check conventional RT data exists
echo ""
echo "[Step 3/4] Checking Conventional RT reference data..."
if [ -f "H2_periodic_20_conventional_rt_rt.data" ]; then
    CONV_LINES=$(grep -c "^[[:space:]]*[0-9]" H2_periodic_20_conventional_rt_rt.data)
    echo "  ✓ Reference data found: $CONV_LINES timesteps"
else
    echo "  ⚠ WARNING: Conventional RT reference data not found"
    echo "    To generate reference: $SALMON_BIN < inputfile_h2_periodic_20_conventional"
fi

# Run Python analysis
echo ""
echo "[Step 4/4] Running physics analysis..."
python3 physics_validation.py

echo ""
echo "========================================================================"
echo "  Validation Complete"
echo "========================================================================"
echo ""
echo "Output files generated:"
ls -lh H2_periodic_20_dg_new_param_rt*.data 2>/dev/null | awk '{print "  " $9 " (" $5 ")"}'
echo ""
echo "For detailed analysis, see: PHYSICS_VALIDATION_REPORT.md"
