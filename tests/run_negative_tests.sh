#!/usr/bin/env bash
# tests/run_negative_tests.sh
# Automated regression test suite asserting expected failure codes and diagnostics for negative test cases.

set -euo pipefail

NF_BIN="/opt/homebrew/bin/nextflow"
if ! command -v "$NF_BIN" &> /dev/null; then
    NF_BIN="nextflow"
fi

echo "=========================================================="
echo "Running LongVarMeth Negative Test Suite"
echo "=========================================================="

FAILED_TESTS=0
TOTAL_TESTS=0

run_neg_test() {
    local test_name="$1"
    local expected_pattern="$2"
    shift 2
    local cmd=("$@")
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    echo -n "Test $TOTAL_TESTS: $test_name ... "
    
    local output
    set +e
    output=$("${cmd[@]}" 2>&1)
    local status=$?
    set -e
    
    if [ "$status" -eq 0 ]; then
        echo "FAILED (Expected non-zero exit code, got 0)"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        return
    fi
    
    if echo "$output" | grep -q "$expected_pattern"; then
        echo "PASSED (exit $status, matched '$expected_pattern')"
    else
        echo "FAILED (Expected pattern '$expected_pattern' not found in output)"
        echo "$output"
        FAILED_TESTS=$((FAILED_TESTS + 1))
    fi
}

# Test 1: Missing --ref_nuc in reference-based mode
run_neg_test \
    "Missing --ref_nuc in reference_based mode" \
    "Reference FASTA file specified by --ref_nuc does not exist" \
    "$NF_BIN" run main.nf -profile test -stub --run_type reference_based --ref_nuc ""

# Test 2: Platform mismatch (PacBio BAM on --platform ont)
cat <<EOF > tests/_tmp_mismatch.csv
run_id,input_path,input_type,kit,model,modified_bases_model
mismatch_sample,tests/data/sample.pacbio.bam,pacbio_bam,,,
EOF

run_neg_test \
    "Platform mismatch rejection (PacBio BAM with --platform ont)" \
    "PacBio format 'pacbio_bam' supplied to ONT workflow" \
    "$NF_BIN" run main.nf -profile test -stub --platform ont --input tests/_tmp_mismatch.csv

rm -f tests/_tmp_mismatch.csv

# Test 3: Missing required column in samplesheet
cat <<EOF > tests/_tmp_missing_col.csv
run_id,input_path
sample1,tests/data/reads.pod5
EOF

run_neg_test \
    "Missing required column in samplesheet" \
    "missing required column 'input_type'" \
    "$NF_BIN" run main.nf -profile test -stub --input tests/_tmp_missing_col.csv

rm -f tests/_tmp_missing_col.csv

echo "=========================================================="
if [ "$FAILED_TESTS" -eq 0 ]; then
    echo "ALL $TOTAL_TESTS NEGATIVE TESTS PASSED SUCCESSFULLY"
    echo "=========================================================="
    exit 0
else
    echo "$FAILED_TESTS / $TOTAL_TESTS TESTS FAILED"
    echo "=========================================================="
    exit 1
fi

