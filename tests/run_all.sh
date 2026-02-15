#!/bin/bash
# Run all Guitar tests
# Usage: micromamba run -n guitar bash tests/run_all.sh

set -e
cd "$(dirname "$0")/.."

PASS=0
FAIL=0

for test in tests/test_basic.R tests/test_multigroup.R tests/test_prebuild.R; do
    echo ""
    echo "========================================"
    echo "Running: $test"
    echo "========================================"
    if Rscript "$test"; then
        PASS=$((PASS + 1))
    else
        FAIL=$((FAIL + 1))
    fi
done

echo ""
echo "========================================"
echo "Test suites: $PASS passed, $FAIL failed"
echo "========================================"
[ $FAIL -eq 0 ]
