#!/bin/bash
#
# Run all tests for the Romanian palatalization pipeline
#
# This script runs pytest tests to validate:
# - Pipeline invariants (TP domain, NDE classification, etc.)
# - TP summary table snapshots (stable counts)
#
# Usage:
#   ./run_tests.sh           # Run all tests
#   ./run_tests.sh -v        # Verbose output
#   ./run_tests.sh -k test_gimpe  # Run specific test pattern

set -euo pipefail

# Check if pytest is installed
if ! python3 -c "import pytest" 2>/dev/null; then
    echo "pytest not found. Installing..."
    pip3 install pytest pandas
fi

# Run tests
echo "Running Romanian palatalization pipeline tests..."
echo "=================================================="
python3 -m pytest tests/ "$@"

echo ""
echo "=================================================="
echo "All tests passed!"
