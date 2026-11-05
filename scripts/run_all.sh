#!/usr/bin/env bash
set -e
cd "$(dirname "$0")/.."
python src/sweep.py
python src/analysis.py
python src/compare.py
