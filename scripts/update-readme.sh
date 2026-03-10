#!/usr/bin/env bash
# scripts/update-readme.sh
# Runs the full reproducible pipeline to generate a fresh README.md

set -e

echo "🚀 Starting robscale README update pipeline..."

# Ensure we are in the project root
cd "$(dirname "$0")/.."

# Run the targets pipeline
# This will compile robscale variationally and run benchmarks
Rscript -e "targets::tar_make()"

echo "✅ README.md has been successfully updated and is honest!"
