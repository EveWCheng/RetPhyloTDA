#!/usr/bin/env bash
# Full sim pipeline: simulate + find cycles, write analysis reports, then plot networks
# (reticulations coloured by the reports). sim_main.py wipes outputs/sim_phylo_outputs first.
# Usage: ./run_all.sh
set -euo pipefail
cd "$(dirname "$0")"

uv run python sim_main.py
uv run python main_result_analysis.py
Rscript ../Rscripts/plot_network.R ../Rscripts/config_main.yaml
