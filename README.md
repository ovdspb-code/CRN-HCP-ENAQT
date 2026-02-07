# CRN-HCP-ENAQT: Noise-Assisted Transport in Human Connectome Subgraphs

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Code and data for:

> Dolgikh O. (2026). Noise-Assisted Transport Windows in Human Connectome Subgraphs: Evidence from Basal Ganglia and Motor Pathways. *PLOS Computational Biology* (submitted).

## Overview

This repository contains the complete simulation pipeline and analysis code for testing noise-assisted transport (ENAQT-type) phenomena in human structural connectomes derived from the Human Connectome Project (HCP). The study examines two anatomical pathways — basal ganglia gating (T2) and thalamocortical motor relay (T3) — using GKSL open-system dynamics as a functional proxy for wave-like transport with tunable dephasing.

## Repository Structure

```
CRN-HCP-ENAQT/
├── README.md
├── LICENSE                          # MIT
├── requirements.txt                 # Python dependencies
│
├── inputs/
│   ├── graphml/                     # HCP connectomes (5 subjects, .graphml)
│   │   ├── 100206_repeated10_scale125.graphml
│   │   ├── 304727_repeated10_scale125.graphml
│   │   ├── 314225_repeated10_scale125.graphml
│   │   ├── 333330_repeated10_scale125.graphml
│   │   └── 519647_repeated10_scale125.graphml
│   └── human_task_battery_v1.json   # Port definitions (T0/T1/T2/T3)
│
├── code/
│   ├── step11_hashed_only_stats.py          # Main GKSL + surrogate pipeline (P0–P4)
│   ├── step12b_t2_motifs_metriclist.py      # Topological metric analysis (M1–M8)
│   ├── step12c_Bconfound_diagnostic.py      # B-freeze confounder diagnostics
│   ├── step13_edge_weight_betweenness_correlation.py  # Sequestration analysis
│   ├── step14_ctrw_surrogate_control.py     # CTRW baseline (2 subjects)
│   └── step15_ctrw_surrogates_remaining3.py # CTRW baseline (3 remaining)
│
├── results/
│   ├── Step11_HashedOnly12Stats/
│   │   ├── all_trials_hashed12_clean.csv    # 4320 trials, all G values
│   │   ├── per_seed_summary.csv             # Per-seed ΔG
│   │   ├── per_subject_labels.csv           # Enhancement/structure-driven labels
│   │   └── group_summary.csv                # Group-level statistics
│   │
│   ├── Step12b_T2Motifs_MetricList/
│   │   ├── t2_surrogate_metrics_deltas.csv  # M1–M8 for 160 surrogates
│   │   ├── t2_metric_correlations_raw.csv   # Aggregate Spearman ρ
│   │   └── t2_metric_correlations_by_B.csv  # B-stratified correlations
│   │
│   ├── Step12c_BconfoundDiagnostic/
│   │   ├── B_stats_all.csv                  # Variance analysis per B
│   │   └── diagnostic_key_results.json      # Key var-ratio numbers
│   │
│   ├── Step13_EdgeWeightBetweenness/
│   │   └── edge_weight_betweenness_correlation.csv
│   │
│   ├── Step14_CTRW_Surrogates/
│   │   └── ctw_surrogates_summary.csv       # 16 cases (2 subjects)
│   │
│   └── Step15_CTRW_Remaining3/
│       ├── ctrw_surrogates_remaining3_summary.csv   # 48 cases
│       ├── ctrw_surrogates_remaining3_curves.csv    # Full R₀(κ) curves
│       └── rk2_vs_rk4_sanitycheck_304727_T2_B0_orig.csv
│
├── figures/
│   ├── Fig1_ENAQT_window/           # Placeholder for main figure panels
│   ├── Fig2_CTRW_control/
│   ├── Fig3_pathway_specific/
│   └── Fig4_topo_mechanisms/
│
└── supplementary/
    └── tables/                      # S1–S7 Tables (CSV format)
```

## Quick Start

### Dependencies

```bash
pip install numpy scipy networkx pandas matplotlib
```

### Reproduce Main Results

```bash
# Step 11: GKSL + surrogates (hashed seeds, 12 per subject)
python code/step11_hashed_only_stats.py \
  --input_dir inputs/graphml \
  --battery inputs/human_task_battery_v1.json \
  --output_dir results/Step11_HashedOnly12Stats

# Step 12b: Topological metrics on T2 surrogates
python code/step12b_t2_motifs_metriclist.py \
  --input_dir inputs/graphml \
  --battery inputs/human_task_battery_v1.json \
  --trials results/Step11_HashedOnly12Stats/all_trials_hashed12_clean.csv \
  --output_dir results/Step12b_T2Motifs_MetricList

# Step 15: CTRW classical baseline
python code/step15_ctrw_surrogates_remaining3.py \
  --input_dir inputs/graphml \
  --battery inputs/human_task_battery_v1.json \
  --output_dir results/Step15_CTRW_Remaining3
```

## Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| J | 0.1 | Coupling constant |
| ε_dis | 3.0 | Disorder amplitude (σ of Gaussian) |
| γ_T | 1.0 | Target sink rate |
| γ_L | 5×10⁻⁴ | Non-target leakage rate |
| dt | 0.05 | Integration timestep |
| T_end | 10.0 | Total integration time |
| K_extra | 15 | Intermediate nodes in subgraph |
| n_shuffles | 8 | Surrogates per (subject, port, B) |
| n_seeds | 12 | Hashed disorder seeds per subject |
| B_list | [0, 5, 10, 20] | Betweenness-freeze percentages |

## Canonical Hamiltonian (Canonical-H)

```
A = log1p(W) / log1p(W_max_full_graph)
L = D − A
H = J · (L + diag(ξ))
```

where ξ ~ N(0, ε²) with subject-specific hashing:
```python
seed_eff = int(SHA256(f"{base_seed}_{subject_id}").hexdigest()[:8], 16)
```

## Citation

```bibtex
@article{dolgikh2026noise,
  title={Noise-Assisted Transport Windows in Human Connectome Subgraphs:
         Evidence from Basal Ganglia and Motor Pathways},
  author={Dolgikh, Oleg},
  journal={PLOS Computational Biology},
  year={2026},
  note={Submitted}
}
```

## Related Work

- CRN Core Theory: [Dolgikh 2024, Zenodo](https://doi.org/10.5281/zenodo.18249250)
- Companion ICC Methods Paper: Dolgikh (in preparation)

## License

MIT License. See [LICENSE](LICENSE) for details.

HCP data are subject to the [HCP Data Use Agreement](https://www.humanconnectome.org/study/hcp-young-adult/data-use-terms).
