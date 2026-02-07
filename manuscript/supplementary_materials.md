# Supplementary Materials

## For: Noise-Assisted Transport Windows in Human Connectome Subgraphs: Evidence from Basal Ganglia and Motor Pathways

Oleg Dolgikh

---

## S1 Table. Per-Subject ENAQT Gain and Surrogate Enhancement (Hashed Seeds)

GKSL dynamics on T2 and T3 subgraphs. 12 hashed disorder seeds per subject, 8 strength-preserving surrogates per (subject, port, B). ΔG = median over seeds of (G_surr − G_orig); negative = original better. CI from hierarchical bootstrap (3000 resamples).

| Subject | Port | B (%) | G_orig (median) | ΔG (median) | CI low | CI high | Label |
|---------|------|-------|-----------------|-------------|--------|---------|-------|
| 100206 | T2 | 0 | 1.93 | −0.69 | −1.64 | −0.21 | Enhancement |
| 304727 | T2 | 0 | 1.93 | −0.68 | −1.04 | −0.28 | Enhancement |
| 314225 | T2 | 0 | 1.93 | −0.15 | −0.62 | +0.22 | Ambiguous |
| 333330 | T2 | 0 | 1.93 | −0.57 | −0.86 | −0.12 | Enhancement |
| 519647 | T2 | 0 | 1.93 | −1.05 | −1.84 | −0.37 | Enhancement |
| 100206 | T2 | 5 | 1.93 | −0.75 | −1.79 | −0.24 | Enhancement |
| 304727 | T2 | 5 | 1.93 | −0.86 | −1.38 | −0.38 | Enhancement |
| 314225 | T2 | 5 | 1.93 | −0.25 | −0.78 | +0.25 | Ambiguous |
| 333330 | T2 | 5 | 1.93 | −0.26 | −0.60 | +0.42 | Ambiguous |
| 519647 | T2 | 5 | 1.93 | −0.41 | −0.82 | −0.13 | Enhancement |
| 100206 | T2 | 10 | 1.93 | −1.03 | −1.64 | −0.40 | Enhancement |
| 304727 | T2 | 10 | 1.93 | −1.22 | −1.76 | −0.67 | Enhancement |
| 314225 | T2 | 10 | 1.93 | −0.19 | −0.74 | +0.18 | Ambiguous |
| 333330 | T2 | 10 | 1.93 | −0.41 | −0.64 | +0.22 | Ambiguous |
| 519647 | T2 | 10 | 1.93 | −0.38 | −0.63 | −0.01 | Enhancement |
| 100206 | T2 | 20 | 1.93 | −0.51 | −1.17 | −0.12 | Enhancement |
| 304727 | T2 | 20 | 1.93 | −0.76 | −1.50 | −0.38 | Enhancement |
| 314225 | T2 | 20 | 1.93 | −0.16 | −0.93 | +0.17 | Ambiguous |
| 333330 | T2 | 20 | 1.93 | −0.47 | −0.76 | −0.18 | Enhancement |
| 519647 | T2 | 20 | 1.93 | −0.39 | −1.11 | −0.10 | Enhancement |
| 100206 | T3 | 0 | 5.00 | −0.64 | −2.02 | +0.93 | Ambiguous |
| 304727 | T3 | 0 | 5.00 | −1.21 | −2.68 | −0.09 | Enhancement |
| 314225 | T3 | 0 | 5.00 | +0.33 | −2.34 | +1.55 | Ambiguous |
| 333330 | T3 | 0 | 5.00 | +0.37 | −0.52 | +1.30 | Ambiguous |
| 519647 | T3 | 0 | 5.00 | +0.40 | −1.16 | +2.23 | Ambiguous |

**Note**: G_orig median is across all subjects and seeds for the given port. "Enhancement" = 95% CI entirely below zero (original outperforms surrogates). "Ambiguous" = CI crosses zero.

---

## S2 Table. CTRW Classical Baseline Summary

Pauli master equation (classical continuous-time random walk) on same subgraphs. 21 log-spaced κ from 10⁻³ to 10. RK2 integration, dt = 0.05, T_end = 10.

| Subject | Port | B | Variant | R₀(κ_min) | R₀(κ_max) | κ* | G_CTRW | Monotonic |
|---------|------|---|---------|-----------|-----------|-----|--------|-----------|
| 100206 | T2 | 0 | orig | 57.8 | 150.5 | 10.0 | 0.0 | Yes |
| 100206 | T2 | 0 | surr0 | 55.8 | 150.3 | 10.0 | 0.0 | Yes |
| 100206 | T3 | 0 | orig | 262.7 | 513.2 | 10.0 | 0.0 | Yes |
| 333330 | T2 | 0 | orig | 56.3 | 150.4 | 10.0 | 0.0 | Yes |
| 333330 | T3 | 0 | orig | 281.3 | 512.9 | 10.0 | 0.0 | Yes |
| 304727 | T2 | 0 | orig | 56.3 | 150.5 | 10.0 | 0.0 | Yes |
| 304727 | T3 | 0 | orig | 279.6 | 513.2 | 10.0 | 0.0 | Yes |
| 314225 | T2 | 0 | orig | 60.3 | 150.3 | 10.0 | 0.0 | Yes |
| 314225 | T3 | 0 | orig | 212.0 | 512.3 | 10.0 | 0.0 | Yes |
| 519647 | T2 | 0 | orig | 55.7 | 150.4 | 10.0 | 0.0 | Yes |
| 519647 | T3 | 0 | orig | 222.8 | 512.4 | 10.0 | 0.0 | Yes |

**Total**: 64 cases tested (5 subjects × 2 ports × {orig + 3 surrogates} × {B = 0, B = 10}). All G_CTRW = 0.0, all curves monotonically increasing. No inverted-U window under classical dynamics.

**RK2 vs RK4 sanity check** (304727, T2, B = 0, orig): max |ΔR₀| = 8.74 × 10⁻⁵, max relative error = 5.84 × 10⁻⁷.

---

## S3 Table. Topological Metrics for T2 Surrogates (M1–M8)

Spearman rank correlation ρ(ΔG, ΔM) across 160 surrogates (5 subjects × 4 B-levels × 8 shuffles). ΔM = M_surr − M_orig; ΔG = G_orig − G_surr (positive = original better).

### S3a. Aggregate (all B pooled, N = 160)

| Metric | Description | ρ | p-value | Expected sign | Observed | Match |
|--------|------------|---|---------|--------------|----------|-------|
| M2 | Effective resistance (S→T) | +0.410 | 7.3 × 10⁻⁸ | + | + | Yes |
| M5 | Max edge BC at targets | +0.279 | 3.5 × 10⁻⁴ | + | + | Yes |
| M3 | Spectral gap (λ₂) | +0.238 | 2.4 × 10⁻³ | ? | + | — |
| M1 | Weighted shortest path | +0.197 | 1.2 × 10⁻² | + | + | Yes |
| M8 | IPR target overlap | +0.138 | 0.083 | − | + | No |
| M7 | Fiedler alignment | +0.057 | 0.478 | − | + | No |
| M6 | Strength Gini | −0.046 | 0.567 | + | − | No |
| M4 | BC Gini | −0.320 | 3.8 × 10⁻⁵ | ? | − | — |

### S3b. Stratified by B-freeze level (N = 40 per B)

| B | M1 ρ | M2 ρ | M3 ρ | M4 ρ | M5 ρ | M6 ρ | M7 ρ | M8 ρ |
|---|------|------|------|------|------|------|------|------|
| 0 | +0.36 | −0.11 | +0.17 | −0.33 | 0.00 | −0.17 | +0.27 | −0.16 |
| 5 | +0.19 | **+0.56** | +0.42 | −0.21 | +0.40 | +0.09 | −0.20 | +0.31 |
| 10 | +0.10 | **+0.64** | +0.33 | **−0.58** | **+0.49** | −0.03 | −0.18 | +0.18 |
| 20 | +0.22 | **+0.45** | +0.06 | −0.17 | +0.24 | −0.16 | +0.29 | +0.24 |

Bold: p < 0.005. Note M2 signal absent at B = 0 (ρ = −0.11, p = 0.50) but strong at B ≥ 5.

---

## S4 Table. Edge Weight–Betweenness Correlation (Sequestration Analysis)

Spearman correlation between edge weight (DTI streamline count) and edge betweenness centrality per subgraph.

| Subject | Port | N_nodes | N_edges | ρ (Spearman) | p-value | Significant? |
|---------|------|---------|---------|-------------|---------|-------------|
| 100206 | T2 | 20 | 152 | −0.325 | 4.3 × 10⁻⁵ | Yes |
| 100206 | T3 | 26 | 254 | −0.532 | 6.2 × 10⁻²⁰ | Yes |
| 304727 | T2 | 20 | 157 | −0.315 | 5.7 × 10⁻⁵ | Yes |
| 304727 | T3 | 26 | 280 | −0.479 | 1.9 × 10⁻¹⁷ | Yes |
| 314225 | T2 | 20 | 135 | −0.131 | 0.129 | No |
| 314225 | T3 | 26 | 226 | −0.385 | 2.0 × 10⁻⁹ | Yes |
| 333330 | T2 | 20 | 158 | −0.235 | 0.003 | Yes |
| 333330 | T3 | 26 | 237 | −0.393 | 3.4 × 10⁻¹⁰ | Yes |
| 519647 | T2 | 20 | 148 | −0.117 | 0.155 | No |
| 519647 | T3 | 26 | 227 | −0.360 | 2.3 × 10⁻⁸ | Yes |

T3: consistently negative (ρ = −0.36 to −0.53), all p < 10⁻⁸. T2: weaker, 3/5 significant. Interpretation: heavy DTI edges are not positioned on high-betweenness bridges (sequestration).

---

## S5 Table. B-Confounder Diagnostic: M2 Variance Analysis

Comparison of ΔM2 and ΔG variance between B = 0 (unconstrained) and B = 10 (backbone frozen).

| Variable | var(B = 0) | var(B = 10) | Ratio (B10/B0) | Levene p |
|----------|-----------|------------|---------------|----------|
| ΔM2 | 4.48 × 10⁻³ | 5.16 × 10⁻⁵ | **0.0115** (87× compression) | 0.012 |
| ΔG | 0.358 | 0.290 | 0.811 (stable) | 0.973 |

Spearman ρ(ΔG, ΔM2) at B = 0: −0.109 (p = 0.50); at B = 10: +0.643 (p = 7.7 × 10⁻⁶).

### Per-subject var(ΔM2) ratio (B = 10 / B = 0)

| Subject | var(ΔM2) B = 0 | var(ΔM2) B = 10 | Ratio |
|---------|---------------|-----------------|-------|
| 100206 | 2.66 × 10⁻⁵ | 1.48 × 10⁻⁵ | 0.556 |
| 304727 | 2.15 × 10⁻⁵ | 1.83 × 10⁻⁵ | 0.855 |
| 314225 | 1.05 × 10⁻² | 1.87 × 10⁻⁵ | 0.002 |
| 333330 | 1.88 × 10⁻³ | 1.02 × 10⁻⁵ | 0.005 |
| 519647 | 3.13 × 10⁻⁵ | 1.18 × 10⁻⁵ | 0.378 |

Subjects 314225 and 333330 show extreme variance compression (100–500×), indicating their bridge topology is highly sensitive to unconstrained shuffling.

---

## S6 Table. T0 Negative Control: ENAQT in Random Cortical Pairs

10 random cortical ROI pairs per subject, K = 15 subgraph, 1 hashed seed per pair. GKSL dynamics with canonical parameters.

| Subject | Median G | frac(G > 0) | Mean G |
|---------|----------|-------------|--------|
| 100206 | 0.000 | 3/10 (30%) | 0.22 |
| 304727 | 0.131 | 3/10 (30%) | 0.14 |
| 314225 | −0.010 | 2/10 (20%) | −0.03 |
| 333330 | 0.000 | 4/10 (40%) | 0.28 |
| 519647 | 0.149 | 2/10 (20%) | 0.19 |
| **Group** | **0.01** | **14/50 (28%)** | **0.16** |

**Comparison**: T2 functional: median G ≈ 1.93, frac > 0 = 98%. T3 functional: median G ≈ 5.0, frac > 0 = 90%. T0 random: median G ≈ 0.0, frac > 0 = 28%. ENAQT is pathway-specific.

---

## S7 Table. T2 Subgraph Topological Features

Per-subject structural properties of T2 subgraphs (K = 15 expansion, N = 20 nodes).

| Subject | N_edges | Density | BC Gini | Avg edge connectivity | Sequestration ρ |
|---------|---------|---------|---------|----------------------|----------------|
| 100206 | 152 | 0.80 | 0.56 | 13.2 | −0.325 |
| 304727 | 157 | 0.83 | 0.52 | 14.7 | −0.315 |
| 314225 | 135 | 0.71 | 0.63 | 10.8 | −0.131 (ns) |
| 333330 | 158 | 0.83 | 0.48 | 14.2 | −0.235 |
| 519647 | 148 | 0.78 | 0.61 | 12.0 | −0.117 (ns) |

Subject 314225 is the boundary case: lowest density, fewest edges, lowest edge connectivity, highest BC Gini, non-significant sequestration. This sparse, less redundant topology provides fewer alternative paths for noise-assisted transport.

---

## S8 Table. Group-Level Summary Statistics

Aggregated across 5 subjects, 12 hashed seeds (60 observations per port).

| Port | B | G_orig (median) | G_orig (mean) | ΔG (median) | frac(G_orig > 0) | frac(ΔG < 0) |
|------|---|----------------|--------------|-------------|------------------|--------------|
| T2 | 0 | 1.93 | 1.82 | −0.62 | 98.3% | 90.0% |
| T2 | 5 | 1.93 | 1.82 | −0.49 | 98.3% | 78.3% |
| T2 | 10 | 1.93 | 1.82 | −0.49 | 98.3% | 85.0% |
| T2 | 20 | 1.93 | 1.82 | −0.47 | 98.3% | 85.0% |
| T3 | 0 | 5.00 | 4.48 | −0.13 | 90.0% | 55.0% |
| T3 | 5 | 5.00 | 4.48 | −0.53 | 90.0% | 58.3% |
| T3 | 10 | 5.00 | 4.48 | −0.19 | 90.0% | 53.3% |
| T3 | 20 | 5.00 | 4.48 | +0.13 | 90.0% | 46.7% |

T2 shows consistent enhancement (ΔG < 0, frac > 78%) across all B. T3 shows mixed pattern (frac ≈ 50%), consistent with structure-driven transport.

---

## Supplementary Figures (Descriptions)

**S1 Fig. κ-Grid probe bias resolution.** Comparison of 5-point vs. 7-point κ-grids for T3. Left: 5-point grid misses peaks near κ = 0.063 in subjects 333330 and 519647. Right: 7-point grid resolves all inverted-U peaks. Per-subject R₀(κ) curves shown with seed-averaged bands.

**S2 Fig. ICC correction via hashed seeding.** Left: G distributions under shared seeding (ICC = 0.74 for T2); visible clustering by subject. Right: Hashed seeding (ICC = 0.03); distributions overlap. Bottom: forest plot showing per-subject ΔG with 95% CI under both seeding modes.

**S3 Fig. CTRW R₀(κ) curves.** All 12 panels for 3 remaining subjects (304727, 314225, 519647) × 2 ports (T2, T3) × 2 B-levels (0, 10). Each panel shows original + 3 surrogates. All curves monotonically increasing, confirming G_CTRW = 0.

**S4 Fig. Topological metric scatter plots.** ΔG vs. ΔM scatter for M1 (weighted shortest path), M5 (target edge betweenness), M4 (BC Gini), colored by B-level. Complement to main text Fig. 4C which shows M2.

**S5 Fig. Per-subject T2 enhancement detail.** Boxplots of per-seed ΔG for each subject separately, stratified by B. Shows that 314225 is consistently near zero while other subjects show robust negative ΔG.

---

## Supplementary Methods

### SM1. Canonical Hamiltonian Construction

The effective Hamiltonian is constructed from the DTI adjacency matrix W as follows:

1. **Weight normalization**: A = log1p(W) / log1p(W_max), where W_max is computed on the full 68-node connectome (not the subgraph). The log1p transform compresses the heavy-tailed streamline count distribution.

2. **Laplacian**: L = D − A, where D = diag(Σⱼ Aᵢⱼ) is the degree matrix of the normalized adjacency.

3. **Disorder**: ξᵢ ~ N(0, ε²) with ε = 3.0, drawn with subject-specific hashed seeds.

4. **Full Hamiltonian**: H = J · (L + diag(ξ)), with coupling constant J = 0.1.

5. **Non-Hermitian sinks**: H_eff = H − (i/2) · diag(γ), where γᵢ = γ_T = 1.0 for target nodes, γᵢ = γ_L = 5 × 10⁻⁴ for all others.

### SM2. Iterative Proportional Fitting (IPF) for Surrogates

Strength-preserving shuffles use IPF to ensure that shuffled weight matrices W_surr have identical node strengths to W_orig:

1. Permute non-frozen edge weights randomly (seeded).
2. Initialize W_surr with permuted weights.
3. Iterate: (a) scale rows to match target row sums; (b) scale columns to match target column sums; (c) symmetrize: W_surr = (W_surr + W_surr^T) / 2.
4. Repeat until max |strength_surr − strength_orig| < 10⁻³.

Convergence typically occurs within 50–100 iterations for 20-node subgraphs.

### SM3. ICC Computation

Intraclass correlation was computed via one-way random-effects ANOVA:

ICC = (MS_between − MS_within) / (MS_between + (k − 1) · MS_within)

where k = number of seeds per subject. Design effect DE = 1 + (k − 1) · ICC. Effective sample size N_eff = N / DE.

---

*End of Supplementary Materials*
