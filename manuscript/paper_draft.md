# Noise-Assisted Transport Windows in Human Connectome Subgraphs: Evidence from Basal Ganglia and Motor Pathways

**Oleg Dolgikh**

Independent Researcher, Sant Cugat del Vallès, Spain

ORCID: 0009-0008-0159-1718

Correspondence: hanlon@occam.world

---

## Abstract

Noise-assisted transport — where intermediate environmental coupling enhances signal delivery beyond both coherent and diffusive limits — has been demonstrated in photosynthetic complexes and small biological networks, but not in human brain connectivity data. Here we test whether human structural connectomes derived from diffusion tensor imaging (DTI) exhibit such transport windows in functionally defined neural pathways. Using an open-system Lindblad master equation as a mechanistic proxy for wave-like dynamics with tunable dephasing, we analyze two anatomical ports in five Human Connectome Project subjects: T2 (basal ganglia gating loop) and T3 (thalamocortical motor relay). We report three principal findings. First, inverted-U transport profiles (transport efficiency vs. dephasing rate κ) appear robustly in both pathways across all five subjects, and are entirely absent under a classical continuous-time random walk (CTRW) control on the same graphs. Second, the two pathways exhibit a functional dissociation: in T2, native DTI edge weights outperform strength-preserving surrogates (median ΔG = +0.62, hierarchical bootstrap p < 0.01), indicating optimization beyond random wiring, whereas T3 transport is structure-driven, with weight perturbations producing negligible change. Third, within-surrogate analysis identifies effective resistance and target-adjacent edge betweenness as the topological predictors of transport loss, contingent on preservation of the routing backbone — suggesting a two-level architecture of skeleton plus conductance tuning. We further report a methodological finding: shared random disorder seeds across simulated subjects inflated intraclass correlation from 0.03 to 0.74, creating pseudoreplication that could affect any connectome simulation study with nested designs (detailed in a companion methods paper). All code and data are publicly available. These results demonstrate that noise-assisted transport windows are a reproducible, pathway-specific property of human DTI connectomes, consistent with the Coherent Resonant Netting (CRN) hypothesis that biological networks exploit transient wave-like dynamics for energy-efficient signal routing.

**Keywords:** noise-assisted transport; open-system dynamics; connectome; dephasing; basal ganglia; transport efficiency; surrogate controls; DTI tractography

---

## 1. Introduction

Neural computation is metabolically expensive. Action potentials, synaptic transmission, and dendritic integration consume a substantial fraction of the brain's energy budget, and any search process that relies on repeated spiking updates incurs costs proportional to the number of hypotheses explored [1,2]. This creates an energy–information bottleneck: as decision spaces grow, purely diffusive exploration becomes prohibitively slow on structured graphs, while brute-force acceleration via increased spiking rates risks exceeding metabolic capacity [3].

A complementary principle has emerged from studies of energy transfer in photosynthetic complexes: environmental noise, rather than degrading transport, can enhance it. In the phenomenon termed Environment-Assisted Quantum Transport (ENAQT), intermediate dephasing — neither too coherent nor too incoherent — produces optimal delivery of excitation energy to reaction centers [4–6]. The key signature is an inverted-U profile: transport efficiency peaks at an intermediate noise level and declines at both extremes. Related phenomena have been observed in engineered photonic lattices [7] and predicted for disordered networks more generally [8].

The Coherent Resonant Netting (CRN) framework hypothesizes that this principle extends beyond molecular systems to neural circuit architecture [9]. CRN proposes a two-stage decision architecture: Stage I performs low-amplitude, wave-like exploration on the structural connectome ("netting"), rapidly pruning hypothesis spaces at low energetic cost, while Stage II commits to a decision via high-gain spiking fixation. Critically, CRN does not claim microscopic quantum coherence in brain tissue. Instead, it uses the Gorini-Kossakowski-Sudarshan-Lindblad (GKSL) master equation as a *functional proxy* for any transient wave-like process with tunable damping: at zero dephasing (κ → 0) the model is fully coherent; at high dephasing (κ → ∞) it converges to classical diffusion. Any linear wave substrate with tunable phase randomization — including dendritic cable oscillations, subthreshold membrane potential fluctuations, or local field potential coupling — could in principle support the same qualitative signatures [9, §2.3].

Prior CRN work demonstrated noise-assisted transport windows and disorder-enhanced selectivity (DES) in C. elegans touch circuits, Drosophila larva mushroom body, and mouse cortical proxy models [9]. The present study asks whether the same signatures appear in *human* DTI connectomes — networks orders of magnitude larger, derived from macroscopic neuroimaging rather than cellular-resolution tracing.

We focus on two anatomically defined pathways: **T2** (basal ganglia → thalamus/pallidum), a gating circuit implicated in action selection and cognitive flexibility [10], and **T3** (thalamus → motor cortex), a relay circuit supporting motor execution [11]. These ports were chosen because they represent functionally distinct computational demands — selective gating (T2) versus reliable broadcast (T3) — and because their anatomy is well-characterized across subjects.

Our analysis proceeds through a systematic validation pipeline: (P0) resolution of κ-grid probe bias; (P1) classical CTRW baseline to confirm that inverted-U profiles require coherent dynamics; (P2) surrogate controls to test whether native DTI weights are specifically optimized for transport; (P3) negative-control random cortical pairs to test pathway specificity; and (P4) elimination of pseudoreplication via hashed disorder seeding. We then characterize the topological predictors of transport loss in surrogates and identify a two-level architecture — routing backbone plus conductance tuning — that explains the pathway-specific findings.

---

## 2. Materials and Methods

### 2.1 Dataset and Functional Ports

Structural connectomes were obtained from five participants in the Human Connectome Project (HCP) [12]: subjects 100206, 304727, 314225, 333330, and 519647. Parcellation followed the Desikan-Killiany atlas (68 cortical and subcortical ROIs per hemisphere); edge weights represent DTI streamline counts.

Two functional ports were defined based on known basal ganglia–thalamocortical circuit anatomy:

- **T2 (Basal Ganglia loop)**: Sources = {left caudate, putamen, accumbens}; Targets = {left pallidum, thalamus}.
- **T3 (Motor output)**: Source = {left thalamus}; Targets = {left precentral gyrus, paracentral lobule}.

A negative-control set **T0** consisted of 10 random cortical ROI pairs per subject.

For each port, subgraphs were constructed by selecting core source and target nodes plus K = 15 highest-weight intermediate nodes from the ipsilateral hemisphere, yielding subgraphs of approximately 20 nodes (T2) and 26 nodes (T3). Intermediate node selection was deterministic: ranked by total edge weight to core nodes, with node ID as tiebreaker.

### 2.2 Open-System Transport Model (GKSL Proxy)

Stage-I dynamics were modeled with the GKSL master equation for density matrix ρ(t):

dρ/dt = −i[H, ρ] + κ·D_deph(ρ) + L_sink(ρ)

where H = J·(L + diag(ξ)) is the effective Hamiltonian constructed from the weighted graph Laplacian L = D − A with normalized adjacency A = log1p(W) / log1p(W_max), diagonal disorder ξ ~ N(0, ε²) with ε = 3.0, coupling constant J = 0.1, dephasing superoperator D_deph(ρ) = −κ·Σᵢ[|i⟩⟨i|, [|i⟩⟨i|, ρ]] implementing pure site dephasing at rate κ, and sink superoperator L_sink with target absorption rate γ_T = 1.0 and non-target leakage rate γ_L = 5 × 10⁻⁴.

Integration used a fourth-order Runge-Kutta (RK4) scheme with dt = 0.05 over T_end = 10.0 time units, with burn-in t = 0.5.

The primary observable was R₀_max(κ) = max_t [P_T(t) − P_D(t)], the cumulative target absorption minus leakage penalty. The ENAQT gain metric was defined as:

G = R₀_max(κ*) − R₀_edge

where κ* is the dephasing rate maximizing R₀_max and R₀_edge = min(R₀|κ_min, R₀|κ_max). An inverted-U profile yields G > 0.

### 2.3 Classical Baseline (CTRW)

To verify that inverted-U profiles require coherent dynamics, we replaced GKSL evolution with a classical continuous-time random walk (Pauli master equation):

dp/dt = −J·L·p + κ·(uniform redistribution − p) − diag(γ)·p

This removes all off-diagonal coherences (no interference, no wave-like transport). Integration parameters, subgraphs, and observables were identical to the GKSL model. A 21-point log-spaced κ-grid from 10⁻³ to 10 was used to ensure no peaks were missed.

### 2.4 Surrogate Controls

Strength-preserving weight shuffles tested whether ENAQT advantage depends on specific edge weight configurations rather than node strength profiles alone. For each subgraph:

1. Edges were ranked by raw weight; the top 5% were frozen (preserving heavy-tail structure).
2. Edges were optionally ranked by unweighted edge betweenness centrality; the top B% (B ∈ {0, 5, 10, 20}) were additionally frozen (preserving routing backbone).
3. Remaining edge weights were permuted and iteratively rescaled via Iterative Proportional Fitting (IPF) to restore original node strengths within tolerance 10⁻³.

Eight surrogates were generated per (subject, port, B) combination using deterministic seeding: RNG seed = MD5("surr" || subject || port || B || shuffle_id). The metric ΔG = G_orig − median(G_surr) quantified enhancement (ΔG > 0 = native weights better).

### 2.5 Disorder Seeding and ICC Control

Diagonal disorder ξ was drawn from N(0, ε²) on the full connectome graph, then mapped to subgraph nodes. To prevent pseudoreplication from shared random seeds, each subject received subject-specific disorder via hashing:

seed_eff = int(SHA256(f"{base_seed}_{subject_id}").hexdigest()[:8], 16)

Twelve base seeds (20250212–20250223) provided 12 independent disorder realizations per subject. Intraclass correlation (ICC) was measured via nested ANOVA. A companion methods paper details this correction and its impact on effective sample size (Dolgikh, in preparation).

### 2.6 κ-Grid Resolution

Preliminary analysis revealed probe bias in T3: a 5-point κ-grid missed inverted-U peaks near κ ≈ 0.063 in 2/5 subjects. All reported results use a 7-point grid for T3: {0.001, 0.010, 0.063, 0.251, 0.631, 1.0, 10.0} and a 5-point grid for T2: {0.001, 0.251, 0.631, 1.0, 10.0}. CTRW controls used a denser 21-point log-spaced grid.

### 2.7 Network Topology Metrics

To identify topological predictors of ENAQT loss in surrogates, eight metrics were computed on each surrogate and original subgraph (T2 only, where enhancement was observed):

- **M1**: Mean weighted shortest path length (source → target), using 1/w as distance.
- **M2**: Effective resistance between source and target sets, computed from the pseudoinverse of the graph Laplacian.
- **M3**: Spectral gap (second-smallest eigenvalue of the weighted Laplacian).
- **M4**: Gini coefficient of edge betweenness centrality distribution.
- **M5**: Maximum edge betweenness among edges incident to target nodes.
- **M6**: Mean node strength Gini coefficient.
- **M7**: Fiedler vector alignment (correlation of |Fiedler vector| with target indicator).
- **M8**: Inverse Participation Ratio of Laplacian eigenvectors projected onto target subspace.

Spearman rank correlations ρ(ΔG, ΔM) were computed across 160 surrogates (5 subjects × 4 B-levels × 8 shuffles), and stratified by B to diagnose confounding.

### 2.8 Statistical Framework

Group-level significance of T2 enhancement used subject as the unit of analysis (N = 5), with per-subject ΔG estimated as the median across 12 hashed seeds. Confidence intervals on per-subject ΔG were obtained via hierarchical bootstrap (3000 resamples, sampling seeds within subjects, then subjects across group) following recommendations for nested neuroscience data [13].

### 2.9 Code and Data Availability

All simulation code, configuration files, and output datasets are available at [GitHub repository URL] under MIT license. Raw HCP connectomes are available through the HCP data use agreement [12].

---

## 3. Results

### 3.1 Inverted-U Transport Profiles Exist in Human Connectomes (P0)

All five subjects exhibited inverted-U R₀_max(κ) profiles in both T2 and T3 pathways under GKSL dynamics (Fig. 1). Transport efficiency peaked at intermediate dephasing (κ* ≈ 0.25–1.0 for T2; κ* ≈ 0.06–0.63 for T3) and declined at both coherent (κ → 0) and diffusive (κ → ∞) extremes.

Group-level ENAQT gain: T2 median G = 1.93 (IQR 0.88–2.40), fraction G > 0: 59/60 seeds (98%). T3 median G = 5.00 (IQR 3.2–6.5), fraction G > 0: 54/60 seeds (90%).

The 7-point κ-grid was critical for T3: two subjects (333330, 519647) showed artificially low G with the initial 5-point grid due to missed peaks near κ = 0.063. After grid correction, all subjects showed positive G.

### 3.2 Classical Transport Shows No Inverted-U (P1)

Under the classical CTRW model, R₀_max(κ) was monotonically increasing in all 64 tested configurations (5 subjects × 2 ports × {original + 3 surrogates} × {B = 0, B = 10}; Fig. 2). Maximum R₀ always occurred at κ_max = 10 (maximum randomization). No inverted-U window was observed (G_CTRW = 0.0 in all cases).

RK2-vs-RK4 sanity check on a representative case confirmed negligible integrator error (max |ΔR₀| = 8.7 × 10⁻⁵, relative error < 10⁻⁶).

These results establish that inverted-U profiles are specific to coherent (GKSL) dynamics and cannot be reproduced by classical stochastic diffusion on the same graphs with the same parameters.

### 3.3 Pathway-Specific Enhancement and Structure-Driven Transport (P2)

#### T2 (Basal Ganglia): Topology-Specific Enhancement

Across all B-levels (unconstrained and constrained surrogates), native DTI weights outperformed strength-preserving surrogates in T2. At B = 0 (strongest null — only top-5% weight edges frozen):

- 4/5 subjects showed significant enhancement (hierarchical bootstrap 95% CI excluding zero).
- 1/5 (subject 314225) was ambiguous (ΔG = −0.15, CI [−0.62, +0.22]).
- Group median ΔG = −0.62 (surrogates worse than original; convention: negative delta = original better).
- Group-level test (subjects as unit): t(4) = 4.2, p = 0.009.

Enhancement was robust across B-levels: median ΔG remained negative (original superior) at B = 5, 10, and 20% (Fig. 3).

#### T3 (Motor Output): Structure-Driven

T3 showed mixed results with no consistent enhancement:

- At B = 0: 1/5 enhancement, 4/5 structure-driven or ambiguous.
- Median ΔG near zero; group CI crossed zero.
- B-freeze modulated effect size but did not reverse the pattern.

#### Negative Control (T0): Random Pairs

Ten random cortical ROI pairs per subject showed minimal ENAQT: median G ≈ 0.0, fraction G > 0 = 28% (vs. 98% for T2, 90% for T3). This confirms that inverted-U profiles are pathway-specific, not an artifact of subgraph extraction or DTI processing.

### 3.4 Topological Predictors of Transport Loss (Step 12b)

Within-surrogate analysis on T2 (160 surrogates, 5 subjects × 4 B × 8 shuffles) identified three significant predictors of ENAQT loss (Fig. 4):

| Metric | Spearman ρ | p-value | Predicted sign | Match |
|--------|-----------|---------|---------------|-------|
| M2 (effective resistance) | +0.41 | 7 × 10⁻⁸ | + | Yes |
| M5 (target edge betweenness) | +0.28 | 3 × 10⁻⁴ | + | Yes |
| M4 (BC Gini) | −0.32 | 4 × 10⁻⁵ | ? | — |
| M1 (weighted shortest path) | +0.20 | 0.012 | + | Yes |

Surrogates with higher effective resistance (worse conductance), more concentrated betweenness at target nodes (bottleneck), or longer weighted paths showed greater ENAQT loss. The negative M4 signal indicates that more uniform betweenness distribution (lower Gini) degraded performance — counterintuitively, directed heterogeneity (flow channelization) supports ENAQT.

Spectral metrics (M7 Fiedler alignment, M8 IPR overlap) were non-significant (p > 0.4), likely due to the high density of 20-node subgraphs (density 0.7–0.8).

#### B-Dependence: Two-Level Architecture

Stratification by betweenness-freeze level revealed that the M2 signal was B-dependent: ρ = −0.11 (ns) at B = 0, rising to ρ = +0.64 (p < 10⁻⁵) at B = 10. Diagnostic analysis confirmed that var(ΔM2) compressed 87-fold between B = 0 and B = 10 (Levene p = 0.012), while var(ΔG) remained stable (ratio 0.81, Levene p = 0.97).

This dissociation — compressed predictor variance without compressed outcome variance — rules out trivial range-compression confounding. Instead, it supports a **two-level architecture**: (i) a routing backbone (bridges and hub connectivity, captured by betweenness-freeze) sets baseline transport capacity; (ii) fine-grained edge weight distribution (conductance, captured by M2) optimizes performance within the preserved skeleton. Effective resistance informs the fine-tuning level; target betweenness (M5) and heterogeneity (M4) remain informative across all B.

### 3.5 Weight–Betweenness Anti-Correlation (Sequestration)

Edge weight and edge betweenness centrality showed consistent negative Spearman correlation in T3 subgraphs (ρ = −0.36 to −0.53, all p < 10⁻⁸ across 5 subjects). In T2, the correlation was weaker and non-significant in 2/5 subjects (ρ = −0.12 to −0.33). This "sequestration" pattern — heavy DTI edges occupying local cores rather than global bridges — is consistent with biological routing strategies that avoid chokepoint formation on high-capacity links.

### 3.6 Hashed Seeding Eliminates Pseudoreplication (P4)

Shared disorder seeds inflated ICC to 0.74 for T2 (design effect DE = 4.7, effective N ≈ 6 out of 30). Hashed (subject-specific) seeding reduced ICC to 0.03 (DE = 1.17, effective N ≈ 25). All reported statistics use hashed mode. The main enhancement claim for T2 survives both seed-level (t(29) = 14.4, p < 10⁻⁴) and conservative subject-level (t(4) = 4.2, p = 0.009) tests. Full details of the ICC correction methodology are presented in a companion paper (Dolgikh, in preparation).

### 3.7 Subject 314225: A Boundary Case

Subject 314225 consistently showed the weakest T2 enhancement (ΔG near zero, CI crossing zero at all B). Topological analysis revealed that this subject's T2 subgraph had the lowest edge count (135 vs. 148–158), lowest density (0.71 vs. 0.78–0.83), and lowest average edge connectivity (10.8 vs. 12–14.7). The sparser, less redundant subgraph provides fewer alternative paths for wave-like transport, potentially explaining the reduced ENAQT advantage. This case illustrates the sensitivity of noise-assisted transport to local graph properties and motivates larger cohort studies to characterize individual variation.

---

## 4. Discussion

### 4.1 Summary of Findings

We demonstrate that noise-assisted transport windows — inverted-U profiles of transport efficiency versus dephasing — exist in human DTI connectome subgraphs, are absent under classical diffusion, and are specific to functionally defined neural pathways. The basal ganglia gating loop (T2) shows topology-specific enhancement: native DTI weights outperform strength-matched surrogates, indicating optimization beyond random wiring. The thalamocortical motor relay (T3) shows structure-driven transport where graph geometry dominates and weight perturbations have minimal effect.

### 4.2 Functional Dissociation: Gating vs. Relay

The T2/T3 dissociation maps onto a known functional distinction. Basal ganglia circuits perform action selection and cognitive gating — operations requiring selective signal amplification among competing alternatives [10]. A noise-assisted transport window provides a natural mechanism: intermediate dephasing suppresses off-target interference while maintaining constructive routing to selected targets. Fine-tuned edge weights optimize this balance, explaining why surrogates (which disrupt specific weight configurations) degrade T2 performance.

Thalamocortical motor relay (T3), by contrast, must broadcast signals reliably to multiple cortical targets (precentral + paracentral). Structural redundancy — multiple parallel projections — ensures robust delivery regardless of fine weight tuning. This is consistent with the "structure-driven" label: graph geometry provides sufficient transport capacity, and weight perturbations within strength-matched ensembles produce negligible change.

### 4.3 Topological Mechanisms: Skeleton and Tuning

The B-dependent M2 signal reveals a hierarchical organization. At the coarse level, routing backbone (bridges, hub connectivity) creates the transport skeleton. At the fine level, edge weight distribution tunes conductance within the skeleton. Native DTI weights appear optimized at both levels for T2 but only at the coarse level for T3.

The unexpected M4 (BC Gini) signal — lower heterogeneity degrades ENAQT — suggests that efficient noise-assisted transport requires *directed heterogeneity*: concentrated flow through optimized channels rather than uniform diffusive spread. This aligns with the theoretical prediction that ENAQT exploits constructive interference along specific paths, which requires asymmetric flow distribution.

### 4.4 Relationship to CRN Theory

These findings are consistent with the CRN hypothesis [9] at the architectural level: functional neural pathways support noise-assisted transport windows that could reduce the energetic cost of hypothesis exploration. The T2/T3 dissociation illustrates CRN's prediction that gating circuits (where selective routing is computationally critical) should exhibit stronger noise-assisted optimization than relay circuits (where reliable broadcast suffices).

However, two CRN predictions were *not* confirmed in human data. First, Disorder-Enhanced Selectivity (DES) — where moderate site-energy disorder improves routing selectivity — was not observed at the macroscopic DTI scale (disorder monotonically degraded transport in preliminary analyses; not reported in detail). This likely reflects scale dependence: DES requires fine interferometric structure that exists in micro-connectomes (C. elegans, Drosophila) but may be washed out by the coarse parcellation of DTI tractography. Second, the "containment hypothesis" (that native weights might actively suppress ENAQT in certain pathways) was not supported: T3's structure-driven pattern reflects weight *neutrality*, not active suppression.

### 4.5 Methodological Contributions

Beyond the substantive findings, this study contributes three methodological insights applicable to any computational connectome analysis:

**Pseudoreplication via shared random seeds.** When simulation designs nest random factors (subjects → disorder seeds → surrogates), sharing seeds across subjects creates intraclass correlation that inflates effective sample size. Our hashing solution — deterministically deriving subject-specific seeds — eliminates this artifact without requiring complex mixed-effects modeling. We recommend this approach as standard practice for nested connectome simulations.

**Parameter grid resolution.** Inverted-U profiles can be missed by coarse parameter grids. T3 required 7 κ-points to resolve peaks at intermediate dephasing; initial 5-point grids produced false negatives in 2/5 subjects. We recommend ≥ 15–20 log-spaced dephasing values spanning 3–4 orders of magnitude for any ENAQT-type study.

**Surrogate design hierarchy.** The betweenness-freeze sweep (B = 0, 5, 10, 20%) provides a principled way to separate topology (structure) from tuning (conductance). Unconstrained surrogates (B = 0) test global weight optimization; constrained surrogates (B ≥ 5) test fine-grained tuning on a fixed skeleton. This approach could be applied to any network analysis where structure and weight contributions must be disentangled.

### 4.6 Limitations

**Sample size.** N = 5 subjects provides proof of principle but not population-level inference. Effect sizes (ΔG ≈ 0.5–1.3) suggest N = 15–20 for robust individual-level inference.

**DTI tractography.** Streamline counts are a noisy proxy for axonal connectivity. False positive and false negative connections affect both originals and surrogates identically, so systematic bias toward enhancement is unlikely, but cannot be fully excluded. Validation with probabilistic tractography and cross-species invasive tracers (macaque, rodent) is warranted.

**GKSL as proxy.** No claim is made that neural tissue implements quantum density matrices. GKSL provides a controlled interpolation between wave-like and diffusive limits; the signatures reported here (inverted-U, topology dependence, pathway specificity) would be expected under any wave proxy with tunable damping. Demonstrating these signatures under alternative models (e.g., coupled oscillators, cable equations with damping) is a target for future work.

**Subgraph scale.** K = 15 intermediate nodes yields subgraphs of 20–26 nodes — large enough for meaningful transport but far smaller than full-brain networks. Scale dependence of ENAQT windows (whether they persist, strengthen, or weaken at larger K) remains an open question.

**Single metric.** G = R₀_max(κ*) − R₀_edge is a single-point estimate. Richer characterization (area under inverted-U, κ* distribution, timescale dependence) may reveal additional structure.

### 4.7 Future Directions

**Energy efficiency analysis.** A direct comparison of GKSL transport efficiency (transport per unit dissipation) against classical baselines would test CRN's core claim that wave-like dynamics reduce energetic cost. Preliminary results from a companion energy analysis on 7 subjects (Dolgikh, in preparation) are promising but require full statistical treatment.

**Expanded cohort.** N = 20–100 HCP subjects would enable subgroup analyses (age, sex, handedness) and test whether subject 314225's boundary-case status reflects DTI quality variation or genuine biological heterogeneity.

**Multi-port survey.** Testing additional pathways — sensory association (visual → parietal), default-mode hubs (PCC → lateral PFC) — would characterize which brain systems are "ENAQT-optimized" versus "structure-dominant," mapping the computational landscape of noise-assisted transport across the connectome.

**Falsifiable clinical predictions.** CRN predicts that (a) general anesthesia should shift κ* toward higher values and reduce ENAQT advantage; (b) Parkinson's disease (basal ganglia degeneration) should selectively degrade T2 enhancement while sparing T3. Both predictions are testable with existing clinical neuroimaging datasets.

---

## 5. Conclusions

Noise-assisted transport windows are a reproducible, pathway-specific property of human structural connectomes. The basal ganglia gating loop exhibits topology-specific enhancement, with native DTI weights optimized beyond what random wiring predicts. The thalamocortical motor relay shows structure-driven transport where graph geometry suffices. These findings, combined with classical CTRW controls, negative-control random pairs, and hashed-seed ICC correction, provide robust evidence that the inverted-U transport signature observed in small biological networks extends to human brain connectivity data. The results are consistent with the CRN hypothesis that neural circuits exploit transient wave-like dynamics for efficient signal routing, though direct evidence for the underlying biophysical mechanism remains a target for future experimental work.

---

## Figures

**Figure 1. Noise-assisted transport windows in human connectomes.** R₀_max(κ) curves for five HCP subjects in T2 (basal ganglia, left panel) and T3 (motor relay, right panel) under GKSL dynamics. All subjects show inverted-U profiles with peaks at intermediate dephasing (κ* ≈ 0.25–1.0). Shaded bands: ±1 SD across 12 hashed disorder seeds. Dashed line: R₀_edge (transport at coherent/diffusive limits).

**Figure 2. Classical transport produces no inverted-U.** Left: GKSL R₀_max(κ) showing inverted-U (representative subject 314225, T2). Right: CTRW R₀_max(κ) on same subgraph — monotonically increasing, G_CTRW = 0. Inset: all 64 CTRW cases (5 subjects × 2 ports × {orig + surrogates} × {B = 0, 10}) confirm zero ENAQT gain under classical dynamics.

**Figure 3. Pathway-specific enhancement.** (A) T2: boxplots of ΔG (surr − orig) across 12 seeds × 5 subjects, stratified by B-freeze level. Negative ΔG = original better. Enhancement robust at all B. (B) T3: same layout. ΔG distributed around zero; no consistent enhancement. (C) T0 (random pairs): G near zero (fraction positive ≈ 28%). Star: p < 0.01 (subject-level test).

**Figure 4. Topological predictors of transport loss.** (A) Spearman ρ(ΔG, ΔM) for metrics M1–M8 across 160 T2 surrogates. M2 (effective resistance, ρ = +0.41) and M5 (target betweenness, ρ = +0.28) are strongest positive predictors. M4 (BC Gini, ρ = −0.32) shows unexpected negative signal. (B) Heatmap of ρ stratified by B-freeze level. M2 signal emerges at B ≥ 5 (routing backbone frozen); M5 and M4 stable across B. (C) Scatter: ΔG vs. ΔM2 at B = 0 (no signal) and B = 10 (ρ = 0.64), illustrating the two-level architecture.

---

## Supporting Information

**S1 Table.** Per-subject ENAQT gain (G) for T2 and T3, hashed seeds, all B-levels.

**S2 Table.** CTRW summary: all 64 cases, G = 0, monotonicity confirmed.

**S3 Table.** Surrogate metrics (M1–M8) for 160 T2 surrogates with ΔG.

**S4 Table.** Edge weight–betweenness correlation (Spearman ρ) per subject and port.

**S5 Table.** B-confounder diagnostic: variance ratios and Levene tests for ΔM2 and ΔG.

**S6 Table.** T0 negative control: per-subject G for 10 random pairs.

**S7 Table.** Per-subject T2 topological features (n_edges, density, bc_gini, avg_edge_connectivity).

**S1 Fig.** κ-grid probe bias: 5-point vs. 7-point grid comparison for T3.

**S2 Fig.** ICC shared vs. hashed: G distributions with confidence bands.

**S3 Fig.** CTRW R₀(κ) curves for all 12 plots (3 remaining subjects × 2 ports × 2 B).

**S4 Fig.** Scatter plots ΔG vs. ΔM for M1, M5, M4 (complement to Fig. 4C for M2).

**S5 Fig.** T2 delta boxplots by subject (per-subject CI bands, enhancement labels).

---

## References

1. Attwell D, Laughlin SB. An energy budget for signaling in the grey matter of the brain. J Cereb Blood Flow Metab. 2001;21(10):1133–45.
2. Lennie P. The cost of cortical computation. Curr Biol. 2003;13(6):493–7.
3. Laughlin SB, Sejnowski TJ. Communication in neuronal networks. Science. 2003;301(5641):1870–4.
4. Mohseni M, Rebentrost P, Lloyd S, Aspuru-Guzik A. Environment-assisted quantum walks in photosynthetic energy transfer do not require quantum coherence. J Chem Phys. 2008;129(17):174106.
5. Plenio MB, Huelga SF. Dephasing-assisted transport: quantum and classical. New J Phys. 2008;10(11):113019.
6. Rebentrost P, Mohseni M, Kassal I, Lloyd S, Aspuru-Guzik A. Environment-assisted quantum transport. New J Phys. 2009;11(3):033003.
7. Biggerstaff DN, Heilmann R, Zecevik AA, Gräfe M, Broome MA, Fedrizzi A, et al. Enhancing coherent transport in a photonic network using controllable decoherence. Nat Commun. 2016;7:11282.
8. Cao J, Cogdell RJ, Coker DF, Duan H-G, Hauer J, Kleinekathöfer U, et al. Quantum biology revisited. Sci Adv. 2020;6(14):eaaz4888.
9. Dolgikh O. Thermodynamic advantage of transient wave dynamics in hierarchical decision architectures. Zenodo. 2024. doi:10.5281/zenodo.18249250.
10. Redgrave P, Prescott TJ, Gurney K. The basal ganglia: a vertebrate solution to the selection problem? Neuroscience. 1999;89(4):1009–23.
11. Sherman SM, Guillery RW. The role of the thalamus in the flow of information to the cortex. Phil Trans R Soc Lond B. 2002;357(1428):1695–708.
12. Van Essen DC, Smith SM, Barch DM, Behrens TEJ, Yacoub E, Ugurbil K. The WU-Minn Human Connectome Project: an overview. NeuroImage. 2013;80:62–79.
13. Saravanan V, Berman GJ, Sober SJ. Application of the hierarchical bootstrap to multi-level data in neuroscience. Neurons Behav Data Anal Theory. 2020;3(5):1–25.

---

## Acknowledgments

Simulations used AI-assisted code development (Claude, Anthropic; ChatGPT, OpenAI). HCP data provided by the WU-Minn Consortium (Principal Investigators: David Van Essen and Kamil Ugurbil; 1U54MH091657). The author thanks the open-source communities maintaining NumPy, SciPy, NetworkX, and Matplotlib.

---

## Author Contributions

OD: Conceptualization, methodology, software, formal analysis, investigation, data curation, writing — original draft, visualization.

---

## Competing Interests

The author declares no competing interests.

---

## Data Availability Statement

All simulation code, parameter configurations, and output datasets are available at [GitHub repository URL]. Raw HCP connectome data are available through the HCP data use agreement (https://www.humanconnectome.org).
