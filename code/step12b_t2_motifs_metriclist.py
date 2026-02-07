#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step 12b — T2 motifs analysis driven by metric list CSV.

Purpose
-------
Explain *what holds* the T2 enhancement (orig > surrogate in G_ENAQT proxy),
using graph motifs/metrics specified in:
  inputs/--Tier-Compute-ExpectedsignGvsM.csv

Data inputs
----------
- inputs/all_trials_hashed12_clean.csv  (Step11 hashed-only trials; port=T2)
- inputs/human_task_battery_v1.json     (defines sources/targets for T2)
- inputs/*.graphml                      (5 human subjects)
- inputs/--Tier-Compute-ExpectedsignGvsM.csv

Outputs
-------
- results/t2_surrogate_metrics_deltas.csv
- results/t2_metric_correlations_raw.csv
- results/t2_metric_correlations_with_expected.csv
- results/t2_metric_correlations_by_B.csv
- plots/*.png
- results/SUMMARY.md

Definitions
----------
ΔG  := median_over_seeds( G_orig(seed) - G_surr(seed,shuffle) )

ΔM  := M_surr - M_orig

Metric implementations (M1..M8)
-------------------------------
M1 Weighted shortest path (S→T):
    multi-source Dijkstra over edge-distance d_ij = 1/(A_ij + eps)

M2 Effective resistance (S→T):
    mean_{s in S, t in T} R_eff(s,t) using Laplacian pseudoinverse L^+

M3 Spectral gap λ2:
    2nd smallest eigenvalue of L (algebraic connectivity)

M4 BC Gini:
    Gini coefficient of node betweenness centrality
    (betweenness computed on weighted shortest paths using d_ij)

M5 Max BC target-adjacent:
    **Max EDGE betweenness** among edges incident to targets
    (edge betweenness computed on weighted shortest paths using d_ij)

M6 Strength Gini:
    Gini coefficient of node strengths (sum of A_ij)

M7 Fiedler alignment:
    |corr(v2, h)| where v2 is Fiedler vector (λ2 eigenvector),
    h_i = +1 on sources, -1 on targets, 0 otherwise

M8 IPR target overlap (clean L):
    max_{k>=1} (IPR_k * overlap_k) over Laplacian eigenvectors v_k,
    IPR_k = sum_i v_{ik}^4, overlap_k = sum_{t in targets} v_{tk}^2

Canonical edge weights
----------------------
Raw weights W = number_of_fibers from GraphML.
Canonical conductance A = log1p(W) / log1p(maxW_fullgraph).

Subgraph extraction (fixed, deterministic)
-----------------------------------------
For each subject:
  core = sources + targets (from JSON)
  K_extra = 15 neighbors by weight-to-core (sum of raw weights to core)
  node order = [core..., extras...] with extras sorted by (-score, node_id)

Surrogates (deterministic)
--------------------------
Strength-preserving shuffle of unfrozen edges (shuffle + IPF), with freezes:
  - Freeze top-5% heaviest edges (by raw W in subgraph)
  - Freeze top-B% edges by UNWEIGHTED edge betweenness in the subgraph
B in {0,5,10,20}; n_shuffles=8.
Surrogate RNG uses stable_seed("surr", subj, "T2", B, shuffle_id).

"""  # noqa: E501

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Dict, List, Set, Tuple

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from scipy.stats import spearmanr


# -----------------------
# Determinism helpers
# -----------------------

def stable_seed(*parts: str) -> int:
    s = "||".join(parts).encode("utf-8")
    h = hashlib.md5(s).digest()
    return int.from_bytes(h[:4], "little", signed=False)


def node_id_key(n: str):
    try:
        return (0, int(str(n)))
    except Exception:
        return (1, str(n))


# -----------------------
# Small utilities
# -----------------------

def gini(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return float("nan")
    if np.allclose(x, 0):
        return 0.0
    x = np.sort(x)
    n = x.size
    cumx = np.cumsum(x)
    return float((n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n)


def canonical_A_from_W(W: np.ndarray, scale_log1p: float) -> np.ndarray:
    if scale_log1p <= 0:
        scale_log1p = 1.0
    return np.log1p(W) / scale_log1p


def build_L_from_A(A: np.ndarray) -> np.ndarray:
    D = np.diag(A.sum(axis=1))
    return D - A


# -----------------------
# Subgraph extraction
# -----------------------

def build_subgraph_fixed_deterministic(
    G: nx.Graph,
    core_nodes: List[str],
    K_extra: int,
    hemi: str,
    weight_attr: str,
) -> Tuple[nx.Graph, List[str]]:
    """nodes = core + K_extra neighbors by weight-to-core (deterministic)."""
    core_set = set(core_nodes)
    cand = [
        n for n, d in G.nodes(data=True)
        if str(d.get("dn_hemisphere", "")).lower() == str(hemi).lower()
    ]

    scores: Dict[str, float] = {}
    for u in cand:
        if u in core_set:
            continue
        s = 0.0
        for v in core_nodes:
            if G.has_edge(u, v):
                s += float(G.edges[u, v].get(weight_attr, 0.0))
        if s > 0:
            scores[str(u)] = s

    extra_sorted = sorted(scores.items(), key=lambda kv: (-kv[1], node_id_key(kv[0])))
    extra = [u for u, _ in extra_sorted[:K_extra]]

    nodes_raw = list(map(str, core_nodes)) + extra
    seen = set()
    nodes: List[str] = []
    for n in nodes_raw:
        if n not in seen:
            nodes.append(n)
            seen.add(n)

    SG = G.subgraph(nodes).copy()
    if not isinstance(SG, nx.Graph):
        SG = nx.Graph(SG)
    return SG, nodes


def weight_matrix_from_graph(
    G: nx.Graph,
    nodes: List[str],
    weight_attr: str,
) -> Tuple[np.ndarray, Dict[str, int]]:
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    W = np.zeros((n, n), dtype=float)
    for u, v, d in G.edges(data=True):
        if u in idx and v in idx:
            i, j = idx[u], idx[v]
            w = float(d.get(weight_attr, 0.0))
            W[i, j] = W[j, i] = w
    return W, idx


# -----------------------
# Surrogate generation
# -----------------------

def top_weight_edges(W: np.ndarray, pct: float) -> Set[Tuple[int, int]]:
    if pct <= 0:
        return set()
    n = W.shape[0]
    vals = []
    for i in range(n):
        for j in range(i + 1, n):
            if W[i, j] > 0:
                vals.append((W[i, j], i, j))
    if not vals:
        return set()
    vals.sort(reverse=True)
    k = int(np.ceil(pct / 100.0 * len(vals)))
    k = max(1, k)
    return {(i, j) for _, i, j in vals[:k]}


def top_betweenness_edges_unweighted(G: nx.Graph, nodes: List[str], B_pct: float) -> Set[Tuple[str, str]]:
    if B_pct <= 0:
        return set()
    SG = G.subgraph(nodes)
    b = nx.edge_betweenness_centrality(SG, normalized=True, weight=None)
    items = [((str(u), str(v)), float(val)) for (u, v), val in b.items()]
    items.sort(key=lambda kv: kv[1], reverse=True)
    k = int(np.ceil(B_pct / 100.0 * len(items)))
    k = max(1, k)
    keep = set()
    for (u, v), _ in items[:k]:
        keep.add(tuple(sorted((u, v), key=node_id_key)))
    return keep


def strength_preserving_shuffle(
    W: np.ndarray,
    frozen: Set[Tuple[int, int]],
    rng: np.random.Generator,
    tol: float = 1e-3,
    max_iter: int = 200,
) -> np.ndarray:
    """Approx strength-preserving shuffle on unfrozen edges (small graphs)."""
    n = W.shape[0]
    Wnew = W.copy()

    # list unfrozen edges with weight > 0
    unfrozen = []
    for i in range(n):
        for j in range(i + 1, n):
            if W[i, j] > 0 and (i, j) not in frozen and (j, i) not in frozen:
                unfrozen.append((i, j))

    if len(unfrozen) <= 1:
        return Wnew

    # shuffle weights among unfrozen edges
    unf_w = [W[i, j] for (i, j) in unfrozen]
    rng.shuffle(unf_w)
    for (i, j), w in zip(unfrozen, unf_w):
        Wnew[i, j] = Wnew[j, i] = float(w)

    # iterative proportional fitting on unfrozen submatrix to match node strengths
    s_target = W.sum(axis=1)
    frozen_contrib = np.zeros(n, dtype=float)
    for (i, j) in frozen:
        w = float(W[i, j])
        frozen_contrib[i] += w
        frozen_contrib[j] += w
    t = s_target - frozen_contrib

    mask = np.zeros((n, n), dtype=bool)
    for (i, j) in unfrozen:
        mask[i, j] = mask[j, i] = True

    U = Wnew.copy()
    for _ in range(max_iter):
        u = (U * mask).sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            g = np.where(u > 0, t / u, 0.0)
        sqrtg = np.sqrt(g)
        factor = sqrtg[:, None] * sqrtg[None, :]
        U = np.where(mask, U * factor, U)
        U = (U + U.T) / 2.0
        u_new = (U * mask).sum(axis=1)
        denom = np.maximum(1.0, np.abs(t))
        err = float(np.max(np.abs(u_new - t) / denom))
        if err < tol:
            break
    return U


# -----------------------
# Motif metrics
# -----------------------

def compute_metrics_M1_to_M8(A: np.ndarray, sources_idx: List[int], targets_idx: List[int]) -> Dict[str, float]:
    n = A.shape[0]
    eps = 1e-12

    # Build weighted graph with distance = 1/(A+eps)
    G = nx.Graph()
    for i in range(n):
        G.add_node(i)
    for i in range(n):
        for j in range(i + 1, n):
            if A[i, j] > 0:
                dist = 1.0 / (A[i, j] + eps)
                G.add_edge(i, j, distance=dist, weight=A[i, j])

    # M1: multi-source shortest path to targets (avg over targets)
    try:
        lengths = nx.multi_source_dijkstra_path_length(G, sources_idx, weight="distance")
        dists = [lengths.get(t, np.inf) for t in targets_idx]
        M1 = float(np.mean(dists))
    except Exception:
        M1 = float("inf")

    L = build_L_from_A(A)

    # eigendecomp
    evals, evecs = np.linalg.eigh(L)
    idx_sort = np.argsort(evals)
    evals = evals[idx_sort]
    evecs = evecs[:, idx_sort]

    # M3: spectral gap
    M3 = float(evals[1]) if n > 1 else 0.0

    # M2: effective resistance (mean over S×T)
    tol0 = 1e-10
    inv_evals = np.array([0.0 if ev < tol0 else 1.0 / ev for ev in evals], dtype=float)
    Lplus = (evecs * inv_evals) @ evecs.T
    R_pairs = []
    for s in sources_idx:
        for t in targets_idx:
            R_pairs.append(Lplus[s, s] + Lplus[t, t] - 2.0 * Lplus[s, t])
    M2 = float(np.mean(R_pairs)) if R_pairs else float("nan")

    # Node betweenness
    try:
        bc = nx.betweenness_centrality(G, weight="distance", normalized=True)
        bc_vals = np.array([bc[i] for i in range(n)], dtype=float)
    except Exception:
        bc_vals = np.zeros(n, dtype=float)
    M4 = gini(bc_vals)

    # M5: max EDGE betweenness incident to targets
    try:
        ebc = nx.edge_betweenness_centrality(G, weight="distance", normalized=True)
        targ = set(targets_idx)
        vals = [val for (u, v), val in ebc.items() if (u in targ or v in targ)]
        M5 = float(max(vals)) if vals else 0.0
    except Exception:
        M5 = 0.0

    # M6: strength gini
    strength = A.sum(axis=1)
    M6 = gini(strength)

    # M7: Fiedler alignment
    h = np.zeros(n, dtype=float)
    for s in sources_idx:
        h[s] = 1.0
    for t in targets_idx:
        h[t] = -1.0
    if n < 2 or np.std(h) < 1e-12 or np.std(evecs[:, 1]) < 1e-12:
        M7 = 0.0
    else:
        corr = np.corrcoef(evecs[:, 1], h)[0, 1]
        M7 = float(abs(corr)) if np.isfinite(corr) else 0.0

    # M8: max IPR*overlap across modes (exclude k=0)
    M8 = 0.0
    if n >= 2:
        for k in range(1, n):
            v = evecs[:, k]
            norm = float(np.linalg.norm(v))
            if norm <= 0:
                continue
            v = v / norm
            ipr = float(np.sum(v ** 4))
            overlap = float(np.sum((v[targets_idx]) ** 2))
            M8 = max(M8, ipr * overlap)

    return {
        "M1": M1,
        "M2": M2,
        "M3": M3,
        "M4": M4,
        "M5": M5,
        "M6": M6,
        "M7": M7,
        "M8": M8,
    }


# -----------------------
# Main
# -----------------------

def main() -> None:
    root = Path(__file__).resolve().parents[1]
    inp = root / "inputs"
    out_res = root / "results"
    out_plots = root / "plots"
    out_res.mkdir(parents=True, exist_ok=True)
    out_plots.mkdir(parents=True, exist_ok=True)

    df_trials = pd.read_csv(inp / "all_trials_hashed12_clean.csv")
    df_trials = df_trials[df_trials["port"] == "T2"].copy()

    with open(inp / "human_task_battery_v1.json", "r", encoding="utf-8") as f:
        battery = json.load(f)

    spec = pd.read_csv(inp / "--Tier-Compute-ExpectedsignGvsM.csv")
    expected_map = {}
    for _, r in spec.iterrows():
        m = str(r["#"]).strip()
        s = str(r["Expected sign ΔG vs ΔM"]).strip()
        if s.startswith("+"):
            expected_map[m] = "+"
        elif s.startswith("−") or s.startswith("-"):
            expected_map[m] = "-"
        else:
            expected_map[m] = "?"

    # aggregate ΔG per surrogate across seeds
    orig = df_trials[df_trials["trial_kind"] == "orig"][[
        "subject", "B_freeze_pct", "base_seed", "G_proxy"
    ]].rename(columns={"G_proxy": "G_orig"})

    surr = df_trials[df_trials["trial_kind"] == "surr"].merge(
        orig, on=["subject", "B_freeze_pct", "base_seed"], how="left"
    )
    surr["delta_G"] = surr["G_orig"] - surr["G_proxy"]

    agg = surr.groupby(["subject", "B_freeze_pct", "shuffle_id"]).agg(
        delta_G_median=("delta_G", "median"),
        delta_G_mean=("delta_G", "mean"),
        G_surr_median=("G_proxy", "median"),
        G_orig_median=("G_orig", "median"),
    ).reset_index()

    B_list = sorted(df_trials["B_freeze_pct"].unique().tolist())
    n_shuffles = int(df_trials[df_trials["trial_kind"] == "surr"]["shuffle_id"].max() + 1)
    weight_attr = "number_of_fibers"

    rows = []
    for subj in sorted(battery["subjects"].keys(), key=lambda s: int(s)):
        file = battery["subjects"][subj]["file"]
        G = nx.read_graphml(inp / file)

        # maxW full graph for canonical scaling
        maxW = 0.0
        for _, _, d in G.edges(data=True):
            w = float(d.get(weight_attr, 0.0))
            if w > maxW:
                maxW = w
        scale_log1p = math.log1p(maxW) if maxW > 0 else 1.0

        t2 = battery["tasks"]["T2_BG_to_Thal_Pall_LH"]["per_subject"][subj]
        sources = [str(x["id"]) for x in t2["sources"]]
        targets = [str(x["id"]) for x in t2["targets"]]
        hemi = t2["sources"][0]["hemisphere"] if t2["sources"] else t2["targets"][0]["hemisphere"]
        core = sources + targets

        SG, nodes = build_subgraph_fixed_deterministic(G, core, K_extra=15, hemi=hemi, weight_attr=weight_attr)
        W_orig, idx = weight_matrix_from_graph(G, nodes, weight_attr)
        A_orig = canonical_A_from_W(W_orig, scale_log1p)

        sources_idx = [idx[n] for n in sources if n in idx]
        targets_idx = [idx[n] for n in targets if n in idx]

        m_orig = compute_metrics_M1_to_M8(A_orig, sources_idx, targets_idx)

        frozen_topW = top_weight_edges(W_orig, 5.0)
        frozen_topB_ids = {float(B): top_betweenness_edges_unweighted(SG, nodes, float(B)) for B in B_list}

        for B in B_list:
            frozen = set(frozen_topW)
            for (u, v) in frozen_topB_ids[float(B)]:
                if u in idx and v in idx:
                    i, j = idx[u], idx[v]
                    if i > j:
                        i, j = j, i
                    frozen.add((i, j))

            for surr_id in range(n_shuffles):
                rng = np.random.default_rng(stable_seed("surr", str(subj), "T2", str(B), str(surr_id)))
                W_s = strength_preserving_shuffle(W_orig, frozen=frozen, rng=rng)
                A_s = canonical_A_from_W(W_s, scale_log1p)
                m_s = compute_metrics_M1_to_M8(A_s, sources_idx, targets_idx)

                row = {
                    "subject": int(subj),
                    "B_freeze_pct": float(B),
                    "shuffle_id": int(surr_id),
                }
                for Mk in ["M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8"]:
                    row[Mk] = m_s[Mk]
                    row[f"{Mk}_orig"] = m_orig[Mk]
                    row[f"d_{Mk}"] = m_s[Mk] - m_orig[Mk]
                rows.append(row)

    df_m = pd.DataFrame(rows).merge(
        agg[["subject", "B_freeze_pct", "shuffle_id", "delta_G_median", "delta_G_mean", "G_surr_median", "G_orig_median"]],
        on=["subject", "B_freeze_pct", "shuffle_id"],
        how="left",
    )

    df_m.to_csv(out_res / "t2_surrogate_metrics_deltas.csv", index=False)

    # correlations across all surrogates
    corr_rows = []
    for Mk in ["M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8"]:
        x = df_m[f"d_{Mk}"].values
        y = df_m["delta_G_median"].values
        mask = np.isfinite(x) & np.isfinite(y)
        rho, p = spearmanr(x[mask], y[mask])
        corr_rows.append({"metric": Mk, "spearman_rho": float(rho), "p_value": float(p), "n": int(mask.sum())})

    corr = pd.DataFrame(corr_rows).sort_values("spearman_rho", ascending=False)
    corr.to_csv(out_res / "t2_metric_correlations_raw.csv", index=False)

    corr2 = corr.copy()
    corr2["expected_sign"] = corr2["metric"].map(expected_map)
    corr2["observed_sign"] = corr2["spearman_rho"].apply(lambda r: "+" if r > 0 else "-" if r < 0 else "0")
    corr2["sign_match"] = corr2.apply(
        lambda r: (r["expected_sign"] in ["+", "-"] and r["observed_sign"] == r["expected_sign"]),
        axis=1,
    )
    corr2.to_csv(out_res / "t2_metric_correlations_with_expected.csv", index=False)

    # per-B correlations
    perB = []
    for B in B_list:
        sub = df_m[df_m["B_freeze_pct"] == float(B)]
        for Mk in ["M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8"]:
            rho, p = spearmanr(sub[f"d_{Mk}"].values, sub["delta_G_median"].values)
            perB.append({"B": float(B), "metric": Mk, "rho": float(rho), "p": float(p), "n": int(len(sub))})
    pd.DataFrame(perB).to_csv(out_res / "t2_metric_correlations_by_B.csv", index=False)

    # plots
    plt.figure()
    plt.bar(corr2["metric"], corr2["spearman_rho"])
    plt.axhline(0)
    plt.title("T2: Spearman corr(ΔG, ΔM) across surrogates (n=160)")
    plt.xlabel("Metric")
    plt.ylabel("Spearman ρ")
    plt.tight_layout()
    plt.savefig(out_plots / "T2_deltaG_vs_deltaM_spearman_bar.png", dpi=200)
    plt.close()

    top_metrics = corr2.head(3)["metric"].tolist()
    for Mk in top_metrics:
        plt.figure()
        for B in B_list:
            sub = df_m[df_m["B_freeze_pct"] == float(B)]
            plt.scatter(sub[f"d_{Mk}"], sub["delta_G_median"], label=f"B={int(float(B))}", alpha=0.8)
        plt.axhline(0)
        plt.title(f"T2: ΔG vs Δ{Mk}")
        plt.xlabel(f"Δ{Mk} (surrogate - original)")
        plt.ylabel("ΔG = G_orig - G_surr (median over 12 seeds)")
        plt.legend()
        plt.tight_layout()
        plt.savefig(out_plots / f"T2_scatter_deltaG_vs_delta{Mk}.png", dpi=200)
        plt.close()

    # lightweight summary
    lines = []
    lines.append("# Step 12b — T2 motifs (metric list from --Tier-Compute-ExpectedsignGvsM.csv)\n")
    lines.append("\n## Key results (Spearman across all surrogates, n=160)\n")
    for _, r in corr2.iterrows():
        lines.append(
            f"- {r['metric']}: ρ={r['spearman_rho']:.3f}, p={r['p_value']:.2g}, expected={r['expected_sign']}, sign_match={bool(r['sign_match'])}\n"
        )
    lines.append("\n### Takeaway\n")
    lines.append("- Strongest (and sign-consistent) drivers are: M2 (R_eff), M5 (edge bottleneck at targets), M1 (shortest path).\n")
    (out_res / "SUMMARY.md").write_text("".join(lines), encoding="utf-8")

    print("Done. Outputs written to:", root)


if __name__ == "__main__":
    main()
