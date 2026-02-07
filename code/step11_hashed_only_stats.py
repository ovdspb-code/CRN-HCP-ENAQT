#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 11 — Final hashed-only statistics (12 hashed seeds) on CanonV2 / Canonical-H
for Human DTI connectomes (5 subjects), ports T2 + T3, B-sweep (B=0/5/10/20).

What this step does
-------------------
1) Uses ONLY hashed (subject-specific) disorder seeds (no shared seeds in N_eff).
   base_seeds = 20250212..20250223 (12 seeds)
   seed_eff = int(sha256(f"{seed}_{subj}").hexdigest()[:8], 16)

2) Uses Canonical-H consistently:
   A = log1p(W) / log1p(maxW_full)
   L = D - A
   H = J*(L + diag(xi_sub))

3) Uses GKSL dynamics with:
   - pure dephasing (kappa)
   - target sink gamma_T on target nodes
   - leakage sink gamma_L on non-target nodes
   Integrator: batched RK2 midpoint, dt=0.05, T_end=10.0, burn_in=0.5

4) Surrogates (controls):
   Strength-preserving weight shuffle with:
     - Freeze top-5% heaviest edges by raw weight
     - Freeze top-B% edges by unweighted edge betweenness (B in {0,5,10,20})
   n_shuffles=8 per (subject,port,B), reused across all disorder seeds.

Kappa grids
-----------
T2 uses 5-point probes (bias=0 confirmed earlier):
  [0.001, 0.251, 0.631, 1.0, 10.0]
T3 uses 7-point probes (fixes known probe bias for small-kappa peaks):
  [0.001, 0.010, 0.063, 0.251, 0.631, 1.0, 10.0]

Outputs
-------
results/all_trials_hashed12.csv
results/per_seed_summary.csv
results/per_subject_labels.csv
results/group_summary.csv
plots/*.png
SUMMARY.md

Reproducibility / determinism
-----------------------------
- node order: sorted(node_id) everywhere
- K-extra selection: weight-to-core then node_id tie-break
- surrogate RNG: stable_seed("surr", subj, port, B, shuffle_id)
- disorder mapping: generate xi on full graph ordered by sorted(node_id), then map by node_id into subgraph
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Dict, List, Set, Tuple

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


# -----------------------
# Determinism helpers
# -----------------------

def stable_seed(*parts: str) -> int:
    s = "||".join(parts).encode("utf-8")
    h = hashlib.md5(s).digest()
    return int.from_bytes(h[:4], "little", signed=False)


def node_id_key(n: str):
    try:
        return (0, int(n))
    except Exception:
        return (1, str(n))


def sha256_seed_int32(base_seed: int, subject_id: str) -> int:
    """EXACT spec: int(sha256(f"{seed}_{subj}").hexdigest()[:8], 16)"""
    h = hashlib.sha256(f"{base_seed}_{subject_id}".encode("utf-8")).hexdigest()
    return int(h[:8], 16)


# -----------------------
# Task battery I/O
# -----------------------

def get_port_nodes_from_battery(battery: dict, subj: str, port: str) -> Tuple[List[str], List[str], str]:
    port_key_map = {
        "T2": "T2_BG_to_Thal_Pall_LH",
        "T3": "T3_Thal_to_Motor_LH",
        "T0": "T0_LOF_to_SF_RH",
        "T1": "T1_Sensory_to_Motor_LH",
    }
    key = port_key_map[port]
    entry = battery["tasks"][key]["per_subject"][subj]

    sources = [str(x["id"]) for x in entry.get("sources", [])]
    targets = [str(x["id"]) for x in entry.get("targets", [])]

    hemi = None
    if entry.get("sources"):
        hemi = entry["sources"][0].get("hemisphere")
    elif entry.get("targets"):
        hemi = entry["targets"][0].get("hemisphere")
    hemi = str(hemi) if hemi is not None else "left"
    return sources, targets, hemi


# -----------------------
# Subgraph construction
# -----------------------

def build_subgraph_fixed_deterministic(
    G: nx.Graph, core_nodes: List[str], K_extra: int, hemi: str, weight_attr: str
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
            scores[u] = s

    extra_sorted = sorted(scores.items(), key=lambda kv: (-kv[1], node_id_key(kv[0])))
    extra = [u for u, _ in extra_sorted[:K_extra]]

    nodes_raw = list(core_nodes) + extra
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


# -----------------------
# Matrices
# -----------------------

def weight_matrix_from_graph(G: nx.Graph, nodes: List[str], weight_attr: str) -> Tuple[np.ndarray, Dict[str, int]]:
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    W = np.zeros((n, n), dtype=float)
    for u, v, d in G.edges(data=True):
        if u in idx and v in idx:
            i, j = idx[u], idx[v]
            w = float(d.get(weight_attr, 0.0))
            W[i, j] = W[j, i] = w
    return W, idx


def build_A_normalized(W: np.ndarray, scale_log1p: float) -> np.ndarray:
    return np.log1p(W) / scale_log1p


def build_L_from_A(A: np.ndarray) -> np.ndarray:
    D = np.diag(A.sum(axis=1))
    return D - A


def build_H_from_L_xi(L: np.ndarray, xi_sub: np.ndarray, J: float) -> np.ndarray:
    return J * (L + np.diag(xi_sub))


# -----------------------
# Surrogates: strength-preserving + freezes
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
# GKSL simulator (batched RK2)
# -----------------------

def simulate_R0max_kappa_batch_RK2(
    H: np.ndarray,
    sources_idx: List[int],
    targets_idx: List[int],
    kappa_grid: np.ndarray,
    gamma_T: float,
    gamma_L: float,
    dt: float,
    T_end: float,
    burn_in: float,
    eps_stab: float,
) -> np.ndarray:
    """Batched RK2 midpoint for GKSL with pure dephasing + sinks."""
    Hc = H.astype(np.complex128)
    n = Hc.shape[0]
    K = int(len(kappa_grid))
    kappa_arr = np.array(kappa_grid, dtype=float)

    gamma = np.full(n, gamma_L, dtype=float)
    gamma[targets_idx] = gamma_T

    Gmat = 0.5 * (gamma[:, None] + gamma[None, :]).astype(np.complex128)

    targets_mask = np.zeros(n, dtype=bool)
    targets_mask[targets_idx] = True

    steps = int(T_end / dt)
    burn_steps = int(burn_in / dt)

    rho = np.zeros((K, n, n), dtype=np.complex128)
    for s in sources_idx:
        rho[:, s, s] += 1.0 / max(1, len(sources_idx))

    P_T = np.zeros(K, dtype=float)
    P_D = np.zeros(K, dtype=float)
    R0max = np.zeros(K, dtype=float)

    diag_idx = np.arange(n)
    offmask = (~np.eye(n, dtype=bool)).astype(np.complex128)
    kappa_b = kappa_arr[:, None, None].astype(np.float64)

    def drho(r: np.ndarray) -> np.ndarray:
        d = -1j * (Hc @ r - r @ Hc)
        d -= kappa_b * (r * offmask[None, :, :])
        d -= Gmat[None, :, :] * r
        return d

    for step in range(steps):
        k1 = drho(rho)
        rho_mid = rho + 0.5 * dt * k1
        k2 = drho(rho_mid)
        rho = rho + dt * k2

        pop = np.real(rho[:, diag_idx, diag_idx])
        pop = np.clip(pop, 0, None)

        P_T += gamma_T * pop[:, targets_mask].sum(axis=1) * dt
        P_D += gamma_L * pop[:, ~targets_mask].sum(axis=1) * dt

        if step >= burn_steps:
            R0 = P_T / (P_D + eps_stab)
            R0max = np.maximum(R0max, R0)

    return R0max


def proxy_from_R0max(R0max: np.ndarray, kappa_probes: np.ndarray) -> Dict[str, float]:
    R_edge = float(max(R0max[0], R0max[-1]))
    mid = R0max[1:-1]
    mid_k = kappa_probes[1:-1]
    j = int(np.argmax(mid))
    R_peak = float(mid[j])
    k_star = float(mid_k[j])
    G = float(R_peak - R_edge)
    return {
        "R_edge": R_edge,
        "R_peak": R_peak,
        "kappa_star": k_star,
        "G_proxy": G,
    }


# -----------------------
# Hierarchical bootstrap for per-subject labels
# -----------------------

def hierarchical_bootstrap_delta(
    deltas_by_seed: Dict[int, np.ndarray],
    stat=np.median,
    n_boot: int = 3000,
    seed: int = 0,
) -> Tuple[float, float]:
    """
    deltas_by_seed: {seed_eff: array of deltas (shuffle-level, e.g. G_surr - G_orig)}
    Outer bootstrap samples seeds with replacement.
    Inner bootstrap samples shuffles within each selected seed with replacement,
    computes seed-level statistic, then aggregates across seeds.
    """
    rng = np.random.default_rng(seed)
    seeds = np.array(list(deltas_by_seed.keys()), dtype=int)
    if len(seeds) == 0:
        return (float("nan"), float("nan"))

    boot_vals = []
    for _ in range(n_boot):
        sampled_seeds = rng.choice(seeds, size=len(seeds), replace=True)
        seed_stats = []
        for s in sampled_seeds:
            arr = np.array(deltas_by_seed[int(s)], dtype=float)
            if len(arr) == 0:
                continue
            arr_s = rng.choice(arr, size=len(arr), replace=True)
            seed_stats.append(stat(arr_s))
        if len(seed_stats) == 0:
            boot_vals.append(np.nan)
        else:
            boot_vals.append(stat(np.array(seed_stats, dtype=float)))
    boot_vals = np.array(boot_vals, dtype=float)
    boot_vals = boot_vals[np.isfinite(boot_vals)]
    if len(boot_vals) == 0:
        return (float("nan"), float("nan"))
    return (float(np.percentile(boot_vals, 2.5)), float(np.percentile(boot_vals, 97.5)))


def label_from_ci(ci_low: float, ci_high: float, eps: float = 1e-9) -> str:
    # delta defined as median(G_surr) - G_orig:
    #  >0  => surrogates better (containment-ish)
    #  <0  => orig better (enhancement)
    if not np.isfinite(ci_low) or not np.isfinite(ci_high):
        return "ambiguous"
    if ci_low > eps:
        return "containment_or_surr_better"
    if ci_high < -eps:
        return "enhancement_or_orig_better"
    return "structure-driven_or_ambiguous"


# -----------------------
# Main
# -----------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", type=str, required=True, help="Folder with graphml + human_task_battery_v1.json")
    ap.add_argument("--outdir", type=str, required=True)
    ap.add_argument("--K_extra", type=int, default=15)
    ap.add_argument("--B_list", type=str, default="0,5,10,20")
    ap.add_argument("--n_shuffles", type=int, default=8)
    ap.add_argument("--eps_dis", type=float, default=3.0)
    ap.add_argument("--J", type=float, default=0.1)
    ap.add_argument("--gamma_T", type=float, default=1.0)
    ap.add_argument("--gamma_L", type=float, default=5e-4)
    ap.add_argument("--dt", type=float, default=0.05)
    ap.add_argument("--T_end", type=float, default=10.0)
    ap.add_argument("--burn_in", type=float, default=0.5)
    ap.add_argument("--eps_stab", type=float, default=0.001)
    args = ap.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "results").mkdir(exist_ok=True)
    (outdir / "plots").mkdir(exist_ok=True)

    battery = json.loads((indir / "human_task_battery_v1.json").read_text(encoding="utf-8"))

    subjects = ["100206", "304727", "314225", "333330", "519647"]
    ports = ["T2", "T3"]
    B_list = [float(x) for x in args.B_list.split(",") if x.strip()]

    # Hashed-only 12 seeds
    base_seeds = list(range(20250212, 20250224))

    kappa_T2 = np.array([0.001, 0.251, 0.631, 1.0, 10.0], dtype=float)
    kappa_T3 = np.array([0.001, 0.010, 0.063, 0.251, 0.631, 1.0, 10.0], dtype=float)

    all_rows = []

    # Cache for subgraphs and surrogates per (subj,port,B)
    cache = {}

    for subj in subjects:
        G = nx.read_graphml(indir / f"{subj}_repeated10_scale125.graphml")

        # Detect weight attr
        sample_edge = next(iter(G.edges(data=True)))[2]
        weight_attr = "number_of_fibers" if "number_of_fibers" in sample_edge else ("weight" if "weight" in sample_edge else list(sample_edge.keys())[0])

        full_nodes_sorted = sorted(G.nodes(), key=node_id_key)
        maxW = max(float(d.get(weight_attr, 0.0)) for _, _, d in G.edges(data=True))
        scale_log1p = float(np.log1p(maxW))

        for port in ports:
            sources, targets, hemi = get_port_nodes_from_battery(battery, subj, port)
            core = sources + targets

            SG, nodes = build_subgraph_fixed_deterministic(G, core, args.K_extra, hemi, weight_attr)
            W_orig, idx = weight_matrix_from_graph(SG, nodes, weight_attr)

            # Map source/target indices inside subgraph
            sources_idx = [idx[n] for n in sources if n in idx]
            targets_idx = [idx[n] for n in targets if n in idx]

            # Determine kappa probes for this port
            if port == "T2":
                kappa_probes = kappa_T2
            else:
                kappa_probes = kappa_T3

            # Precompute frozen edges sets (by index)
            frozen_topW = top_weight_edges(W_orig, 5.0)
            frozen_topB_ids = {b: top_betweenness_edges_unweighted(SG, nodes, b) for b in B_list}

            for B in B_list:
                key = (subj, port, B)
                if key not in cache:
                    # convert frozen_topB_ids (node pairs) to index pairs
                    frozen = set(frozen_topW)
                    for (u, v) in frozen_topB_ids[B]:
                        if u in idx and v in idx:
                            i, j = idx[u], idx[v]
                            if i > j:
                                i, j = j, i
                            frozen.add((i, j))

                    # build surrogate weight matrices (raw W) deterministically
                    surrogates = []
                    for surr_id in range(args.n_shuffles):
                        rng = np.random.default_rng(stable_seed("surr", subj, port, str(B), str(surr_id)))
                        W_s = strength_preserving_shuffle(W_orig, frozen=frozen, rng=rng)
                        surrogates.append(W_s)
                    cache[key] = {
                        "nodes": nodes,
                        "idx": idx,
                        "W_orig": W_orig,
                        "W_surr": surrogates,
                        "sources_idx": sources_idx,
                        "targets_idx": targets_idx,
                        "scale_log1p": scale_log1p,
                        "full_nodes_sorted": full_nodes_sorted,
                        "kappa_probes": kappa_probes,
                    }

                c = cache[key]

                # Precompute L for original and each surrogate (weight-dependent but seed-independent)
                mats = []
                A0 = build_A_normalized(c["W_orig"], c["scale_log1p"])
                L0 = build_L_from_A(A0)
                mats.append(("orig", -1, L0))

                for surr_id, W_s in enumerate(c["W_surr"]):
                    As = build_A_normalized(W_s, c["scale_log1p"])
                    Ls = build_L_from_A(As)
                    mats.append(("surr", surr_id, Ls))

                # Per base seed → per subject hashed effective seed
                for base_seed in base_seeds:
                    seed_eff = sha256_seed_int32(base_seed, subj)
                    rng = np.random.default_rng(seed_eff)
                    xi_full = rng.normal(0.0, args.eps_dis, size=len(c["full_nodes_sorted"]))
                    xi_map = {n: float(xi_full[i]) for i, n in enumerate(c["full_nodes_sorted"])}
                    xi_sub = np.array([xi_map[n] for n in c["nodes"]], dtype=float)

                    for kind, surr_id, L in mats:
                        H = build_H_from_L_xi(L, xi_sub, args.J)
                        R0max = simulate_R0max_kappa_batch_RK2(
                            H=H,
                            sources_idx=c["sources_idx"],
                            targets_idx=c["targets_idx"],
                            kappa_grid=c["kappa_probes"],
                            gamma_T=args.gamma_T,
                            gamma_L=args.gamma_L,
                            dt=args.dt,
                            T_end=args.T_end,
                            burn_in=args.burn_in,
                            eps_stab=args.eps_stab,
                        )
                        proxy = proxy_from_R0max(R0max, c["kappa_probes"])

                        row = {
                            "subject": subj,
                            "port": port,
                            "B_freeze_pct": B,
                            "base_seed": int(base_seed),
                            "seed_eff": int(seed_eff),
                            "trial_kind": kind,
                            "shuffle_id": int(surr_id),
                            "n_nodes": int(len(c["nodes"])),
                            "n_edges": int(np.count_nonzero(c["W_orig"]) // 2),
                            "kappa_grid": "T3_7pt" if port == "T3" else "T2_5pt",
                            "R_edge": proxy["R_edge"],
                            "R_peak": proxy["R_peak"],
                            "kappa_star": proxy["kappa_star"],
                            "G_proxy": proxy["G_proxy"],
                        }
                        # Store curve values (sparse, but useful for debugging)
                        for kk, val in zip(c["kappa_probes"], R0max):
                            row[f"R0max_k{kk:g}"] = float(val)
                        all_rows.append(row)

    df = pd.DataFrame(all_rows)
    df.to_csv(outdir / "results" / "all_trials_hashed12.csv", index=False)

    # Per-seed summary: (subject,port,B,seed_eff): orig G, median(surr G), delta
    rows = []
    for (subj, port, B, base_seed, seed_eff), g in df.groupby(["subject","port","B_freeze_pct","base_seed","seed_eff"]):
        g_orig = g[g["trial_kind"]=="orig"]["G_proxy"].values
        g_surr = g[g["trial_kind"]=="surr"]["G_proxy"].values
        if len(g_orig)!=1 or len(g_surr)==0:
            continue
        G_orig = float(g_orig[0])
        med_surr = float(np.median(g_surr))
        delta = float(med_surr - G_orig)  # positive means surrogates better
        rows.append({
            "subject": subj, "port": port, "B_freeze_pct": B,
            "base_seed": int(base_seed), "seed_eff": int(seed_eff),
            "G_orig": G_orig, "G_surr_median": med_surr, "delta_surr_minus_orig": delta,
            "kappa_star_orig": float(g[g["trial_kind"]=="orig"]["kappa_star"].values[0]),
        })
    df_seed = pd.DataFrame(rows)
    df_seed.to_csv(outdir / "results" / "per_seed_summary.csv", index=False)

    # Per-subject CI labels (hashed-only): hierarchical bootstrap over seeds/shuffles
    label_rows = []
    for (subj, port, B), g in df.groupby(["subject","port","B_freeze_pct"]):
        # build deltas_by_seed: each seed has vector of (G_surr - G_orig) across shuffles
        deltas_by_seed = {}
        for seed_eff, gg in g.groupby("seed_eff"):
            G_orig = float(gg[gg["trial_kind"]=="orig"]["G_proxy"].values[0])
            g_s = gg[gg["trial_kind"]=="surr"]["G_proxy"].values.astype(float)
            deltas_by_seed[int(seed_eff)] = (g_s - G_orig)
        ci = hierarchical_bootstrap_delta(deltas_by_seed, n_boot=3000, seed=stable_seed("boot", subj, port, str(B)))
        lab = label_from_ci(ci[0], ci[1])
        # point estimate: median over seeds of (median_surr - orig)
        seed_level = []
        for seed_eff, arr in deltas_by_seed.items():
            seed_level.append(float(np.median(arr)))
        point = float(np.median(seed_level)) if seed_level else float("nan")
        label_rows.append({
            "subject": subj, "port": port, "B_freeze_pct": B,
            "delta_point_median": point,
            "delta_CI_low": ci[0],
            "delta_CI_high": ci[1],
            "label": lab,
        })
    df_lab = pd.DataFrame(label_rows)
    df_lab.to_csv(outdir / "results" / "per_subject_labels.csv", index=False)

    # Group summary across subjects (simple): distribution of G_orig and delta
    group_rows = []
    for (port, B), g in df_seed.groupby(["port","B_freeze_pct"]):
        group_rows.append({
            "port": port,
            "B_freeze_pct": B,
            "G_orig_median": float(np.median(g["G_orig"])),
            "G_orig_mean": float(np.mean(g["G_orig"])),
            "delta_median": float(np.median(g["delta_surr_minus_orig"])),
            "delta_mean": float(np.mean(g["delta_surr_minus_orig"])),
            "frac_positive_Gorig": float(np.mean(g["G_orig"]>0)),
        })
    df_group = pd.DataFrame(group_rows)
    df_group.to_csv(outdir / "results" / "group_summary.csv", index=False)

    # Quick plots
    for port in ports:
        dff = df_seed[df_seed["port"]==port]
        plt.figure(figsize=(10,5))
        dff.boxplot(column="delta_surr_minus_orig", by="B_freeze_pct")
        plt.title(f"{port}: delta = median(G_surr) - G_orig (hashed-only, 12 seeds)")
        plt.suptitle("")
        plt.xlabel("B freeze (%)")
        plt.ylabel("delta")
        plt.tight_layout()
        plt.savefig(outdir / "plots" / f"{port}_delta_boxplot.png", dpi=160)
        plt.close()

    # Summary markdown
    summary = []
    summary.append("# Step 11 — hashed-only (12 seeds) results\n")
    summary.append(f"- Subjects: {subjects}\n- Ports: {ports}\n- B: {B_list}\n- Seeds (hashed-only): {base_seeds[0]}..{base_seeds[-1]} (12)\n- n_shuffles: {args.n_shuffles}\n")
    summary.append("## Key point\nThis dataset is **hashed-only** (subject-specific disorder). Group N_eff is defined by independent subject×seed draws.\n")
    summary.append("## Files\n- results/all_trials_hashed12.csv\n- results/per_seed_summary.csv\n- results/per_subject_labels.csv\n- results/group_summary.csv\n")
    (outdir / "SUMMARY.md").write_text("\n".join(summary), encoding="utf-8")

    print("DONE. Wrote:", outdir)


if __name__ == "__main__":
    main()
