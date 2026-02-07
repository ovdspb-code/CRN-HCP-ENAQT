#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step15: CTRW on surrogates (remaining 3 subjects) on Canonical-v2 pipeline.

Goal:
  Show that classical Pauli master equation / CTRW on the same subgraphs and surrogate controls
  does NOT produce an ENAQT inverted-U (G_ENAQT == 0), even when weights are shuffled under
  FreezeTop5% + FreezeBetweenness(B%) constraints.

Subjects: 304727, 314225, 519647
Ports:    T2 (BG->Thal/Pall), T3 (Thal->Motor)
B list:   {0, 10}
K_extra:  15
Surrogates per case: 3

CTRW model:
  dp/dt = Q p - (gamma_T*I_T + gamma_L*I_notT) p + kappa*(mean(p)-p)
  dP_T/dt = gamma_T * sum_{i in targets} p_i
  dP_D/dt = gamma_L * sum_{i not in targets} p_i

Rates:
  W_rate = J * A_norm
  where A_norm = log1p(W_raw)/scale_log1p(full_graph)

Integration:
  RK2 midpoint, dt=0.05, T_end=10.0, burnin=0.5
Metric:
  R0(t) = P_T(t) / (P_D(t) + eps_stab)
  R0_max(kappa) = max_{t>=burnin} R0(t)
  G_ENAQT = max_k R0_max(k) - max(R0_max(k_min), R0_max(k_max))
"""
from __future__ import annotations

import os, math, json
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

WEIGHT_KEY = "number_of_fibers"

def canonical_node_sort(nodes):
    return sorted(nodes, key=lambda x: int(x))

def build_full_weight_matrix(G: nx.Graph, weight_key: str = WEIGHT_KEY):
    nodes = canonical_node_sort(G.nodes())
    idx = {n:i for i,n in enumerate(nodes)}
    N=len(nodes)
    W=np.zeros((N,N), dtype=float)
    for u,v,dat in G.edges(data=True):
        w=float(dat.get(weight_key,0.0))
        if w<=0:
            continue
        i,j=idx[u], idx[v]
        W[i,j]=W[j,i]=w
    return W, nodes, idx

def get_port_nodes(battery: dict, task_key: str, subj_id: str):
    task = battery["tasks"][task_key]["per_subject"][str(subj_id)]
    sources = [d["id"] for d in task["sources"]]
    targets = [d["id"] for d in task["targets"]]
    return sources, targets

def select_subgraph_nodes(G: nx.Graph, sources, targets, K_extra=15, weight_key: str = WEIGHT_KEY):
    core = list(dict.fromkeys([*sources, *targets]))  # preserve order
    core_set=set(core)
    scores={}
    for c in core:
        for nbr in G.neighbors(c):
            if nbr in core_set:
                continue
            w=float(G.edges[c,nbr].get(weight_key,0.0))
            scores[nbr]=scores.get(nbr,0.0)+w
    ranked=sorted(scores.items(), key=lambda kv: (-kv[1], int(kv[0])))
    extra=[nid for nid,_ in ranked[:K_extra]]
    nodes = core + extra
    return nodes

def build_sub_weight_matrix(G: nx.Graph, nodes, weight_key: str = WEIGHT_KEY):
    idx={n:i for i,n in enumerate(nodes)}
    N=len(nodes)
    W=np.zeros((N,N), dtype=float)
    for u,v,dat in G.edges(data=True):
        if u in idx and v in idx:
            w=float(dat.get(weight_key,0.0))
            i,j=idx[u], idx[v]
            W[i,j]=W[j,i]=w
    return W, idx

def weight_freeze_set(W: np.ndarray, frac: float=0.05):
    edges=[(i,j,W[i,j]) for i in range(W.shape[0]) for j in range(i+1,W.shape[1]) if W[i,j]>0]
    if not edges:
        return set()
    edges_sorted=sorted(edges, key=lambda x: -x[2])
    k=max(1, int(math.ceil(frac*len(edges_sorted))))
    return set((i,j) for i,j,_ in edges_sorted[:k])

def betweenness_freeze_set(nodes, W: np.ndarray, percent: float=10.0):
    if percent<=0:
        return set()
    SG=nx.Graph()
    SG.add_nodes_from(nodes)
    for i,u in enumerate(nodes):
        for j in range(i+1,len(nodes)):
            if W[i,j]>0:
                SG.add_edge(u, nodes[j])
    if SG.number_of_edges()==0:
        return set()
    bc=nx.edge_betweenness_centrality(SG, normalized=True, weight=None)
    edges_sorted=sorted(bc.items(), key=lambda kv: (-kv[1], int(min(kv[0])), int(max(kv[0]))))
    m=len(edges_sorted)
    k=max(1, int(math.ceil((percent/100.0)*m)))
    frozen=set()
    node_to_i={n:i for i,n in enumerate(nodes)}
    for (u,v),_ in edges_sorted[:k]:
        i=node_to_i[u]; j=node_to_i[v]
        frozen.add((i,j) if i<j else (j,i))
    return frozen

def symmetric_ipf_scale(U: np.ndarray, r: np.ndarray, max_iter=5000, tol=1e-8):
    N=U.shape[0]
    U=np.maximum(U,0.0)
    np.fill_diagonal(U,0.0)
    x=np.ones(N, dtype=float)
    x[r<=0]=0.0
    for it in range(max_iter):
        Ux=U.dot(x)
        s=x*Ux
        mask=r>0
        if not mask.any():
            break
        ratio=np.ones_like(r)
        ok=mask & (s>0)
        ratio[ok]=r[ok]/s[ok]
        x_new=x*np.sqrt(ratio)
        if it%10==0:
            Ux_new=U.dot(x_new)
            s_new=x_new*Ux_new
            err=np.max(np.abs(s_new[mask]-r[mask])/(np.maximum(r[mask],1e-12)))
            if err<tol:
                x=x_new
                break
        x=x_new
    return (x[:,None]*U)*x[None,:]

def strength_preserving_shuffle(W_orig: np.ndarray, frozen_set: set[tuple[int,int]], rng: np.random.Generator):
    N=W_orig.shape[0]
    unfrozen=[(i,j) for i in range(N) for j in range(i+1,N) if W_orig[i,j]>0 and (i,j) not in frozen_set]
    weights=[W_orig[i,j] for (i,j) in unfrozen]
    perm=weights.copy()
    rng.shuffle(perm)
    Wtmp=W_orig.copy()
    for w,(i,j) in zip(perm,unfrozen):
        Wtmp[i,j]=Wtmp[j,i]=w
    frozen_mat=np.zeros_like(W_orig)
    for (i,j) in frozen_set:
        frozen_mat[i,j]=frozen_mat[j,i]=W_orig[i,j]  # keep original weights for frozen edges
    s_orig=W_orig.sum(axis=1)
    r=s_orig - frozen_mat.sum(axis=1)
    r[r<0]=0.0
    U=np.maximum(Wtmp - frozen_mat, 0.0)
    np.fill_diagonal(U,0.0)
    U_scaled=symmetric_ipf_scale(U, r)
    return frozen_mat + U_scaled

def ctrw_curve_R0max(A_norm: np.ndarray, sources_idx, targets_idx, kappa_list,
                     dt=0.05, T_end=10.0, burnin=0.5, gamma_T=1.0, gamma_L=5e-4, eps_stab=0.001, J=0.1):
    n=A_norm.shape[0]
    W_rate = J*A_norm.copy()
    np.fill_diagonal(W_rate, 0.0)
    out_rate=W_rate.sum(axis=0)
    Q=W_rate.copy()
    for j in range(n):
        Q[j,j]=-out_rate[j]
    targets_mask=np.zeros(n, dtype=float); targets_mask[targets_idx]=1.0
    leak_mask=1.0-targets_mask
    sink_vec=gamma_T*targets_mask
    leak_vec=gamma_L*leak_mask
    p0=np.zeros(n, dtype=float)
    p0[sources_idx]=1.0/len(sources_idx)
    steps=int(round(T_end/dt))
    burn_steps=int(round(burnin/dt))
    R0max=[]
    for kappa in kappa_list:
        p=p0.copy(); PT=0.0; PD=0.0; r0_max=0.0
        leak_bool=leak_mask.astype(bool)
        for step in range(steps):
            def deriv(p, PT, PD):
                dp=Q.dot(p) - (sink_vec+leak_vec)*p
                if kappa>0:
                    m=p.sum()/n
                    dp += kappa*(m-p)
                dPT=gamma_T*p[targets_idx].sum()
                dPD=gamma_L*p[leak_bool].sum()
                return dp,dPT,dPD
            dp1,dpt1,dpd1=deriv(p,PT,PD)
            p_mid=p+0.5*dt*dp1
            PT_mid=PT+0.5*dt*dpt1
            PD_mid=PD+0.5*dt*dpd1
            dp2,dpt2,dpd2=deriv(p_mid,PT_mid,PD_mid)
            p=p+dt*dp2
            PT=PT+dt*dpt2
            PD=PD+dt*dpd2
            p=np.maximum(p,0.0)
            if step+1>=burn_steps:
                r0=PT/(PD+eps_stab)
                if r0>r0_max:
                    r0_max=r0
        R0max.append(r0_max)
    return np.array(R0max)

def compute_G_ENAQT(curve: np.ndarray):
    m=float(curve.max())
    edge=float(max(curve[0], curve[-1]))
    return m-edge, int(np.argmax(curve))

def main():
    subjects=["304727","314225","519647"]
    ports={"T2":"T2_BG_to_Thal_Pall_LH","T3":"T3_Thal_to_Motor_LH"}
    B_list=[0,10]
    K_extra=15
    n_surr=3
    kappa_list=np.logspace(-3,1,21)
    dt=0.05; T_end=10.0; burnin=0.5
    gamma_T=1.0; gamma_L=5e-4; eps_stab=0.001; J=0.1

    with open("human_task_battery_v1.json","r") as f:
        battery=json.load(f)

    out_rows=[]
    curves=[]
    os.makedirs("plots", exist_ok=True)

    for subj in subjects:
        G=nx.read_graphml(f"{subj}_repeated10_scale125.graphml")
        Wfull,_,_=build_full_weight_matrix(G)
        scale_log1p=float(np.log1p(Wfull.max())) if Wfull.max()>0 else 1.0

        for pshort, task_key in ports.items():
            sources, targets=get_port_nodes(battery, task_key, subj)
            nodes_sub=select_subgraph_nodes(G, sources, targets, K_extra=K_extra)
            Wraw, idx=build_sub_weight_matrix(G, nodes_sub)
            A_orig=np.log1p(Wraw)/scale_log1p
            s_idx=[idx[s] for s in sources]
            t_idx=[idx[t] for t in targets]

            for B in B_list:
                frozen=weight_freeze_set(Wraw, 0.05) | betweenness_freeze_set(nodes_sub, Wraw, float(B))
                # orig
                curve=ctrw_curve_R0max(A_orig, s_idx, t_idx, kappa_list, dt,T_end,burnin,gamma_T,gamma_L,eps_stab,J)
                Gval,argm=compute_G_ENAQT(curve)
                out_rows.append(dict(subject=subj,port=pshort,B=B,variant="orig",kappa_star=float(kappa_list[argm]),
                                     R0_star=float(curve.max()),R0_kmin=float(curve[0]),R0_kmax=float(curve[-1]),G_ENAQT=float(Gval)))
                for k,r0 in zip(kappa_list, curve):
                    curves.append(dict(subject=subj,port=pshort,B=B,variant="orig",kappa=float(k),R0max=float(r0)))

                # surrogates
                base_rng=np.random.default_rng((int(subj)*100 + (0 if pshort=="T2" else 1)*10 + B) % (2**32-1))
                for si in range(n_surr):
                    rng=np.random.default_rng(base_rng.integers(0,2**32-1))
                    Wsh=strength_preserving_shuffle(Wraw, frozen, rng)
                    A_sh=np.log1p(Wsh)/scale_log1p
                    curve_sh=ctrw_curve_R0max(A_sh, s_idx, t_idx, kappa_list, dt,T_end,burnin,gamma_T,gamma_L,eps_stab,J)
                    Gs,argms=compute_G_ENAQT(curve_sh)
                    out_rows.append(dict(subject=subj,port=pshort,B=B,variant=f"surr{si}",kappa_star=float(kappa_list[argms]),
                                         R0_star=float(curve_sh.max()),R0_kmin=float(curve_sh[0]),R0_kmax=float(curve_sh[-1]),G_ENAQT=float(Gs)))
                    for k,r0 in zip(kappa_list, curve_sh):
                        curves.append(dict(subject=subj,port=pshort,B=B,variant=f"surr{si}",kappa=float(k),R0max=float(r0)))

                # plot
                df=pd.DataFrame([c for c in curves if c["subject"]==subj and c["port"]==pshort and c["B"]==B])
                plt.figure()
                for var,g in df.groupby("variant"):
                    g=g.sort_values("kappa")
                    plt.plot(g["kappa"],g["R0max"],marker='o',linewidth=1,markersize=3,label=var)
                plt.xscale("log")
                plt.xlabel("κ (CTRW mixing-rate)")
                plt.ylabel("R0_max")
                plt.title(f"CTRW R0max vs κ | subj {subj} | {pshort} | B={B} | K=15")
                plt.legend(fontsize=7)
                plt.tight_layout()
                plt.savefig(f"plots/CTRW_{subj}_{pshort}_B{B}.png",dpi=200)
                plt.close()

    pd.DataFrame(out_rows).to_csv("ctrw_summary.csv",index=False)
    pd.DataFrame(curves).to_csv("ctrw_curves.csv",index=False)
    print("Done. All G_ENAQT are expected to be 0 (no interior peak).")

if __name__=="__main__":
    main()
