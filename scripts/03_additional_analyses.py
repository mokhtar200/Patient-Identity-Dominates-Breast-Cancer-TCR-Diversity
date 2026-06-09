#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Additional Analyses  —  TCRb Repertoire Pipeline
=================================================
All analyses save tables + figures into:
    D:\\PRJNA301507\\analysis_june_2026\\additional_analyses
"""

import os, sys, re, warnings, itertools, json
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import Counter

warnings.filterwarnings("ignore")
RNG = np.random.default_rng(20240517)
np.random.seed(20240517)

DATA_DIR = r"D:\PRJNA301507\data"
OUT = r"D:\PRJNA301507\analysis_june_2026\additional_analyses"
META_FILE = os.path.join(DATA_DIR, "sample_metadata.csv")
SUBSETS = ["Treg", "CD4_Memory", "Other"]
SCOL = {"Treg": "#E64B35", "CD4_Memory": "#4DBBD5", "Other": "#00A087"}
COMPARTMENTS = ["Blood", "Tumor", "Lymph node"]
CCOL = {"Blood": "#E64B35", "Tumor": "#4DBBD5", "Lymph node": "#00A087"}
AA20 = set("ACDEFGHIKLMNPQRSTVWY")

os.makedirs(OUT, exist_ok=True)
os.makedirs(os.path.join(OUT, "figures"), exist_ok=True)
os.makedirs(os.path.join(OUT, "tables"), exist_ok=True)

def log(msg): print(f"[ADD] {msg}", flush=True)
def FIG(name): return os.path.join(OUT, "figures", name)
def TBL(name): return os.path.join(OUT, "tables", name)

# ---------------------------------------------------------------------------
# 1. DATA LOADING
# ---------------------------------------------------------------------------
COLUMN_ALIASES = {
    "cdr3aa": ["cdr3aa","aaseqcdr3","amino_acid","cdr3_amino_acid","cdr3","aaseq",
               "cdr3_aa","junction_aa","cdr3.aa","aaseq"],
    "v":      ["v","vgene","bestvgene","v_gene","vgenename","bestvhit","v_call","vgene_name","vregion"],
    "j":      ["j","jgene","bestjgene","j_gene","jgenename","bestjhit","j_call","jgene_name","jregion"],
    "count":  ["count","clonecount","templates","reads","productive_templates",
               "duplicate_count","clones","cloneCount","frequencycount","readcount"],
    "freq":   ["freq","frequency","productive_frequency","clonefraction","clonefreq","clone_fraction"],
}

def _match_col(columns, aliases):
    low = {c.lower().strip(): c for c in columns}
    for a in aliases:
        if a in low: return low[a]
    for a in aliases:
        for lc, orig in low.items():
            if a in lc: return orig
    return None

def read_clonotype_table(path):
    sep = "\t" if path.lower().endswith((".tsv", ".txt")) else None
    try:
        df = pd.read_csv(path, sep=sep, engine="python")
    except Exception:
        return None
    if df.shape[1] == 1:
        try: df = pd.read_csv(path, sep="\t", engine="python")
        except: pass
    cols = list(df.columns)
    chain_col = _match_col(cols, ["chain"])
    if chain_col is not None:
        df = df[df[chain_col].astype(str).str.upper().str.contains("TRB")]
    c_cdr3 = _match_col(cols, COLUMN_ALIASES["cdr3aa"])
    c_v = _match_col(cols, COLUMN_ALIASES["v"])
    c_j = _match_col(cols, COLUMN_ALIASES["j"])
    c_cnt = _match_col(cols, COLUMN_ALIASES["count"])
    if any(x is None for x in (c_cdr3, c_v, c_j, c_cnt)):
        return None
    out = pd.DataFrame({
        "cdr3aa": df[c_cdr3].astype(str),
        "v": df[c_v].astype(str),
        "j": df[c_j].astype(str),
        "count": pd.to_numeric(df[c_cnt], errors="coerce")
    })
    return out

def standardise(df):
    if df is None or len(df) == 0: return None
    df = df.copy()
    df["v"] = df["v"].str.replace(r"\*.*$", "", regex=True)
    df["j"] = df["j"].str.replace(r"\*.*$", "", regex=True)
    df = df[df["cdr3aa"].notna() & (df["cdr3aa"] != "")]
    df = df[~df["cdr3aa"].str.contains(r"[*._]", regex=True, na=True)]
    df = df[df["cdr3aa"].apply(lambda s: len(s)>=5 and len(s)<=25 and set(s) <= AA20)]
    df = df[df["v"].str.startswith("TRBV") & df["j"].str.startswith("TRBJ")]
    df = df[pd.to_numeric(df["count"], errors="coerce").fillna(0) > 0]
    if len(df) == 0: return None
    df = df.groupby(["cdr3aa","v","j"], as_index=False)["count"].sum()
    df["freq"] = df["count"] / df["count"].sum()
    return df

log("Loading data...")
meta = pd.read_csv(META_FILE)
meta = meta[meta["subset"].isin(SUBSETS)].reset_index(drop=True)
samples = {}
for _, row in meta.iterrows():
    fp = os.path.join(DATA_DIR, str(row["file_name"]))
    if not os.path.exists(fp): fp = str(row["file_name"])
    s = standardise(read_clonotype_table(fp))
    if s is not None and len(s) >= 100 and s["count"].sum() >= 1000:
        samples[row["sample_id"]] = s
meta_qc = meta[meta["sample_id"].isin(samples)].reset_index(drop=True)
log(f"Loaded {len(samples)} samples, {meta_qc['patient_id'].nunique()} patients")

def pool_subset(sub):
    ids = meta_qc[meta_qc["subset"] == sub]["sample_id"].tolist()
    parts = [samples[i] for i in ids if i in samples]
    if not parts: return None
    bb = pd.concat(parts, ignore_index=True)
    bb = bb.groupby(["cdr3aa","v","j"], as_index=False)["count"].sum()
    bb["freq"] = bb["count"] / bb["count"].sum()
    return bb

subset_pools = {}
for s in SUBSETS:
    p = pool_subset(s)
    if p is not None: subset_pools[s] = p
log("Data ready.")

# =========================================================================
# A1: V/J Gene Usage
# =========================================================================
log("--- A1: V/J Gene Usage ---")
vg = {}
jg = {}
for s in SUBSETS:
    vg[s] = subset_pools[s].groupby("v")["count"].sum()
    vg[s] = vg[s] / vg[s].sum()
    jg[s] = subset_pools[s].groupby("j")["count"].sum()
    jg[s] = jg[s] / jg[s].sum()
all_v = sorted(set().union(*[set(v.index) for v in vg.values()]))
all_j = sorted(set().union(*[set(j.index) for j in jg.values()]))
v_df = pd.DataFrame({s: [vg[s].get(v, 0) for v in all_v] for s in SUBSETS}, index=all_v)
j_df = pd.DataFrame({s: [jg[s].get(j, 0) for j in all_j] for s in SUBSETS}, index=all_j)
v_df.to_csv(TBL("A1_V_gene_usage.csv")); j_df.to_csv(TBL("A1_J_gene_usage.csv"))

# Plot V genes
fig, ax = plt.subplots(figsize=(14, 5))
x = np.arange(len(all_v)); w = 0.25
for i, s in enumerate(SUBSETS):
    ax.bar(x + i*w, [vg[s].get(v,0)*100 for v in all_v], w, label=s, color=SCOL[s])
ax.set_xticks(x + w); ax.set_xticklabels(all_v, rotation=90, fontsize=7)
ax.set_ylabel("Frequency (%)"); ax.set_title("TRBV Gene Usage by Subset")
ax.legend(fontsize=8); plt.tight_layout(); fig.savefig(FIG("A1_V_gene_usage.png"), dpi=200)
plt.close(fig)

fig, ax = plt.subplots(figsize=(10, 4))
x = np.arange(len(all_j)); w = 0.25
for i, s in enumerate(SUBSETS):
    ax.bar(x + i*w, [jg[s].get(j,0)*100 for j in all_j], w, label=s, color=SCOL[s])
ax.set_xticks(x + w); ax.set_xticklabels(all_j, rotation=90, fontsize=8)
ax.set_ylabel("Frequency (%)"); ax.set_title("TRBJ Gene Usage by Subset")
ax.legend(fontsize=8); plt.tight_layout(); fig.savefig(FIG("A1_J_gene_usage.png"), dpi=200)
plt.close(fig)

# Differential V (chi-square)
v_rows = []
for v in all_v:
    tbl = np.array([[vg[s].get(v,0)*100, 100-vg[s].get(v,0)*100] for s in SUBSETS])
    try:
        _, p = stats.chi2_contingency(tbl)[:2]
        v_rows.append(dict(V=v, p=round(p, 5), **{s: round(vg[s].get(v,0)*100, 2) for s in SUBSETS}))
    except: pass
vdiff = pd.DataFrame(v_rows).sort_values("p")
vdiff.to_csv(TBL("A1_V_gene_differential.csv"), index=False)
log(f"  V genes: {len(all_v)}, most diff: {vdiff.iloc[0]['V'] if len(vdiff) else 'none'} (p={vdiff.iloc[0]['p']:.5f})")

# =========================================================================
# A2: CDR3 Sequence Properties
# =========================================================================
log("--- A2: CDR3 Properties ---")
def cdr3_props(df):
    lens = df["cdr3aa"].str.len()
    seqs = "".join(df["cdr3aa"].tolist())
    aa_cnt = Counter(seqs)
    total = sum(aa_cnt.values())
    return {"mean_len": lens.mean(), "median_len": lens.median(), "std_len": lens.std(),
            "min_len": lens.min(), "max_len": lens.max(),
            **{aa: aa_cnt.get(aa,0)/total for aa in sorted(AA20)}}

prop_rows = [{**cdr3_props(subset_pools[s]), "subset": s} for s in SUBSETS]
cdr3_props_df = pd.DataFrame(prop_rows).set_index("subset")
cdr3_props_df.to_csv(TBL("A2_CDR3_properties.csv"))
log(f"  Mean lengths: {cdr3_props_df['mean_len'].to_dict()}")

# Length distribution
fig, ax = plt.subplots(figsize=(8, 5))
for s in SUBSETS:
    lens = subset_pools[s]["cdr3aa"].str.len()
    cnt = lens.value_counts().sort_index()
    ax.plot(cnt.index, cnt.values/cnt.sum()*100, marker="o", label=s, color=SCOL[s])
ax.set_xlabel("CDR3 length (aa)"); ax.set_ylabel("Frequency (%)")
ax.set_title("CDR3 Length Distribution by Subset"); ax.legend()
plt.tight_layout(); fig.savefig(FIG("A2_CDR3_length_dist.png"), dpi=200)
plt.close(fig)

# AA composition
fig, ax = plt.subplots(figsize=(10, 4))
aa_chart = cdr3_props_df[[aa for aa in sorted(AA20)]].T
im = ax.imshow(aa_chart.values, cmap="Reds", aspect="auto")
ax.set_xticks(range(len(SUBSETS))); ax.set_xticklabels(SUBSETS)
ax.set_yticks(range(len(aa_chart.index))); ax.set_yticklabels(aa_chart.index)
ax.set_title("Amino Acid Composition by Subset")
plt.colorbar(im, label="Frequency"); plt.tight_layout()
fig.savefig(FIG("A2_AA_composition.png"), dpi=200); plt.close(fig)

# =========================================================================
# A3: Differential Clonal Abundance (optimized - only well-represented clones)
# =========================================================================
log("--- A3: Differential Clonal Abundance ---")

# Build clone x subset matrix efficiently
all_clones = []
for s in SUBSETS:
    d = subset_pools[s][["cdr3aa","v","j","count"]].copy()
    d["subset"] = s
    all_clones.append(d)
clone_master = pd.concat(all_clones, ignore_index=True)
total_by_subset = clone_master.groupby("subset")["count"].sum()

# Pivot: clone -> count per subset
clone_pivot = clone_master.pivot_table(index=["cdr3aa","v","j"], columns="subset",
                                        values="count", aggfunc="sum", fill_value=0)
for s in SUBSETS:
    if s not in clone_pivot.columns: clone_pivot[s] = 0.0
clone_pivot = clone_pivot[SUBSETS]

# Filter to clones with minimum total count (speed + relevance)
min_count = max(10, int(total_by_subset.min() * 0.0001))
mask = (clone_pivot["Treg"] + clone_pivot["CD4_Memory"]) >= min_count
clone_pivot_filt = clone_pivot[mask]
log(f"  Testing {len(clone_pivot_filt)}/{len(clone_pivot)} clones (min total count={min_count})")

# Vectorized approximate test: log2 FC + z-test
treg_total = total_by_subset["Treg"]
mem_total = total_by_subset["CD4_Memory"]

treg_p = clone_pivot_filt["Treg"].values / treg_total
mem_p = clone_pivot_filt["CD4_Memory"].values / mem_total
# Add pseudocount
treg_p = (clone_pivot_filt["Treg"].values + 1) / (treg_total + 1)
mem_p = (clone_pivot_filt["CD4_Memory"].values + 1) / (mem_total + 1)
log2fc = np.log2(treg_p / mem_p)

# z-test for two proportions
p_pool = (clone_pivot_filt["Treg"].values + clone_pivot_filt["CD4_Memory"].values) / (treg_total + mem_total)
se = np.sqrt(p_pool * (1-p_pool) * (1/treg_total + 1/mem_total))
z = (treg_p - mem_p) / (se + 1e-15)
pvals = 2 * (1 - stats.norm.cdf(np.abs(z)))
p_fdr = stats.false_discovery_control(pvals, method="bh")

fd = pd.DataFrame({
    "cdr3aa": clone_pivot_filt.index.get_level_values(0),
    "v": clone_pivot_filt.index.get_level_values(1),
    "j": clone_pivot_filt.index.get_level_values(2),
    "count_Treg": clone_pivot_filt["Treg"].values,
    "count_Memory": clone_pivot_filt["CD4_Memory"].values,
    "log2FC": np.round(log2fc, 3),
    "p": np.round(pvals, 6),
    "p_fdr": np.round(p_fdr, 6),
})
fd = fd.sort_values("p_fdr")
fd.to_csv(TBL("A3_differential_clones.csv"), index=False)
sig = fd[fd["p_fdr"] < 0.05]
log(f"  Clones tested: {len(fd)}, FDR<0.05: {len(sig)} "
    f"(Treg-enriched: {(sig.log2FC>0).sum()}, Mem-enriched: {(sig.log2FC<0).sum()})")

# Volcano plot
if fd["log2FC"].nunique() > 1:
    fig, ax = plt.subplots(figsize=(7, 6))
    fd["neg_log10_p"] = -np.log10(fd["p_fdr"].clip(1e-15))
    colors = ["#E64B35" if r.log2FC > 1 and r.p_fdr < 0.05 else
              "#4DBBD5" if r.log2FC < -1 and r.p_fdr < 0.05 else "grey"
              for _, r in fd.iterrows()]
    ax.scatter(fd["log2FC"], fd["neg_log10_p"], c=colors, s=3, alpha=0.6)
    ax.axhline(-np.log10(0.05), ls="--", color="grey", alpha=0.5)
    ax.axvline(1, ls="--", color="grey", alpha=0.3)
    ax.axvline(-1, ls="--", color="grey", alpha=0.3)
    ax.set_xlabel("log2 FC (Treg / CD4_Memory)"); ax.set_ylabel("-log10(FDR)")
    ax.set_title("Differential Clonal Abundance")
    plt.tight_layout(); fig.savefig(FIG("A3_volcano.png"), dpi=200); plt.close(fig)

# =========================================================================
# A4: Compartment-Specific Analyses
# =========================================================================
log("--- A4: Compartment-Specific ---")
comp_rows = []
for _, r in meta_qc.iterrows():
    s = samples.get(r["sample_id"])
    if s is None: continue
    p = s["freq"].values; H = -np.sum(p * np.log(p))
    comp_rows.append(dict(sample_id=r["sample_id"], patient_id=r["patient_id"],
                          subset=r["subset"], compartment=r["compartment"],
                          Shannon=H, Sobs=len(s)))
cd = pd.DataFrame(comp_rows)
cd.to_csv(TBL("A4_compartment_diversity.csv"), index=False)

# Boxplots
fig, axes = plt.subplots(1, 3, figsize=(14, 4))
for i, s in enumerate(SUBSETS):
    ax = axes[i]
    data = [cd[(cd.subset==s)&(cd.compartment==c)]["Shannon"].dropna() for c in COMPARTMENTS]
    bp = ax.boxplot(data, labels=COMPARTMENTS, patch_artist=True)
    for patch, c in zip(bp["boxes"], COMPARTMENTS): patch.set_facecolor(CCOL[c])
    ax.set_title(f"{s}"); ax.set_ylabel("Shannon H'" if i==0 else "")
plt.tight_layout(); fig.savefig(FIG("A4_compartment_shannon.png"), dpi=200)
plt.close(fig)

# Compartment overlap
overlap_rows = []
for pid in meta_qc["patient_id"].unique():
    pid_rows = meta_qc[meta_qc["patient_id"] == pid]
    for c1, c2 in itertools.combinations(COMPARTMENTS, 2):
        ids1 = pid_rows[pid_rows.compartment==c1]["sample_id"].tolist()
        ids2 = pid_rows[pid_rows.compartment==c2]["sample_id"].tolist()
        if not ids1 or not ids2: continue
        parts1 = [samples[i]["cdr3aa"] for i in ids1 if i in samples]
        parts2 = [samples[i]["cdr3aa"] for i in ids2 if i in samples]
        if not parts1 or not parts2: continue
        cc1 = set(pd.concat(parts1).unique())
        cc2 = set(pd.concat(parts2).unique())
        if len(cc1) and len(cc2):
            overlap_rows.append(dict(patient_id=pid, c1=c1, c2=c2,
                                     shared=len(cc1&cc2), jaccard=round(len(cc1&cc2)/len(cc1|cc2), 4)))
odf = pd.DataFrame(overlap_rows)
odf.to_csv(TBL("A4_compartment_overlap.csv"), index=False)
log(f"  Compartment pairs: {len(odf)}, mean Jaccard={odf['jaccard'].mean():.3f}" if len(odf) else "")

# Blood vs Tumor Treg
dom_rows = []
for pid in meta_qc["patient_id"].unique():
    for comp in ["Blood", "Tumor"]:
        ids = meta_qc[(meta_qc.patient_id==pid)&(meta_qc.subset=="Treg")&(meta_qc.compartment==comp)]["sample_id"].tolist()
        if not ids: continue
        parts = [samples[i] for i in ids if i in samples]
        if not parts: continue
        cc = pd.concat(parts, ignore_index=True)
        cc = cc.groupby(["cdr3aa","v","j"], as_index=False)["count"].sum()
        dom_rows.append(dict(patient_id=pid, compartment=comp,
                             Treg_top_freq=cc["count"].max()/cc["count"].sum(), n=len(cc)))
tdf = pd.DataFrame(dom_rows)
tdf.to_csv(TBL("A4_Treg_dominance_by_compartment.csv"), index=False)
if len(tdf):
    fig, ax = plt.subplots(figsize=(6, 5))
    data = [tdf[tdf.compartment==c]["Treg_top_freq"].dropna() for c in ["Blood", "Tumor"]]
    bp = ax.boxplot(data, labels=["Blood","Tumor"], patch_artist=True)
    for patch, c in zip(bp["boxes"], ["Blood","Tumor"]): patch.set_facecolor(CCOL[c])
    ax.set_ylabel("Top clone frequency (Treg)"); ax.set_title("Treg Dominance: Blood vs Tumor")
    plt.tight_layout(); fig.savefig(FIG("A4_Treg_dominance_blood_vs_tumor.png"), dpi=200); plt.close(fig)
    if len(data[0])>2 and len(data[1])>2:
        _, p = stats.mannwhitneyu(data[0], data[1])
        log(f"  Treg top freq: Blood={np.median(data[0]):.3f}, Tumor={np.median(data[1]):.3f}, p={p:.4f}")

# =========================================================================
# A5: Enhanced Beta Diversity
# =========================================================================
log("--- A5: Beta Diversity ---")
all_clones_beta = pd.concat([samples[sid] for sid in samples], ignore_index=True)
all_clones_beta["clone"] = all_clones_beta["cdr3aa"]+"|"+all_clones_beta["v"]+"|"+all_clones_beta["j"]
top_clones = all_clones_beta.groupby("clone")["count"].sum().sort_values(ascending=False).head(5000).index.tolist()
sid = list(samples.keys())
M = np.zeros((len(sid), len(top_clones)))
cidx = {c:i for i,c in enumerate(top_clones)}
for ri, sid_i in enumerate(sid):
    s = samples[sid_i].copy()
    s["clone"] = s["cdr3aa"]+"|"+s["v"]+"|"+s["j"]
    s = s[s["clone"].isin(cidx)]
    for _, row in s.iterrows(): M[ri, cidx[row["clone"]]] = row["count"]
kr = M.sum(1) > 0; M = M[kr]; sid = [s for s,k in zip(sid, kr) if k]
lab = meta_qc.set_index("sample_id").loc[sid, "subset"].values
D = squareform(pdist(M, metric="braycurtis"))

# PERMDISP
def permdisp(D, labels, n_perm=999):
    labels = np.asarray(labels); groups = np.unique(labels); n = D.shape[0]
    dists = np.zeros(n)
    for g in groups:
        idx = np.where(labels == g)[0]
        if len(idx) > 0:
            for i in idx: dists[i] = np.sqrt(D[i, idx].mean())
    obs_f, _ = stats.f_oneway(*[dists[labels==g] for g in groups])
    perm_f = np.empty(n_perm)
    for k in range(n_perm):
        lp = RNG.permutation(labels)
        pdists = np.zeros(n)
        for g in groups:
            idx = np.where(lp == g)[0]
            if len(idx) > 0:
                for i in range(n):
                    if lp[i] == g: pdists[i] = np.sqrt(D[i, idx].mean())
        perm_f[k], _ = stats.f_oneway(*[pdists[lp==g] for g in groups])
    p = (np.sum(perm_f >= obs_f) + 1) / (n_perm + 1)
    return obs_f, p

f_disp, p_disp = permdisp(D, lab, 999)
log(f"  PERMDISP: F={f_disp:.3f}, p={p_disp:.3f}")
pd.DataFrame([dict(test="PERMDISP", F=f_disp, p=p_disp)]).to_csv(TBL("A5_PERMDISP.csv"), index=False)

# ANOSIM
def anosim(D, labels, n_perm=999):
    labels = np.asarray(labels); n = D.shape[0]
    r_between, r_within = [], []
    for i in range(n):
        for j in range(i+1, n):
            if labels[i] == labels[j]: r_within.append(D[i,j])
            else: r_between.append(D[i,j])
    R = (np.mean(r_between) - np.mean(r_within)) * 2 / n if n > 1 else 0
    perm_R = []
    for _ in range(n_perm):
        lp = RNG.permutation(labels)
        pw, pb = [], []
        for i in range(n):
            for j in range(i+1, n):
                if lp[i] == lp[j]: pw.append(D[i,j])
                else: pb.append(D[i,j])
        perm_R.append((np.mean(pb)-np.mean(pw))*2/n if n > 1 else 0)
    p = (np.sum(np.array(perm_R) >= R) + 1) / (n_perm + 1)
    return R, p

R_ano, p_ano = anosim(D, lab, 999)
log(f"  ANOSIM: R={R_ano:.3f}, p={p_ano:.3f}")
pd.DataFrame([dict(test="ANOSIM", R=R_ano, p=p_ano)]).to_csv(TBL("A5_ANOSIM.csv"), index=False)

# =========================================================================
# A6: Clonal Dominance
# =========================================================================
log("--- A6: Clonal Dominance ---")
dom_rows = []
for sid_i in samples:
    s = samples[sid_i]
    cnt = s["count"].values / s["count"].sum()
    dom_rows.append(dict(sample_id=sid_i, S=len(s), top1=cnt[0], top5=cnt[:5].sum(),
                         top10=cnt[:10].sum(), top20=cnt[:20].sum(),
                         clones_gt_1pc=int(np.sum(cnt > 0.01)),
                         clones_gt_5pc=int(np.sum(cnt > 0.05))))
ddf = pd.DataFrame(dom_rows)
ddf = ddf.merge(meta_qc[["sample_id","subset","patient_id","compartment"]], on="sample_id")
ddf.to_csv(TBL("A6_clonal_dominance.csv"), index=False)

fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for i, metric in enumerate(["top1", "top5", "top10"]):
    ax = axes[i]
    data = [ddf[ddf.subset==s][metric].dropna() for s in SUBSETS]
    bp = ax.boxplot(data, labels=SUBSETS, patch_artist=True)
    for patch, s in zip(bp["boxes"], SUBSETS): patch.set_facecolor(SCOL[s])
    ax.set_title(metric)
plt.tight_layout(); fig.savefig(FIG("A6_clonal_dominance.png"), dpi=200); plt.close(fig)

# =========================================================================
# A7: Public Clonotype Analysis
# =========================================================================
log("--- A7: Public Clonotypes ---")
patient_cdr3s = {}
for pid in meta_qc["patient_id"].unique():
    ids = meta_qc[meta_qc["patient_id"]==pid]["sample_id"].tolist()
    parts = [samples[i] for i in ids if i in samples]
    if not parts: continue
    cc = pd.concat(parts, ignore_index=True)
    patient_cdr3s[pid] = set(cc["cdr3aa"].unique())
sharing = pd.Series([c for s in patient_cdr3s.values() for c in s]).value_counts()
sharing.to_csv(TBL("A7_clone_sharing_freq.csv"))
log(f"  CDR3s: {len(sharing)} total, {(sharing>=2).sum()} public, {(sharing==1).sum()} private")

fig, ax = plt.subplots(figsize=(8, 5))
hist = sharing.value_counts().sort_index()
ax.bar(hist.index, hist.values, width=0.8, color="steelblue", edgecolor="k")
ax.set_xlabel("Patients sharing a CDR3"); ax.set_ylabel("Number of CDR3s")
ax.set_title("CDR3 Sharing Frequency"); ax.set_yscale("log")
plt.tight_layout(); fig.savefig(FIG("A7_sharing_distribution.png"), dpi=200); plt.close(fig)

# =========================================================================
# A8: ML Classification
# =========================================================================
log("--- A8: ML Classification ---")
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler

feat_rows = []
for _, r in meta_qc.iterrows():
    s = samples.get(r["sample_id"])
    if s is None: continue
    p = s["freq"].values; cnt = s["count"].values
    H = -np.sum(p * np.log(p + 1e-15))
    feat = {"Shannon": H, "invSimpson": 1/np.sum(p**2), "S": len(s),
            "top1": cnt[0]/cnt.sum(), "top5": cnt[:5].sum()/cnt.sum(),
            "mean_CDR3_len": s["cdr3aa"].str.len().mean()}
    v_f = s.groupby("v")["count"].sum(); v_f = v_f / v_f.sum()
    feat["V_entropy"] = -np.sum(v_f * np.log(v_f + 1e-15))
    feat["sample_id"] = r["sample_id"]; feat["subset"] = r["subset"]
    feat_rows.append(feat)
fdf = pd.DataFrame(feat_rows)

bin_df = fdf[fdf["subset"].isin(["Treg","CD4_Memory"])].copy()
if len(bin_df) > 10:
    X = StandardScaler().fit_transform(bin_df.drop(columns=["sample_id","subset"]).values)
    y = (bin_df["subset"] == "Treg").astype(int)
    rf = RandomForestClassifier(n_estimators=200, random_state=2024, class_weight="balanced")
    aucs = cross_val_score(rf, X, y, cv=StratifiedKFold(5, shuffle=True, random_state=2024), scoring="roc_auc")
    log(f"  RF Treg vs CD4_Memory: AUC={np.mean(aucs):.3f} +/- {np.std(aucs):.3f}")
    pd.DataFrame([dict(mean_AUC=np.mean(aucs), std_AUC=np.std(aucs))]).to_csv(TBL("A8_RF_classification.csv"), index=False)
    rf.fit(X, y)
    imp = pd.DataFrame({"feature": bin_df.drop(columns=["sample_id","subset"]).columns,
                        "importance": rf.feature_importances_}).sort_values("importance", ascending=False)
    imp.to_csv(TBL("A8_RF_feature_importance.csv"), index=False)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.barh(imp["feature"], imp["importance"], color="steelblue")
    ax.set_xlabel("Importance"); ax.set_title("RF Feature Importance")
    plt.tight_layout(); fig.savefig(FIG("A8_RF_importance.png"), dpi=200); plt.close(fig)

# =========================================================================
# A9: Mixed Model Extensions
# =========================================================================
log("--- A9: Mixed Models ---")
import statsmodels.formula.api as smf

cd_model = cd.dropna(subset=["Shannon"])
try:
    md = smf.mixedlm("Shannon ~ subset * compartment", cd_model, groups=cd_model["patient_id"])
    mf = md.fit(reml=True, method="lbfgs")
    with open(TBL("A9_mixed_model_interaction.txt"), "w") as f:
        f.write(str(mf.summary()))
    vp = float(mf.cov_re.iloc[0,0]); vr = float(mf.scale)
    log(f"  Subset*Compartment: ICC={vp/(vp+vr):.3f}")
    fe = pd.DataFrame({"coef": mf.fe_params, "p": mf.pvalues})
    fe.to_csv(TBL("A9_mixed_model_fixed_effects.csv"))
except Exception as e:
    log(f"  Mixed model error: {e}")

# Variance partitioning
vp_rows = []
for formula in ["Shannon ~ 1", "Shannon ~ subset", "Shannon ~ compartment", "Shannon ~ subset + compartment"]:
    try:
        md = smf.mixedlm(formula, cd_model, groups=cd_model["patient_id"])
        mf = md.fit(reml=True, method="lbfgs")
        vp = float(mf.cov_re.iloc[0,0]); vr = float(mf.scale)
        vp_rows.append(dict(model=formula, var_patient=vp, var_residual=vr,
                            ICC=vp/(vp+vr), AIC=mf.aic, BIC=mf.bic))
    except: pass
pd.DataFrame(vp_rows).to_csv(TBL("A9_variance_partitioning.csv"), index=False)

# =========================================================================
# A10: Sequence Similarity Network
# =========================================================================
log("--- A10: Sequence Network ---")
treg_seqs = subset_pools["Treg"]["cdr3aa"].unique().tolist()
if len(treg_seqs) > 3000:
    treg_seqs = list(RNG.choice(treg_seqs, 3000, replace=False))
log(f"  Treg CDR3s sampled: {len(treg_seqs)}")

def l1_variants(s):
    n = len(s)
    return [s[:i] + aa + s[i+1:] for i in range(n) for aa in "ACDEFGHIKLMNPQRSTVWY" if aa != s[i]]

seq_set = set(treg_seqs)
edges = []
for seq in treg_seqs:
    for var in l1_variants(seq):
        if var in seq_set and var != seq:
            edges.append((seq, var))
edges = list(set(edges))
log(f"  Levenshtein-1 edges: {len(edges)}")
pd.DataFrame(edges, columns=["cdr3_1","cdr3_2"]).to_csv(TBL("A10_L1_network_edges.csv"), index=False)

if len(edges):
    import igraph as ig
    all_seqs = list(set([e[0] for e in edges] + [e[1] for e in edges]))
    si = {s:i for i,s in enumerate(all_seqs)}
    g = ig.Graph(n=len(all_seqs))
    g.add_edges([(si[e[0]], si[e[1]]) for e in edges])
    csizes = [len(c) for c in g.components()]
    log(f"  Components: {len(csizes)}, max={max(csizes)}, median={np.median(csizes):.0f}")
    pd.DataFrame({"component_size": csizes}).to_csv(TBL("A10_L1_component_sizes.csv"), index=False)
    fig, ax = plt.subplots(figsize=(7, 5))
    sc, sv = zip(*sorted(Counter(csizes).items()))
    ax.bar(sc, sv, color="steelblue", edgecolor="k")
    ax.set_xlabel("Component size"); ax.set_ylabel("Count")
    ax.set_title("CDR3 Similarity Network: Component Sizes")
    plt.tight_layout(); fig.savefig(FIG("A10_L1_component_sizes.png"), dpi=200); plt.close(fig)

# =========================================================================
# A11: Subset Overlap (within patient)
# =========================================================================
log("--- A11: Subset Overlap ---")
turn_rows = []
for pid in meta_qc["patient_id"].unique():
    pd_ = meta_qc[meta_qc["patient_id"] == pid]
    for s1, s2 in itertools.combinations(SUBSETS, 2):
        ids1 = pd_[pd_.subset==s1]["sample_id"].tolist()
        ids2 = pd_[pd_.subset==s2]["sample_id"].tolist()
        if not ids1 or not ids2: continue
        parts1 = [samples[i]["cdr3aa"] for i in ids1 if i in samples]
        parts2 = [samples[i]["cdr3aa"] for i in ids2 if i in samples]
        if not parts1 or not parts2: continue
        cc1 = set(pd.concat(parts1).unique())
        cc2 = set(pd.concat(parts2).unique())
        if len(cc1) and len(cc2):
            turn_rows.append(dict(patient_id=pid, s1=s1, s2=s2,
                                  shared=len(cc1&cc2), unique_s1=len(cc1-cc2),
                                  unique_s2=len(cc2-cc1),
                                  jaccard=round(len(cc1&cc2)/len(cc1|cc2), 4)))
tdf = pd.DataFrame(turn_rows)
tdf.to_csv(TBL("A11_subset_overlap.csv"), index=False)
log(f"  Subset pairs: {len(tdf)}, mean Jaccard={tdf['jaccard'].mean():.3f}" if len(tdf) else "")

if len(tdf):
    fig, ax = plt.subplots(figsize=(6, 5))
    jmat = tdf.pivot_table(index="s1", columns="s2", values="jaccard", aggfunc="mean")
    jmat = jmat.reindex(index=SUBSETS, columns=SUBSETS).fillna(0)
    im = ax.imshow(jmat, cmap="Blues", vmin=0, vmax=1)
    for i in range(3):
        for j in range(3):
            v = jmat.iloc[i,j]
            if v > 0:
                ax.text(j, i, f"{v:.3f}", ha="center", va="center", fontsize=12,
                       color="white" if v > 0.5 else "black")
    ax.set_xticks(range(3)); ax.set_yticks(range(3))
    ax.set_xticklabels(SUBSETS); ax.set_yticklabels(SUBSETS)
    ax.set_title("Mean Jaccard Overlap"); plt.colorbar(im)
    plt.tight_layout(); fig.savefig(FIG("A11_subset_overlap_heatmap.png"), dpi=200); plt.close(fig)

# =========================================================================
# A12: Jensen-Shannon Divergence
# =========================================================================
log("--- A12: JSD ---")
def jsd(p, q):
    p = np.asarray(p, float); p = p/p.sum()
    q = np.asarray(q, float); q = q/q.sum()
    m = (p + q) / 2
    # handle 0*log(0) = 0
    def kl(x, y):
        mask = x > 0
        return np.sum(x[mask] * np.log(x[mask] / y[mask]))
    val = 0.5 * (kl(p, m) + kl(q, m))
    return round(np.sqrt(max(0, val)), 4)

jsd_res = {}
for s1, s2 in itertools.combinations(SUBSETS, 2):
    keys = set(subset_pools[s1]["cdr3aa"]) | set(subset_pools[s2]["cdr3aa"])
    d1 = dict(zip(subset_pools[s1]["cdr3aa"], subset_pools[s1]["freq"]))
    d2 = dict(zip(subset_pools[s2]["cdr3aa"], subset_pools[s2]["freq"]))
    v1 = np.array([d1.get(k, 0) for k in keys])
    v2 = np.array([d2.get(k, 0) for k in keys])
    jsd_res[f"{s1}_vs_{s2}"] = jsd(v1, v2)
pd.DataFrame([jsd_res]).to_csv(TBL("A12_JSD.csv"), index=False)
log(f"  JSD: {jsd_res}")

# =========================================================================
# SUMMARY
# =========================================================================
summary_path = os.path.join(OUT, "additional_analyses_summary.txt")
lines = [
    "=" * 60,
    "ADDITIONAL ANALYSES SUMMARY",
    f"Data: {len(samples)} samples, {meta_qc['patient_id'].nunique()} patients",
    "=" * 60,
    f"A1: {len(all_v)} V genes, {len(all_j)} J genes. Most diff: {vdiff.iloc[0]['V']} p={vdiff.iloc[0]['p']:.4f}",
    f"A2: Mean CDR3 lengths {cdr3_props_df['mean_len'].to_dict()}",
    f"A3: {len(fd)} clones tested, {len(sig)} FDR<0.05",
    f"A4: {len(odf)} comp. pairs (mean J={odf['jaccard'].mean():.3f})" if len(odf) else "A4: N/A",
    f"A5: PERMDISP F={f_disp:.3f} p={p_disp:.3f}, ANOSIM R={R_ano:.3f} p={p_ano:.3f}",
    f"A6: Top1 per subset: {ddf.groupby('subset')['top1'].mean().to_dict()}",
    f"A7: {len(sharing)} CDR3s, {(sharing>=2).sum()} public, {(sharing==1).sum()} private",
    f"A8: RF AUC={np.mean(aucs):.3f}" if len(bin_df)>10 else "A8: Insufficient data",
    f"A9: Variance partitioning ({len(vp_rows)} models)",
    f"A10: {len(edges)} L1 edges, {len(csizes)} components" if len(edges) else "A10: No edges",
    f"A11: {len(tdf)} subset-overlap pairs" if len(tdf) else "A11: N/A",
    f"A12: JSD = {jsd_res}",
    "=" * 60,
]
with open(summary_path, "w") as f: f.write("\n".join(lines))
log("=== DONE ===")
