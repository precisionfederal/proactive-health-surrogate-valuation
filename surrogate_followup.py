"""Run 2 (declared in BASELINE-REGISTRY APPEND 2026-11-03): seed-mean C interval for A4 vs A1, and Freedman PTE
for current smoking. Reuses surrogate_bench.py functions and its cohort cache. Writes results/surrogate_followup.json"""
import json, os, sys, time
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import surrogate_bench as sb
from sklearn.ensemble import HistGradientBoostingClassifier

t0 = time.time()
df, chem_cols, cbc_cols = pd.read_pickle(os.path.join(sb.D, "_cohort_cache.pkl"))
need = ["LBDSALSI", "LBDSCRSI", "LBDSGLSI", "LBXCRP", "LBXLYPCT", "LBXMCVSI", "LBXRDW", "LBXSAPSI", "LBXWBCSI"]
df = df[(df.RIDAGEYR >= 20) & (df.RIDAGEYR <= 84)].dropna(subset=need)
df = df[df.LBXCRP > 0].copy(); df["sex"] = (df.RIAGENDR == 2).astype(float); df["phenoage"] = sb.phenoage(df)
df = df[np.isfinite(df.phenoage)].copy()
df["t10"] = np.minimum(df.permth_exm, sb.HORIZON).astype(float)
df["e10"] = ((df.mortstat == 1) & (df.permth_exm <= sb.HORIZON)).astype(int)
panel = list(dict.fromkeys(["RIDAGEYR", "sex", "LBXCRP"] + list(chem_cols) + list(cbc_cols)))
dev = df[df.cycle <= 3].copy(); val = df[df.cycle >= 4].copy()
assert len(val) == 15189, len(val)
preds = []
for s in sb.SEEDS:
    m = HistGradientBoostingClassifier(learning_rate=0.05, max_iter=400, max_leaf_nodes=15, l2_regularization=1.0,
                                       min_samples_leaf=40, early_stopping=True, validation_fraction=0.15,
                                       n_iter_no_change=30, random_state=s).fit(dev[panel].values, dev.e10.values)
    preds.append(m.predict_proba(val[panel].values)[:, 1])
tv = val.t10.values; ev = val.e10.values
ev0, comp, concP = sb.cmats(tv, ev, val.phenoage.values.astype(float))
concS = [sb.cmats(tv, ev, p)[2] for p in preds]
cP = sb.cindex(ev0, comp, concP); cS = [sb.cindex(ev0, comp, c) for c in concS]
point = float(np.mean(cS) - cP)
rng = np.random.default_rng(20261103); n = len(val); boot = []
for _ in range(sb.NBOOT):
    w = rng.multinomial(n, np.ones(n) / n).astype(np.float32); we = w[ev0]; den = we @ (comp @ w)
    base = we @ (concP @ w) / den
    boot.append(float(np.mean([we @ (c @ w) / den for c in concS]) - base))
lo, hi = np.percentile(boot, [2.5, 97.5])
del concS, concP, comp

# ---- Freedman PTE for current smoking (validation cycles) ----
smq = pd.concat([sb.rd(f)[["SEQN", "SMQ020", "SMQ040"]] for f in ["SMQ", "SMQ_B", "SMQ_C", "SMQ_D", "SMQ_E", "SMQ_F"]])
v = val.merge(smq, on="SEQN", how="left")
v["lo_ours"] = np.log(np.clip(np.mean(preds, 0), 1e-6, 1 - 1e-6))
v["lo_ours"] = np.log(np.clip(np.mean(preds, 0), 1e-6, 1 - 1e-6)) - np.log(1 - np.clip(np.mean(preds, 0), 1e-6, 1 - 1e-6))
cur = (v.SMQ020 == 1) & v.SMQ040.isin([1, 2]); nev = (v.SMQ020 == 2)
v = v[cur | nev].copy(); v["smoke"] = ((v.SMQ020 == 1) & v.SMQ040.isin([1, 2])).astype(float)
v["pa_z"] = (v.phenoage - v.phenoage.mean()) / v.phenoage.std(); v["so_z"] = (v.lo_ours - v.lo_ours.mean()) / v.lo_ours.std()
def pte(d, sur):
    t = d.t10.values; e = d.e10.values
    X1 = np.column_stack([d.smoke, d.RIDAGEYR, d.sex]).astype(float); X1 -= X1.mean(0)
    X2 = np.column_stack([d.smoke, d.RIDAGEYR, d.sex, d[sur]]).astype(float); X2 -= X2.mean(0)
    b1, _ = sb.cox_hr(X1, t, e); b2, _ = sb.cox_hr(X2, t, e)
    return float(b1[0]), float(b2[0]), float(1 - b2[0] / b1[0])
res = {}
for sur in ["pa_z", "so_z"]:
    b1, b2, p = pte(v, sur); bs = []
    r2 = np.random.default_rng(7)
    for _ in range(200):
        ix = r2.integers(0, len(v), len(v)); bs.append(pte(v.iloc[ix], sur)[2])
    l2, h2 = np.percentile(bs, [2.5, 97.5])
    res[sur] = {"beta_smoking_without_surrogate": b1, "beta_smoking_with_surrogate": b2, "hr_smoking_without": float(np.exp(b1)),
                "hr_smoking_with": float(np.exp(b2)), "pte": p, "ci95_low": float(l2), "ci95_high": float(h2), "n_boot": 200}
d_pte = [None]
out = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"),
       "headline_seed_mean": {"c_phenoage": cP, "c_ours_per_seed": cS, "c_ours_seed_mean": float(np.mean(cS)),
                              "delta": point, "ci95_low": float(lo), "ci95_high": float(hi), "n": int(n), "events_10y": int(ev.sum()),
                              "seeds": sb.SEEDS, "method": "paired bootstrap, 1,000 resamples; per resample mean over 5 seeds of C(ours) minus C(PhenoAge); Harrell C, 10-year horizon",
                              "separates": bool(lo > 0)},
       "pte_current_smoking": {"n": int(len(v)), "n_current": int(v.smoke.sum()), "deaths_10y": int(v.e10.sum()),
                               "method": "Freedman, Graubard, Schatzkin, Stat Med 1992;11(2):167-178; Cox (Breslow), 10-year horizon, adjusted for age and sex; observational, individual-level",
                               "phenoage": res["pa_z"], "ours": res["so_z"]},
       "runtime_s": round(time.time() - t0, 1)}
json.dump(out, open(os.path.join(HERE, "results", "surrogate_followup.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
