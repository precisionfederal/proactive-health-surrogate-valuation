"""Run 3: paired PTE difference, ours (A4 seed-mean) minus PhenoAge, current vs never smoking, validation cycles."""
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
dev = df[df.cycle <= 3].copy(); val = df[df.cycle >= 4].copy(); assert len(val) == 15189
preds = [HistGradientBoostingClassifier(learning_rate=0.05, max_iter=400, max_leaf_nodes=15, l2_regularization=1.0,
         min_samples_leaf=40, early_stopping=True, validation_fraction=0.15, n_iter_no_change=30, random_state=s)
         .fit(dev[panel].values, dev.e10.values).predict_proba(val[panel].values)[:, 1] for s in sb.SEEDS]
np.savez(os.path.join(HERE, "results", "val_preds.npz"), seqn=val.SEQN.values, preds=np.array(preds), phenoage=val.phenoage.values)
smq = pd.concat([sb.rd(f)[["SEQN", "SMQ020", "SMQ040"]] for f in ["SMQ", "SMQ_B", "SMQ_C", "SMQ_D", "SMQ_E", "SMQ_F"]])
v = val.merge(smq, on="SEQN", how="left"); pm = np.clip(np.mean(preds, 0), 1e-6, 1 - 1e-6)
v["lo_ours"] = np.log(pm) - np.log(1 - pm)
cur = (v.SMQ020 == 1) & v.SMQ040.isin([1, 2]); nev = (v.SMQ020 == 2)
v = v[cur | nev].copy(); v["smoke"] = cur[cur | nev].astype(float).values
v["pa_z"] = (v.phenoage - v.phenoage.mean()) / v.phenoage.std(); v["so_z"] = (v.lo_ours - v.lo_ours.mean()) / v.lo_ours.std()
def pte(d, sur):
    t = d.t10.values; e = d.e10.values
    X1 = np.column_stack([d.smoke, d.RIDAGEYR, d.sex]).astype(float); X1 -= X1.mean(0)
    X2 = np.column_stack([d.smoke, d.RIDAGEYR, d.sex, d[sur]]).astype(float); X2 -= X2.mean(0)
    b1, _ = sb.cox_hr(X1, t, e); b2, _ = sb.cox_hr(X2, t, e); return float(1 - b2[0] / b1[0])
pO, pP = pte(v, "so_z"), pte(v, "pa_z"); r = np.random.default_rng(7); diffs = []
for _ in range(200):
    ix = r.integers(0, len(v), len(v)); d = v.iloc[ix]; diffs.append(pte(d, "so_z") - pte(d, "pa_z"))
lo, hi = np.percentile(diffs, [2.5, 97.5])
out = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"), "pte_ours": pO, "pte_phenoage": pP, "delta": pO - pP,
       "ci95_low": float(lo), "ci95_high": float(hi), "n": int(len(v)), "n_boot": 200, "seeds": sb.SEEDS,
       "method": "paired bootstrap (identical resamples) of Freedman PTE difference; Cox 10-year, age and sex adjusted",
       "separates": bool(lo > 0 or hi < 0), "runtime_s": round(time.time() - t0, 1)}
json.dump(out, open(os.path.join(HERE, "results", "pte_paired.json"), "w"), indent=1); print(json.dumps(out, indent=1))
