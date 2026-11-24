"""Surrogate bench (declared in BASELINE-REGISTRY.md before the first run).

NHANES 1999-2010 + NCHS 2019 public-use Linked Mortality Files. Temporal split: develop on 1999-2004,
validate on 2005-2010. Arms A0-A4 as registered. Writes results/surrogate_bench.json.
Run: python3 surrogate_bench.py
"""
import json, os, sys, time
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)

CYC = [  # (SDDSRVYR, demo, chem, cbc, crp, mort)
    (1, "DEMO", "LAB18", "LAB25", "LAB11", "MORT_1999_2000"),
    (2, "DEMO_B", "L40_B", "L25_B", "L11_B", "MORT_2001_2002"),
    (3, "DEMO_C", "L40_C", "L25_C", "L11_C", "MORT_2003_2004"),
    (4, "DEMO_D", "BIOPRO_D", "CBC_D", "CRP_D", "MORT_2005_2006"),
    (5, "DEMO_E", "BIOPRO_E", "CBC_E", "CRP_E", "MORT_2007_2008"),
    (6, "DEMO_F", "BIOPRO_F", "CBC_F", "CRP_F", "MORT_2009_2010"),
]
SEEDS = [11, 23, 37, 41, 59]
NBOOT = 1000
HORIZON = 120  # months


def rd(name):
    return pd.read_sas(os.path.join(D, name + ".xpt"))


def mort(name):
    rows = []
    with open(os.path.join(D, name + ".dat")) as f:
        for line in f:
            line = line.rstrip("\n")
            if len(line) < 48:
                line = line.ljust(48)
            seqn = line[0:6].strip(); el = line[14:15].strip(); ms = line[15:16].strip(); pm = line[45:48].strip()
            rows.append((int(seqn), el, ms, pm))
    m = pd.DataFrame(rows, columns=["SEQN", "eligstat", "mortstat", "permth_exm"])
    m = m[m.eligstat == "1"].copy()
    m["mortstat"] = pd.to_numeric(m.mortstat, errors="coerce")
    m["permth_exm"] = pd.to_numeric(m.permth_exm, errors="coerce")
    return m.dropna(subset=["mortstat", "permth_exm"])


def load():
    chem_cols = None; cbc_cols = None; parts = []
    chems = {c[0]: rd(c[2]) for c in CYC}; cbcs = {c[0]: rd(c[3]) for c in CYC}
    # 2001-2002 files name some analytes LBD... where every other cycle uses LBX...; align to the 1999 names
    for group in (chems, cbcs):
        ref = set(group[1].columns)
        for k in group:
            ren = {c: "LBX" + c[3:] for c in group[k].columns
                   if c.startswith("LBD") and ("LBX" + c[3:]) in ref and ("LBX" + c[3:]) not in group[k].columns
                   and c not in ref}
            group[k] = group[k].rename(columns=ren)
    for k in chems:
        chem_cols = set(chems[k].columns) if chem_cols is None else chem_cols & set(chems[k].columns)
        cbc_cols = set(cbcs[k].columns) if cbc_cols is None else cbc_cols & set(cbcs[k].columns)
    # one unit per analyte: drop the conventional column where an SI twin exists
    chem_cols = sorted(c for c in chem_cols if c != "SEQN" and not (c.startswith("LBX") and ("LBD" + c[3:] + "SI") in chem_cols))
    cbc_cols = sorted(c for c in cbc_cols if c != "SEQN")
    for (cy, demo, chem, cbc, crp, mo) in CYC:
        d = rd(demo)[["SEQN", "SDDSRVYR", "RIAGENDR", "RIDAGEYR"]]
        x = d.merge(chems[cy][["SEQN"] + chem_cols], on="SEQN", how="left")
        x = x.merge(cbcs[cy][["SEQN"] + cbc_cols], on="SEQN", how="left")
        x = x.merge(rd(crp)[["SEQN", "LBXCRP"]], on="SEQN", how="left")
        x = x.merge(mort(mo), on="SEQN", how="inner")
        x["cycle"] = cy
        parts.append(x)
    df = pd.concat(parts, ignore_index=True)
    return df, chem_cols, cbc_cols


def phenoage(df):
    # Levine 2018 Table 1 / Liu 2018 equation; full-precision constants as published in the BioAge package
    xb = (-19.90667 - 0.03359355 * df.LBDSALSI + 0.009506491 * df.LBDSCRSI + 0.1953192 * df.LBDSGLSI
          + 0.09536762 * np.log(df.LBXCRP) - 0.01199984 * df.LBXLYPCT + 0.02676401 * df.LBXMCVSI
          + 0.3306156 * df.LBXRDW + 0.001868778 * df.LBXSAPSI + 0.05542406 * df.LBXWBCSI + 0.08035356 * df.RIDAGEYR)
    m = 1 - np.exp((-1.51714 * np.exp(xb)) / 0.007692696)
    return (np.log(-0.0055305 * np.log(1 - m)) / 0.090165) + 141.50225


def cox_hr(X, t, e):
    """Cox partial likelihood, Breslow ties, Newton-Raphson. Returns beta, se."""
    o = np.argsort(-t, kind="mergesort"); X = X[o]; t = t[o]; e = e[o]
    b = np.zeros(X.shape[1])
    for _ in range(50):
        eta = X @ b; w = np.exp(eta - eta.max())
        S0 = np.cumsum(w); S1 = np.cumsum(w[:, None] * X, 0)
        S2 = np.cumsum(w[:, None, None] * X[:, :, None] * X[:, None, :], 0)
        # risk set for time t_i = everyone with t >= t_i -> last index of tie block
        idx = np.searchsorted(-t, -t, side="right") - 1
        ev = np.where(e == 1)[0]; j = idx[ev]
        mu = S1[j] / S0[j][:, None]
        g = X[ev].sum(0) - mu.sum(0)
        H = -((S2[j] / S0[j][:, None, None]).sum(0) - np.einsum('ni,nj->ij', mu, mu))
        step = np.linalg.solve(H, g); b = b - step
        if np.max(np.abs(step)) < 1e-9:
            break
    se = np.sqrt(np.diag(np.linalg.inv(-H)))
    return b, se


def cmats(t, e, r):
    ev = np.where((e == 1) & (t <= HORIZON))[0]
    ti = t[ev][:, None]; tj = t[None, :]
    comp = (tj > ti) | ((tj == ti) & (e[None, :] == 0))
    ri = r[ev][:, None]; rj = r[None, :]
    conc = np.where(comp, (ri > rj) + 0.5 * (ri == rj), 0.0).astype(np.float32)
    return ev, comp.astype(np.float32), conc


def cindex(ev, comp, conc, w=None):
    if w is None:
        return float(conc.sum() / comp.sum())
    we = w[ev]
    return float(we @ (conc @ w) / (we @ (comp @ w)))


def main():
    t0 = time.time()
    cache = os.path.join(D, "_cohort_cache.pkl")
    if os.path.exists(cache):
        df, chem_cols, cbc_cols = pd.read_pickle(cache)
    else:
        df, chem_cols, cbc_cols = load(); pd.to_pickle((df, chem_cols, cbc_cols), cache)
    need = ["LBDSALSI", "LBDSCRSI", "LBDSGLSI", "LBXCRP", "LBXLYPCT", "LBXMCVSI", "LBXRDW", "LBXSAPSI", "LBXWBCSI"]
    n_all = len(df)
    df = df[(df.RIDAGEYR >= 20) & (df.RIDAGEYR <= 84)]
    df = df.dropna(subset=need)
    df = df[df.LBXCRP > 0].copy()
    df["sex"] = (df.RIAGENDR == 2).astype(float)
    df["phenoage"] = phenoage(df)
    df = df[np.isfinite(df.phenoage)].copy()
    df["t10"] = np.minimum(df.permth_exm, HORIZON).astype(float)
    df["e10"] = ((df.mortstat == 1) & (df.permth_exm <= HORIZON)).astype(int)
    df["known10"] = (df.e10 == 1) | (df.permth_exm >= HORIZON)

    dev = df[df.cycle <= 3].copy(); val = df[df.cycle >= 4].copy()
    nine = need[:-1] + ["LBXWBCSI"]
    feats9 = ["RIDAGEYR", "sex"] + need
    panel = ["RIDAGEYR", "sex", "LBXCRP"] + [c for c in chem_cols] + [c for c in cbc_cols]
    panel = list(dict.fromkeys(panel))

    # ---------------- reproduction anchor (whole 1999-2010 cohort) ----------------
    k = df[df.known10]
    auc_pheno_all = roc_auc_score(k.e10, k.phenoage); auc_age_all = roc_auc_score(k.e10, k.RIDAGEYR)
    X = np.column_stack([df.phenoage.values, df.RIDAGEYR.values, df.sex.values]).astype(float)
    X = X - X.mean(0)
    b, se = cox_hr(X, df.permth_exm.values.astype(float), (df.mortstat.values == 1).astype(int))
    hr = float(np.exp(b[0])); hr_lo = float(np.exp(b[0] - 1.96 * se[0])); hr_hi = float(np.exp(b[0] + 1.96 * se[0]))

    # ---------------- arms on the validation cycles ----------------
    ydev = dev.e10.values  # every dev participant has >=15 y of follow-up through 2019
    assert dev.known10.all()
    risks = {"A0_age": val.RIDAGEYR.values.astype(float), "A1_phenoage": val.phenoage.values.astype(float)}
    med = dev[panel].median()
    Xd = dev[panel].fillna(med).values; Xv = val[panel].fillna(med).values
    mu, sd = Xd.mean(0), Xd.std(0) + 1e-9
    lr = LogisticRegression(max_iter=5000, C=1.0).fit((Xd - mu) / sd, ydev)
    risks["A2_logit_panel"] = lr.predict_proba((Xv - mu) / sd)[:, 1]
    per_seed = {"A3_gbm9": [], "A4_gbm_panel": []}
    preds = {"A3_gbm9": [], "A4_gbm_panel": []}
    for s in SEEDS:
        for arm, cols in (("A3_gbm9", feats9), ("A4_gbm_panel", panel)):
            m = HistGradientBoostingClassifier(learning_rate=0.05, max_iter=400, max_leaf_nodes=15, l2_regularization=1.0,
                                               min_samples_leaf=40, early_stopping=True, validation_fraction=0.15,
                                               n_iter_no_change=30, random_state=s)
            m.fit(dev[cols].values, ydev)
            preds[arm].append(m.predict_proba(val[cols].values)[:, 1])
    tv = val.t10.values; evv = val.e10.values
    mats = {}
    for a in list(risks) + list(preds):
        pass
    res_arms = {}
    for a, r in risks.items():
        mats[a] = cmats(tv, evv, r)
    for a in preds:
        for i, p in enumerate(preds[a]):
            per_seed[a].append(cindex(*cmats(tv, evv, p)))
        risks[a] = np.mean(preds[a], axis=0)  # seed-mean ensemble prediction for the bootstrap
        mats[a] = cmats(tv, evv, risks[a])
    kv = val.known10.values
    for a, r in risks.items():
        res_arms[a] = {"c_index": cindex(*mats[a]), "auc10_known": float(roc_auc_score(evv[kv], r[kv]))}
    for a in per_seed:
        res_arms[a]["c_index_per_seed"] = per_seed[a]
        res_arms[a]["c_index_seed_mean"] = float(np.mean(per_seed[a]))
        assert abs(res_arms[a]["c_index_seed_mean"] - np.mean(per_seed[a])) < 1e-12

    rng = np.random.default_rng(20261103); n = len(val)
    deltas = {f"{a}_minus_A1_phenoage": [] for a in ["A0_age", "A2_logit_panel", "A3_gbm9", "A4_gbm_panel"]}
    for _ in range(NBOOT):
        w = rng.multinomial(n, np.ones(n) / n).astype(np.float32)
        ev0, comp0, _ = mats["A1_phenoage"]; we = w[ev0]; den = we @ (comp0 @ w)
        cw = {a: float(we @ (mats[a][2] @ w) / den) for a in mats}
        for a in ["A0_age", "A2_logit_panel", "A3_gbm9", "A4_gbm_panel"]:
            deltas[f"{a}_minus_A1_phenoage"].append(cw[a] - cw["A1_phenoage"])
    dres = {}
    for kname, v in deltas.items():
        a = kname.split("_minus_")[0]
        if a in per_seed:
            point = res_arms[a]["c_index_seed_mean"] - res_arms["A1_phenoage"]["c_index"]
        else:
            point = res_arms[a]["c_index"] - res_arms["A1_phenoage"]["c_index"]
        lo, hi = np.percentile(v, [2.5, 97.5])
        dres[kname] = {"delta": float(point), "ci95_low": float(lo), "ci95_high": float(hi), "n": int(n),
                       "events_10y": int(evv.sum()), "seeds": SEEDS if a in per_seed else "n/a — deterministic arm",
                       "method": f"paired bootstrap, {NBOOT} resamples of validation participants; Harrell C, 10-year horizon",
                       "separates": bool(lo > 0 or hi < 0)}
    out = {
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "data": {"nhanes_cycles": "1999-2010 (6 continuous cycles)", "linkage": "NCHS 2019 public-use Linked Mortality Files, follow-up through 2019-12-31",
                 "n_records_linked_eligible": int(n_all), "n_cohort": int(len(df)), "deaths_total_2019": int((df.mortstat == 1).sum()),
                 "deaths_within_10y": int(df.e10.sum()), "n_dev": int(len(dev)), "n_val": int(len(val)),
                 "dev_cycles": "1999-2004", "val_cycles": "2005-2010", "panel_features": panel, "n_panel_features": len(panel)},
        "reproduction": {
            "phenoage_auc10": {"published_value": 0.88, "reproduced_value": round(float(auc_pheno_all), 4),
                               "abs_diff": round(abs(float(auc_pheno_all) - 0.88), 4), "pct_diff": round(100 * abs(float(auc_pheno_all) - 0.88) / 0.88, 2),
                               "n": int(len(k)), "source_citation": "Liu Z et al. PLoS Med 2018;15(12):e1002718, Fig 7",
                               "verified_date": "2026-11-03"},
            "age_auc10": {"published_value": 0.86, "reproduced_value": round(float(auc_age_all), 4),
                          "abs_diff": round(abs(float(auc_age_all) - 0.86), 4), "source_citation": "Liu 2018 Fig 7", "verified_date": "2026-11-03"},
            "phenoage_hr_per_year_adj_age_sex": {"published_value": 1.09, "published_ci": [1.08, 1.10], "reproduced_value": round(hr, 4),
                                                 "reproduced_ci95": [round(hr_lo, 4), round(hr_hi, 4)], "abs_diff": round(abs(hr - 1.09), 4),
                                                 "method": "Cox (Breslow) on 2019 follow-up; Liu used Gompertz on follow-up through 2011",
                                                 "source_citation": "Liu 2018 Table 1", "verified_date": "2026-11-03"}},
        "arms_validation": res_arms,
        "deltas_vs_phenoage": dres,
        "runtime_s": round(time.time() - t0, 1),
    }
    with open(os.path.join(OUT, "surrogate_bench.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: out[k] for k in ["data", "reproduction", "deltas_vs_phenoage"]}, indent=1)[:6000])
    print({a: {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in r.items()} for a, r in res_arms.items()})


if __name__ == "__main__":
    main()
