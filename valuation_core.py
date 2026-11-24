"""Valuation engine life-table core: rebuild NCHS's U.S. life table, 2019 (Arias & Xu, NVSR 70(19), 2022) from its
printed death probabilities qx and compare life expectancy at birth and at 65 with the printed values.
Conventions fixed before the run: a0 = 0.10 (infant deaths early in the year), ax = 0.5 for ages 1-99,
open interval 100+ closed with a constant hazard mu = -ln(1 - q99). Printed lx/Lx/Tx/ex are NOT used as inputs.
Run: python3 valuation_core.py  -> results/valuation_core.json"""
import json, os, math
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
t = pd.read_excel(os.path.join(HERE, "sources", "lt2019_t1.xlsx"), header=None).iloc[3:104]
q = t[1].astype(float).tolist(); printed_e = t[6].astype(float).tolist()
l = [100000.0]; L = []
for x in range(100):
    d = l[x] * q[x]; a = 0.10 if x == 0 else 0.5
    L.append(l[x] - (1 - a) * d); l.append(l[x] - d)
mu = -math.log(1 - q[99]); L.append(l[100] / mu)
T = [sum(L[x:]) for x in range(101)]; e = [T[x] / l[x] for x in range(101)]
out = {"source_citation": "Arias E, Xu J. United States Life Tables, 2019. National Vital Statistics Reports 70(19), March 22, 2022; Table 1",
       "verified_date": "2026-11-03",
       "e0": {"published_value": 78.8, "published_table_value": round(printed_e[0], 4), "reproduced_value": round(e[0], 4),
              "abs_diff": round(abs(e[0] - printed_e[0]), 4), "pct_diff": round(100 * abs(e[0] - printed_e[0]) / printed_e[0], 3)},
       "e65": {"published_table_value": round(printed_e[65], 4), "reproduced_value": round(e[65], 4),
               "abs_diff": round(abs(e[65] - printed_e[65]), 4)},
       "interval": "n/a — deterministic computation on the printed qx column (101 ages)",
       "conventions": "a0=0.10; ax=0.5 ages 1-99; 100+ closed with constant hazard -ln(1-q99)"}
os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
json.dump(out, open(os.path.join(HERE, "results", "valuation_core.json"), "w"), indent=1); print(json.dumps(out, indent=1))
