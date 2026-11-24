"""Valuation engine dollar core: reproduce HHS/ASPE Appendix D (April 2021), "Updating Value per Statistical Life (VSL)
Estimates for Inflation and Changes in Real Income", from its printed inputs only, and compare every printed value.
Declared before the run in BASELINE-REGISTRY.md (append of 2026-11-03). Printed tables are parsed from the source text
(sources/aspe-2021-appendix-d-vsl-update.txt, pdftotext -layout), never retyped.
Run: python3 hhs_valuation_repro.py  -> results/hhs_valuation_repro.json"""
import json, os, re
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "sources", "aspe-2021-appendix-d-vsl-update.txt")
txt = open(SRC, encoding="utf-8").read()

# ---- printed inputs (Appendix D text, D.3.1 bullets, D.3.2, D.4) ----
VSL13 = {"low": 4.2e6, "central": 9.0e6, "high": 13.7e6}
CPI13, CPI20 = 232.957, 258.811
EARN13, EARN20 = 333.0, 380.0          # "about $333 to $380 per week"
ELAST, G = 1.0, 0.008                   # elasticity 1.0; CBO 2021 real earnings growth 0.8%/yr
PVQ40 = {0.03: 19.6, 0.07: 11.8}        # PV future QALYs, age 40 (printed, rounded)
PVQ75_3 = 7.94
for s in ["232.957", "258.811", "$333 to $380", "14.1 percent", "0.8 percent per year", "19.6 QALYs", "11.8 QALYs", "7.94 QALYs"]:
    assert s in txt, s

vsl20 = {k: v * (CPI20 / CPI13) * (EARN20 / EARN13) ** ELAST for k, v in VSL13.items()}
def vsl(k, y): return vsl20[k] * (1 + G) ** (ELAST * (y - 2020))
def vq(k, y, r): return vsl(k, y) / PVQ40[r]

def m(x): return f"${x/1e6:.1f} million"
def q(x):
    v = int(round(x / 1e4)) * 10000
    return f"${v:,}"
def b(x): return f"${x/1e9:.2f} billion"

rows = []

# ---- Table D.1 and D.2 (parsed) ----
results = {}
def check(table, label, ours_str, printed_str, ours_raw):
    rows.append({"table": table, "item": label, "ours_raw": round(ours_raw, 2), "ours_printed_form": ours_str,
                 "printed": printed_str, "match": ours_str == printed_str})

a = txt.index("Table D.2. VSL estimates"); seg = txt[a:txt.index("D.3.2 Constant value", a)]
for line in seg.splitlines():
    mm = re.match(r"\s*(20[2-4]\d)\s+(\$[\d.]+ million)\s+(\$[\d.]+ million)\s+(\$[\d.]+ million)", line)
    if mm:
        y = int(mm.group(1))
        for k, p in zip(["low", "central", "high"], mm.groups()[1:]):
            check("D.2" if y > 2020 else "D.1/D.2", f"VSL {k} {y}", m(vsl(k, y)), p, vsl(k, y))

# ---- Table D.3 / D.4a / D.4b (parsed) ----
for tag, r, start, end in [("D.4a", 0.03, "Table D.4a.", "Table D.4b."), ("D.4b", 0.07, "Table D.4b.", "D.4 Examples")]:
    a = txt.index(start); seg = txt[a:txt.index(end, a)]
    for line in seg.splitlines():
        mm = re.match(r"\s*(20[2-4]\d)\s+(\$[\d,]+)\s+(\$[\d,]+)\s+(\$[\d,]+)", line)
        if mm:
            y = int(mm.group(1))
            for k, p in zip(["low", "central", "high"], mm.groups()[1:]):
                check(tag if y > 2020 else f"D.3/{tag}", f"value per QALY {k} {y} at {int(r*100)}%", q(vq(k, y, r)), p, vq(k, y, r))

# ---- Examples D.5-D.9 (3%, base year 2023, initial year undiscounted) ----
yrs = range(2023, 2028)
pv = lambda vals: sum(v / 1.03 ** i for i, v in enumerate(vals))
d5 = pv([100 * vsl("central", y) for y in yrs])
d6 = pv([100 * vsl("central", y) / 1.03 ** 3 for y in yrs])
d7 = pv([100 * PVQ75_3 * vq("central", y, 0.03) for y in yrs])
d8 = pv([50 * 0.2 * vq("central", y, 0.03) for y in yrs])
d9 = d5 + d8
ex_printed = {"D.5": "$5.60 billion", "D.6": "$5.12 billion", "D.7": "$2.25 billion", "D.8": "$28.3 million", "D.9": "$5.63 billion"}
for tag, ours, form in [("D.5", d5, b), ("D.6", d6, b), ("D.7", d7, b), ("D.8", d8, m), ("D.9", d9, b)]:
    assert ex_printed[tag] in txt, ex_printed[tag]
    check(tag, f"example net present value {tag}", form(ours), ex_printed[tag], ours)

n = len(rows); ok = sum(r["match"] for r in rows)
out = {"source_citation": "U.S. Department of Health and Human Services, ASPE. Appendix D: Updating Value per Statistical Life (VSL) "
                          "Estimates for Inflation and Changes in Real Income. April 2021. Supplement to the HHS Guidelines for Regulatory Impact Analysis (2016).",
       "source_file": "sources/aspe-2021-appendix-d-vsl-update.pdf", "verified_date": "2026-11-03",
       "interval": "n/a — deterministic arithmetic on printed inputs",
       "headline": {"vsl_central_2020": {"published_value": "$11.4 million", "reproduced_value": round(vsl20["central"], 0)},
                    "value_per_qaly_central_2020_3pct": {"published_value": "$580,000", "reproduced_value": round(vq("central", 2020, 0.03), 0)},
                    "example_D5_npv": {"published_value": "$5.60 billion", "reproduced_value": round(d5, 0)}},
       "values_compared": n, "values_matched": ok,
       "mismatches": [r for r in rows if not r["match"]], "rows": rows}
os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
json.dump(out, open(os.path.join(HERE, "results", "hhs_valuation_repro.json"), "w"), indent=1)
print(json.dumps({k: out[k] for k in ["headline", "values_compared", "values_matched", "mismatches"]}, indent=1))
