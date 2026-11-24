# Baseline registry (written before the first run; append only)

Written 2026-11-03, before surrogate_bench.py was run for the first time.

## Cohort (declared)
NHANES 1999–2010 (six continuous cycles, public-use files from wwwn.cdc.gov), linked to the NCHS 2019 public-use
Linked Mortality Files (ftp.cdc.gov, follow-up through 2019-12-31). Adults 20–84 at exam, eligible for linkage,
complete on the nine Phenotypic Age biomarkers, CRP > 0. Outcome: all-cause death within 120 months of the exam
(permth_exm), censored at 120 months.

## Split (declared)
Temporal: development = cycles 1999–2000, 2001–2002, 2003–2004; validation = 2005–2006, 2007–2008, 2009–2010.
Nothing fitted touches the validation cycles.

## Metric (declared)
Headline: Harrell's concordance (C) for 10-year all-cause mortality on the validation cycles. Secondary: AUC for
death within 10 years among validation participants with known 10-year status. Delta vs baseline with a paired
bootstrap 95% interval (1,000 resamples of validation participants). Unweighted.

## Arms (declared before any result)
| Arm | What | Why it is here |
|---|---|---|
| A0 age | chronological age alone | the no-mechanism arm |
| A1 PhenoAge | Levine 2018 / Liu 2018 published formula, coefficients fixed, no refit | THE published surrogate to beat |
| A2 logit-panel | logistic regression (scikit-learn, standardized) on the full routine panel below | the stock-library arm anyone can fit |
| A3 GBM-9 | gradient boosting on PhenoAge's own nine biomarkers + age + sex | same inputs, our learner (isolates method) |
| A4 GBM-panel (HEADLINE) | gradient boosting on the routine panel: every chemistry-panel, CBC and CRP analyte reported in all six cycles + age + sex | our surrogate |

Seeds: 5 for every stochastic arm; the headline is the mean across seeds; per-seed values ship beside it.
Verb rule: "beats" only if the 95% interval on (A4 − A1) excludes zero; otherwise TIE and the S5 ladder runs.

## Published anchors to reproduce (declared)
- Liu et al., PLoS Med 2018;15(12):e1002718, Fig 7: 10-year mortality AUC 0.88 for Phenotypic Age, 0.86 for
  chronological age (NHANES IV 1999–2010, n = 11,432, deaths through 2011-12-31).
- Liu et al. 2018 Table 1: HR 1.09 (95% CI 1.08–1.10) per 1-year Phenotypic Age, adjusted for chronological age and
  sex (Gompertz). We fit Cox (Breslow); stated as a method difference.
Known deviations from Liu: 2019 linkage (longer follow-up), no fasting-hours exclusion, unweighted.

## APPEND 2026-11-03 (after run 1, before run 2) — declared before surrogate_followup.py runs
Run 1 bootstrapped the five-seed ENSEMBLE prediction while the point delta was the MEAN of per-seed C, so point and
interval described two different objects. Run 2 recomputes the headline interval on the declared quantity: per
resample, mean over the five seeds of C(A4, seed) minus C(A1). Same seeds, same RNG seed, same 1,000 resamples.
Also declared now: an individual-level surrogacy statistic on observational data, the proportion of a modifiable
exposure's mortality hazard explained (Freedman, Graubard and Schatzkin, Stat Med 1992). Exposure: current smoking
(SMQ020 = 1 and SMQ040 in 1, 2) vs never smoking (SMQ020 = 2); former smokers excluded. Validation cycles only.
Cox (Breslow) on 10-year mortality: beta1 from [smoking, age, sex]; beta2 adding the surrogate's log-odds;
PTE = 1 − beta2/beta1, 200-resample bootstrap 95% interval. Computed for A1 PhenoAge (standardized) and A4 (seed-mean
log-odds). It is an observational individual-level statistic, never written as trial-level surrogacy.

## APPEND 2026-11-03 (after run 2, before run 3)
Run 3 adds one declared quantity: the paired bootstrap difference PTE(A4) − PTE(A1) on identical resamples (200,
seed 7), so a comparison may be printed only if its 95% interval excludes zero. Saves validation predictions to
results/val_preds.npz so no later step refits.

## APPEND 2026-11-03 (before hhs_valuation_repro.py runs) — the valuation engine's dollar core vs HHS
Published evaluation to reproduce: HHS/ASPE, "Appendix D: Updating Value per Statistical Life (VSL) Estimates for
Inflation and Changes in Real Income", April 2021, supplement to the HHS 2016 Guidelines for Regulatory Impact
Analysis (sources/aspe-2021-appendix-d-vsl-update.pdf, downloaded 2026-11-03 from aspe.hhs.gov).
Inputs, ONLY as printed in Appendix D: VSL 2013 low/central/high $4.2M/$9.0M/$13.7M; CPI-U 232.957 (2013) and
258.811 (2020); real weekly earnings "about $333 to $380" (14.1 percent); income elasticity 1.0; real income growth
0.8 percent per year (CBO 2021); discount rates 3 and 7 percent, no discounting in the initial year; present value of
future QALYs at age 40 = 19.6 (3%) and 11.8 (7%), at age 75 = 7.94 (3%); example regulation effective 2023-2027,
100 deaths averted per year, 50 nonfatal cases at 0.2 QALY each, 3-year cessation lag.
Quantities declared: Table D.1 (3 values), D.2 (90), D.3 (6), D.4a and D.4b (180), and the five example totals
D.5 NPV $5.60 billion, D.6 $5.12 billion, D.7 $2.25 billion, D.8 $28.3 million, D.9 $5.63 billion.
Match rule, fixed now: a value matches when our unrounded result, rounded the way the table prints it (0.1 million
for VSL; $10,000 for value per QALY; 0.01 billion; 0.1 million), equals the printed string. Every mismatch is
reported with both values; none is tuned away. Deterministic: interval "n/a — deterministic arithmetic".
The QALY denominators (19.6, 11.8, 7.94) are taken as printed; their HRQL input (Hanmer 2006 EQ-5D) is not openly
available, so this run reproduces the valuation arithmetic, not the QALY denominator.
