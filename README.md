# Proactive health surrogate and valuation core: measured on public data

Apache License 2.0. Precision Federal, 2026.

What this repository reproduces, from public files only:

- A ten-year all-cause mortality surrogate learned from the routine blood panel (chemistry panel, complete blood
  count, CRP, age, sex) on NHANES 1999–2004 and tested on NHANES 2005–2010, linked to the NCHS 2019 public-use
  Linked Mortality Files.
- Phenotypic Age (Levine et al., Aging 2018) computed with its published coefficients and compared on the same
  held-out adults; our code reproduces the ten-year mortality AUC printed by Liu et al., PLoS Medicine 2018.
- The proportion of current smoking's ten-year mortality hazard explained by each score (Freedman, Graubard and
  Schatzkin, Statistics in Medicine 1992).
- The valuation engine's life-table core, which rebuilds the NCHS United States Life Tables, 2019 (Arias and Xu,
  NVSR 70(19), 2022) from the printed death probabilities alone.
- The valuation engine's dollar core, which reproduces HHS's published valuation of mortality risk (ASPE, Appendix D:
  Updating Value per Statistical Life Estimates, April 2021) from the inputs that document prints, and reports every
  printed value it does and does not match.

## Data (downloaded by you, never redistributed here)

- NHANES public-use files, 1999–2010: DEMO, chemistry panel (LAB18, L40_B, L40_C, BIOPRO_D–F), CBC (LAB25, L25_B,
  L25_C, CBC_D–F), CRP (LAB11, L11_B, L11_C, CRP_D–F), smoking (SMQ, SMQ_B–F), from https://wwwn.cdc.gov/nchs/nhanes/
- NCHS 2019 public-use Linked Mortality Files, NHANES 1999–2010, from
  https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/
- NVSR 70(19) Table 1, from https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/NVSR/70-19/Table01.xlsx
- ASPE Appendix D (April 2021), from https://aspe.hhs.gov/reports/updating-vsl-estimates; save its text with
  `pdftotext -layout` as sources/aspe-2021-appendix-d-vsl-update.txt

Place the NHANES .xpt and mortality .dat files in data/.

## Run

    python3 surrogate_bench.py      # arms, reproduction, bootstrap intervals -> results/surrogate_bench.json
    python3 surrogate_followup.py   # seed-mean interval and smoking hazard explained -> results/surrogate_followup.json
    python3 pte_paired.py           # paired comparison of hazard explained -> results/pte_paired.json
    python3 valuation_core.py       # life-table core -> results/valuation_core.json
    python3 hhs_valuation_repro.py  # HHS VSL and worked-example reproduction -> results/hhs_valuation_repro.json

Requires Python 3, numpy, pandas, openpyxl, scikit-learn. Baselines and every comparison were declared before each run.

The two small public-domain source files the valuation cores read are included in sources/ (NVSR 70(19) Table 1 saved
as lt2019_t1.xlsx, and the text of ASPE Appendix D), so those two scripts run as cloned. The results/ folder holds the
outputs of every script as run; rerun any script and compare.
