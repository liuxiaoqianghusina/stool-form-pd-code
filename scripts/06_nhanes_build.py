# External validation dataset: NHANES 2005-2010 bowel health questionnaire (BHQ) + NCHS public-use linked mortality (to 2019).
# Exposure categories and covariates are harmonised with the UK Biobank analysis where the data allow.
import glob
import re

import numpy as np
import pandas as pd

RAW = "/path/to/project/data/nhanes/raw/"
OUT = "/path/to/project/data/nhanes/nhanes_bowel.csv"
CYCLES = {"D": 1, "E": 2, "F": 3}


def read(name, cols=None):
    parts = []
    for s in CYCLES:
        d = pd.read_sas(f"{RAW}{name}_{s}.xpt", format="xport")
        if cols:
            d = d[[c for c in cols if c in d.columns]]
        d["cycle"] = s
        parts.append(d)
    return pd.concat(parts, ignore_index=True)


demo = read("DEMO", ["SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH1", "DMDEDUC2", "INDFMPIR", "WTMEC2YR", "SDMVPSU", "SDMVSTRA"])
bhq = read("BHQ", ["SEQN", "BHD050", "BHQ060"]).drop(columns="cycle")
smq = read("SMQ", ["SEQN", "SMQ020", "SMQ040"]).drop(columns="cycle")
alq = read("ALQ", ["SEQN", "ALQ101", "ALQ130"]).drop(columns="cycle")
bmx = read("BMX", ["SEQN", "BMXBMI"]).drop(columns="cycle")
diq = read("DIQ", ["SEQN", "DIQ010"]).drop(columns="cycle")
ghb = read("GHB", ["SEQN", "LBXGH"]).drop(columns="cycle")
bpq = read("BPQ", ["SEQN", "BPQ020"]).drop(columns="cycle")
mcq = read("MCQ", ["SEQN", "MCQ160B", "MCQ160C", "MCQ160E", "MCQ160F", "MCQ220"]).drop(columns="cycle")
dpq = read("DPQ").drop(columns="cycle")
dr1 = read("DR1TOT", ["SEQN", "DR1TFIBE", "DR1TKCAL", "DR1_320Z", "DR1TCAFF"]).drop(columns="cycle")
paq = read("PAQ")
rx = read("RXQ_RX", ["SEQN", "RXDDRUG"]).drop(columns="cycle")

d = demo.merge(bhq, on="SEQN", how="inner")
for t in [smq, alq, bmx, diq, ghb, bpq, mcq, dpq, dr1, rx.groupby("SEQN").RXDDRUG.apply(lambda s: "|".join(
        x.decode() if isinstance(x, bytes) else str(x) for x in s)).reset_index()]:
    d = d.merge(t, on="SEQN", how="left")
d = d.merge(paq[[c for c in paq.columns if c in ("SEQN", "PAD200", "PAD320", "PAQ650", "PAQ665")]], on="SEQN", how="left")

# ------------------------------------------------------------------ linked mortality (fixed-width public-use file)
mort = []
for f in sorted(glob.glob(RAW + "NHANES_*_MORT_2019_PUBLIC.dat")):
    m = pd.read_fwf(f, colspecs=[(0, 6), (14, 15), (15, 16), (16, 19), (19, 20), (20, 21), (42, 45), (45, 48)],
                    names=["SEQN", "ELIGSTAT", "MORTSTAT", "UCOD_LEADING", "DIABETES_MC", "HYPERTEN_MC", "PERMTH_INT", "PERMTH_EXM"],
                    na_values=["."])
    mort.append(m)
mort = pd.concat(mort, ignore_index=True)
d = d.merge(mort, on="SEQN", how="left")

# ------------------------------------------------------------------ analytic variables
out = pd.DataFrame({"SEQN": d.SEQN, "cycle": d.cycle})
out["age"] = d.RIDAGEYR
out["female"] = (d.RIAGENDR == 2).astype(int)
out["race"] = d.RIDRETH1.map({1: "Mexican American", 2: "Other Hispanic", 3: "Non-Hispanic White", 4: "Non-Hispanic Black", 5: "Other"})
out["education"] = d.DMDEDUC2.map({1: "<High school", 2: "<High school", 3: "High school", 4: "Some college", 5: "College+"})
out["pir"] = d.INDFMPIR
out["weight"] = d.WTMEC2YR / 3  # three combined 2-year cycles
out["psu"], out["strata"] = d.SDMVPSU, d.SDMVSTRA
bm = d.BHD050.where(d.BHD050 < 777)
out["bm_week"] = bm
out["bmf_cat"] = pd.cut(bm, [0, 6.5, 13.5, 20.5, 1e9], labels=["<1/day", "1/day", "2/day", ">=3/day"]).astype(str).where(bm.notna())
out["bmf_low"] = np.where(bm <= 2, "<=2/week", np.where(bm < 7, "3-6/week", None))
bristol = d.BHQ060.where(d.BHQ060.between(1, 7))
out["bristol"] = bristol
out["stool_form"] = np.select([bristol <= 2, bristol >= 6], ["Hard (BSFS 1-2)", "Loose (BSFS 6-7)"], default="Normal (BSFS 3-5)")
out.loc[bristol.isna(), "stool_form"] = None
out["smoking"] = np.select([d.SMQ020 == 2, (d.SMQ020 == 1) & (d.SMQ040 == 3), (d.SMQ020 == 1) & d.SMQ040.isin([1, 2])],
                           ["Never", "Former", "Current"], default="Missing")
out["drinks_day"] = d.ALQ130.where(d.ALQ130 < 777).fillna(0).where(d.ALQ101.notna() | d.ALQ130.notna())
out["bmi"] = d.BMXBMI
out["active"] = ((d.get("PAD200") == 1) | (d.get("PAD320") == 1) | (d.get("PAQ650") == 1) | (d.get("PAQ665") == 1)).astype(int)
out["fibre_g"] = d.DR1TFIBE
out["energy_kcal"] = d.DR1TKCAL
out["water_g"] = d.DR1_320Z
out["caffeine_mg"] = d.DR1TCAFF
out["diabetes"] = ((d.DIQ010 == 1) | (d.LBXGH >= 6.5)).astype(int)
out["hypertension"] = (d.BPQ020 == 1).astype(int)
out["cvd_history"] = d[["MCQ160B", "MCQ160C", "MCQ160E", "MCQ160F"]].eq(1).any(axis=1).astype(int)
out["cancer_history"] = (d.MCQ220 == 1).astype(int)
phq = d[[f"DPQ0{i}0" for i in range(1, 10)]].where(lambda x: x <= 3)
out["phq9"] = phq.sum(axis=1, min_count=9)
out["depression"] = (out.phq9 >= 10).astype(float).where(out.phq9.notna())
drugs = d.RXDDRUG.fillna("").str.upper()
flag = lambda rx: drugs.str.contains(rx, regex=True).astype(int)
out["metformin"] = flag(r"METFORMIN")
out["opioid"] = flag(r"CODEINE|HYDROCODONE|OXYCODONE|MORPHINE|TRAMADOL|FENTANYL|METHADONE|HYDROMORPHONE|BUPRENORPHINE|PROPOXYPHENE")
out["laxative"] = flag(r"DOCUSATE|SENNA|POLYETHYLENE GLYCOL|LACTULOSE|BISACODYL|PSYLLIUM|METHYLCELLULOSE|LAXATIVE|MAGNESIUM HYDROXIDE")
out["antidepressant"] = flag(r"SERTRALINE|CITALOPRAM|ESCITALOPRAM|FLUOXETINE|PAROXETINE|VENLAFAXINE|DULOXETINE|BUPROPION|TRAZODONE|MIRTAZAPINE|AMITRIPTYLINE|NORTRIPTYLINE")
out["ppi"] = flag(r"OMEPRAZOLE|ESOMEPRAZOLE|LANSOPRAZOLE|PANTOPRAZOLE|RABEPRAZOLE")
# treated Parkinson's disease: PD-specific drugs (primary); + dopamine agonists, also used for restless legs (sensitivity)
out["pd_specific_rx"] = flag(r"LEVODOPA|CARBIDOPA|RASAGILINE|SELEGILINE|ENTACAPONE|TOLCAPONE")
out["pd_any_rx"] = ((out.pd_specific_rx == 1) | (flag(r"PRAMIPEXOLE|ROPINIROLE|ROTIGOTINE") == 1)).astype(int)
out["eligible_mort"] = (d.ELIGSTAT == 1).astype(int)
out["death"] = d.MORTSTAT
out["fu_months"] = d.PERMTH_EXM
ucod = d.UCOD_LEADING
out["death_cvd"] = ((d.MORTSTAT == 1) & ucod.isin([1, 5])).astype(int)
out["death_cancer"] = ((d.MORTSTAT == 1) & (ucod == 2)).astype(int)
out["death_other_noncvd_noncancer"] = ((d.MORTSTAT == 1) & ~ucod.isin([1, 2, 5])).astype(int)
out.to_csv(OUT, index=False)
a = out[(out.age >= 20) & out.bmf_cat.notna() & (out.eligible_mort == 1)]
print("BHQ respondents:", len(out), "| adults with BMF and mortality linkage:", len(a))
print("deaths:", int(a.death.sum()), "| CVD deaths:", int(a.death_cvd.sum()), "| cancer deaths:", int(a.death_cancer.sum()),
      "| median FU (y):", round(a.fu_months.median() / 12, 1))
print(a.bmf_cat.value_counts().to_dict()); print(a.stool_form.value_counts().to_dict())
