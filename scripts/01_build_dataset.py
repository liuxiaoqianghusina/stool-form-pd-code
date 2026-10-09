# Build the analytic dataset for: bowel habits (UKB digestive health questionnaire, 2017-2020) and multisystem health.
# Outputs (data/):
#   cohort.csv    one row per respondent: exposures, covariates, exclusion flags, follow-up end (days from questionnaire)
#   outcomes.csv  eid + one column per outcome = days from questionnaire to first diagnosis
#                 (negative/zero = prevalent, -99999 = prevalent with unknown date, empty = never)
#   outcome_dictionary.csv  outcome code, label, chapter, source
import csv
import re

import numpy as np
import pandas as pd

D = "/path/to/ukb_data/"
OUT = "/path/to/project/data/"
SPECIAL = {"1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07"}

rows = list(csv.DictReader(open(D + "app45551_20251118060954.dataset.data_dictionary.csv", encoding="utf-8")))
codings = {}
for r in csv.reader(open(D + "app45551_20251118060954.dataset.codings.csv", encoding="utf-8")):
    if r:
        codings.setdefault(r[0], {})[r[1]] = r[2]

# ------------------------------------------------------------------ exposures (digestive health questionnaire)
dig = [f"p{i}" for i in list(range(21023, 21045)) + [21064, 21065, 21066, 21067, 21068]]
of = pd.read_csv(D + "UKB619657_online_followup_2025Nov.csv", usecols=["eid"] + dig, dtype=str)
of["q_date"] = pd.to_datetime(of.p21023, errors="coerce")
of = of[of.q_date.notna()].copy()
num = lambda s: pd.to_numeric(s, errors="coerce")
of["bmf_day"] = num(of.p21044).where(lambda s: s >= 0)
of["bmf_least_week"] = num(of.p21043).where(lambda s: s >= 0)
of["bmf_most_day"] = num(of.p21042).where(lambda s: s >= 0)
FREQ = {"-500": 0, "-501": 1, "-502": 2, "-503": 3, "-504": 4}  # never, sometimes, often, most of the time, always
of["hard_stool"] = of.p21033.map(FREQ)
of["loose_stool"] = of.p21034.map(FREQ)
of["ibs_ever_dx_q"] = (num(of.p21024) == 1).astype(float).where(num(of.p21024).isin([0, 1]))
of["sensitive_stomach"] = (num(of.p21064) == 1).astype(float).where(num(of.p21064).isin([0, 1]))
of["bowel_satisfaction"] = num(of.p21040).where(lambda s: s >= 0)
of = of[["eid", "q_date", "bmf_day", "bmf_least_week", "bmf_most_day", "hard_stool", "loose_stool",
         "ibs_ever_dx_q", "sensitive_stomach", "bowel_satisfaction"]]

# ------------------------------------------------------------------ baseline covariates (instance 0)
header = open(D + "UKB619657_assessment_center_2025Nov.csv").readline().replace('"', "").strip().split(",")
want = [21003, 53, 54, 21000, 6138, 20116, 1558, 22040, 21001, 1309, 1289, 1299, 1438, 1448, 1458, 1468, 1528, 1488,
        1498, 137, 6154, 6150, 2443, 20127]
ac_cols = ["eid"] + [c for c in header if re.match(r"p(\d+)_i0(_a\d+)?$", c)
                     and int(re.match(r"p(\d+)", c).group(1)) in want]
ac_cols += [c for c in header if re.match(r"p20003_i0_a\d+$", c)]
ac = pd.read_csv(D + "UKB619657_assessment_center_2025Nov.csv", usecols=ac_cols, dtype=str)
pc = pd.read_csv(D + "UKB619657_population_characteristics_2025Nov.csv", usecols=["eid", "p31", "p34", "p52", "p22189"], dtype=str)
cv = ac.merge(pc, on="eid", how="left")
n = lambda c: pd.to_numeric(cv[c], errors="coerce")
cov = pd.DataFrame({"eid": cv.eid})
cov["sex"] = n("p31")
cov["birth_year"] = n("p34")
cov["birth_month"] = n("p52")
cov["centre"] = n("p54_i0")
cov["baseline_date"] = pd.to_datetime(cv.p53_i0, errors="coerce")
eth = n("p21000_i0")
cov["ethnicity"] = np.select([eth.isin([1, 1001, 1002, 1003]), eth.isin([3, 3001, 3002, 3003, 3004, 5]),
                              eth.isin([4, 4001, 4002, 4003]), eth.isin([2, 2001, 2002, 2003, 2004, 6])],
                             ["White", "Asian", "Black", "Mixed/other"], default="Unknown")
cov["tdi"] = n("p22189")
q = cv[[c for c in cv.columns if c.startswith("p6138_i0")]].astype(str).agg("|".join, axis=1)
cov["degree"] = q.str.contains(r"(^|\|)1(\||$)").astype(float).where(q.str.contains(r"(^|\|)(-7|[1-6])(\||$)"))
cov["smoking"] = n("p20116_i0").where(lambda s: s >= 0)
cov["alcohol_freq"] = n("p1558_i0").where(lambda s: s > 0)
cov["met_min"] = n("p22040_i0")
cov["bmi"] = n("p21001_i0")
def portions(c):  # "-10 less than one" -> 0.5 ; negatives (do not know / prefer not) -> NaN
    v = n(c)
    return v.where(v >= 0, np.where(v == -10, 0.5, np.nan))
cov["fruit"] = portions("p1309_i0")
cov["veg"] = portions("p1289_i0").fillna(0) + portions("p1299_i0").fillna(0)
cov["veg"] = cov["veg"].where(cv.p1289_i0.notna() | cv.p1299_i0.notna())
cov["bread_slices"] = portions("p1438_i0")
cov["wholemeal_bread"] = (n("p1448_i0") == 3).astype(float).where(n("p1448_i0") > 0)
cov["cereal_bowls"] = portions("p1458_i0")
cov["fibre_cereal"] = n("p1468_i0").isin([1, 2, 3, 4]).astype(float).where(n("p1468_i0") > 0)
cov["water_glasses"] = portions("p1528_i0")
cov["tea_cups"] = portions("p1488_i0")
cov["coffee_cups"] = portions("p1498_i0")
cov["n_medications"] = n("p137_i0")
cov["neuroticism"] = n("p20127_i0")
otc = cv[[c for c in cv.columns if c.startswith("p6154_i0")]].astype(str).agg("|".join, axis=1)
lax_code = [k for k, v in codings[[r for r in rows if r["name"] == "p6154_i0"][0]["coding_name"]].items() if v.startswith("Laxatives")][0]
ome_code = [k for k, v in codings[[r for r in rows if r["name"] == "p6154_i0"][0]["coding_name"]].items() if v.startswith("Omeprazole")][0]
cov["laxative_otc"] = otc.str.contains(rf"(^|\|){lax_code}(\||$)").astype(float)
cov["omeprazole_otc"] = otc.str.contains(rf"(^|\|){ome_code}(\||$)").astype(float)
cod4 = codings["data_coding_4"]
meds = cv[[c for c in cv.columns if c.startswith("p20003_i0")]].apply(pd.to_numeric, errors="coerce")
def medflag(rx):
    codes = {int(k) for k, v in cod4.items() if re.search(rx, v, re.I)}
    return meds.isin(codes).any(axis=1).astype(float)
cov["laxative_rx"] = medflag(r"lactulose|senna|bisacodyl|macrogol|movicol|ispaghula|fybogel|docusate|laxative|methylcellulose|sterculia|picosulfate")
cov["laxative_any"] = ((cov.laxative_rx == 1) | (cov.laxative_otc == 1)).astype(float)
cov["metformin"] = medflag(r"metformin")
cov["ppi"] = ((medflag(r"omeprazole|lansoprazole|pantoprazole|rabeprazole|esomeprazole") == 1) | (cov.omeprazole_otc == 1)).astype(float)
cov["opioid"] = medflag(r"codeine|co-codamol|co-dydramol|dihydrocodeine|tramadol|(?<!apo)morphine|oxycodone|fentanyl|buprenorphine|pethidine|tapentadol|meptazinol|diamorphine")
cov["antidepressant"] = medflag(r"amitriptyline|citalopram|escitalopram|fluoxetine|sertraline|paroxetine|venlafaxine|duloxetine|mirtazapine|dosulepin|trazodone|nortriptyline|imipramine|clomipramine|lofepramine")
cov["antibiotic"] = medflag(r"amoxicillin|doxycycline|clarithromycin|erythromycin|ciprofloxacin|trimethoprim|nitrofurantoin|flucloxacillin|cefalexin|metronidazole|penicillin")
cov["iron_supp"] = medflag(r"ferrous|iron")
cov["calcium_channel_blocker"] = medflag(r"amlodipine|nifedipine|felodipine|diltiazem|verapamil|lercanidipine")

# ------------------------------------------------------------------ outcomes: first occurrences (3-character ICD-10)
ho_header = open(D + "UKB619657_health_outcomes_2025Nov.csv").readline().replace('"', "").strip().split(",")
fo = {}
label = {}
for r in rows:
    m = re.match(r"Date ([A-Z]\d\d) first reported \((.*)\)", r["title"])
    if m and "First occurrences" in r["folder_path"] and r["name"] in ho_header:
        fo[m.group(1)] = r["name"]
        label[m.group(1)] = m.group(2)
fo = {k: v for k, v in fo.items() if k[0] not in "OPQRSTUVWXYZ"}  # drop pregnancy, perinatal, congenital, symptoms, injuries
cancer_date_cols = [c for c in ho_header if c.startswith("p40005_i")]
cancer_type_cols = [c for c in ho_header if c.startswith("p40006_i")]
op_cols = ["p41272"] + [c for c in ho_header if c.startswith("p41282_a")]
ho = pd.read_csv(D + "UKB619657_health_outcomes_2025Nov.csv",
                 usecols=["eid", "p40000_i0", "p40001_i0"] + list(fo.values()) + cancer_date_cols + cancer_type_cols + op_cols,
                 dtype=str)
ho = ho[ho.eid.isin(of.eid)].reset_index(drop=True)

# administrative end of follow-up by UK nation of assessment centre. Hospital-data cut-offs were read off the monthly
# counts of first-occurrence events (sharp drop afterwards): England 2023-03-31, Scotland 2022-09-30, Wales 2022-05-31.
# The death register is complete to 2024-09-30.
SCOT, WALES = {11004, 11005}, {11003, 11022, 11023}
cen = cov.set_index("eid").centre.reindex(ho.eid).values
nation = np.where(np.isin(cen, list(SCOT)), "Scotland", np.where(np.isin(cen, list(WALES)), "Wales", "England"))
alld = ho[list(fo.values())].apply(pd.to_datetime, errors="coerce")
end_by_nation = {"England": pd.Timestamp("2023-03-31"), "Scotland": pd.Timestamp("2022-09-30"), "Wales": pd.Timestamp("2022-05-31")}
DEATH_END = pd.Timestamp("2024-09-30")
print("Administrative censoring dates:", {k: str(v.date()) for k, v in end_by_nation.items()})

qd = of.set_index("eid").q_date.reindex(ho.eid).values
qd = pd.to_datetime(pd.Series(qd))
outcomes = pd.DataFrame({"eid": ho.eid})
odict = []
for code, col in sorted(fo.items()):
    raw = ho[col]
    days = (alld[col] - qd).dt.days.astype("float")
    days[raw.isin(SPECIAL)] = -99999
    outcomes[code] = days
    odict.append(dict(code=code, label=label[code], chapter=code[0], source="first_occurrence"))

# cancer registry (3-character C codes), earliest date per code
canc = {}
for dc, tc in zip(sorted(cancer_date_cols, key=lambda c: int(c.split("_i")[1])), sorted(cancer_type_cols, key=lambda c: int(c.split("_i")[1]))):
    t = ho[tc].fillna("").str[:3]
    d = pd.to_datetime(ho[dc], errors="coerce")
    for code in t[t.str.match(r"^C\d\d$")].unique():
        sel = t == code
        cur = canc.get(code, pd.Series(pd.NaT, index=ho.index))
        cur = cur.where(~sel | (cur.notna() & (cur <= d)), d.where(sel))
        canc[code] = cur
grouped = {"C18-C20 colorectal": ["C18", "C19", "C20"]}
for code, d in sorted(canc.items()):
    outcomes["CR_" + code] = (d - qd).dt.days.astype("float")
    odict.append(dict(code="CR_" + code, label=f"Cancer {code}", chapter="C", source="cancer_registry"))
for name, codes in grouped.items():
    stack = pd.concat([canc[c] for c in codes if c in canc], axis=1).min(axis=1)
    outcomes["CR_" + name.split()[0]] = (stack - qd).dt.days.astype("float")
    odict.append(dict(code="CR_" + name.split()[0], label=name, chapter="C", source="cancer_registry"))

# deaths (all-cause and by underlying-cause chapter)
death = pd.to_datetime(ho.p40000_i0, errors="coerce")
cause = ho.p40001_i0.fillna("")
outcomes["DEATH_all"] = (death - qd).dt.days.astype("float")
odict.append(dict(code="DEATH_all", label="All-cause mortality", chapter="death", source="death_register"))
for nm, rx in {"CVD": r"^I", "cancer": r"^C", "respiratory": r"^J", "digestive": r"^K", "neuro_mental": r"^[FG]"}.items():
    outcomes["DEATH_" + nm] = ((death - qd).dt.days.astype("float")).where(cause.str.match(rx))
    odict.append(dict(code="DEATH_" + nm, label=f"{nm} mortality", chapter="death", source="death_register"))

# ------------------------------------------------------------------ cohort file
cohort = of.merge(cov, on="eid", how="left")
cohort = cohort.set_index("eid").reindex(ho.eid).reset_index()
cohort["age_q"] = (cohort.q_date - pd.to_datetime(dict(year=cohort.birth_year, month=cohort.birth_month, day=15))).dt.days / 365.25
cohort["nation"] = nation
end = pd.Series([end_by_nation[x] for x in nation])
cohort["fu_end_days"] = np.minimum((end - cohort.q_date).dt.days, (death - cohort.q_date).dt.days.fillna(1e9))
cohort["fu_end_death_days"] = np.minimum((DEATH_END - cohort.q_date).dt.days, (death - cohort.q_date).dt.days.fillna(1e9))
# exclusions before questionnaire: IBD (K50, K51), coeliac/malabsorption (K90), colorectal cancer, colectomy/rectal excision
def prevalent(code):
    v = outcomes[code]
    return (v <= 0) | (v == -99999)
cohort["excl_ibd"] = (prevalent("K50") | prevalent("K51")).astype(int).values
cohort["excl_coeliac"] = prevalent("K90").astype(int).values
crc = outcomes["CR_C18-C20"]
cohort["excl_crc"] = (crc <= 0).astype(int).values
ops = ho.p41272.fillna("").str.split("|")
op_dates = ho[[c for c in ho.columns if c.startswith("p41282_a")]]
colect = []
for i, codes in enumerate(ops):
    hit = [j for j, c in enumerate(codes) if re.match(r"H0[4-9]|H1[01]|H33", c)]
    ds = [pd.to_datetime(op_dates.iat[i, j], errors="coerce") for j in hit if j < op_dates.shape[1]]
    colect.append(any(pd.notna(x) and x <= cohort.q_date.iat[i] for x in ds))
cohort["excl_colectomy"] = np.array(colect, dtype=int)
# prevalent conditions at questionnaire (for adjustment)
for nm, codes in {"prev_diabetes": ["E10", "E11", "E12", "E13", "E14"], "prev_hypertension": ["I10"],
                  "prev_ihd": ["I20", "I21", "I22", "I23", "I24", "I25"], "prev_stroke": ["I60", "I61", "I63", "I64"],
                  "prev_depression": ["F32", "F33"], "prev_anxiety": ["F40", "F41"], "prev_ibs": ["K58"],
                  "prev_ckd": ["N18"], "prev_copd": ["J44"], "prev_hypothyroid": ["E03"]}.items():
    cohort[nm] = np.column_stack([prevalent(c).values for c in codes]).any(axis=1).astype(int)
cohort["prev_cancer"] = (pd.concat([outcomes[c] for c in outcomes.columns if c.startswith("CR_C")], axis=1) <= 0).any(axis=1).astype(int).values

cohort.drop(columns=["baseline_date"]).to_csv(OUT + "cohort.csv", index=False)
outcomes.to_csv(OUT + "outcomes.csv", index=False)
pd.DataFrame(odict).to_csv(OUT + "outcome_dictionary.csv", index=False)
print("cohort", cohort.shape, "| outcomes", outcomes.shape)
print("Exclusions:", {c: int(cohort[c].sum()) for c in cohort.columns if c.startswith("excl_")})
print("Follow-up (years) median:", round(cohort.fu_end_days.median() / 365.25, 2))
