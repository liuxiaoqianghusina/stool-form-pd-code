# Brain iron (quantitative susceptibility mapping, QSM) of the substantia nigra and bowel habits.
# Cross-sectional: imaging visit (instance 2 or 3) closest in time to the digestive health questionnaire, within 3 years.
# Outcomes: mean left/right median QSM and T2* of the substantia nigra (higher QSM / lower T2* = more iron), plus other
# subcortical nuclei as comparison regions. Linear models, outcome standardised (SD units).
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

D = "/path/to/ukb_data/"
P = "/path/to/project/"
REGIONS = {"Substantia nigra": (24481, 24482), "Putamen": (24471, 24472), "Caudate": (24469, 24470),
           "Pallidum": (24473, 24474), "Thalamus": (24467, 24468), "Hippocampus": (24475, 24476),
           "Amygdala": (24477, 24478), "Accumbens": (24479, 24480)}
T2S_SN = (24483, 24484)
fields = [f for pair in list(REGIONS.values()) + [T2S_SN] for f in pair] + [53, 54, 25000, 21003]
cols = ["eid"] + [f"p{f}_i{i}" for f in fields for i in (2, 3)]
import os
if os.path.exists(P + "data/imaging_qsm_raw.pkl"):  # cached extract of the 37 GB assessment-centre file
    img = pd.read_pickle(P + "data/imaging_qsm_raw.pkl")
else:
    header = open(D + "UKB619657_assessment_center_2025Nov.csv").readline().replace('"', "").strip().split(",")
    img = pd.read_csv(D + "UKB619657_assessment_center_2025Nov.csv", usecols=[c for c in cols if c in header], dtype=str)
    img.to_pickle(P + "data/imaging_qsm_raw.pkl")
img["eid"] = img.eid.astype("int64")

co = pd.read_csv(P + "data/cohort.csv", parse_dates=["q_date"])
oc = pd.read_csv(P + "data/outcomes.csv", usecols=["eid", "G20"])
df = co.merge(img, on="eid", how="inner").merge(oc, on="eid", how="left")
num = lambda c: pd.to_numeric(df.get(c), errors="coerce")
# choose the imaging instance closest to the questionnaire
gaps = {}
for i in (2, 3):
    gaps[i] = (pd.to_datetime(df[f"p53_i{i}"], errors="coerce") - df.q_date).dt.days.abs()
use3 = gaps[3].notna() & (gaps[2].isna() | (gaps[3] < gaps[2]))
inst = np.where(use3, 3, 2)
pick = lambda f: pd.Series(np.where(use3, num(f"p{f}_i3"), num(f"p{f}_i2")), index=df.index)
df["gap_days"] = np.where(use3, gaps[3], gaps[2])
df["img_date"] = np.where(use3, pd.to_datetime(df.p53_i3, errors="coerce"), pd.to_datetime(df.p53_i2, errors="coerce"))
df["age_img"] = pick(21003)
df["img_centre"] = pick(54).astype("Int64").astype(str)
df["head_size"] = pick(25000)
for name, (l, r) in REGIONS.items():
    df[f"qsm_{name}"] = pd.concat([pick(l), pick(r)], axis=1).mean(axis=1)
df["t2s_Substantia nigra"] = pd.concat([pick(T2S_SN[0]), pick(T2S_SN[1])], axis=1).mean(axis=1)

ok = (df.gap_days <= 3 * 365.25) & df["qsm_Substantia nigra"].notna() & df.bmf_day.notna()
ok &= ~(pd.to_numeric(df.G20, errors="coerce") <= (pd.to_datetime(df.img_date) - df.q_date).dt.days)  # no PD at imaging
ok &= (df.excl_ibd == 0) & (df.excl_coeliac == 0) & (df.excl_crc == 0) & (df.excl_colectomy == 0)
d = df[ok].copy()
d["bmf_cat"] = pd.Categorical(np.select([d.bmf_day == 0, d.bmf_day == 1, d.bmf_day == 2], ["<1/day", "1/day", "2/day"], ">=3/day"),
                              categories=["1/day", "<1/day", "2/day", ">=3/day"])
d["hard_cat"] = pd.Categorical(np.select([d.hard_stool == 0, d.hard_stool == 1], ["Never", "Sometimes"], "Often+"),
                               categories=["Never", "Sometimes", "Often+"])
d["loose_cat"] = pd.Categorical(np.select([d.loose_stool == 0, d.loose_stool == 1], ["Never", "Sometimes"], "Often+"),
                                categories=["Never", "Sometimes", "Often+"])
for c in ["tdi", "bmi", "met_min"]:
    d[c + "_m"] = d[c].fillna(d[c].median())
d["smoking_c"] = d.smoking.fillna(-1).astype(int).astype(str)
d["alcohol_c"] = d.alcohol_freq.fillna(-1).astype(int).astype(str)
COV1 = "age_img + I(age_img**2) + C(sex) + C(img_centre) + head_size"
COV2 = COV1 + " + tdi_m + bmi_m + np.log1p(met_min_m) + C(smoking_c) + C(alcohol_c) + laxative_any + metformin + opioid + prev_diabetes + prev_hypertension"
print(f"N with QSM within 3 y of questionnaire (no PD, exclusions applied): {len(d)}; median gap {d.gap_days.median()/365.25:.2f} y")
print(d.hard_cat.value_counts().to_dict(), d.loose_cat.value_counts().to_dict())
rows = []
outs = [f"qsm_{n}" for n in REGIONS] + ["t2s_Substantia nigra"]
for y in outs:
    t = d.dropna(subset=[y, "head_size", "age_img"]).copy()
    lo, hi = t[y].quantile([0.005, 0.995])
    t = t[(t[y] >= lo) & (t[y] <= hi)]  # trim extreme values (segmentation failures)
    t["yz"] = (t[y] - t[y].mean()) / t[y].std()
    for ex in ["hard_cat", "loose_cat", "bmf_cat"]:
        for mname, cov in [("Model 1", COV1), ("Model 2", COV2)]:
            m = smf.ols(f"yz ~ {ex} + {cov}", data=t).fit(cov_type="HC1")
            for k in m.params.index:
                if k.startswith(ex):
                    rows.append(dict(outcome=y, exposure=ex, level=k.split("T.")[1].rstrip("]"), model=mname, beta_sd=m.params[k],
                                     lo=m.conf_int().loc[k, 0], hi=m.conf_int().loc[k, 1], p=m.pvalues[k], n=int(m.nobs)))
        # linear trend across the five-level frequency (0-4)
        if ex != "bmf_cat":
            v = "hard_stool" if ex == "hard_cat" else "loose_stool"
            m = smf.ols(f"yz ~ {v} + {COV2}", data=t).fit(cov_type="HC1")
            rows.append(dict(outcome=y, exposure=v + " (per level, 0-4)", level="trend", model="Model 2", beta_sd=m.params[v],
                             lo=m.conf_int().loc[v, 0], hi=m.conf_int().loc[v, 1], p=m.pvalues[v], n=int(m.nobs)))
res = pd.DataFrame(rows)
res.to_csv(P + "results/qsm_nigra.csv", index=False)
pd.set_option("display.width", 220)
show = res[(res.model == "Model 2") & (res.level.isin(["Often+", "trend"]))].copy()
show["est"] = show.apply(lambda x: f"{x.beta_sd:+.3f} ({x.lo:+.3f},{x.hi:+.3f}) p={x.p:.2g}", axis=1)
print(show.pivot_table(index="outcome", columns="exposure", values="est", aggfunc="first").to_string())
