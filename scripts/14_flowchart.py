# Study design and participant flow (Supplementary Fig. 1): UK Biobank cohort, brain-imaging subsample,
# NHANES external sample and the genetic (summary-statistics) analyses. Counts are recomputed from the data files
# with the same definitions as prep.R, 07_nhanes_analysis.R and 09_qsm_nigra.py.
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from fig_tiff import save_tiff
import pandas as pd
from matplotlib.patches import FancyBboxPatch

P = "/path/to/project/"

# ------------------------------------------------------------------ UK Biobank counts
co = pd.read_csv(P + "data/cohort.csv", parse_dates=["q_date"])
oc = pd.read_csv(P + "data/outcomes.csv", usecols=["eid", "G20", "F00", "F01", "F02", "F03", "G30"])
n_all = len(co)
miss_bmf = co.bmf_day.isna()
s1 = co[~miss_bmf]
EX = ["excl_ibd", "excl_coeliac", "excl_crc", "excl_colectomy"]
ex_n = {c: int(s1[c].sum()) for c in EX}
any_ex = s1[EX].sum(axis=1) > 0
coh = s1[~any_ex & (s1.fu_end_days > 0)].merge(oc, on="eid")
prev = lambda cols: ((coh[cols] <= 0) | (coh[cols] == -99999)).any(axis=1)
pd_risk, dem_risk = int((~prev(["G20"])).sum()), int((~prev(["F00", "F01", "F02", "F03", "G30"])).sum())
k = pd.read_csv(P + "results/key_outcomes.csv")
ev = lambda o: int(k[(k.outcome == o) & k.analysis.str.startswith("Model 2")].events.iat[0])
n_outcomes = pd.read_csv(P + "results/phewas_model2.csv").code.nunique()
fu_d = np.median(coh.fu_end_days) / 365.25
fu_m = np.median(coh.fu_end_death_days) / 365.25

# imaging subsample (closest imaging visit to the questionnaire, as in 09_qsm_nigra.py)
img = pd.read_pickle(P + "data/imaging_qsm_raw.pkl"); img["eid"] = img.eid.astype("int64")
b = coh.merge(img, on="eid", how="left")
num = lambda c: pd.to_numeric(b.get(c), errors="coerce")
gaps = {i: (pd.to_datetime(b[f"p53_i{i}"], errors="coerce") - b.q_date).dt.days.abs() for i in (2, 3)}
use3 = gaps[3].notna() & (gaps[2].isna() | (gaps[3] < gaps[2]))
pick = lambda f: pd.Series(np.where(use3, num(f"p{f}_i3"), num(f"p{f}_i2")))
sn = pd.concat([pick(24481), pick(24482)], axis=1).mean(axis=1)
gap = pd.Series(np.where(use3, gaps[3], gaps[2]))
img_date = pd.Series(np.where(use3, pd.to_datetime(b.p53_i3, errors="coerce"), pd.to_datetime(b.p53_i2, errors="coerce")))
has = sn.notna()
within = has & (gap <= 3 * 365.25)
pd_img = pd.to_numeric(b.G20, errors="coerce") <= (img_date - b.q_date).dt.days
img_final = within & ~pd_img
q = pd.read_csv(P + "results/qsm_nigra.csv")
n_qsm = int(q[(q.outcome == "qsm_Substantia nigra") & (q.model == "Model 2")].n.iat[0])

# ------------------------------------------------------------------ NHANES counts
nh = pd.read_csv(P + "data/nhanes/nhanes_bowel.csv")
nh_bmf = nh[nh.bmf_cat.notna()]
nh_final = nh_bmf[(nh_bmf.eligible_mort == 1) & nh_bmf.fu_months.notna()]
nh_pd = nh[(nh.age >= 40) & nh.stool_form.notna()]

# ------------------------------------------------------------------ drawing
# Boxes are placed by their top edge; height follows the number of text lines, and every text line is checked
# against its box after rendering so nothing spills over.
plt.rcParams.update({"font.family": "Arial", "font.size": 7})
fig, ax = plt.subplots(figsize=(7.2, 8.6))
ax.set_xlim(0, 100); ax.set_ylim(12, 99); ax.axis("off")
MAIN, EXCL, SUB, EDGE = "#eaf2fb", "#f4f4f2", "#fdf3e7", "#4a4a4a"
LH, PAD = 1.55, 1.3
BOXES, TEXTS = [], []


def box(x, top, w, lines, fc=MAIN, bold_first=False, center=False, fs=6.4, bold_all=False):
    h = len(lines) * LH + 2 * PAD
    ax.add_patch(FancyBboxPatch((x, top - h), w, h, boxstyle="round,pad=0,rounding_size=0.9", fc=fc, ec=EDGE, lw=0.7))
    for i, ln in enumerate(lines):
        t = ax.text(x + w / 2 if center else x + 1.3, top - PAD - (i + 0.5) * LH, ln, ha="center" if center else "left",
                    va="center", fontsize=fs, fontweight="bold" if (bold_all or (bold_first and i == 0)) else "normal")
        TEXTS.append((t, x, w))
    BOXES.append((x, top - h, w, h))
    return top - h


def arrow(x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", lw=0.7, color=EDGE, mutation_scale=8,
                                                                    shrinkA=0, shrinkB=0))


def side(x0, x1, y):
    ax.plot([x0, x1], [y, y], lw=0.7, color=EDGE)


def label(x, y, text):
    ax.text(x, y, text, fontsize=8, fontweight="bold", ha="left", va="bottom", color="#1f1f1f")


# --- a: UK Biobank (main flow in column 1, exclusions and the imaging branch in column 2)
label(0, 96.6, "a  UK Biobank")
MX, MW, CX = 0, 33, 16.5
EXX, EXW = 36, 31.5
b1 = box(MX, 95.5, MW, ["Completed the Digestive Health", "Questionnaire (Mar 2017 – Jun 2020)", f"n = {n_all:,}"], center=True)
e1_top = b1 - 1.5
e1 = box(EXX, e1_top, EXW, ["Excluded: missing bowel movement", f"frequency, n = {int(miss_bmf.sum()):,}"], fc=EXCL)
side(CX, EXX, (e1_top + e1) / 2)
b2_top = e1 - 2.2
arrow(CX, b1, CX, b2_top)
b2 = box(MX, b2_top, MW, [f"n = {n_all - int(miss_bmf.sum()):,}"], center=True)
e2_top = b2 - 1.5
e2 = box(EXX, e2_top, EXW, ["Excluded: diagnosed before the", f"questionnaire, n = {int(any_ex.sum()):,}*",
                             f"  IBD  {ex_n['excl_ibd']:,}", f"  Coeliac/malabsorption  {ex_n['excl_coeliac']:,}",
                             f"  Colorectal cancer  {ex_n['excl_crc']:,}", f"  Colectomy/rectal excision  {ex_n['excl_colectomy']:,}"],
         fc=EXCL)
side(CX, EXX, (e2_top + e2) / 2)
b3_top = e2 - 2.2
arrow(CX, b2, CX, b3_top)
b3 = box(MX, b3_top, MW, ["Analytic cohort", f"n = {len(coh):,}"], center=True, bold_all=True)
br_top = b3 - 4.0
arrow(CX, b3, CX, br_top)
ax.annotate("", xy=(EXX + EXW / 2, br_top), xytext=(MX + MW, (b3_top + b3) / 2),
            arrowprops=dict(arrowstyle="-|>", lw=0.7, color=EDGE, mutation_scale=8, shrinkA=0, shrinkB=0,
                            connectionstyle="angle,angleA=0,angleB=90,rad=0"))
pdc, dmc = ev("Parkinson's disease"), ev("Dementia")
box(MX, br_top, MW, ["Prospective analyses", f"Phenome-wide scan: {n_outcomes} outcomes", "  (≥100 incident events each)", "",
                     "Parkinson's disease", f"  at risk {pd_risk:,}; incident cases {pdc}", "Dementia",
                     f"  at risk {dem_risk:,}; incident cases {dmc}", "", "Median follow-up",
                     f"  {fu_d:.1f} y (disease); {fu_m:.1f} y (death)"], fc=SUB, bold_first=True)
i1 = box(EXX, br_top, EXW, ["Brain imaging subsample", "Nigral QSM at the imaging visit", "closest to the questionnaire",
                            f"n = {int(has.sum()):,}"], fc=SUB, bold_first=True)
i2_top = i1 - 2.6
arrow(EXX + EXW / 2, i1, EXX + EXW / 2, i2_top)
i2 = box(EXX, i2_top, EXW, ["Excluded", f"  >3 y from questionnaire  {int((has & ~(gap <= 3 * 365.25)).sum()):,}",
                            f"  PD before imaging  {int((within & pd_img).sum()):,}"], fc=EXCL)
i3_top = i2 - 2.6
arrow(EXX + EXW / 2, i2, EXX + EXW / 2, i3_top)
box(EXX, i3_top, EXW, ["Imaging analysis", f"n = {int(img_final.sum()):,} (median interval {np.median(gap[img_final]) / 365.25:.1f} y)",
                       f"{n_qsm:,} after removing extreme values"], fc=SUB, bold_first=True)

# --- b: NHANES (column 3)
RX, RW = 70.5, 29.5
label(RX, 96.6, "b  NHANES 2005–2010")
n1 = box(RX, 95.5, RW, ["Adults ≥20 y who completed the", "bowel health questionnaire", f"n = {len(nh):,}"], center=True)
n2_top = n1 - 2.6
arrow(RX + RW / 2, n1, RX + RW / 2, n2_top)
n2 = box(RX, n2_top, RW, ["Excluded", f"  missing frequency data  {len(nh) - len(nh_bmf):,}",
                          f"  no mortality linkage  {len(nh_bmf) - len(nh_final):,}"], fc=EXCL)
n3_top = n2 - 2.6
arrow(RX + RW / 2, n2, RX + RW / 2, n3_top)
n3 = box(RX, n3_top, RW, ["Mortality analysis", f"n = {len(nh_final):,}; {int(nh_final.death.sum()):,} deaths",
                          f"median follow-up {np.median(nh_final.fu_months) / 12:.1f} y"], fc=SUB, bold_first=True)
n4_top = n3 - 2.6
n4 = box(RX, n4_top, RW, ["Antiparkinsonian drug use", "cross-sectional, age ≥40 y", f"n = {len(nh_pd):,}; {int(nh_pd.pd_any_rx.sum())} users of any drug"],
         fc=SUB, bold_first=True)
# elbow from the survey sample to the drug-use subsample (drawn along the right margin)
yr0, yr1 = (95.5 + n1) / 2, (n4_top + n4) / 2
ax.plot([RX + RW, RX + RW + 1.6, RX + RW + 1.6], [yr0, yr0, yr1], lw=0.7, color=EDGE, clip_on=False)
ax.annotate("", xy=(RX + RW, yr1), xytext=(RX + RW + 1.6, yr1), arrowprops=dict(arrowstyle="-|>", lw=0.7, color=EDGE,
            mutation_scale=8, shrinkA=0, shrinkB=0), annotation_clip=False)

# --- c: genetic analyses (column 3)
c_lab = n4 - 4.0
label(RX, c_lab, "c  Genetic analyses")
box(RX, c_lab - 1.1, RW, ["Summary-statistics MR and coloc", "Stool frequency GWAS (n = 167,966)", "  12 lead variants (F = 29–96)", "",
                          "Forward MR → FinnGen R12", "  PD 5,861 cases / 494,487 controls", "",
                          "Reverse MR: PD GWAS (Nalls 2019)", "  20 variants; FinnGen R12 PD, 3", "",
                          "Colocalisation: ±500 kb, 12 loci"], fc=SUB, bold_first=True)
foot_y = min(y for _, y, _, _ in BOXES) - 3.2
ax.set_ylim(foot_y - 0.5, 99)
ax.text(0, foot_y, "*Exclusion categories overlap. PD, Parkinson's disease; QSM, quantitative susceptibility mapping; IBD, inflammatory "
                 "bowel disease;\nMR, Mendelian randomization; GWAS, genome-wide association study; NHANES, US National Health and "
                 "Nutrition Examination Survey.", fontsize=6, va="bottom", color="#555555")

# overflow check
fig.canvas.draw()
inv = ax.transData.inverted()
for t, x, w in TEXTS:
    bb = t.get_window_extent()
    (x0, _), (x1, _) = inv.transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
    if x0 < x + 0.3 or x1 > x + w - 0.3:
        print(f"OVERFLOW: {t.get_text()!r} spans {x0:.1f}-{x1:.1f}, box {x:.1f}-{x + w:.1f}")

for ext in ("pdf", "png"):
    fig.savefig(P + f"figures/FigS1_flowchart.{ext}", dpi=300, bbox_inches="tight")
save_tiff(fig, P + "figures/FigS1_flowchart.tif")
print("UKB", n_all, int(miss_bmf.sum()), int(any_ex.sum()), ex_n, len(coh), pd_risk, dem_risk)
print("imaging", int(has.sum()), int((has & ~(gap <= 3 * 365.25)).sum()), int((within & pd_img).sum()), int(img_final.sum()), n_qsm)
print("NHANES", len(nh), len(nh) - len(nh_bmf), len(nh_bmf) - len(nh_final), len(nh_final), len(nh_pd))
