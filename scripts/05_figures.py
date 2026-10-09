# Draft main figures for the bowel-habits paper (static PDF + PNG).
# Colour encodes one job only, direction of association: red = higher risk, blue = lower risk, grey = not FDR-significant.
import textwrap

import matplotlib.pyplot as plt
import numpy as np
from fig_tiff import save_tiff
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, NullFormatter

P = "/path/to/project/"
RED, BLUE, GREY = "#e34948", "#2a78d6", "#c9c7c1"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e4df"
plt.rcParams.update({"font.size": 8, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "font.family": "DejaVu Sans"})
CHAPTERS = {"A": "Infections", "B": "Infections", "D": "Blood/immune", "E": "Endocrine/metabolic", "F": "Mental",
            "G": "Nervous", "H": "Eye/ear", "I": "Circulatory", "J": "Respiratory", "K": "Digestive", "L": "Skin",
            "M": "Musculoskeletal", "N": "Genitourinary", "C": "Cancer", "death": "Mortality"}
CH_ORDER = ["Infections", "Blood/immune", "Endocrine/metabolic", "Mental", "Nervous", "Eye/ear", "Circulatory",
            "Respiratory", "Skin", "Musculoskeletal", "Genitourinary", "Cancer", "Mortality"]
SHORT = {"G20": "Parkinson's disease", "F32": "Depression", "F41": "Anxiety", "E11": "Type 2 diabetes",
         "I10": "Hypertension", "I50": "Heart failure", "J45": "Asthma", "A09": "Infectious gastroenteritis",
         "N17": "Acute kidney injury", "F02": "Dementia (other diseases)", "D50": "Iron-deficiency anaemia",
         "E55": "Vitamin D deficiency", "DEATH_all": "All-cause mortality", "M47": "Spondylosis",
         "I95": "Hypotension", "F10": "Alcohol use disorder", "E66": "Obesity", "G47": "Sleep disorders",
         "A49": "Bacterial infection", "N32": "Bladder disorders", "DEATH_neuro_mental": "Neuro/mental mortality"}


def save(fig, name):
    fig.savefig(P + f"figures/{name}.pdf", bbox_inches="tight")
    fig.savefig(P + f"figures/{name}.png", dpi=300, bbox_inches="tight")
    save_tiff(fig, P + f"figures/{name}.tif")
    plt.close(fig)


# ------------------------------------------------------------------ Figure 2: phenome-wide signed Manhattan
r = pd.read_csv(P + "results/phewas_model2.csv")
r["chap"] = r.chapter.map(CHAPTERS)
r = r[r.chap.notna() & (r.chap != "Digestive")].copy()  # digestive outcomes go to the supplement
r["chap_i"] = r.chap.map({c: i for i, c in enumerate(CH_ORDER)})
order = r.drop_duplicates("code").sort_values(["chap_i", "code"]).reset_index(drop=True)
xpos = dict(zip(order.code, np.arange(len(order))))
r["x"] = r.code.map(xpos)
r["y"] = -np.log10(r.p) * np.sign(np.log(r.HR))
# stool form first (Parkinson's-disease-focused manuscript), then bowel movement frequency
PANELS = [("hard_cat", "Often+", "a  Hard/lumpy stools often or more (vs never)"),
          ("loose_cat", "Often+", "b  Loose/watery stools often or more (vs never)"),
          ("bmf_cat", ">=3/day", "c  Bowel movements ≥3 per day (vs 1/day)"),
          ("bmf_cat", "<1/day", "d  Bowel movements <1 per day (vs 1/day)")]
fig, axes = plt.subplots(4, 1, figsize=(7.2, 9.0), sharex=True)
for ax, (ex, lvl, title) in zip(axes, PANELS):
    d = r[(r.exposure == ex) & (r.level == lvl)]
    sig = d.fdr < 0.05
    thr = -np.log10(d.loc[sig, "p"].max()) if sig.any() else None
    ax.scatter(d.x[~sig], d.y[~sig], s=9, color=GREY, linewidths=0, zorder=2)
    ax.scatter(d.x[sig & (d.HR > 1)], d.y[sig & (d.HR > 1)], s=14, color=RED, edgecolor="white", linewidths=0.4, zorder=3)
    ax.scatter(d.x[sig & (d.HR < 1)], d.y[sig & (d.HR < 1)], s=14, color=BLUE, edgecolor="white", linewidths=0.4, zorder=3)
    ax.axhline(0, color=INK2, lw=0.6)
    if thr:
        for s in (1, -1):
            ax.axhline(s * thr, color=INK2, lw=0.6, ls=(0, (3, 3)))
    lab = pd.concat([d[sig & (d.HR < 1)], d[sig].sort_values("p").head(6)]).drop_duplicates("code")  # every inverse association labelled
    placed = []
    fig.canvas.draw()
    axbox = ax.get_window_extent()
    for _, row in lab.iterrows():
        name = SHORT.get(row.code, textwrap.shorten(str(row.label).capitalize(), 26, placeholder="…"))
        up = row.y > 0
        cands = [(4, 3 if up else -9, "left"), (-4, 3 if up else -9, "right"), (4, -9 if up else 3, "left"),
                 (-4, -9 if up else 3, "right")] + [(10, (3 + 8 * k) * (1 if up else -1), "left") for k in range(1, 4)]
        for dx, dy, ha in cands:  # first candidate position that overlaps no placed label and stays inside the axes
            t = ax.annotate(name, (row.x, row.y), xytext=(dx, dy), textcoords="offset points", fontsize=6.5,
                            color=INK, ha=ha,
                            arrowprops=None if abs(dy) <= 9 else dict(arrowstyle="-", color=INK2, lw=0.4))
            bb = t.get_window_extent()
            inside = axbox.x0 <= bb.x0 and bb.x1 <= axbox.x1 and axbox.y0 <= bb.y0 and bb.y1 <= axbox.y1
            if inside and not any(bb.overlaps(p) for p in placed):
                placed.append(bb)
                break
            t.remove()
    ax.set_title(title, loc="left", fontsize=8.5, color=INK, fontweight="bold")
    ax.set_ylabel("Signed −log$_{10}$ P", fontsize=7.5)
    ax.grid(axis="y", color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    for i in range(1, len(CH_ORDER)):
        b = order[order.chap_i == i].index.min()
        if not np.isnan(b):
            ax.axvline(b - 0.5, color=GRID, lw=0.5, zorder=1)
mids = order.groupby("chap_i").apply(lambda g: g.index.to_numpy().mean())
axes[-1].set_xticks(mids.values)
axes[-1].set_xticklabels([CH_ORDER[i] for i in mids.index], rotation=40, ha="right", fontsize=7)
axes[0].legend(handles=[Line2D([], [], marker="o", ls="", color=RED, label="Higher risk (FDR<0.05)"),
                        Line2D([], [], marker="o", ls="", color=BLUE, label="Lower risk (FDR<0.05)"),
                        Line2D([], [], marker="o", ls="", color=GREY, label="Not significant")],
               loc="upper right", fontsize=6.5, frameon=False, ncol=3, bbox_to_anchor=(1, 1.32))
fig.tight_layout(h_pad=1.2)
save(fig, "Fig2_phenome_wide")


# ------------------------------------------------------------------ forest helper
def forest(ax, rows, title, xlim=(0.2, 5)):
    y = np.arange(len(rows))[::-1]
    for yi, (lab, hr, lo, hi, p) in zip(y, rows):
        col = (RED if hr > 1 else BLUE) if p < 0.05 else GREY
        ax.plot([lo, hi], [yi, yi], color=col, lw=1.6, solid_capstyle="round")
        ax.plot(hr, yi, "o", color=col, ms=5, mec="white", mew=0.6)
        ax.text(xlim[1] * 1.06, yi, f"{hr:.2f} ({lo:.2f}–{hi:.2f})", va="center", fontsize=6.5, color=INK)
    ax.axvline(1, color=INK2, lw=0.7)
    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    ticks = [t for t in [0.25, 0.5, 1, 2, 4] if xlim[0] <= t <= xlim[1]]
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=7)
    ax.set_xlabel("Hazard ratio (95% CI)", fontsize=7.5)
    ax.set_title(title, loc="left", fontsize=8.5, color=INK, fontweight="bold")
    ax.grid(axis="x", color=GRID, lw=0.5)
    ax.set_axisbelow(True)


k = pd.read_csv(P + "results/key_outcomes.csv")
ANALYSES = [("Model 1: age + sex", "Model 1 (age, sex)"), ("Model 2: fully adjusted (main)", "Model 2 (main)"),
            ("Model 3: + coffee and tea", "+ coffee and tea"), ("Excluding events in first 2 years", "Excl. events <2 y"),
            ("Excluding events in first 4 years", "Excl. events <4 y"), ("Excluding IBS", "Excl. IBS"),
            ("Excluding laxative, metformin and opioid users", "Excl. laxative/metformin/opioid"),
            ("Follow-up window 0-2 y", "Follow-up 0–2 y"), ("Follow-up window 2-4 y", "Follow-up 2–4 y"),
            ("Follow-up window 4-end y", "Follow-up ≥4 y")]


def rows_for(outcome, ex, lvl):
    d = k[(k.outcome == outcome) & (k.exposure == ex) & (k.level == lvl)].set_index("analysis")
    return [(lab, d.at[a, "HR"], d.at[a, "lo"], d.at[a, "hi"], d.at[a, "p"]) for a, lab in ANALYSES if a in d.index]


# ------------------------------------------------------------------ Figure 3: Parkinson's disease, bidirectional
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6), sharey=True)
forest(axes[0], rows_for("Parkinson's disease", "hard_cat", "Often+"), "a  Hard stools often+ and Parkinson's disease", (0.2, 5))
forest(axes[1], rows_for("Parkinson's disease", "loose_cat", "Often+"), "b  Loose stools often+ and Parkinson's disease", (0.2, 5))
fig.tight_layout(w_pad=7)
save(fig, "Fig3_parkinson_bidirectional")

# ------------------------------------------------------------------ Figure 4: bowel frequency and cardiometabolic / mental outcomes
OUTS = ["Type 2 diabetes", "Heart failure", "Hypertension", "Stroke", "Ischaemic heart disease", "Atrial fibrillation",
        "Chronic kidney disease", "Depression", "Anxiety", "All-cause mortality"]
main = k[k.analysis == "Model 2: fully adjusted (main)"]
fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.6), sharey=True)
for ax, (lvl, title) in zip(axes, [("<1/day", "a  <1 per day"), ("2/day", "b  2 per day"), (">=3/day", "c  ≥3 per day")]):
    d = main[(main.exposure == "bmf_cat") & (main.level == lvl)].set_index("outcome")
    forest(ax, [(o, d.at[o, "HR"], d.at[o, "lo"], d.at[o, "hi"], d.at[o, "p"]) for o in OUTS], title, (0.5, 2.5))
    for t in ax.texts:
        t.set_visible(lvl == ">=3/day")
fig.suptitle("Bowel movement frequency (reference: 1 per day), fully adjusted", x=0.01, ha="left", fontsize=8, color=INK2)
fig.tight_layout(w_pad=1.0)
save(fig, "Fig4_frequency_cardiometabolic")

# ------------------------------------------------------------------ Figure 5: dose-response
dr = pd.read_csv(P + "results/dose_response.csv")
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.7))
specs = [("hard5", ["Parkinson's disease", "Depression"], "a  Hard/lumpy stool frequency"),
         ("loose5", ["Parkinson's disease", "Dementia", "Depression"], "b  Loose/watery stool frequency"),
         ("bmf6", ["Type 2 diabetes", "Heart failure", "Depression", "All-cause mortality"], "c  Bowel movement frequency")]
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]  # categorical slots 1-4, fixed order
LEVELS = {"hard5": ["Never", "Sometimes", "Often", "Most of the time", "Always"],
          "loose5": ["Never", "Sometimes", "Often", "Most of the time", "Always"],
          "bmf6": ["<=2/week", "3-6/week", "1/day", "2/day", "3/day", ">=4/day"]}
TICKS = {"hard5": ["Never", "Some-\ntimes", "Often", "Mostly", "Always"], "loose5": ["Never", "Some-\ntimes", "Often", "Mostly", "Always"],
         "bmf6": ["≤2/\nweek", "3–6/\nweek", "1/\nday", "2/\nday", "3/\nday", "≥4/\nday"]}
for ax, (ex, outs, title) in zip(axes, specs):
    lv = LEVELS[ex]
    for i, o in enumerate(outs):
        d = dr[(dr.exposure == ex) & (dr.outcome == o)].set_index("level").reindex(lv)
        x = np.arange(len(lv)) + (i - (len(outs) - 1) / 2) * 0.12
        ax.errorbar(x, d.HR, yerr=[d.HR - d.lo, d.hi - d.HR], fmt="o-", color=SERIES[i], ms=4, lw=1.2,
                    elinewidth=1, capsize=0, mec="white", mew=0.5, label=o)
    ax.axhline(1, color=INK2, lw=0.7)
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_yticks([0.1, 0.25, 0.5, 1, 2, 4, 8] if ex != "bmf6" else [0.5, 0.75, 1, 1.5, 2])
    ax.set_xticks(np.arange(len(lv)))
    ax.set_xticklabels(TICKS[ex], fontsize=6.5)
    ax.set_title(title, loc="left", fontsize=8, color=INK, fontweight="bold")
    ax.grid(axis="y", color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    if ex == "bmf6":
        ax.legend(fontsize=6, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=2)
    else:
        ax.legend(fontsize=6, frameon=False, loc="upper left")
axes[0].set_ylabel("Hazard ratio (95% CI)", fontsize=7.5)
fig.tight_layout(w_pad=1.2)
save(fig, "Fig5_dose_response")
print("figures saved")
