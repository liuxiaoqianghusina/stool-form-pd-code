# Main and supplementary figures for the Parkinson's-disease-focused manuscript (npj Parkinson's Disease version).
#   Fig3_pd_transit_spectrum : dose-response for hard and loose stools + transit spectrum score (spline and categories)
#   Fig4_pd_robustness       : hard and loose stools across models, time windows, sensitivity analyses, outcome
#                              definitions, missing-data handling; negative-control exposures
#   Fig5_pd_prediction       : incremental discrimination, case enrichment, joint stool form x other prodromal markers
#   FigS2_pd_subgroups       : subgroup forest plot with interaction P values
# Colour encodes direction only: red = hard stools / higher risk, blue = loose stools / lower risk, grey = neutral.
import sys

import matplotlib.pyplot as plt
import numpy as np
from fig_tiff import save_tiff
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

P = "/path/to/project/"
R = P + "results/"
RED, BLUE, GREY, DGREY = "#e34948", "#2a78d6", "#a9a7a1", "#6b6a66"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e4df"
plt.rcParams.update({"font.size": 7.5, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "font.family": "Arial"})
fmt_x = FuncFormatter(lambda v, _: f"{v:g}")


def save(fig, name):
    fig.savefig(P + f"figures/{name}.pdf", bbox_inches="tight")
    fig.savefig(P + f"figures/{name}.png", dpi=300, bbox_inches="tight")
    save_tiff(fig, P + f"figures/{name}.tif")


def panel(ax, letter, title):
    ax.set_title(f"{letter}  {title}", loc="left", fontsize=8.5, color=INK, fontweight="bold")


def logy(ax, ticks):
    ax.set_yscale("log"); ax.yaxis.set_major_locator(FixedLocator(ticks)); ax.yaxis.set_major_formatter(fmt_x)
    ax.yaxis.set_minor_formatter(NullFormatter()); ax.axhline(1, color=INK2, lw=0.8, zorder=1)
    for t in ticks:
        ax.axhline(t, color=GRID, lw=0.5, zorder=0)


# ================================================================== Fig 3: dose-response and transit spectrum
def fig3():
    dr = pd.read_csv(R + "dose_response.csv")
    spl = pd.read_csv(R + "tier1_transit_spline.csv")
    sc = pd.read_csv(R + "tier1_transit_score.csv")
    cnt = pd.read_csv(R + "tier1_transit_score_counts.csv")
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.55), gridspec_kw={"width_ratios": [1, 1, 1.45], "wspace": 0.38})
    LV = ["Never", "Sometimes", "Often", "Most of the time", "Always"]
    LAB = ["Never", "Sometimes", "Often", "Most of the time", "Always"]
    for ax, ex, col, letter, title in [(axes[0], "hard5", RED, "a", "Hard/lumpy stools"),
                                       (axes[1], "loose5", BLUE, "b", "Loose/watery stools")]:
        for outc, c, mk, lab, off in [("Parkinson's disease", col, "o", "Parkinson's disease", -0.08),
                                      ("Dementia", DGREY, "s", "All-cause dementia", 0.08)]:
            d = dr[(dr.exposure == ex) & (dr.outcome == outc)].set_index("level").reindex(LV)
            x = np.arange(5) + off
            ax.errorbar(x, d.HR, yerr=[d.HR - d.lo, d.hi - d.HR], fmt=mk + "-", color=c, ms=4.2, lw=1.4, elinewidth=1,
                        capsize=0, mec="white", mew=0.5, label=lab, zorder=3)
        logy(ax, [0.1, 0.25, 0.5, 1, 2, 4, 8, 16])
        ax.set_ylim(0.06, 18) if ex == "loose5" else ax.set_ylim(0.5, 18)
        ax.set_xticks(range(5)); ax.set_xticklabels(LAB, fontsize=6.3, rotation=35, ha="right", rotation_mode="anchor")
        ax.set_xlim(-0.4, 4.4)
        ev = d.events.astype(int)
        panel(ax, letter, title)
        p_tr = dr[(dr.exposure == ex) & (dr.outcome == "Parkinson's disease")].p_trend.iat[0]
        per = dr[(dr.exposure == ex) & (dr.outcome == "Parkinson's disease")].HR_per_level.iat[0]
        mant, expo = f"{p_tr:.0e}".split("e")
        ax.text(0.03, 0.97 if ex == "hard5" else 0.05, f"PD: HR {per:.2f} per level\n" + rf"$P$-trend = ${mant}\times10^{{{int(expo)}}}$",
                transform=ax.transAxes, fontsize=6.3, color=col, va="top" if ex == "hard5" else "bottom")
    axes[0].set_ylabel("Hazard ratio (95% CI)")
    axes[1].legend(frameon=False, fontsize=6.3, loc="upper right", handlelength=1.6)
    # c: transit spectrum
    ax = axes[2]
    ax.fill_between(spl.score, spl.lo, spl.hi, color="#d9d7d2", lw=0, zorder=1)
    ax.plot(spl.score, spl.HR, color=INK, lw=1.4, zorder=2)
    cat = sc[sc.analysis == "categorical"].copy()
    pos = {"score_cat<=-3": -3.5, "score_cat-2": -2, "score_cat-1": -1, "score_cat0": 0, "score_cat1": 1,
           "score_cat2": 2, "score_cat>=3": 3.5}
    for _, r in cat.iterrows():
        x = pos[r.term]; c = BLUE if x < 0 else (RED if x > 0 else INK2)
        ax.errorbar(x, r.HR, yerr=[[r.HR - r.lo], [r.hi - r.HR]], fmt="o", color=c, ms=4.5, elinewidth=1, mec="white",
                    mew=0.5, zorder=4)
    logy(ax, [0.125, 0.25, 0.5, 1, 2, 4, 8])
    ax.set_ylim(0.09, 9)
    ax.set_xticks([-3.5, -2, -1, 0, 1, 2, 3.5]); ax.set_xticklabels(["≤−3", "−2", "−1", "0", "1", "2", "≥3"])
    ax.set_xlim(-4.3, 4.3)
    ax.set_xlabel("Transit spectrum score (hard-stool level − loose-stool level)", fontsize=6.8)
    ax.text(-4.2, 0.16, "← looser", color=BLUE, fontsize=6.5, va="bottom", bbox=dict(fc="white", ec="none", pad=0.5)); ax.text(4.2, 0.16, "harder →", color=RED, fontsize=6.5, ha="right", va="bottom")
    lin = sc[sc.analysis == "per 1-point increase"].iloc[0]
    pnl = pd.read_csv(R + "tier1_transit_info.csv").p_nonlinearity.iat[0]
    ax.text(-4.1, 7.2, f"HR per point {lin.HR:.2f} ({lin.lo:.2f}–{lin.hi:.2f})\n" + rf"$P$-nonlinearity = {pnl:.2f}", fontsize=6.3, va="top")
    n_by = dict(zip(cnt.score, cnt.n))
    ax2 = ax.inset_axes([0.0, -0.005, 1, 0.13], sharex=ax)
    ax2.bar(cnt.score, cnt.n, width=0.6, color=["#9cc0ee" if s < 0 else ("#f3a9a8" if s > 0 else "#cfcdc8") for s in cnt.score])
    ax2.axis("off")
    panel(ax, "c", "Transit spectrum and Parkinson's disease")
    save(fig, "Fig3_pd_transit_spectrum")


# ================================================================== Fig 4: robustness
def fig4():
    k = pd.read_csv(R + "key_outcomes.csv"); k = k[k.outcome == "Parkinson's disease"]
    fg = pd.read_csv(R + "key_outcomes_finegray.csv"); fg = fg[fg.outcome == "Parkinson's disease"]
    ma = pd.read_csv(R + "pd_mutual_adjustment.csv")
    de = pd.read_csv(R + "tier1_pd_definitions.csv")
    med = pd.read_csv(R + "tier2_medication.csv") if __import__("os").path.exists(R + "tier2_medication.csv") else None
    mi = pd.read_csv(R + "tier2_missing_data.csv") if __import__("os").path.exists(R + "tier2_missing_data.csv") else None
    neg = pd.read_csv(R + "tier1_negative_control.csv")
    slc = pd.read_csv(R + "sens_loose_causes.csv")

    def kk(analysis, ex):
        r = k[(k.analysis == analysis) & (k.exposure == ex) & (k.level == "Often+")].iloc[0]; return r.HR, r.lo, r.hi
    def term(df, mask):
        r = df[mask].iloc[0]; return r.HR, r.lo, r.hi
    groups = [
        ("Model", [("Age and sex", lambda e: kk("Model 1: age + sex", e)),
                   ("Fully adjusted (main)", lambda e: kk("Model 2: fully adjusted (main)", e)),
                   ("+ coffee and tea", lambda e: kk("Model 3: + coffee and tea", e)),
                   ("+ other stool form and frequency", lambda e: term(ma, (ma.analysis == "Mutually adjusted: hard + loose + frequency") & (ma.term == f"{e}Often+")))]),
        ("Timing", [("Events <2 y excluded", lambda e: kk("Excluding events in first 2 years", e)),
                    ("Follow-up 0–2 y", lambda e: kk("Follow-up window 0-2 y", e)),
                    ("Follow-up 2–4 y", lambda e: kk("Follow-up window 2-4 y", e)),
                    ("Follow-up ≥4 y", lambda e: kk("Follow-up window 4-end y", e))]),
        ("Confounding and selection", [("IBS excluded", lambda e: kk("Excluding IBS", e)),
                                       ("Laxative, metformin, opioid users excluded", lambda e: kk("Excluding laxative, metformin and opioid users", e))]
         + ([("Drugs affecting parkinsonism excluded",
              lambda e: term(med, (med.analysis == "Excluding all three drug groups") & (med.term == f"{e}Often+")))] if med is not None else [])
         + [("Causes of loose stools excluded", lambda e: term(slc, slc.term == f"{e}Often+"))]
         + ([("Multiple imputation", lambda e: term(mi, (mi.analysis.str.startswith("Multiple")) & (mi.term == f"{e}Often+"))),
             ("Complete cases", lambda e: term(mi, (mi.analysis == "Complete-case analysis") & (mi.term == f"{e}Often+")))] if mi is not None else [])
         + [("Death as competing risk (Fine–Gray)", lambda e: term(fg, (fg.exposure == e) & (fg.level == "Often+")))]),
        ("Outcome definition", [("UK Biobank algorithmic PD", lambda e: term(de, (de.analysis.str.startswith("UK Biobank")) & (de.term == f"{e}Often+"))),
                                ("Hospital or death records only", lambda e: term(de, (de.analysis.str.startswith("First report")) & (de.term == f"{e}Often+"))),
                                ("Secondary/atypical parkinsonism excluded", lambda e: term(de, (de.analysis.str.startswith("Excluding cases")) & (de.term == f"{e}Often+")))]),
    ]
    rows = []
    for g, items in groups:
        rows.append(("header", g, None, None))
        for lab, fn in items:
            rows.append(("row", lab, fn("hard_cat"), fn("loose_cat")))
    n = len(rows)
    fig = plt.figure(figsize=(7.2, 0.17 * n + 1.9))
    gs = fig.add_gridspec(2, 3, height_ratios=[n, 4.2], width_ratios=[1.35, 1, 1], hspace=0.35, wspace=0.08)
    axl = fig.add_subplot(gs[0, 0]); axh = fig.add_subplot(gs[0, 1]); axs = fig.add_subplot(gs[0, 2])
    y = np.arange(n)[::-1]
    axl.set_ylim(-0.6, n - 0.4); axl.axis("off")
    for yi, (kind, lab, h, s) in zip(y, rows):
        axl.text(0 if kind == "header" else 0.04, yi, lab, fontsize=6.6 if kind == "row" else 7,
                 fontweight="bold" if kind == "header" else "normal", va="center", color=INK if kind == "header" else INK2,
                 transform=axl.get_yaxis_transform())
    for ax, idx, col, title in [(axh, 2, RED, "Hard stools, often or more"), (axs, 3, BLUE, "Loose stools, often or more")]:
        for yi, r in zip(y, rows):
            if r[0] != "row":
                continue
            hr, lo, hi = r[idx]
            main = r[1] == "Fully adjusted (main)"
            ax.plot([lo, hi], [yi, yi], color=col, lw=1.1, solid_capstyle="round")
            ax.plot(hr, yi, marker="D" if main else "o", color=col, ms=4.6 if main else 3.6, mec="white", mew=0.5)
            ax.text(5.6 if col == RED else 1.45, yi, f"{hr:.2f} ({lo:.2f}–{hi:.2f})", fontsize=5.9, va="center", color=INK2)
        ax.set_xscale("log"); ax.axvline(1, color=INK2, lw=0.8)
        ax.set_xlim((0.8, 24) if col == RED else (0.12, 4.2))
        ax.xaxis.set_major_locator(FixedLocator([1, 2, 4] if col == RED else [0.125, 0.25, 0.5, 1]))
        ax.xaxis.set_major_formatter(fmt_x); ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_ylim(-0.6, n - 0.4); ax.set_yticks([]); ax.spines["left"].set_visible(False)
        ax.set_xlabel("Hazard ratio (95% CI)", fontsize=6.8)
        ax.set_title(title, fontsize=7.5, color=col, fontweight="bold")
    fig.text(0.125, 0.955, "a  Hard and loose stools and Parkinson's disease across analyses", fontsize=8.5, fontweight="bold", va="top")
    # b: negative-control exposures alongside stool form
    axn = fig.add_subplot(gs[1, :])
    nn = neg[neg.analysis == "Single exposure"]
    k2 = k[k.analysis == "Model 2: fully adjusted (main)"]
    items = [("Hard stools, often or more", *k2[(k2.exposure == "hard_cat") & (k2.level == "Often+")][["HR", "lo", "hi"]].iloc[0], RED),
             ("Loose stools, often or more", *k2[(k2.exposure == "loose_cat") & (k2.level == "Often+")][["HR", "lo", "hi"]].iloc[0], BLUE),
             ("Abdominal pain ≥1 day/week", *nn[nn.term == "abd_pain_freqWeekly or more"][["HR", "lo", "hi"]].iloc[0], GREY),
             ("Current abdominal pain", *nn[nn.term == "abd_pain_nowYes"][["HR", "lo", "hi"]].iloc[0], GREY),
             ("Current abdominal distension", *nn[nn.term == "distension_nowYes"][["HR", "lo", "hi"]].iloc[0], GREY)]
    xs = np.arange(len(items))
    for x, (lab, hr, lo, hi, c) in zip(xs, items):
        axn.plot([x, x], [lo, hi], color=c, lw=1.4); axn.plot(x, hr, "o", color=c, ms=5, mec="white", mew=0.5)
        axn.text(x + 0.08, hr, f"{hr:.2f}", fontsize=6.2, va="center", color=INK2)
    logy(axn, [0.25, 0.5, 1, 2, 4]); axn.set_ylim(0.2, 5)
    axn.set_xticks(xs); axn.set_xticklabels([i[0].replace(", ", ",\n").replace("pain ≥", "pain\n≥").replace("Current abdominal ", "Current abdominal\n")
                                             for i in items], fontsize=6.4); axn.set_xlim(-0.5, len(items) - 0.5)
    axn.set_ylabel("HR (95% CI)", fontsize=6.8)
    panel(axn, "b", "Stool form versus negative-control symptoms from the same questionnaire")
    save(fig, "Fig4_pd_robustness")


# ================================================================== Fig 5: prediction and prodromal markers
def fig5():
    pf = pd.read_csv(R + "pd_prediction_perf.csv").set_index("model")
    cv = pd.read_csv(R + "pd_prediction_cv.csv").set_index("model").C_cv
    en = pd.read_csv(R + "pd_prediction_enrich.csv")
    jc = R + "tier2_joint_cuminc.csv"
    has_joint = __import__("os").path.exists(jc)
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.0), gridspec_kw={"wspace": 0.32, "hspace": 0.62})
    axes = axs.ravel()
    ax = axes[0]
    labs = ["Conventional\nmodel", "+ bowel\nfrequency", "+ stool\nform", "+ both"]
    app = [pf.C_base.iat[0], pf.loc["freq", "C_model"], pf.loc["form", "C_model"], pf.loc["all", "C_model"]]
    cvs = [cv["base"], np.nan, cv["form"], cv["all"]]
    cols = [GREY, DGREY, INK, INK]
    for i, (a, c, col) in enumerate(zip(app, cvs, cols)):
        ax.plot(i, a, "o", color=col, ms=5.5, mec="white", mew=0.5, zorder=3)
        if not np.isnan(c):
            ax.plot(i, c, "o", mfc="white", mec=col, ms=5, mew=1.1, zorder=3)
    for i, kname in [(1, "freq"), (2, "form")]:
        ax.text(i, app[i] + 0.004, f"ΔC {pf.loc[kname, 'dC']:+.3f}\n({pf.loc[kname, 'dC_lo']:.3f} to {pf.loc[kname, 'dC_hi']:.3f})",
                ha="center", fontsize=5.8, color=INK2, va="bottom")
    ax.set_xticks(range(4)); ax.set_xticklabels(labs, fontsize=6.3); ax.set_xlim(-0.5, 3.5); ax.set_ylim(0.735, 0.80)
    ax.set_ylabel("C index"); ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    ax.plot([], [], "o", color=INK, ms=4.5, label="apparent"); ax.plot([], [], "o", mfc="white", mec=INK, ms=4.5, label="10-fold CV")
    ax.legend(frameon=False, fontsize=6, loc="lower right", handletextpad=0.2)
    panel(ax, "a", "Discrimination")
    ax = axes[1]
    for j, (q, lab) in enumerate([(0.05, "Top 5%"), (0.10, "Top 10%")]):
        for i, (m, col) in enumerate([("base", GREY), ("form", INK)]):
            r = en[(en.model == m) & (en.top_fraction == q)].iloc[0]
            x = j + (i - 0.5) * 0.36
            ax.bar(x, r.pct_cases_captured, width=0.34, color=col)
            ax.text(x, r.pct_cases_captured + 0.6, f"{r.pct_cases_captured:.1f}", ha="center", fontsize=6, color=INK2)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Top 5% of\npredicted risk", "Top 10% of\npredicted risk"], fontsize=6.3)
    ax.set_ylabel("PD cases captured (%)"); ax.set_ylim(0, 46)
    ax.legend(handles=[Patch(color=GREY, label="Conventional model"), Patch(color=INK, label="+ stool form")],
              frameon=False, fontsize=6, loc="upper left")
    panel(ax, "b", "Case enrichment")
    if has_joint:
        ax = axes[2]
        j = pd.read_csv(jc)
        order = [("Loose often+", BLUE), ("Neither", GREY), ("Hard often+", RED)]
        for gi, (form, col) in enumerate(order):
            for pi, (pm, alpha) in enumerate([("No other prodromal marker", 0.45), ("Any other prodromal marker", 1.0)]):
                r = j[j.joint == f"{form} | {pm}"]
                if r.empty:
                    continue
                r = r.iloc[0]; x = gi + (pi - 0.5) * 0.38
                ax.bar(x, r.cuminc5, width=0.36, color=col, alpha=alpha)
                ax.plot([x, x], [r.lo, r.hi], color=INK2, lw=0.8)
                ax.text(x, r.hi + 0.04, f"{r.cuminc5:.2f}", ha="center", fontsize=5.8, color=INK2)
        ax.set_xticks(range(3)); ax.set_xticklabels(["Loose stools\noften+", "Neither", "Hard stools\noften+"], fontsize=6.3)
        ax.set_ylabel("5-year cumulative incidence of PD (%)")
        ax.legend(handles=[Patch(color=DGREY, alpha=0.45, label="No other prodromal marker"),
                           Patch(color=DGREY, alpha=1, label="≥1 other prodromal marker")], frameon=False, fontsize=6, loc="upper left")
        panel(ax, "c", "Stool form and other prodromal markers")
    # d: accelerometry sub-study, joint low acceleration x hard stools
    ax = axes[3]
    aj = pd.read_csv(R + "tier2_accel_joint.csv"); ajc = pd.read_csv(R + "tier2_accel_joint_counts.csv")
    info = pd.read_csv(R + "tier2_accel_info.csv").iloc[0]
    groups = [("Higher activity\nhard stools −", None, GREY, "acc_low=0, hard_often=False"),
              ("Higher activity\nhard stools +", "interaction(acc_low, hard_stool >= 2)0.TRUE", "#f3a9a8", "acc_low=0, hard_often=True"),
              ("Low activity\nhard stools −", "interaction(acc_low, hard_stool >= 2)1.FALSE", DGREY, "acc_low=1, hard_often=False"),
              ("Low activity\nhard stools +", "interaction(acc_low, hard_stool >= 2)1.TRUE", RED, "acc_low=1, hard_often=True")]
    for i, (lab, t, col, ck) in enumerate(groups):
        if t is None:
            hr, lo, hi = 1, 1, 1
        else:
            r = aj[aj.term == t].iloc[0]; hr, lo, hi = r.HR, r.lo, r.hi
        ax.plot([i, i], [lo, hi], color=col, lw=1.4); ax.plot(i, hr, "o", color=col, ms=6, mec="white", mew=0.5)
        c = ajc[ajc.term.str.lower() == ck.lower()]
        ev = f"{int(c.events.iat[0])}/{int(c.n.iat[0]):,}" if not c.empty else ""
        ax.text(i + 0.12, hr, f"{hr:.2f}" + (f"\n{ev}" if ev else ""), fontsize=5.8, va="center", color=INK2)
    logy(ax, [0.5, 1, 2, 4, 8, 16]); ax.set_ylim(0.5, 16)
    ax.set_xticks(range(4)); ax.set_xticklabels([g[0] for g in groups], fontsize=6.2); ax.set_xlim(-0.4, 3.75)
    ax.text(0.03, 0.97, "Low activity: lowest quintile of average acceleration\nHard stools +: often or more", transform=ax.transAxes,
            fontsize=5.6, va="top", color=INK2)
    ax.set_ylabel("HR for Parkinson's disease (95% CI)", fontsize=6.8)
    panel(ax, "d", f"Wrist accelerometry sub-study (n={int(info.n):,})")
    save(fig, "Fig5_pd_prediction")


# ================================================================== Supplementary: subgroups
def figS2():
    s = pd.read_csv(R + "tier1_subgroups.csv"); s = s[s.term.str.contains("Often")]
    SGL = {"sg_sex": "Sex", "sg_age": "Age", "sg_smoke": "Smoking", "sg_coffee": "Coffee", "sg_bmi": "Body mass index",
           "sg_lax": "Laxative use", "sg_dep": "Depression history"}
    rows = []
    for g in SGL:
        h = s[(s.subgroup == g) & (s.exposure == "hard_cat")]; l = s[(s.subgroup == g) & (s.exposure == "loose_cat")]
        rows.append(("header", SGL[g], h.p_interaction.iat[0], l.p_interaction.iat[0], None))
        for lv in h.lvl:
            rows.append(("row", lv.replace(">=", "≥").replace("kg/m2", "kg/m²"), h[h.lvl == lv].iloc[0], l[l.lvl == lv].iloc[0], None))
    n = len(rows); y = np.arange(n)[::-1]
    fig, (axl, axh, axs) = plt.subplots(1, 3, figsize=(7.2, 0.17 * n + 0.7), gridspec_kw={"width_ratios": [1.1, 1, 1], "wspace": 0.06})
    axl.axis("off"); axl.set_ylim(-0.6, n - 0.4)
    for yi, r in zip(y, rows):
        if r[0] == "header":
            axl.text(0, yi, r[1], fontsize=7, fontweight="bold", va="center")
        else:
            axl.text(0.05, yi, f"{r[1]}  (n={int(r[2].n):,}; {int(r[2].events)} cases)", fontsize=6.4, va="center", color=INK2)
    for ax, idx, col, title, xl, ticks in [(axh, 2, RED, "Hard stools, often or more", (0.15, 12), [0.25, 0.5, 1, 2, 4, 8]),
                                           (axs, 3, BLUE, "Loose stools, often or more", (0.1, 4), [0.125, 0.25, 0.5, 1, 2])]:
        for yi, r in zip(y, rows):
            if r[0] == "header":
                ax.text(xl[1], yi, f"P-int = {r[idx]:.2g}", fontsize=6, ha="right", va="center", color=INK2)
                continue
            d = r[idx]
            ax.plot([d.lo, d.hi], [yi, yi], color=col, lw=1.1); ax.plot(d.HR, yi, "o", color=col, ms=3.6, mec="white", mew=0.5)
        ax.set_xscale("log"); ax.axvline(1, color=INK2, lw=0.8); ax.set_xlim(*xl)
        ax.xaxis.set_major_locator(FixedLocator(ticks)); ax.xaxis.set_major_formatter(fmt_x); ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_ylim(-0.6, n - 0.4); ax.set_yticks([]); ax.spines["left"].set_visible(False)
        ax.set_xlabel("Hazard ratio (95% CI)", fontsize=6.8); ax.set_title(title, fontsize=7.5, color=col, fontweight="bold")
    save(fig, "FigS2_pd_subgroups")


# ================================================================== Supplementary: genetics and imaging
def figS3():
    ss = pd.read_csv(R + "mr_singlesnp.csv"); ss = ss[(ss.outcome == "Parkinson's disease") & ss.SNP.str.startswith("rs")]
    mr = pd.read_csv(R + "mr_stoolfreq_finngen.csv"); mr = mr[mr.outcome == "Parkinson's disease"]
    co = pd.read_csv(R + "coloc_stoolfreq_pd.csv")
    q = pd.read_csv(R + "qsm_nigra.csv"); q = q[(q.model == "Model 2") & (q.level == "Often+")]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.7), gridspec_kw={"width_ratios": [1.15, 1, 1.15], "wspace": 0.55})
    ax = axes[0]
    ss = ss.sort_values("b")
    y = np.arange(len(ss))
    for yi, (_, r) in zip(y, ss.iterrows()):
        c = RED if r.SNP == "rs2732706" else INK2
        ax.plot([r.b - 1.96 * r.se, r.b + 1.96 * r.se], [yi, yi], color=c, lw=1); ax.plot(r.b, yi, "o", color=c, ms=3.5)
    ax.set_yticks(y); ax.set_yticklabels([s_ + (" (MAPT)" if s_ == "rs2732706" else "") for s_ in ss.SNP], fontsize=5.8)
    ivw = [("IVW, all", "Inverse variance weighted", "All instruments"), ("IVW, MAPT removed", "Inverse variance weighted", "Excluding MAPT-region variant"),
           ("MR-PRESSO corrected", "MR-PRESSO Outlier-corrected", "All instruments")]
    for k_, (lab, m_, set_) in enumerate(ivw):
        r = mr[(mr.method == m_) & (mr.instrument_set == set_)].iloc[0]; yy = len(ss) + 0.8 + k_
        ax.plot([r.b - 1.96 * r.se, r.b + 1.96 * r.se], [yy, yy], color=BLUE, lw=1.3); ax.plot(r.b, yy, "D", color=BLUE, ms=4)
    ax.set_yticks(list(y) + [len(ss) + 0.8 + k_ for k_ in range(3)])
    ax.set_yticklabels([s_ + (" (MAPT)" if s_ == "rs2732706" else "") for s_ in ss.SNP] + [i[0] for i in ivw], fontsize=5.8)
    ax.axvline(0, color=INK2, lw=0.8); ax.set_xlabel("log OR for PD per unit stool frequency", fontsize=6.5)
    panel(ax, "a", "Mendelian randomization")
    ax = axes[1]
    co = co.sort_values("PP_H4")
    ax.barh(range(len(co)), co.PP_H4, color=[RED if l == "rs2732706" else GREY for l in co.lead])
    ax.barh(range(len(co)), co.PP_H3, left=co.PP_H4, color="#e7e5e0")
    ax.set_yticks(range(len(co))); ax.set_yticklabels([l + (" (MAPT)" if l == "rs2732706" else "") for l in co.lead], fontsize=5.8)
    ax.axvline(0.8, color=INK2, lw=0.6, ls="--"); ax.set_xlim(0, 1); ax.set_xlabel("Posterior probability", fontsize=6.5)
    ax.legend(handles=[Patch(color=GREY, label="PP.H4 (shared variant)"), Patch(color="#e7e5e0", label="PP.H3 (distinct variants)")],
              frameon=False, fontsize=5.6, loc="center right")
    panel(ax, "b", "Colocalisation")
    ax = axes[2]
    REG = ["Substantia nigra", "Putamen", "Caudate", "Pallidum", "Thalamus", "Hippocampus", "Amygdala", "Accumbens"]
    for j, (ex, col, off) in enumerate([("hard_cat", RED, -0.15), ("loose_cat", BLUE, 0.15)]):
        for i, reg in enumerate(REG):
            r = q[(q.outcome == f"qsm_{reg}") & (q.exposure == ex)]
            if r.empty:
                continue
            r = r.iloc[0]; yy = len(REG) - 1 - i + off
            ax.plot([r.lo, r.hi], [yy, yy], color=col, lw=1); ax.plot(r.beta_sd, yy, "o", color=col, ms=3.3)
    ax.set_yticks(range(len(REG))); ax.set_yticklabels(REG[::-1], fontsize=6)
    ax.axvline(0, color=INK2, lw=0.8); ax.set_xlabel("QSM difference (SD), often+ vs never", fontsize=6.3)
    ax.plot([], [], "o", color=RED, label="Hard stools"); ax.plot([], [], "o", color=BLUE, label="Loose stools")
    ax.legend(frameon=False, fontsize=5.8, loc="center left", bbox_to_anchor=(1.0, 0.5))
    panel(ax, "c", "Brain iron (QSM)")
    save(fig, "FigS3_genetics_imaging")


if __name__ == "__main__":
    which = sys.argv[1:] or ["fig3", "fig4", "fig5", "figS2", "figS3"]
    for w in which:
        globals()[w]()
        print("saved", w)
