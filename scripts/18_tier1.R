# Tier-1 robustness and interpretation analyses for Parkinson's disease (PD), all with Model 2 covariates:
#  (1) negative-control exposures from the same questionnaire (abdominal pain and distension);
#  (2) alternative PD definitions (UK Biobank algorithmic PD; first reports from hospital/death records only;
#      cases with secondary or atypical parkinsonism codes not counted);
#  (3) a "transit spectrum" score (hard-stool level minus loose-stool level, -4 to +4), categorical and spline;
#  (4) subgroup analyses with interaction tests;
#  (5) absolute risk: crude Aalen-Johansen cumulative incidence (death as competing risk) and Cox-standardised
#      5-year risks with bootstrap CIs for risk differences.
# Usage: Rscript 18_tier1.R [B]  (bootstrap replicates for (5), default 200)
suppressPackageStartupMessages({library(parallel); library(splines)})
source("/path/to/project/scripts/prep.R")
B <- if (length(commandArgs(TRUE))) as.integer(commandArgs(TRUE)[1]) else 200L
set.seed(20261009)
R_ <- paste0(P, "results/")

ex <- fread(paste0(P, "data/extra_fields.csv"), colClasses = "character")
mx <- match(co$eid, as.integer(ex$eid))
num <- function(v) suppressWarnings(as.numeric(ex[[v]][mx]))
pain_f <- num("p21025")
co[, abd_pain_freq := factor(fcase(pain_f == 0, "Never", pain_f %in% 1:3, "Less than weekly", pain_f %in% 4:6, "Weekly or more"),
                             levels = c("Never", "Less than weekly", "Weekly or more"))]
co[, abd_pain_now := factor(fcase(num("p21035") == 0, "No", num("p21035") == 1, "Yes"), levels = c("No", "Yes"))]
co[, distension_now := factor(fcase(num("p21038") == 0, "No", num("p21038") == 1, "Yes"), levels = c("No", "Yes"))]

hr_rows <- function(m, pattern, label) {
  s <- summary(m); r <- grep(pattern, rownames(s$coefficients))
  data.table(analysis = label, term = rownames(s$coefficients)[r], HR = s$conf.int[r, 1], lo = s$conf.int[r, 3],
             hi = s$conf.int[r, 4], p = s$coefficients[r, 5], events = m$nevent, n = m$n)
}
pd_data <- function(d_days) {                      # analysis table for a PD definition given as days from questionnaire
  tt <- outcome_tt(d_days, FALSE)
  dd <- co[tt$sel]; dd[, `:=`(time = tt$time[tt$sel], event = tt$event[tt$sel])]; dd
}
G20 <- oc[["G20"]]
dd <- pd_data(G20)

# ---------------------------------------------------------------- (1) negative-control exposures
neg <- list()
for (e in c("abd_pain_freq", "abd_pain_now", "distension_now")) {
  m1 <- coxph(as.formula(paste("Surv(time, event) ~", e, "+", M2)), data = dd)
  m2 <- coxph(as.formula(paste("Surv(time, event) ~", e, "+ hard_cat + loose_cat +", M2)), data = dd)
  neg[[e]] <- rbind(hr_rows(m1, paste0("^", e), "Single exposure"), hr_rows(m2, paste0("^", e), "Adjusted for stool form"))
}
neg_tab <- rbindlist(neg)
cnt_neg <- rbindlist(lapply(c("abd_pain_freq", "abd_pain_now", "distension_now"), function(e)
  dd[!is.na(get(e)), .(n = .N, events = sum(event)), by = .(level = get(e))][, exposure := e]))
fwrite(neg_tab, paste0(R_, "tier1_negative_control.csv")); fwrite(cnt_neg, paste0(R_, "tier1_negative_control_counts.csv"))

# ---------------------------------------------------------------- (2) alternative PD definitions
def <- list()
fitdef <- function(d_days, label) {
  x <- pd_data(d_days)
  rbindlist(lapply(c("hard_cat", "loose_cat"), function(e)
    hr_rows(coxph(as.formula(paste("Surv(time, event) ~", e, "+", M2)), data = x), paste0("^", e), label)))
}
def[["main"]] <- fitdef(G20, "First occurrence G20 (main)")
# (a) UK Biobank algorithmically defined PD; algorithm data run to their last update, so later time is censored there
alg_date <- as.IDate(ex$p42032[mx])
special <- !is.na(alg_date) & year(alg_date) %in% c(1900, 1901, 1902, 1903, 2037)
alg_days <- as.numeric(alg_date - as.IDate(co$q_date)); alg_days[special] <- -99999
alg_end <- max(alg_date[!special], na.rm = TRUE)
fu_save <- co$fu_end_days
co[, fu_end_days := pmin(fu_end_days, as.numeric(alg_end - as.IDate(q_date)))]
def[["alg"]] <- fitdef(alg_days, "UK Biobank algorithmic PD (42032)")
co[, fu_end_days := fu_save]
# (b) first report from hospital admissions or death register (source codes 20, 21, 40, 41); other first reports are
#     not counted as events (follow-up censored at that report)
src <- ex$p131023[mx]
hosp <- src %in% c("20", "21", "40", "41")
x <- pd_data(G20)
keep_rows <- which(outcome_tt(G20, FALSE)$sel)
x[, event := as.integer(event == 1 & hosp[keep_rows])]
def[["hosp"]] <- rbindlist(lapply(c("hard_cat", "loose_cat"), function(e)
  hr_rows(coxph(as.formula(paste("Surv(time, event) ~", e, "+", M2)), data = x), paste0("^", e),
          "First report from hospital or death records")))
# (c) cases with any secondary parkinsonism (G21) or other degenerative basal-ganglia disease (G23) code not counted
atyp <- (!is.na(oc$G21) & oc$G21 <= co$fu_end_days) | (!is.na(oc$G23) & oc$G23 <= co$fu_end_days)
x <- pd_data(G20)
x[, event := as.integer(event == 1 & !atyp[keep_rows])]
def[["typ"]] <- rbindlist(lapply(c("hard_cat", "loose_cat"), function(e)
  hr_rows(coxph(as.formula(paste("Surv(time, event) ~", e, "+", M2)), data = x), paste0("^", e),
          "Excluding cases with G21/G23 codes")))
def_tab <- rbindlist(def)
def_info <- data.table(alg_cases_total = sum(!is.na(alg_date)), alg_last_date = as.character(alg_end),
                       src_hosp_or_death = sum(hosp & !is.na(G20) & G20 > 0), atypical_any = sum(atyp & !is.na(G20) & G20 > 0))
fwrite(def_tab, paste0(R_, "tier1_pd_definitions.csv")); fwrite(def_info, paste0(R_, "tier1_pd_definitions_info.csv"))

# ---------------------------------------------------------------- (3) transit spectrum score
dd[, score := hard_stool - loose_stool]
dd[, score_cat := factor(fcase(score <= -3, "<=-3", score == -2, "-2", score == -1, "-1", score == 0, "0", score == 1, "1",
                               score == 2, "2", score >= 3, ">=3"), levels = c("0", "<=-3", "-2", "-1", "1", "2", ">=3"))]
ds <- dd[!is.na(score)]
m_cat <- coxph(as.formula(paste("Surv(time, event) ~ score_cat +", M2)), data = ds)
m_lin <- coxph(as.formula(paste("Surv(time, event) ~ score +", M2)), data = ds)
nsb <- ns(ds$score, df = 3)
ds[, paste0("ns", 1:3) := as.data.table(unclass(nsb)[, 1:3])]
m_ns <- coxph(as.formula(paste("Surv(time, event) ~ ns1 + ns2 + ns3 +", M2)), data = ds)
grid <- seq(-4, 4, by = 0.1)
Xg <- predict(nsb, grid) - matrix(predict(nsb, 0), nrow = length(grid), ncol = 3, byrow = TRUE)
bb <- coef(m_ns)[paste0("ns", 1:3)]; V <- vcov(m_ns)[paste0("ns", 1:3), paste0("ns", 1:3)]
lp <- drop(Xg %*% bb); se <- sqrt(rowSums((Xg %*% V) * Xg))
spl <- data.table(score = grid, HR = exp(lp), lo = exp(lp - 1.96 * se), hi = exp(lp + 1.96 * se))
p_nonlin <- anova(m_lin, m_ns)$`Pr(>|Chi|)`[2]
cat_tab <- rbind(data.table(analysis = "categorical", term = "score_cat0", HR = 1, lo = 1, hi = 1, p = NA, events = m_cat$nevent, n = m_cat$n),
                 hr_rows(m_cat, "^score_cat", "categorical"), hr_rows(m_lin, "^score$", "per 1-point increase"))
cnt_score <- ds[, .(n = .N, events = sum(event)), by = score][order(score)]
# sensitivity: mixed pattern (often hard and often loose) removed
m_lin_nomix <- coxph(as.formula(paste("Surv(time, event) ~ score +", M2)), data = ds[!(hard_stool >= 2 & loose_stool >= 2)])
cat_tab <- rbind(cat_tab, hr_rows(m_lin_nomix, "^score$", "per 1-point increase, mixed pattern excluded"))
fwrite(cat_tab, paste0(R_, "tier1_transit_score.csv")); fwrite(spl, paste0(R_, "tier1_transit_spline.csv"))
fwrite(cnt_score, paste0(R_, "tier1_transit_score_counts.csv"))
fwrite(data.table(p_nonlinearity = p_nonlin, n_mixed = sum(ds$hard_stool >= 2 & ds$loose_stool >= 2)),
       paste0(R_, "tier1_transit_info.csv"))

# ---------------------------------------------------------------- (4) subgroups with interaction tests
dd[, `:=`(sg_sex = factor(fifelse(sex == 1, "Men", "Women"), levels = c("Women", "Men")),
          sg_age = factor(fifelse(age_q < 65, "<65 years", ">=65 years")),
          sg_smoke = factor(fifelse(as.character(smoking) == "0", "Never smoker", "Ever smoker"), levels = c("Never smoker", "Ever smoker")),
          sg_coffee = factor(fifelse(coffee_cups >= 2, ">=2 cups/day", "<2 cups/day"), levels = c("<2 cups/day", ">=2 cups/day")),
          sg_bmi = factor(fifelse(bmi >= 25, ">=25 kg/m2", "<25 kg/m2"), levels = c("<25 kg/m2", ">=25 kg/m2")),
          sg_lax = factor(fifelse(laxative_any == 1, "Laxative use", "No laxative use"), levels = c("No laxative use", "Laxative use")),
          sg_dep = factor(fifelse(prev_depression == 1, "Depression history", "No depression history"),
                          levels = c("No depression history", "Depression history")))]
SG <- c(sg_sex = "sex", sg_age = NA, sg_smoke = "smoking", sg_coffee = NA, sg_bmi = NA, sg_lax = "laxative_any", sg_dep = "prev_depression")
drop_term <- function(f, term) if (is.na(term)) f else gsub(paste0("(^| )", term, " \\+ "), "\\1", gsub(paste0(" \\+ ", term, "( |$)"), "\\1", f))
sub_rows <- list()
for (g in names(SG)) {
  cov <- drop_term(M2, SG[[g]])
  for (e in c("hard_cat", "loose_cat")) {
    m0 <- coxph(as.formula(paste("Surv(time, event) ~", g, "+", e, "+", cov)), data = dd)
    m1 <- coxph(as.formula(paste("Surv(time, event) ~", g, "+", g, ":", e, "+", cov)), data = dd)
    p_int <- anova(m0, m1)$`Pr(>|Chi|)`[2]
    s <- summary(m1); r <- grep(paste0("^", g, ".*:", e), rownames(s$coefficients))
    ev <- dd[!is.na(get(e)), .(n = .N, events = sum(event)), by = .(lvl = as.character(get(g)))]
    sub_rows[[paste(g, e)]] <- data.table(subgroup = g, exposure = e, term = rownames(s$coefficients)[r],
                                          HR = s$conf.int[r, 1], lo = s$conf.int[r, 3], hi = s$conf.int[r, 4],
                                          p = s$coefficients[r, 5], p_interaction = p_int)[, lvl := sub(paste0("^", g), "", sub(":.*", "", term))][
                                            ev, on = "lvl", `:=`(n = i.n, events = i.events)]
  }
}
sub_tab <- rbindlist(sub_rows)
fwrite(sub_tab, paste0(R_, "tier1_subgroups.csv"))

# ---------------------------------------------------------------- (5) absolute risk
dth <- oc$DEATH_all[keep_rows]
dd[, status3 := factor(fifelse(event == 1, 1L, fifelse(!is.na(dth) & dth > 0 & dth / 365.25 <= time + 1e-9, 2L, 0L)),
                       levels = 0:2, labels = c("censor", "PD", "death"))]
aj <- rbindlist(lapply(c("hard_cat", "loose_cat"), function(e) {
  f <- survfit(as.formula(paste("Surv(time, status3) ~", e)), data = dd[!is.na(get(e))])
  s <- summary(f, times = 5, extend = TRUE)
  k <- which(f$states == "PD")
  data.table(exposure = e, level = sub(".*=", "", s$strata), cuminc5 = 100 * s$pstate[, k],
             lo = 100 * s$lower[, k], hi = 100 * s$upper[, k])
}))
slim <- dd[, c("time", "event", "hard_cat", "loose_cat", "age_q", all.vars(as.formula(paste("~", M2)))), with = FALSE]
slim <- slim[, unique(names(slim)), with = FALSE]
rm(oc, x); invisible(gc())
std_risk <- function(data, e, sub = NULL) {
  m <- coxph(as.formula(paste("Surv(time, event) ~", e, "+", M2)), data = data)
  bh <- basehaz(m, centered = TRUE); H <- approx(bh$time, bh$hazard, xout = 5, rule = 2)$y
  tgt <- if (is.null(sub)) data else data[eval(sub)]
  sapply(levels(data[[e]]), function(l) { nd <- copy(tgt); set(nd, j = e, value = factor(l, levels = levels(data[[e]])))
    mean(1 - exp(-H * exp(predict(m, newdata = nd, type = "lp")))) })
}
std_all <- function(data) c(std_risk(data, "hard_cat"), std_risk(data, "loose_cat"),
                            std_risk(data, "hard_cat", quote(age_q >= 65)), std_risk(data, "loose_cat", quote(age_q >= 65)))
est <- std_all(slim[!is.na(hard_cat) & !is.na(loose_cat)])
labs <- c(paste("hard", levels(slim$hard_cat)), paste("loose", levels(slim$loose_cat)),
          paste("age65 hard", levels(slim$hard_cat)), paste("age65 loose", levels(slim$loose_cat)))
base_dat <- slim[!is.na(hard_cat) & !is.na(loose_cat)]
bt <- do.call(rbind, mclapply(seq_len(B), function(b) tryCatch(std_all(base_dat[sample.int(nrow(base_dat), replace = TRUE)]),
                                                                error = function(e) rep(NA, length(est))), mc.cores = 4))
rd <- function(v, ref) 1000 * (v - ref)        # risk difference per 1,000 over 5 years
ref_idx <- c(1, 1, 1, 4, 4, 4, 7, 7, 7, 10, 10, 10)
abs_tab <- data.table(group = labs, risk5_pct = 100 * est,
                      risk5_lo = 100 * apply(bt, 2, quantile, 0.025, na.rm = TRUE), risk5_hi = 100 * apply(bt, 2, quantile, 0.975, na.rm = TRUE),
                      rd_per1000 = rd(est, est[ref_idx]),
                      rd_lo = sapply(seq_along(est), function(i) quantile(rd(bt[, i], bt[, ref_idx[i]]), 0.025, na.rm = TRUE)),
                      rd_hi = sapply(seq_along(est), function(i) quantile(rd(bt[, i], bt[, ref_idx[i]]), 0.975, na.rm = TRUE)))
fwrite(aj, paste0(R_, "tier1_cuminc_aj.csv")); fwrite(abs_tab, paste0(R_, "tier1_absolute_risk.csv"))

options(width = 200)
print(neg_tab); print(cnt_neg); print(def_tab[grepl("Often|Sometimes", term)]); print(def_info)
print(cat_tab); print(cnt_score); cat("p nonlinearity", p_nonlin, "\n")
print(sub_tab[grepl("Often", term), .(subgroup, exposure, lvl, n, events, HR = round(HR, 2), lo = round(lo, 2), hi = round(hi, 2), p_interaction = signif(p_interaction, 2))])
print(aj); print(abs_tab)
