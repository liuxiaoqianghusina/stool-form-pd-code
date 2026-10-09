# Tier-2 analyses for Parkinson's disease (PD), Model 2 covariates unless stated:
#  (6) stool form jointly with other (non-gastrointestinal) prodromal markers recorded before the questionnaire
#      (hypotension, urinary dysfunction, sleep disorders, depression/anxiety): joint HRs, Aalen-Johansen 5-year
#      cumulative incidence, and additive (RERI) and multiplicative interaction;
#  (7) accelerometry sub-study (2013-2015): stool form and average acceleration as joint predictors of PD;
#  (8) medication sensitivity: excluding baseline users of drugs that can cause parkinsonism, anticholinergics or
#      dopaminergic drugs;
#  (9) multiple imputation of missing covariates (mice, m = 5, PMM / logistic / polytomous models with the event
#      indicator and Nelson-Aalen cumulative hazard), pooled with Rubin's rules; plus complete-case analysis.
suppressPackageStartupMessages({library(parallel); library(mice)})
source("/path/to/project/scripts/prep.R")
set.seed(20261010)
R_ <- paste0(P, "results/")
hr_rows <- function(m, pattern, label) {
  s <- summary(m); r <- grep(pattern, rownames(s$coefficients))
  data.table(analysis = label, term = rownames(s$coefficients)[r], HR = s$conf.int[r, 1], lo = s$conf.int[r, 3],
             hi = s$conf.int[r, 4], p = s$coefficients[r, 5], events = m$nevent, n = m$n)
}
prev <- function(codes) { m <- as.matrix(oc[, ..codes]); as.integer(rowSums(!is.na(m) & (m <= 0)) > 0) }
co[, `:=`(hypotension = prev("I95"), urinary = prev(c("N31", "N32", "N39")), sleep = prev("G47"),
          mood = as.integer(prev_depression == 1 | prev_anxiety == 1))]
co[, n_prodromal := hypotension + urinary + sleep + mood]
tt <- outcome_tt(oc[["G20"]], FALSE); keep_rows <- which(tt$sel)
dd <- co[keep_rows]; dd[, `:=`(time = tt$time[keep_rows], event = tt$event[keep_rows])]
dth <- oc$DEATH_all[keep_rows]
dd[, status3 := factor(fifelse(event == 1, 1L, fifelse(!is.na(dth) & dth > 0 & dth / 365.25 <= time + 1e-9, 2L, 0L)),
                       levels = 0:2, labels = c("censor", "PD", "death"))]

# ---------------------------------------------------------------- (6) joint with other prodromal markers
dd[, form3 := factor(fcase(hard_stool >= 2 & loose_stool < 2, "Hard often+",
                           loose_stool >= 2 & hard_stool < 2, "Loose often+",
                           hard_stool < 2 & loose_stool < 2, "Neither", default = NA_character_),
                     levels = c("Neither", "Hard often+", "Loose often+"))]
dd[, prod2 := factor(fifelse(n_prodromal >= 1, "Any other prodromal marker", "No other prodromal marker"),
                     levels = c("No other prodromal marker", "Any other prodromal marker"))]
dd[, joint := factor(paste(form3, prod2, sep = " | "))]
dd[, joint := relevel(joint, ref = "Neither | No other prodromal marker")]
M2p <- gsub(" \\+ prev_depression| \\+ prev_anxiety", "", M2)       # mood is part of the marker count
jt <- dd[!is.na(form3)]
m_joint <- coxph(as.formula(paste("Surv(time, event) ~ joint +", M2p)), data = jt)
joint_tab <- rbind(data.table(analysis = "joint", term = "jointNeither | No other prodromal marker", HR = 1, lo = 1, hi = 1, p = NA,
                              events = m_joint$nevent, n = m_joint$n), hr_rows(m_joint, "^joint", "joint"))
cnt_joint <- jt[, .(n = .N, events = sum(event)), by = joint]
f <- survfit(Surv(time, status3) ~ joint, data = jt); s <- summary(f, times = 5, extend = TRUE); k <- which(f$states == "PD")
aj_joint <- data.table(joint = sub("joint=", "", s$strata), cuminc5 = 100 * s$pstate[, k], lo = 100 * s$lower[, k], hi = 100 * s$upper[, k])
# stool form within strata of other prodromal markers; multiplicative interaction; RERI for hard often+ x any marker
within <- rbindlist(lapply(levels(dd$prod2), function(l) rbind(
  hr_rows(coxph(as.formula(paste("Surv(time, event) ~ hard_cat +", M2p)), data = dd[prod2 == l]), "^hard_cat", l),
  hr_rows(coxph(as.formula(paste("Surv(time, event) ~ loose_cat +", M2p)), data = dd[prod2 == l]), "^loose_cat", l))))
p_mult <- sapply(c("hard_cat", "loose_cat"), function(e) {
  m0 <- coxph(as.formula(paste("Surv(time, event) ~ prod2 +", e, "+", M2p)), data = dd)
  m1 <- coxph(as.formula(paste("Surv(time, event) ~ prod2 *", e, "+", M2p)), data = dd)
  anova(m0, m1)$`Pr(>|Chi|)`[2] })
dd[, hard2 := as.integer(hard_stool >= 2)]; dd[, prodA := as.integer(n_prodromal >= 1)]
m_reri <- coxph(as.formula(paste("Surv(time, event) ~ hard2 * prodA +", M2p)), data = dd[!is.na(hard2)])
b <- coef(m_reri)[c("hard2", "prodA", "hard2:prodA")]; V <- vcov(m_reri)[names(b), names(b)]
reri <- exp(sum(b)) - exp(b[1]) - exp(b[2]) + 1
g <- c(exp(sum(b)) - exp(b[1]), exp(sum(b)) - exp(b[2]), exp(sum(b)))      # delta-method gradient
se_reri <- sqrt(drop(t(g) %*% V %*% g))
inter <- data.table(p_mult_hard = p_mult[["hard_cat"]], p_mult_loose = p_mult[["loose_cat"]], RERI = reri,
                    RERI_lo = reri - 1.96 * se_reri, RERI_hi = reri + 1.96 * se_reri,
                    n_any_marker = sum(dd$n_prodromal >= 1), pct_any_marker = 100 * mean(dd$n_prodromal >= 1))
fwrite(joint_tab, paste0(R_, "tier2_joint.csv")); fwrite(cnt_joint, paste0(R_, "tier2_joint_counts.csv"))
fwrite(aj_joint, paste0(R_, "tier2_joint_cuminc.csv")); fwrite(within, paste0(R_, "tier2_within_prodromal.csv"))
fwrite(inter, paste0(R_, "tier2_interaction.csv"))

# ---------------------------------------------------------------- (7) accelerometry
acc <- fread(paste0(P, "data/accelerometry.csv"), colClasses = "character")
ma <- match(dd$eid, as.integer(acc$eid))
dd[, `:=`(acc_mg = suppressWarnings(as.numeric(acc$p90012[ma])), acc_wear = acc$p90015[ma], acc_cal = acc$p90016[ma],
          acc_date = as.IDate(substr(acc$p90010[ma], 1, 10)))]
ac <- dd[!is.na(acc_mg) & acc_wear == "1" & acc_cal == "1" & acc_mg < 100 & !is.na(hard_cat) & !is.na(loose_cat)]
ac[, acc_sd := (acc_mg - mean(acc_mg)) / sd(acc_mg)]
ac[, acc_low := as.integer(acc_mg <= quantile(acc_mg, 0.2))]
acc_rows <- list()
m_a <- coxph(as.formula(paste("Surv(time, event) ~ acc_sd +", M2)), data = ac)
m_f <- coxph(as.formula(paste("Surv(time, event) ~ hard_cat + loose_cat +", M2)), data = ac)
m_af <- coxph(as.formula(paste("Surv(time, event) ~ acc_sd + hard_cat + loose_cat +", M2)), data = ac)
acc_rows <- rbind(hr_rows(m_a, "^acc_sd", "acceleration only"), hr_rows(m_f, "^hard_cat|^loose_cat", "stool form only"),
                  hr_rows(m_af, "^acc_sd|^hard_cat|^loose_cat", "mutually adjusted"))
# discrimination in the sub-study: age + sex base, + acceleration, + stool form, + both
BASEa <- "age_q + I(age_q^2) + sex"
cfit <- function(f) coxph(as.formula(paste("Surv(time, event) ~", f)), data = ac)
mods <- list(base = cfit(BASEa), acc = cfit(paste(BASEa, "+ acc_sd")), form = cfit(paste(BASEa, "+ hard_cat + loose_cat")),
             both = cfit(paste(BASEa, "+ acc_sd + hard_cat + loose_cat")))
cc <- concordance(mods$base, mods$acc, mods$form, mods$both)
dC <- function(i, j) { v <- numeric(4); v[i] <- -1; v[j] <- 1; d <- sum(v * cc$concordance); se <- sqrt(drop(t(v) %*% cc$var %*% v))
  c(dC = d, lo = d - 1.96 * se, hi = d + 1.96 * se) }
acc_disc <- data.table(model = c("age+sex", "+acceleration", "+stool form", "+both"), C = cc$concordance,
                       rbind(c(NA, NA, NA), dC(1, 2), dC(1, 3), dC(2, 4)))
setnames(acc_disc, c("model", "C", "dC_vs_previous", "lo", "hi"))
acc_disc[, comparison := c("", "vs age+sex", "vs age+sex", "+both vs +acceleration")]
acc_joint <- ac[, .(n = .N, events = sum(event)), by = .(acc_low, hard2 = hard_stool >= 2)]
m_aj <- coxph(as.formula(paste("Surv(time, event) ~ interaction(acc_low, hard_stool >= 2) +", M2)), data = ac)
acc_info <- data.table(n = nrow(ac), events = sum(ac$event), median_gap_y = median(as.numeric(as.IDate(ac$q_date) - ac$acc_date), na.rm = TRUE) / 365.25,
                       mean_mg = mean(ac$acc_mg), sd_mg = sd(ac$acc_mg))
fwrite(acc_rows, paste0(R_, "tier2_accel_hr.csv")); fwrite(acc_disc, paste0(R_, "tier2_accel_discrimination.csv"))
fwrite(acc_info, paste0(R_, "tier2_accel_info.csv"))
fwrite(rbind(acc_joint[, .(term = paste0("acc_low=", acc_low, ", hard_often=", hard2), n, events)]), paste0(R_, "tier2_accel_joint_counts.csv"))
fwrite(hr_rows(m_aj, "^interaction", "joint low acceleration x hard often+"), paste0(R_, "tier2_accel_joint.csv"))

# ---------------------------------------------------------------- (8) medication sensitivity
mf <- fread(paste0(P, "data/meds_flags.csv"))
mm <- match(dd$eid, mf$eid)
dd[, `:=`(dip_drug = mf$dip_drug[mm] %in% 1, anticholinergic = mf$anticholinergic[mm] %in% 1, dopaminergic = mf$dopaminergic[mm] %in% 1)]
med_rows <- rbindlist(lapply(list(
  list("Excluding drugs that can cause parkinsonism", quote(dip_drug == FALSE)),
  list("Excluding anticholinergic drugs", quote(anticholinergic == FALSE)),
  list("Excluding dopaminergic drugs", quote(dopaminergic == FALSE)),
  list("Excluding all three drug groups", quote(dip_drug == FALSE & anticholinergic == FALSE & dopaminergic == FALSE))), function(a)
  rbind(hr_rows(coxph(as.formula(paste("Surv(time, event) ~ hard_cat +", M2)), data = dd[eval(a[[2]])]), "^hard_cat", a[[1]]),
        hr_rows(coxph(as.formula(paste("Surv(time, event) ~ loose_cat +", M2)), data = dd[eval(a[[2]])]), "^loose_cat", a[[1]]))))
med_cnt <- dd[, .(dip_drug = sum(dip_drug), anticholinergic = sum(anticholinergic), dopaminergic = sum(dopaminergic),
                  pd_events_any_drug = sum(event[dip_drug | anticholinergic | dopaminergic]))]
fwrite(med_rows, paste0(R_, "tier2_medication.csv")); fwrite(med_cnt, paste0(R_, "tier2_medication_counts.csv"))

# ---------------------------------------------------------------- (9) multiple imputation and complete-case analysis
raw <- fread(paste0(P, "data/cohort.csv"))[match(dd$eid, eid)]
NUM <- c("tdi", "met_min", "bmi", "fruit", "veg", "water_glasses", "n_medications")
CAT <- c("ethnicity", "degree", "smoking", "alcohol_freq", "wholemeal_bread", "fibre_cereal")
FIX <- c("age_q", "sex", "laxative_any", "metformin", "ppi", "opioid", "antidepressant", "prev_diabetes", "prev_hypertension",
         "prev_ihd", "prev_stroke", "prev_depression", "prev_anxiety", "prev_ibs", "prev_ckd", "prev_cancer")
mi <- data.table(dd[, ..FIX], raw[, ..NUM], raw[, ..CAT], hard_cat = dd$hard_cat, loose_cat = dd$loose_cat,
                 time = dd$time, event = dd$event)
mi[, met_min := log1p(met_min)]
for (v in c("ethnicity", "smoking", "alcohol_freq")) set(mi, j = v, value = factor(mi[[v]]))
for (v in c("degree", "wholemeal_bread", "fibre_cereal")) set(mi, j = v, value = factor(mi[[v]]))
mi[, nelson := nelsonaalen(as.data.frame(mi), time, event)]
cc_n <- sum(complete.cases(mi))
meth <- make.method(mi); pm <- make.predictorMatrix(mi); pm[, "time"] <- 0
imp <- mice(mi, m = 5, maxit = 5, method = meth, predictorMatrix = pm, printFlag = FALSE, seed = 20261010)
COVMI <- paste(c(FIX[-(1:2)], NUM, CAT), collapse = " + ")
pool_rows <- rbindlist(lapply(c("hard_cat", "loose_cat"), function(e) {
  fit <- with(imp, coxph(as.formula(paste("Surv(time, event) ~", e, "+ age_q + sex +", COVMI))))
  sm <- summary(pool(fit), conf.int = TRUE); sm <- sm[grepl(paste0("^", e), sm$term), ]
  data.table(analysis = "Multiple imputation (m = 5)", term = as.character(sm$term), HR = exp(sm$estimate),
             lo = exp(sm$`2.5 %`), hi = exp(sm$`97.5 %`), p = sm$p.value, events = sum(mi$event), n = nrow(mi))
}))
cc_rows <- rbindlist(lapply(c("hard_cat", "loose_cat"), function(e)
  hr_rows(coxph(as.formula(paste("Surv(time, event) ~", e, "+ age_q + sex +", COVMI)), data = mi[complete.cases(mi)]),
          paste0("^", e), "Complete-case analysis")))
fwrite(rbind(pool_rows, cc_rows), paste0(R_, "tier2_missing_data.csv"))
fwrite(data.table(n = nrow(mi), complete_cases = cc_n, pct_complete = 100 * cc_n / nrow(mi)), paste0(R_, "tier2_missing_info.csv"))

options(width = 220)
print(joint_tab); print(cnt_joint); print(aj_joint); print(within[grepl("Often", term)]); print(inter)
print(acc_info); print(acc_rows); print(acc_disc); print(med_cnt); print(med_rows[grepl("Often", term)])
print(rbind(pool_rows, cc_rows)); cat("complete cases:", cc_n, "\n")
