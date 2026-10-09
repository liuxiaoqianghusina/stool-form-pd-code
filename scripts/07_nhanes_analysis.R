# External validation in NHANES 2005-2010: survey-weighted Cox models for mortality (to 2019) and survey-weighted
# logistic models for cross-sectional depression, diabetes and hypertension, with exposures harmonised to UKB.
suppressPackageStartupMessages({library(data.table); library(survey); library(survival)})
options(survey.lonely.psu = "adjust")
P <- "/path/to/project/"
d <- fread(paste0(P, "data/nhanes/nhanes_bowel.csv"), na.strings = c("", "NA"))

d[, in_analysis := age >= 20 & !is.na(bmf_cat) & eligible_mort == 1 & !is.na(fu_months)]
d[, bmf_cat := factor(bmf_cat, levels = c("1/day", "<1/day", "2/day", ">=3/day"))]
d[, stool_form := factor(stool_form, levels = c("Normal (BSFS 3-5)", "Hard (BSFS 1-2)", "Loose (BSFS 6-7)"))]
d[, bmf6 := factor(fcase(bm_week <= 2, "<=2/week", bm_week < 7, "3-6/week", bm_week < 14, "1/day",
                         bm_week < 21, "2/day", bm_week < 28, "3/day", bm_week >= 28, ">=4/day"),
                   levels = c("1/day", "<=2/week", "3-6/week", "2/day", "3/day", ">=4/day"))]
d[, time_y := fu_months / 12]
fill <- function(x) { m <- is.na(x); x[m] <- median(x, na.rm = TRUE); list(x, as.integer(m)) }
d[, diet_miss := as.integer(is.na(fibre_g))]  # fibre, energy and caffeine come from the same dietary recall
for (v in c("pir", "bmi", "fibre_g", "energy_kcal", "caffeine_mg")) {
  r <- fill(d[[v]]); set(d, j = v, value = r[[1]]); set(d, j = paste0(v, "_miss"), value = r[[2]])
}
d[, alcohol := factor(fcase(is.na(drinks_day), "Missing", drinks_day == 0, "0", drinks_day <= 2, "1-2", default = "3+"))]
for (v in c("race", "education")) set(d, j = v, value = factor(fifelse(is.na(d[[v]]) | d[[v]] == "", "Missing", d[[v]])))
d[, depression_f := factor(fifelse(is.na(depression), "Missing", as.character(depression)))]
M1 <- "age + female"
M2 <- paste(M1, "+ race + education + pir + pir_miss + smoking + alcohol + bmi + bmi_miss + active +",
            "fibre_g + energy_kcal + diet_miss + diabetes + hypertension + cvd_history +",
            "cancer_history + depression_f + metformin + opioid + laxative + antidepressant + ppi")
M3 <- paste(M2, "+ caffeine_mg")

# Rows outside the analytic subpopulation stay in the design (for correct variance) but carry zero weight after subset();
# give them reference covariate values so no model-matrix column varies only within those rows.
covs <- unique(unlist(strsplit(gsub("\\s", "", M3), "\\+")))
ref <- d[in_analysis == TRUE][1]
for (v in covs) set(d, i = which(!d$in_analysis), j = v, value = ref[[v]])
d[in_analysis == FALSE, `:=`(time_y = 1, death = 0L, death_cvd = 0L, death_cancer = 0L)]

des_all <- svydesign(id = ~psu, strata = ~strata, weights = ~weight, nest = TRUE, data = d)
des <- subset(des_all, in_analysis)

tidy <- function(m, ex, outcome, model) {
  s <- summary(m)$coefficients; ci <- confint(m)
  rows <- grep(paste0("^", ex), rownames(s))
  data.table(outcome = outcome, exposure = ex, level = sub(ex, "", rownames(s)[rows], fixed = TRUE), model = model,
             estimate = exp(s[rows, 1]), lo = exp(ci[rows, 1]), hi = exp(ci[rows, 2]), p = s[rows, ncol(s)])
}
# drop covariate terms that are aliased (constant or collinear) in the analytic rows, using an unweighted fit
dealias <- function(lhs, ex, cov, data) {
  terms <- strsplit(gsub("\\s", "", cov), "\\+")[[1]]
  repeat {
    mm <- model.matrix(as.formula(paste("~", ex, "+", paste(terms, collapse = "+"))), data)
    q <- qr(mm)
    if (q$rank == ncol(mm)) break
    bad <- colnames(mm)[q$pivot[(q$rank + 1):ncol(mm)]]
    drop <- terms[sapply(terms, function(t) any(startsWith(bad, t)))]
    terms <- setdiff(terms, drop[length(drop)])
  }
  as.formula(paste(lhs, "~", ex, "+", paste(terms, collapse = " + ")))
}
an <- d[in_analysis == TRUE]
res <- list()
# ---------------------------------------------------------------- mortality (survey-weighted Cox)
for (ev in c("death", "death_cvd", "death_cancer")) {
  for (ex in c("bmf_cat", "stool_form", "bmf6")) {
    for (mod in list(c("Model 1", M1), c("Model 2", M2), c("Model 3 (+ caffeine)", M3))) {
      f <- dealias(sprintf("Surv(time_y, %s)", ev), ex, mod[2], an[!is.na(get(ex))])
      m <- svycoxph(f, design = subset(des, !is.na(get(ex))))
      res[[length(res) + 1]] <- tidy(m, ex, ev, mod[1])
    }
    # exclude deaths in the first 2 years
    f <- dealias(sprintf("Surv(time_y, %s)", ev), ex, M2, an[!is.na(get(ex))])
    m <- svycoxph(f, design = subset(des, !is.na(get(ex)) & !(get(ev) == 1 & time_y < 2)))
    res[[length(res) + 1]] <- tidy(m, ex, ev, "Model 2, excl. deaths <2 y")
  }
}
# ---------------------------------------------------------------- cross-sectional outcomes (survey-weighted logistic)
M2_noDep <- sub("depression_f \\+ ", "", M2)
for (y in list(c("depression", sub(" \\+ antidepressant", "", M2_noDep)),
               c("diabetes", sub("diabetes \\+ ", "", sub(" \\+ metformin", "", M2))),
               c("hypertension", sub("hypertension \\+ ", "", M2)))) {
  for (ex in c("bmf_cat", "stool_form")) {
    m <- svyglm(dealias(y[1], ex, y[2], an[!is.na(get(ex)) & !is.na(get(y[1]))]), design = subset(des, !is.na(get(ex)) & !is.na(get(y[1]))),
                family = quasibinomial())
    res[[length(res) + 1]] <- tidy(m, ex, paste0(y[1], " (cross-sectional, OR)"), "Model 2")
  }
}
res <- rbindlist(res)
n_tab <- d[in_analysis == TRUE, .(n = .N, deaths = sum(death), cvd_deaths = sum(death_cvd)), by = bmf_cat]
fwrite(res, paste0(P, "results/nhanes_validation.csv"))
fwrite(n_tab, paste0(P, "results/nhanes_counts.csv"))
print(res[model %in% c("Model 2", "Model 2, excl. deaths <2 y")], digits = 3)
