# Shared data preparation: analytic sample, exposure categories, covariates and model formulas.
# Sourced by 02_phewas_cox.R and 03_key_outcomes.R so both use identical definitions.
suppressPackageStartupMessages({library(data.table); library(survival)})
P <- "/path/to/project/"

co  <- fread(paste0(P, "data/cohort.csv"))
oc  <- fread(paste0(P, "data/outcomes.csv"))
dic <- fread(paste0(P, "data/outcome_dictionary.csv"))
stopifnot(identical(co$eid, oc$eid))

# ---------------------------------------------------------------- analytic sample and exposures
keep <- !is.na(co$bmf_day) & co$excl_ibd == 0 & co$excl_coeliac == 0 & co$excl_crc == 0 & co$excl_colectomy == 0 &
        co$fu_end_days > 0
co <- co[keep]; oc <- oc[keep]
co[, bmf_cat := factor(fcase(bmf_day == 0, "<1/day", bmf_day == 1, "1/day", bmf_day == 2, "2/day", bmf_day >= 3, ">=3/day"),
                       levels = c("1/day", "<1/day", "2/day", ">=3/day"))]
co[, hard_cat  := factor(fcase(hard_stool == 0, "Never", hard_stool == 1, "Sometimes", hard_stool >= 2, "Often+"),
                         levels = c("Never", "Sometimes", "Often+"))]
co[, loose_cat := factor(fcase(loose_stool == 0, "Never", loose_stool == 1, "Sometimes", loose_stool >= 2, "Often+"),
                         levels = c("Never", "Sometimes", "Often+"))]
co[, constip := (bmf_day == 0 | hard_stool >= 2)]
co[, diarr   := (bmf_day >= 3 | loose_stool >= 2)]
co[, phenotype := factor(fcase(constip & diarr, "Mixed", constip, "Constipation-type", diarr, "Diarrhoea-type",
                               !is.na(constip) & !is.na(diarr), "Normal"),
                         levels = c("Normal", "Constipation-type", "Diarrhoea-type", "Mixed"))]

# ---------------------------------------------------------------- covariates (missing -> indicator / category)
fill_num <- function(x) { m <- is.na(x); x[m] <- median(x, na.rm = TRUE); list(x, as.integer(m)) }
for (v in c("tdi", "met_min", "bmi", "fruit", "veg", "water_glasses", "n_medications", "coffee_cups", "tea_cups")) {
  r <- fill_num(co[[v]]); set(co, j = v, value = r[[1]]); set(co, j = paste0(v, "_miss"), value = r[[2]])
}
co[, met_min := log1p(met_min)]
fac <- function(x) { x <- as.character(x); x[is.na(x)] <- "Missing"; factor(x) }
for (v in c("ethnicity", "degree", "smoking", "alcohol_freq", "wholemeal_bread", "fibre_cereal"))
  set(co, j = v, value = fac(co[[v]]))

M1 <- "age_q + sex"
M2 <- paste(M1, "+ ethnicity + tdi + tdi_miss + degree + smoking + alcohol_freq + met_min + met_min_miss + bmi + bmi_miss +",
            "fruit + fruit_miss + veg + veg_miss + wholemeal_bread + fibre_cereal + water_glasses + water_glasses_miss +",
            "n_medications + n_medications_miss + laxative_any + metformin + ppi + opioid + antidepressant +",
            "prev_diabetes + prev_hypertension + prev_ihd + prev_stroke + prev_depression + prev_anxiety + prev_ibs +",
            "prev_ckd + prev_cancer")
M3 <- paste(M2, "+ coffee_cups + coffee_cups_miss + tea_cups + tea_cups_miss")  # M2 + caffeine sources
EXPOSURES <- c("bmf_cat", "hard_cat", "loose_cat", "phenotype")

# event / time for one outcome column (days from questionnaire); deaths use the death-register follow-up end
outcome_tt <- function(d, is_death) {
  fu <- if (is_death) co$fu_end_death_days else co$fu_end_days
  prevalent <- if (is_death) rep(FALSE, length(d)) else (!is.na(d) & d <= 0)
  event <- as.integer(!is.na(d) & d > 0 & d <= fu)
  time  <- ifelse(event == 1, d, fu) / 365.25
  list(event = event, time = time, sel = !prevalent & time > 0)
}
