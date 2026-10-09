# Dose-response: five-level hard/loose stool frequency and six-level bowel movement frequency, fully adjusted (Model 2).
suppressPackageStartupMessages(library(parallel))
source("/path/to/project/scripts/prep.R")

FREQ5 <- c("Never", "Sometimes", "Often", "Most of the time", "Always")
co[, hard5  := factor(FREQ5[hard_stool + 1],  levels = FREQ5)]
co[, loose5 := factor(FREQ5[loose_stool + 1], levels = FREQ5)]
co[, bmf6 := factor(fcase(bmf_day == 0 & bmf_least_week <= 2, "<=2/week",
                          bmf_day == 0, "3-6/week",
                          bmf_day == 1, "1/day", bmf_day == 2, "2/day", bmf_day == 3, "3/day", bmf_day >= 4, ">=4/day"),
                    levels = c("1/day", "<=2/week", "3-6/week", "2/day", "3/day", ">=4/day"))]
co[, bmf6_score := fcase(bmf6 == "<=2/week", 1, bmf6 == "3-6/week", 2, bmf6 == "1/day", 3, bmf6 == "2/day", 4,
                         bmf6 == "3/day", 5, bmf6 == ">=4/day", 6)]

PLAN <- list(
  list(ex = "hard5",  trend = "hard_stool",  outcomes = c("Parkinson's disease", "Dementia", "Depression", "Anxiety")),
  list(ex = "loose5", trend = "loose_stool", outcomes = c("Parkinson's disease", "Dementia", "Depression", "Anxiety", "Type 2 diabetes")),
  list(ex = "bmf6",   trend = NA,            outcomes = c("Type 2 diabetes", "Heart failure", "Hypertension", "Stroke",
                                                          "Ischaemic heart disease", "Atrial fibrillation", "Depression",
                                                          "All-cause mortality", "Parkinson's disease")))
KEY <- list("All-cause mortality" = "DEATH_all", "Ischaemic heart disease" = c("I20", "I21", "I22", "I23", "I24", "I25"),
            "Heart failure" = "I50", "Atrial fibrillation" = "I48", "Stroke" = c("I60", "I61", "I63", "I64"),
            "Hypertension" = "I10", "Type 2 diabetes" = "E11", "Depression" = c("F32", "F33"), "Anxiety" = c("F40", "F41"),
            "Dementia" = c("F00", "F01", "F02", "F03", "G30"), "Parkinson's disease" = "G20")

combine_days <- function(codes) {
  m <- as.matrix(oc[, ..codes])
  prev <- rowSums(!is.na(m) & m <= 0) > 0
  d <- suppressWarnings(apply(m, 1, function(r) if (all(is.na(r))) NA_real_ else min(r, na.rm = TRUE)))
  d[prev] <- -1
  d
}

jobs <- rbindlist(lapply(PLAN, function(p) data.table(ex = p$ex, trend = p$trend, outcome = p$outcomes)))
run <- function(i) {
  j <- jobs[i]; codes <- KEY[[j$outcome]]
  d <- if (length(codes) == 1) oc[[codes]] else combine_days(codes)
  tt <- outcome_tt(d, startsWith(codes[1], "DEATH_"))
  dd <- co[tt$sel]; dd[, `:=`(time = tt$time[tt$sel], event = tt$event[tt$sel])]
  m <- coxph(as.formula(paste("Surv(time, event) ~", j$ex, "+", M2)), data = dd)
  s <- summary(m); rows <- grep(paste0("^", j$ex), rownames(s$coefficients))
  ev <- dd[, .(events = sum(event), n = .N), by = c(j$ex)]; setnames(ev, j$ex, "level")
  ref <- levels(dd[[j$ex]])[1]
  out <- rbind(data.table(level = ref, HR = 1, lo = 1, hi = 1, p = NA_real_),
               data.table(level = sub(j$ex, "", rownames(s$coefficients)[rows], fixed = TRUE),
                          HR = s$conf.int[rows, "exp(coef)"], lo = s$conf.int[rows, "lower .95"],
                          hi = s$conf.int[rows, "upper .95"], p = s$coefficients[rows, "Pr(>|z|)"]))
  out <- merge(out, ev, by = "level", all.x = TRUE)
  # trend: per one-level increase (ordinal score); for bmf6 trends are fitted separately below and above 1/day
  if (!is.na(j$trend)) {
    mt <- coxph(as.formula(paste("Surv(time, event) ~", j$trend, "+", M2)), data = dd)
    out[, p_trend := summary(mt)$coefficients[j$trend, "Pr(>|z|)"]]
    out[, HR_per_level := exp(coef(mt)[j$trend])]
  } else {
    lo_side <- coxph(as.formula(paste("Surv(time, event) ~ bmf6_score +", M2)), data = dd[bmf6_score <= 3])
    hi_side <- coxph(as.formula(paste("Surv(time, event) ~ bmf6_score +", M2)), data = dd[bmf6_score >= 3])
    out[, `:=`(p_trend_below_1day = summary(lo_side)$coefficients["bmf6_score", "Pr(>|z|)"],
               p_trend_above_1day = summary(hi_side)$coefficients["bmf6_score", "Pr(>|z|)"])]
  }
  out[, `:=`(exposure = j$ex, outcome = j$outcome, level = factor(level, levels = levels(dd[[j$ex]])))]
  setorder(out, level)
  out
}
t0 <- Sys.time()
res <- rbindlist(mclapply(seq_len(nrow(jobs)), run, mc.cores = 3, mc.preschedule = FALSE), fill = TRUE)
setcolorder(res, c("exposure", "outcome", "level"))
fwrite(res, paste0(P, "results/dose_response.csv"))
cat(sprintf("Done in %.1f min, %d rows\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(res)))
