# Time-origin comparison for Parkinson's disease. Bowel habits were reported on the Digestive Health Questionnaire
# (2017-2020), so follow-up should start at questionnaire completion (main analysis). A previous UK Biobank analysis
# (Zhang et al., Neuroepidemiology 2025) started follow-up at recruitment (2006-2010) and excluded only PD present at
# recruitment, so PD diagnosed between recruitment and the questionnaire was counted as incident although the bowel
# habits were reported after diagnosis. Here both designs are fitted in the same analytic sample with Model 2.
source("/path/to/project/scripts/prep.R")
bd <- fread(paste0(P, "data/baseline_date.csv"))
co[, baseline_date := as.IDate(bd$p53_i0[match(co$eid, bd$eid)])]
co[, offset := as.numeric(as.IDate(q_date) - baseline_date)]  # days from recruitment to questionnaire
d <- oc[["G20"]]                                              # days from questionnaire; -99999 = prevalent, date unknown
fu <- co$fu_end_days

fit <- function(dd, label) {
  rbindlist(lapply(c("hard_cat", "loose_cat"), function(ex) {
    m <- coxph(as.formula(paste("Surv(time, event) ~", ex, "+", M2)), data = dd)
    s <- summary(m); r <- grep(paste0("^", ex), rownames(s$coefficients))
    data.table(design = label, exposure = ex, level = sub(ex, "", rownames(s$coefficients)[r], fixed = TRUE),
               HR = s$conf.int[r, 1], lo = s$conf.int[r, 3], hi = s$conf.int[r, 4], p = s$coefficients[r, 5],
               events = m$nevent, n = m$n)
  }))
}

# A. main design: origin = questionnaire, PD before the questionnaire excluded
a <- copy(co); a[, `:=`(event = as.integer(!is.na(d) & d > 0 & d <= fu), keep = !(!is.na(d) & d <= 0))]
a[, time := ifelse(event == 1, d, fu) / 365.25]
A <- fit(a[keep & time > 0], "A: origin at questionnaire (main analysis)")

# B. recruitment-origin design: origin = recruitment, only PD before recruitment (or with unknown date) excluded
b <- copy(co)
b[, pre_recruit := !is.na(d) & (d == -99999 | d <= -offset)]
b[, event := as.integer(!is.na(d) & !pre_recruit & d <= fu)]
b[, time := ifelse(event == 1, d + offset, fu + offset) / 365.25]
b[, before_q := event == 1 & d <= 0]                           # "incident" cases diagnosed before bowel habits were reported
B <- fit(b[!pre_recruit & !is.na(offset) & time > 0], "B: origin at recruitment (as in Zhang et al.)")
# B8: recruitment origin, excluding cases in the first 8 years (their longest lag)
b8 <- b[!pre_recruit & !is.na(offset) & time > 0 & !(event == 1 & time < 8)]
B8 <- fit(b8, "B8: origin at recruitment, events in first 8 years excluded")
# C. only the cases diagnosed between recruitment and the questionnaire (follow-up censored at the questionnaire)
# (computed on the full table so that the outcome vector d stays aligned with the rows)
cc <- copy(b)
cc[, `:=`(event = as.integer(before_q), time = ifelse(before_q, d + offset, offset) / 365.25)]
C <- fit(cc[!pre_recruit & !is.na(offset) & time > 0], "C: cases diagnosed between recruitment and questionnaire only")

res <- rbind(A, B, B8, C)
fwrite(res, paste0(P, "results/pd_time_origin_comparison.csv"))
cnt <- data.table(
  n_main = nrow(a[keep & time > 0]), events_main = sum(a[keep & time > 0]$event),
  n_recruit = nrow(b[!pre_recruit & !is.na(offset) & time > 0]),
  events_recruit = sum(b[!pre_recruit & !is.na(offset) & time > 0]$event),
  events_before_questionnaire = sum(b$before_q), events_before_q_after_8y = sum(b8$event == 1 & b8$before_q),
  median_fu_recruit = median(b[!pre_recruit & !is.na(offset) & time > 0]$time),
  median_offset_y = median(co$offset, na.rm = TRUE) / 365.25)
fwrite(cnt, paste0(P, "results/pd_time_origin_counts.csv"))
print(cnt)
print(res[level == "Often+" | level == "Sometimes", .(design, exposure, level, HR = round(HR, 2), lo = round(lo, 2),
                                                      hi = round(hi, 2), p = signif(p, 2), events, n)])
