# Sensitivity analysis for the loose-stool finding: exclude participants with recorded conditions or baseline drugs that
# commonly cause loose stools (hyperthyroidism E05, other non-infective gastroenteritis/colitis including microscopic
# colitis K52, chronic pancreatic disease K86, and users of PPIs or metformin), then refit Model 2 for PD.
source("/path/to/project/scripts/prep.R")
prev <- function(codes) { m <- as.matrix(oc[, ..codes]); as.integer(rowSums(!is.na(m) & (m <= 0)) > 0) }
co[, loose_cause := as.integer(prev(c("E05", "K52", "K86")) == 1 | ppi == 1 | metformin == 1)]
tt <- outcome_tt(oc[["G20"]], FALSE)
dd <- co[tt$sel]; dd[, `:=`(time = tt$time[tt$sel], event = tt$event[tt$sel])]
x <- dd[loose_cause == 0]
res <- rbindlist(lapply(c("hard_cat", "loose_cat"), function(e) {
  m <- coxph(as.formula(paste("Surv(time, event) ~", e, "+", M2)), data = x)
  s <- summary(m); r <- grep(paste0("^", e), rownames(s$coefficients))
  data.table(analysis = "Excluding causes of loose stools (E05, K52, K86, PPI, metformin)", term = rownames(s$coefficients)[r],
             HR = s$conf.int[r, 1], lo = s$conf.int[r, 3], hi = s$conf.int[r, 4], p = s$coefficients[r, 5], events = m$nevent, n = m$n)
}))
fwrite(res, paste0(P, "results/sens_loose_causes.csv"))
print(dd[, .(excluded = sum(loose_cause), pct = 100 * mean(loose_cause), pd_events_excluded = sum(event[loose_cause == 1]))])
print(res)
