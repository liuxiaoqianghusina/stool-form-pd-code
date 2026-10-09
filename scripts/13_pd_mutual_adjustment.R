# Parkinson's disease: (1) hard stools, loose stools and bowel movement frequency fitted jointly (mutual adjustment),
# (2) proportional-hazards check (Schoenfeld residuals) for the main Model 2 fits.
source("/path/to/project/scripts/prep.R")
tt <- outcome_tt(oc[["G20"]], FALSE)
dd <- co[tt$sel]; dd[, `:=`(time = tt$time[tt$sel], event = tt$event[tt$sel])]
dd <- dd[!is.na(hard_cat) & !is.na(loose_cat)]
tidy <- function(m, analysis) {
  s <- summary(m); rows <- grep("^(hard_cat|loose_cat|bmf_cat)", rownames(s$coefficients))
  data.table(analysis = analysis, term = rownames(s$coefficients)[rows], HR = s$conf.int[rows, "exp(coef)"],
             lo = s$conf.int[rows, "lower .95"], hi = s$conf.int[rows, "upper .95"], p = s$coefficients[rows, "Pr(>|z|)"],
             events = m$nevent, n = m$n)
}
res <- list(); zph <- list()
for (ex in c("hard_cat", "loose_cat", "bmf_cat")) {
  m <- coxph(as.formula(paste("Surv(time, event) ~", ex, "+", M2)), data = dd)
  res[[ex]] <- tidy(m, paste("Single exposure (common sample):", ex))
  # Schoenfeld test; the full Model 2 information matrix is singular for a few sparse covariate levels,
  # so fall back to Model 1 (age + sex) when cox.zph cannot invert it
  z <- tryCatch(list(cox.zph(m)$table, "Model 2"), error = function(e)
    list(cox.zph(coxph(as.formula(paste("Surv(time, event) ~", ex, "+", M1)), data = dd))$table, "Model 1"))
  zph[[ex]] <- data.table(model = z[[2]], exposure = ex, term = c(ex, "GLOBAL"), chisq = z[[1]][c(ex, "GLOBAL"), "chisq"],
                          p = z[[1]][c(ex, "GLOBAL"), "p"])
}
m <- coxph(as.formula(paste("Surv(time, event) ~ hard_cat + loose_cat +", M2)), data = dd)
res[["hl"]] <- tidy(m, "Mutually adjusted: hard + loose")
m <- coxph(as.formula(paste("Surv(time, event) ~ hard_cat + loose_cat + bmf_cat +", M2)), data = dd)
res[["hlb"]] <- tidy(m, "Mutually adjusted: hard + loose + frequency")
out <- rbindlist(res)
fwrite(out, paste0(P, "results/pd_mutual_adjustment.csv"))
fwrite(rbindlist(zph), paste0(P, "results/pd_schoenfeld.csv"))
print(out[, .(analysis, term, HR = round(HR, 2), lo = round(lo, 2), hi = round(hi, 2), p = signif(p, 2), events, n)])
print(rbindlist(zph))
