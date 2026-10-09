# Key-outcome analyses: model ladder, sensitivity analyses, competing risk of death, time-stratified hazard ratios
# (prodromal/reverse-causation check) and E-values.
# Usage: Rscript 03_key_outcomes.R [main|finegray]  (Fine-Gray models are memory-heavy, so they run as a separate step)
suppressPackageStartupMessages(library(parallel))
STEP <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else "main"
CORES <- if (STEP == "main") 6L else 3L
source("/path/to/project/scripts/prep.R")

KEY <- list(
  "All-cause mortality" = "DEATH_all", "CVD mortality" = "DEATH_CVD",
  "Ischaemic heart disease" = c("I20", "I21", "I22", "I23", "I24", "I25"), "Heart failure" = "I50",
  "Atrial fibrillation" = "I48", "Stroke" = c("I60", "I61", "I63", "I64"), "Hypertension" = "I10",
  "Type 2 diabetes" = "E11", "Chronic kidney disease" = "N18", "Depression" = c("F32", "F33"),
  "Anxiety" = c("F40", "F41"), "Dementia" = c("F00", "F01", "F02", "F03", "G30"), "Parkinson's disease" = "G20",
  "Colorectal cancer" = "CR_C18-C20")
KEY_EXPOSURES <- c("bmf_cat", "hard_cat", "loose_cat")

combine_days <- function(codes) {
  m <- as.matrix(oc[, ..codes])
  prev <- rowSums(!is.na(m) & m <= 0) > 0
  d <- suppressWarnings(apply(m, 1, function(r) if (all(is.na(r))) NA_real_ else min(r, na.rm = TRUE)))
  d[prev] <- -1
  d
}

evalue <- function(hr, lo, hi) {
  rr <- ifelse(hr < 1, 1 / hr, hr)
  e  <- rr + sqrt(rr * (rr - 1))
  lim <- ifelse(hr < 1, ifelse(hi >= 1, 1, 1 / hi), ifelse(lo <= 1, 1, lo))
  c(e, lim + sqrt(lim * (lim - 1)))
}

tidy <- function(m, ex, analysis) {
  s <- summary(m)$coefficients; ci <- summary(m)$conf.int
  rows <- grep(paste0("^", ex), rownames(s))
  rows <- rows[!grepl(":", rownames(s)[rows])]
  data.table(exposure = ex, level = sub(ex, "", rownames(s)[rows], fixed = TRUE), analysis = analysis,
             HR = ci[rows, "exp(coef)"], lo = ci[rows, "lower .95"], hi = ci[rows, "upper .95"], p = s[rows, "Pr(>|z|)"],
             events = m$nevent, n = m$n)
}

run_key <- function(name) {
  codes <- KEY[[name]]
  is_death <- startsWith(codes[1], "DEATH_")
  d <- if (length(codes) == 1) oc[[codes]] else combine_days(codes)
  tt <- outcome_tt(d, is_death)
  dd <- co[tt$sel]; dd[, `:=`(time = tt$time[tt$sel], event = tt$event[tt$sel])]
  # death from other causes, for the competing-risk model
  dth <- oc$DEATH_all[tt$sel]
  dd[, status3 := fifelse(event == 1, 1L, fifelse(!is.na(dth) & dth > 0 & dth / 365.25 <= time + 1e-9, 2L, 0L))]
  out <- list()
  for (ex in KEY_EXPOSURES) {
    fit <- function(cov, data, analysis) {
      m <- tryCatch(coxph(as.formula(paste("Surv(time, event) ~", ex, "+", cov)), data = data), error = function(e) NULL)
      if (!is.null(m)) out[[length(out) + 1]] <<- tidy(m, ex, analysis)
    }
    if (STEP == "main") {
      fit(M1, dd, "Model 1: age + sex")
      fit(M2, dd, "Model 2: fully adjusted (main)")
      fit(M3, dd, "Model 3: + coffee and tea")
      fit(M2, dd[!(event == 1 & time < 2)], "Excluding events in first 2 years")
      fit(M2, dd[!(event == 1 & time < 4)], "Excluding events in first 4 years")
      fit(M2, dd[prev_ibs == 0 & (is.na(ibs_ever_dx_q) | ibs_ever_dx_q == 0)], "Excluding IBS")
      fit(M2, dd[laxative_any == 0 & metformin == 0 & opioid == 0], "Excluding laxative, metformin and opioid users")
      # time-stratified HRs: separate counting-process fits within follow-up windows
      for (w in list(c(0, 2), c(2, 4), c(4, 99))) {
        dw <- dd[time > w[1]]
        dw[, `:=`(t0 = w[1], t1 = pmin(time, w[2]), ev = as.integer(event == 1 & time <= w[2]))]
        m <- tryCatch(coxph(as.formula(paste("Surv(t0, t1, ev) ~", ex, "+", M2)), data = dw), error = function(e) NULL)
        if (!is.null(m) && m$nevent >= 20) out[[length(out) + 1]] <- tidy(m, ex, sprintf("Follow-up window %g-%s y", w[1], ifelse(w[2] == 99, "end", w[2])))
      }
    }
    if (STEP == "finegray" && !is_death) {
      fg <- tryCatch({
        fgd <- finegray(as.formula(paste("Surv(time, factor(status3)) ~", ex, "+", M2)), data = dd, etype = "1")
        coxph(as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~", ex, "+", M2)), data = fgd, weights = fgwt)
      }, error = function(e) NULL)
      if (!is.null(fg)) out[[length(out) + 1]] <- tidy(fg, ex, "Fine-Gray (competing risk of death)")
    }
  }
  res <- rbindlist(out); res[, outcome := name]
  res
}

t0 <- Sys.time()
res <- rbindlist(mclapply(names(KEY), run_key, mc.cores = CORES, mc.preschedule = FALSE), fill = TRUE)
if (STEP == "finegray") {
  fwrite(res, paste0(P, "results/key_outcomes_finegray.csv"))
  quit(save = "no")
}
main <- res[analysis == "Model 2: fully adjusted (main)"]
ev <- t(mapply(evalue, main$HR, main$lo, main$hi))
main[, `:=`(evalue_point = ev[, 1], evalue_ci = ev[, 2])]
res <- merge(res, main[, .(outcome, exposure, level, evalue_point, evalue_ci)], by = c("outcome", "exposure", "level"), all.x = TRUE)
setcolorder(res, c("outcome", "exposure", "level", "analysis"))
fwrite(res, paste0(P, "results/key_outcomes.csv"))
cat(sprintf("Done in %.1f min, %d rows\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(res)))
