# Phenome-wide Cox analysis: bowel habits (digestive health questionnaire) and incident disease / mortality.
# Usage: Rscript 02_phewas_cox.R [model]   model = 1 (age + sex), 2 (main, fully adjusted; default) or 3 (2 + caffeine)
suppressPackageStartupMessages(library(parallel))
args  <- commandArgs(trailingOnly = TRUE)
MODEL <- if (length(args)) as.integer(args[1]) else 2L
MIN_EVENTS <- 100
source("/path/to/project/scripts/prep.R")
COV <- switch(MODEL, M1, M2, M3)

# ---------------------------------------------------------------- one outcome
fit_outcome <- function(code) {
  tt <- outcome_tt(oc[[code]], startsWith(code, "DEATH_"))
  if (sum(tt$event[tt$sel]) < MIN_EVENTS) return(NULL)
  dd <- co[tt$sel]; dd[, `:=`(time = tt$time[tt$sel], event = tt$event[tt$sel])]
  out <- list()
  for (ex in EXPOSURES) {
    f <- as.formula(paste("Surv(time, event) ~", ex, "+", COV))
    m <- tryCatch(coxph(f, data = dd, ties = "efron"), error = function(e) NULL)
    if (is.null(m)) next
    s <- summary(m)$coefficients; ci <- summary(m)$conf.int
    rows <- grep(paste0("^", ex), rownames(s))
    ev_by <- dd[, .(events = sum(event), n = .N), by = ex]
    for (r in rows) {
      lvl <- sub(ex, "", rownames(s)[r], fixed = TRUE)
      out[[length(out) + 1]] <- data.table(code = code, exposure = ex, level = lvl,
        HR = ci[r, "exp(coef)"], lo = ci[r, "lower .95"], hi = ci[r, "upper .95"], p = s[r, "Pr(>|z|)"],
        events_level = ev_by[get(ex) == lvl, events], n_level = ev_by[get(ex) == lvl, n],
        events_total = sum(dd$event), n_total = nrow(dd))
    }
  }
  rbindlist(out)
}

codes <- setdiff(names(oc), "eid")
t0 <- Sys.time()
res <- rbindlist(mclapply(codes, fit_outcome, mc.cores = 11), fill = TRUE)
res <- merge(res, dic, by = "code", all.x = TRUE)
res[, fdr := p.adjust(p, method = "BH"), by = exposure]
setorder(res, exposure, p)
fwrite(res, paste0(P, sprintf("results/phewas_model%d.csv", MODEL)))
cat(sprintf("Model %d: %d outcomes analysed, %d rows, %.1f min\n", MODEL, uniqueN(res$code), nrow(res),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("FDR<0.05 by exposure level:\n"); print(res[fdr < 0.05, .N, by = .(exposure, level)])
