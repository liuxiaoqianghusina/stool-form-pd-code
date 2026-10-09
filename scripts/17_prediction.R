# Incremental predictive value of stool form for incident Parkinson's disease (PD) over a conventional risk model.
# Conventional model: risk and prodromal markers of the MDS research criteria that are available in UK Biobank
# (age, sex, smoking, coffee, physical activity, family history of PD, diabetes, depression, anxiety, hypotension,
# urinary dysfunction, sleep disorders) plus the existing constipation measures (diagnosed functional bowel disorder
# K59 and laxative use). Metrics: Harrell's C (apparent and 10-fold cross-validated), likelihood-ratio test,
# continuous NRI at 5 years, and enrichment of cases among the highest predicted risks; bootstrap CIs.
# Usage: Rscript 17_prediction.R [B]   (B = bootstrap replicates, default 200)
suppressPackageStartupMessages(library(parallel))
source("/path/to/project/scripts/prep.R")
B <- if (length(commandArgs(TRUE))) as.integer(commandArgs(TRUE)[1]) else 200L
set.seed(20261008)
T5 <- 5

# ---------------------------------------------------------------- predictors
prev <- function(codes) { m <- as.matrix(oc[, ..codes]); as.integer(rowSums(!is.na(m) & (m <= 0)) > 0) }
co[, `:=`(hypotension = prev("I95"), urinary = prev(c("N31", "N32", "N39")), sleep = prev("G47"),
          constip_dx = prev("K59"))]
fh <- fread(paste0(P, "data/family_history_raw.csv"), colClasses = "character")
fh_cols <- grep("^p2011[01]_i|^p20107_i", names(fh), value = TRUE)
pd_fam <- Reduce(`|`, lapply(fh_cols, function(v) grepl("(^|[^0-9-])11([^0-9]|$)", fh[[v]])))
co[, family_pd := as.integer(pd_fam[match(co$eid, as.integer(fh$eid))] %in% TRUE)]
co[, smoke3 := factor(fifelse(as.character(smoking) %in% c("0", "1", "2"), as.character(smoking), "Missing"))]

tt <- outcome_tt(oc[["G20"]], FALSE)
d <- co[tt$sel, .(eid, age_q, sex, smoke3, coffee_cups, coffee_cups_miss, met_min, met_min_miss, family_pd,
                  prev_diabetes, prev_depression, prev_anxiety, hypotension, urinary, sleep, constip_dx, laxative_any,
                  hard_cat, loose_cat, bmf_cat)]
d[, `:=`(time = tt$time[tt$sel], event = tt$event[tt$sel])]
d <- d[!is.na(hard_cat) & !is.na(loose_cat)]
rm(oc); invisible(gc())

# (no missing-coffee indicator: that group had no PD events, which gives an infinite coefficient)
BASE <- paste("age_q + I(age_q^2) + sex + smoke3 + coffee_cups + met_min + met_min_miss + family_pd +",
              "prev_diabetes + prev_depression + prev_anxiety + hypotension + urinary + sleep + constip_dx + laxative_any")
MODELS <- c(base = BASE, freq = paste(BASE, "+ bmf_cat"), form = paste(BASE, "+ hard_cat + loose_cat"),
            all = paste(BASE, "+ hard_cat + loose_cat + bmf_cat"))
fitm <- function(f, data) coxph(as.formula(paste("Surv(time, event) ~", f)), data = data, x = FALSE)
risk5 <- function(m, newdata) {                       # predicted 5-year risk (centred lp with centred baseline hazard)
  bh <- basehaz(m, centered = TRUE); H <- approx(bh$time, bh$hazard, xout = T5, rule = 2)$y
  1 - exp(-H * exp(predict(m, newdata = newdata, type = "lp")))
}
km_risk <- function(time, event) { s <- summary(survfit(Surv(time, event) ~ 1), times = T5, extend = TRUE)$surv; 1 - s }
nri5 <- function(r_old, r_new, time, event) {         # continuous NRI at 5 years, KM-based (Pencina et al. 2011)
  up <- r_new > r_old; down <- r_new < r_old
  pD <- km_risk(time, event); pU <- mean(up); pDn <- mean(down)
  pD_up <- if (any(up)) km_risk(time[up], event[up]) else 0
  pD_dn <- if (any(down)) km_risk(time[down], event[down]) else 0
  ev <- (pD_up * pU - pD_dn * pDn) / pD
  ne <- ((1 - pD_dn) * pDn - (1 - pD_up) * pU) / (1 - pD)
  c(nri = ev + ne, nri_events = ev, nri_nonevents = ne)
}
cstat <- function(m) unname(concordance(m)$concordance)

# ---------------------------------------------------------------- apparent performance
fits <- lapply(MODELS, fitm, data = d)
cc <- concordance(fits$base, fits$freq, fits$form, fits$all)
ctr <- function(i) { v <- c(-1, rep(0, 3)); v[i] <- 1; v }
perf <- rbindlist(lapply(2:4, function(i) {
  k <- ctr(i); dC <- sum(k * cc$concordance); se <- sqrt(drop(t(k) %*% cc$var %*% k))
  lr <- anova(fits$base, fits[[i]])
  data.table(model = names(MODELS)[i], C_base = cc$concordance[1], C_model = cc$concordance[i], dC = dC,
             dC_lo = dC - 1.96 * se, dC_hi = dC + 1.96 * se, dC_p = 2 * pnorm(-abs(dC / se)),
             LR_chisq = lr$Chisq[2], LR_df = lr$Df[2], LR_p = lr$`Pr(>|Chi|)`[2])
}))
r_base <- risk5(fits$base, d); r_form <- risk5(fits$form, d); r_all <- risk5(fits$all, d); r_freq <- risk5(fits$freq, d)
nri_tab <- rbind(data.table(model = "freq", t(nri5(r_base, r_freq, d$time, d$event))),
                 data.table(model = "form", t(nri5(r_base, r_form, d$time, d$event))),
                 data.table(model = "all", t(nri5(r_base, r_all, d$time, d$event))))
perf <- merge(perf, nri_tab, by = "model")

# ---------------------------------------------------------------- 10-fold cross-validated C
fold <- sample(rep(1:10, length.out = nrow(d)))
cv_lp <- sapply(c("base", "form", "all"), function(k) {
  lp <- numeric(nrow(d))
  for (f in 1:10) { m <- fitm(MODELS[[k]], d[fold != f]); lp[fold == f] <- predict(m, newdata = d[fold == f], type = "lp") }
  lp
})
cv <- sapply(c("base", "form", "all"), function(k) concordance(Surv(d$time, d$event) ~ cv_lp[, k], reverse = TRUE)$concordance)

# ---------------------------------------------------------------- enrichment: top decile / top 5% of predicted risk
enrich <- rbindlist(lapply(list(base = r_base, form = r_form, all = r_all), function(r) {
  rbindlist(lapply(c(0.10, 0.05), function(q) {
    top <- r >= quantile(r, 1 - q)
    data.table(top_fraction = q, cases_captured = sum(d$event[top]), cases_total = sum(d$event),
               pct_cases_captured = 100 * mean(d$event[top] == 1) * sum(top) / sum(d$event),
               risk5_top_km = 100 * km_risk(d$time[top], d$event[top]), risk5_rest_km = 100 * km_risk(d$time[!top], d$event[!top]))
  }))
}), idcol = "model")
enrich[, pct_cases_captured := 100 * cases_captured / cases_total]

# ---------------------------------------------------------------- bootstrap CIs (delta C and NRI for the stool-form model)
boot1 <- function(b) {
  i <- sample.int(nrow(d), replace = TRUE); db <- d[i]
  m0 <- fitm(MODELS[["base"]], db); m1 <- fitm(MODELS[["form"]], db)
  c(dC = cstat(m1) - cstat(m0), nri5(risk5(m0, db), risk5(m1, db), db$time, db$event))
}
bt <- do.call(rbind, mclapply(seq_len(B), function(b) tryCatch(boot1(b), error = function(e) rep(NA, 4)), mc.cores = 4))
bci <- apply(bt, 2, quantile, c(0.025, 0.975), na.rm = TRUE)

# ---------------------------------------------------------------- output
out <- list(perf = perf, cv = data.table(model = names(cv), C_cv = cv), enrich = enrich,
            boot = data.table(stat = colnames(bci), lo = bci[1, ], hi = bci[2, ], B_ok = sum(complete.cases(bt))),
            sample = data.table(n = nrow(d), events = sum(d$event), family_pd = sum(d$family_pd),
                                hypotension = sum(d$hypotension), urinary = sum(d$urinary), sleep = sum(d$sleep),
                                constip_dx = sum(d$constip_dx), km_risk5 = 100 * km_risk(d$time, d$event)),
            base_hr = as.data.table(summary(fits$base)$conf.int[, c(1, 3, 4)], keep.rownames = "term"))
for (k in names(out)) fwrite(out[[k]], paste0(P, "results/pd_prediction_", k, ".csv"))
print(out$sample); print(out$perf); print(out$cv); print(out$boot); print(out$enrich)
