# Two-sample Mendelian randomization: genetically predicted stool frequency (Diaz-Munoz 2026 Gut, European meta-analysis,
# GCST90565484, n = 167,966) -> FinnGen R12 outcomes (no UK Biobank overlap), plus reverse MR (PD liability -> stool frequency).
# Instruments: genome-wide significant lead variants >= 10 Mb apart (independent by distance; avoids an LD reference panel).
suppressPackageStartupMessages({library(data.table); library(TwoSampleMR); library(MRPRESSO)})
P  <- "/path/to/project/"
MR <- paste0(P, "data/mr/")
OUTCOMES <- c(G6_PARKINSON = "Parkinson's disease", T2D = "Type 2 diabetes", F5_DEPRESSIO = "Depression",
              G6_ALZHEIMER = "Alzheimer's disease", F5_DEMENTIA = "Dementia", F5_ALLANXIOUS = "Anxiety",
              I9_HYPTENS = "Hypertension", I9_HEARTFAIL = "Heart failure")
FLAG <- c(rs2732706 = "17q21.31 MAPT region (PD risk locus)", rs12273363 = "11p14 BDNF region")

lead_by_distance <- function(d, chr = "chr", pos = "pos", p = "p", window = 1e7) {
  o <- order(d[[p]]); d <- d[o]; keep <- logical(nrow(d))  # index computed outside [ ] (data.table scoping)
  for (i in seq_len(nrow(d))) {
    prev <- d[keep]
    keep[i] <- !any(prev[[chr]] == d[[chr]][i] & abs(prev[[pos]] - d[[pos]][i]) < window)
  }
  d[keep]
}

# ---------------------------------------------------------------- exposure instruments
ex <- fread(paste0(MR, "exposure_gws.tsv"))
setnames(ex, c("chromosome", "base_pair_location", "p_value"), c("chr", "pos", "p"))
ex <- lead_by_distance(ex)
ex[, F := (beta / standard_error)^2]
exp_dat <- format_data(as.data.frame(ex), type = "exposure", snp_col = "rsid", beta_col = "beta", se_col = "standard_error",
                       eaf_col = "effect_allele_frequency", effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                       pval_col = "p", chr_col = "chr", pos_col = "pos")
exp_dat$samplesize.exposure <- 167966; exp_dat$exposure <- "Stool frequency"

read_finngen_at <- function(code, chr, pos) {
  f <- paste0(MR, "finngen_r12/finngen_R12_", code, ".gz")
  key <- paste(chr, pos)
  x <- fread(cmd = sprintf("gzcat '%s'", f), select = c("#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"))
  setnames(x, "#chrom", "chrom")
  x[paste(chrom, pos) %in% key]
}

results <- list(); loo <- list(); snpwise <- list()
for (code in names(OUTCOMES)) {
  og <- read_finngen_at(code, ex$chr, ex$pos)
  og <- merge(og, ex[, .(chr, pos, rsid, ea = effect_allele, oa = other_allele)], by.x = c("chrom", "pos"), by.y = c("chr", "pos"))
  og <- og[(ref == oa & alt == ea) | (ref == ea & alt == oa)]  # same variant, either orientation
  out_dat <- format_data(as.data.frame(og), type = "outcome", snp_col = "rsid", beta_col = "beta", se_col = "sebeta",
                         eaf_col = "af_alt", effect_allele_col = "alt", other_allele_col = "ref", pval_col = "pval")
  out_dat$outcome <- OUTCOMES[[code]]
  h <- harmonise_data(exp_dat, out_dat, action = 2)
  h <- h[h$mr_keep, ]
  for (set in c("All instruments", "Excluding MAPT-region variant", "Excluding MAPT and BDNF variants")) {
    hh <- switch(set, "All instruments" = h, "Excluding MAPT-region variant" = h[h$SNP != "rs2732706", ],
                 "Excluding MAPT and BDNF variants" = h[!h$SNP %in% names(FLAG), ])
    r <- mr(hh, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
    r$instrument_set <- set
    het <- mr_heterogeneity(hh, method_list = "mr_ivw"); ple <- mr_pleiotropy_test(hh)
    r$Q_pval_ivw <- het$Q_pval[1]; r$egger_intercept <- ple$egger_intercept; r$egger_intercept_p <- ple$pval
    results[[length(results) + 1]] <- r
  }
  pr <- tryCatch(mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome",
                           SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = h,
                           NbDistribution = 3000, SignifThreshold = 0.05), error = function(e) NULL)
  if (!is.null(pr)) {
    mt <- pr$`Main MR results`
    results[[length(results) + 1]] <- data.frame(outcome = OUTCOMES[[code]], exposure = "Stool frequency",
      method = paste("MR-PRESSO", mt$`MR Analysis`), nsnp = nrow(h), b = mt$`Causal Estimate`, se = mt$Sd, pval = mt$`P-value`,
      instrument_set = "All instruments", Q_pval_ivw = NA, egger_intercept = NA,
      egger_intercept_p = pr$`MR-PRESSO results`$`Global Test`$Pvalue)
  }
  l <- mr_leaveoneout(h); l$outcome <- OUTCOMES[[code]]; loo[[code]] <- l
  s <- mr_singlesnp(h); s$outcome <- OUTCOMES[[code]]; snpwise[[code]] <- s
}
res <- rbindlist(results, fill = TRUE)
res[, `:=`(OR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se))]
fwrite(res, paste0(P, "results/mr_stoolfreq_finngen.csv"))
fwrite(rbindlist(loo, fill = TRUE), paste0(P, "results/mr_leaveoneout.csv"))
fwrite(rbindlist(snpwise, fill = TRUE), paste0(P, "results/mr_singlesnp.csv"))
fwrite(ex[, .(rsid, chr, pos, effect_allele, other_allele, beta, standard_error, effect_allele_frequency, p, F)],
       paste0(P, "results/mr_instruments.csv"))

# ---------------------------------------------------------------- reverse MR: PD liability (FinnGen) -> stool frequency
pd <- fread(cmd = sprintf("gzcat '%s'", paste0(MR, "finngen_r12/finngen_R12_G6_PARKINSON.gz")),
            select = c("#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"))
setnames(pd, "#chrom", "chrom")
pd_iv <- lead_by_distance(pd[pval < 5e-8], chr = "chrom", pos = "pos", p = "pval")
pd_iv[, rsid := sub(",.*", "", rsids)]
exp_pd <- format_data(as.data.frame(pd_iv), type = "exposure", snp_col = "rsid", beta_col = "beta", se_col = "sebeta",
                      eaf_col = "af_alt", effect_allele_col = "alt", other_allele_col = "ref", pval_col = "pval")
exp_pd$exposure <- "Parkinson's disease liability"
sf <- fread(cmd = sprintf("gzcat '%s'", paste0(MR, "GCST90565484.h.tsv.gz")),
            select = c("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error",
                       "effect_allele_frequency", "p_value", "rsid"))
sf <- merge(sf, pd_iv[, .(chrom, pos, ref, alt, rsid_pd = rsid)], by.x = c("chromosome", "base_pair_location"), by.y = c("chrom", "pos"))
sf <- sf[(effect_allele == alt & other_allele == ref) | (effect_allele == ref & other_allele == alt)]
out_sf <- format_data(as.data.frame(sf), type = "outcome", snp_col = "rsid_pd", beta_col = "beta", se_col = "standard_error",
                      eaf_col = "effect_allele_frequency", effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                      pval_col = "p_value")
out_sf$outcome <- "Stool frequency"
hr <- harmonise_data(exp_pd, out_sf, action = 2); hr <- hr[hr$mr_keep, ]
rev <- mr(hr, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
rev$Q_pval_ivw <- mr_heterogeneity(hr, method_list = "mr_ivw")$Q_pval[1]
fwrite(as.data.table(rev), paste0(P, "results/mr_reverse_pd_to_stoolfreq.csv"))
cat("Forward MR rows:", nrow(res), "| reverse MR instruments:", nrow(hr), "\n")
print(res[instrument_set == "All instruments" & method %in% c("Inverse variance weighted", "Weighted median", "MR Egger"),
          .(outcome, method, nsnp, OR = round(OR, 3), lo = round(lo, 3), hi = round(hi, 3), p = signif(pval, 3))])
print(as.data.table(rev)[, .(method, nsnp, b = round(b, 4), se = round(se, 4), p = signif(pval, 3))])
