# Reverse Mendelian randomization: genetic liability to Parkinson's disease -> stool frequency.
# Exposure: Nalls et al. 2019 PD GWAS (GCST009325, harmonised GRCh38); replication with FinnGen R12 PD.
# Outcome: stool frequency (Diaz-Munoz 2026 Gut, GCST90565484). Instruments: genome-wide significant lead variants
# >= 1 Mb apart, estimated in the full meta-analysis (Nalls: N cases > 10,000). Sensitivity: excluding 17q21.31 (MAPT).
suppressPackageStartupMessages({library(data.table); library(TwoSampleMR); library(MRPRESSO)})
P  <- "/path/to/project/"
MR <- paste0(P, "data/mr/")
MAPT <- c(chr = 17, start = 45.4e6, end = 46.8e6)  # GRCh38 inversion region

lead_by_distance <- function(d, window = 1e6) {
  o <- order(d$p); d <- d[o]; keep <- logical(nrow(d))
  for (i in seq_len(nrow(d))) keep[i] <- !any(d$chr[keep] == d$chr[i] & abs(d$pos[keep] - d$pos[i]) < window)
  d[keep]
}

sf <- fread(cmd = sprintf("gzcat '%s'", paste0(MR, "GCST90565484.h.tsv.gz")),
            select = c("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error",
                       "effect_allele_frequency", "p_value", "rsid"))
setnames(sf, c("chromosome", "base_pair_location", "p_value"), c("chr", "pos", "p"))
sf[, chr := as.character(chr)]

run_reverse <- function(iv, label) {
  iv[, chr := as.character(chr)]
  o <- merge(sf, iv[, .(chr, pos, ea_x = ea, oa_x = oa, snp_id = snp)], by = c("chr", "pos"))
  o <- o[(effect_allele == ea_x & other_allele == oa_x) | (effect_allele == oa_x & other_allele == ea_x)]
  exp_dat <- format_data(as.data.frame(iv), type = "exposure", snp_col = "snp", beta_col = "beta", se_col = "se",
                         eaf_col = "eaf", effect_allele_col = "ea", other_allele_col = "oa", pval_col = "p")
  exp_dat$exposure <- label
  out_dat <- format_data(as.data.frame(o), type = "outcome", snp_col = "snp_id", beta_col = "beta", se_col = "standard_error",
                         eaf_col = "effect_allele_frequency", effect_allele_col = "effect_allele",
                         other_allele_col = "other_allele", pval_col = "p")
  out_dat$outcome <- "Stool frequency"; out_dat$samplesize.outcome <- 167966
  h <- harmonise_data(exp_dat, out_dat, action = 2); h <- h[h$mr_keep, ]
  res <- list()
  for (set in c("All instruments", "Excluding MAPT region")) {
    hh <- if (set == "All instruments") h else h[!(h$SNP %in% iv[chr == "17" & pos >= MAPT["start"] & pos <= MAPT["end"], snp]), ]
    r <- mr(hh, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
    r$instrument_set <- set
    r$Q_pval_ivw <- mr_heterogeneity(hh, method_list = "mr_ivw")$Q_pval[1]
    r$egger_intercept_p <- mr_pleiotropy_test(hh)$pval
    res[[set]] <- r
  }
  pr <- tryCatch(mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome",
                           SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = h,
                           NbDistribution = 5000, SignifThreshold = 0.05), error = function(e) NULL)
  out <- rbindlist(res, fill = TRUE)
  if (!is.null(pr)) {
    mt <- as.data.table(pr$`Main MR results`)
    out <- rbind(out, data.table(exposure = label, outcome = "Stool frequency", method = paste("MR-PRESSO", mt$`MR Analysis`),
                                 nsnp = nrow(h), b = mt$`Causal Estimate`, se = mt$Sd, pval = mt$`P-value`,
                                 instrument_set = "All instruments",
                                 egger_intercept_p = pr$`MR-PRESSO results`$`Global Test`$Pvalue), fill = TRUE)
  }
  st <- tryCatch(directionality_test(transform(h, samplesize.exposure = iv$n[1], ncase.exposure = iv$ncase[1],
                                               ncontrol.exposure = iv$ncontrol[1], units.exposure = "log odds",
                                               prevalence.exposure = 0.005, units.outcome = "SD")),
                 error = function(e) NULL)
  if (!is.null(st)) out[, steiger_correct_direction := st$correct_causal_direction[1]]
  out
}

# ---------------------------------------------------------------- Nalls 2019
na <- fread(cmd = sprintf("gzcat '%s'", paste0(MR, "GCST009325_nalls2019.h.tsv.gz")),
            select = c("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error",
                       "effect_allele_frequency", "p_value", "rsid", "N_cases", "N_controls"))
na <- na[p_value < 5e-8 & N_cases > 10000 & !is.na(rsid)]
setnames(na, c("chromosome", "base_pair_location", "effect_allele", "other_allele", "standard_error", "effect_allele_frequency",
               "p_value", "rsid"), c("chr", "pos", "ea", "oa", "se", "eaf", "p", "snp"))
iv_nalls <- lead_by_distance(na)
iv_nalls[, `:=`(ncase = N_cases, ncontrol = N_controls, n = N_cases + N_controls)]
res_nalls <- run_reverse(iv_nalls, "PD liability (Nalls 2019)")

# ---------------------------------------------------------------- FinnGen R12 replication
fg <- fread(cmd = sprintf("gzcat '%s'", paste0(MR, "finngen_r12/finngen_R12_G6_PARKINSON.gz")),
            select = c("#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"))
setnames(fg, c("#chrom", "ref", "alt", "pval", "sebeta", "af_alt"), c("chr", "oa", "ea", "p", "se", "eaf"))
fg <- fg[p < 5e-8]
fg[, snp := sub(",.*", "", rsids)]
iv_fg <- lead_by_distance(fg)
iv_fg[, `:=`(ncase = 5861, ncontrol = 494487, n = 500348)]
res_fg <- run_reverse(iv_fg, "PD liability (FinnGen R12)")

res <- rbind(res_nalls, res_fg, fill = TRUE)
fwrite(res, paste0(P, "results/mr_reverse_pd_to_stoolfreq.csv"))
fwrite(rbind(iv_nalls[, .(source = "Nalls 2019", snp, chr, pos, ea, oa, beta, se, p)],
             iv_fg[, .(source = "FinnGen R12", snp, chr, pos, ea, oa, beta, se, p)]),
       paste0(P, "results/mr_reverse_instruments.csv"))
print(res[, .(exposure, instrument_set, method, nsnp, b = signif(b, 3), se = signif(se, 3), p = signif(pval, 3),
              Q_p = signif(Q_pval_ivw, 3), egger_p = signif(egger_intercept_p, 3))])
if ("steiger_correct_direction" %in% names(res)) print(unique(res[, .(exposure, steiger_correct_direction)]))
