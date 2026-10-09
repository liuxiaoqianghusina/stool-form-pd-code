# NHANES 2005-2010: treated Parkinson's disease (antiparkinsonian prescriptions) and usual stool form (Bristol), adults >= 40.
suppressPackageStartupMessages({library(data.table); library(survey)}); options(survey.lonely.psu = "adjust")
P <- "/path/to/project/"
d <- fread(paste0(P, "data/nhanes/nhanes_bowel.csv"), na.strings = c("", "NA"))
d[, ok := age >= 40 & !is.na(stool_form)]
d[, stool_form := factor(stool_form, levels = c("Normal (BSFS 3-5)", "Hard (BSFS 1-2)", "Loose (BSFS 6-7)"))]
d[ok == FALSE, `:=`(age = 60, female = 0L, stool_form = "Normal (BSFS 3-5)")]
des <- subset(svydesign(id = ~psu, strata = ~strata, weights = ~weight, nest = TRUE, data = d), ok)
out <- rbindlist(lapply(c(pd_specific_rx = "PD-specific drugs (levodopa/carbidopa, MAO-B inhibitors, COMT inhibitors)",
                          pd_any_rx = "Any antiparkinsonian drug (also dopamine agonists)"), function(lab) {
  y <- names(which(c(pd_specific_rx = "PD-specific drugs (levodopa/carbidopa, MAO-B inhibitors, COMT inhibitors)",
                     pd_any_rx = "Any antiparkinsonian drug (also dopamine agonists)") == lab))
  m <- svyglm(as.formula(paste(y, "~ stool_form + age + female")), design = des, family = quasibinomial())
  s <- summary(m)$coefficients[2:3, ]; ci <- confint(m)[2:3, ]
  cases <- d[ok == TRUE, .(cases = sum(get(y)), n = .N), by = stool_form]
  data.table(definition = lab, level = sub("stool_form", "", rownames(s)), OR = exp(s[, 1]), lo = exp(ci[, 1]),
             hi = exp(ci[, 2]), p = s[, 4], cases_total = d[ok == TRUE, sum(get(y))])
}))
fwrite(out, paste0(P, "results/nhanes_pd_meds.csv")); print(out)
