# Bayesian colocalisation (approximate Bayes factor method of Giambartolomei et al. 2014, as in coloc.abf with default
# priors p1 = p2 = 1e-4, p12 = 1e-5) between stool frequency (Diaz-Munoz 2026, GRCh38) and Parkinson's disease
# (FinnGen R12, GRCh38) in a +/-500 kb window around each stool-frequency lead variant.
import gzip

import numpy as np
import pandas as pd
from scipy.special import logsumexp

P = "/path/to/project/"
MR = P + "data/mr/"
N_SF, N_CASE, N_CTRL = 167966, 5861, 494487
P1, P2, P12 = 1e-4, 1e-4, 1e-5
WINDOW = 500_000

leads = pd.read_csv(P + "results/mr_instruments.csv")


def read_region(path, chrom_col, pos_col, usecols, regions, sep="\t"):
    """Stream a gzipped summary-statistics file and keep rows inside any (chrom, start, end) window."""
    keep = []
    for chunk in pd.read_csv(path, sep=sep, usecols=usecols, chunksize=2_000_000, dtype={chrom_col: str}):
        chunk[chrom_col] = chunk[chrom_col].astype(str)
        m = np.zeros(len(chunk), dtype=bool)
        for c, s, e in regions:
            m |= (chunk[chrom_col] == str(c)) & chunk[pos_col].between(s, e)
        if m.any():
            keep.append(chunk[m])
    return pd.concat(keep, ignore_index=True)


regions = [(int(r.chr), int(r.pos) - WINDOW, int(r.pos) + WINDOW) for r in leads.itertuples()]
sf = read_region(MR + "GCST90565484.h.tsv.gz", "chromosome", "base_pair_location",
                 ["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error",
                  "effect_allele_frequency", "rsid"], regions)
pdg = read_region(MR + "finngen_r12/finngen_R12_G6_PARKINSON.gz", "#chrom", "pos",
                  ["#chrom", "pos", "ref", "alt", "beta", "sebeta", "af_alt"], regions)
pdg = pdg.rename(columns={"#chrom": "chromosome", "pos": "base_pair_location"})
m = sf.merge(pdg, on=["chromosome", "base_pair_location"])
same = (m.effect_allele == m.alt) & (m.other_allele == m.ref)
flip = (m.effect_allele == m.ref) & (m.other_allele == m.alt)
flipped = flip[same | flip].values
m = m[same | flip].copy()
m["b_pd"] = np.where(flipped, -m.beta_y, m.beta_y)  # PD log-OR aligned to the stool-frequency effect allele (beta_x)
m["maf"] = np.minimum(m.effect_allele_frequency, 1 - m.effect_allele_frequency)
m = m[(m.maf > 0.01) & m.standard_error.notna() & m.sebeta.notna()]

# sdY for the quantitative trait, estimated as in coloc::sdY.est: varbeta ~ sdY^2 / (2 * N * maf * (1 - maf))
sdY = np.sqrt(np.median(m.standard_error ** 2 * 2 * N_SF * m.maf * (1 - m.maf)))


def labf(beta, se, W):
    V = se ** 2
    r = W / (W + V)
    z = beta / se
    return 0.5 * (np.log(1 - r) + r * z ** 2)


rows = []
for r in leads.itertuples():
    reg = m[(m.chromosome.astype(str) == str(int(r.chr))) & m.base_pair_location.between(int(r.pos) - WINDOW, int(r.pos) + WINDOW)]
    if len(reg) < 50:
        continue
    l1 = labf(reg.beta_x.values, reg.standard_error.values, (0.15 * sdY) ** 2)   # stool frequency (quantitative)
    l2 = labf(reg.b_pd.values, reg.sebeta.values, 0.2 ** 2)                       # PD (case-control, log-OR scale)
    s1, s2, s12 = logsumexp(l1), logsumexp(l2), logsumexp(l1 + l2)
    lH = np.array([0.0, np.log(P1) + s1, np.log(P2) + s2,
                   np.log(P1) + np.log(P2) + s1 + s2 + np.log1p(-np.exp(s12 - s1 - s2)),
                   np.log(P12) + s12])
    pp = np.exp(lH - logsumexp(lH))
    top = reg.iloc[np.argmax(l1 + l2)]
    rows.append(dict(lead=r.rsid, chr=r.chr, pos=r.pos, n_snps=len(reg), PP_H0=pp[0], PP_H1=pp[1], PP_H2=pp[2], PP_H3=pp[3],
                     PP_H4=pp[4], top_shared_snp=top.rsid, top_shared_pos=int(top.base_pair_location),
                     min_p_pd=float(2 * __import__("scipy").stats.norm.sf(np.abs(reg.b_pd / reg.sebeta)).min())))
res = pd.DataFrame(rows)
res.to_csv(P + "results/coloc_stoolfreq_pd.csv", index=False)
pd.set_option("display.width", 200)
print(f"sdY (stool frequency) = {sdY:.3f}; variants matched = {len(m)}")
print(res.round(3).to_string(index=False))
