# Baseline (2006-2010) medication flags from field 20003 for PD sensitivity analyses: drugs that can cause
# drug-induced parkinsonism and anticholinergic (antimuscarinic) drugs. Output: data/meds_flags.csv.
import csv
import re

import pandas as pd

D = "/path/to/ukb_data/"
P = "/path/to/project/"
cod4 = {r[1]: r[2] for r in csv.reader(open(D + "app45551_20251118060954.dataset.codings.csv", encoding="utf-8"))
        if r and r[0] == "data_coding_4"}
GROUPS = {
    "dip_drug": r"haloperidol|chlorpromazine|trifluoperazine|fluphenazine|flupentixol|flupenthixol|zuclopenthixol|sulpiride|"
                r"amisulpride|risperidone|olanzapine|quetiapine|aripiprazole|clozapine|pericyazine|levomepromazine|pimozide|"
                r"perphenazine|promazine|metoclopramide|prochlorperazine|cinnarizine|flunarizine|lithium|valpro",
    "anticholinergic": r"oxybutynin|tolterodine|solifenacin|trospium|darifenacin|fesoterodine|propiverine|hyoscine|procyclidine|"
                       r"orphenadrine|benzatropine|benztropine|trihexyphenidyl|benzhexol|dicycloverine|dicyclomine|propantheline",
    "dopaminergic": r"levodopa|co-careldopa|co-beneldopa|sinemet|madopar|stalevo|ropinirole|pramipexole|rotigotine|rasagiline|"
                    r"selegiline|entacapone|amantadine",
}
codes = {g: {k for k, v in cod4.items() if re.search(rx, v, re.I)} for g, rx in GROUPS.items()}
med = pd.read_csv(P + "data/meds_baseline_raw.csv", dtype=str)
cols = [c for c in med.columns if c.startswith("p20003_")]
vals = med[cols]
out = pd.DataFrame({"eid": med.eid})
for g, cs in codes.items():
    out[g] = vals.isin(cs).any(axis=1).astype(int)
    print(g, len(cs), "codes;", int(out[g].sum()), "users;", sorted(cod4[c] for c in list(cs))[:8])
out.to_csv(P + "data/meds_flags.csv", index=False)
