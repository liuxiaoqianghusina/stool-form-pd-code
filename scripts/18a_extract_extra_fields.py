# Extra fields for the tier-1 analyses: abdominal pain/distension items of the Digestive Health Questionnaire
# (negative-control exposures) and alternative Parkinson's disease definitions (algorithmic outcome 42032/42033 and the
# source of the first-occurrence G20 report 131023). Output: data/extra_fields.csv (eid + raw values).
import pandas as pd

D = "/path/to/ukb_data/"
P = "/path/to/project/"
a = pd.read_csv(D + "UKB619657_online_followup_2025Nov.csv",
                usecols=["eid", "p21025", "p21035", "p21036", "p21038", "p21039"], dtype=str)
b = pd.read_csv(D + "UKB619657_health_outcomes_2025Nov.csv", usecols=["eid", "p42032", "p42033", "p131023"], dtype=str)
d = a.merge(b, on="eid", how="outer")
d.to_csv(P + "data/extra_fields.csv", index=False)
print(len(d)); print(d.notna().sum().to_string())
print(d.p42032.dropna().sort_values().tail(3).tolist())
