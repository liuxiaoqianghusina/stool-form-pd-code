# Accelerometry (UK Biobank wrist accelerometer sub-study, 2013-2015) and baseline medication codes (field 20003,
# instance 0) for the tier-2 analyses. Output: data/accelerometry.csv, data/meds_baseline_raw.csv.
import pandas as pd

D = "/path/to/ukb_data/"
P = "/path/to/project/"
acc = pd.read_csv(D + "UKB619657_additional_exposures_2025Nov.csv",
                  usecols=["eid", "p90010", "p90012", "p90013", "p90015", "p90016", "p90051"], dtype=str)
acc.to_csv(P + "data/accelerometry.csv", index=False)
print("accelerometry rows with average acceleration:", acc.p90012.notna().sum())
hdr = open(D + "UKB619657_assessment_center_2025Nov.csv").readline().replace('"', "").strip().split(",")
med = pd.read_csv(D + "UKB619657_assessment_center_2025Nov.csv",
                  usecols=["eid"] + [c for c in hdr if c.startswith("p20003_i0_a")], dtype=str)
med.to_csv(P + "data/meds_baseline_raw.csv", index=False)
print("medication rows:", len(med))
