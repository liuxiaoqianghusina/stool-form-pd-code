# One-off renumbering of supplementary items in the manuscript sources (Chinese and English) so that they are numbered
# in order of first citation and split into Supplementary Tables (in the Supplementary Information file) and
# Supplementary Data (separate Excel workbook). The mapping is saved to manuscript/supplement_map.json and reused by
# 12_tables.py to build the submission files. A marker file prevents the text from being renumbered twice.
import json
import os
import re

P = "/path/to/project/manuscript/"
TABLE = {"15": 1, "21": 2, "23": 3, "20": 4, "19": 5, "17": 6, "22": 7, "18": 8, "24": 9, "25": 10, "16": 11, "9": 12,
         "10": 13, "13": 14, "12": 15, "1": 16}
DATA = {"2": 1, "3": 2, "4a": 3, "4b": 4, "7": 5, "5": 6, "6": 7, "8": 8, "11": 9, "14": 10}
FIG = {"1": 2, "2": 1, "3": 3}
json.dump({"table": TABLE, "data": DATA, "fig": FIG}, open(P + "supplement_map.json", "w"), indent=1)
MARK = P + ".supplement_renumbered"
if os.path.exists(MARK):
    raise SystemExit("already renumbered; mapping file refreshed only")


def expand(spec):
    out = []
    for part in re.split(r"\s*(?:,|、|和| and )\s*", spec):
        if "–" in part:
            a, b = part.split("–"); out += [str(i) for i in range(int(a), int(b) + 1)]
        else:
            out.append(part)
    return [y for x in out for y in (["4a", "4b"] if x == "4" else [x])]


def compress(nums, lang):
    nums = sorted(nums); runs = []; i = 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        runs.append(f"{nums[i]}–{nums[j]}" if j > i else str(nums[i])); i = j + 1
    if lang == "zh":
        return "、".join(runs)
    return runs[0] if len(runs) == 1 else ", ".join(runs[:-1]) + " and " + runs[-1]


def relabel(spec, lang):
    items = expand(spec)
    t = [TABLE[x] for x in items if x in TABLE]; d = [DATA[x] for x in items if x in DATA]
    assert len(t) + len(d) == len(items), (spec, items)
    plural_t = len(t) > 1
    parts = []
    if t:
        parts.append(("补充表 " if lang == "zh" else ("Supplementary Tables " if plural_t else "Supplementary Table ")) + compress(t, lang))
    if d:
        parts.append(("补充数据 " if lang == "zh" else "Supplementary Data ") + compress(d, lang))
    return ("；" if lang == "zh" else "; ").join(parts)


NUM = r"\d+(?:\s*(?:–|、|和| and |, )\s*\d+)*"
for lang, files in [("zh", ["src/结果.md", "src/讨论.md", "src/方法.md", "src/图注.md"]),
                    ("en", ["src_en/results.md", "src_en/discussion.md", "src_en/methods.md", "src_en/legends.md"])]:
    for f in files:
        s = open(P + f, encoding="utf-8").read()
        if lang == "zh":
            s = re.sub(r"补充表 (" + NUM + ")", lambda m: relabel(m.group(1), "zh"), s)
            s = re.sub(r"补充图 (\d)", lambda m: f"补充图 {FIG[m.group(1)]}", s)
        else:
            s = re.sub(r"Supplementary Tables? (" + NUM + ")", lambda m: relabel(m.group(1), "en"), s)
            s = re.sub(r"Supplementary Fig\. (\d)", lambda m: f"Supplementary Fig. {FIG[m.group(1)]}", s)
        if f.endswith(("图注.md", "legends.md")):      # put supplementary figure legends back in numerical order
            paras = s.split("\n\n")
            key = (lambda p: re.match(r"\*\*(?:补充图|Supplementary Fig\.) (\d)", p))
            supp = sorted([p for p in paras if key(p)], key=lambda p: int(key(p).group(1)))
            it = iter(supp); s = "\n\n".join(next(it) if key(p) else p for p in paras)
        open(P + f, "w", encoding="utf-8").write(s)
open(MARK, "w").write("done\n")
print("renumbered")
