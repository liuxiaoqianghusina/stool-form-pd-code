# Assemble the Chinese (default) or English (argument 'en') full-text manuscript (npj Parkinson's Disease order: Introduction, Results, Discussion, Methods),
# number citations by first appearance, and build the reference list from PubMed metadata so every reference is real.
# Sources: manuscript/src/{题目与摘要,引言,结果,讨论,方法}.md. Citations are written as {PMID,PMID,...} or {REF:key}
# for the few non-PubMed items.
# Output: manuscript/全文_中文初稿_npjPD.md and .docx (Word conversion via scratch md2docx-style writer below).
import json
import re
import subprocess
import time
import xml.etree.ElementTree as ET

import sys

P = "/path/to/project/manuscript/"
LANG = sys.argv[1] if len(sys.argv) > 1 else "zh"          # "zh" (Chinese draft) or "en" (English manuscript)
rd = lambda f: open(P + f, encoding="utf-8").read().strip()
if LANG == "en":
    intro, results, disc, methods = rd("src_en/introduction.md"), rd("src_en/results.md"), rd("src_en/discussion.md"), rd("src_en/methods.md")
    legends, ta, tables = rd("src_en/legends.md"), rd("src_en/title_abstract.md"), rd("src_en/tables.md")
    # supplementary figure legends belong to the Supplementary Information file, not the main manuscript
    legends = "\n\n".join(x for x in legends.split("\n\n") if not x.startswith("**Supplementary Fig."))
    title, abstract = ta.split("\n", 1)[0], "# Abstract" + ta.split("# Abstract", 1)[1]
else:
    intro, results, disc, methods = rd("src/引言.md"), rd("src/结果.md"), rd("src/讨论.md"), rd("src/方法.md")
    legends, ta, tables = rd("src/图注.md"), rd("src/题目与摘要.md"), rd("src/表格.md")
    title, abstract = ta.split("\n", 1)[0], "# 摘要" + ta.split("# 摘要", 1)[1]

body = "\n\n".join([intro, results, disc, methods])

# ---------------------------------------------------------------- number citations by first appearance
order = []
for m in re.finditer(r"\{([^{}]+)\}", body):
    for k in m.group(1).split(","):
        k = k.strip()
        if k not in order:
            order.append(k)
num = {k: i + 1 for i, k in enumerate(order)}


def compress(ns):
    ns = sorted(set(ns)); out = []; i = 0
    while i < len(ns):
        j = i
        while j + 1 < len(ns) and ns[j + 1] == ns[j] + 1:
            j += 1
        out.append(f"{ns[i]}–{ns[j]}" if j - i >= 2 else ",".join(str(x) for x in ns[i:j + 1]))
        i = j + 1
    return ",".join(out)


body = re.sub(r"\{([^{}]+)\}", lambda m: "^" + compress([num[k.strip()] for k in m.group(1).split(",")]) + "^", body)

# ---------------------------------------------------------------- reference metadata from PubMed (efetch XML)
pmids = [k for k in order if not k.startswith("REF:")]
xml = None
for a in range(6):
    r = subprocess.run(["curl", "-s", "--max-time", "120", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                        "--data-urlencode", "db=pubmed", "--data-urlencode", "id=" + ",".join(pmids),
                        "--data-urlencode", "retmode=xml"], capture_output=True, text=True)
    try:
        xml = ET.fromstring(r.stdout); break
    except ET.ParseError:
        time.sleep(3 + 3 * a)
meta = {}


def pages(pg):
    if not pg:
        return ""
    m = re.fullmatch(r"(\w*?)(\d+)-(\w*?)(\d+)", pg)
    if m and len(m.group(4)) < len(m.group(2)):
        end = m.group(2)[:len(m.group(2)) - len(m.group(4))] + m.group(4)
        return f"{m.group(1)}{m.group(2)}–{m.group(1)}{end}"
    return pg.replace("-", "–")


for art in xml.findall(".//PubmedArticle"):
    pmid = art.findtext(".//MedlineCitation/PMID")
    names = []
    for au in art.findall(".//AuthorList/Author"):
        if au.findtext("CollectiveName"):
            names.append(au.findtext("CollectiveName")); continue
        ini = " ".join(c + "." for c in (au.findtext("Initials") or ""))
        names.append(f"{au.findtext('LastName')}, {ini}".strip(", "))
    au_s = (names[0] + " et al.") if len(names) > 5 else (", ".join(names[:-1]) + " & " + names[-1] if len(names) > 1 else (names[0] if names else ""))
    au_s = au_s if au_s.endswith(".") else au_s + "."
    t = "".join(art.find(".//ArticleTitle").itertext()).strip().rstrip(".")
    j = art.findtext(".//Journal/ISOAbbreviation") or art.findtext(".//Journal/Title")
    y = art.findtext(".//JournalIssue/PubDate/Year") or (art.findtext(".//JournalIssue/PubDate/MedlineDate") or "")[:4]
    vol = art.findtext(".//JournalIssue/Volume") or ""
    pg = pages(art.findtext(".//Pagination/MedlinePgn"))
    if not pg:
        el = [e.text for e in art.findall(".//ELocationID") if e.get("EIdType") != "doi"]
        pg = el[0] if el else ""
    doi = next((e.text for e in art.findall(".//ArticleIdList/ArticleId") if e.get("IdType") == "doi"), "")
    loc = f" **{vol}**, {pg}" if vol else (f" {pg}" if pg else "")
    t = t if t.endswith(("?", "!")) else t + "."
    meta[pmid] = f"{au_s} {t} *{j}*{loc} ({y}). doi:{doi}" if doi else f"{au_s} {t} *{j}*{loc} ({y})."
MANUAL = {
    "REF:finegray": "Fine, J. P. & Gray, R. J. A proportional hazards model for the subdistribution of a competing risk. "
                    "*J. Am. Stat. Assoc.* **94**, 496–509 (1999).",
    "REF:bh": "Benjamini, Y. & Hochberg, Y. Controlling the false discovery rate: a practical and powerful approach to multiple "
              "testing. *J. R. Stat. Soc. B* **57**, 289–300 (1995).",
    "REF:nchs": "National Center for Health Statistics. 2019 Public-Use Linked Mortality Files. "
                "https://www.cdc.gov/nchs/data-linkage/mortality-public.htm",
}
missing = [k for k in order if k not in meta and k not in MANUAL]
assert not missing, f"PubMed metadata missing for {missing}"
refs = "\n".join(f"{num[k]}. {meta.get(k) or MANUAL[k]}" for k in order)

if LANG == "en":
    words = lambda t: len(re.findall(r"[A-Za-z0-9][A-Za-z0-9'’.\-–]*", re.sub(r"\^[^^]+\^", "", t)))
    main_words = sum(words(x) for x in (intro, results, disc))
    note = (f"(Draft for npj Parkinson's Disease. Abstract {words(abstract.split(chr(10), 1)[1])} words; main text "
            f"[Introduction, Results, Discussion] {main_words} words; Methods {words(methods)} words. Highlighted 【】 items "
            "are to be completed or verified. References were generated from PubMed and will be formatted with reference software.)")
    full = f"{title}\n\n{note}\n\n{abstract.strip()}\n\n{body}\n\n# References\n\n{refs}\n\n{tables}\n\n\\pagebreak\n\n{legends}\n"
    open(P + "Manuscript_npjPD_EN.md", "w", encoding="utf-8").write(full)
    print(note)
else:
    note = ("（中文全文初稿，按 npj Parkinson's Disease 的 Article 结构排列：摘要—引言—结果—讨论—方法。引用编号按首次出现顺序排列；"
            "参考文献信息从 PubMed 自动生成，投稿前用 EndNote/Zotero 统一为期刊格式。黄色高亮的【】内容是待补充或待核实之处。）")
    full = f"{title}\n\n{note}\n\n{abstract.strip()}\n\n{body}\n\n# 参考文献\n\n{refs}\n\n{tables}\n\n\\pagebreak\n\n{legends}\n"
    open(P + "全文_中文初稿_npjPD.md", "w", encoding="utf-8").write(full)
    parts = re.split(r"(?m)^(?=# )", body)
    rm = "".join(x for x in parts if x.startswith(("# 结果", "# 方法")))
    open(P + "方法与结果_中文初稿_npjPD.md", "w", encoding="utf-8").write(f"{title}\n\n{note}\n\n{rm.strip()}\n\n# 参考文献\n\n{refs}\n")
print(f"{len(order)} references ({len(pmids)} PubMed, {len(order) - len(pmids)} manual)")
print(re.sub(r"\n.*", "", refs[:0]))
for k in order[:3]:
    print(num[k], meta.get(k) or MANUAL[k])
