# Main-text tables (Table 1 and Table 2) as Markdown for the manuscripts, placed after the references (Nature style).
# Reads the tables written by 12_tables.py and writes manuscript/src_en/tables.md (English) and manuscript/src/表格.md
# (Chinese labels). Sub-rows are indented with an em space, which scripts/md2docx.py keeps.
import pandas as pd

P = "/path/to/project/manuscript/"
x = pd.ExcelFile(P + "Supplementary_Tables.xlsx")


def load(sheet):
    d = pd.read_excel(x, sheet, dtype=str).fillna("")
    first = d.columns[0]
    d[first] = d[first].str.strip()          # sub-rows carry leading spaces in the source table; indent is re-added below
    return d[(d[first] != "") | (d.iloc[:, 1:] != "").any(axis=1)].loc[lambda t: ~t[first].str.startswith("Abbreviations")]


t1, t2 = load("Table 1"), load("Table 2")
SUB = {"Never", "Previous", "Current", "Laxatives", "Proton pump inhibitors", "Metformin", "Opioids", "Antidepressants",
       "Diabetes", "Hypertension", "Ischaemic heart disease", "Stroke", "Depression", "Anxiety", "Irritable bowel syndrome",
       "Chronic kidney disease", "Cancer", "Hard/lumpy stools often or more", "Loose/watery stools often or more"}
EM = " "

NOTE_EN = {
    1: ("Values are n (%), mean (SD) or median (IQR). Lifestyle and medication data are from the UK Biobank baseline assessment "
        "(2006–2010); conditions are those recorded before the Digestive Health Questionnaire (2017–2020). Participants with IBD, "
        "coeliac disease or malabsorption, colorectal cancer or prior colectomy were excluded. IBD, inflammatory bowel disease; "
        "IQR, interquartile range; MET, metabolic equivalent of task; SD, standard deviation."),
    2: ("Model 1: age and sex. Model 2: additionally adjusted for ethnicity, Townsend index, education, smoking, alcohol, physical "
        "activity, BMI, diet (fruit, vegetables, wholemeal bread, fibre cereal), water intake, number of medications, laxatives, "
        "metformin, PPIs, opioids, antidepressants and prevalent diabetes, hypertension, ischaemic heart disease, stroke, "
        "depression, anxiety, IBS, CKD and cancer. †Linear trend across the five response levels (never, sometimes, often, most "
        "of the time, always), model 2. BMI, body mass index; CI, confidence interval; CKD, chronic kidney disease; HR, hazard "
        "ratio; IBS, irritable bowel syndrome; PPI, proton pump inhibitor; PY, person-years."),
}
TITLE_EN = {1: "Table 1 Characteristics of participants by average bowel movement frequency",
            2: "Table 2 Stool form and incident Parkinson's disease and dementia"}

ZH = {
    "Characteristic": "特征", "Overall": "合计", "<1 per day": "<1 次/天", "1 per day": "1 次/天", "2 per day": "2 次/天",
    "≥3 per day": "≥3 次/天", "Participants, n": "人数，n", "Age at questionnaire, years, mean (SD)": "问卷时年龄，岁，均值（SD）",
    "Female, n (%)": "女性，n（%）", "White ethnicity, n (%)": "白人，n（%）",
    "Townsend deprivation index, median (IQR)": "Townsend 剥夺指数，中位数（IQR）",
    "College or university degree, n (%)": "大学学历，n（%）", "Smoking status, n (%)": "吸烟状态，n（%）",
    "Never": "从不", "Previous": "既往", "Current": "当前",
    "Alcohol daily or almost daily, n (%)": "每天或几乎每天饮酒，n（%）",
    "Physical activity, MET-min/week, median (IQR)": "体力活动，MET-分钟/周，中位数（IQR）",
    "Body mass index, kg/m², mean (SD)": "体重指数，kg/m²，均值（SD）",
    "Fresh fruit, pieces/day, median (IQR)": "新鲜水果，份/天，中位数（IQR）",
    "Vegetables, tablespoons/day, median (IQR)": "蔬菜，汤匙/天，中位数（IQR）", "Wholemeal bread, n (%)": "食用全麦面包，n（%）",
    "Water, glasses/day, median (IQR)": "饮水，杯/天，中位数（IQR）", "Coffee, cups/day, median (IQR)": "咖啡，杯/天，中位数（IQR）",
    "Medication use, n (%)": "用药，n（%）", "Laxatives": "泻药", "Proton pump inhibitors": "质子泵抑制剂", "Metformin": "二甲双胍",
    "Opioids": "阿片类药物", "Antidepressants": "抗抑郁药", "Conditions before the questionnaire, n (%)": "问卷前已有疾病，n（%）",
    "Diabetes": "糖尿病", "Hypertension": "高血压", "Ischaemic heart disease": "缺血性心脏病", "Stroke": "卒中",
    "Depression": "抑郁", "Anxiety": "焦虑", "Irritable bowel syndrome": "肠易激综合征", "Chronic kidney disease": "慢性肾病",
    "Cancer": "癌症", "Stool form in the past 3 months, n (%)": "过去 3 个月的大便性状，n（%）",
    "Hard/lumpy stools often or more": "经常及以上出现硬便或块状便", "Loose/watery stools often or more": "经常及以上出现稀便",
    "Outcome": "结局", "Exposure": "暴露", "Frequency": "频率", "Events/N": "事件数/人数", "Rate per 1,000 PY": "发病率（/1,000 人年）",
    "HR (95% CI), model 1": "HR（95% CI），模型 1", "HR (95% CI), model 2": "HR（95% CI），模型 2", "P for trend†": "*P* 趋势†",
    "Parkinson's disease": "帕金森病", "Dementia": "痴呆", "Hard/lumpy stools": "硬便或块状便", "Loose/watery stools": "稀便",
    "Sometimes": "有时", "Often or more": "经常及以上", "1.00 (reference)": "1.00（参照）",
}
NOTE_ZH = {
    1: ("数值为 n（%）、均值（SD）或中位数（IQR）。生活方式和用药信息来自 UK Biobank 基线评估（2006–2010 年），已有疾病为消化健康问卷"
        "（2017–2020 年）前已记录者。已排除炎症性肠病、乳糜泻或吸收不良、结直肠癌或既往结直肠切除者。缩写：IQR，四分位距；MET，代谢当量；"
        "SD，标准差。"),
    2: ("模型 1 调整年龄和性别；模型 2 进一步调整种族、Townsend 剥夺指数、教育程度、吸烟、饮酒、体力活动、体重指数、饮食（水果、蔬菜、"
        "全麦面包、高纤维谷物）、饮水量、用药总数、泻药、二甲双胍、质子泵抑制剂、阿片类药物、抗抑郁药，以及既往糖尿病、高血压、缺血性"
        "心脏病、卒中、抑郁、焦虑、肠易激综合征、慢性肾病和癌症。†为模型 2 中跨五个应答等级（从不、有时、经常、大多数时候、总是）的线性"
        "趋势检验。缩写：CI，置信区间；HR，风险比。"),
}
TITLE_ZH = {1: "表 1 按平均每天排便次数分组的参与者特征", 2: "表 2 大便性状与新发帕金森病和痴呆的关联"}


def to_md(d, lang, indent_first_col):
    tr = (lambda v: ZH.get(v, v).replace(" to ", " 至 ")) if lang == "zh" else (lambda v: v)
    head = "| " + " | ".join(tr(c) for c in d.columns) + " |"
    sep = "| " + " | ".join("---" for _ in d.columns) + " |"
    rows = []
    for r in d.itertuples(index=False):
        cells = [tr(str(v)) for v in r]
        if indent_first_col and r[0] in SUB:
            cells[0] = EM + cells[0]
        rows.append("| " + " | ".join(cells) + " |")
    return "\n".join([head, sep] + rows)


for lang, out, titles, notes, hdr in [("en", "src_en/tables.md", TITLE_EN, NOTE_EN, "# Tables"),
                                      ("zh", "src/表格.md", TITLE_ZH, NOTE_ZH, "# 表格")]:
    parts = [hdr]
    for k, d in [(1, t1), (2, t2)]:
        parts += ["\\pagebreak", f"**{titles[k]}**", to_md(d, lang, k == 1), notes[k]]
    open(P + out, "w", encoding="utf-8").write("\n\n".join(parts) + "\n")
    print("written", out)
