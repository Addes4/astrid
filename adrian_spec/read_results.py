#!/usr/bin/env python
"""First-pass readout of an ASTRID result CSV.

Usage: python adrian_spec/read_results.py adrian_spec/test_run/p2_ASTRID_Result.csv
"""
import sys, re, collections
import pandas as pd

f = sys.argv[1] if len(sys.argv) > 1 else "adrian_spec/test_run/p2_ASTRID_Result.csv"
d = pd.read_csv(f)
key = d.columns[0]                      # the final clustering level
n = d.CellCount.sum()
print(f"{f}\n{len(d)} clusters, {n} cells, final level = {key}\n")

# 1. composition, weighted by cluster size
g = d.groupby("SingleR_CellType").agg(clusters=("CellCount", "size"), cells=("CellCount", "sum"))
g["pct"] = (100 * g.cells / n).round(1)
print("COMPOSITION (SingleR call)")
print(g.sort_values("cells", ascending=False).to_string(), "\n")

# 2. did the marker logic confirm the call?
passed = d.FormulaPassed.astype(str) == "True"
na = d.FormulaPassed.astype(str) == "not applicable"
print(f"VALIDATION  confirmed {passed.sum()}/{len(d)} clusters "
      f"({d.loc[passed,'CellCount'].sum()}/{n} cells, "
      f"{100*d.loc[passed,'CellCount'].sum()/n:.1f}%)"
      + (f"; {na.sum()} clusters have no formula" if na.sum() else ""))

# 3. contradictions: call not among the types the markers support
def other(r):
    return [] if pd.isna(r.OtherCellTypesPassed) else str(r.OtherCellTypesPassed).split(";")
d["contradicted"] = d.apply(lambda r: (not pd.isna(r.OtherCellTypesPassed)
                                       and str(r.OtherCellTypesPassed) != ""
                                       and r.SingleR_CellType not in other(r)), axis=1)
con = d[d.contradicted].sort_values("CellCount", ascending=False)
print(f"CONTRADICTED  {len(con)} clusters / {con.CellCount.sum()} cells "
      f"where markers support something other than the SingleR call")
if len(con):
    print(con.head(8)[[key, "SingleR_CellType", "CellCount",
                       "OtherCellTypesPassed", "TopOddsRatioGenes"]].to_string(index=False))
print()

# 4. unassigned
un = d[d[key].astype(str) == "unassigned"]
if len(un):
    c = int(un.CellCount.iloc[0])
    print(f"UNASSIGNED  {c} cells ({100*c/n:.1f}%) never reached the final level "
          f"(skipped by validation, ASTRID_v0.01.py:446)\n")

# 5. CNV — relative spread only; absolute values are not comparable across runs
print("CNV (compare clusters within this run, not across runs)")
big = d[d.CellCount >= 20]
print(big[["newCNVScoreABS", "newCNVScoreSQR", "CNVCorrelation"]].describe()
        .loc[["min", "50%", "max"]].round(4).to_string())
top = big.nlargest(5, "newCNVScoreABS")
print("\nmost CNV-damaged clusters:")
print(top[[key, "SingleR_CellType", "CellCount", "newCNVScoreABS", "CNVCorrelation"]].to_string(index=False))

# 6. ground truth, if author labels were supplied
if "paperCellTypeAbundance" in d and d.paperCellTypeAbundance.notna().any():
    tot = collections.Counter()
    for _, r in d.iterrows():
        for m in re.finditer(r"([^;()]+)\(([\d.]+)%\)", str(r.paperCellTypeAbundance)):
            tot[m.group(1).strip()] += r.CellCount * float(m.group(2)) / 100
    t = pd.Series(tot).sort_values(ascending=False)
    print("\nAUTHOR LABELS (approx cells, top 10)")
    print(t.head(10).round(0).to_string())

if "cancer_percentage" in d and d.cancer_percentage.max() == 0:
    print("\nNOTE: cancer_percentage is all zero. This is the AUTHOR label matching the "
          "string 'Cancer' (ASTRID_v0.01.py:399), not an ASTRID prediction.")
