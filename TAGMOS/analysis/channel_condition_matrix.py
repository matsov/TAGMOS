#!/usr/bin/env python3
"""Channel-by-condition effect-size matrix and the transdiagnostic ordering.

For each axis in the panel and each condition, the covariate-adjusted
within-study effect (age, sex, body mass, log read depth and functional
richness), taken as the median across the studies contributing that condition.
Conditions are then ordered by the difference between the mean level of their
danger axes and the mean level of their protective axes, the valences being
those declared in the panel file.

Writes cMD_disease_channel_matrix_V42_CORRECTED.tsv and
cMD_transdiagnostic_ordering_CORRECTED.tsv.
"""
import tagmos_io
import pandas as pd, numpy as np, json
from scipy import stats
_a  = tagmos_io.cli("Channel-by-condition effect-size matrix and transdiagnostic ordering.",
              ec=True, panel=True)
OUT = _a.out
ECF, META = _a.ec, _a.meta
AX  = tagmos_io.load_panel(_a.panel)
AXES={n:v["ecs"] for n,v in AX.items() if v.get("ecs")}
# valences come from the panel file, not from this script
VAL={n:v.get("valence") for n,v in AX.items() if v.get("ecs")}
allec=sorted({e for ecs in AXES.values() for e in ecs})
hdr=pd.read_csv(ECF,sep="\t",nrows=0).columns.tolist(); use=[e for e in allec if e in hdr]
ec=pd.read_csv(ECF,sep="\t",usecols=["sample_id"]+use,low_memory=False).drop_duplicates("sample_id").set_index("sample_id")
for e in use: ec[e]=pd.to_numeric(ec[e],errors="coerce").fillna(0.0)
rich={}
for ch in pd.read_csv(ECF,sep="\t",low_memory=False,chunksize=1500):
    ch=ch.set_index("sample_id"); rich.update(((ch.apply(pd.to_numeric,errors="coerce").fillna(0)>0).sum(1)).to_dict())
ec["n_ec"]=pd.Series(rich).reindex(ec.index)
def z(s): s=pd.to_numeric(s,errors="coerce"); return (s-s.mean())/(s.std()+1e-12)
CH={n:pd.concat([z(ec[e]) for e in ecs if e in ec.columns],axis=1).mean(1) for n,ecs in AXES.items() if any(e in ec.columns for e in ecs)}
chan=pd.DataFrame(CH)
meta=pd.read_csv(META,sep="\t",dtype=str,low_memory=False).drop_duplicates("sample_id").set_index("sample_id")
X=chan.join(ec[["n_ec"]]).join(meta[["study_name","study_condition","disease_subtype","age","age_category","gender","BMI","antibiotics_current_use","number_reads","body_site"]],how="inner")
X["age_n"]=pd.to_numeric(X.age,errors="coerce");X["BMI_n"]=pd.to_numeric(X.BMI,errors="coerce")
X["logdepth"]=np.log10(pd.to_numeric(X.number_reads,errors="coerce").clip(lower=1)); X["logrich"]=np.log10(X.n_ec.clip(lower=1))
X["sex01"]=X.gender.map({"male":1,"female":0});X["abx01"]=X.antibiotics_current_use.map({"yes":1,"no":0})
X=X[X.body_site.fillna("stool")=="stool"]; X=X[~((X.age_category.isin(["newborn","child","schoolage"]))|(X.age_n<18)).fillna(False)]
CTRL={"control","healthy"};COV=["age_n","BMI_n","sex01","abx01","logdepth","logrich"]  # + logrich (ERRATA #3)
def eff(sub,case,col):
    u=[c for c in COV if sub[c].notna().mean()>=0.8 and sub[c].nunique(dropna=True)>1]
    cols={"p":z(sub[col]).values};[cols.update({c:z(sub[c]).values}) for c in u]
    Dm=pd.DataFrame(cols,index=sub.index);ok=Dm.notna().all(1)&np.isfinite(case);Dm=Dm[ok];y=case[ok.values]
    if len(y)<20 or y.sum()<6 or (1-y).sum()<6:return None
    Xm=np.column_stack([np.ones(len(y))]+[Dm[c].values for c in Dm.columns if np.std(Dm[c])>0])
    b,_,rk,_=np.linalg.lstsq(Xm,y,rcond=None); return b[1] if rk==Xm.shape[1] else None
CONF=["CRC","IBD","CD","UC","T2D","cirrhosis"]
recs={}
for lab in CONF:
    col="disease_subtype" if lab in ("CD","UC") else "study_condition"
    vec={}
    for ch in AXES:
        bs=[]
        for st,g in X.groupby("study_name"):
            ct=g[g.study_condition.isin(CTRL)]; ids=g.index[g[col]==lab]
            if len(ct)<8 or len(ids)<8: continue
            sub=pd.concat([g.loc[ids],ct]); case=np.r_[np.ones(len(ids)),np.zeros(len(ct))].astype(float)
            b=eff(sub,case,ch)
            if b is not None: bs.append(b)
        vec[ch]=np.median(bs) if bs else np.nan
    recs[lab]=vec
M=pd.DataFrame(recs).round(3); M["valence"]=[VAL[a] for a in M.index]
M.to_csv(f"{OUT}/cMD_disease_channel_matrix_V42_CORRECTED.tsv",sep="\t")
# transdiagnostic ordering with corrected valences: danger vs protective (exclude flags/context/research)
dang=[a for a in M.index if M.loc[a,"valence"]=="danger"]
prot=[a for a in M.index if M.loc[a,"valence"]=="protective"]
rows=[]
for p in CONF:
    dm=M.loc[dang,p].astype(float).mean(); pm=M.loc[prot,p].astype(float).mean()
    rows.append(dict(phenotype=p,danger_mean=round(dm,3),protective_mean=round(pm,3),transdiag_score=round(dm-pm,3)))
O=pd.DataFrame(rows).sort_values("transdiag_score",ascending=False)
O.to_csv(f"{OUT}/cMD_transdiagnostic_ordering_CORRECTED.tsv",sep="\t",index=False)
print("\ndanger axes(%d)=%s"%(len(dang),dang)); print("protective axes(%d)=%s"%(len(prot),prot))
print("\n=== transdiagnostic ordering ===")
print(O.to_string(index=False))
# compare vs old
try:
    old=pd.read_csv(f"{OUT}/cMD_transdiagnostic_ordering_V42.tsv",sep="\t").set_index("phenotype")["transdiag_score"]
    print("\ndelta versus a previous ordering in the same output directory:")
    for _,r in O.iterrows(): print("  %-10s new=%.3f old=%.3f"%(r.phenotype,r.transdiag_score,old.get(r.phenotype,np.nan)))
except Exception as e: print("old compare:",e)
