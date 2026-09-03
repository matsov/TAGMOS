#!/usr/bin/env python3
"""Evidence grade for every channel-by-condition association.

For each axis in the panel and each condition of the compendium: odds ratio per
standard deviation from a within-study logistic model adjusted for age, sex,
body mass, log read depth and functional richness; a DerSimonian-Laird
random-effects pool across studies; the number of independent studies; and an
evidence grade (STRONG when at least three studies agree at q < 0.05, then
ASSOCIATION, WEAK, and MECHANISM for a single study).

Every comparison is made within its own study and never pooled at the sample
level. Writes channel_evidence_grade_cMD.tsv.
"""
import tagmos_io
import pandas as pd, numpy as np, json
from datetime import datetime
from scipy import stats
def log(*a): print(f"[{datetime.now():%H:%M:%S}]",*a,flush=True)
_a  = tagmos_io.cli("Evidence grade per channel and condition.", ec=True, panel=True, richness=True)
OUT = _a.out
AX  = tagmos_io.load_panel(_a.panel)
ec   = tagmos_io.load_ec(_a.ec)
nec  = tagmos_io.load_richness(_a.richness, ec)
meta = tagmos_io.load_meta(_a.meta)
AXES={n:v["ecs"] for n,v in AX.items() if v.get("ecs")}

log("Loading inputs...")
for c in ec.columns:
    if c!="study": ec[c]=pd.to_numeric(ec[c],errors="coerce").fillna(0.0)*1e6
X=ec.join(meta[["study_name","study_condition","disease_subtype","age","age_category","gender","BMI","number_reads","body_site"]],how="inner")
X["age_n"]=pd.to_numeric(X.age,errors="coerce");X["BMI_n"]=pd.to_numeric(X.BMI,errors="coerce")
X["sex01"]=X.gender.map({"male":1,"female":0});X["logreads"]=np.log10(pd.to_numeric(X.number_reads,errors="coerce").clip(lower=1))
X["lognec"]=np.log10(nec.reindex(X.index).clip(lower=1))
X=X[X.body_site.fillna("stool")=="stool"]; X=X[~((X.age_category.isin(["newborn","child","schoolage"]))|(X.age_n<18)).fillna(False)]
def z(s): s=pd.to_numeric(s,errors="coerce"); sd=s.std(); return (s-s.mean())/(sd if sd>0 else 1)
# axis value = mean of z(log10(EC_cpm + 0.1)) over the enzymes present
for ch,ecs in AXES.items():
    cols=[e for e in ecs if e in ec.columns]
    if not cols: continue
    X[f"CH_{ch}"]=pd.concat([z(np.log10(X[e]+0.1)) for e in cols],axis=1).mean(1)
CH=[c for c in X.columns if c.startswith("CH_")]
COV=["age_n","BMI_n","sex01","logreads","lognec"]; CTRL={"control","healthy"}
CONDS={"CRC":("study_condition","CRC"),"adenoma":("study_condition","adenoma"),"IBD":("study_condition","IBD"),
 "UC":("disease_subtype","UC"),"CD":("disease_subtype","CD"),"cirrhosis":("study_condition","cirrhosis"),
 "STH":("study_condition","STH"),"T2D":("study_condition","T2D"),"IGT":("study_condition","IGT"),
 "ACVD":("study_condition","ACVD"),"hypertension":("study_condition","hypertension"),
 "schizophrenia":("study_condition","schizophrenia"),"ME/CFS":("study_condition","ME/CFS")}
def logit_beta(g,ic,ict,xcol):
    order=list(ic)+list(ict); u=[c for c in COV if g.loc[order,c].notna().mean()>=0.8 and g.loc[order,c].nunique(dropna=True)>1]
    yv=np.r_[np.ones(len(ic)),np.zeros(len(ict))]
    M=np.column_stack([np.ones(len(order)),z(g.loc[order,xcol]).values]+[z(g.loc[order,c]).values for c in u])
    ok=np.all(np.isfinite(M),1); M,yv=M[ok],yv[ok]
    if len(yv)<20 or yv.sum()<6 or (len(yv)-yv.sum())<6: return None
    b=np.zeros(M.shape[1])
    for _ in range(60):
        p=1/(1+np.exp(-M@b));W=np.clip(p*(1-p),1e-6,None);zz=M@b+(yv-p)/W
        try:bn=np.linalg.solve((M*W[:,None]).T@M,(M*W[:,None]).T@zz)
        except np.linalg.LinAlgError:return None
        if np.max(np.abs(bn-b))<1e-7:b=bn;break
        b=bn
    p=1/(1+np.exp(-M@b));W=np.clip(p*(1-p),1e-6,None)
    try:cov=np.linalg.inv((M*W[:,None]).T@M)
    except np.linalg.LinAlgError:return None
    return b[1],np.sqrt(cov[1,1])
def dl(bs):
    b=np.array([x[0] for x in bs]);s=np.clip(np.array([x[1] for x in bs]),1e-4,None);v=s**2;w=1/v;k=len(b)
    m=np.sum(w*b)/np.sum(w);Q=np.sum(w*(b-m)**2);den=np.sum(w)-np.sum(w**2)/np.sum(w)
    t2=max(0,(Q-(k-1))/den) if (k>1 and den>0) else 0;wr=1/(v+t2);mr=np.sum(wr*b)/np.sum(wr);se=np.sqrt(1/np.sum(wr))
    return mr,2*stats.norm.sf(abs(mr/se)),k
NSTRONG={"CRC":11,"adenoma":5}
rows=[]
for cond,(col,val) in CONDS.items():
    S=[]
    for st,g in X.groupby("study_name"):
        c=g.index[g[col]==val]; ct=g.index[g.study_condition.isin(CTRL)]
        if len(c)>=6 and len(ct)>=6: S.append((c,ct,g))
    for ch in CH:
        bs=[logit_beta(g,c,ct,ch) for c,ct,g in S]; bs=[b for b in bs if b]
        if len(bs)<1: continue
        mr,p,k=dl(bs); OR=np.exp(mr)
        rows.append(dict(channel=ch[3:],condition=cond,OR=round(OR,3),p=round(p,4),n_studies=k))
M=pd.DataFrame(rows)
M["q"]=np.nan
for cond,sub in M.groupby("condition"):
    m=sub.p.between(0,1); M.loc[sub.index[m],"q"]=stats.false_discovery_control(sub.loc[m,"p"].values,method="bh")
def grade(r):
    if r.n_studies>=3 and r.q<0.05: return "FORTE"
    if (r.n_studies>=3 and r.q<0.10) or (r.n_studies==2 and r.p<0.05): return "ASSOCIAZIONE"
    if r.n_studies==1 and r.p<0.05: return "MECCANISMO"
    return "DEBOLE"
M["grade"]=M.apply(grade,axis=1)
FEATURED={"ENG","SACCHAROLYTIC","RIBOFLAVIN","TMA","PROT_stickland","LPS_hexaacyl"}
CONTEXT={"POLYAMINES","IRON_pathobiont","IRON","BA_deconj","BA_secondary","MUC_z","MUC","HIS","UREM_pcresol","UREM_tyramine","ETU","LACTATE"}
def roster(ch): return "featured_tierA" if ch in FEATURED else ("contesto/descrittivo" if ch in CONTEXT else "altro")
M["roster"]=M.channel.map(roster)
M=M[["channel","condition","OR","p","q","n_studies","grade","roster"]].sort_values(["condition","q"])
M.to_csv(f"{OUT}/channel_evidence_grade_cMD.tsv",sep="\t",index=False)
log("rows:",len(M),"| FORTE:",int((M.grade=='FORTE').sum()))
print("\n=== STRONG channels (>=3 studies, q<0.05) ===")
print(M[M.grade=='FORTE'].sort_values(["channel","condition"])[["channel","condition","OR","q","n_studies","roster"]].to_string(index=False))
print("\n=== per-condition summary ===")
feat=M[(M.channel.isin(FEATURED))&(M.condition.isin(["CRC","adenoma","IBD","cirrhosis","T2D","STH"]))]
print(feat.pivot_table(index="channel",columns="condition",values="OR").round(2).to_string())
