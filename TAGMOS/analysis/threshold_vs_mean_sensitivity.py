#!/usr/bin/env python3
"""Sensitivity of the threshold-only fraction to the gate and the FDR cutoff.

Repeats the analysis of threshold_vs_mean.py across gate thresholds
(z in 1.0, 1.5, 2.0 against the study's own controls) and false-discovery
cutoffs (q in 0.05, 0.10), to show that the fraction of associations detectable
only through the threshold is not an artefact of either choice. The mean
reading is threshold-independent and is computed once.

Writes tail_vs_mean_sensitivity_cMD.tsv and its summary grid.
"""
import tagmos_io
import pandas as pd, numpy as np, json
from datetime import datetime
from scipy import stats
def log(*a): print(f"[{datetime.now():%H:%M:%S}]",*a,flush=True)
_a  = tagmos_io.cli("Sensitivity of the threshold-only fraction to gate and FDR cutoffs.", ec=True, panel=True, richness=True)
OUT = _a.out
AX  = tagmos_io.load_panel(_a.panel)
ec   = tagmos_io.load_ec(_a.ec)
nec  = tagmos_io.load_richness(_a.richness, ec)
meta = tagmos_io.load_meta(_a.meta)
AXES={n:v for n,v in AX.items() if v.get("ecs")}
for c in ec.columns:
    if c!="study": ec[c]=pd.to_numeric(ec[c],errors="coerce").fillna(0.0)*1e6
X=ec.join(meta[["study_name","study_condition","disease_subtype","age","age_category","gender","BMI","number_reads","body_site"]],how="inner")
X["age_n"]=pd.to_numeric(X.age,errors="coerce");X["BMI_n"]=pd.to_numeric(X.BMI,errors="coerce")
X["sex01"]=X.gender.map({"male":1,"female":0});X["logreads"]=np.log10(pd.to_numeric(X.number_reads,errors="coerce").clip(lower=1))
X["logrich"]=np.log10(nec.reindex(X.index).clip(lower=1))
X=X[X.body_site.fillna("stool")=="stool"];X=X[~((X.age_category.isin(["newborn","child","schoolage"]))|(X.age_n<18)).fillna(False)]
def z(s): s=pd.to_numeric(s,errors="coerce");sd=s.std();return (s-s.mean())/(sd if sd>0 else 1)
CH={};VAL={}
for ch,v in AXES.items():
    cols=[e for e in v["ecs"] if e in ec.columns]
    if cols: CH[ch]=pd.concat([np.log10(X[e]+0.1) for e in cols],axis=1).mean(1);VAL[ch]=v.get("valence","danger")
CHdf=pd.DataFrame(CH)
CTRL={"control","healthy"};COV=["age_n","BMI_n","sex01","logreads","logrich"]
CONDS={"CRC":("study_condition","CRC"),"adenoma":("study_condition","adenoma"),"IBD":("study_condition","IBD"),
 "UC":("disease_subtype","UC"),"CD":("disease_subtype","CD"),"cirrhosis":("study_condition","cirrhosis"),
 "STH":("study_condition","STH"),"T2D":("study_condition","T2D"),"IGT":("study_condition","IGT"),
 "ACVD":("study_condition","ACVD"),"hypertension":("study_condition","hypertension")}
THRS=[1.0,1.5,2.0]
def logit(y,cols):
    M=np.column_stack([np.ones(len(y))]+cols);ok=np.all(np.isfinite(M),1)&np.isfinite(y);M,y=M[ok],y[ok]
    if len(y)<20 or y.sum()<5 or (len(y)-y.sum())<5 or np.std(M[:,1])==0: return None
    b=np.zeros(M.shape[1])
    for _ in range(60):
        p=1/(1+np.exp(-np.clip(M@b,-30,30)));W=np.clip(p*(1-p),1e-6,None)
        try: bn=np.linalg.solve((M*W[:,None]).T@M,(M*W[:,None]).T@(M@b+(y-p)/W))
        except np.linalg.LinAlgError: return None
        if np.max(np.abs(bn-b))<1e-7: b=bn;break
        b=bn
    p=1/(1+np.exp(-np.clip(M@b,-30,30)));W=np.clip(p*(1-p),1e-6,None)
    try: cov=np.linalg.inv((M*W[:,None]).T@M)
    except np.linalg.LinAlgError: return None
    return b[1],np.sqrt(cov[1,1])
def dl(bs):
    b=np.array([x[0] for x in bs]);s=np.clip(np.array([x[1] for x in bs]),1e-4,None);v=s**2;w=1/v;k=len(b)
    m=np.sum(w*b)/np.sum(w);Q=np.sum(w*(b-m)**2);den=np.sum(w)-np.sum(w**2)/np.sum(w)
    t2=max(0,(Q-(k-1))/den) if (k>1 and den>0) else 0;wr=1/(v+t2);mr=np.sum(wr*b)/np.sum(wr);se=np.sqrt(1/np.sum(wr))
    return mr,2*stats.norm.sf(abs(mr/se)),k
log("Sweep: mean reading plus threshold reading at three gates...")
rows=[]
for cond,(col,val) in CONDS.items():
    for ch in CH:
        prot=VAL[ch] in ("protective","protective_redox_flag")
        mb=[];tb={t:[] for t in THRS};cM=[];cT={t:[] for t in THRS}
        for st,g in X.groupby("study_name"):
            case=g.index[g[col]==val];ctrl=g.index[g.study_condition.isin(CTRL)]
            if len(case)<5 or len(ctrl)<8: continue
            cv=CHdf.loc[g.index,ch];mu=cv.loc[ctrl].mean();sd=cv.loc[ctrl].std()
            if not sd or sd==0: continue
            zc=(cv-mu)/sd;order=list(case)+list(ctrl);y=np.r_[np.ones(len(case)),np.zeros(len(ctrl))]
            covs=[z(g.loc[order,c]).values for c in COV if g.loc[order,c].notna().mean()>=0.8 and g.loc[order,c].nunique()>1]
            rm=logit(y,[zc.loc[order].values]+covs)
            if rm: mb.append(rm);cM.append(np.sign(rm[0]))
            for t in THRS:
                gate=(zc>=t).astype(float) if not prot else (zc<=-t).astype(float)
                gv=gate.loc[order].values
                if 0<gv.sum()<len(gv):
                    rt=logit(y,[gv]+covs)
                    if rt: tb[t].append(rt);cT[t].append(np.sign(rt[0]))
        rec=dict(condition=cond,channel=ch,valence=VAL[ch])
        if len(mb)>=2: mm,mp,mk=dl(mb);rec.update(mean_OR=round(np.exp(mm),3),mean_p=mp,mean_concord=round(np.mean(np.array(cM)==np.sign(mm)),2))
        for t in THRS:
            if len(tb[t])>=2: tm,tp,tk=dl(tb[t]);rec[f"tail_p_{t}"]=tp;rec[f"tail_OR_{t}"]=round(np.exp(tm),3);rec[f"tail_concord_{t}"]=round(np.mean(np.array(cT[t])==np.sign(tm)),2)
        rows.append(rec)
R=pd.DataFrame(rows)
# FDR
R["mean_q"]=np.nan;m=R.mean_p.between(0,1);R.loc[m,"mean_q"]=stats.false_discovery_control(R.loc[m,"mean_p"].values,method="bh")
for t in THRS:
    col=f"tail_p_{t}";R[f"tail_q_{t}"]=np.nan;m=R[col].between(0,1);R.loc[m,f"tail_q_{t}"]=stats.false_discovery_control(R.loc[m,col].values,method="bh")
R.to_csv(f"{OUT}/tail_vs_mean_sensitivity_cMD.tsv",sep="\t",index=False)

print("\n=== SENSITIVITY GRID: threshold-only fraction among detectable associations ===")
print(f"{'soglia':>7} {'q':>6} {'detect':>7} {'BOTH':>5} {'MEANonly':>9} {'TAILonly':>9} {'%TAILonly':>10} {'concT':>6} {'concM':>6}")
grid=[]
for t in THRS:
    for q in [0.05,0.10]:
        ms=R.mean_q<q; ts=R[f"tail_q_{t}"]<q
        both=(ms&ts).sum();mo=(ms&~ts).sum();to=(ts&~ms).sum();det=both+mo+to
        cT=R.loc[ts,f"tail_concord_{t}"].dropna().mean();cM=R.loc[ms,"mean_concord"].dropna().mean()
        pct=100*to/max(det,1)
        print(f"{t:7.1f} {q:6.2f} {det:7d} {both:5d} {mo:9d} {to:9d} {pct:9.0f}% {cT:6.2f} {cM:6.2f}")
        grid.append(dict(gate_threshold=t,q_cutoff=q,detectable=det,both=both,mean_only=mo,tail_only=to,pct_tail_only=round(pct,1),concord_tail=round(cT,2),concord_mean=round(cM,2)))
pd.DataFrame(grid).to_csv(f"{OUT}/tail_vs_mean_sensitivity_grid_cMD.tsv",sep="\t",index=False)
log("DONE.")
