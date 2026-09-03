#!/usr/bin/env python3
"""Threshold versus mean detection, for every channel-by-condition association.

The claim this script tests is that a systematic class of disease associations
is visible only as a threshold and is missed entirely by a comparison of means.

For each axis in the panel and each condition, two readings of the same data are
computed within study and adjusted for age, sex, body mass, log read depth and
functional richness, pooled by random effects and corrected by Benjamini-Hochberg:

  MEAN       odds ratio per standard deviation of the axis, z-scored against the
             study's own controls -- the "compare the means" reading;
  THRESHOLD  odds ratio for occupying the dysbiotic tail (z >= +1 against the
             study's controls for a danger axis, z <= -1 for a protective one)
             -- the "gate" reading.

Each cell is classified BOTH, TAIL_ONLY, MEAN_ONLY or neither, and the
cross-study sign concordance of the two readings is compared.

Writes tail_vs_mean_channel_condition_cMD.tsv.
"""
import tagmos_io
import pandas as pd, numpy as np, json
from datetime import datetime
from scipy import stats
def log(*a): print(f"[{datetime.now():%H:%M:%S}]",*a,flush=True)
_a  = tagmos_io.cli("Threshold versus mean detection, per channel and condition.", ec=True, panel=True, richness=True)
OUT = _a.out
AX  = tagmos_io.load_panel(_a.panel)
ec   = tagmos_io.load_ec(_a.ec)
nec  = tagmos_io.load_richness(_a.richness, ec)
meta = tagmos_io.load_meta(_a.meta)
AXES={n:v for n,v in AX.items() if v.get("ecs")}

log("Loading inputs...")
for c in ec.columns:
    if c!="study": ec[c]=pd.to_numeric(ec[c],errors="coerce").fillna(0.0)*1e6
X=ec.join(meta[["study_name","study_condition","disease_subtype","age","age_category","gender","BMI","number_reads","body_site"]],how="inner")
X["age_n"]=pd.to_numeric(X.age,errors="coerce");X["BMI_n"]=pd.to_numeric(X.BMI,errors="coerce")
X["sex01"]=X.gender.map({"male":1,"female":0});X["logreads"]=np.log10(pd.to_numeric(X.number_reads,errors="coerce").clip(lower=1))
X["logrich"]=np.log10(nec.reindex(X.index).clip(lower=1))
X=X[X.body_site.fillna("stool")=="stool"];X=X[~((X.age_category.isin(["newborn","child","schoolage"]))|(X.age_n<18)).fillna(False)]
def z(s): s=pd.to_numeric(s,errors="coerce");sd=s.std();return (s-s.mean())/(sd if sd>0 else 1)
# raw axis value = mean of log10(EC_cpm + 0.1) over the enzymes
CH={}; VAL={}
for ch,v in AXES.items():
    cols=[e for e in v["ecs"] if e in ec.columns]
    if cols: CH[ch]=pd.concat([np.log10(X[e]+0.1) for e in cols],axis=1).mean(1); VAL[ch]=v.get("valence","danger")
CHdf=pd.DataFrame(CH)
CTRL={"control","healthy"};COV=["age_n","BMI_n","sex01","logreads","logrich"]
CONDS={"CRC":("study_condition","CRC"),"adenoma":("study_condition","adenoma"),"IBD":("study_condition","IBD"),
 "UC":("disease_subtype","UC"),"CD":("disease_subtype","CD"),"cirrhosis":("study_condition","cirrhosis"),
 "STH":("study_condition","STH"),"T2D":("study_condition","T2D"),"IGT":("study_condition","IGT"),
 "ACVD":("study_condition","ACVD"),"hypertension":("study_condition","hypertension")}
def logit(y,cols):
    M=np.column_stack([np.ones(len(y))]+cols); ok=np.all(np.isfinite(M),1)&np.isfinite(y);M,y=M[ok],y[ok]
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
log("Sweep: mean versus threshold reading...")
rows=[]
for cond,(col,val) in CONDS.items():
    for ch in CH:
        prot = VAL[ch] in ("protective","protective_redox_flag")
        mb=[];tb=[];signM=[];signT=[]
        for st,g in X.groupby("study_name"):
            case=g.index[g[col]==val]; ctrl=g.index[g.study_condition.isin(CTRL)]
            if len(case)<5 or len(ctrl)<8: continue
            cv=CHdf.loc[g.index,ch]; mu=cv.loc[ctrl].mean(); sd=cv.loc[ctrl].std()
            if not sd or sd==0: continue
            zc=(cv-mu)/sd
            idx=case.union(ctrl); y=g.loc[idx,col].eq(val).astype(float).values if col=="disease_subtype" else np.r_[np.ones(len(case)),np.zeros(len(ctrl))]
            order=list(case)+list(ctrl); y=np.r_[np.ones(len(case)),np.zeros(len(ctrl))]
            covs=[z(g.loc[order,c]).values for c in COV if g.loc[order,c].notna().mean()>=0.8 and g.loc[order,c].nunique()>1]
            # MEDIA: z continuo
            rm=logit(y,[zc.loc[order].values]+covs)
            if rm: mb.append(rm); signM.append(np.sign(rm[0]))
            # CODA: gate binario riferito ai controlli
            gate=(zc>=1).astype(float) if not prot else (zc<=-1).astype(float)
            gv=gate.loc[order].values
            if 0<gv.sum()<len(gv):
                rt=logit(y,[gv]+covs)
                if rt: tb.append(rt); signT.append(np.sign(rt[0]))
        if len(mb)<2 and len(tb)<2: continue
        rec=dict(condition=cond,channel=ch,valence=VAL[ch])
        if len(mb)>=2:
            mm,mp,mk=dl(mb); rec.update(mean_OR=round(np.exp(mm),3),mean_p=mp,mean_nstudies=mk,mean_concord=round(np.mean(np.array(signM)==np.sign(mm)),2))
        if len(tb)>=2:
            tm,tp,tk=dl(tb); rec.update(tail_OR=round(np.exp(tm),3),tail_p=tp,tail_nstudies=tk,tail_concord=round(np.mean(np.array(signT)==np.sign(tm)),2))
        rows.append(rec)
R=pd.DataFrame(rows)
for c in ["mean_p","tail_p"]:
    R[c.replace("_p","_q")]=np.nan; m=R[c].between(0,1); R.loc[m,c.replace("_p","_q")]=stats.false_discovery_control(R.loc[m,c].values,method="bh")
def cls(r):
    ms=(r.get("mean_q",1)<0.10); ts=(r.get("tail_q",1)<0.10)
    return "BOTH" if ms and ts else ("TAIL_ONLY" if ts and not ms else ("MEAN_ONLY" if ms and not ts else "neither"))
R["class"]=R.apply(cls,axis=1)
R.to_csv(f"{OUT}/tail_vs_mean_channel_condition_cMD.tsv",sep="\t",index=False)
det=R[R["class"]!="neither"]
log("celle testate:",len(R),"| rilevabili (coda o media, q<0.10):",len(det))
print("\n=== classification of the detectable associations ===")
print(R["class"].value_counts().to_string())
n_tail=(R["class"]=="TAIL_ONLY").sum(); n_both=(R["class"]=="BOTH").sum(); n_mean=(R["class"]=="MEAN_ONLY").sum()
print(f"\nTHRESHOLD-ONLY FRACTION among detectable associations: {n_tail}/{len(det)} = {100*n_tail/max(len(det),1):.0f}%")
print(f"  (BOTH={n_both}, MEAN_ONLY={n_mean}, TAIL_ONLY={n_tail})")
print("\n=== reproducibility: cross-study sign concordance (mean where detectable) ===")
print(f"  concordanza CODA:  {det['tail_concord'].dropna().mean():.2f}")
print(f"  concordanza MEDIA: {det['mean_concord'].dropna().mean():.2f}")
print("\n=== threshold-only examples (readable as a gate, invisible to the mean) ===")
ex=R[R["class"]=="TAIL_ONLY"].sort_values("tail_q").head(15)
print(ex[["condition","channel","valence","mean_OR","mean_q","tail_OR","tail_q"]].round(3).to_string(index=False))
log("DONE.")
