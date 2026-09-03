#!/usr/bin/env python3
"""A taxonomic guild as a negative control against random species panels.

Scores a species guild supplied with --guild across the colorectal-cancer
studies of the compendium: per-study odds ratio from a within-study logistic
model adjusted for age, sex, body mass, log read depth and species richness, a
random-effects pool, and leave-one-study-out. The control is 150 null panels of
the same size drawn at matched detection prevalence, giving a pooled null odds
ratio and an empirical P.

Writes butyrate_tax_CRC_negcontrol.tsv.
"""
import tagmos_io
import pandas as pd, numpy as np
from datetime import datetime
from scipy import stats
def log(*a): print(f"[{datetime.now():%H:%M:%S}]",*a,flush=True)
_a  = tagmos_io.cli("Taxonomic guild as a negative control against prevalence-matched random species panels.",
              species=True, guild=True)
CONDITION = _a.condition
OUT = _a.out
SPF = _a.species
gr  = pd.read_csv(_a.guild, sep="\t")
meta = tagmos_io.load_meta(_a.meta)
rng=np.random.default_rng(20260712)
# studies contributing the target condition, derived from the metadata rather than
# fixed in the source: any study with at least MIN_N cases and MIN_N controls.
MIN_N=8
_c=meta[meta["study_condition"].isin([CONDITION,"control"])]
CRC_STUDIES=sorted(st for st,g in _c.groupby("study_name")
                   if (g["study_condition"]==CONDITION).sum()>=MIN_N
                   and (g["study_condition"]=="control").sum()>=MIN_N)
log(f"{CONDITION}: {len(CRC_STUDIES)} studies with >={MIN_N} cases and >={MIN_N} controls")
# guild panel supplied by the user; Akkermansia is excluded as in the paper
BUT=[s for s in gr.iloc[:,0].astype(str) if "akkermansia" not in s.lower() and s!="feature"]
log("guild panel n=",len(BUT))

hdr=pd.read_csv(SPF,sep="\t",nrows=0).columns.tolist()
crc_cols=[c for c in hdr[1:] if c.split("|")[0] in CRC_STUDIES]
log("sample columns from those studies:",len(crc_cols))
sp=pd.read_csv(SPF,sep="\t",usecols=["feature"]+crc_cols,low_memory=False).set_index("feature")
sp=sp.apply(pd.to_numeric,errors="coerce").fillna(0.0)
S=sp.T  # samples x species
S.index=[c.split("|")[-1] for c in S.index]   # -> sample_id
M=meta.reindex(S.index)
df=pd.DataFrame(index=S.index)
df["study"]=M.study_name.values; df["cond"]=M.study_condition.values
df["age"]=pd.to_numeric(M.age,errors="coerce").values; df["bmi"]=pd.to_numeric(M.BMI,errors="coerce").values
df["sex01"]=M.gender.map({"male":1,"female":0}).values
df["logreads"]=np.log10(pd.to_numeric(M.number_reads,errors="coerce").clip(lower=1)).values
df["srich"]=np.log10((S>0).sum(1).clip(lower=1)).values
df=df[df.cond.isin(["CRC","control","healthy"])]
S=S.reindex(df.index)
butcols=[s for s in BUT if s in S.columns]
log("BUT specie presenti nel master:",len(butcols))
df["BUT_tax"]=S[butcols].sum(1)

def z(s): s=pd.to_numeric(s,errors="coerce"); sd=s.std(); return (s-s.mean())/(sd if sd>0 else 1)
COV=["age","bmi","sex01","logreads","srich"]
def study_OR(dd,score):
    dd=dd.copy(); dd["_s"]=score.reindex(dd.index)
    outs=[]
    for st,g in dd.groupby("study"):
        y=(g.cond=="CRC").astype(float).values
        if y.sum()<6 or (len(y)-y.sum())<6: continue
        u=[c for c in COV if g[c].notna().mean()>=0.8 and g[c].nunique()>1]
        X=np.column_stack([np.ones(len(g)),z(g["_s"]).values]+[z(g[c]).values for c in u])
        ok=np.all(np.isfinite(X),1); X2,y2=X[ok],y[ok]
        if len(y2)<20: continue
        b=np.zeros(X2.shape[1])
        for _ in range(50):
            p=1/(1+np.exp(-X2@b));W=np.clip(p*(1-p),1e-6,None);zz=X2@b+(y2-p)/W
            try:bn=np.linalg.solve((X2*W[:,None]).T@X2,(X2*W[:,None]).T@zz)
            except np.linalg.LinAlgError:bn=None;break
            if bn is None or np.max(np.abs(bn-b))<1e-7:b=bn if bn is not None else b;break
            b=bn
        if b is None: continue
        p=1/(1+np.exp(-X2@b));W=np.clip(p*(1-p),1e-6,None)
        try:cov=np.linalg.inv((X2*W[:,None]).T@X2)
        except np.linalg.LinAlgError:continue
        outs.append((st,b[1],np.sqrt(cov[1,1])))
    return outs
def dl(bs):
    b=np.array([x[1] for x in bs]);s=np.clip(np.array([x[2] for x in bs]),1e-4,None);v=s**2;w=1/v;k=len(b)
    m=np.sum(w*b)/np.sum(w);Q=np.sum(w*(b-m)**2);den=np.sum(w)-np.sum(w**2)/np.sum(w)
    t2=max(0,(Q-(k-1))/den) if (k>1 and den>0) else 0;wr=1/(v+t2);mr=np.sum(wr*b)/np.sum(wr);se=np.sqrt(1/np.sum(wr))
    return mr,se,2*stats.norm.sf(abs(mr/se)),k

log("=== per-study OR (guild panel) ===")
bs=study_OR(df,df["BUT_tax"]);
perstudy=[(st,round(np.exp(b),3)) for st,b,se in bs]
mr,se,p,k=dl(bs); pooled_OR=np.exp(mr)
log("pooled OR=%.3f p=%.2g (%d studies)"%(pooled_OR,p,k))
# LOSO
loso=[]
for drop in set(x[0] for x in bs):
    sub=[x for x in bs if x[0]!=drop]; m2,_,_,_=dl(sub); loso.append(np.exp(m2))
loso_lo,loso_hi=min(loso),max(loso)
log("LOSO OR range: %.3f - %.3f"%(loso_lo,loso_hi))

# === CONTROLLO NEGATIVO: 150 pannelli nulli di 40 specie appaiate per prevalenza ===
log("=== negative control (150 null panels, prevalence-matched species) ===")
prev=(S>0).mean(0)
nonbut=[c for c in S.columns if c not in butcols and prev[c]>0.02]
but_prev=prev[butcols].values
bins=np.quantile(prev[nonbut],np.linspace(0,1,11))
def match_panel():
    chosen=[]
    for pv in rng.choice(but_prev,40,replace=True):
        # candidati con prevalenza simile
        lo,hi=pv*0.5,pv*1.5+0.02
        cand=[c for c in nonbut if lo<=prev[c]<=hi]
        if not cand: cand=nonbut
        chosen.append(rng.choice(cand))
    return list(dict.fromkeys(chosen))
null_ORs=[]
for i in range(150):
    panel=match_panel()
    score=S[panel].sum(1)
    nb=study_OR(df,score)
    if len(nb)>=3:
        m2,_,_,_=dl(nb); null_ORs.append(np.exp(m2))
null_ORs=np.array(null_ORs); null_med=np.median(null_ORs)
# p empirico: quanti pannelli nulli altrettanto o piu protettivi dell'osservato
p_emp=(np.sum(null_ORs<=pooled_OR)+1)/(len(null_ORs)+1)
log("null median OR=%.3f | n_null=%d | p_empirical=%.4f"%(null_med,len(null_ORs),p_emp))

rows=[dict(study=st,OR=orr) for st,orr in perstudy]
rows.append(dict(study="__POOLED__",OR=round(pooled_OR,3)))
rows.append(dict(study="__pooled_p__",OR=round(p,5)))
rows.append(dict(study="__LOSO_range__",OR=f"{loso_lo:.3f}-{loso_hi:.3f}"))
rows.append(dict(study="__null_median_OR__",OR=round(null_med,3)))
rows.append(dict(study="__p_empirical__",OR=round(p_emp,4)))
rows.append(dict(study="__n_studies__",OR=k))
rows.append(dict(study="__n_BUT_species__",OR=len(butcols)))
pd.DataFrame(rows).to_csv(f"{OUT}/butyrate_tax_CRC_negcontrol.tsv",sep="\t",index=False)
print("\n=== per-study OR ==="); [print(f"  {st:20s} {orr}") for st,orr in perstudy]
print(f"\nPOOLED OR={pooled_OR:.3f} (p={p:.2g}, {k} studies) | LOSO {loso_lo:.3f}-{loso_hi:.3f} | null med {null_med:.3f} | p_emp {p_emp:.4f}")
log("DONE.")
