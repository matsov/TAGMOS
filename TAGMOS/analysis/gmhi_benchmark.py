#!/usr/bin/env python3
"""
GMHI (Gupta et al., Nat Commun 2020) sul compendio cMD 93-studi,
passato per la stessa pipeline di validazione del composito curato:
  M0  AUC entro-studio (media sui contrasti malattia-vs-controllo)
  M1  concordanza di segno cross-studio (media sulle 4 malattie)
  eta2 paese / eta2 studio fra i controlli sani
  null di 300 pannelli specie-random appaiati per dimensione e prevalenza
"""
import re, json, sys
import numpy as np, pandas as pd
import tagmos_io

RNG = np.random.default_rng(20260903)
_a = tagmos_io.cli('GMHI (Gupta et al. 2020) benchmarked through the same validation pipeline as the curated composite.',
             species=True)
OUT = _a.out
META_F, SPECIES_F = _a.meta, _a.species
N_NULL = int(__import__('os').environ.get('N_NULL', 300))
DISEASES = ['CRC', 'IBD', 'T2D', 'adenoma']
MIN_N = 8
MIN_HEALTHY_COUNTRY = 30

# ---------------------------------------------------------------- pannelli GMHI
MH_NOMINAL = ['Alistipes senegalensis', 'Bacteroidales bacterium ph8',
              'Bifidobacterium adolescentis', 'Bifidobacterium angulatum',
              'Bifidobacterium catenulatum', 'Lachnospiraceae bacterium 8_1_57FAA',
              'Sutterella wadsworthensis']
MN_NOMINAL = ['Anaerotruncus colihominis','Atopobium parvulum','Bifidobacterium dentium',
 'Blautia producta','candidate division TM7 single cell isolate TM7c',
 'Clostridiales bacterium 1_7_47FAA','Clostridium asparagiforme','Clostridium bolteae',
 'Clostridium citroniae','Clostridium clostridioforme','Clostridium hathewayi',
 'Clostridium nexile','Clostridium ramosum','Clostridium symbiosum','Eggerthella lenta',
 'Erysipelotrichaceae bacterium 2_2_44A','Flavonifractor plautii','Fusobacterium nucleatum',
 'Gemella morbillorum','Gemella sanguinis','Granulicatella adiacens','Holdemania filiformis',
 'Klebsiella pneumoniae','Lachnospiraceae bacterium 1_4_56FAA','Lachnospiraceae bacterium 2_1_58FAA',
 'Lachnospiraceae bacterium 3_1_57FAA_CT1','Lachnospiraceae bacterium 5_1_57FAA',
 'Lachnospiraceae bacterium 9_1_43BFAA','Lactobacillus salivarius','Peptostreptococcus stomatis',
 'Ruminococcaceae bacterium D16','Ruminococcus gnavus','Solobacterium moorei',
 'Streptococcus anginosus','Streptococcus australis','Streptococcus gordonii',
 'Streptococcus infantis','Streptococcus mitis_oralis_pneumoniae','Streptococcus sanguinis',
 'Streptococcus vestibularis','Subdoligranulum sp 4_3_54A2FAA','Subdoligranulum variabile',
 'Veillonella atypica']

# sinonimi di nomenclatura post-MetaPhlAn2 (LPSN / GTDB riclassificazioni)
SYNONYMS = {
 'Atopobium parvulum': ['Lancefieldella parvula'],
 'Clostridium asparagiforme': ['Enterocloster asparagiformis'],
 'Clostridium bolteae': ['Enterocloster bolteae'],
 'Clostridium citroniae': ['Enterocloster citroniae'],
 'Clostridium clostridioforme': ['Enterocloster clostridioformis'],
 'Clostridium hathewayi': ['Hungatella hathewayi'],
 'Clostridium nexile': ['Tyzzerella nexilis'],
 'Clostridium ramosum': ['Erysipelatoclostridium ramosum'],
 'Clostridium symbiosum': ['[Clostridium] symbiosum'],
 'Lactobacillus salivarius': ['Ligilactobacillus salivarius'],
 'Ruminococcus gnavus': ['[Ruminococcus] gnavus'],
 'Streptococcus anginosus': ['Streptococcus milleri'],           # anginosus-group
 'Streptococcus mitis_oralis_pneumoniae': ['Streptococcus mitis', 'Streptococcus oralis',
                                           'Streptococcus pneumoniae'],  # clade MP2 ricomposto
}

def norm(s):
    s = re.sub(r'^[a-z]__', '', str(s).strip())
    s = s.replace('_', ' ')
    return re.sub(r'\s+', ' ', s).strip().lower()

# ---------------------------------------------------------------- metadati
meta_all = pd.read_csv(META_F, sep='\t', low_memory=False,
                   usecols=['study_name','sample_id','subject_id','body_site','study_condition',
                            'disease','age','age_category','country','number_reads'])
# le colonne della matrice specie usano DUE convenzioni: "study|sample_id" e "sample_id" nudo
hdr = open(SPECIES_F).readline().rstrip('\n').split('\t')[1:]
k1 = {s + '|' + str(i): j for j, (s, i) in
      enumerate(zip(meta_all.study_name, meta_all.sample_id))}
uniq = meta_all.drop_duplicates('sample_id', keep=False)
k2 = {str(i): j for i, j in zip(uniq.sample_id, uniq.index)}
row_for_col, cols_ok = [], []
for c in hdr:
    j = k1.get(c, k2.get(c))
    if j is not None:
        row_for_col.append(j); cols_ok.append(c)
print(f'[key] colonne mappate su metadati: {len(cols_ok)}/{len(hdr)}', file=sys.stderr)
meta = meta_all.iloc[row_for_col].copy()
meta['col'] = cols_ok
meta['age_num'] = pd.to_numeric(meta.age, errors='coerce')

keep = (meta.body_site == 'stool')
keep &= ~meta.age_category.isin(['newborn', 'child', 'schoolage'])
keep &= ~(meta.age_num < 18)
meta = meta[keep].copy()
print(f'[meta] campioni stool adulti: {len(meta)}', file=sys.stderr)

# ---------------------------------------------------------------- matrice specie
sp = pd.read_csv(SPECIES_F, sep='\t',
                 index_col=0, dtype={'feature': str})
sp = sp.astype(np.float32)
print(f'[species] {sp.shape[0]} specie x {sp.shape[1]} campioni', file=sys.stderr)

keep_cols = meta.col.tolist()
sp = sp[keep_cols]
meta = meta.reset_index(drop=True)
print(f'[join] campioni analizzabili (stool, adulti): {sp.shape[1]}', file=sys.stderr)

# --- per-sample renormalisation and 1e-5 floor (GMHI_2020 reference implementation)
X = sp.to_numpy(dtype=np.float64)
tot = X.sum(axis=0)
bad = tot <= 0
X = X / np.where(bad, 1.0, tot)
X[X < 1e-5] = 0.0
keep_s = ~bad
X = X[:, keep_s]
meta = meta[keep_s].reset_index(drop=True)
species = list(sp.index)
del sp
print(f'[renorm] campioni con profilo valido: {X.shape[1]}', file=sys.stderr)

richness = (X > 0).sum(axis=0)          # species richness after the floor
prevalence = (X > 0).mean(axis=1)       # prevalenza per specie sul set analizzato

# ---------------------------------------------------------------- mapping pannelli
idx_by_norm = {}
for i, s in enumerate(species):
    idx_by_norm.setdefault(norm(s), []).append(i)

def resolve(name):
    """-> (lista indici, etichetta risoluzione)"""
    if norm(name) in idx_by_norm:
        return idx_by_norm[norm(name)], 'exact'
    for alt in SYNONYMS.get(name, []):
        if norm(alt) in idx_by_norm:
            hits = []
            for a in SYNONYMS[name]:
                hits += idx_by_norm.get(norm(a), [])
            return hits, 'synonym:' + '+'.join(SYNONYMS[name])
    return [], 'absent'

mapping = []
MH_idx, MN_idx = [], []
for lab, panel, store in [('MH', MH_NOMINAL, MH_idx), ('MN', MN_NOMINAL, MN_idx)]:
    for nm in panel:
        hits, how = resolve(nm)
        mapping.append({'panel': lab, 'gmhi_species': nm, 'resolution': how,
                        'cmd_features': ';'.join(species[i] for i in hits),
                        'prevalence': float(prevalence[hits].max()) if hits else np.nan})
        if hits:
            store.append(hits)          # lista di liste: clade ricomposti = 1 pseudo-specie
mapdf = pd.DataFrame(mapping)
mapdf.to_csv(f'{OUT}/gmhi_species_mapping_cMD.tsv', sep='\t', index=False)
print(f'[panel] MH recuperate {len(MH_idx)}/7 · MN recuperate {len(MN_idx)}/43', file=sys.stderr)

def collapse(groups):
    """somma le feature dentro ogni gruppo -> matrice (n_pseudo x n_sample)"""
    if not groups:
        return np.zeros((0, X.shape[1]))
    return np.vstack([X[g, :].sum(axis=0) if len(g) > 1 else X[g[0], :] for g in groups])

# ---------------------------------------------------------------- formula GMHI
def gmhi(mh_mat, mn_mat):
    def psi(M, const):
        pres = M > 0
        R = pres.sum(axis=0)
        with np.errstate(divide='ignore', invalid='ignore'):
            L = np.where(pres, M * np.log(np.where(pres, M, 1.0)), 0.0)
        shannon = -L.sum(axis=0)
        return (R / const) * shannon
    return np.log10((psi(mh_mat, 7.0) + 1e-5) / (psi(mn_mat, 31.0) + 1e-5))

GMHI = gmhi(collapse(MH_idx), collapse(MN_idx))
meta['GMHI'] = GMHI
meta['richness'] = richness

# ---------------------------------------------------------------- metriche
def auc_ctrl_gt_case(ctrl, case):
    """P(score_controllo > score_caso) — >0.5 = direzione attesa (score alto = sano)"""
    n1, n2 = len(ctrl), len(case)
    allv = np.concatenate([ctrl, case])
    r = pd.Series(allv).rank().to_numpy()
    return (r[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n2)

def residualize_within_study(score, study, logrich):
    """residui di score ~ log(ricchezza) stimati entro-studio (§12 del resoconto)"""
    out = np.array(score, dtype=float).copy()
    for s in pd.unique(study):
        m = (study == s).to_numpy()
        if m.sum() < 5:
            out[m] = score[m] - np.mean(score[m]); continue
        x = logrich[m]; y = np.asarray(score)[m]
        if np.std(x) < 1e-12:
            out[m] = y - y.mean(); continue
        b = np.polyfit(x, y, 1)
        out[m] = y - np.polyval(b, x)
    return out

is_ctrl = (meta.study_condition == 'control')
contrasts = []
for d in DISEASES:
    for st, g in meta.groupby('study_name'):
        cs = g[g.study_condition == d]
        ct = g[g.study_condition == 'control']
        if len(cs) >= MIN_N and len(ct) >= MIN_N:
            contrasts.append((d, st, ct.index.to_numpy(), cs.index.to_numpy()))

healthy = meta[is_ctrl].copy()
cc = healthy.country.value_counts()
ok_countries = set(cc[cc >= MIN_HEALTHY_COUNTRY].index)
hidx = healthy.index[healthy.country.isin(ok_countries)].to_numpy()

def eta2(values, groups):
    df = pd.DataFrame({'v': values, 'g': groups}).dropna()
    gm = df.v.mean()
    ssb = df.groupby('g').v.agg(['mean', 'count']).eval('count * (mean - @gm) ** 2').sum()
    sst = ((df.v - gm) ** 2).sum()
    return float(ssb / sst) if sst > 0 else np.nan

def all_metrics(score):
    score = np.asarray(score, dtype=float)
    logr = np.log(np.maximum(richness, 1))
    adj = residualize_within_study(score, meta.study_name, logr)
    res = {}
    for tag, sc in [('raw', score), ('adj', adj)]:
        aucs, per_dis = [], {}
        for d, st, ci, si in contrasts:
            a = auc_ctrl_gt_case(sc[ci], sc[si])
            aucs.append(a); per_dis.setdefault(d, []).append(a)
        res[f'M0_{tag}'] = float(np.mean(aucs))
        conc = []
        for d, v in per_dis.items():
            signs = np.sign(np.array(v) - 0.5)
            signs = signs[signs != 0]
            if len(signs) == 0:
                conc.append(0.5); continue
            conc.append(max((signs > 0).mean(), (signs < 0).mean()))
        res[f'M1_{tag}'] = float(np.mean(conc))
        res[f'M0_per_disease_{tag}'] = {d: float(np.mean(v)) for d, v in per_dis.items()}
    res['eta2_country'] = eta2(score[hidx], meta.country.to_numpy()[hidx])
    res['eta2_study'] = eta2(score[healthy.index.to_numpy()],
                             meta.study_name.to_numpy()[healthy.index.to_numpy()])
    # eta2 paese dopo aggiustamento per studio (§5)
    hs = healthy.index.to_numpy()
    adj_study = residualize_within_study(score[hs], meta.study_name.iloc[hs],
                                         np.zeros(len(hs)))
    sel = np.isin(hs, hidx)
    res['eta2_country_studyadj'] = eta2(adj_study[sel], meta.country.to_numpy()[hs][sel])
    return res

real = all_metrics(GMHI)
print(json.dumps({k: v for k, v in real.items()}, indent=1), file=sys.stderr)

# ---------------------------------------------------------------- null appaiato
panel_flat = set(i for g in MH_idx + MN_idx for i in g)
order = np.argsort(prevalence)
rank_of = np.empty(len(prevalence), dtype=int); rank_of[order] = np.arange(len(prevalence))

def matched_draw(target_idx_groups, taken):
    """pesca una specie con prevalenza simile a quella del gruppo target"""
    prev_t = max(prevalence[i] for i in target_idx_groups)
    j = int(np.searchsorted(prevalence[order], prev_t))
    for w in (25, 50, 100, 200, 400, len(order)):
        lo, hi = max(0, j - w), min(len(order), j + w)
        cand = [c for c in order[lo:hi] if c not in taken and c not in panel_flat]
        if cand:
            return int(RNG.choice(cand))
    raise RuntimeError('no candidate')

def null_panel(n_mh_groups, n_mn_groups):
    taken = set()
    mh, mn = [], []
    for g in n_mh_groups:
        c = matched_draw(g, taken); taken.add(c); mh.append([c])
    for g in n_mn_groups:
        c = matched_draw(g, taken); taken.add(c); mn.append([c])
    return mh, mn

def run_null(mh_targets, mn_targets, n, tag):
    rows = []
    for b in range(n):
        mh, mn = null_panel(mh_targets, mn_targets)
        s = gmhi(collapse(mh), collapse(mn))
        m = all_metrics(s)
        rows.append({'panel': b, 'M0_raw': m['M0_raw'], 'M0_adj': m['M0_adj'],
                     'M1_raw': m['M1_raw'], 'M1_adj': m['M1_adj'],
                     'eta2_country': m['eta2_country'], 'eta2_study': m['eta2_study'],
                     'eta2_country_studyadj': m['eta2_country_studyadj']})
        if (b + 1) % 25 == 0:
            print(f'  [{tag}] {b+1}/{n}', file=sys.stderr)
    return pd.DataFrame(rows)

# null A: appaiato alla dimensione EFFETTIVA del pannello recuperato
nullA = run_null(MH_idx, MN_idx, N_NULL, 'effective')
nullA.to_csv(f'{OUT}/null_effective.tsv', sep='\t', index=False)

# null B: appaiato alla dimensione NOMINALE 7+43 (prevalenze delle specie recuperate riciclate)
mh_t = (MH_idx * 3)[:7]
mn_t = (MN_idx * 2)[:43]
nullB = run_null(mh_t, mn_t, N_NULL, 'nominal')
nullB.to_csv(f'{OUT}/null_nominal_7_43.tsv', sep='\t', index=False)

meta[['study_name','sample_id','country','study_condition','age_category',
      'GMHI','richness']].to_csv(f'{OUT}/gmhi_persample_cMD.tsv', sep='\t', index=False)

summary = {
    'n_samples': int(X.shape[1]),
    'n_studies': int(meta.study_name.nunique()),
    'n_healthy': int(is_ctrl.sum()),
    'n_healthy_countries_ge30': len(ok_countries),
    'n_healthy_in_country_analysis': int(len(hidx)),
    'MH_recovered': len(MH_idx), 'MN_recovered': len(MN_idx),
    'contrasts': [{'disease': d, 'study': s, 'n_ctrl': len(c), 'n_case': len(k)}
                  for d, s, c, k in contrasts],
    'real': real,
}
json.dump(summary, open(f'{OUT}/summary_real.json', 'w'), indent=1, default=float)
print('DONE', file=sys.stderr)
