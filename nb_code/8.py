# Cell 2
import os
import os.path as op
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.stats as stats
from collections import defaultdict
from statsmodels.stats.multitest import multipletests
from matplotlib.patches import Patch
from scipy.stats import gaussian_kde

# ---

# Cell 4
MAGNET_CONDITIONS = ['DCM', 'HCM', 'NF']
GTEX_TISSUES = [
    'Bladder', 'Heart - Left Ventricle', 'Kidney - Cortex',
    'Liver', 'Pancreas', 'Spleen', 'Stomach',
]
ALL_DATASETS = [('MAGNet', c) for c in MAGNET_CONDITIONS] + \
               [('GTEx', t) for t in GTEX_TISSUES]

PART3_RESULTS = {
    'MAGNet': './results_magnet_plau_filter/',
    'GTEx':   './results_gtex_plau_filter/',
}
APPRIS_PATH  = 'appris_data.appris.txt'
DIGGER_PATH  = 'digger_data.csv'
BIOMART_PATH = 'biomart.txt'
TF_LIST_PATH = 'allTFs_hg38.txt'
SF_LIST_PATH = 'comprehensive_sfs.csv'

FIG_DIR  = './sa/figures_annotation/'
TBL_DIR  = './sa/tables_annotation/'
for d in [FIG_DIR, TBL_DIR,
          op.join(FIG_DIR, 'main'), op.join(FIG_DIR, 'supplementary')]:
    os.makedirs(d, exist_ok=True)

N_PERM         = 1000
SEED           = 42
SAVE_FIGURES   = True
SAVE_TABLES    = True
RUN_PERMUTATION = True

# ---

# Cell 6
DS_COLOR = {
    'DCM':                    '#90A4C4',
    'HCM':                    '#ADC9E4',
    'NF':                     '#C0D8E8',
    'Bladder':                '#DBDDEF',
    'Heart - Left Ventricle': '#EBF0FC',
    'Kidney - Cortex':        '#FCEDE9',
    'Liver':                  '#F7DAD4',
    'Pancreas':               '#F5C9B6',
    'Spleen':                 '#F3B8AC',
    'Stomach':                '#F29F95',
}

APPRIS_COLORS_ISO = {'PRINCIPAL': '#44935B', 'ALTERNATIVE': '#7DBA7F', 'MINOR': '#B3D7AE'}
APPRIS_COLORS_BG  = {'PRINCIPAL': '#4884AF', 'ALTERNATIVE': '#78AAC8', 'MINOR': '#B3CEDE'}

TRIFID_BINS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01]
TRIFID_LABELS = ['0.0-0.1','0.1-0.2','0.2-0.3','0.3-0.4','0.4-0.5',
                 '0.5-0.6','0.6-0.7','0.7-0.8','0.8-0.9','0.9-1.0']
TRIFID_COLORS = ['#2171B5','#4292C6','#6BAED6','#9ECAE1','#C6DBEF',
                 '#FEE0D2','#FCBBA1','#FC9272','#FB6A4A','#EF3B2C']

COL_OBS = '#C44E52'
COL_BG  = '#B0B0B0'

TARGET_ORDER = ['target_gene_specific','target_equivalent',
                'target_ambiguous','target_isoform_specific']
SF_ORDER = ['sf_splicing_supported_specific','sf_expression_associated','sf_ambiguous']

def dataset_color(dtype, group):
    return DS_COLOR.get(group, '#888888')

def dataset_label(dtype, group):
    return group

plt.rcParams.update({
    'figure.dpi': 150, 'figure.facecolor': 'white',
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial','Helvetica','DejaVu Sans'],
    'font.size': 13, 'axes.titlesize': 14, 'axes.labelsize': 13,
    'axes.linewidth': 0.6,
    'axes.spines.top': True, 'axes.spines.right': True,
    'axes.spines.bottom': True, 'axes.spines.left': True,
    'axes.grid': False,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.labelsize': 12,  'ytick.labelsize': 12,
    'legend.fontsize': 11,  'legend.frameon': True,
    'savefig.bbox': 'tight', 'savefig.dpi': 300,
})

# ---

# Cell 8
appris = pd.read_csv(APPRIS_PATH, sep='\t')
appris['transcript_id'] = appris['Transcript ID'].str.split('.').str[0]
tx_appris_label = dict(zip(appris['transcript_id'], appris['APPRIS Annotation']))
tx_trifid = dict(zip(appris['transcript_id'], appris['Trifid Score']))

biomart = pd.read_csv(BIOMART_PATH, sep='\t')
gene2txs    = biomart.groupby('Gene stable ID')['Transcript stable ID'].apply(set).to_dict()
tx2gene     = dict(zip(biomart['Transcript stable ID'], biomart['Gene stable ID']))
gene2symbol = dict(zip(biomart['Gene stable ID'], biomart['Gene name']))
symbol2ensg = dict(zip(biomart['Gene name'], biomart['Gene stable ID']))

digger = pd.read_csv(DIGGER_PATH, sep=',')
digger['transcript_id'] = digger['Transcript stable ID'].str.split('.').str[0]
tx_exons = digger.groupby('transcript_id')['Exon stable ID'].apply(
    lambda x: set(x.dropna())).to_dict()
tx_pfam = (digger[digger['Pfam ID'].notna()]
           .groupby('transcript_id')['Pfam ID']
           .apply(lambda x: set(x.dropna())).to_dict())

tf_names = set(pd.read_csv(TF_LIST_PATH, header=None)[0])
sf_df = pd.read_csv(SF_LIST_PATH)
sf_names = set(sf_df.iloc[:, 0])

tf_ensg = {symbol2ensg[s] for s in tf_names if s in symbol2ensg}
sf_ensg = {symbol2ensg[s] for s in sf_names if s in symbol2ensg}

# ---

# Cell 10
def sanitize_prefix(dtype, group):
    return group

def load_dataset_tables(dtype, group):
    rdir   = PART3_RESULTS[dtype]
    prefix = sanitize_prefix(dtype, group)
    d = {}
    for tname, fname in [('t1', 'set_a'), ('t2', 'set_b'), ('t3', 'set_c')]:
        path = op.join(rdir, f'{prefix}_{fname}_plausible.tsv')
        d[tname] = pd.read_csv(path, sep='\t')
    return d

DATA = {}
for dtype, group in ALL_DATASETS:
    DATA[(dtype, group)] = load_dataset_tables(dtype, group)

# ---

# Cell 12
def classify_appris(label):
    if pd.isna(label) or label == '':
        return 'UNKNOWN'
    s = str(label).upper()
    if 'PRINCIPAL' in s:   return 'PRINCIPAL'
    if 'ALTERNATIVE' in s: return 'ALTERNATIVE'
    if 'MINOR' in s:       return 'MINOR'
    return 'UNKNOWN'

def find_principal_txs(gene_id):
    return {tx for tx in gene2txs.get(gene_id, set())
            if 'PRINCIPAL' in str(tx_appris_label.get(tx, '')).upper()}

def compute_unique_missing_domains(tx_id, gene_id):
    label = classify_appris(tx_appris_label.get(tx_id, ''))
    if label != 'PRINCIPAL':
        ref_txs = find_principal_txs(gene_id) or (gene2txs.get(gene_id, set()) - {tx_id})
    else:
        ref_txs = gene2txs.get(gene_id, set()) - {tx_id}
    ref_pfam = set().union(*(tx_pfam.get(t, set()) for t in ref_txs)) if ref_txs else set()
    my_pfam  = tx_pfam.get(tx_id, set())
    return {
        'has_unique_domain':  len(my_pfam - ref_pfam) > 0,
        'has_missing_domain': len(ref_pfam - my_pfam) > 0,
    }

def perm_pvalue(observed, null_dist):
    return (1 + np.sum(null_dist >= observed)) / (1 + len(null_dist))

def plot_perm_hist(ax, null_dist, obs_val, p_val, xlabel, title):
    ax.hist(null_dist, bins=35, color=COL_BG, edgecolor='white',
            linewidth=0.4, alpha=0.75)
    ax.axvline(obs_val, color=COL_OBS, linewidth=2, linestyle='--',
               label=f'Observed = {obs_val:.2f}')
    sig = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'
    ax.text(0.95, 0.92, f'p = {p_val:.4f} ({sig})',
            transform=ax.transAxes, ha='right', va='top', fontsize=8,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='#CCCCCC', alpha=0.9))
    ax.set_xlabel(xlabel); ax.set_ylabel('Frequency')
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=7, loc='upper left')


# ---

# Cell 14
def compute_section_A(key, data):
    t2 = data['t2'].copy()
    st, pl = {}, {}
    st['n_edges'] = len(t2)
    n_iso = (t2['target_category'] == 'target_isoform_specific').sum()
    st['n_isospec'] = int(n_iso)
    st['frac_isospec'] = n_iso / max(len(t2), 1) * 100
    iso = t2[t2['target_category'] == 'target_isoform_specific']

    pl['ntgt_isospec'] = iso.groupby('target_gene')['tgt_n_isoforms'].first().dropna().values
    st['ntgt_median_isospec'] = float(np.median(pl['ntgt_isospec'])) if len(pl['ntgt_isospec']) else np.nan

    pl['dom_isospec'] = iso['tgt_dominance'].dropna().values
    st['dom_median_isospec'] = float(np.median(pl['dom_isospec'])) if len(pl['dom_isospec']) else np.nan

    if 'target_tx_set' in t2.columns:
        iso_with_tx = iso[iso['target_tx_set'].str.len() > 0].copy()
        if len(iso_with_tx) > 0:
            iso_with_tx['target_tx_list'] = iso_with_tx['target_tx_set'].str.split(',')
            exploded = iso_with_tx.explode('target_tx_list')
            exploded = exploded[exploded['target_tx_list'].str.strip().str.len() > 0]
            regs_per_tx = exploded.groupby('target_tx_list')['regulator_gene'].nunique()
            pl['regs_per_tgt_tx'] = regs_per_tx.values
            st['regs_per_tgt_tx_median'] = float(regs_per_tx.median())
    return st, pl

RESULTS_A = {}
for key in DATA:
    st, pl = compute_section_A(key, DATA[key])
    RESULTS_A[key] = {'stats': st, 'plot': pl}

# ---

# Cell 16

all_keys = list(DATA.keys())
all_labels = [dataset_label(*k) for k in all_keys]
all_colors = [DS_COLOR.get(k[1], '#888') for k in all_keys]
n_mag = sum(1 for k in all_keys if k[0] == 'MAGNet')

def add_divider(ax, n_mag, n_total, zero_indexed=True):
    if 0 < n_mag < n_total:
        xpos = (n_mag - 0.5) if zero_indexed else (n_mag + 0.5)
        ax.axvline(x=xpos, color='#333', lw=0.8, ls=':', alpha=0.6)

def raincloud(ax, data_list, labels, colors, ylabel):
    positions = np.arange(1, len(data_list) + 1)
    for i, (vals, col) in enumerate(zip(data_list, colors)):
        pos = positions[i]
        vals = np.asarray(vals, dtype=float)
        vals = vals[np.isfinite(vals)]
        if len(vals) < 3:
            continue
        try:
            kde = gaussian_kde(vals, bw_method=0.3)
            y_range = np.linspace(vals.min(), vals.max(), 200)
            density = kde(y_range)
            density = density / density.max() * 0.35
            ax.fill_betweenx(y_range, pos, pos + density, color=col, alpha=0.5)
        except:
            pass
        jitter = np.random.default_rng(42).uniform(-0.15, 0, size=len(vals))
        ax.scatter(pos + jitter, vals, s=3, alpha=0.2, color=col, edgecolors='none')
        bp = ax.boxplot([vals], positions=[pos], widths=0.12, vert=True,
                        patch_artist=True, showfliers=False,
                        medianprops=dict(color='#333', lw=1),
                        whiskerprops=dict(color='#666', lw=0.5),
                        capprops=dict(color='#666', lw=0.5))
        bp['boxes'][0].set_facecolor(col)
        bp['boxes'][0].set_alpha(0.7)
        bp['boxes'][0].set_edgecolor('#666')
        bp['boxes'][0].set_linewidth(0.5)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=12, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=13, fontweight='bold')
    add_divider(ax, n_mag, len(data_list), zero_indexed=False)

fig, ax = plt.subplots(figsize=(6, 4.5))
x = np.arange(len(all_keys))
frac_vals = [RESULTS_A[k]['stats'].get('frac_isospec', 0) for k in all_keys]
ax.bar(x, frac_vals, color=all_colors, edgecolor='white', linewidth=0.5, width=0.65)
ax.set_xticks(x)
ax.set_xticklabels(all_labels, rotation=45, ha='right', fontsize=12, fontweight='bold')
ax.set_ylabel('Fraction', fontsize=13, fontweight='bold')
add_divider(ax, n_mag, len(all_keys), zero_indexed=True)
for label in ax.get_yticklabels(): label.set_fontweight('bold')
ax.yaxis.get_offset_text().set_fontweight('bold')
ax.yaxis.get_offset_text().set_fontsize(12)
plt.tight_layout()
if SAVE_FIGURES: plt.savefig(op.join(FIG_DIR, 'main', 'A1_isospec_fraction.pdf'))
plt.show()

a2_data = [RESULTS_A[k]['plot'].get('ntgt_isospec', np.array([])) for k in all_keys]
a2_valid = [(d, l, c) for d, l, c in zip(a2_data, all_labels, all_colors) if len(d) > 0]
if a2_valid:
    fig, ax = plt.subplots(figsize=(6, 4.5))
    raincloud(ax, [v[0] for v in a2_valid], [v[1] for v in a2_valid],
              [v[2] for v in a2_valid], 'Expressed isoforms')
    for label in ax.get_yticklabels(): label.set_fontweight('bold')
    ax.yaxis.get_offset_text().set_fontweight('bold')
    plt.tight_layout()
    if SAVE_FIGURES: plt.savefig(op.join(FIG_DIR, 'main', 'A2_ntgt_raincloud.pdf'))
    plt.show()

a3_data = [RESULTS_A[k]['plot'].get('dom_isospec', np.array([])) for k in all_keys]
a3_valid = [(d, l, c) for d, l, c in zip(a3_data, all_labels, all_colors) if len(d) > 0]
if a3_valid:
    fig, ax = plt.subplots(figsize=(6, 4.5))
    raincloud(ax, [v[0] for v in a3_valid], [v[1] for v in a3_valid],
              [v[2] for v in a3_valid], 'Dominance')
    for label in ax.get_yticklabels(): label.set_fontweight('bold')
    ax.yaxis.get_offset_text().set_fontweight('bold')
    plt.tight_layout()
    if SAVE_FIGURES: plt.savefig(op.join(FIG_DIR, 'main', 'A3_dominance_raincloud.pdf'))
    plt.show()

a4_data = [RESULTS_A[k]['plot'].get('regs_per_tgt_tx', np.array([])) for k in all_keys]
a4_valid = [(d, l, c) for d, l, c in zip(a4_data, all_labels, all_colors) if len(d) > 0]
if a4_valid:
    fig, ax = plt.subplots(figsize=(6, 4.5))
    raincloud(ax, [v[0] for v in a4_valid], [v[1] for v in a4_valid],
              [v[2] for v in a4_valid], 'Unique regulators')
    for label in ax.get_yticklabels(): label.set_fontweight('bold')
    ax.yaxis.get_offset_text().set_fontweight('bold')
    plt.tight_layout()
    if SAVE_FIGURES: plt.savefig(op.join(FIG_DIR, 'main', 'A4_regs_raincloud.pdf'))
    plt.show()

# ---

# Cell 18
def compute_section_B(key, data):
    t3 = data['t3'].copy()
    st, pl = {}, {}
    st['n_edges'] = len(t3)
    n_sup = (t3['sf_category'] == 'sf_splicing_supported_specific').sum()
    st['n_supported'] = int(n_sup)
    st['frac_supported'] = n_sup / max(len(t3), 1) * 100
    sup = t3[t3['sf_category'] == 'sf_splicing_supported_specific']

    pl['ntgt_sup'] = sup.groupby('target_gene')['tgt_n_isoforms'].first().dropna().values
    st['ntgt_median_sup'] = float(np.median(pl['ntgt_sup'])) if len(pl['ntgt_sup']) else np.nan

    pl['dom_sup'] = sup['dominance'].dropna().values
    st['dom_median_sup'] = float(np.median(pl['dom_sup'])) if len(pl['dom_sup']) else np.nan

    if 'target_tx_set' in t3.columns:
        sup_with_tx = sup[sup['target_tx_set'].str.len() > 0].copy()
        if len(sup_with_tx) > 0:
            sup_with_tx['target_tx_list'] = sup_with_tx['target_tx_set'].str.split(',')
            exploded = sup_with_tx.explode('target_tx_list')
            exploded = exploded[exploded['target_tx_list'].str.strip().str.len() > 0]
            sfs_per_tx = exploded.groupby('target_tx_list')['sf_gene'].nunique()
            pl['sfs_per_tgt_tx'] = sfs_per_tx.values
            st['sfs_per_tgt_tx_median'] = float(sfs_per_tx.median())
    return st, pl

RESULTS_B = {}
for key in DATA:
    st, pl = compute_section_B(key, DATA[key])
    RESULTS_B[key] = {'stats': st, 'plot': pl}

# ---

# Cell 20
fig, ax = plt.subplots(figsize=(6, 4.5))
x = np.arange(len(all_keys))
sup_fracs = [RESULTS_B[k]['stats'].get('frac_supported', 0) for k in all_keys]
ax.bar(x, sup_fracs, color=all_colors, edgecolor='white', linewidth=0.5, width=0.65)
ax.set_xticks(x)
ax.set_xticklabels(all_labels, rotation=45, ha='right', fontsize=12, fontweight='bold')
ax.set_ylabel('Fraction', fontsize=13, fontweight='bold')
add_divider(ax, n_mag, len(all_keys), zero_indexed=True)
for label in ax.get_yticklabels(): label.set_fontweight('bold')
ax.yaxis.get_offset_text().set_fontweight('bold')
ax.yaxis.get_offset_text().set_fontsize(12)
plt.tight_layout()
if SAVE_FIGURES: plt.savefig(op.join(FIG_DIR, 'main', 'B1_supported_fraction.pdf'))
plt.show()

b2_data = [RESULTS_B[k]['plot'].get('ntgt_sup', np.array([])) for k in all_keys]
b2_valid = [(d, l, c) for d, l, c in zip(b2_data, all_labels, all_colors) if len(d) > 0]
if b2_valid:
    fig, ax = plt.subplots(figsize=(6, 4.5))
    raincloud(ax, [v[0] for v in b2_valid], [v[1] for v in b2_valid],
              [v[2] for v in b2_valid], 'Expressed isoforms')
    for label in ax.get_yticklabels(): label.set_fontweight('bold')
    ax.yaxis.get_offset_text().set_fontweight('bold')
    plt.tight_layout()
    if SAVE_FIGURES: plt.savefig(op.join(FIG_DIR, 'main', 'B2_ntgt_raincloud.pdf'))
    plt.show()

b3_data = [RESULTS_B[k]['plot'].get('dom_sup', np.array([])) for k in all_keys]
b3_valid = [(d, l, c) for d, l, c in zip(b3_data, all_labels, all_colors) if len(d) > 0]
if b3_valid:
    fig, ax = plt.subplots(figsize=(6, 4.5))
    raincloud(ax, [v[0] for v in b3_valid], [v[1] for v in b3_valid],
              [v[2] for v in b3_valid], 'Dominance')
    for label in ax.get_yticklabels(): label.set_fontweight('bold')
    ax.yaxis.get_offset_text().set_fontweight('bold')
    plt.tight_layout()
    if SAVE_FIGURES: plt.savefig(op.join(FIG_DIR, 'main', 'B3_dominance_raincloud.pdf'))
    plt.show()

b4_data = [RESULTS_B[k]['plot'].get('sfs_per_tgt_tx', np.array([])) for k in all_keys]
b4_valid = [(d, l, c) for d, l, c in zip(b4_data, all_labels, all_colors) if len(d) > 0]
if b4_valid:
    fig, ax = plt.subplots(figsize=(6, 4.5))
    raincloud(ax, [v[0] for v in b4_valid], [v[1] for v in b4_valid],
              [v[2] for v in b4_valid], 'Unique SF regulators')
    for label in ax.get_yticklabels(): label.set_fontweight('bold')
    ax.yaxis.get_offset_text().set_fontweight('bold')
    plt.tight_layout()
    if SAVE_FIGURES: plt.savefig(op.join(FIG_DIR, 'main', 'B4_sfs_raincloud.pdf'))
    plt.show()

# ---

# Cell 22
def compute_section_C(key, data):
    rng = np.random.default_rng(SEED)
    t1 = data['t1'].copy()
    st, pl = {}, {}
    t1['appris_coarse'] = t1['reg_appris'].apply(classify_appris)
    t1['trifid_score']  = t1['best_tx'].map(tx_trifid)

    tx_pairs = t1[['best_tx','regulator_gene']].drop_duplicates(subset=['best_tx'])
    recs = []
    for _, row in tx_pairs.iterrows():
        if pd.isna(row['best_tx']): continue
        r = compute_unique_missing_domains(row['best_tx'], row['regulator_gene'])
        r['best_tx'] = row['best_tx']
        recs.append(r)
    if recs:
        adf = pd.DataFrame(recs)
        t1 = t1.merge(adf[['best_tx','has_unique_domain','has_missing_domain']],
                       on='best_tx', how='left')

    iso_t1 = t1[t1['source_category']=='source_isoform_specific']
    iso_tx = iso_t1[['best_tx','appris_coarse','trifid_score']].drop_duplicates(subset=['best_tx'])
    bg_tx  = t1[['best_tx','appris_coarse','trifid_score']].drop_duplicates(subset=['best_tx'])

    st['n_iso_tx'] = len(iso_tx)
    st['n_bg_tx']  = len(bg_tx)
    st['n_edges']  = len(t1)

    # APPRIS: only PRINCIPAL/ALTERNATIVE/MINOR (exclude UNKNOWN)
    for prefix, tx_df in [('iso', iso_tx), ('bg', bg_tx)]:
        known = tx_df[tx_df['appris_coarse'].isin(['PRINCIPAL','ALTERNATIVE','MINOR'])]
        total_known = max(len(known), 1)
        for cat in ['PRINCIPAL','ALTERNATIVE','MINOR']:
            st[f'{prefix}_appris_{cat}_pct'] = (known['appris_coarse']==cat).sum() / total_known * 100

    # TRIFID bins
    for prefix, tx_df in [('iso', iso_tx), ('bg', bg_tx)]:
        scores = tx_df['trifid_score'].dropna().values
        counts, _ = np.histogram(scores, bins=TRIFID_BINS)
        total = max(counts.sum(), 1)
        for j, lbl in enumerate(TRIFID_LABELS):
            st[f'{prefix}_trifid_bin_{lbl}'] = counts[j] / total * 100

    # Permutation
    if RUN_PERMUTATION and 'has_missing_domain' in t1.columns:
        def _count_np(edf, feat):
            txdf = edf[['best_tx','appris_coarse',feat]].drop_duplicates(subset=['best_tx'])
            return txdf[txdf['appris_coarse'] != 'PRINCIPAL'][feat].sum()
        obs_miss = _count_np(iso_t1, 'has_missing_domain')
        obs_uniq = _count_np(iso_t1, 'has_unique_domain')
        n_iso_e = len(iso_t1)
        bg_miss, bg_uniq = np.zeros(N_PERM), np.zeros(N_PERM)
        for i in range(N_PERM):
            idx = rng.choice(len(t1), size=min(n_iso_e, len(t1)), replace=False)
            s = t1.iloc[idx]
            bg_miss[i] = _count_np(s, 'has_missing_domain')
            bg_uniq[i] = _count_np(s, 'has_unique_domain')
        st['perm_missing_p'] = perm_pvalue(obs_miss, bg_miss)
        st['perm_unique_p']  = perm_pvalue(obs_uniq, bg_uniq)
        pl['perm_miss'] = {'obs': float(obs_miss), 'null': bg_miss}
        pl['perm_uniq'] = {'obs': float(obs_uniq), 'null': bg_uniq}

    return st, pl

RESULTS_C = {}
for key in DATA:
    st, pl = compute_section_C(key, DATA[key])
    RESULTS_C[key] = {'stats': st, 'plot': pl}

# ---

# Cell 24
appris_cats = ['PRINCIPAL', 'ALTERNATIVE', 'MINOR']
mag_keys = [k for k in DATA if k[0]=='MAGNet']
gtex_keys = [k for k in DATA if k[0]=='GTEx']

def _appris_sidebyside(ax, keys, results):
    x = np.arange(len(keys))
    w = 0.35
    bot_iso = np.zeros(len(keys)); bot_bg = np.zeros(len(keys))
    handles = []
    for cat in appris_cats:
        iso_vals = np.array([results[k]['stats'].get(f'iso_appris_{cat}_pct', 0) for k in keys])
        bg_vals  = np.array([results[k]['stats'].get(f'bg_appris_{cat}_pct', 0) for k in keys])
        b1 = ax.bar(x - w/2, iso_vals, w, bottom=bot_iso,
                     color=APPRIS_COLORS_ISO[cat], edgecolor='white', linewidth=0.3)
        b2 = ax.bar(x + w/2, bg_vals, w, bottom=bot_bg,
                     color=APPRIS_COLORS_BG[cat], edgecolor='white', linewidth=0.3)
        handles.append((b1, f'{cat} (s)'))
        handles.append((b2, f'{cat} (f)'))
        bot_iso += iso_vals
        bot_bg  += bg_vals
    ax.set_xticks(x)
    ax.set_xticklabels([dataset_label(*k) for k in keys], rotation=45, ha='right', fontsize=12, fontweight='bold')
    ax.set_ylabel('Fraction', fontsize=13, fontweight='bold')
    ax.set_ylim(0, 105)
    for label in ax.get_yticklabels(): label.set_fontweight('bold')
    # Legend: s_i_s group first, then f_s_a group
    iso_handles = [h for h, l in handles if '(s)' in l]
    iso_labels  = [l for h, l in handles if '(s)' in l]
    bg_handles  = [h for h, l in handles if '(f)' in l]
    bg_labels   = [l for h, l in handles if '(f)' in l]
    ax.legend(iso_handles + bg_handles, iso_labels + bg_labels,
              loc='upper left', bbox_to_anchor=(1.01, 1), ncol=1,
              prop={'weight': 'bold', 'size': 9})

# C1: APPRIS - MAGNet
fig, ax = plt.subplots(figsize=(6, 4))
_appris_sidebyside(ax, mag_keys, RESULTS_C)
plt.tight_layout()
if SAVE_FIGURES: plt.savefig(op.join(FIG_DIR, 'main', 'C1_appris_magnet.pdf'))
plt.show()

# C2: APPRIS - GTEx
fig, ax = plt.subplots(figsize=(10, 4))
_appris_sidebyside(ax, gtex_keys, RESULTS_C)
plt.tight_layout()
if SAVE_FIGURES: plt.savefig(op.join(FIG_DIR, 'main', 'C2_appris_gtex.pdf'))
plt.show()

def _trifid_sidebyside(ax, keys, results):
    x = np.arange(len(keys))
    w = 0.35
    bot_iso = np.zeros(len(keys)); bot_bg = np.zeros(len(keys))
    for j, lbl in enumerate(TRIFID_LABELS):
        iso_vals = np.array([results[k]['stats'].get(f'iso_trifid_bin_{lbl}', 0) for k in keys])
        bg_vals  = np.array([results[k]['stats'].get(f'bg_trifid_bin_{lbl}', 0) for k in keys])
        ax.bar(x - w/2, iso_vals, w, bottom=bot_iso,
               color=TRIFID_COLORS[j], edgecolor='white', linewidth=0.2,
               label=lbl if True else '')
        ax.bar(x + w/2, bg_vals, w, bottom=bot_bg,
               color=TRIFID_COLORS[j], edgecolor='white', linewidth=0.2)
        bot_iso += iso_vals
        bot_bg  += bg_vals
    ax.set_xticks(x)
    ax.set_xticklabels([dataset_label(*k) for k in keys], rotation=45, ha='right', fontsize=12, fontweight='bold')
    ax.set_ylabel('Fraction', fontsize=13, fontweight='bold')
    ax.set_ylim(0, 105)
    for label in ax.get_yticklabels(): label.set_fontweight('bold')
    # Single legend (10 bins, one color each)
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), ncol=1,
              prop={'weight': 'bold', 'size': 7})

# C3: TRIFID - MAGNet
fig, ax = plt.subplots(figsize=(6, 4))
_trifid_sidebyside(ax, mag_keys, RESULTS_C)
plt.tight_layout()
if SAVE_FIGURES: plt.savefig(op.join(FIG_DIR, 'main', 'C3_trifid_magnet.pdf'))
plt.show()

# C4: TRIFID - GTEx
fig, ax = plt.subplots(figsize=(10, 4))
_trifid_sidebyside(ax, gtex_keys, RESULTS_C)
plt.tight_layout()
if SAVE_FIGURES: plt.savefig(op.join(FIG_DIR, 'main', 'C4_trifid_gtex.pdf'))
plt.show()

for key in DATA:
    pm = RESULTS_C[key]['plot'].get('perm_miss')
    pu = RESULTS_C[key]['plot'].get('perm_uniq')
    if pm is None or pu is None:
        continue
    label = dataset_label(*key)

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10, 3.5))

    for ax, pdata, pval_key, xlabel_short in [
        (ax0, pu, 'perm_unique_p',  'Non-principal with unique domain'),
        (ax1, pm, 'perm_missing_p', 'Non-principal with missing domain'),
    ]:
        null = pdata['null']
        obs = pdata['obs']
        p_val = RESULTS_C[key]['stats'][pval_key]

        all_vals = np.concatenate([null, [obs]])
        xmin = max(0, all_vals.min() - 3)
        xmax = all_vals.max() + 3

        ax.hist(null, bins=np.arange(xmin, xmax + 1, 1), color='#AEC7E8', alpha=0.6,
                edgecolor='white', linewidth=0.3, label='Null distribution')
        ax.axvline(obs, color='#C44E52', lw=2, ls='--',
                   label=f'Observed = {obs:.0f}')

        sig = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'
        ax.text(0.95, 0.92, f'p = {p_val:.4f} ({sig})',
                transform=ax.transAxes, ha='right', va='top', fontsize=11, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          edgecolor='#CCC', alpha=0.9))
        ax.set_xlabel(xlabel_short, fontsize=13, fontweight='bold')
        ax.set_ylabel('Count', fontsize=13, fontweight='bold')
        ax.set_xlim(xmin, xmax)
        for label in ax.get_yticklabels(): label.set_fontweight('bold')
        for label in ax.get_xticklabels(): label.set_fontweight('bold')
        ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1),
                  prop={'weight': 'bold', 'size': 9})

    plt.tight_layout()
    if SAVE_FIGURES:
        safe = key[1].replace(' ', '_').replace('-', '')
        plt.savefig(op.join(FIG_DIR, 'main', f'C5C6_perm_{safe}.pdf'))
    plt.show()

# ---

# Cell 26
summary_rows = []
for key in DATA:
    dtype, group = key
    A = RESULTS_A[key]['stats']
    B = RESULTS_B[key]['stats']
    C = RESULTS_C[key]['stats']
    summary_rows.append({
        'dataset_type': dtype, 'group': group,
        'SetB_edges': A.get('n_edges'),
        'SetB_frac_isospec': A.get('frac_isospec'),
        'SetB_ntgt_median': A.get('ntgt_median_isospec'),
        'SetB_dom_median': A.get('dom_median_isospec'),
        'SetC_edges': B.get('n_edges'),
        'SetC_frac_supported': B.get('frac_supported'),
        'SetC_ntgt_median': B.get('ntgt_median_sup'),
        'SetC_dom_median': B.get('dom_median_sup'),
        'SetA_edges': C.get('n_edges'),
        'SetA_iso_MINOR_pct': C.get('iso_appris_MINOR_pct'),
        'SetA_perm_missing_p': C.get('perm_missing_p'),
        'SetA_perm_unique_p': C.get('perm_unique_p'),
    })
summary_df = pd.DataFrame(summary_rows)
display(summary_df)
if SAVE_TABLES:
    summary_df.to_csv(op.join(TBL_DIR, 'summary_all_datasets.csv'), index=False)

# ---