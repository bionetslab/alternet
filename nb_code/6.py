# Cell 2
import os, os.path as op, re, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

plt.rcParams.update({
    'figure.dpi': 150, 'figure.facecolor': 'white',
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 13, 'axes.titlesize': 14, 'axes.labelsize': 13,
    'axes.linewidth': 0.6,
    'axes.spines.top': True, 'axes.spines.right': True,
    'axes.spines.bottom': True, 'axes.spines.left': True,
    'axes.grid': False,
    'xtick.labelsize': 12, 'ytick.labelsize': 12,
    'legend.fontsize': 11, 'legend.frameon': True,
    'savefig.bbox': 'tight', 'savefig.dpi': 300,
})

ORA_DIR = './ora/'
os.makedirs(ORA_DIR, exist_ok=True)
P_THRESHOLD = 0.05

# Colors for L1/L2/L3
L_COLORS_C = {'L1': '#B3D7AE', 'L2': '#7DBA7F', 'L3': '#44935B'}  # C_hits (green)
L_COLORS_B = {'L1': '#B3CEDE', 'L2': '#78AAC8', 'L3': '#4884AF'}  # B_hits (blue)

DATASETS = pd.DataFrame([
    {'dataset_id': 'DCM',                    'cohort': 'MAGNet', 'label': 'DCM',                    'n': 166, 'csv': 'DCM_gprofiler.csv'},
    {'dataset_id': 'HCM',                    'cohort': 'MAGNet', 'label': 'HCM',                    'n': 28,  'csv': 'HCM_gprofiler.csv'},
    {'dataset_id': 'NF',                     'cohort': 'MAGNet', 'label': 'NF',                     'n': 166, 'csv': 'NF_gprofiler.csv'},
    {'dataset_id': 'Bladder',                'cohort': 'GTEx',   'label': 'Bladder',                'n': 77,  'csv': 'Bladder.csv'},
    {'dataset_id': 'Heart - Left Ventricle', 'cohort': 'GTEx',   'label': 'Heart - Left Ventricle', 'n': 452, 'csv': 'Heart_LV.csv'},
    {'dataset_id': 'Kidney - Cortex',        'cohort': 'GTEx',   'label': 'Kidney - Cortex',        'n': 104, 'csv': 'Kidney.csv'},
    {'dataset_id': 'Liver',                  'cohort': 'GTEx',   'label': 'Liver',                  'n': 262, 'csv': 'Liver.csv'},
    {'dataset_id': 'Pancreas',               'cohort': 'GTEx',   'label': 'Pancreas',               'n': 362, 'csv': 'Pancreas.csv'},
    {'dataset_id': 'Spleen',                 'cohort': 'GTEx',   'label': 'Spleen',                 'n': 277, 'csv': 'Spleen.csv'},
    {'dataset_id': 'Stomach',                'cohort': 'GTEx',   'label': 'Stomach',                'n': 407, 'csv': 'Stomach.csv'},
])
DATASETS = DATASETS.set_index('dataset_id')
DS_ORDER = list(DATASETS.index)
DS_LABELS = DATASETS['label'].to_dict()


# ---

# Cell 4
SET_A = [
    "mitochond", "mitochondrial",
    "oxidative phosphorylation", "respiration", "electron transport",
    "ribosome", "translation",
    "proteasome", "ubiquitin", "protein catabolic"
]

SET_B = [
    "transcription", "rna processing", "mrna", "splicing",
    "chromatin", "nuclear", "polymerase",
    "gene expression",
]

SET_C = {
    "Liver": [
        "xenobiotic", "bile acid", "bile salt", "bile",
        "lipoprotein", "cholesterol", "gluconeogenesis",
        "peroxisom", "coagulation",
    ],
    "Pancreas": [
        "secretory granule", "regulated exocytosis",
        "endoplasmic reticulum", "unfolded protein",
        "insulin", "glucose homeostasis",
        "digestive", "zymogen", "pancreatic",
    ],
    "Stomach": [
        "gastric acid", "acid secretion", "proton transport",
        "digestion", "mucin", "mucus",
        "smooth muscle", "muscle contraction",
    ],
    "Heart - Left Ventricle": [
        "cardiac muscle contraction", "cardiac conduction",
        "action potential", "ion channel",
        "calcium", "sarcoplasmic reticulum",
        "sarcomere", "myofibril",
        "heart development", "heart rate",
    ],
    "Bladder": [
        "urothel", "urothelial", "bladder",
        "tight junction", "keratinization",
        "detrusor", "micturition", "urinary",
    ],
    "Spleen": [
        "antigen presentation", "mhc",
        "b cell", "t cell", "lymphocyte",
        "interferon", "complement",
        "phagocytosis", "innate immune",
    ],
    "Kidney - Cortex": [
        "glomerul", "filtration", "nephron", 
        "renal", "tubule", "podocyte",
        "solute", "sodium ion transport",
    ],
    "DCM": [
        "extracellular matrix", "collagen", "fibrosis",
        "cardiac muscle contraction", "sarcomere", "myofibril",
        "calcium", "sarcoplasmic reticulum",
        "heart failure", "dilation", "ventricular",
    ],
    "HCM": [
        "sarcomere", "myofibril", "myosin", "troponin",
        "hypertrophy", "hypertrophic",
        "calcium", "sarcoplasmic reticulum",
        "mechanical stimulus", "z disc",
    ],
    "NF": [
        "cardiac muscle contraction", "cardiac conduction",
        "action potential", "ion channel",
        "calcium", "sarcoplasmic reticulum",
        "sarcomere", "myofibril",
        "heart development", "heart rate",
    ],
}


# ---

# Cell 6
def load_gprofiler(filepath, dataset_id):
    """Load one g:Profiler CSV, return standardized long-form DataFrame."""
    raw = pd.read_csv(filepath)

    # Parse highlighted
    raw['hl_L1'] = raw['highlighted'].str.split(',').str[0].str.strip() == 'true'
    raw['hl_L2'] = raw['highlighted'].str.split(',').str[1].str.strip() == 'true'
    raw['hl_L3'] = raw['highlighted'].str.split(',').str[2].str.strip() == 'true'

    col_map = {
        'adjusted_p_value__L1_Canonical': 'pval_L1',
        'adjusted_p_value__L2_Source_Isoform': 'pval_L2',
        'adjusted_p_value__L3_Target_Isoform': 'pval_L3',
    }
    raw = raw.rename(columns=col_map)

    # Melt to long form: one row per (term, query)
    rows = []
    for _, r in raw.iterrows():
        for L in ['L1', 'L2', 'L3']:
            rows.append({
                'dataset_id': dataset_id,
                'query': L,
                'term_id': r.get('term_id', r.get('native', '')),
                'term_name': r['term_name'],
                'source': r.get('source', ''),
                'p_adj': r[f'pval_{L}'],
                'is_highlighted': r[f'hl_{L}'],
            })
    return pd.DataFrame(rows)

frames = []
for ds_id in DS_ORDER:
    csv_path = DATASETS.loc[ds_id, 'csv']
    try:
        df = load_gprofiler(csv_path, ds_id)
        frames.append(df)
        n_terms = df['term_name'].nunique()
    except FileNotFoundError:

ALL_TERMS = pd.concat(frames, ignore_index=True)
loaded_ids = ALL_TERMS['dataset_id'].unique().tolist()
loaded_order = [ds for ds in DS_ORDER if ds in loaded_ids]

# ---

# Cell 8
def normalize_term(s):
    if pd.isna(s):
        return ''
    s = str(s).lower()
    s = re.sub(r'[-_/]', ' ', s)
    s = re.sub(r'\s+', ' ', s).strip()
    return s

ALL_TERMS['term_norm'] = ALL_TERMS['term_name'].apply(normalize_term)

# ---

# Cell 10
def compile_pat(keyword_list):
    parts = []
    for kw in keyword_list:
        escaped = re.escape(kw.lower())
        if len(kw) <= 3:
            parts.append(r'\b' + escaped + r'\b')
        else:
            parts.append(escaped)
    return re.compile('|'.join(parts), re.IGNORECASE)

pat_A = compile_pat(SET_A)
pat_B = compile_pat(SET_B)
pat_C = {ds: compile_pat(kws) for ds, kws in SET_C.items()}

# Filter to significant terms, deduplicate per (dataset, query, term_id)
SIG = ALL_TERMS[ALL_TERMS['p_adj'] < P_THRESHOLD].copy()
SIG = SIG.drop_duplicates(subset=['dataset_id', 'query', 'term_id'])

# Step 1: Independent boolean flags (non-disjoint)
SIG['is_A'] = SIG['term_norm'].str.contains(pat_A, regex=True)
SIG['is_B'] = SIG['term_norm'].str.contains(pat_B, regex=True)
SIG['is_C'] = SIG.apply(
    lambda r: bool(pat_C.get(r['dataset_id']) and pat_C[r['dataset_id']].search(r['term_norm'])),
    axis=1)

# Step 2: Disjoint main_set for B vs C only (priority C > B)
def assign_main(term_norm, dataset_id):
    pc = pat_C.get(dataset_id)
    if pc and pc.search(term_norm):
        return 'C'
    if pat_B.search(term_norm):
        return 'B'
    return None

SIG['main_set'] = SIG.apply(lambda r: assign_main(r['term_norm'], r['dataset_id']), axis=1)

# Backward-compatible kw_set: uses main_set for B/C, falls back to is_A for A
SIG['kw_set'] = SIG['main_set']
SIG.loc[(SIG['kw_set'].isna()) & (SIG['is_A']), 'kw_set'] = 'A'


# ---

# Cell 12
summary_rows = []
summary_hl_rows = []

for ds_id in loaded_order:
    for L in ['L1', 'L2', 'L3']:
        # --- Significant terms ---
        sub = SIG[(SIG['dataset_id'] == ds_id) & (SIG['query'] == L)]
        total = len(sub)
        a = sub['is_A'].sum()              # non-disjoint
        b = (sub['main_set'] == 'B').sum()  # disjoint B vs C
        c = (sub['main_set'] == 'C').sum()  # disjoint B vs C
        summary_rows.append({
            'dataset_id': ds_id,
            'label': DS_LABELS[ds_id],
            'cohort': DATASETS.loc[ds_id, 'cohort'],
            'query': L,
            'total_sig': total,
            'A_hits': int(a), 'B_hits': int(b), 'C_hits': int(c),
            'A_frac': a / max(total, 1),
            'B_frac': b / max(total, 1),
            'C_frac': c / max(total, 1),
        })

        # --- Highlighted terms ---
        sub_hl = sub[sub['is_highlighted'] == True]
        total_hl = len(sub_hl)
        a_hl = (sub_hl['kw_set'] == 'A').sum()
        b_hl = (sub_hl['kw_set'] == 'B').sum()
        c_hl = (sub_hl['kw_set'] == 'C').sum()
        summary_hl_rows.append({
            'dataset_id': ds_id,
            'label': DS_LABELS[ds_id],
            'cohort': DATASETS.loc[ds_id, 'cohort'],
            'query': L,
            'total_sig': total_hl,
            'A_hits': int(a_hl), 'B_hits': int(b_hl), 'C_hits': int(c_hl),
            'A_frac': a_hl / max(total_hl, 1),
            'B_frac': b_hl / max(total_hl, 1),
            'C_frac': c_hl / max(total_hl, 1),
        })

summary_df = pd.DataFrame(summary_rows)
summary_hl_df = pd.DataFrame(summary_hl_rows)

for metric in ['C_hits', 'B_hits', 'total_sig']:
    piv = summary_df.pivot(index='label', columns='query', values=metric)
    piv = piv.reindex([DS_LABELS[d] for d in loaded_order])[['L1', 'L2', 'L3']]
    display(piv)

# ---

# Cell 14
C_L3_ONLY = {}
for ds_id in loaded_order:
    c_terms = {}
    for L in ['L1', 'L2', 'L3']:
        sub = SIG[(SIG['dataset_id']==ds_id) & (SIG['query']==L) & (SIG['kw_set']=='C')]
        c_terms[L] = set(sub['term_norm'])
    C_L3_ONLY[ds_id] = c_terms['L3'] - c_terms['L1'] - c_terms['L2']

idx_rows = []
for ds_id in loaded_order:
    s = summary_df[summary_df['dataset_id'] == ds_id].set_index('query')
    r = {
        'dataset_id': ds_id,
        'label': DS_LABELS[ds_id],
        'cohort': DATASETS.loc[ds_id, 'cohort'],
        'RI':     round(s.loc['L2', 'B_frac'] - s.loc['L1', 'B_frac'], 4),
        'TG':     round(s.loc['L3', 'C_frac'] - s.loc['L1', 'C_frac'], 4),
        'TG32':   round(s.loc['L3', 'C_frac'] - s.loc['L2', 'C_frac'], 4),
        'DC_hits':    int(s.loc['L3', 'C_hits'] - s.loc['L1', 'C_hits']),
        'DC_pp':  round((s.loc['L3', 'C_frac'] - s.loc['L1', 'C_frac']) * 100, 2),
        'C_L3_only':  len(C_L3_ONLY[ds_id]),
        'L2_drift': (s.loc['L2', 'B_frac'] > s.loc['L1', 'B_frac'] and
                     s.loc['L2', 'C_frac'] < s.loc['L1', 'C_frac']),
    }
    idx_rows.append(r)

indices_df = pd.DataFrame(idx_rows)
display(indices_df)

indices_df.to_csv(op.join(ORA_DIR, 'indices.csv'), index=False)
summary_df.to_csv(op.join(ORA_DIR, 'summary_per_level.csv'), index=False)


# ---

# Cell 16
pw_rows = []
for ds_id in loaded_order:
    s = summary_df[summary_df['dataset_id'] == ds_id].set_index('query')
    label = DS_LABELS[ds_id]

    # L2 vs L1: L2 captures more regulation (Set B)
    b_l1, b_l2 = s.loc['L1','B_hits'], s.loc['L2','B_hits']
    pw_rows.append({
        'Dataset': label, 'Comparison': 'L2 vs L1',
        'Hypothesis': 'L2 captures more regulation',
        'Metric': 'B_hits',
        'L_low': int(b_l1), 'L_high': int(b_l2),
        'Delta': int(b_l2 - b_l1),
        'Confirmed': b_l2 > b_l1,
    })

    # L3 vs L2: L3 captures more tissue (Set C)
    c_l2, c_l3 = s.loc['L2','C_hits'], s.loc['L3','C_hits']
    pw_rows.append({
        'Dataset': label, 'Comparison': 'L3 vs L2',
        'Hypothesis': 'L3 captures more tissue biology',
        'Metric': 'C_hits',
        'L_low': int(c_l2), 'L_high': int(c_l3),
        'Delta': int(c_l3 - c_l2),
        'Confirmed': c_l3 > c_l2,
    })

    # L3 vs L1: primary = tissue gain, secondary = regulation gain
    b_l1, b_l3 = s.loc['L1','B_hits'], s.loc['L3','B_hits']
    c_l1, c_l3 = s.loc['L1','C_hits'], s.loc['L3','C_hits']
    dc_31 = int(c_l3 - c_l1)
    db_31 = int(b_l3 - b_l1)
    pw_rows.append({
        'Dataset': label, 'Comparison': 'L3 vs L1',
        'Hypothesis': 'L3 captures more tissue biology',
        'Metric': 'ΔC_31 (primary) + ΔB_31 (secondary)',
        'L_low': int(c_l1), 'L_high': int(c_l3),
        'Delta': dc_31,
        'Delta_B': db_31,
        'Confirmed': dc_31 > 0,
        'B_direction': '↑' if db_31 > 0 else ('↓' if db_31 < 0 else '='),
    })

pw_df = pd.DataFrame(pw_rows)

# Summary counts
for comp in ['L2 vs L1', 'L3 vs L2', 'L3 vs L1']:
    sub = pw_df[pw_df['Comparison'] == comp]
    n_yes = sub['Confirmed'].sum()

piv = pw_df.pivot_table(
    index='Dataset', columns='Comparison',
    values='Delta', aggfunc='first',
)[['L2 vs L1', 'L3 vs L2', 'L3 vs L1']]

# Add ΔB_31 annotation column
db31 = pw_df[pw_df['Comparison']=='L3 vs L1'].set_index('Dataset')
piv['ΔB(L3−L1)'] = db31['Delta_B']
piv['B dir.'] = db31['B_direction']
display(piv)

pw_df.to_csv(op.join(ORA_DIR, 'pairwise_comparisons.csv'), index=False)

# ---

# Cell 18

labels = [DS_LABELS[d] for d in loaded_order]
n = len(loaded_order)
y = np.arange(n)
h = 0.6

fig, ax = plt.subplots(figsize=(10, max(5, n * 0.55)))

c_cumul = {}
b_cumul = {}

# C_hits: stacked leftward
left_c = np.zeros(n)
for L in ['L1', 'L2', 'L3']:
    c_vals = np.array([summary_df[(summary_df['dataset_id']==d)&(summary_df['query']==L)]['C_hits'].values[0]
                        for d in loaded_order], dtype=float)
    ax.barh(y, -c_vals, h, left=-left_c, color=L_COLORS_C[L], edgecolor='none',
            label=f'{L} (C hits)')
    left_c += c_vals
    c_cumul[L] = -left_c.copy()

# B_hits: stacked rightward
left_b = np.zeros(n)
for L in ['L1', 'L2', 'L3']:
    b_vals = np.array([summary_df[(summary_df['dataset_id']==d)&(summary_df['query']==L)]['B_hits'].values[0]
                        for d in loaded_order], dtype=float)
    ax.barh(y, b_vals, h, left=left_b, color=L_COLORS_B[L], edgecolor='none',
            label=f'{L} (B hits)')
    left_b += b_vals
    b_cumul[L] = left_b.copy()

def ribbon_flow(ax, y_centers, x_outer, x_inner, bar_h, color, alpha=0.15):
    """Ribbon that exactly covers one bar segment and smoothly flows between them."""
    hh = bar_h / 2
    nn = len(y_centers)
    y_pts = []; x_pts_outer = []; x_pts_inner = []
    for i in range(nn):
        y_pts.append(y_centers[i] - hh)
        x_pts_outer.append(x_outer[i])
        x_pts_inner.append(x_inner[i])
        y_pts.append(y_centers[i] + hh)
        x_pts_outer.append(x_outer[i])
        x_pts_inner.append(x_inner[i])
        if i < nn - 1:
            gap_y = np.linspace(y_centers[i] + hh, y_centers[i+1] - hh, 40)
            t = np.linspace(0, 1, 40)
            smooth_t = t * t * (3 - 2 * t)
            gap_outer = x_outer[i] + (x_outer[i+1] - x_outer[i]) * smooth_t
            gap_inner = x_inner[i] + (x_inner[i+1] - x_inner[i]) * smooth_t
            y_pts.extend(gap_y.tolist())
            x_pts_outer.extend(gap_outer.tolist())
            x_pts_inner.extend(gap_inner.tolist())
    y_pts = np.array(y_pts)
    x_pts_outer = np.array(x_pts_outer)
    x_pts_inner = np.array(x_pts_inner)
    ax.fill_betweenx(y_pts, x_pts_inner, x_pts_outer, color=color, alpha=alpha, zorder=0)
    
# L3 ribbon on left: from L2 inner edge to L3 outer edge
ribbon_flow(ax, y, c_cumul['L3'], c_cumul['L2'], h, L_COLORS_C['L3'], alpha=0.15)
# L2 ribbon on right: from L1 inner edge to L2 outer edge
ribbon_flow(ax, y, b_cumul['L2'], b_cumul['L1'], h, L_COLORS_B['L2'], alpha=0.15)

ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=12, fontweight='bold')

xmin_abs = abs(ax.get_xlim()[0])
xmax_abs = abs(ax.get_xlim()[1])
xbound = max(xmin_abs, xmax_abs) * 1.08
ax.set_xlim(-xbound, xbound)

ax.text(-xbound * 0.5, -1.2, 'C hits (tissue)', ha='center', fontsize=13, fontweight='bold', fontstyle='italic')
ax.text(xbound * 0.5, -1.2, 'B hits (regulation)', ha='center', fontsize=13, fontweight='bold', fontstyle='italic')

xticks = ax.get_xticks()
ax.set_xticklabels([f'{abs(int(t))}' for t in xticks])
ax.set_xlabel('Number of enriched terms', fontsize=13, fontweight='bold')

ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), ncol=1,
          prop={'weight': 'bold', 'size': 11})
for label in ax.get_yticklabels(): label.set_fontweight('bold')
for label in ax.get_xticklabels(): label.set_fontweight('bold')
ax.invert_yaxis()
plt.tight_layout()
plt.savefig(op.join(ORA_DIR, 'plot1_bidirectional_CB.pdf'))
plt.show()
# ---

# Cell 19
# Mild colormap: trim deep red/blue ends
full_cmap = plt.cm.get_cmap('RdBu_r', 256)
CMAP = LinearSegmentedColormap.from_list(
    'RdBu_r_mild',
    full_cmap(np.linspace(0.15, 0.85, 256))
)

hm_data = []
for ds_id in loaded_order:
    s = summary_df[summary_df['dataset_id'] == ds_id].set_index('query')
    hm_data.append({
        'Dataset': DS_LABELS[ds_id],
        'dB_21':    int(s.loc['L2','B_hits'] - s.loc['L1','B_hits']),
        'dC_32':    int(s.loc['L3','C_hits'] - s.loc['L2','C_hits']),
        'dC_31':    int(s.loc['L3','C_hits'] - s.loc['L1','C_hits']),
        'dB_pp_21': round((s.loc['L2','B_frac'] - s.loc['L1','B_frac']) * 100, 2),
        'dC_pp_32': round((s.loc['L3','C_frac'] - s.loc['L2','C_frac']) * 100, 2),
        'dC_pp_31': round((s.loc['L3','C_frac'] - s.loc['L1','C_frac']) * 100, 2),
    })
hm_df = pd.DataFrame(hm_data).set_index('Dataset')

count_cols = ['dB_21', 'dC_32', 'dC_31']
pp_cols = ['dB_pp_21', 'dC_pp_32', 'dC_pp_31']
all_counts = hm_df[count_cols].values.flatten()
all_pp = hm_df[pp_cols].values.flatten()
vabs_counts = max(np.abs(all_counts).max(), 1)
vabs_pp = max(np.abs(all_pp).max(), 0.5)

LINTHRESH_COUNTS = 10
LINTHRESH_PP = 2

fig, axes = plt.subplots(2, 3, figsize=(15, max(7, len(loaded_order)*0.8)))

panels = [
    (0, 0, 'dB_21',    '(A) L2 vs L1: ΔB hits\n(regulation gain)',   '.0f', vabs_counts, LINTHRESH_COUNTS),
    (0, 1, 'dC_32',    '(B) L3 vs L2: ΔC hits\n(tissue gain)',       '.0f', vabs_counts, LINTHRESH_COUNTS),
    (0, 2, 'dC_31',    '(C) L3 vs L1: ΔC hits\n(tissue gain)',       '.0f', vabs_counts, LINTHRESH_COUNTS),
    (1, 0, 'dB_pp_21', '(A) L2 vs L1: ΔB pp',                       '.1f', vabs_pp, LINTHRESH_PP),
    (1, 1, 'dC_pp_32', '(B) L3 vs L2: ΔC pp',                       '.1f', vabs_pp, LINTHRESH_PP),
    (1, 2, 'dC_pp_31', '(C) L3 vs L1: ΔC pp',                       '.1f', vabs_pp, LINTHRESH_PP),
]

for row, col, colname, title, fmt, vabs, linthresh in panels:
    ax = axes[row, col]
    data = hm_df[[colname]]
    norm = SymLogNorm(linthresh=linthresh, linscale=1, vmin=-vabs, vmax=vabs)
    sns.heatmap(data, annot=True, annot_kws={'fontsize': 11, 'fontweight': 'bold', 'color': 'black'}, fmt=fmt, center=0,
                cmap=CMAP, norm=norm,
                linewidths=0.5, linecolor='white',
                ax=ax, cbar_kws={'shrink': 0.5})
    ax.set_title(title, fontsize=11, fontweight='bold')
    for sp in ax.spines.values(): sp.set_visible(True)
    ax.set_ylabel('')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=10, fontweight='bold')
    if row == 0:
        ax.set_xticklabels([])

plt.tight_layout()
plt.savefig(op.join(ORA_DIR, 'plot3_pairwise_heatmaps.pdf'))
plt.show()
# ---

# Cell 20

raw_data = {s: [] for s in ['A', 'B', 'C']}
for ds_id in loaded_order:
    s = summary_df[summary_df['dataset_id'] == ds_id].set_index('query')
    for setname in ['A', 'B', 'C']:
        raw_data[setname].append({
            'Dataset': DS_LABELS[ds_id],
            'L1': int(s.loc['L1', f'{setname}_hits']),
            'L2': int(s.loc['L2', f'{setname}_hits']),
            'L3': int(s.loc['L3', f'{setname}_hits']),
        })

def mild_cmap(base_name):
    full = plt.cm.get_cmap(base_name, 256)
    return LinearSegmentedColormap.from_list(
        base_name + '_mild', full(np.linspace(0.05, 0.75, 256)))

fig, axes = plt.subplots(1, 3, figsize=(14, max(4, len(loaded_order) * 0.45)))

set_info = [
    ('A', 'Set A: Generic Cellular', mild_cmap('Reds')),
    ('B', 'Set B: Regulation', mild_cmap('Blues')),
    ('C', 'Set C: Tissue / Condition', mild_cmap('Greens')),
]

LINTHRESH_RAW = 15

for ax, (setname, title, cmap) in zip(axes, set_info):
    df = pd.DataFrame(raw_data[setname]).set_index('Dataset')
    vmax = max(df.values.max(), 1)
    norm = SymLogNorm(linthresh=LINTHRESH_RAW, linscale=1, vmin=0, vmax=vmax)
    sns.heatmap(df, annot=True, annot_kws={'fontsize': 11, 'fontweight': 'bold', 'color': 'black'}, fmt='.0f', cmap=cmap,
                norm=norm,
                linewidths=0.5, linecolor='white',
                ax=ax, cbar_kws={'shrink': 0.5})
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.set_ylabel('')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=10, fontweight='bold')
    ax.set_xticklabels(['L1', 'L2', 'L3'], fontsize=12, fontweight='bold')
    for sp in ax.spines.values(): sp.set_visible(True)

plt.tight_layout()
plt.savefig(op.join(ORA_DIR, 'plot3a_raw_ABC_hits.pdf'))
plt.show()

# ---

# Cell 22
DS_LV = 'Heart - Left Ventricle'
DS_NF = 'NF'


# E0: Input gene list overlap
LV_GENE_FILE = './results_gprofiler/Heart - Left Ventricle_all_lists_combined.txt'
NF_GENE_FILE = './results_magnet_gprofiler/NF_all_lists_combined.txt'

def parse_gene_lists(filepath):
    lists = {}; current = None
    with open(filepath, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                current = line[1:].strip(); lists[current] = []
            elif current is not None:
                lists[current].extend(line.split())
    return {k: set(v) for k, v in lists.items()}

def jaccard(a, b):
    inter = len(a & b); union = len(a | b)
    return inter, union, inter / union if union else np.nan

try:
    lv_lists = parse_gene_lists(LV_GENE_FILE)
    nf_lists = parse_gene_lists(NF_GENE_FILE)
    gene_rows = []
    for name in ['L1_Canonical', 'L2_Source_Isoform', 'L3_Target_Isoform']:
        lv_set = lv_lists.get(name, set()); nf_set = nf_lists.get(name, set())
        inter, union, j = jaccard(lv_set, nf_set)
        gene_rows.append({'list': name, 'LV_n': len(lv_set), 'NF_n': len(nf_set),
                          'intersection': inter, 'union': union, 'Jaccard': round(j, 4)})
    gene_overlap_df = pd.DataFrame(gene_rows)
    display(gene_overlap_df)
    gene_overlap_df.to_csv(op.join(ORA_DIR, 'E0_gene_list_overlap.csv'), index=False)
except FileNotFoundError as e:

# E1: g:Profiler term overlap
SCOPES_PER_LEVEL = {'L1': ['all', 'Set_A'], 'L2': ['all', 'Set_B'], 'L3': ['all', 'Set_C']}
SCOPE_FILTERS = {
    'all':   lambda df: df,
    'Set_A': lambda df: df[df['kw_set'] == 'A'],
    'Set_B': lambda df: df[df['kw_set'] == 'B'],
    'Set_C': lambda df: df[df['kw_set'] == 'C'],
}
overlap_rows = []
for L in ['L1', 'L2', 'L3']:
    lv_base = SIG[(SIG['dataset_id']==DS_LV)&(SIG['query']==L)]
    nf_base = SIG[(SIG['dataset_id']==DS_NF)&(SIG['query']==L)]
    for scope_name in SCOPES_PER_LEVEL[L]:
        lv = SCOPE_FILTERS[scope_name](lv_base); nf = SCOPE_FILTERS[scope_name](nf_base)
        s_lv = set(lv['term_id']); s_nf = set(nf['term_id'])
        inter, union, j = jaccard(s_lv, s_nf)
        overlap_rows.append({'list': f'{L} / {scope_name}', 'LV_n': len(s_lv),
                             'NF_n': len(s_nf), 'intersection': inter, 'union': union,
                             'Jaccard': round(j, 4)})
overlap_df = pd.DataFrame(overlap_rows)
display(overlap_df)
overlap_df.to_csv(op.join(ORA_DIR, 'E1_term_overlap.csv'), index=False)

# E2: Shared enriched terms
FOCUS = [('L1', 'A', 'Set A (Generic)'), ('L2', 'B', 'Set B (Regulation)'),
         ('L3', 'C', 'Set C (Tissue)')]
for L, setlabel, setname in FOCUS:
    lv_terms = SIG[(SIG['dataset_id']==DS_LV)&(SIG['query']==L)&(SIG['kw_set']==setlabel)].copy()
    nf_terms = SIG[(SIG['dataset_id']==DS_NF)&(SIG['query']==L)&(SIG['kw_set']==setlabel)].copy()
    shared_ids = set(lv_terms['term_id']) & set(nf_terms['term_id'])
    term_rows = []
    for tid in shared_ids:
        lr = lv_terms[lv_terms['term_id']==tid]; nr = nf_terms[nf_terms['term_id']==tid]
        term_rows.append({'term_id': tid, 'term_name': lr.iloc[0]['term_name'],
                          'p_adj_LV': lr.iloc[0]['p_adj'], 'p_adj_NF': nr.iloc[0]['p_adj']})
    if term_rows:
        term_df = pd.DataFrame(term_rows).sort_values('p_adj_LV')
        display(term_df)
        term_df.to_csv(op.join(ORA_DIR, f'E2_{L}_{setlabel}_terms.csv'), index=False)

# E3: Rank analysis
from scipy.stats import spearmanr as _spearmanr
import matplotlib.cm as cm

LEVEL_KEYS = [('L1_Canonical','L1'), ('L2_Source_Isoform','L2'), ('L3_Target_Isoform','L3')]
rank_dfs = {}

if 'lv_lists' in dir() and 'nf_lists' in dir():
    for list_key, level_label in LEVEL_KEYS:
        lv_ordered = list(lv_lists.get(list_key, set()))
        nf_ordered = list(nf_lists.get(list_key, set()))
        lv_ordered.sort()  # deterministic order
        nf_ordered.sort()
        shared_genes = set(lv_ordered) & set(nf_ordered)
        if len(shared_genes) == 0:
            continue

        lv_rank = {g: i+1 for i, g in enumerate(lv_ordered)}
        nf_rank = {g: i+1 for i, g in enumerate(nf_ordered)}
        rows = [{'gene': g, 'rank_LV': lv_rank[g], 'rank_NF': nf_rank[g]} for g in shared_genes]
        rdf = pd.DataFrame(rows).sort_values('rank_LV')
        rdf['mean_rank'] = (rdf['rank_LV'] + rdf['rank_NF']) / 2
        rank_dfs[level_label] = rdf
        rdf.to_csv(op.join(ORA_DIR, f'E3_shared_{level_label}_gene_ranks.csv'), index=False)

        pass  # individual bump charts replaced by combined panel below

    if rank_dfs:
        bump_colors = {'L1': '#CB181D', 'L2': '#2171B5', 'L3': '#238B45'}
        
        n_panels = len(rank_dfs)
        fig, axes = plt.subplots(1, n_panels, figsize=(4 * n_panels, 6))
        if n_panels == 1: axes = [axes]

        for ax, (ll, rdf) in zip(axes, rank_dfs.items()):
            max_rank = max(rdf['rank_LV'].max(), rdf['rank_NF'].max())
            col_dark = bump_colors[ll]

            for _, row in rdf.iterrows():
                ax.plot([0, 1], [row['rank_LV'], row['rank_NF']],
                        color=col_dark, alpha=0.6, lw=0.8, solid_capstyle='round')

            ax.set_xlim(-0.05, 1.05)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(['Heart - Left Ventricle', 'NF'], fontsize=12, fontweight='bold')
            ax.set_ylabel('Rank' if ax == axes[0] else '', fontsize=13, fontweight='bold')
            ax.invert_yaxis()
            for sp in ax.spines.values(): sp.set_visible(True)
            for label in ax.get_yticklabels(): label.set_fontweight('bold')

            corr = _spearmanr(rdf['rank_LV'], rdf['rank_NF']).correlation
            ax.text(0.5, 1.02, f'{ll} (n={len(rdf)}, \u03C1={corr:.2f})',
                    transform=ax.transAxes, ha='center', fontsize=12, fontweight='bold')

        plt.tight_layout()
        plt.savefig(op.join(ORA_DIR, 'E3a_slope_panel.pdf'))
        plt.show()

            # ---- E3b: Combined rank concordance scatter ----
    if rank_dfs:
        fig, ax = plt.subplots(figsize=(6, 6))
        level_colors = {'L1': '#7788B4', 'L2': '#B0C2E4', 'L3': '#D9BCCD'}
        level_markers = {'L1': 'o', 'L2': 's', 'L3': '^'}

        max_rank = 0
        for ll, rdf in rank_dfs.items():
            corr = _spearmanr(rdf['rank_LV'], rdf['rank_NF']).correlation
            ax.scatter(rdf['rank_LV'], rdf['rank_NF'],
                       s=12, alpha=0.4, c=level_colors[ll], marker=level_markers[ll],
                       edgecolors='#333', linewidths=0.2,
                       label=f'{ll} (n={len(rdf)}, ρ={corr:.2f})')
            max_rank = max(max_rank, rdf['rank_LV'].max(), rdf['rank_NF'].max())

        ax.plot([0, max_rank], [0, max_rank], color='#999', ls='--', lw=0.8)
        ax.set_xlabel('Rank in Heart - Left Ventricle', fontsize=13, fontweight='bold')
        ax.set_ylabel('Rank in NF', fontsize=13, fontweight='bold')
        ax.set_aspect('equal')
        for label in ax.get_xticklabels(): label.set_fontweight('bold')
        for label in ax.get_yticklabels(): label.set_fontweight('bold')
        ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1),
                  prop={'weight': 'bold', 'size': 11})
        plt.tight_layout()
        plt.savefig(op.join(ORA_DIR, 'E3b_rank_concordance_all.pdf'))
        plt.show()

    for list_key, level_label in LEVEL_KEYS:
        if level_label not in rank_dfs:
            continue
        rdf = rank_dfs[level_label]
        lv_ordered = list(lv_lists.get(list_key, []))
        nf_ordered = list(nf_lists.get(list_key, []))
        lv_rank = {g: i+1 for i, g in enumerate(lv_ordered)}
        nf_rank = {g: i+1 for i, g in enumerate(nf_ordered)}
        shared_genes = set(rdf['gene'])

        n_perm = 1000
        rng = np.random.default_rng(42)
        obs_lv = np.array([lv_rank[g] for g in shared_genes])
        obs_nf = np.array([nf_rank[g] for g in shared_genes])
        obs_r = _spearmanr(obs_lv, obs_nf).correlation

        null_corrs = []
        for _ in range(n_perm):
            shuffled_nf = rng.permutation(obs_nf)
            null_corrs.append(_spearmanr(obs_lv, shuffled_nf).correlation)
        null_corrs = np.array(null_corrs)
        pval = (np.sum(null_corrs >= obs_r) + 1) / (len(null_corrs) + 1)

        fig, ax = plt.subplots(figsize=(6, 3.5))
        ax.hist(null_corrs, bins=40, color='#AEC7E8', alpha=0.6,
                edgecolor='white', linewidth=0.3, label='Null distribution')
        ax.axvline(obs_r, color='#C44E52', lw=2, ls='--',
                   label=f'Observed ρ={obs_r:.2f} (p={pval:.3f})')
        ax.set_xlabel('')
        ax.set_ylabel('Count', fontsize=13, fontweight='bold')
        sig = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'
        ax.text(0.95, 0.92, f'p = {pval:.4f} ({sig})',
                transform=ax.transAxes, ha='right', va='top', fontsize=11, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          edgecolor='#CCC', alpha=0.9))
        for label in ax.get_xticklabels(): label.set_fontweight('bold')
        for label in ax.get_yticklabels(): label.set_fontweight('bold')
        ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1),
                  prop={'weight': 'bold', 'size': 11})
        plt.tight_layout()
        plt.savefig(op.join(ORA_DIR, f'E3c_null_{level_label}.pdf'))
        plt.show()


# ---

# Cell 24
for ds_id in loaded_order:
    terms = sorted(C_L3_ONLY.get(ds_id, set()))
    if not terms:
        continue
    for t in terms[:12]:
    if len(terms) > 12:

# ---

# Cell 26
# L3-only tissue terms per dataset
for ds_id in loaded_order:
    terms = sorted(C_L3_ONLY.get(ds_id, set()))
    if terms:
        # Get full info from SIG
        t_df = SIG[(SIG['dataset_id']==ds_id) & (SIG['query']=='L3') &
                    (SIG['kw_set']=='C') &
                    (SIG['term_norm'].isin(C_L3_ONLY[ds_id]))][
            ['term_id','term_name','p_adj','is_highlighted']
        ].drop_duplicates('term_name').sort_values('p_adj')
        safe = ds_id.replace(' ','_').replace('-','')
        t_df.to_csv(op.join(ORA_DIR, f'{safe}_L3_only_C_terms.tsv'),
                    sep='\t', index=False)

# Keywords
kw_rows = []
for kw in SET_A:
    kw_rows.append({'set': 'A_generic', 'dataset': 'ALL', 'keyword': kw})
for kw in SET_B:
    kw_rows.append({'set': 'B_regulation', 'dataset': 'ALL', 'keyword': kw})
for ds, kws in SET_C.items():
    for kw in kws:
        kw_rows.append({'set': 'C_tissue', 'dataset': ds, 'keyword': kw})
pd.DataFrame(kw_rows).to_csv(op.join(ORA_DIR, 'keywords_3set.csv'), index=False)

for f in sorted(glob.glob(op.join(ORA_DIR, '*'))):

# ---

# Cell 27
out_path = op.join(ORA_DIR, 'highlighted_top15_per_level.xlsx')

with pd.ExcelWriter(out_path, engine='openpyxl') as writer:
    for ds_id in loaded_order:
        for L in ['L1', 'L2', 'L3']:
            sub = SIG[(SIG['dataset_id']==ds_id) & (SIG['query']==L) & (SIG['is_highlighted']==True)]
            top = sub.sort_values('p_adj').head(15)[['term_id','term_norm','p_adj','kw_set']]
            top = top.reset_index(drop=True)
            top.index += 1
            top.columns = ['Term ID', 'Term Name', 'Adj. p-value', 'Keyword Set']

            safe = DS_LABELS[ds_id].replace(' ','_').replace('-','')
            sheet_name = f'{safe}_{L}'[:31]
            top.to_excel(writer, sheet_name=sheet_name)


# ---