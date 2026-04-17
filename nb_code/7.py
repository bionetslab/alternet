# Cell 1
import os, os.path as op
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
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.labelsize': 12, 'ytick.labelsize': 12,
    'legend.fontsize': 11, 'legend.frameon': True,
    'savefig.bbox': 'tight', 'savefig.dpi': 300,
})

FIG_DIR = './edge_count_figs/'
os.makedirs(FIG_DIR, exist_ok=True)

DATASETS = [
    ('MAGNet', 'DCM'), ('MAGNet', 'HCM'), ('MAGNet', 'NF'),
    ('GTEx', 'Bladder'), ('GTEx', 'Heart - Left Ventricle'),
    ('GTEx', 'Kidney - Cortex'), ('GTEx', 'Liver'),
    ('GTEx', 'Pancreas'), ('GTEx', 'Spleen'), ('GTEx', 'Stomach'),
]
DS_LABELS = [g for _, g in DATASETS]
n_ds = len(DATASETS)
n_mag = 3

PART1 = {'MAGNet': './results_magnet_net_infer/', 'GTEx': './results_gtex_net_infer/'}
PART2 = {'MAGNet': './results_magnet_edge_cate/', 'GTEx': './results_gtex_edge_cate/'}
PART3 = {'MAGNet': './results_magnet_plau_filter/', 'GTEx': './results_gtex_plau_filter/'}

def prefix(dtype, group):
    return group

def count_rows(filepath):
    return len(pd.read_csv(filepath, sep='\t'))

def count_by_category(filepath, cat_col):
    df = pd.read_csv(filepath, sep='\t')
    return df[cat_col].value_counts().to_dict()

# ---

# Cell 3
a1_data = {'Net1': [], 'Net2': [], 'Net3': []}
for dtype, group in DATASETS:
    p = prefix(dtype, group)
    r = PART1[dtype]
    a1_data['Net1'].append(count_rows(op.join(r, f'{p}_canonical_raw.tsv')))
    a1_data['Net2'].append(count_rows(op.join(r, f'{p}_as_aware_source_raw.tsv')))
    a1_data['Net3'].append(count_rows(op.join(r, f'{p}_fully_as_aware_raw.tsv')))

a2_data = {'Net1': [], 'Net2': [], 'Net3': []}
for dtype, group in DATASETS:
    p = prefix(dtype, group)
    r = PART2[dtype]
    a2_data['Net1'].append(count_rows(op.join(r, f'{p}_net1_filtered.tsv')))
    a2_data['Net2'].append(count_rows(op.join(r, f'{p}_net2_filtered.tsv')))
    a2_data['Net3'].append(count_rows(op.join(r, f'{p}_net3_filtered.tsv')))

B_CONFIG = {
    'B1': {
        'cat_file': '{p}_set_a_source_resolution.tsv',
        'plaus_file': '{p}_set_a_plausible.tsv',
        'cat_col': 'source_category',
        'categories': ['source_isoform_specific', 'source_gene_specific',
                       'source_equivalent', 'source_ambiguous'],
        'colors': {'source_isoform_specific': '#b7282e', 'source_gene_specific': '#d16d5b',
                   'source_equivalent': '#dc917b', 'source_ambiguous': '#eabaa1'},
        'label': 'Set A',
    },
    'B2': {
        'cat_file': '{p}_set_b_target_resolution.tsv',
        'plaus_file': '{p}_set_b_plausible.tsv',
        'cat_col': 'target_category',
        'categories': ['target_isoform_specific', 'target_gene_specific',
                       'target_equivalent', 'target_ambiguous'],
        'colors': {'target_isoform_specific': '#313772', 'target_gene_specific': '#326db6',
                   'target_equivalent': '#478ecc', 'target_ambiguous': '#75b5dc'},
        'label': 'Set B',
    },
    'B3': {
        'cat_file': '{p}_set_c_sf_splicing.tsv',
        'plaus_file': '{p}_set_c_plausible.tsv',
        'cat_col': 'sf_category',
        'categories': ['sf_splicing_supported_specific', 'sf_splicing_supported_broad',
                       'sf_expression_associated', 'sf_ambiguous'],
        'colors': {'sf_splicing_supported_specific': '#376439', 'sf_splicing_supported_broad': '#669877',
                   'sf_expression_associated': '#81b095', 'sf_ambiguous': '#a4cbb7'},
        'label': 'Set C',
    },
    'B4': {
        'cat_file': '{p}_set_d_tfsf_joint_ambiguous.tsv',
        'plaus_file': '{p}_set_d_plausible.tsv',
        'cat_col': 'tfsf_category',
        'categories': ['tfsf_joint', 'tfsf_ambiguous'],
        'colors': {'tfsf_joint': '#4e4e4e', 'tfsf_ambiguous': '#8d8d8d'},
        'label': 'Set D',
    },
}

B_DATA = {}
for bname, cfg in B_CONFIG.items():
    cat_counts = []
    plaus_counts = []
    for dtype, group in DATASETS:
        p = prefix(dtype, group)
        r2 = PART2[dtype]
        r3 = PART3[dtype]
        cat_counts.append(count_by_category(op.join(r2, cfg['cat_file'].format(p=p)), cfg['cat_col']))
        plaus_counts.append(count_by_category(op.join(r3, cfg['plaus_file'].format(p=p)), cfg['cat_col']))
    B_DATA[bname] = {'cat': cat_counts, 'plaus': plaus_counts}

# ---

# Cell 5
NET_COLORS = {'Net1': '#7788B4', 'Net2': '#B0C2E4', 'Net3': '#D9BCCD'}

def add_divider(ax, n_mag, n_total, zero_indexed=True):
    if 0 < n_mag < n_total:
        xpos = (n_mag - 0.5) if zero_indexed else (n_mag + 0.5)
        ax.axvline(x=xpos, color='#333', lw=0.8, ls=':', alpha=0.6)

for plot_label, data, ylabel_extra in [
    ('A1', a1_data, '(raw aggregated, Part 1)'),
    ('A2', a2_data, '(filtered, Part 2)'),
]:
    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(n_ds)
    w = 0.25
    for j, (net, col) in enumerate(NET_COLORS.items()):
        vals = data[net]
        ax.bar(x + (j - 1) * w, vals, w, color=col, edgecolor='none', label=net)
    ax.set_xticks(x)
    ax.set_xticklabels(DS_LABELS, rotation=45, ha='right', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of edges', fontsize=13, fontweight='bold')
    add_divider(ax, n_mag, n_ds)
    ax.legend(fontsize=9, loc='upper left', bbox_to_anchor=(1.01, 1),
              prop={'weight': 'bold', 'size': 11})
    for label in ax.get_yticklabels():
        label.set_fontweight('bold')
    ax.yaxis.get_offset_text().set_fontweight('bold')
    ax.yaxis.get_offset_text().set_fontsize(12)
    plt.tight_layout()
    plt.savefig(op.join(FIG_DIR, f'{plot_label}_edge_counts.pdf'))
    plt.show()

# ---

# Cell 6
for bname, cfg in B_CONFIG.items():
    bd = B_DATA[bname]
    cats = cfg['categories']
    colors = cfg['colors']

    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(n_ds)
    w = 0.3
    gap = 0.03

    bot_c = np.zeros(n_ds)
    for cat in cats:
        vals = np.array([bd['cat'][i].get(cat, 0) for i in range(n_ds)], dtype=float)
        ax.bar(x - w/2 - gap/2, vals, w, bottom=bot_c, color=colors.get(cat, '#999'),
               edgecolor='none')
        bot_c += vals

    bot_p = np.zeros(n_ds)
    for cat in cats:
        vals = np.array([bd['plaus'][i].get(cat, 0) for i in range(n_ds)], dtype=float)
        ax.bar(x + w/2 + gap/2, vals, w, bottom=bot_p, color=colors.get(cat, '#999'),
               edgecolor='none')
        bot_p += vals

    ax.set_xticks(x)
    ax.set_xticklabels(DS_LABELS, rotation=45, ha='right', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of edges', fontsize=13, fontweight='bold')
    add_divider(ax, n_mag, n_ds)

    handles = [Patch(fc=colors[c], label=c) for c in cats]
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.01, 1),
              prop={'weight': 'bold', 'size': 11}, handlelength=1.5, handletextpad=0.5)
    for label in ax.get_yticklabels():
        label.set_fontweight('bold')
    plt.tight_layout()
    plt.savefig(op.join(FIG_DIR, f'{bname}_stacked.pdf'))
    plt.show()

# ---