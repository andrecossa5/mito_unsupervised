"""
Visualization of clones and samples classification performances.
"""

# Code
import sys
import os
import re
import pickle
from mito_utils.preprocessing import *
from mito_utils.diagnostic_plots import *
from mito_utils.heatmaps_plots import *
from mito_utils.utils import *
from mito_utils.plotting_base import *
from matplotlib.gridspec import GridSpec
import matplotlib


##

# Args
path_main = '/Users/IEO5505/Desktop/mito_bench/'

# Set paths
path_data = os.path.join(path_main, 'data')
path_output = os.path.join(path_main, 'results', 'unsupervised_clones', 'output')
path_report = os.path.join(path_main, 'results', 'unsupervised_clones', 'reports')
path_viz = '/Users/IEO5505/Desktop/MI_TO/images_DIPA/images'

# Create report
L = []
pickles = [ x for x in os.listdir(path_output) if bool(re.search('.pickle', x)) ]
for p in pickles:
    with open(os.path.join(path_output, p), 'rb') as f:
        d = pickle.load(f)
    L.append(d['df_performance'])

# Concat and save
pd.concat(L).to_csv(os.path.join(path_report, 'report_unsupervised.csv'))


##


# Load report
df = pd.read_csv(os.path.join(path_report, 'report_unsupervised.csv'), index_col=0)


##


############## Extended summary

# Full summary
df.describe().to_excel(
    os.path.join(path_report, 'full_aggregate_f1.xlsx')
)

############## Ranking and technical summary

# Plot
metric_colors = {'Adjusted Rand Index (ARI)':'#D8D3CA', 'Normalized Mutual Information (NMI)':'#9A9588'}

# Data
df_ = (
    df.loc[:, ['NMI', 'ARI', 'model', 'filtering']]
    .assign(
        model_short=lambda x: np.where(x['model'].str.startswith('leiden'), 'leiden', 'vireoSNP'),
        method=lambda x: x['filtering'] + '_' + x['model_short']
    )
    .drop(columns=['model', 'filtering', 'model_short'])
)
method_order = (
    df_.groupby('method')
    .mean()
    .sort_values('NMI', ascending=False)
    .index
)
df_ = df_.melt(id_vars='method', var_name='metric')
df_['metric'] = np.select(
    [df_['metric'] == 'ARI', df_['metric'] == 'NMI'],
    ['Adjusted Rand Index (ARI)', 'Normalized Mutual Information (NMI)']
)

# Fig
fig, ax = plt.subplots(figsize=(7, 7))
# gs = GridSpec(3, 2, figure=fig, width_ratios=[2.5,3], height_ratios=[0.5,1,2])

# Ax
sns.barplot(
    data=df_, x='value', y='method', hue='metric', 
    order=method_order,
    orient='h', 
    palette=metric_colors.values(), 
    errcolor='k',
    errwidth=.8,
    capsize=.05,
    #width=0.1,
    saturation=.7,
    dodge=1,
    ax=ax,
    edgecolor='k',
    linewidth=.5
)
format_ax(xlabel='Metric value', ylabel='', ax=ax, title='Method')
ax.spines[['right', 'top', 'left']].set_visible(False)
ax.legend([], [], frameon=False)
add_legend(
    label='Metric',
    colors=metric_colors, 
    ax=ax,
    loc='lower right',
    bbox_to_anchor=(.98, .02),
    ncols=1,
    artists_size=9,
    label_size=10,
    ticks_size=9
)

# Save
fig.tight_layout()
fig.savefig(os.path.join(path_viz, 'methods_ranking.png'))

##############


##


############## Top 3 models by sample

# Data
df_ = (
    df.loc[:, ['NMI', 'ARI', 'model', 'filtering', 'sample']]
    .assign(
        model_short=lambda x: np.where(x['model'].str.startswith('leiden'), 'leiden', 'vireoSNP'),
        method=lambda x: x['filtering'] + '_' + x['model_short']
    )
    .drop(columns=['model', 'filtering', 'model_short'])
    .groupby('sample')
    .apply(lambda x: x.sort_values('NMI', ascending=False).head(3))
    .reset_index(drop=True)
)

# Method order
method_order = (
    df_.groupby('method')
    .median()
    .sort_values('NMI', ascending=False)
    .index
)

# Sample order
sample_order = (
    df_.groupby('sample')
    .median()
    .sort_values('NMI', ascending=False)
    .index
)

# Colors
_ = sns.color_palette('inferno', n_colors=10)
_ = [_[9], _[5], _[2]]
methods_color = { m : _[i] for i, m in enumerate(method_order)}

# Plot
fig, axs = plt.subplots(1,2,figsize=(10, 5))

# NMI
strip(df_, 'sample', 'NMI', by='method', c=methods_color, order=sample_order, s=7, ax=axs[0])
axs[0].set_ylim((0.3, 1))
axs[0].spines[['right', 'top']].set_visible(False)
axs[0].hlines(0.5, xmin=0, xmax=df_['sample'].unique().size-0.8, linestyles='dashed', color='k')
medians = df_.groupby('sample')['NMI'].median()
std = df_.groupby('sample')['NMI'].std()
for j, m in enumerate(sample_order):
    axs[0].errorbar(j, medians[m], yerr=std[m], c='k')

# ARI
strip(df_, 'sample', 'ARI', by='method', c=methods_color, order=sample_order, s=8, ax=axs[1])
axs[1].set_ylim((0, 1))
axs[1].spines[['right', 'top']].set_visible(False)
add_legend(
    ax=axs[1],
    label='Method',
    colors=methods_color,
    loc='upper right',
    artists_size=8,
    ticks_size=8,
    label_size=10,
    bbox_to_anchor=(1,1),
    ncols=1
)
axs[1].hlines(0.5, xmin=0, xmax=df_['sample'].unique().size-0.8, linestyles='dashed', color='k')
medians = df_.groupby('sample')['ARI'].median()
std = df_.groupby('sample')['ARI'].std()
for j, m in enumerate(sample_order):
    axs[1].errorbar(j, medians[m], yerr=std[m], c='k')
    
# Save
fig.suptitle('Top 3 methods')
fig.savefig(os.path.join(path_viz, 'top_3_unsupervised.png'))

############## 