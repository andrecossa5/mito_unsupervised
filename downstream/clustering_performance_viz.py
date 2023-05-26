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
path_viz = os.path.join(path_main, 'results', 'unsupervised_clones', 'viz_overall_performance')


os.listdir(path_viz)

# Create report
L = []
pickles = [ x for x in os.listdir(path_output) if bool(re.search('.pickle', x)) ]
for p in pickles:
    with open(os.path.join(path_output, p), 'rb') as f:
        d = pickle.load(f)
    L.append(d['performance_df'])

# Concat and save
pd.concat(L).to_csv(os.path.join(path_report, 'report_f1.csv'))


##


# Load report
clones = pd.read_csv(os.path.join(path_report, 'report_f1.csv'), index_col=0)


##


############## Extended summary, aggregate by sample

# Full summary
clones.describe().to_excel(
    os.path.join(path_report, 'full_aggregate_f1.xlsx')
)

# Agg by sample and run summary
(
    clones.assign(job=lambda x: x['filtering'] + '_' + x['dimred'] + '_' + x['model'] + '_' + x['tuning'])
    .groupby(['sample', 'job'])
    .agg('mean')
    .reset_index(level=1)
    .groupby('sample')
    .apply(lambda x: x.sort_values('f1', ascending=False))
).to_excel(
    os.path.join(path_report, 'sample_job_aggregate_f1.xlsx')
)

# Plot
fig, axs = plt.subplots(1,2, figsize=(13.5, 5), constrained_layout=True)

df_ = (
    clones
    .assign(job=lambda x: x['filtering'] + '_' + x['dimred'] + '_' + x['model'] + '_' + x['tuning'])
    .groupby(['sample', 'job'])
    .agg('mean')
)

# Top
box(df_.iloc[:,:6].melt(), 'variable', 'value', c='#E9E7E7', order=None, ax=axs[0])
strip(df_.query('AUCPR>=0.7').iloc[:,:6].melt(), 'variable', 'value', c='r', s=3, ax=axs[0])
strip(df_.query('AUCPR<0.7').iloc[:,:6].melt(), 'variable', 'value', c='#3b3838', s=2.5, ax=axs[0])

n_good = df_.query('AUCPR>=0.7').shape[0]
format_ax(axs[0], ylabel='value', xticks_size=10, rotx=90, title='Feature type')
add_legend(
    label='Model',
    colors={'AUCPR >= 0.7':'r', 'AUCPR < 0.7':'#3b3838'}, 
    ax=axs[0],
    loc='upper left',
    bbox_to_anchor=(1,1),
    ncols=1,
    label_size=12,
    ticks_size=10
)
axs[0].spines[['top', 'right']].set_visible(False)

# Bottom

# Colors
df_['feature_type'] = df_.reset_index(level=1)['job'].map(lambda x: '_'.join(x.split('_')[:-2])).values
feat_type_colors = create_palette(df_, 'feature_type', 'Spectral', saturation=.85)

# Axes
box(df_.iloc[:,:6].melt(), 'variable', 'value', c='#E9E7E7', order=None, ax=axs[1])
strip(
    df_.loc[:,
        [ 'accuracy', 'balanced_accuracy', 'precision',
        'recall', 'f1', 'AUCPR', 'feature_type']
    ]
    .melt(id_vars='feature_type'),
    'variable', 'value', by='feature_type', c=feat_type_colors, s=2.5, ax=axs[1]
)
format_ax(axs[1], xticks_size=10, rotx=90)
add_legend(
    label='Feature type',
    colors=feat_type_colors,
    ax=axs[1],
    loc='upper left',
    bbox_to_anchor=(1,1),
    ncols=1,
    label_size=12,
    ticks_size=10
)
axs[1].spines[['top', 'right']].set_visible(False)

# Save
fig.suptitle(f"{n_good}/{df_.shape[0]} ({n_good/df_.shape[0]*100:.2f}%) 'overall good' models (>=0.7 mean AUCPR, over all sample clones)")
fig.savefig(os.path.join(path_viz, 'overall_performance_f1.png'))

##############


##


############## Ranking and technical summary

# Plot
metric_colors = {'precision':'#F0D290', 'recall':'#DE834D', 'f1':'#A3423C', 'AUCPR':'#5C4040'}

fig = plt.figure(figsize=(8, 8))
# gs = GridSpec(3, 2, figure=fig, width_ratios=[2.5,3], height_ratios=[0.5,1,2])
# 
df_ = (
    clones
    .assign(job=lambda x: x['filtering'] + '|' + x['dimred'] + '|' + x['model'] + '|' + x['tuning'])
    .groupby(['sample', 'job'])
    .agg('mean')
    .reset_index(level=1)
)

df_[['filtering', 'dimred', 'model', 'tuning']] = df_['job'].str.split('|', expand=True)
df_ = df_.drop(columns=['job']).reset_index(drop=True)

# Feat type ranking plot
df_1 = (
    df_.loc[:, ['f1', 'precision', 'recall', 'AUCPR', 'filtering', 'dimred']]
    .assign(feature_type=lambda x: x['filtering'] + '_' + x['dimred'])
    .drop(columns=['dimred', 'filtering'])
)
feat_order = (
    df_1.groupby('feature_type')
    .mean()
    .sort_values('AUCPR', ascending=False)
    .index
)
df_1 = df_1.melt(id_vars='feature_type', var_name='metric')

# # Ax
sns.barplot(
    data=df_1, x='value', y='feature_type', hue='metric', 
    order=feat_order,
    orient='h', 
    palette=metric_colors.values(), 
    errcolor='k',
    errwidth=.8,
    capsize=.05,
    saturation=1,
    ax=ax,
    edgecolor='k',
    linewidth=.5
)
format_ax(xlabel='Metric value', ylabel='', ax=ax, title='Feature type')
ax.legend([], [], frameon=False)
ax.spines[['right', 'top', 'left']].set_visible(False)

##

# Model
# ax2 = fig.add_subplot(gs[:2, 1])
# box(
#     (
#         df_.loc[:, ['f1', 'precision', 'recall', 'AUCPR', 'model']]
#         .melt(id_vars='model', var_name='metric', value_name='score')
#     ), 
#     'model', 
#     'score',
#     by='metric', 
#     c=metric_colors, 
#     hue_order=metric_colors.keys(),
#     ax=ax2
# )
# format_ax(ax2, xticks_size=10, rotx=0, title='Model', ylabel='Metric value')
# add_legend(
#     label='Metric',
#     colors=metric_colors, 
#     ax=ax2,
#     loc='upper center',
#     bbox_to_anchor=(0.5, -0.15),
#     ncols=4,
#     artists_size=9,
#     label_size=10,
#     ticks_size=9
# )
# ax2.spines[['top', 'right']].set_visible(False)


# Save
fig.tight_layout()
fig.savefig(os.path.join(path_viz, 'feature_type_ranking.png'))
plt.show()


##############


##


############## Overall performance by clone

# Summary
df_ = clones.assign(feature_type=lambda x: x['filtering'] + '_' + x['dimred'])
clones_order = (
    df_.groupby('comparison')['f1']
    .agg('mean')
    .sort_values(ascending=False)
    .index
)
n_clones = df_['comparison'].unique().size
good_clones = df_.query('f1>=0.7')['comparison'].unique()
n_good_clones = good_clones.size
df_['clone_status'] = np.where(df_['comparison'].isin(good_clones), '>=1 good models', 'No good models')

# Good clones per sample
(
    df_.loc[:, ['comparison', 'clone_status', 'sample']]
    .drop_duplicates()
    .groupby('sample')
    .apply(lambda x: x['clone_status'].value_counts())
    .reset_index(level=1, name='count')
    .rename(columns={'level_1':'clone_status'})
    .assign(freq=lambda x: 
        x['count'] / x.reset_index().groupby('sample')['count'].transform('sum').values
    )
).to_excel(
    os.path.join(path_report, 'clones_status_aggregate_f1.xlsx')
)

##

# Plot
fig, axs = plt.subplots(2, 1, figsize=(10, 5), sharex=True, sharey=True)

# Topz
box(df_, 'comparison', 'f1', c='#E9E7E7', ax=axs[0], order=clones_order)
strip(df_.query('f1>=0.7'), 'comparison', 'f1', c='r', s=3, ax=axs[0], order=clones_order)
strip(df_.query('f1<0.7'), 'comparison', 'f1', c='#3b3838', s=1.5, ax=axs[0], order=clones_order)
format_ax(
    axs[0], 
    title=f"{n_good_clones}/{n_clones} ({n_good_clones/n_clones*100:.2f}%) clones with at least 1 'good' model", 
    ylabel='f1 score', 
    rotx=90, 
    xticks_size=5,
    xticks=[]
    )
add_legend(
    label='Model',
    colors={'f1 >= 0.7':'r', 'f1 < 0.7':'#3b3838'}, 
    ax=axs[0],
    loc='upper left',
    bbox_to_anchor=(1,1),
    ncols=1,
    label_size=10,
    ticks_size=8
)
axs[0].spines[['top', 'bottom', 'right']].set_visible(False)

# Bottom:
feat_type_colors = create_palette(df_, 'feature_type', 'Spectral', saturation=.85)
box(df_, 'comparison', 'f1', c='#E9E7E7', ax=axs[1], order=clones_order)
strip(df_, 'comparison', 'f1', by='feature_type', c=feat_type_colors, s=2, ax=axs[1], order=clones_order)
format_ax(
    axs[1],
    ylabel='f1 score', 
    rotx=90,
    xticks_size=5
)
add_legend(
    label='Feature type',
    colors=feat_type_colors, 
    ax=axs[1],
    loc='upper left',
    bbox_to_anchor=(1,1),
    ncols=1,
    label_size=10,
    ticks_size=8
)
axs[1].spines[['top', 'right']].set_visible(False)

# Save
fig.tight_layout()
fig.savefig(os.path.join(path_viz, 'performance_by_clone_f1.png'))

##############


##


############## Overall performance by clone 2: relationship with clonal prevalence

fig, axs = plt.subplots(1,3,figsize=(15,5))

# Right
meta = pd.read_csv(os.path.join(path_data, 'cells_meta.csv'), index_col=0)
df_ = (
    meta.assign(clonal=lambda x: x['sample'] + '_' + x['GBC'])
    .groupby('clonal')
    .size()
    .reset_index(name='n')
    .assign(GBC=lambda x: x['clonal'].apply(lambda x: x.split('_')[-1]))
    .assign(sample=lambda x: x['clonal'].apply(lambda x: '_'.join(x.split('_')[:-1])))
    .assign(prevalence=lambda x: x['n'] / x.groupby('sample')['n'].transform('sum'))
    .loc[:, ['GBC', 'prevalence']]
    .set_index('GBC')
)

df_ = (
    clones
    .assign(job=lambda x: x['filtering'] + '_' + x['dimred'] + '_' + x['model'] + '_' + x['tuning'])
    .assign(GBC=lambda x: clones['comparison'].apply(lambda x: x.split('_')[0]))
    .groupby(['GBC', 'job'])
    .agg('mean')
    .reset_index()
    .loc[:, ['GBC', 'f1', 'precision', 'recall']]
    .set_index('GBC')
    .join(df_)
    .sort_values('f1', ascending=False)
)

# Fq
x = df_['prevalence']
y = df_['f1']
corr, p_value = pearsonr(x, y)

axs[0].plot(x, y, 'ko', markersize=2)
sns.regplot(x=x, y=y, ax=axs[0], scatter=False)
format_ax(
    axs[0],
    title=f"Pearson's rho {corr:.2f} (p={p_value:.2e})", 
    xlabel='Clonal prevalence', 
    ylabel='f1'
)
axs[0].spines[['right', 'top']].set_visible(False)

x = df_['prevalence']
y = df_['precision']
corr, p_value = pearsonr(x, y)

axs[1].plot(x, y, 'ko', markersize=2)
sns.regplot(x=x, y=y, ax=axs[1], scatter=False)
format_ax(
    axs[1],
    title=f"Pearson's rho {corr:.2f} (p={p_value:.2e})", 
    xlabel='Clonal Prevalence', 
    ylabel='Precision'
)
axs[1].spines[['right', 'top']].set_visible(False)

x = df_['prevalence']
y = df_['recall']
corr, p_value = pearsonr(x, y)

axs[2].plot(x, y, 'ko', markersize=2)
sns.regplot(x=x, y=y, ax=axs[2], scatter=False)
format_ax(
    axs[2],
    title=f"Pearson's rho {corr:.2f} (p={p_value:.2e})", 
    xlabel='Clonal Prevalence', 
    ylabel='Recall'
)
axs[2].spines[['right', 'top']].set_visible(False)

# Save
fig.tight_layout()
fig.savefig(os.path.join(path_viz, 'clonal_prevalence_performance_relationship.png'))

############## Overall performance by sample


##


############## Overall performance by sample
fig, axs = plt.subplots(1,3,figsize=(15, 5))

# Left
df_ = (
    clones.loc[:, ['sample','ncells']]
    .drop_duplicates()
    .set_index('sample')
    .loc[['AML_clones', 'MDA_clones', 'MDA_lung', 'MDA_PT']]
    .reset_index()
)

bar(df_, x='sample', y='ncells', c='lightgrey', edgecolor='k', ax=axs[0], s=0.75)
format_ax(axs[0], title='n cells', ylabel='n', xticks=df_['sample'])
axs[0].spines[['top', 'right', 'left']].set_visible(False)

# Center
df_ = (
    clones.loc[:, ['sample','n_clones_analyzed']]
    .drop_duplicates()
    .set_index('sample')
    .loc[['AML_clones', 'MDA_clones', 'MDA_lung', 'MDA_PT']]
    .reset_index()
)

bar(df_, x='sample', y='n_clones_analyzed', c='lightgrey', edgecolor='k', ax=axs[1], s=0.75)
format_ax(axs[1], title='n clones', xticks=df_['sample'])
axs[1].spines[['top', 'right', 'left']].set_visible(False)

# Right
df_ = (
    clones
    .assign(job=lambda x: x['filtering'] + '_' + x['dimred'] + '_' + x['model'] + '_' + x['tuning'])
    .groupby(['sample', 'job'])
    .agg('mean')
    .reset_index()
    .assign(feature_type=lambda x: x['job'].map(lambda x: '_'.join(x.split('_')[:-2])))
    .drop(columns=['job'])
)
samples_order = ['AML_clones', 'MDA_clones', 'MDA_lung', 'MDA_PT']
feat_type_colors = create_palette(df_, 'feature_type', 'Spectral', saturation=.85)
box(df_, 'sample', 'f1', c='#E9E7E7', ax=axs[2], order=samples_order)
strip(df_, 'sample', 'f1', by='feature_type', c=feat_type_colors, s=4, ax=axs[2], order=samples_order)
format_ax(axs[2], title='Classification performance', ylabel='f1 score', xticks_size=10)
add_legend(
    label='Feature type',
    colors=feat_type_colors, 
    ax=axs[2],
    loc='upper left',
    bbox_to_anchor=(1,1),
    ncols=1,
    label_size=10,
    ticks_size=8
)
axs[2].spines[['top', 'right']].set_visible(False)

# Save
fig.tight_layout()
fig.savefig(os.path.join(path_viz, 'performance_by_sample_f1.png'))

##############


##


############## Overall performance by sample 2: relationship with clonal complexity

fig, axs = plt.subplots(1,3,figsize=(15,5))

# Right
meta = pd.read_csv(os.path.join(path_data, 'cells_meta.csv'), index_col=0)
df_ = (
    meta.groupby(['sample', 'GBC'])
    .size()
    .reset_index(name='n')
    .assign(prevalence=lambda x: x['n'] / x.groupby('sample')['n'].transform('sum'))
)

def SH(df, sample):
    freqs = df.query(f'sample == "{sample}"')['prevalence'].values
    return -np.sum(freqs * np.log10(freqs))

df_ = pd.Series({ x : SH(df_, x) for x in df_['sample'].unique() }).to_frame('sh')
df_ = (
    clones
    .assign(job=lambda x: x['filtering'] + '_' + x['dimred'] + '_' + x['model'] + '_' + x['tuning'])
    .groupby(['sample', 'job'])
    .agg('mean')
    .reset_index()
    .loc[:, ['sample', 'f1', 'precision', 'recall']]
    .join(df_, on='sample')
    .sort_values('f1')
)

# Fq
x = df_['sh']
y = df_['f1']
corr, p_value = pearsonr(x, y)

axs[0].plot(x, y, 'ko', markersize=4)
sns.regplot(x=x, y=y, ax=axs[0], scatter=False)
format_ax(
    axs[0],
    title=f"Pearson's rho {corr:.2f} (p={p_value:.2e})", 
    xlabel='Shannon Entropy', 
    ylabel='f1'
)
axs[0].spines[['top', 'right']].set_visible(False)

x = df_['sh']
y = df_['precision']
corr, p_value = pearsonr(x, y)

axs[1].plot(x, y, 'ko', markersize=4)
sns.regplot(x=x, y=y, ax=axs[1], scatter=False)
format_ax(
    axs[1],
    title=f"Pearson's rho {corr:.2f} (p={p_value:.2e})", 
    xlabel='Shannon Entropy', 
    ylabel='Precision'
)
axs[1].spines[['top', 'right']].set_visible(False)

x = df_['sh']
y = df_['recall']
corr, p_value = pearsonr(x, y)

axs[2].plot(x, y, 'ko', markersize=4)
sns.regplot(x=x, y=y, ax=axs[2], scatter=False)
format_ax(
    axs[2],
    title=f"Pearson's rho {corr:.2f} (p={p_value:.2e})", 
    xlabel='Shannon Entropy', 
    ylabel='Recall'
)
axs[2].spines[['top', 'right']].set_visible(False)

# Save
fig.tight_layout()
fig.savefig(os.path.join(path_viz, 'sample_SH_performance_relationship.png'))

##############


##


############## Top 3 models

# Data
df_ = (
    
    clones
    .groupby(['sample', 'comparison'])
    .apply(lambda x: 
        x
        .iloc[:,:5]
        .sort_values('f1', ascending=False)
        .head(3)
    )
    .droplevel(2)
    .sort_values(by='f1', ascending=False)
    .reset_index()
    
)

df_['GBC'] = df_['comparison'].str.split('_').map(lambda x: f'{x[0][:6]}...')

# Vix
samples_order = ['AML_clones', 'MDA_clones', 'MDA_lung', 'MDA_PT']
colors = { s:c for s,c in zip(samples_order, sc.pl.palettes.vega_10) }

# Three good
fig, axs = plt.subplots(1,3,figsize=(11,3.5), sharex=True)
for (ax, sample) in zip(axs, samples_order[:-1]):
    n_clones = df_.query('sample == @sample')['GBC'].unique().size
    strip(
        df_
        .query('sample == @sample'), 
        'GBC',
        'f1',
        by='feature_type',
        c=colors[sample], 
        s=4, 
        ax=ax
    )
    ax.set(ylim=(-0.1,1.1))
    ax.hlines(0.7, 0, n_clones, 'k', 'dashed')
    
    medians = (
        df_.query('sample == @sample')
        .groupby('GBC')['f1']
        .agg('median')
        .sort_values(ascending=False)
    )
    for i, y in enumerate(medians):
        ax.hlines(y, i-.25, i+.25, '#4f4e4a', zorder=2)
    
    format_ax(
        ax, 
        title=sample, 
        xlabel='Ranked clones',
        ylabel='f1',
        rotx=0,
        xticks=[],
        xticks_size=6
    )
    ax.spines[['top', 'right']].set_visible(False)
    
    median_f1 = np.median(df_.query('sample == @sample')['f1'])
    std_f1 = np.std(df_.query('sample == @sample')['f1'])
    median_precision = np.median(df_.query('sample == @sample')['precision'])
    std_precision = np.std(df_.query('sample == @sample')['precision'])
    median_recall = np.median(df_.query('sample == @sample')['recall'])
    std_recall = np.std(df_.query('sample == @sample')['recall'])
    
    ax.text(0.05, 0.71, '0.7', transform=ax.transAxes)
    ax.text(0.05, 0.19, 
        f'precision: {median_precision:.2f} +-{std_precision:.2f}', 
        transform=ax.transAxes
    )
    ax.text(0.05, 0.12, 
        f'recall: {median_recall:.2f} +-{std_recall:.2f}', 
        transform=ax.transAxes
    )
    ax.text(0.05, 0.05, 
        f'f1: {median_f1:.2f} +-{std_f1:.2f}', 
        transform=ax.transAxes
    )

# Save
fig.suptitle('Top 3 models')
fig.tight_layout()
fig.savefig(os.path.join(path_viz, 'top3_models_performances_good_samples.png'))

##
    
# MDA_PT    
sample = 'MDA_PT'
n_clones = df_.query('sample == @sample')['GBC'].unique().size

fig, ax = plt.subplots(figsize=(11,3))
strip(
    df_
    .query('sample == @sample'), 
    'GBC',
    'f1',
    by='feature_type',
    c=colors[sample], 
    s=4, 
    ax=ax
)
ax.set(ylim=(-0.1,1.1))
ax.hlines(0.7, 0, n_clones, 'k', 'dashed')
format_ax(
    ax, 
    title=f'{sample}',
    xlabel='Ranked clones', 
    ylabel='f1',
    xticks=[],
    rotx=90
)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

medians = (
    df_.query('sample == @sample')
    .groupby('GBC')['f1']
    .agg('median')
    .sort_values(ascending=False)
)
for i, y in enumerate(medians):
    ax.hlines(y, i-.25, i+.25, '#4f4e4a', zorder=2)
    
median_f1 = np.median(df_.query('sample == @sample')['f1'])
std_f1 = np.std(df_.query('sample == @sample')['f1'])
median_precision = np.median(df_.query('sample == @sample')['precision'])
std_precision = np.std(df_.query('sample == @sample')['precision'])
median_recall = np.median(df_.query('sample == @sample')['recall'])
std_recall = np.std(df_.query('sample == @sample')['recall'])

ax.text(0.05, 0.71, '0.7', transform=ax.transAxes)
ax.text(0.75, 0.94, 
    f'precision: {median_precision:.2f} +-{std_precision:.2f}', 
    transform=ax.transAxes
)
ax.text(0.75, 0.87, 
    f'recall: {median_recall:.2f} +-{std_recall:.2f}', 
    transform=ax.transAxes
)
ax.text(0.75, 0.80, 
    f'f1: {median_f1:.2f} +-{std_f1:.2f}', 
    transform=ax.transAxes
)

# Save
fig.tight_layout()
fig.savefig(os.path.join(path_viz, 'top3_models_performances_MDA_PT.png'))

##############


##


############## 