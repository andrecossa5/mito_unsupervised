"""
Prep MDA_PT data for iterative_scheme.
"""

# Code
import os
import warnings
from sklearn.metrics import normalized_mutual_info_score
from itertools import chain
warnings.filterwarnings('ignore')

from mito_utils.preprocessing import *
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.clustering import *
from mito_utils.iterative import *
from mito_utils.plotting_base import *
from mito_utils.heatmaps_plots import plot_heatmap
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/mito_bench/'
path_data = os.path.join(path_main, 'PT_subsampled') # Here subsampled
path_output = os.path.join(path_main, 'results', 'unsupervised_clones', 'output')
path_viz = os.path.join(path_main, 'results', 'unsupervised_clones', 'visualization')
path_tmp = os.path.join(path_main, 'results', 'unsupervised_clones', 'downstream_files')

##

# AFM first step, no filtering of low numebr of cells clones + variants filtering
afm = read_one_sample(path_data, 'PT_subsampled', with_GBC=True)

##

# Here we go
AFM_list = [afm]

i = 0
go = True
variants = {}

while go:
    
    if i == 0:
        L = [ 
                one_iteration(
                    x, mode='iterative_splitting', 
                    rank_by='custom_perc_tresholds',
                    min_clone_perc=.9
                ) \
                for x in AFM_list 
        ]
        labels, tests, var, AFM_list = zip(*L)
        variants[f'it_{i}'] = pool_variants(var)
        update_labels(afm, labels, i=i)
    else:
        AFM_list = [ subset_afm(x) for x in AFM_list ]
        AFM_list = list(chain.from_iterable(AFM_list))
        L = [ 
                one_iteration(
                    x, mode='iterative_splitting', 
                    rank_by='custom_perc_tresholds',
                    min_clone_perc=.9
                ) \
                for x in AFM_list if x.shape[0]>10
        ]
        labels, tests, var, AFM_list = zip(*[ x for x in L if x[3] is not None ])
        variants[f'it_{i}'] = pool_variants(var)
        update_labels(afm, labels, i=i)

    if i != 0:
        go = np.sum(tests)>1
    else:
        go = np.sum(tests)>0
    
    i+=1


##


# Test concordance and visualize
iteration = 'it_1'
retain = []
for x in afm.obs[iteration].value_counts().index:
    if bool(re.search('[.]', str(x))):
        x = x.split('.')
        if (x[-1] != 'unassigned') and (x[-1] != 'poor_quality'):
            retain.append('.'.join(x))
    else:
        if x != 'unassigned' and x != 'poor_quality':
            retain.append(x)

df = afm.obs.loc[afm.obs[iteration].isin(retain)]
df[iteration].value_counts()


# All   
fig, axs = plt.subplots(1,2, figsize=(15,8))

ari_all = custom_ARI(df['GBC'].astype('str'), df[iteration].astype('str'))
nmi_all = normalized_mutual_info_score(df['GBC'].astype('str'), df[iteration].astype('str'))
df_ = pd.crosstab(df['GBC'].astype('str'), df[iteration].astype('str'))
plot_heatmap(df_, ax=axs[0], y_names_size=5, rank_diagonal=True, label='n cells')
format_ax(title=f'NMI: {nmi_all:.2f}, ARI: {ari_all:.2f}', ax=axs[0], rotx=90)

# > 10 cells
top_clones = df.groupby('GBC').size().loc[lambda x: x>=10].index
df = df.query('GBC in @top_clones')
ari_top = custom_ARI(df['GBC'].astype('str'), df[iteration].astype('str'))
nmi_top = normalized_mutual_info_score(df['GBC'].astype('str'), df[iteration].astype('str'))

df_ = pd.crosstab(df['GBC'].astype('str'), df[iteration].astype('str'))
plot_heatmap(df_, ax=axs[1], y_names_size=5, rank_diagonal=True, label='n cells')
format_ax(title=f'NMI: {nmi_top:.2f}, ARI: {ari_top:.2f}', ax=axs[1], rotx=90)
fig.tight_layout()
plt.show()





