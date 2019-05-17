import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
import os
import sys

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
res_dir = "/Users/kriemo/Projects/sc_repos/vibhakar/results/2019-05-03_qc"

sc.settings.set_figure_params(dpi=300)

samples = [
  "1_831",
  "1_925_010819",
  "2_1008",
  "2_934",
  "2_945_010819",
  "3_1130_010819",
  "4_1167_010819",
  "4_1355",
  "Foreman_1",
  "Foreman_1125",
  "Foreman_1128",
  "Foreman_2",
  "Foreman_3",
  "Foreman_4"
]

dat_dir = "/Users/kriemo/Projects/sc_repos/vibhakar/data/cellranger/results/"

adatas = [sc.read_10x_mtx(os.path.join(dat_dir, sample, "outs", "filtered_feature_bc_matrix"),
                        var_names='gene_symbols', make_unique = True) for sample in samples]

for idx, sample in enumerate(samples):
  adatas[idx].obs["sample"] = sample

adata = AnnData.concatenate(*adatas)

adata.var = adata.var.iloc[:, 0:2]
sc.pl.highest_expr_genes(adata, n_top=20)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=5)

mito_genes = adata.var_names.str.startswith('MT-')
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
adata.obs['percent_mito'] = np.sum(
    adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1


sc.pl.violin(adata, ['n_genes', 'percent_mito'],
               groupby = "sample",
               jitter=0.4, multi_panel=False,
               rotation = 90.0,stripplot = False,
               save = "_qc_metrics.pdf")


sc.pl.scatter(adata, x='n_counts', y='percent_mito', color= 'sample')
sc.pl.scatter(adata, x='n_counts', y='n_genes', color = 'sample')


adata = adata[adata.obs['n_genes'] < 7000, :]
adata = adata[adata.obs['percent_mito'] < 0.35, :]

sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)

adata.raw = adata

sc.pp.highly_variable_genes(adata, min_mean=0.005, max_mean=3, min_disp=0.35)
sc.pl.highly_variable_genes(adata)


adata = adata[:, adata.var['highly_variable']]
sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

sc.pl.pca_variance_ratio(adata, log=True)

#adata.write(os.path.join(res_dir, "sc.h5ad"))


#tl.paga(adata)
#pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
#tl.umap(adata, init_pos='paga')

sc.pp.neighbors(adata, n_neighbors = 25, n_pcs = 30)
sc.tl.umap(adata)

sc.pl.umap(adata, color=['sample'], save = "_by_sample.png")

genes = [
  'CD3D',
  'PSAP',
  'TOP2A',
  'MYCN',
  'OTX2',
  'CDK6',
  'ACVR1',
  'GFI1B',
  'TERT',
  'SNCAIP',
  'GLI2',
  'YAP1' ]

sc.pl.umap(adata, color=genes, ncols = 4, save = "_various_genes.png")


import bbknn
bbknn.bbknn(adata, batch_key = "batch")
sc.tl.umap(adata)
sc.pl.umap(adata, color=['sample'], save = "_by_sample_bbknn.png")


sc.pl.umap(adata, color=genes, ncols = 4, save = "_various_genes_bbknn.png")

sc.tl.leiden(adata)

sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')


sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

pd.DataFrame(adata.uns['names']).head(5)

]:

result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
a = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names']}).head(5)

marker_genes = list(a[0:4].values.flatten('F'))

ax = sc.pl.stacked_violin(adata, marker_genes, groupby='leiden', rotation=90)


ax = sc.pl.dotplot(adata, marker_genes, groupby='leiden')

sc.pl.umap(adata, color=['CD74', 'FTH1', 'LINC00461'])



adata.write(os.path.join(dat_dir, "sc.h5ad"), compression = 'gzip')


