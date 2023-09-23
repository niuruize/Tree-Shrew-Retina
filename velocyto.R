

################################################################################
# velocyto
################################################################################

# STEP 1: download TS_genomeannotation.gtf and TS_rmsk.gff
# STEP 2: MAKE LOOM FILE
cd /home/fingerstyle/sampleID/outs
nohup samtools sort -@ 10  -t CB -O BAM -o cellsorted_possorted_genome_bam.bam possorted_genome_bam.bam &
nohup velocyto run10x -m /home/fingerstyle/TS_3.0_genome/TS_rmsk.gff  /home/fingerstyle/sampleID /home/fingerstyle/TS_3.0_genome/genes/genes.gtf & 

# STEP 3: anglysis in python
import os
os.chdir("/home/xll-zxcv/TS/loom/Seurat_dataset/")
os.getcwd()

import scanpy as sc
import anndata 
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

type="Allcell"

#(1) load sparse matrix:
X = io.mmread(f"/home/xll-zxcv/TS/loom/Seurat_dataset/Allcell/counts.mtx")

# create anndata object
adata = anndata.AnnData(
  X=X.transpose().tocsr()
)

#(2) load cell metadata:
cell_meta = pd.read_csv(f"/home/xll-zxcv/TS/loom/Seurat_dataset/Allcell/metadata.csv")

# load gene names:
with open(f"/home/xll-zxcv/TS/loom/Seurat_dataset/Allcell/gene_names.csv", 'r') as f:
  gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']

adata.var.index = gene_names

#(3) load dimensional reduction:
pca = pd.read_csv(f"/home/xll-zxcv/TS/loom/Seurat_dataset/Allcell/pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by sampleID to test:
sc.pl.umap(adata, color='seurat_clusters', frameon=False, save=True)

# (4) load loom data
import scvelo as scv
import scanpy as sc
#import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

import matplotlib.pyplot as plt

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, fontsize=6, color_map='viridis',frameon=False)

#cr.settings.verbosity = 2
#adata = sc.read_h5ad('my_data.h5ad')
#adata=sc.read_loom(f'out/{type}.my_data.loom')

#adata.obsm['X_pca'] = pca.to_numpy()
#adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

#adata
# load loom files for spliced/unspliced matrices for each sample:
sample_loom_1="/home/xll-zxcv/TS/loom/Child1_count.loom"
ldata1 = scv.read( sample_loom_1, cache=True)
#ldata1.obs['batch'] = '1'
#ldata1.obs['BATCH'] = '1'
ldata1

sample_loom_2="/home/xll-zxcv/TS/loom/Child2_count.loom"
ldata2 = scv.read( sample_loom_2, cache=True)
#ldata2.obs['batch'] = '2'
#ldata2.obs['BATCH'] = '2'
ldata2

sample_loom_3="/home/xll-zxcv/TS/loom/Child3_count.loom"
ldata3 = scv.read( sample_loom_3, cache=True)
#ldata3.obs['batch'] = '3'
#ldata3.obs['BATCH'] = '3'
ldata3

sample_loom_4="/home/xll-zxcv/TS/loom/Child4_count.loom"
ldata4 = scv.read( sample_loom_4, cache=True)
#ldata4.obs['batch'] = '4'
#ldata4.obs['BATCH'] = '4'
ldata4

sample_loom_5="/home/xll-zxcv/TS/loom/Child5_count.loom"
ldata5 = scv.read( sample_loom_5, cache=True)
#ldata5.obs['batch'] = '5'
#ldata5.obs['BATCH'] = '5'
ldata5

sample_loom_6="/home/xll-zxcv/TS/loom/Adult1_count.loom"
ldata6 = scv.read( sample_loom_6, cache=True)
#ldata6.obs['batch'] = '6'
#ldata6.obs['BATCH'] = '6'
ldata6

sample_loom_7="/home/xll-zxcv/TS/loom/Adult2_count.loom"
ldata7 = scv.read( sample_loom_7, cache=True)
#ldata7.obs['batch'] = '7'
#ldata7.obs['BATCH'] = '7'
ldata7

sample_loom_8="/home/xll-zxcv/TS/loom/Adult3_count.loom"
ldata8 = scv.read( sample_loom_8, cache=True)
#ldata8.obs['batch'] = '8'
#ldata8.obs['BATCH'] = '8'
ldata8

sample_loom_9="/home/xll-zxcv/TS/loom/Adult4_count.loom"
ldata9 = scv.read( sample_loom_9, cache=True)
#ldata9.obs['batch'] = '9'
#ldata9.obs['BATCH'] = '9'
ldata9

sample_loom_10="/home/xll-zxcv/TS/loom/Adult5_count.loom"
ldata10 = scv.read( sample_loom_10, cache=True)
#ldata10.obs['batch'] = '10'
#ldata10.obs['BATCH'] = '10'
ldata10

sample_loom_11="/home/xll-zxcv/TS/loom/Old1_count.loom"
ldata11 = scv.read( sample_loom_11, cache=True)
#ldata11.obs['batch'] = '11'
#ldata11.obs['BATCH'] = '11'
ldata11

sample_loom_12="/home/xll-zxcv/TS/loom/Old2_count.loom"
ldata12 = scv.read( sample_loom_12, cache=True)
#ldata12.obs['batch'] = '12'
#ldata12.obs['BATCH'] = '12'
ldata12

sample_loom_13="/home/xll-zxcv/TS/loom/Old3_count.loom"
ldata13 = scv.read( sample_loom_13, cache=True)
#ldata13.obs['batch'] = '13'
#ldata13.obs['BATCH'] = '13'
ldata13

sample_loom_14="/home/xll-zxcv/TS/loom/Old4_count.loom"
ldata14 = scv.read( sample_loom_14, cache=True)
#ldata14.obs['batch'] = '14'
#ldata14.obs['BATCH'] = '14'
ldata14

sample_loom_15="/home/xll-zxcv/TS/loom/Old5_count.loom"
ldata15 = scv.read( sample_loom_15, cache=True)
#ldata15.obs['batch'] = '15'
#ldata15.obs['BATCH'] = '15'
ldata15

# (5) merge loom data
ldata1.obs.index[0:2]
ldata2.obs.index[0:2]
ldata3.obs.index[0:2]
ldata4.obs.index[0:2]
ldata5.obs.index[0:2]
ldata6.obs.index[0:2]
ldata7.obs.index[0:2]
ldata8.obs.index[0:2]
ldata9.obs.index[0:2]
ldata10.obs.index[0:2]
ldata11.obs.index[0:2]
ldata12.obs.index[0:2]
ldata13.obs.index[0:2]
ldata14.obs.index[0:2]
ldata15.obs.index[0:2]

# add cell id

# ldata1
barcodes = [bc.split(':')[1] for bc in ldata1.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata1.obs.index = barcodes
ldata1.obs.index[0:5]

# ldata2
barcodes = [bc.split(':')[1] for bc in ldata2.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata2.obs.index = barcodes
ldata2.obs.index[0:5]

# ldata3
barcodes = [bc.split(':')[1] for bc in ldata3.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata3.obs.index = barcodes
ldata3.obs.index[0:5]

# ldata4
barcodes = [bc.split(':')[1] for bc in ldata4.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata4.obs.index = barcodes
ldata4.obs.index[0:5]

# ldata5
barcodes = [bc.split(':')[1] for bc in ldata5.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata5.obs.index = barcodes
ldata5.obs.index[0:5]

# ldata6
barcodes = [bc.split(':')[1] for bc in ldata6.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata6.obs.index = barcodes
ldata6.obs.index[0:5]

# ldata7
barcodes = [bc.split(':')[1] for bc in ldata7.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata7.obs.index = barcodes
ldata7.obs.index[0:5]

# ldata8
barcodes = [bc.split(':')[1] for bc in ldata8.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata8.obs.index = barcodes
ldata8.obs.index[0:5]

# ldata9
barcodes = [bc.split(':')[1] for bc in ldata9.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata9.obs.index = barcodes
ldata9.obs.index[0:5]

# ldata10
barcodes = [bc.split(':')[1] for bc in ldata10.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata10.obs.index = barcodes
ldata10.obs.index[0:5]

# ldata11
barcodes = [bc.split(':')[1] for bc in ldata11.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata11.obs.index = barcodes
ldata11.obs.index[0:5]

# ldata12
barcodes = [bc.split(':')[1] for bc in ldata12.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata12.obs.index = barcodes
ldata12.obs.index[0:5]

# ldata13
barcodes = [bc.split(':')[1] for bc in ldata13.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata13.obs.index = barcodes
ldata13.obs.index[0:5]

# ldata14
barcodes = [bc.split(':')[1] for bc in ldata14.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata14.obs.index = barcodes
ldata14.obs.index[0:5]

# ldata15
barcodes = [bc.split(':')[1] for bc in ldata15.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata15.obs.index = barcodes
ldata15.obs.index[0:5]

ldata1.var.head()
ldata2.var.head()
ldata3.var.head()
ldata4.var.head()
ldata5.var.head()
ldata6.var.head()
ldata7.var.head()
ldata8.var.head()
ldata9.var.head()
ldata10.var.head()
ldata11.var.head()
ldata12.var.head()
ldata13.var.head()
ldata14.var.head()
ldata15.var.head()

ldata1.var_names_make_unique()
ldata2.var_names_make_unique()
ldata3.var_names_make_unique()
ldata4.var_names_make_unique()
ldata5.var_names_make_unique()
ldata6.var_names_make_unique()
ldata7.var_names_make_unique()
ldata8.var_names_make_unique()
ldata9.var_names_make_unique()
ldata10.var_names_make_unique()
ldata11.var_names_make_unique()
ldata12.var_names_make_unique()
ldata13.var_names_make_unique()
ldata14.var_names_make_unique()
ldata15.var_names_make_unique()

# merge
ldata = ldata1.concatenate([ldata2,ldata3,ldata4,ldata5,ldata6,ldata7,ldata8,ldata9,ldata10,ldata11,ldata12,ldata13,ldata14,ldata15])
ldata = sc.AnnData.concatenate(ldata1,ldata2,ldata3,ldata4,ldata5,ldata6,ldata7,ldata8,ldata9,ldata10,ldata11,ldata12,ldata13,ldata14,ldata15,batch_key = 'BATCH')
ldata.obs

# (6) merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)
adata.obs

#save data
adata.write_h5ad(f'/home/xll-zxcv/TS/loom/Seurat_dataset/Allcell/Allcell.adata_ldata.h5ad')

# (7) scVelo
adata = sc.read(f'/home/xll-zxcv/TS/loom/Seurat_dataset/Allcell/Allcell.adata_ldata.h5ad')
#scv.pp.filter_and_normalize(adata, min_shared_counts=5, min_shared_cells=3, log=True)
scv.pp.filter_and_normalize(adata)

scv.pp.moments(adata, n_neighbors=30, n_pcs=30)

import gc
gc.collect()
#
temp_pre= f"Allcell_nue.in_process2" 
if False==os.path.exists(f"{temp_pre}.velo.gz.h5ad"):
    scv.tl.recover_dynamics(adata, var_names='all', n_jobs=18)
    scv.tl.velocity(adata, mode='dynamical')
    adata.write(f"{temp_pre}.velo.gz.h5ad", compression='gzip')
    print(">>Write to file")
else:
    adata = sc.read(f"{temp_pre}.velo.gz.h5ad", compression='gzip', ext="h5ad")
    print(">>read from file")

scv.tl.velocity_graph(adata, n_jobs=18)


scv.pl.velocity_embedding_stream(adata, basis="umap")
scv.pl.velocity_embedding(adata, basis="umap", save='embedding_Allcell.pdf')

scv.pl.velocity_embedding_grid(adata, basis='umap', color='celltype', figsize=
(12,12), arrow_size=1.8,  arrow_length=2.5, size=100, alpha=0.1,
save='embedding_grid_Allcell.pdf', title='', scale=0.1)

scv.pl.velocity_embedding_stream(adata, basis='umap', color=
['seurat_clusters'], figsize=(8,8), palette =
("#8cb883","#643d2c","#FF8066",'#c94052',"#caa795",'#80d1c8',"#F7ab08","#2C73D2","#845EC2"),
arrow_size=2, linewidth=1.5, legend_fontsize=25, dpi=900,
save='embedding_stream_seurat_clusters_Allcell.svg', title='')

scv.pl.velocity_embedding_stream(adata, basis='umap', color=
['celltype'], figsize=(8,8), palette =
("#8cb883","#643d2c","#FF8066",'#c94052',"#caa795",'#80d1c8',"#F7ab08","#2C73D2","#845EC2"),
arrow_size=2, linewidth=1.2, legend_fontsize=25, dpi=900,
save='embedding_stream_celltype_Allcell.svg', title='')

scv.pl.velocity_embedding(adata, arrow_length=2, arrow_size=1.5, dpi=900,
figsize=(8,8), save='embedding_stream2_Allcell.pdf')

#(8) visualization
scv.pl.velocity(adata,
["OPN1LW","OPN1SW","ARR3","GRK7"],
figsize=(5,5), dpi=600, save='Interprete velocity_Allcell.png', ncols=2)

scv.pl.velocity(adata,
["PDE6H","PDE6C","GNAT2","RCVRN"],
figsize=(5,5), dpi=600, save='Interprete velocity_Allcell2.png', ncols=2)

scv.pl.velocity(adata,
["CRX","PDE6B","PDE6A","ROM1"], figsize=(5,5), dpi=600,
save='Interprete velocity_Allcell3.png', ncols=2)

scv.pl.velocity(adata,
["SAG","CA10","GRIK1","GRIK1li1"],
figsize=(5,5), dpi=600, save='Interprete velocity_Allcell4.png', ncols=2)

scv.pl.velocity(adata,
["ISL1","TRPM1","NETO1","ONECUT1"],
figsize=(5,5), dpi=600, save='Interprete velocity_Allcell5.png', ncols=2)

scv.pl.velocity(adata,
["ONECUT2","NDRG1","RLBP1","WIF1"], figsize=
(5,5), dpi=600, save='Interprete velocity_Allcell6.png', ncols=2)

scv.pl.velocity(adata,
["CHRDL1","TF","GLUL","C1QA"], figsize=(5,5), dpi=600,
save='Interprete velocity_Allcell7.png', ncols=2)

scv.pl.velocity(adata,
["C1QB","CSF1R","CX3CR1","PTPRC"],
figsize=(5,5), dpi=600, save='Interprete velocity_Allcell8.png', ncols=2)

scv.pl.velocity(adata,
["AQP4","SLC1A3","SLC14A1","NRN1"],
figsize=(5,5), dpi=600, save='Interprete velocity_Allcell9.png', ncols=2)

scv.pl.velocity(adata,
["SLC17A6","RBPMS","THY1","PAX6"],
figsize=(5,5), dpi=600, save='Interprete velocity_Allcell10.png', ncols=2)

scv.pl.velocity(adata,
["SYNPR","NRXN1","MEIS2","TFAP2A"],
figsize=(5,5), dpi=600, save='Interprete velocity_Allcell11.png', ncols=2)

scv.pl.velocity(adata,
["TFAP2B","C1QL2","SLC6A1","TCF4"],
figsize=(5,5), dpi=600, save='Interprete velocity_Allcell12.png', ncols=2)

scv.pl.velocity(adata,
["SLC6A9","SLC17A8","NFIA","NFIB"],
figsize=(5,5), dpi=600, save='Interprete velocity_Allcell13.png', ncols=2)

scv.pl.velocity(adata,
["SLC18A3","MEGF11","RPE65","CHAT"],
figsize=(5,5), dpi=600, save='Interprete velocity_Allcell14.png', ncols=2)

scv.pl.velocity(adata,
["GRM8","PAX2","SPARCL1","INPP5D"],
figsize=(5,5), dpi=600, save='Interprete velocity_Allcell15.png', ncols=2)

scv.pl.scatter(adata, 'ARR3', color=
['celltype', 'velocity'], save='Interprete velocity2_Allcell.pdf')

#(8)velocity and consistency
scv.tl.velocity_confidence(adata)
keys = 'velocity_length','velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], dpi=900, figsize=(6,5), save='Speed and coherence_Allcell.pdf')

df = adata.obs.groupby('celltype')[keys].mean().T
df.style.background_gradient(cmap='coolwarm', axis=1)

scv.pl.velocity_graph(adata, threshold=.1, save='velocity_graph_Allcell.pdf')

#(9) Graph of velocity and pseudotime
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', figsize=(5,5), dpi=300 ,save='velocity_pseudotime_Allcell.pdf')
#
adata_subset = adata[:40000]
scv.tl.velocity_pseudotime(adata_subset)
scv.pl.scatter(adata_subset, color='velocity_pseudotime', cmap='gnuplot', figsize=(5,5), dpi=300 ,save='velocity_pseudotime_Allcell.pdf')

###PAGA
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='celltype')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')

scv.pl.paga(adata, basis='umap', size=80, alpha=.15,
            min_edge_width=2, node_size_scale=1.5, figsize=(6,5),arrowsize=20, node_size_power=0.5,dpi=900, save='paga_Allcell.svg')
#PAGA2
adata_subset = adata[:40000]
adata_subset.uns['neighbors']['distances'] = adata_subset.obsp['distances']
adata_subset.uns['neighbors']['connectivities'] = adata_subset.obsp['connectivities']

scv.tl.paga(adata_subset, groups='celltype')
df = scv.get_df(adata_subset, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')

scv.pl.paga(adata_subset, basis='umap', size=40, alpha=.1,
            min_edge_width=1.5, node_size_scale=1.2, figsize=(6,5),arrowsize=20, node_size_power=0.2,dpi=300, save='paga_Allcell.svg')
#
##Latent time
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, dpi=900, save='latent_time_Allcell.pdf')
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]

scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='celltype', n_convolve=100, save='heatmap_Allcell.pdf')






