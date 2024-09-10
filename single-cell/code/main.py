###
# Written by Takwon Yoo
# email : djinissrh1@snu.ac.kr
###


from base import *
from plotting import sc_umap, sc_matrixplot, plot_cell_frac_dpi, sc_paga, sc_draw_graph
from pipeline import filtering, normalize, check_min_max_adata, finding_hvg, scaling, batch_correction


CELL_TYPE_COLOR = ['#9467BD', '#1D6FA8', '#FF7F0E', '#2CA02C', '#D62728']
# Quiescent MuSC, Activated MuSC, Myoblast/Myocyte, Myonuclei(IIb), Myonuclei(IIa)
BATCH_COLORS = ['#ef476f', '#06d6a0', '#CCCCCC']  # Geri, Mid, public
DPI_COLORS = ['#0005FB', '#4E4FF9', '#AAACFF', '#CCCCCC', '#FFA8B0', '#FA554F', '#F70000']


def concat_BD_data():
    """
    1. concat middle-aged and geriatric single-cell data (BD rhapsody)
    """
    sample_list = ['SampleTag06_mm', 'SampleTag07_mm']
    adata = load_BD_csvs(sample_list)
    # AnnData object with n_obs × n_vars = 1410 × 25235

    adata.obs['is.FACS'] = True
    adata.obs['tissue'] = 'Tibialis Anterior'
    adata.obs['batch'] = None
    adata.obs.loc[adata.obs['dataset'] == 'SampleTag06_mm', 'batch'] = 'Mid'
    adata.obs.loc[adata.obs['dataset'] == 'SampleTag07_mm', 'batch'] = 'Geri'
    adata.obs['dataset'] = None
    adata.obs['dataset'] = 'YFP+_MuSCs'
    adata.obs['cell_type'] = 'Our_data'
    adata.obs['sample'] = None
    adata.obs['sample'] = 'Our_data'
    adata.obs['DPI'] = None
    adata.obs['DPI'] = '0'
    adata.write('./adata_combined.h5ad')
    return None


def preprocess_BD_data():
    """
    2. data preprocessing (filter & normalize, BD)
    """
    os.chdir('./')
    adata = sc.read_h5ad('./adata_combined.h5ad')
    adata.var['mt'] = adata.var_names.str.startswith(f'mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True)
    adata = filtering(adata, plotdir=PLOT_DIR)  # 60 cells filtered
    adata = normalize(adata, target_sum=10 ** 4)
    adata.obs['percent.mt'] = adata.obs['pct_counts_mt']
    check_min_max_adata(adata)  # (7.06, 0.0)
    adata.write('./adata_combined_norm.h5ad')
    return None


def load_public_data():
    """
    1. load public data & select GSE162172, GSE138826 and MuSC-lineage mononuncleated cells
    """
    os.chdir('./')
    adata = sc.read_h5ad(f'./data/scMuscle_mm10_slim_v1-1.h5ad')

    check_min_max_adata(adata)    # (8.9, 0): already normalized
    adata.obsm['X_umap'] = adata.obsm['X_umap_harmony']
    sc_umap(adata, save='_harmony_res.1.2_IDs.png', color='harmony_res.1.2_IDs')    # harmony_res.1.2_IDs: clusters in original paper

    # select GSE162172 and GSE138826
    adata_sub = adata[adata.obs['sample'].isin(['Old_D0_A', 'Old_D0_B', 'Old_D0_C', 'Old_D0_D', 'Old_D1_A', 'Old_D1_B',
                                            'Old_D2_A', 'Old_D2_B', 'Old_D2_C', 'Old_D3_5_A', 'Old_D3_5_B', 'Old_D5_A',
                                            'Old_D5_B', 'Old_D5_C', 'Old_D7_A', 'Old_D7_B', 'Old_D7_C', 'Old_D1_C',
                                            'Old_D1_D', 'Old_D3_5_C', 'Old_D3_5_D', 'oprescu_D0_5', 'oprescu_D2',
                                            'oprescu_D3_5', 'oprescu_D5', 'oprescu_D10', 'oprescu_D21', 'oprescu_D0'])]
    sc_umap(adata_sub, save='_harmony_res.1.2_IDs_subset.png', color='harmony_res.1.2_IDs')

    object_cells = ['MuSCs', 'Myoblasts/Progenitors', 'Myonuclei (Type IIb)', 'Myonuclei (Type IIx)']
    adata_muscle = adata_sub[adata_sub.obs['harmony_res.1.2_IDs'].isin(object_cells)]
    sc_umap(adata_muscle, save='_public_MuSC_lineage.png', color='harmony_res.1.2_IDs')

    adata_muscle.obs['dataset'] = None
    adata_muscle.obs.loc[adata_muscle.obs['sample'].str.startswith('Old'), 'dataset'] = 'McKellar'
    adata_muscle.obs.loc[adata_muscle.obs['sample'].str.startswith('oprescu'), 'dataset'] = 'Oprescu'

    adata_muscle.obs['DPI'] = None
    adata_muscle.obs.loc[adata_muscle.obs['sample'].str.contains('D0'), 'DPI'] = '0'
    adata_muscle.obs.loc[adata_muscle.obs['sample'].str.contains('D0_5'), 'DPI'] = '0_5'
    adata_muscle.obs.loc[adata_muscle.obs['sample'].str.contains('D1'), 'DPI'] = '1'
    adata_muscle.obs.loc[adata_muscle.obs['sample'].str.contains('D2'), 'DPI'] = '2'
    adata_muscle.obs.loc[adata_muscle.obs['sample'].str.contains('D3_5'), 'DPI'] = '3_5'
    adata_muscle.obs.loc[adata_muscle.obs['sample'].str.contains('D5'), 'DPI'] = '5'
    adata_muscle.obs.loc[adata_muscle.obs['sample'].str.contains('D7'), 'DPI'] = '7'
    adata_muscle.obs.loc[adata_muscle.obs['sample'].str.contains('D10'), 'DPI'] = '10'
    adata_muscle.obs.loc[adata_muscle.obs['sample'].str.contains('D21'), 'DPI'] = '21'

    adata_muscle.obs['batch'] = None
    adata_muscle.obs['batch'] = 'Public'
    adata_muscle.obs['total_counts'] = adata_muscle.obs['nCount_RNA']
    adata_muscle.obs['n_genes_by_counts'] = adata_muscle.obs['nFeature_RNA']
    adata_muscle.obs['cell_type'] = adata_muscle.obs['harmony_res.1.2_IDs']
    adata_muscle.write('./adata_public_muscle.h5ad')
    return None


def concat_BD_and_10x():
    """
    2. concat our data (Mid & Geri) to public data (GSE138826 & GSE162172)
    """
    adata_public = sc.read_h5ad('./adata_public_muscle.h5ad')
    adata = sc.read_h5ad('./adata_combined_norm.h5ad')

    common_vars = set(adata.var_names) & set(adata_public.var_names)    # 20,499
    adata = adata[:, adata.var_names.isin(common_vars)]
    adata_public = adata_public[:, adata_public.var_names.isin(common_vars)]
    adata_com = ad.concat([adata_public, adata])
    check_min_max_adata(adata_com)
    adata_com.write('./adata_BD_10x_combined.h5ad')
    return None


def BD_10x_processing():
    """
    3. data processing and integration (our data & public data)
    """
    os.chdir('./')
    adata = sc.read_h5ad('./adata_BD_10x_combined.h5ad')
    check_min_max_adata(adata)  # (8.6, 0.0)
    adata = finding_hvg(adata)  # 2000 HVGs
    adata = scaling(adata)
    check_min_max_adata(adata)  # (10.0, -7.2)
    adata = batch_correction(adata, batch_key='sample', n_pcs=50, save='_combined')
    adata.write('./adata_BD_10x_integrated.h5ad')
    return None


def BD_10x_unbiased_clustering():
    """
    4. unbiased clustering + gene expression
    """
    os.chdir('./')
    adata = sc.read_h5ad('./adata_BD_10x_integrated.h5ad')
    adata.obsm['X_umap'] *= -1

    sc.tl.louvain(adata, resolution=0.5, key_added=f'louvain_0.5')
    sc_umap(adata, color=f'louvain_0.5', ncols=1, save=f'_louvain_0.5.png')
    var_names = {'MuSC Quiescence': ['Cd34', 'Spry1', 'Sdc4'], 'MuSC Marker': ['Pax7', 'Itga7', 'Vcam1', 'Myf5'],
                 'Myoblast/Myocyte': ['Myod1', 'Myog', 'Myh3'], 'Myofiber': ['Myh1', 'Myh4', 'Ckm']}
    sc_matrixplot(adata, var_names=var_names, groupby='louvain_0.5', swap_axes=False, save='raw.png',
                  dendrogram=False, vmax=2)
    annot_dict = {'0': 'Myonuclei (Fast fiber)', '3': 'Activated MuSC', '2': 'Myonuclei (Slow fiber)',
                  '1': 'Quiescent MuSC', '4': 'Myoblast/Myocyte', '9': 'del', '6': 'Myonuclei (Fast fiber)',
                  '5': 'del', '7': 'del', '8': 'del', '10': 'del'}
    adata.obs['annot'] = adata.obs[f'louvain_0.5'].map(annot_dict).astype('category')
    adata = adata[adata.obs['annot'] != 'del']
    desired_order = ['Quiescent MuSC', 'Activated MuSC', 'Myoblast/Myocyte',
                     'Myonuclei (Slow fiber)', 'Myonuclei (Fast fiber)']
    adata.obs['annot'] = adata.obs['annot'].astype(pd.CategoricalDtype(categories=desired_order, ordered=True))
    sc_matrixplot(adata, var_names=var_names, groupby='annot', swap_axes=False, save='filtered.png',
                  dendrogram=False, vmax=2)
    sc_umap(adata, color='annot', save='_new_cluster.png', palette=CELL_TYPE_COLOR, size=50)
    adata.write('./adata_BD_10x_clustered.h5ad')
    return None


def visualize_clustered_BD_10x():
    """
    5-1 load pre-clustered adata and generates various plots based on metadata
    """
    os.chdir('./')
    adata = sc.read_h5ad('./adata_BD_10x_clustered.h5ad')

    desired_order = ['Mid', 'Geri', 'Public']
    adata.obs['batch'] = adata.obs['batch'].astype(pd.CategoricalDtype(categories=desired_order, ordered=True))
    adata_public = adata[adata.obs['batch'] == 'Public'].copy()
    adata_wo_mid = adata[adata.obs['batch'] != 'Mid'].copy()
    adata_wo_geri = adata[adata.obs['batch'] != 'Geri'].copy()
    sc_umap(adata_public, color='batch', palette=[BATCH_COLORS[2]], save='_batch_all.png', size=50)
    sc_umap(adata_wo_mid, color='batch', palette=[BATCH_COLORS[0], BATCH_COLORS[2]], save='_batch_wo_mid.png', size=50)
    sc_umap(adata_wo_geri, color='batch', palette=[BATCH_COLORS[1], BATCH_COLORS[2]], save='_batch_wo_geri.png', size=50)

    cell_list = ['Quiescent MuSC', 'Activated MuSC', 'Myoblast/Myocyte']
    adata = adata[adata.obs['annot'].isin(cell_list)]
    column_list = ['Mid', 'Geri', '0', '0_5', '1', '2', '3_5', '5', '7', '10', '21']
    adata.obs['DPI'] = adata.obs['DPI'].astype(pd.CategoricalDtype(categories=column_list, ordered=True))
    adata.obs.loc[adata.obs['batch'] == 'Mid', 'DPI'] = 'Mid'
    adata.obs.loc[adata.obs['batch'] == 'Geri', 'DPI'] = 'Geri'
    cell_type_palette = CELL_TYPE_COLOR[:3]
    figname = 'stacked_bar_dpi'
    plot_cell_frac_dpi(cell_list, adata, column_list, cell_type_palette, figname)
    return None


def BD_10x_trajectory_inference_analysis():
    """
    5-2 trajectory inference
    """
    os.chdir('./')
    adata = sc.read_h5ad('./adata_BD_10x_clustered.h5ad')

    color = 'annot'
    gene_list = ['Pax7', 'Vcam1', 'Itga7', 'Cd63', 'Cd200', 'Cd34', 'Myod1', 'Myog', 'Spry1',
                 'Sdc4', 'Myl9', 'Myh1', 'Myh3', 'Myh4', 'Myh7']

    sc.tl.draw_graph(adata)
    sc.tl.diffmap(adata)
    sc.pp.neighbors(adata, use_rep='X_diffmap')
    sc.tl.draw_graph(adata)
    sc.tl.paga(adata, groups=color)
    sc.pl.paga(adata, color=color, save=f'_harmony_{color}.png', frameon=False, threshold=0)
    sc.tl.draw_graph(adata, init_pos='paga')

    adata.obsm['X_draw_graph_fa'][:, 0] = adata.obsm['X_draw_graph_fa'][:, 0] * -1
    sc_paga(adata, color=color, save=f'_harmony_{color}.png', frameon=False, threshold=0)

    sc_draw_graph(adata, color=color, save=f'_harmony_{color}.png')
    sc_draw_graph(adata, color=gene_list, save='_harmony_markers.png', color_map='magma')

    adata_public = adata[adata.obs['batch'] == 'Public'].copy()
    adata_wo_mid = adata[adata.obs['batch'] != 'Mid'].copy()
    adata_wo_geri = adata[adata.obs['batch'] != 'Geri'].copy()
    sc_draw_graph(adata_public, color='batch', palette=[BATCH_COLORS[2]], save='_batch_all.png', size=50)
    sc_draw_graph(adata_wo_mid, color='batch', palette=[BATCH_COLORS[0], BATCH_COLORS[2]], save='_batch_wo_mid.png', size=50)
    sc_draw_graph(adata_wo_geri, color='batch', palette=[BATCH_COLORS[1], BATCH_COLORS[2]], save='_batch_wo_geri.png', size=50)

    adata_dpi = adata[adata.obs['DPI'].isin(['0', '0_5', '3_5', '5', '7', '10', '21'])]
    sc_draw_graph(adata_dpi, color='DPI', palette=DPI_COLORS, save='_dpi.png')

    sc.pl.paga_compare(
        adata, threshold=0, title='', right_margin=0.2, size=10, edge_width_scale=0.5,
        legend_fontsize=12, fontsize=12, frameon=False, edges=True, save=True)
    adata.write('./adata_BD_10x_trajectory.h5ad')
    return None