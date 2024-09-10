from base import *
from plotting import sc_violin, sc_scatter_with_line, sc_pca_variance_ratio, sc_umap, sc_hvg


def check_min_max_adata(adata):
    print((np.max(adata.X), np.min(adata.X)))
    return None


def filtering(adata, plotdir='', min_genes=200, mt_pct=20):

    # Filtering
    sc.pp.filter_cells(adata, min_genes=min_genes)
    max_counts = int(np.mean(adata.obs['total_counts']) + 2 * np.std(adata.obs['total_counts']))
    sc.pp.filter_cells(adata, max_counts=max_counts)
    adata = adata[adata.obs['pct_counts_mt'] < mt_pct, :]

    # Plotting
    sc_violin(adata)
    sc_scatter_with_line(adata.obs['log1p_total_counts'], adata.obs['log1p_n_genes_by_counts'],
                            adata.obs['pct_counts_mt'], x_label='log(total RNA counts)',
                            y_label='log(number of genes)', log=True,
                            plotdir=f'{plotdir}scatter_genes_by_counts.png', max_x=max_counts, min_y=min_genes)
    sc_scatter_with_line(adata.obs['n_genes_by_counts'], adata.obs['pct_counts_mt'], adata.obs['log1p_total_counts'],
                        x_label='Number of genes', y_label='% mitochondrial RNA counts',
                        plotdir=f'{plotdir}scatter_pct_mt_by_genes.png', min_x=min_genes, min_y=mt_pct)
    return adata


def normalize(adata, target_sum=1e6):
    print('before log1p norm')
    print('max: ', np.max(adata.X))
    print('min: ', np.min(adata.X))
    print('mean: ', np.mean(adata.X))
    print("=" * 100)

    # Normalization and Log

    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.raw = adata

    print('after log1p norm')
    print('max: ', np.max(adata.X))
    print('min: ', np.min(adata.X))
    print('mean: ', np.mean(adata.X))
    return adata


def finding_hvg(adata, hvg_num=2000, flavor='seurat', save='.png'):
    # n_top_genes(hvg_num): specific number or None(default)
    sc.pp.highly_variable_genes(adata, n_top_genes=hvg_num, flavor=flavor)
    print(f"detected {np.sum(adata.var['highly_variable'])} highly variable genes")
    print("=" * 100)
    sc_hvg(adata, save=save)
    return adata


def scaling(adata, regress_out=False):
    if regress_out:
        sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    # remove the influence of sequencing depth per cell & mitochondrial gene percent
    sc.pp.scale(adata, max_value=10)
    print('max: ', np.max(adata.X))
    print('min: ', np.min(adata.X))
    print('mean: ', np.mean(adata.X))
    print("=" * 100)
    return adata


def batch_correction(adata, batch_key, n_pcs, save=''):
    sc.tl.pca(adata)
    sc_pca_variance_ratio(adata, log=True, save=f'{save}.png')

    sc.pp.neighbors(adata, n_pcs=n_pcs)
    sc.tl.umap(adata)
    sc_umap(adata, color=batch_key, save=f'{save}_uncorrected.png')
    adata.obsm['X_umap_uncorrected'] = adata.obsm['X_umap'].copy()
    # end uncorrected

    sc.external.pp.harmony_integrate(adata, key=batch_key, max_iter_harmony=50)
    sc.pp.neighbors(adata, n_pcs=n_pcs, use_rep='X_pca_harmony')
    sc.tl.umap(adata)
    sc_umap(adata, color=batch_key, save=f'{save}_harmony.png')
    adata.obsm['X_umap_harmony'] = adata.obsm['X_umap'].copy()
    # end harmony integration
    return adata