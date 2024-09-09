import matplotlib.pyplot as plt
from base import *


def sc_umap(adata, color=None, save=None, ncols=4, palette=None, legend=False, color_map=None, size=None,
            frameon=False, edges=False, legend_fontsize=None, title=None, use_raw=None):
    if legend:
        sc.pl.umap(adata, color=color, ncols=ncols, save=save, legend_loc='on data', title=title,
                   legend_fontsize=legend_fontsize, legend_fontoutline=True, palette=palette, color_map=color_map,
                   size=size, frameon=frameon, edges=edges, use_raw=use_raw)
    else:
        sc.pl.umap(adata, color=color, ncols=ncols, palette=palette, save=save, color_map=color_map, size=size,
                   frameon=frameon, legend_fontsize=legend_fontsize, edges=edges, title=title, use_raw=use_raw)
    return None


def sc_paga(adata, color=None, save=None, frameon=False, threshold=0):
    sc.pl.paga(adata, color=color, save=save, frameon=frameon, threshold=threshold)
    return None


def sc_draw_graph(adata, color=None, save=None, color_map='magma', palette=None, size=50, legend_loc='right margin'):
    sc.pl.draw_graph(adata, color=color, save=save, color_map=color_map, palette=palette, size=size, legend_loc=legend_loc)
    return None


def sc_violin(adata, key=None, multi_panal=True, groupby=None, save='.png'):
    if key is None:
        key = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
    sc.pl.violin(adata, keys=key, groupby=groupby,
                 jitter=0.4, multi_panel=multi_panal, save=save)
    return None


def sc_matrixplot(adata, var_names, groupby='louvain_0.5', swap_axes=False, save='.png', dendrogram=False, vmax=2):
    sc.pl.matrixplot(adata, var_names=var_names, groupby=groupby, swap_axes=swap_axes, save=save,
                     dendrogram=dendrogram, vmax=vmax)
    return None


def sc_pca_variance_ratio(adata, log, save):
    sc.pl.pca_variance_ratio(adata, log=log, save=save)
    return None


def sc_scatter_with_line(x, y, c, x_label, y_label, plotdir, min_x=None, max_x=None, min_y=None, max_y=None, log=False):
    if log:
        if min_x is not None:
            min_x = np.log1p(min_x)
        if max_x is not None:
            max_x = np.log1p(max_x)
        if min_y is not None:
            min_y = np.log1p(min_y)
        if max_y is not None:
            max_y = np.log1p(max_y)
    plt.scatter(x, y, c=c, s=0.01)

    if min_x is not None:
        plt.axvline(x=min_x, color='r', linestyle='--', linewidth=1)  # min_counts
    if max_x is not None:
        plt.axvline(x=max_x, color='r', linestyle='--', linewidth=1)  # max_counts
    if min_y is not None:
        plt.axhline(y=min_y, color='r', linestyle='--', linewidth=1)  # min_genes
    if max_y is not None:
        plt.axhline(y=max_y, color='r', linestyle='--', linewidth=1)  # max_genes
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.colorbar()
    plt.savefig(plotdir)
    plt.close()
    return None


def plot_cell_frac_dpi(cell_list, adata, column_list, cell_type_palette, figname):
    df_list = []
    adata.obs['DPI'] = adata.obs['DPI'].astype(pd.CategoricalDtype(categories=column_list, ordered=True))
    new_column = []
    for dpi in column_list:
        fraction_list = []
        adata_sub = adata[adata.obs['DPI'] == dpi].copy()
        total_obs, total_vars = adata_sub.shape
        if total_obs <= 50:
            continue
        else:
            new_column.append(dpi)
        for cell_type in cell_list:
            adata_sub_cell = adata_sub[adata_sub.obs['annot'] == cell_type]
            cluster_obs, cluster_vars = adata_sub_cell.shape
            fraction_list.append(cluster_obs / total_obs * 100)
        df_list.append(fraction_list)
    df = pd.DataFrame(data=df_list)
    df.index = new_column
    df.columns = cell_list
    fig = df.plot.bar(stacked=True, color=cell_type_palette, grid=False, legend=False,
                      figsize=(6, 8), ylim=(0, 105)).get_figure()
    fig.tight_layout()
    fig.savefig(f'./figures/{figname}.png')
    plt.close()
    return None