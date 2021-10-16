# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 23:22:49 2020

@author: Administrator
"""
from typing import Union, Sequence, Mapping
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import matplotlib.cm
import seaborn as sns

#import utils_plot as uplt



# In[]

def drop_rows(df, names_drop):
    ''' temp functioin (drop a dataset)
    '''
    if isinstance(names_drop, str):
        names_drop = [names_drop]
    tmp = df.copy()
    for name in names_drop:
        if name in df.index.get_level_values(0):
            tmp = tmp.drop(name, level=0)
        if name in df.index.get_level_values(1):
            tmp = tmp.drop(name, level=1)
    return tmp



def sort_columns(df, by='mean'):
    col_ords = df.mean(0).sort_values().index
    return df[col_ords]
    
    
def sort_df(df):
    df = df.sort_values(by = list(df.columns)[::-1])
    return df

def get_rank(lst):
    order = np.argsort(lst)
    ranks = order.argsort()
    return ranks


def map_names_fuzzy(names, d: dict, fmt='{}'):
    def foo(nm):
        for k, v in d.items():
            if k in nm.lower():
                return fmt.format(v)
        return nm
    return foo(names) if isinstance(names, str) else list(map(foo, names))

def melt_df(df, cols=None, key_value='value', key_column='Method'):
    '''
    Un-pivot a dataframe.
    the output dataframne will have two columns: `key_value` and `key_column`
    
    Example:
    >>> ubp.melt_df(df, method_orders, 'Accuracy', 'Method')
    '''
    if cols is None:
        cols = df.columns
        
    dfs = []
    for c in cols:
        subdf = df[c].to_frame(key_value)
        subdf[key_column] = c
        dfs.append(subdf)
    df_melt = pd.concat(dfs)
    
    return df_melt

    
def merge_dfs(dfs: Union[Sequence[pd.DataFrame], Mapping], 
              idents=None, key_idents='ident'):
    ''' concatenate dataframes by indexes (vstack), and add a column named
    `key_ident` as the identity of each of th dataframes `dfs`
    '''
    
    if isinstance(dfs, (list, tuple)):
        if idents is None:
            idents = np.arange(len(dfs))
    
        df_dct = dict(zip(idents, dfs))
    elif isinstance(dfs, Mapping):
        df_dct = dfs
        
    df_list = []
    for idn, df in df_dct.items():
        tmp = df.copy()
        tmp[key_idents] = idn
        df_list.append(tmp)
        
    df_merge = pd.concat(df_list)
    return df_merge
        


# In[]
def _save_with_adjust(fig, fpath=None, figsize=None, **kwds):
    
    if figsize is not None:
        fig.set_size_inches(*figsize)
    if fpath is not None:
        fig.savefig(fpath, bbox_inches = 'tight', **kwds)
        print(f'figure has been saved into:\n\t{fpath}')
    else:
        fig.show()
    plt.close()


# In[]

'''
    visualization of the performace comparisons
'''


def plot_heatmap(df, 
                 norm_mtd = 'rank', norm_axis = 1,
                 cmap = 'RdBu_r',
                 figsize=(6, 3),
                 fn = None,
                 ylabel = 'pairs of datasets',
                 tt = None,
                 **kwds):
    
    if norm_mtd == 'rank':
        foo = get_rank
    elif norm_mtd == 'max':
        foo = lambda x: x / x.max()
    elif norm_mtd == 'minmax':
        foo = lambda x: (x - x.min()) / (x.max() - x.min())
    
    df = df.apply(foo, axis=norm_axis)
    
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(df, cmap = cmap, yticklabels='', 
                linewidths=0, 
    #            alpha = 0.75,
                cbar=False, ax=ax, **kwds)
    ax.set_ylabel(ylabel)
    if tt is not None:
        ax.set_title(tt)
    if fn is not None:
        fig.savefig(fn, bbox_inches = 'tight')
    return ax

        
def plot_boxes(df, 
               figsize=(6, 3),
               cmap='RdBu_r', 
               xrotation = 0, 
               saturation=1, 
               linewidth=1.5,
               ax = None,
               tt = None,
               ylabel=None,
               fn = None,
               **kwds):
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    #ax1 = sns.violinplot(data=tbl, palette=cmap_box, vmin=0, vmax=1,)
    sns.boxplot(data=df, palette=cmap, ax=ax, 
                saturation=saturation, linewidth=linewidth, **kwds)
    ax.set_yticks([0., 0.5, 1])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=xrotation, ha='right')
    if ylabel:
        ax.set_ylabel(ylabel)
    if tt is not None:
        ax.set_title(tt)
    if fn is not None:
        fig.savefig(fn, bbox_inches = 'tight')
        return ax
    
def plot_strips(df, 
               figsize=(6, 3),
               cmap='RdBu_r', 
               xrotation = 0, 
               alpha=0.5,
               marker='.',
               linewidth=1.5,
               tt = None,
               fn = None,
               **kwds):
    fig, ax = plt.subplots(figsize=figsize)
    #ax1 = sns.violinplot(data=tbl, palette=cmap_box, vmin=0, vmax=1,)
    sns.stripplot(data=df, palette=cmap, ax=ax, 
                alpha=alpha, marker=marker, 
                **kwds)
    ax.set_yticks([0., 0.5, 1])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=xrotation)
    if tt is not None:
        ax.set_title(tt)
    if fn is not None:
        fig.savefig(fn, bbox_inches = 'tight')
        return ax
    
    
def scatter_compare(df, 
                    x: str, y: str, 
                    hue=None, #'unknown_rate', 
                    ax = None,
                    figsize = (3, 3), # ignored if ax is not None
                    xlabel = None,
                    ylabel = None,
                    linewidth=0,
                    legend_loc = (1.03, 0),
                    text_loc = (0.5, 0.03), 
                    text_fmt = '{}',
                    text=None,
                    tt = None,
                    fn=None, 
                    **kwds):
    '''
    legend_loc: None for no legend
    '''
    if ax is None:
        fig, ax = plt.subplots(figsize = figsize)
    ax.plot([0, 1], [0, 1], linestyle='dashed', c = '#222831')
    sns.scatterplot(x = x, y = y, 
                    hue = hue,
                    data = df, ax = ax, 
                    linewidth=linewidth,
                    legend = 'auto', # 
                    **kwds)
    if hue:
        if legend_loc is None:
            ax.legend_.remove()
        else:
            ax.legend(*ax.get_legend_handles_labels(), loc = legend_loc)
    ax.set_xticks([0., 0.5, 1])
    ax.set_yticks([0., 0.5, 1])
    if text is not None:
        ax.text(*text_loc, text_fmt.format(text))
    
    if ylabel is not None:
        ax.set_ylabel(ylabel) #, size = 12
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    tt = '' if tt is None else tt
    ax.set_title(tt)    
    if fn is not None:
        ax.figure.savefig(fn, bbox_inches = 'tight')
    return ax
    
def scatter_compare_grid(df, x, y, 
                         hue=None,
                         tts = None,
                         figsize = (15, 3),
                         texts = None,
                         nrows = 1,
                         fn = None, 
                         gridspec_kw=None,
                         **kwds):
    '''
    gridspec_kw={'wspace': grid_wspace, 'hspace': grid_hspace}
    wspace : float, optional
    The amount of width reserved for space between subplots, expressed as a fraction of the average axis width.
    
    hspace : float, optional
    The amount of height reserved for space between subplots, expressed as a fraction of the average axis height.
    '''
    baselines = x if isinstance(y, str) else y
    n = len(baselines)
    if isinstance(tts, str) or tts is None:
        tts = [tts] * n
    if isinstance(texts, str) or texts is None:
        texts = [texts] * n
    ncols = n // nrows + min(n % nrows, 1)
    fig, axs = plt.subplots(
            nrows, ncols, figsize=figsize, sharey=True, 
            gridspec_kw=gridspec_kw)
    axs_flat = axs.flatten()
    for i, bm, txt, tt in zip(range(len(baselines)), baselines, texts, tts):
#        hue = 'unknown_rate' 
        ax = axs_flat[i]
        if hue not in df.columns: 
            hue = None
        scatter_compare(df, x = bm, y = y,
                        hue = hue,
                        ax = ax, tt = tt,
                        legend_loc=None if i < n - 1 else (1.03, 0),
                        text = txt,
                        **kwds)  
    if fn is not None:
        ax.figure.savefig(fn, bbox_inches = 'tight')
    return fig, axs


#def plot_umap




