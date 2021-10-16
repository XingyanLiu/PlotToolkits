# -*- coding: UTF-8 -*-
"""
@Author: Xingyan Liu
@CreateDate: 2021-10-17
@File: base.py
@Project: PlotToolkits
"""
import os
from pathlib import Path
from typing import Union, Optional, Sequence, Mapping
import time
import logging
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


def _save_with_adjust(fig, fpath=None, figsize=None, **kwds):
    if figsize is not None:
        fig.set_size_inches(*figsize)
    if fpath is not None:
        fig.savefig(fpath, bbox_inches='tight', **kwds)
        print(f'figure has been saved into:\n\t{fpath}')
    else:
        fig.show()
    plt.close()


def rotate_xticklabels(ax, angle=45, **kwargs):
    ax.set_xticklabels(ax.get_xticklabels(), rotation=angle,
                       **kwargs)


def rotate_yticklabels(ax, angle=45, **kwargs):
    ax.set_yticklabels(ax.get_yticklabels(), rotation=angle,
                       **kwargs)


# In[]
# colors
def view_color_map(cmap='viridis', n=None, figsize=(6, 2), s=150, k=20,
                   ax=None,
                   grid=False, **kwds):
    """
    n: total number of colors
    k: number of colors to be plotted on each line.

    Examples
    --------
    import funx as fx
    colors = ['Set1', 'viridis', 'Spectral']
    fx.view_color_map(colors[-1], n=20)

    cmap = sc.pl.palettes.zeileis_26
    cmap = sc.pl.palettes.default_64
    fx.view_color_map(cmap, k=16)
    """
    if not isinstance(cmap, (np.ndarray, list)):
        #        from matplotlib import cm
        cmp = mpl.cm.get_cmap(cmap, n)
        n = 20 if n is None else n
        cmp = [cmp(i) for i in range(n)]
    else:
        n = len(cmap) if n is None else n
        cmp = cmap
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    for i in range(n):
        ax.scatter(i % k, i // k, color=cmp[i], s=s)
    plt.grid(b=grid, )  # axis='y')
    plt.show()
    return ax


def get_colors(cmap='Spectral', n=5, to_hex=True):
    """
    fx.get_colors('Reds', 4)
    fx.view_color_map('Reds', 4)
    """
    #    import matplotlib.colors as mcolors
    cmp = plt.cm.get_cmap(cmap, n)
    colors = [cmp(i) for i in range(n)]
    if to_hex:
        return [mcolors.to_hex(c) for c in colors]
    else:
        return colors


def __test__():
    pass


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(filename)s-%(lineno)d-%(funcName)s(): '
               '%(levelname)s\n %(message)s')
    t = time.time()

    __test__()

    print('Done running file: {}\nTime: {}'.format(
        os.path.abspath(__file__), time.time() - t,
    ))
