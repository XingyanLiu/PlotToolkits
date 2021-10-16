# -*- coding: UTF-8 -*-
"""
@Author: Xingyan Liu
@CreateDate: 2021-10-16
@File: compute.py
@Project: PlotToolkits
"""
import os
from pathlib import Path
from typing import Union, Optional, Sequence, Mapping
import time
import logging
import numpy as np
import pandas as pd
from scipy import sparse


def wrapper_confus_mat(y_true, y_pred, classes_on=None,
                       normalize='true', as_df=True):
    """
    normalize: 'true', 'pred', 'all', None
        by default, normalized by row (true classes)
    """
    from sklearn import metrics
    if classes_on is None:
        classes_on = np.unique(list(y_true) + list(y_pred))
    try:
        mat = metrics.confusion_matrix(y_true, y_pred, labels=classes_on,
                                       normalize=normalize)
    except:
        logging.warning(
            'The argument `normalize` may not be accepted by '
            'the previous version of scikit-learn')
        mat = metrics.confusion_matrix(y_true, y_pred, labels=classes_on, )
    if as_df:
        mat = pd.DataFrame(mat, index=classes_on, columns=classes_on)
    return mat


def wrapper_contingency_mat(y_true, y_pred,
                            order_rows=True,
                            order_cols=False,
                            normalize_axis=None,
                            as_df=True,
                            eps=None,
                            assparse=False
                            ):
    """
    Modified and wrapped function from `sklearn`:
    >>> mat = sklearn.metrics.cluster.contingency_matrix(
    ...        y_true, y_pred, eps=eps, sparse=assparse)
    """
    if eps is not None and sparse:
        raise ValueError("Cannot set 'eps' when sparse=True")

    classes, class_idx = np.unique(y_true, return_inverse=True)
    clusters, cluster_idx = np.unique(y_pred, return_inverse=True)
    n_classes = classes.shape[0]
    n_clusters = clusters.shape[0]
    # Using coo_matrix to accelerate simple histogram calculation,
    # i.e. bins are consecutive integers
    # Currently, coo_matrix is faster than histogram2d for simple cases
    mat = sparse.coo_matrix(
        (np.ones(class_idx.shape[0]), (class_idx, cluster_idx)),
        shape=(n_classes, n_clusters), dtype=np.int
    )
    if assparse:
        mat = mat.tocsr()
        mat.sum_duplicates()
        if normalize_axis is not None:  # 0 for columns
            mat = normalize_norms(mat, axis=normalize_axis)
    else:
        mat = mat.toarray()
        if eps is not None:
            # don't use += as mat is integer
            mat = mat + eps
        if normalize_axis is not None:  # 0 for columns
            mat = normalize_norms(mat, axis=normalize_axis)

        if as_df:
            mat = pd.DataFrame(mat, index=classes, columns=clusters)
        # reorder to make clusters and classes matching each other as possible
        if order_cols:
            mat = order_contingency_mat(mat, 0)
        if order_rows:
            mat = order_contingency_mat(mat, 1)
    return mat


def normalize_col(X, scale_factor=1., by='sum'):
    """
    make the column elements of X to unit sum

    Parameters
    ----------
    X:
        a (sparse) matrix
    scale_factor: float or None
        if None, use the median of sum level as the scaling factor.
    by: str, {'sum', 'max'}

    """
    if by == 'sum':
        norms = X.sum(axis=0)
    elif by == 'max':
        norms = X.max(axis=0)
    else:
        raise ValueError(f'`by` should be either "sum" or "max", got {by}')

    if scale_factor is None:
        is_zero = norms == 0
        scale_factor = np.median(norms[~ is_zero])
    norms = norms / scale_factor
    # for those rows or columns that summed to 0, just do nothing
    if hasattr(norms, 'A'):
        norms = norms.A.flatten()
    norms[norms == 0] = 1

    norm_ = 1 / norms

    if sparse.isspmatrix(X):
        logging.info('sparse normalization')
        X_new = X.dot(sparse.diags(norm_))
    else:
        logging.info('dense normalization')
        X_new = X.dot(np.diag(norm_))

    if isinstance(X, pd.DataFrame):
        X_new.columns = X.columns
    return X_new


def normalize_row(X, scale_factor=1, by='sum'):
    """
    make the row elements of X to unit sum

    Parameters
    ----------
    X:
        a (sparse) matrix
    scale_factor: float or None
        if None, use the median of sum level as the scaling factor.
    by: str, {'sum', 'max'}

    """
    if by == 'sum':
        norms = X.sum(axis=1)
    elif by == 'max':
        norms = X.max(axis=1)
    else:
        raise ValueError(f'`by` should be either "sum" or "max", got {by}')

    if scale_factor is None:
        is_zero = norms == 0
        scale_factor = np.median(norms[~ is_zero])
    norms = norms / scale_factor
    # for those rows or columns that summed to 0, just do nothing
    if hasattr(norms, 'A'): norms = norms.A.flatten()
    norms[norms == 0] = 1
    norm_ = 1 / norms

    if sparse.isspmatrix(X):
        logging.info('sparse normalization')
        X_new = sparse.diags(norm_).dot(X)
    else:
        logging.info('dense normalization')
        X_new = np.diag(norm_).dot(X)

    if isinstance(X, pd.DataFrame):
        X_new.columns = X.columns
    return X_new


def normalize_norms(X, scale_factor=1, axis=0, by='sum'):
    """ wrapper of `normalize_colsum` and `normalize_rowsum`

    Parameters
    ----------
    X:
        a (sparse) matrix
    scale_factor: numeric, None
        if None, use the median of sum level as the scaling factor.
    axis: int, {0, 1}
        if axis = 0, apply to each column;
        if axis = 1, apply to each row.
    by: str, {'sum', 'max'}
        normalization method

    """
    foo = normalize_col if axis == 0 else normalize_row
    return foo(X, scale_factor=scale_factor, by=by)


def zscore(X, with_mean=True, scale=True, ):
    """ For each column of X, do centering (z-scoring)
    """
    # code borrowed from `scanpy.pp._simple`
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler(with_mean=with_mean, copy=True).partial_fit(X)
    if scale:
        # user R convention (unbiased estimator)
        e_adjust = np.sqrt(X.shape[0] / (X.shape[0] - 1))
        scaler.scale_ *= e_adjust
    else:
        scaler.scale_ = np.array([1] * X.shape[1])
    X_new = scaler.transform(X)
    if isinstance(X, pd.DataFrame):
        X_new = pd.DataFrame(X_new, index=X.index, columns=X.columns)
    return X_new

def _order_contingency_array(mat, axis=1):
    """
    axis = 0: re-order the columns
    axis = 1: re-order the rows
    """
    order = np.argsort(np.argmax(mat, axis=axis))
    if axis == 1:
        print('Re-order the rows')
        return mat[order, :]
    else:
        print('Re-order the columns')
        return mat[:, order]


def _order_contingency_df(df: pd.DataFrame, axis=1):
    """
    axis = 0: re-order the columns
    axis = 1: re-order the rows
    """
    order = np.argsort(np.argmax(df.values, axis=axis))
    if axis == 1:
        print('Re-order the rows')
        return df.iloc[order, :]
    else:
        print('Re-order the columns')
        return df.iloc[:, order]


def order_contingency_mat(mat, axis=1):
    if isinstance(mat, pd.DataFrame):
        return _order_contingency_df(mat, axis=axis)
    else:
        return _order_contingency_array(mat, axis=axis)


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
