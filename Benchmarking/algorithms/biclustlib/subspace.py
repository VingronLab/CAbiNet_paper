"""
    biclustlib: A Python library of biclustering algorithms and evaluation measures.
    Copyright (C) 2017  Victor Alexandre Padilha

    This file is part of biclustlib.

    biclustlib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    biclustlib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import sys
import numpy as np

from munkres import Munkres
from itertools import product


def clustering_error(predicted_biclustering, reference_biclustering, num_rows, num_cols):
    """The Clustering Error (CE) external evaluation measure.

    CE computes the similarity between two subspace clusterings. This measure was originally
    introduced in (Patrikainen and Meila, 2006) as a dissimilarity measure. In this package, it
    was implemented as a similarity measure as presented in (Horta and Campello, 2014). This measure
    lies in the interval [0, 1], where values close to 1 indicate better biclustering solutions.

    Reference
    ---------
    Patrikainen, A., & Meila, M. (2006). Comparing subspace clusterings. IEEE Transactions on
    Knowledge and Data Engineering, 18(7), 902-916.

    Horta, D., & Campello, R. J. G. B. (2014). Similarity measures for comparing biclusterings.
    IEEE/ACM Transactions on Computational Biology and Bioinformatics, 11(5), 942-954.

    Parameters
    ----------
    predicted_biclustering : biclustlib.model.Biclustering
        Predicted biclustering solution.

    reference_biclustering : biclustlib.model.Biclustering
        Reference biclustering solution.

    num_rows : int
        Number of rows of the dataset.

    num_cols : int
        Number of columns of the dataset.

    Returns
    -------
    ce : float
        Similarity score between 0.0 and 1.0.
    """
    

    check = check_biclusterings(predicted_biclustering, reference_biclustering)

    if isinstance(check, float):
        return check

    union_size = _calculate_size(predicted_biclustering, reference_biclustering, num_rows, num_cols, 'union')
    dmax = _calculate_dmax(predicted_biclustering, reference_biclustering)

    return float(dmax) / union_size

def relative_non_intersecting_area(predicted_biclustering, reference_biclustering, num_rows, num_cols):
    """The Relative Non-Intersecting Area (RNIA) external evaluation measure.

    RNIA computes the similarity between two subspace clusterings. This measure was originally
    introduced in (Patrikainen and Meila, 2006) as a dissimilarity measure. In this package, it
    was implemented as a similarity measure as presented in (Horta and Campello, 2014). This measure
    lies in the interval [0, 1], where values close to 1 indicate better biclustering solutions.

    Reference
    ---------
    Patrikainen, A., & Meila, M. (2006). Comparing subspace clusterings. IEEE Transactions on
    Knowledge and Data Engineering, 18(7), 902-916.

    Horta, D., & Campello, R. J. G. B. (2014). Similarity measures for comparing biclusterings.
    IEEE/ACM Transactions on Computational Biology and Bioinformatics, 11(5), 942-954.

    Parameters
    ----------
    predicted_biclustering : biclustlib.model.Biclustering
        Predicted biclustering solution.

    reference_biclustering : biclustlib.model.Biclustering
        Reference biclustering solution.

    num_rows : int
        Number of rows of the dataset.

    num_cols : int
        Number of columns of the dataset.

    Returns
    -------
    rnia : float
        Similarity score between 0.0 and 1.0.
    """
    check = check_biclusterings(predicted_biclustering, reference_biclustering)

    if isinstance(check, float):
        return check

    union_size = _calculate_size(predicted_biclustering, reference_biclustering, num_rows, num_cols, 'union')
    intersection_size = _calculate_size(predicted_biclustering, reference_biclustering, num_rows, num_cols, 'intersection')

    return float(intersection_size) / union_size

def _calculate_size(predicted_biclustering, reference_biclustering, num_rows, num_cols, operation):
    pred_count = _count_biclusters(predicted_biclustering, num_rows, num_cols)
    true_count = _count_biclusters(reference_biclustering, num_rows, num_cols)

    if operation == 'union':
        return np.sum(np.maximum(pred_count, true_count))
    elif operation == 'intersection':
        return np.sum(np.minimum(pred_count, true_count))

    valid_operations = ('union', 'intersection')

    raise ValueError("operation must be one of {0}, got {1}".format(valid_operations, operation))

def _calculate_dmax(predicted_biclustering, reference_biclustering):
    pred_sets = _bic2sets(predicted_biclustering)
    true_sets = _bic2sets(reference_biclustering)
    cost_matrix = [[sys.maxsize - len(b.intersection(g)) for g in true_sets] for b in pred_sets]
    indices = Munkres().compute(cost_matrix)
    return sum(sys.maxsize - cost_matrix[i][j] for i, j in indices)

def _count_biclusters(biclustering, num_rows, num_cols):
    count = np.zeros((num_rows, num_cols), dtype=np.int)

    for b in biclustering.biclusters:
        if b.rows.ndim == 0 or b.cols.ndim == 0:
            count[b.rows, b.cols] += 1
        else:
            count[np.ix_(b.rows, b.cols)] += 1

    return count

# def _bic2sets(biclust):
#     return [set(product(b.rows, b.cols)) for b in biclust.biclusters]
def _bic2sets(biclust):
    prd = []
    for b in biclust.biclusters:
        if b.cols.ndim == 0:
            b.cols = np.array([b.cols])
        if b.rows.ndim == 0:
            b.rows = np.array([b.rows])
            
        prd.append(set(product(b.rows, b.cols)))

    return prd
    
 ############################# check.py
    

def check_biclusterings(b1, b2):
    if not isinstance(b1, Biclustering) or not isinstance(b2, Biclustering):
        raise ValueError('b1 and b2 must be Biclustering instances')

    if len(b1.biclusters) == 0 and len(b2.biclusters) == 0:
        return 1.0

    if len(b1.biclusters) == 0 or len(b2.biclusters) == 0:
        return 0.0

    return None
    
    
    
  ############ models.py
  
  

class Bicluster:
    """This class models a bicluster.

    Parameters
    ----------
    rows : numpy.array
        Rows of the bicluster (assumes that row indexing starts at 0).

    cols : numpy.array
        Columns of the bicluster (assumes that column indexing starts at 0).

    data : numpy.ndarray
        bla
    """

    def __init__(self, rows, cols, data=None):
        if isinstance(rows, np.ndarray) and rows.dtype == np.bool and cols.dtype == np.bool:
            self.rows = np.nonzero(rows)[0]
            self.cols = np.nonzero(cols)[0]
        elif isinstance(cols, np.ndarray) and rows.dtype == np.int and cols.dtype == np.int:
            self.rows = rows
            self.cols = cols
        else:
            raise ValueError("rows and cols must be bool or int numpy.arrays")

        if data is not None:
            n, m = len(self.rows), len(self.cols)

            if isinstance(data, np.ndarray) and (data.shape == (n, m) or (len(data) == 0 and n == 0)):
                self.data = data
            else:
                raise ValueError("")

    def intersection(self, other):
        """Returns a bicluster that represents the area of overlap between two biclusters."""
        rows_intersec = np.intersect1d(self.rows, other.rows)
        cols_intersec = np.intersect1d(self.cols, other.cols)
        return Bicluster(rows_intersec, cols_intersec)

    def union(self, other):
        rows_union = np.union1d(self.rows, other.rows)
        cols_union = np.union1d(self.cols, other.cols)
        return Bicluster(rows_union, cols_union)

    def overlap(self, other):
        min_area = min(self.area, other.area)
        return self.intersection(other).area / min_area

    @property
    def area(self):
        """Calculates the number of matrix elements of the bicluster."""
        return len(self.rows) * len(self.cols)

    def sort(self):
        """Sorts the array of row and the array of column indices of the bicluster."""
        self.rows.sort()
        self.cols.sort()

    def __str__(self):
        return 'Bicluster(rows={0}, cols={1})'.format(self.rows, self.cols)


class Biclustering:
    """This class models a biclustering.

    Parameters
    ----------
    biclusters : list
        A list of instances from the Bicluster class.
    """

    def __init__(self, biclusters):
        if all(isinstance(b, Bicluster) for b in biclusters):
            self.biclusters = biclusters
        else:
            raise ValueError("biclusters list contains an element that is not a Bicluster instance")

    def __str__(self):
        return '\n'.join(str(b) for b in self.biclusters)
        
        
def convert2biclustlib(biclustering):

    biclusters = []
    # f = open("biclustering", "wb")
    # pickle.dump(biclustering, f)

    # print("biclustering", flush = True)

    for i in range(len(biclustering)):

        row_indices, col_indices = biclustering[i]

        row_indices = np.asarray(row_indices)
        col_indices = np.asarray(col_indices)

        if row_indices.size == 0 or col_indices.size == 0:
            break

        biclusters.append(Bicluster(row_indices, col_indices))
    
    biclustering = Biclustering(biclusters)
    
    return(biclustering)
