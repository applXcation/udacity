# -*- coding: utf-8 -*-

"""
Misc functions
"""
# Facilities en provenance de Python 3
from __future__ import print_function, division#, unicode_literals
# Statistiques
import random



def createFoldIndexes(index, K, seed):
  """
  Folds generation from an index.

  Parameters
  ----------
  index : 1-d array / list
    List of indexes to be splitted into folds.

  K : int
    Number of folds.

  seed : int
    Seeds.


  Returns
  -------
  foldsIndex : list of length K
    For each item k, contains a dictionnary with item "in" being the in-fold
    indexes (1/K-th of the data), and "out" the out-fold indexes (the rest).
  """

  # EXTRACTION D'INFO
  n = len(index)
  chunkSize = int(n / K)

  # RANDOMISATION
  random.seed(seed)
  index = random.sample(index, n)

  # BOUCLE
  foldsIndex = []
  i = 0
  for k in range(1, K + 1):
    if k != K:
      kthFold = {"out" : index[i:(i + chunkSize)]}
      kthFold["in"] = index[:i] + index[(i + chunkSize):]
    else: 
      kthFold = {"out" : index[i:]}
      kthFold["in"] = index[:i]
    foldsIndex.append(kthFold)
    i += chunkSize

  # RETOUR
  return foldsIndex


