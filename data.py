# -*- coding: utf-8 -*-

"""
Synthetic data
"""
# Facilities en provenance de Python 3
from __future__ import print_function, division, unicode_literals
# Classes d'objets
import numpy as np
# Stats
import scipy.stats as stats



def sugyData(N, muTrain = 1, sigmaTrain = .5, muTest = 2, sigmaTest = .25,
    seed = None):

  """
  Generates synthetic data.
  
      
  Parameters
  ----------
  muTrain, muTest : float
      E(x) in train and test respectively.

  sigmaTrain, sigmaTest : float
      sd(x) in train and test respectively.

  N : integer
      Sample size.

  seed: integer
      seed.


  Returns
  -------
  Dictionary containing the data (xTrain, xTest, yTrain, yTest) and the
  theoretical ratio f_test(x)/f_train(x) for every observation in the train and
  test samples (wTheoriqueTrain, wTheoriqueTest).

  """

  # SEEDS
  np.random.seed(seed)


  # X
  # Donnees
  xTrain = np.random.normal(muTrain, sigmaTrain, size = N)
  xTest = np.random.normal(muTest, sigmaTest, size = N)

  # Reshape vers matrice-colonne
  xTrain = xTrain.reshape(-1, 1) 
  xTest = xTest.reshape(-1, 1) 


  # RATIO THEORIQUE
  x = np.concatenate((xTrain, xTest), axis = 0)

  # Calcul du ratio théorique f_test(x)/f_train(x)
  distTrain = stats.norm(muTrain, sigmaTrain)
  distTest = stats.norm(muTest, sigmaTest)
  fTrain = distTrain.pdf(x)
  fTest = distTest.pdf(x)

  # Ratios
  wTheorique = fTest / fTrain
  wTheoriqueTrain = wTheorique[:N]
  wTheoriqueTest = wTheorique[N:]


  # Y
  # Epsilon
  epsTrain = np.random.normal(0, 1/4, size = N).reshape(-1, 1)
  epsTest = np.random.normal(0, 1/4, size = N).reshape(-1, 1)

  # y
  yTrain = np.sinc(xTrain) + epsTrain
  yTest = np.sinc(xTest) + epsTest


  # RETOUR
  retour = {
    'xTrain' : xTrain,
    'yTrain' : yTrain,
    'wTheoriqueTrain' : wTheoriqueTrain, 
    'xTest' : xTest,
    'yTest' : yTest,
    'wTheoriqueTest' : wTheoriqueTest
    }

  return retour
