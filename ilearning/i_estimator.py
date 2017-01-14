# -*- coding: utf-8 -*-

"""
Estimation of f_test(x) / f_train(x)
"""
# Facilities en provenance de Python 3
from __future__ import print_function, division#, unicode_literals
# Reference interne
from ilearning.misc import createFoldIndexes
# Classes d'objet
import numpy as np
import pandas as pd
from pandas import Series, DataFrame
# Statistiques
import scipy.stats as stats
# ML
from sklearn.grid_search import GridSearchCV
from sklearn.calibration import CalibratedClassifierCV
from sklearn.base import clone
# Misc
import time
import copy



class importanceRatioEstimator(object):
  """
  Estimate f_test(x) / f_train(x).

  Parameters
  ----------
  scorer : classifier
      The estimator used for scoring the probability of belonging to the test set.
      
  K : int
      Number of folds for spliting the data into probability predictor sets / prediction sets.

  percCal : float
      Percentage of scorer training folds used for calibrating the scorers.

  calMethod : character
    Calibration method ("None", "isotonic", "sigmoid").

  grid : dict
    Grid for cross-validating the scorers' metaparameters.

  kCVScorer : int
    Number of folds for cross-validating the scorers' metaparameters.

  nthreads : int
    Number of threads.

  seed : None or int
    Seed.


  Attributes
  ----------
  Input arguments are attached as attributes to the object. In addition:

  wTrain, wTest: weights for xTrain, ordered as the xTrain fed in the .fit
  method (resp xTest)

  w : pandas (len(xTrain) + len(xTest), 2)
    Concatenation of xTrain and xTest, with an additional column stating
    whether the weight comes from the train or test sample.

  grid_scores_ : dict
    Grid of scores obtained from cross-validating the scorer.

  best_params_ : dict
    best parameters obtained from cross-validating the scorer.
  

  TODO
  ----
  * Grid Search CV for the scorer parameter is mandatory ; should be an option.
    Also, gridsearch should be done separately as a method.

  * w is deprecated ; to be dropped.
  """

  def __init__(self, scorer, K, percCal, calMethod, grid, kCVScorer, nthreads, seed):

    # Input
    self.scorer = scorer
    self.K = K
    self.percCal = percCal
    self.calMethod = calMethod
    self.grid = grid
    self.kCVScorer = kCVScorer
    self.nthreads = nthreads
    self.seed = seed
    # Attributs privés
    self._scorerCalibration = None
    self._scorers = []
    self._probaFuns = []
    self._correcIntercept = []
    self._preds = []
    # Output
    self.wTrain = None
    self.wTest = None
    self.w = None
    self.grid_scores_ = None
    self.best_params_ = None



  def fitAndPred(self, xTrain, xTest):
    """
    Learns weight from and for xTrain and xTest. Sets values for self.wTrain,
    self.wTest and self.w.

    Arugments
    ---------
    xTrain, xTest : datasets that may be fed to sklearn .fit methods
      Train and test samples.


    """
    # INITIALISATION
    # Parametres
    scorer = self.scorer
    K = self.K
    percCal = self.percCal
    calMethod = self.calMethod
    grid = self.grid
    kCVScorer = self.kCVScorer
    nthreads = self.nthreads
    seed = self.seed

    # Seed setting
    np.random.seed(seed)  # utilisé par GridSearchCV


    # PREPARATION DONNEES
    # concatenation de xTrain et xTest, et creation d'un vecteur indiquant pour
    # chaque ligne de quel jeu elle provient ("train", "test")
    x, y = self._dataPrep(xTrain = xTrain, xTest = xTest)


    # PREPARATION FOLDS
    # Creation indexes de folds pour apprendre fonction de proba / prédire les
    # proba
    foldIds = createFoldIndexes(index = range(x.shape[0]), K = K, seed = seed)

    # Subsplit : pour chaque fold d'apprentissage, on split en un fold
    # d'apprentissage de la fonction de scoring, et un de la fonction de
    # calibration.
    foldIds = self._subsplitFoldIndexes(K = K, foldIds = foldIds, calMethod =
        calMethod, percCal = percCal)


    # ENTRAINEMENT PREMIER SCOREUR
    # Entraînement du premier score avec selection des metapara par crossval
    ids = foldIds[0]["in"]["scoreLearning"]
    scorerCV = GridSearchCV(estimator = scorer, param_grid = grid, 
      cv = kCVScorer, n_jobs = nthreads)
    cvResult = scorerCV.fit(X = x[ids, :], y = y[ids])


    # STOCKAGE DANS SELF
    self.grid_scores_ = cvResult.grid_scores_
    self.best_params_ = cvResult.best_params_
    bestScorer = cvResult.best_estimator_
    self._scorers.append(bestScorer)
    self._scorerCalibration = cvResult


    # ENTRAINEMENT K-1 SCOREURS SUIVANTS
    for k in range(1, K):
      kthScorer = clone(bestScorer)
      ids = foldIds[k]["in"]["scoreLearning"]
      kthScorer.fit(X = x[ids, :], y = y[ids])
      self._scorers.append(kthScorer)


    # CALIBRATION DES K SCOREURS
    # Calibration
    for k in range(K):
      probaFun, correcIntercept = self._calibration(foldIds = foldIds[k],
        scorer = self._scorers[k], calMethod = calMethod, x = x, y = y)
      self._probaFuns.append(probaFun)
      self._correcIntercept.append(correcIntercept)

    # Nettoyage
    del probaFun, correcIntercept


    # PREDICTION
    # Récupération des prédictions pour chaque jeu
    for k in range(K):
      self._preds.append(self._ratioPredict(foldIds = foldIds[k],
        probaFun = self._probaFuns[k], correc = self._correcIntercept[k],
        x = x))

    # Mise bout-à-bout dans le même ordre que l'index
    w = []
    ids = []
    for k in range(K):
      w.extend(self._preds[k])
      ids.extend(foldIds[k]["out"])
    w = DataFrame({"wPred" : w, "ids" : ids})
    w.sort_values(by = "ids", inplace = True)
    w["y"] = y


    # CORRECTION VALEURS EXTRÊMES
    # Pour les valeurs de wPred = Inf, on les remplace par la valeur maximale 
    # de wPred (dans le train et dans le test, respectivement).
    # Et inversement pour wPred <= 0 (on remplace par min)
    w = self._extremeRatioCorrec(w)

    # Stockage
    self.w = w
    self.wTrain = w.loc[w.y == "train", "wPred"]
    self.wTest = w.loc[w.y == "test", "wPred"]



  def _dataPrep(self, xTrain, xTest):
    """
    From a training and a test datasets, returns x (concatenation) and y
    (whether the observation in x belongs to train or test data) together in a
    tuple (x, y)
    """
    # On concatene train et test
    x = np.concatenate((xTrain, xTest), axis = 0)

    # Si x a une dimension seulement, on le passe en deux dimension
    # (necessaire pour sklearn)
    if len(x.shape) <= 1:
      x = x.reshape(-1, 1) 

    # On génère un vecteur binaire, 0 si train, 1 si test
    y = np.append(
      np.repeat("train", xTrain.shape[0]),
      np.repeat("test", xTest.shape[0]))

    return (x, y)



  def _subsplitFoldIndexes(self, K, foldIds, calMethod, percCal):
    """
    Subplits fold Ids into scoring / calibration folds.

    """
    for k in range(0, K):
      kthFoldInOld = foldIds[k]["in"]

      # Si calibration quelconque :
      if calMethod != 'None':
        nScorer = int(len(kthFoldInOld) * (1 - percCal))

        # Resampling des anciens ids
        index = np.random.choice(xrange(len(kthFoldInOld)), len(kthFoldInOld),
          replace = False)
        kthFoldInOldRan = [kthFoldInOld[i] for i in index]

        # Sélection
        kthFoldInNew = {"scoreLearning" : kthFoldInOldRan[:nScorer]}
        kthFoldInNew["calLearning"] = kthFoldInOldRan[nScorer:]
        foldIds[k]["in"] = kthFoldInNew

      # Si pas de calibration : 
      else:
        kthFoldInNew = {"scoreLearning" : kthFoldInOld, "calLearning" : []}
        foldIds[k]["in"] = kthFoldInNew

    return foldIds



  def _calibration(self, foldIds, scorer, calMethod, x, y):
    """
    Returns calibrated scorer (initial scorer if calMethod == "None") and
    correction intercept as as tuple (probaFun, correcIntercept)
    """
    # Calibration
    if calMethod != 'None':
      ids = foldIds["in"]["calLearning"]
      calib = CalibratedClassifierCV(
        base_estimator = scorer,
        method = calMethod,
        cv = "prefit")
      calib.fit(
        X = x[ids, :],
        y = y[ids])
      probaFun = calib
    else:
      probaFun = scorer

    # Apprentissage de la constante de correction
    if calMethod != 'None':
      ids = foldIds["in"]["calLearning"]
    else:
      ids = foldIds["in"]["scoreLearning"]
    correc = sum(y[ids] == "train") / sum(y[ids] == "test")
    correcIntercept = correc

    return (probaFun, correcIntercept)



  def _ratioPredict(self, foldIds, probaFun, correc, x):
    """
    Predicts probability ratio.
    """
    ids = foldIds["out"]
    probaFun = probaFun
    correc = correc

    # Quelle colonne des prédictions contient la proba de test ?
    col = probaFun.classes_.tolist().index("test")

    # Prédiction: proba d'appartenir au test sachant x
    proba = probaFun.predict_proba(X = x[ids])[:, col]

    # Calcul w(x)
    w = proba / (1 - proba) * correc

    # Retour
    return w



  def _extremeRatioCorrec(self, w):
    """
    Correcting for extreme values
    """
    # wPred == Inf
    locExtr = (w.wPred == np.inf)
    locTrain = (w.y == "train")
    w.loc[locExtr & locTrain, "wPred"] = w.loc[(~locExtr) & locTrain, "wPred"] .max()
    w.loc[locExtr & (~locTrain), "wPred"] = w.loc[(~locExtr) & (~locTrain), "wPred"].max()

    ## wPred <= 0
    locExtr = (w.wPred <= 0)
    locTrain = (w.y == "train")
    w.loc[locExtr & locTrain, "wPred"] = w.loc[(~locExtr) & locTrain, "wPred"].min()
    w.loc[locExtr & (~locTrain), "wPred"] = w.loc[(~locExtr) & (~locTrain), "wPred"].min()

    # Retour
    return w


