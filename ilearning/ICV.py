# -*- coding: utf-8 -*-

"""
ICV
"""
# Facilities en provenance de Python 3
from __future__ import print_function, division#, unicode_literals
# Classes d'objet
import numpy as np
# ML
from skl_scorerweights.grid_search import GridSearchCV
from skl_scorerweights.base import clone



class importanceWeightedCV(object):
  """
  I-cross validation.

  Parameters
  ----------
  weightExpList: list of numeric values, default [1]
    List of exponents for the the training sample weights.

  weightFit: boolean, default True
    Should the `sample_weight` be used when fitting models?

  weightEval: boolean, default True
    Should the `sample_weight` be used when evaluatting models?

  **kwargs: arguments to be passed to sklearn.grid_search.GridSearchCV. 


  Attributes
  ----------
  Input arguments are attached as attributes to the object. In addition:

  best_estimator_: estimator from sklearn
    Best estimator trained.

  best_gamma_: numeric value
    Best exponent for the training weights.

  best_params_: dict
    Metaparameters used to learn the best estmator.

  best_score_: numeric value
    Best scored obtained among all exponent / model parameter values tested.

  best_scores_: list
    For each value of weightExpList, best score among all parameters
    value for this exponent value.

  cv_list_: list
    For each value of weightExpList, stores the GridSearchCV object.

  grid_scores_: list
    For each value of weightExpList, gives the grid_scores_ associated to
    the cross-validation for this exponent value.
  """

  def __init__(self, weightExpList = [1], weightFit = True,
    weightEval = True, **kwargs):

    self.weightExpList = weightExpList
    self.weightFit = weightFit
    self.weightEval = weightEval
    self.kwargs = kwargs if kwargs is not None else {}
    # Output
    best_estimator_ = None
    best_gamma_ = None
    best_score_ = None
    best_params = None
    grid_scores_ = []
    best_scores_ = []
    cv_list_ = []



  def fit(self, X, y, sample_weight):
    """
    Fit method for importance_weighted cross-validation

    Parameters
    ----------
    
    X : array-like, shape = [n_samples, n_features]
      Training vector, where n_samples is the number of samples and
      n_features is the number of features.

    y : array-like, shape = [n_samples] or [n_samples, n_output], optional
      Target relative to X for classification or regression;
      None for unsupervised learning.

    sample_weight: numpy array of shape [n_samples]
      Sample weights, as provided by `importanceRatioEstimator.wTrain`

    TODO
    ----
    
    * Sanity check: vérifier qu'aucun objet type refit n'est passé dans le
      kwargs

    """

    # INITIALISATIONS
    # Pointeurs rapides
    gammaList = self.weightExpList
    
    # Objets supplémentaires
    ## Listes d'output à renvoyer au self
    grid_scores_ = []
    best_scores_ = []
    cv_list_ = []


    # CV POUR CHAQUE GAMMA
    for i in range(len(gammaList)):
      scorerCV = self._cvOverGamma(gamma = gammaList[i], sample_weight = 
        sample_weight, weightFit = self.weightFit, weightEval = 
        self.weightEval, kwargs = self.kwargs, X = X, y = y)
      
      # Stockage
      grid_scores_.append(scorerCV.grid_scores_)
      best_scores_.append(scorerCV.best_score_)
      cv_list_.append(scorerCV)


    # INFERENCE
    # Quelle est l'exposant associé au meilleur estimateur ?
    # Si on a plusieurs fois le même score maximal, alors on prendra le plus
    # petit exposant possible (==> max sparsity)
    bestIndex = best_scores_.index(np.max(best_scores_))
    best_gamma_ = gammaList[bestIndex]

    # Stockage résultats
    best_score_ = best_scores_[bestIndex]
    best_params_ = cv_list_[bestIndex].best_params_


    # FIT
    # Fit du modele avec meilleur gamma et meilleurs metapara 
    # sur toutes les donnees

    # On recupere le modele donné en argument, et on lui donne les meilleurs
    # metapara
    model = clone(self.kwargs["estimator"])
    model.set_params(**best_params_)

    # On entraine le modele avec la meilleure exponentiation
    ## Les poids ne sont utilisés que si besoin, pour une raison déjà évoquée
    ## lors de la cross-validation des modèles.
    wTrain = sample_weight**best_gamma_ 
    if self.weightFit and best_gamma_ != 0:
      best_estimator_ = model.fit(X = X, y = y, sample_weight = wTrain)
    else:
      best_estimator_ = model.fit(X = X, y = y)


    # STOCKAGE DANS SELF
    self.best_estimator_ = best_estimator_
    self.best_gamma_ = best_gamma_
    self.best_score_ = best_score_
    self.grid_scores_ = grid_scores_
    self.best_scores_ = best_scores_
    self.best_params_ = best_params_
    self.cv_list_ = cv_list_



  def _cvOverGamma(self, gamma, sample_weight, weightFit, weightEval, kwargs,
      X, y):
    """
    Performs a gridsearch with given weight sample_weight and exponent gamma
    """
    # Pour l'apprentissage, on prendra soin de ne pas mettre de poids du tout
    # si gamma == 0. En
    # effet, sklearn ne renvoie pas le même résultat avec des poids égaux à
    # 1, et s'il n'y a pas de poids. La facon de proceder choisir ici permet
    # donc d'avoir une coherence entre gamma = 0 et l'absence de poids.

    # paramètres de poids à envoyer au fit et au scoring
    fit_params = {"sample_weight" : sample_weight**gamma} if (weightFit
      and gamma != 0) else {}
    scorer_params = {"sample_weight" : sample_weight} if (weightEval) else {}

    # Cross-validation
    scorerCV = GridSearchCV(
      refit = False,
      fit_params = fit_params,
      scorer_params = scorer_params,
      **kwargs)
    scorerCV.fit(X = X, y = y)

    # Retour
    return scorerCV
