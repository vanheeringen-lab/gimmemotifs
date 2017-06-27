# Copyright (c) 2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
""" Module for motif activity prediction """
from __future__ import print_function
def warn(*args, **kwargs):
        pass
import warnings
warnings.warn = warn

import os
import sys
from functools import partial
try: 
    from itertools import izip
except:
    izip = zip
import itertools
import logging
from multiprocessing import Pool

import pandas as pd 
import numpy as np
from scipy.stats import ks_2samp, hypergeom,mannwhitneyu
from scipy.cluster.hierarchy import linkage, fcluster
from statsmodels.sandbox.stats.multicomp import multipletests

# scikit-learn
from sklearn.model_selection import train_test_split,GridSearchCV
from sklearn.ensemble import BaggingClassifier,RandomForestClassifier
from sklearn.multiclass import OneVsRestClassifier
from sklearn.linear_model import Ridge,MultiTaskLasso,BayesianRidge
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.preprocessing import scale, LabelEncoder
from lightning.classification import CDClassifier
from lightning.regression import CDRegressor

import xgboost

from gimmemotifs.motif import read_motifs
from gimmemotifs.scanner import Scanner
from gimmemotifs.config import MotifConfig, GM_VERSION

logger = logging.getLogger("gimme.maelstrom")

class Moap(object):
    """Moap base class.

    Motif activity prediction.
    """
    _predictors = {}
    name = None

    @classmethod
    def create(cls, name):
        """Create a Moap instance based on the predictor name.

        Parameters
        ----------
        name : str
            Name of the predictor (eg. Xgboost, BayesianRidge, ...)
        
        Returns
        -------
        moap : Moap instance
            moap instance.
        """
        try:
            return cls._predictors[name.lower()]()
        except KeyError:
            raise Exception("Unknown class")

    @classmethod
    def register_predictor(cls, name):
        """Register method to keep list of predictors."""
        def decorator(subclass):
            """Register as decorator function."""
            cls._predictors[name.lower()] = subclass
            subclass.name = name.lower()
            return subclass
        return decorator

    @classmethod
    def list_predictors(self):
        """List available predictors."""
        return list(self._predictors.keys())

    @classmethod
    def list_classification_predictors(self):
        """List available classification predictors."""
        preds = [self.create(x) for x in self._predictors.keys()]
        return [x.name for x in preds if x.ptype == "classification"]

    @classmethod
    def list_regression_predictors(self):
        """List available regression predictors."""
        preds = [self.create(x) for x in self._predictors.keys()]
        return [x.name for x in preds if x.ptype == "regression"]

register_predictor = Moap.register_predictor

def br_fit(X, y):
    model = BayesianRidge()
    model.fit(X, y)
    return model.coef_

def br_fit_star(args):
    return br_fit(*args)

@register_predictor('BayesianRidge')
class BayesianRidgeMoap(Moap):
    def __init__(self, scale=True):
        """Predict motif activities using Bayesian Ridge Regression.

        Parameters
        ----------
        scale : boolean, optional, default True
            If ``True``, the motif scores will be scaled 
            before classification.
       
        Attributes
        ----------
        act_ : DataFrame, shape (n_motifs, n_clusters)
            Coefficients of the regression model. 
        """
        
        self.act_description = ("activity values: coefficients of the"
                               "regression model")
        
        self.ncpus = int(MotifConfig().get_default_params().get("ncpus", 2))
        self.scale = scale
        self.act_ = None
        self.pref_table = "score"
        self.supported_tables = ["score", "count"]
        self.ptype = "regression"

    def fit(self, df_X, df_y):
        logger.info("Fitting BayesianRidge")
        
        if not df_y.shape[0] == df_X.shape[0]:
            raise ValueError("number of regions is not equal")
        
        if self.scale:
            logger.debug("Scaling motif scores")
            # Scale motif scores
            df_X = df_X.apply(scale)
        
        #logger.debug("Scaling y")
        
        # Normalize across samples and features
        #y = df_y.apply(scale, 1).apply(scale, 0)
        y = df_y

        X = df_X.loc[y.index]

        logger.debug("Fitting model")
        pool = Pool(self.ncpus) 
        coefs = pool.map(br_fit_star, izip(itertools.repeat(X), [y[col] for col in y.columns]))
        logger.info("Done")

        self.act_ = pd.DataFrame(coefs, columns=X.columns, index=y.columns).T

@register_predictor('Xgboost')
class XgboostRegressionMoap(Moap):
    def __init__(self, scale=True):
        """Predict motif activities using XGBoost.

        Parameters
        ----------
        scale : boolean, optional, default True
            If ``True``, the motif scores will be scaled 
            before classification
       
        Attributes
        ----------
        act_ : DataFrame, shape (n_motifs, n_clusters)
            Feature scores. 
        """
        
        self.act_description = ("activity values: feature scores from"
                               "fitted model")
        
        self.ncpus = int(MotifConfig().get_default_params().get("ncpus", 2))
        self.scale = scale
        
        self.act_ = None
        self.pref_table = "score"
        self.supported_tables = ["score", "count"]
        self.ptype = "regression"
    
    def fit(self, df_X, df_y):
        logger.info("Fitting XGBoostRegression")
        
        if not df_y.shape[0] == df_X.shape[0]:
            raise ValueError("number of regions is not equal")
        
        if self.scale:
            # Scale motif scores
            df_X = df_X.apply(scale)
        
        # Normalize across samples and features
        #y = df_y.apply(scale, 1).apply(scale, 0)
        y = df_y
        X = df_X.loc[y.index]

        # Define model
        xgb = xgboost.XGBRegressor(
                n_estimators=100,
                learning_rate=0.1,
                nthread=self.ncpus,
                min_child_weight=2,
                max_depth=3,
                subsample=0.75,
                colsample_bytree=0.75,
                objective='reg:linear')

        logger.debug("xgb: 0%")
        
        self.act_ = pd.DataFrame(index=X.columns)
        
        # Fit model
        for i,col in enumerate(y.columns):
            xgb.fit(X, y[col].values)
            d = xgb.booster().get_fscore()
            self.act_[col] = [d.get(m, 0) for m in X.columns]
            
            for motif in self.act_.index:
                if self.act_.loc[motif, col] != 0:
                    high = df_y.loc[df_X[motif] >= df_X[motif].quantile(0.75), col].mean()
                    low = df_y.loc[df_X[motif] <= df_X[motif].quantile(0.25), col].mean()
                    if low > high:
                        self.act_.loc[motif, col] *= -1
            
            logger.debug("..{}%".format(int(float(i + 1)/ len(y.columns) * 100)))
        logger.info("Done")

@register_predictor('LightningRegressor')
class LightningRegressionMoap(Moap):
    def __init__(self, scale=True, cv=3):
        """Predict motif activities using lightning CDRegressor 

        Parameters
        ----------
        scale : boolean, optional, default True
            If ``True``, the motif scores will be scaled 
            before classification
       
        cv : int, optional, default 3
            Cross-validation k-fold parameter.

        Attributes
        ----------
        act_ : DataFrame, shape (n_motifs, n_clusters)
            fitted coefficients

        sig_ : DataFrame, shape (n_motifs,)
            boolean values, if coefficients are higher/lower than
            the 1%t from random permutation
        """
        
        self.act_description = ("activity values: coefficients from "
                               "fitted model")
        
        self.ncpus = int(MotifConfig().get_default_params().get("ncpus", 2))
        self.kfolds = cv
        self.scale = scale
        
        self.act_ = None
        self.pref_table = "score"
        self.supported_tables = ["score", "count"]
        self.ptype = "regression"
    
    def fit(self, df_X, df_y):
        logger.info("Fitting LightningRegression")
        
        if self.scale:
            # Scale motif scores
            df_X = df_X.apply(scale)
        
        # Normalize across samples and features
        #y = df_y.apply(scale, 1).apply(scale, 0)
        y = df_y
        X = df_X.loc[y.index]

        if not y.shape[0] == X.shape[0]:
            raise ValueError("number of regions is not equal")
        
        # Define model
        cd = CDRegressor(penalty="l1/l2", C=1.0/X.shape[0])
        parameters = {
            "alpha": [np.exp(-x) for x in np.arange(0, 8, 1/2.5)],
        }
        clf = GridSearchCV(cd, parameters, n_jobs=self.ncpus)

        # Fit model
        clf.fit(X.values, y.values)
        
        # Get coefficients
        self.act_ = pd.DataFrame(clf.best_estimator_.coef_.T, 
                index=X.columns, columns = y.columns)
        
        logger.info("Done")

@register_predictor('LightningClassification')
class LightningClassificationMoap(Moap):
    def __init__(self, scale=True, permute=False):
        """Predict motif activities using lightning CDClassifier 

        Parameters
        ----------
        scale : boolean, optional, default True
            If ``True``, the motif scores will be scaled 
            before classification
       
        Attributes
        ----------
        act_ : DataFrame, shape (n_motifs, n_clusters)
            fitted coefficients

        sig_ : DataFrame, shape (n_motifs,)
            boolean values, if coefficients are higher/lower than
            the 1%t from random permutation
        """
        
        self.act_description = ("activity values: coefficients from "
                               "fitted model")
        
        #self.cdc = CDClassifier(random_state=args.seed)
        self.cdc = CDClassifier()
        
        self.parameters = {
            "penalty": ["l1/l2"],
            "loss": ["squared_hinge"],
            "multiclass":[True],
            "max_iter":[20],
            "alpha": [np.exp(-x) for x in np.arange(0, 10, 1/3.0)],
            "C":[0.001, 0.01, 0.1, 0.5, 1.0],
            "tol":[1e-3]
        }

        self.kfolds = 10
        
        self.clf = GridSearchCV(self.cdc, self.parameters, 
                cv=self.kfolds, n_jobs=-1)

        self.scale = scale
        self.permute = permute

        self.act_ = None
        self.sig_ = None
        self.pref_table = "score"
        self.supported_tables = ["score", "count"]
        self.ptype = "classification"

    def fit(self, df_X, df_y):
        logger.info("Fitting LightningClassification")
        
        if not df_y.shape[0] == df_X.shape[0]:
            raise ValueError("number of regions is not equal")
        if df_y.shape[1] != 1:
            raise ValueError("y needs to have 1 label column")
        
        if self.scale:
            # Scale motif scores
            df_X = df_X.apply(scale)
        
        idx = list(range(df_y.shape[0]))

        y = df_y.iloc[idx]
        X = df_X.loc[y.index].values
        y = y.values.flatten()
            
        # Convert (putative) string labels 
        l = LabelEncoder()
        y = l.fit_transform(y)

        # Split data 
        X_train,X_test,y_train,y_test = train_test_split(X,y)

        logger.debug("Setting parameters through cross-validation")
        # Determine best parameters based on CV
        self.clf.fit(X_train,y_train)
    
        logger.debug("Average score ({} fold CV): {}".format(
                self.kfolds,
                self.clf.score(X_test, y_test)
                ))

        logger.debug("Estimate coefficients using bootstrapping")

        # Estimate coefficients using bootstrappig
        #b = BaggingClassifier(self.clf.best_estimator_, 
        #        max_samples=0.75, n_jobs=-1, random_state=state)
        b = BaggingClassifier(self.clf.best_estimator_, 
                max_samples=0.75, n_jobs=-1)
        b.fit(X,y)
        
        # Get mean coefficients
        coeffs = np.array([e.coef_ for e in b.estimators_]).mean(axis=0)
        
        # Create dataframe of predicted coefficients 
        if len(l.classes_) == 2:
            self.act_ = pd.DataFrame(np.hstack((-coeffs.T, coeffs.T)))
        else:
            self.act_ = pd.DataFrame(coeffs.T)

        # Convert labels back to original names
        self.act_.columns = l.inverse_transform(range(len(l.classes_)))
        self.act_.index = df_X.columns
        
        if self.permute:
            # Permutations
            logger.debug("Permutations")
            random_dfs = []
            for _ in range(10):
                y_random = np.random.permutation(y)
                b.fit(X,y_random)
                coeffs = np.array([e.coef_ for e in b.estimators_]).mean(axis=0)
            
                if len(l.classes_) == 2:
                    random_dfs.append(pd.DataFrame(np.hstack((-coeffs.T, coeffs.T))))
                else:
                    random_dfs.append(pd.DataFrame(coeffs.T))
            random_df = pd.concat(random_dfs)
    
            # Select cutoff based on percentile
            high_cutoffs = random_df.quantile(0.99)
            low_cutoffs = random_df.quantile(0.01)
    
            # Set significance
            self.sig_ = pd.DataFrame(index=df_X.columns)
            self.sig_["sig"] = False

            for col,c_high,c_low in zip(
                    self.act_.columns, high_cutoffs, low_cutoffs):
                self.sig_["sig"].loc[self.act_[col] >= c_high] = True
                self.sig_["sig"].loc[self.act_[col] <= c_low] = True
        logger.info("Done")

@register_predictor('MWU')
class MWUMoap(Moap):
    def __init__(self):
        """Predict motif activities using Mann-Whitney U p-value
    
        This method compares the motif score distribution of each 
        cluster versus the motif score distribution of all other 
        clusters.
        
        Parameters
        ----------
       
        Attributes
        ----------
        act_ : DataFrame, shape (n_motifs, n_clusters)
            -log10 of the Mann-Whitney U p-value, corrected for multiple
            testing using the Benjamini-Hochberg correction
        """
        self.act_ = None
        self.act_description = ("activity values: BH-corrected "
                               "-log10 Mann-Whitney U p-value")
        self.pref_table = "score"
        self.supported_tables = ["score"]
        self.ptype = "classification"
    
    def fit(self, df_X, df_y):
        logger.info("Fitting MWU")

        if not df_y.shape[0] == df_X.shape[0]:
            raise ValueError("number of regions is not equal")
        if df_y.shape[1] != 1:
            raise ValueError("y needs to have 1 label column")
        
        # calculate Mann-Whitney U p-values
        pvals = []
        clusters  =  df_y[df_y.columns[0]].unique()
        for cluster in clusters:
            pos = df_X[df_y.iloc[:,0] == cluster]
            neg = df_X[df_y.iloc[:,0] != cluster]
            p = []
            for m in pos:
                try:
                    p.append(mannwhitneyu(pos[m], neg[m], alternative="greater")[1])
                except Exception as e:
                    sys.stderr.write(str(e) + "\n")
                    sys.stderr.write("motif {} failed, setting to p = 1\n".format(m))
                    p.append(1)
            pvals.append(p)
        
        # correct for multipe testing
        pvals = np.array(pvals)
        fpr = multipletests(pvals.flatten(), 
                method="fdr_bh")[1].reshape(pvals.shape)
        
        # create output DataFrame
        self.act_ = pd.DataFrame(-np.log10(fpr.T), 
                columns=clusters, index=df_X.columns)
        
        logger.info("Done")

@register_predictor('Hypergeom')
class HypergeomMoap(Moap):
    def __init__(self):
        """Predict motif activities using hypergeometric p-value

        Parameters
        ----------
       
        Attributes
        ----------
        act_ : DataFrame, shape (n_motifs, n_clusters)
            -log10 of the hypergeometric p-value, corrected for multiple
            testing using the Benjamini-Hochberg correction
        """
        self.act_ = None
        self.act_description = ("activity values: BH-corrected "
                                "hypergeometric p-values")
        self.pref_table = "count"
        self.supported_tables = ["count"]
        self.ptype = "classification"
    
    def fit(self, df_X, df_y):
        logger.info("Fitting Hypergeom")
        if not df_y.shape[0] == df_X.shape[0]:
            raise ValueError("number of regions is not equal")
        if df_y.shape[1] != 1:
            raise ValueError("y needs to have 1 label column")
       
        if set(df_X.dtypes) != set([np.dtype(int)]):
            raise ValueError("need motif counts, not scores")

        # calculate hypergeometric p-values
        pvals = []
        clusters = df_y[df_y.columns[0]].unique()
        M = df_X.shape[0]
        for cluster in clusters:
            pos = df_X[df_y.iloc[:,0] == cluster]
            neg = df_X[df_y.iloc[:,0] != cluster]
            
            pos_true = (pos > 0).sum(0)
            pos_false = (pos == 0).sum(0)
            neg_true = (neg > 0).sum(0)
            
            p = []
            for pt, pf, nt in zip(pos_true, pos_false, neg_true):
                n = pt + nt
                N = pt + pf
                x = pt - 1
                p.append(hypergeom.sf(x, M, n, N))
            
            pvals.append(p)
        
        # correct for multipe testing
        pvals = np.array(pvals)
        fpr = multipletests(pvals.flatten(), 
                method="fdr_bh")[1].reshape(pvals.shape)
        
        # create output DataFrame
        self.act_ = pd.DataFrame(-np.log10(fpr.T), 
                columns=clusters, index=df_X.columns)

        logger.info("Done")

@register_predictor('RF')
class RFMoap(Moap):
    def __init__(self):
        """Predict motif activities using a random forest classifier

        Parameters
        ----------
       
        Attributes
        ----------
        act_ : DataFrame, shape (n_motifs, n_clusters)
            feature importances from the model
        
        """
        self.act_ = None
        self.ncpus = int(MotifConfig().get_default_params().get("ncpus", 2))
        self.act_description = ("activity values: feature importances "
                               "from fitted Random Forest model")
        self.pref_table = "score"
        self.supported_tables = ["score", "count"]
        self.ptype = "classification"

    def fit(self, df_X, df_y):
        logger.info("Fitting RF")
        if not df_y.shape[0] == df_X.shape[0]:
            raise ValueError("number of regions is not equal")
        if df_y.shape[1] != 1:
            raise ValueError("y needs to have 1 label column")

        le = LabelEncoder()
        y = le.fit_transform(df_y.iloc[:,0].values)

        clf = RandomForestClassifier(n_estimators=100, n_jobs=self.ncpus)
        
        # Multiclass
        if len(le.classes_) > 2:
            orc = OneVsRestClassifier(clf)
            orc.fit(df_X.values, y)

            importances = np.array([c.feature_importances_ for c in orc.estimators_]).T
        else: # Only two classes
            clf.fit(df_X.values, y)
            importances = np.array([
                clf.feature_importances_,
                clf.feature_importances_
                ]).T
        
        for i,c in enumerate(le.classes_):
            diff = df_X.loc[y == i].quantile(q=0.75) - df_X.loc[y != i].quantile(q=0.75)
            sign = (diff >= 0) * 2 - 1
            importances[:,i] *= sign
        
        # create output DataFrame
        self.act_ = pd.DataFrame(importances,
                columns=le.inverse_transform(range(len(le.classes_))),
                index=df_X.columns)
        logger.info("Done")

@register_predictor('Lasso')
class LassoMoap(Moap):
    def __init__(self, scale=True, kfolds=5, alpha_stepsize=1.0):
        """Predict motif activities using Lasso MultiTask regression

        Parameters
        ----------
        scale : boolean, optional, default True
            If ``True``, the motif scores will be scaled 
            before classification
 
        kfolds : integer, optional, default 5
            number of kfolds for parameter search
        
        alpha_stepsize : float, optional, default 1.0
            stepsize for use in alpha gridsearch

        Attributes
        ----------
        act_ : DataFrame, shape (n_motifs, n_clusters)
            fitted motif activities
    
        sig_ : DataFrame, shape (n_motifs,)
            boolean values, if coefficients are higher/lower than
            the 1%t from random permutation
        """
        
        self.kfolds = kfolds
        self.act_description = ("activity values: coefficients from "
                               "fitted model")

        self.scale = scale
        self.ncpus = int(MotifConfig().get_default_params().get("ncpus", 2))

        # initialize attributes
        self.act_ = None 
        self.sig_ = None 
    
        mtk = MultiTaskLasso()
        parameters = {
            "alpha": [np.exp(-x) for x in np.arange(0, 10, alpha_stepsize)],
        }
        self.clf = GridSearchCV(mtk, parameters, cv=kfolds, n_jobs=self.ncpus)
        self.pref_table = "score"
        self.supported_tables = ["score", "count"]
        self.ptype = "regression"

    def fit(self, df_X, df_y):
        logger.info("Fitting Lasso")
        if not df_y.shape[0] == df_X.shape[0]:
            raise ValueError("number of regions is not equal")
       
        if self.scale:
            # Scale motif scores
            df_X = df_X.apply(scale)

        idx = list(range(df_y.shape[0]))
        y = df_y.iloc[idx]
        X = df_X.loc[y.index].values
        y = y.values
       
        X_train,X_test,y_train,y_test = train_test_split(X,y)
        sys.stderr.write("set alpha through cross-validation\n")
        # Determine best parameters based on CV
        self.clf.fit(X_train, y_train)
        
        sys.stdout.write("average score ({} fold CV): {}\n".format(
                    self.kfolds,
                    self.clf.score(X_test, y_test)
                    ))

        sys.stderr.write("Estimate coefficients using bootstrapping\n")

        # fit coefficients
        coefs = self._get_coefs(X, y)
        self.act_ = pd.DataFrame(coefs.T)
        
        # convert labels back to original names
        self.act_.columns = df_y.columns
        self.act_.index = df_X.columns

        # Permutations
        sys.stderr.write("permutations\n")
        random_dfs = []
        for _ in range(10):
            y_random = y[np.random.permutation(range(y.shape[0]))]
            coefs = self._get_coefs(X, y_random)
            random_dfs.append(pd.DataFrame(coefs.T))
        random_df = pd.concat(random_dfs)

        # Select cutoff based on percentile
        high_cutoffs = random_df.quantile(0.99)
        low_cutoffs = random_df.quantile(0.01)
        
        # Set significance
        self.sig_ = pd.DataFrame(index=df_X.columns)
        self.sig_["sig"] = False

        for col,c_high,c_low in zip(
                self.act_.columns, high_cutoffs, low_cutoffs):
            self.sig_["sig"].loc[self.act_[col] >= c_high] = True
            self.sig_["sig"].loc[self.act_[col] <= c_low] = True
        
        logger.info("Done")

    def _get_coefs(self, X, y):
        n_samples = 0.75 * X.shape[0]
        max_samples = X.shape[0]
        m = self.clf.best_estimator_
        coefs = []
        for _ in range(10):
            idx = np.random.randint(0, n_samples, max_samples)
            m.fit(X[idx], y[idx])
            coefs.append(m.coef_)
        coefs = np.array(coefs).mean(axis=0)
        return coefs

def moap(inputfile, method="hypergeom", scoring=None, outfile=None, motiffile=None, pwmfile=None, genome=None, fpr=0.01):
    """Run a single motif activity prediction algorithm.
    
    Parameters
    ----------
    inputfile : str
        File with regions (chr:start-end) in first column and either cluster 
        name in second column or a table with values.
    
    method : str, optional
        Motif activity method to use. Any of 'hypergeom', 'lasso', 
        'lightningclassification', 'lightningregressor', 'bayesianridge', 
        'rf', 'xgboost'. Default is 'hypergeom'. 
    
    scoring:  str, optional
        Either 'score' or 'count'
    
    outfile : str, optional
        Name of outputfile to save the fitted activity values.
    
    motiffile : str, optional
        Table with motif scan results. First column should be exactly the same
        regions as in the inputfile.
    
    pwmfile : str, optional
        File with motifs in pwm format. Required when motiffile is not 
        supplied.
    
    genome : str, optional
        Genome name, as indexed by gimme. Required when motiffile is not
        supplied
    
    fpr : float, optional
        FPR for motif scanning
    
    Returns
    -------
    pandas DataFrame with motif activity
    """

    if scoring and scoring not in ['score', 'count']:
        raise ValueError("valid values are 'score' and 'count'")
    
    config = MotifConfig()

    m2f = None
    
    # read data
    df = pd.read_table(inputfile, index_col=0)
    
    clf = Moap.create(method)

    if clf.ptype == "classification":
        if df.shape[1] != 1:
            raise ValueError("1 column expected for {}".format(method))
    else:
        if np.dtype('object') in set(df.dtypes):
            raise ValueError(
                    "columns should all be numeric for {}".format(method))

    if motiffile is None:
        if genome is None:
            raise ValueError("need a genome")
        # check pwmfile
        if pwmfile is None:
            pwmfile = config.get_default_params().get("motif_db", None)
            if pwmfile is not None:
                pwmfile = os.path.join(config.get_motif_dir(), pwmfile)
        
        if pwmfile is None:
            raise ValueError("no pwmfile given and no default database specified")

        if not os.path.exists(pwmfile):
            raise ValueError("{} does not exist".format(pwmfile))

        try:
            motifs = read_motifs(open(pwmfile))
        except:
            sys.stderr.write("can't read motifs from {}".format(pwmfile))
            raise

        base = os.path.splitext(pwmfile)[0]
        map_file = base + ".motif2factors.txt"
        if os.path.exists(map_file):
            m2f = pd.read_table(map_file, index_col=0)

        # initialize scanner
        s = Scanner()
        sys.stderr.write(pwmfile + "\n")
        s.set_motifs(pwmfile)
        s.set_genome(genome)

        # scan for motifs
        sys.stderr.write("scanning for motifs\n")
        motif_names = [m.id for m in read_motifs(open(pwmfile))]
        scores = []
        if method == 'classic' or scoring == "count":
            s.set_threshold(fpr=fpr)
            for row in s.count(list(df.index)):
                scores.append(row)
        else:
            for row in s.best_score(list(df.index)):
                scores.append(row)

        motifs = pd.DataFrame(scores, index=df.index, columns=motif_names)
    else:
        motifs = pd.read_table(motiffile, index_col=0)   

    if outfile and os.path.exists(outfile):
        out = pd.read_table(outfile, index_col=0, comment="#")
        ncols = df.shape[1]
        if ncols == 1:
            ncols = len(df.iloc[:,0].unique())
        
        if out.shape[0] == motifs.shape[1] and out.shape[1] == ncols:
            logger.warn("%s output already exists... skipping", method)
            return out
    
    motifs = motifs.loc[df.index]
    
    clf.fit(motifs, df)
    
    if outfile:
        with open(outfile, "w") as f:
            f.write("# maelstrom - GimmeMotifs version {}\n".format(GM_VERSION))
            f.write("# method: {} with motif {}\n".format(method, scoring))
            if genome:
                f.write("# genome: {}\n".format(genome))
            if motiffile:
                f.write("# motif table: {}\n".format(motiffile))
            f.write("# {}\n".format(clf.act_description))
        
        with open(outfile, "a") as f:
            clf.act_.to_csv(f, sep="\t")

    return clf.act_
