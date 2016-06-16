# Copyright (c) 2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
""" Module for motif activity prediction """

import os
import sys
import argparse
from functools import partial

from multiprocessing import Pool

import pandas as pd 
import numpy as np
from scipy.stats import scoreatpercentile, ks_2samp, hypergeom,mannwhitneyu
from scipy.cluster.hierarchy import linkage, fcluster
from statsmodels.sandbox.stats.multicomp import multipletests

# scikit-learn
from sklearn.cluster import AgglomerativeClustering
from sklearn.cross_validation import train_test_split
from sklearn.grid_search import GridSearchCV
from sklearn.ensemble import BaggingClassifier,RandomForestClassifier
from sklearn.multiclass import OneVsRestClassifier
from sklearn.linear_model import Ridge,MultiTaskLasso
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.preprocessing import scale, LabelEncoder
from sklearn.multiclass import OneVsRestClassifier

from lightning.classification import CDClassifier

import pymc as pm

from gimmemotifs.motif import read_motifs
from gimmemotifs.scanner import Scanner
from gimmemotifs.mara import make_model
from gimmemotifs.config import MotifConfig, GM_VERSION

CLUSTER_METHODS = ["classic", "ks", "lightning", "rf", "mwu"]
VALUE_METHODS = ["lasso", "mara"]

class LightningMoap(object):
    def __init__(self, scale=True):
        """Predict motif activities using lighting CDClassifier

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
        
        self.act_ = None
        self.sig_ = None

    def fit(self, df_X, df_y):
        
        if not df_y.shape[0] == df_X.shape[0]:
            raise ValueError("number of regions is not equal")
        if df_y.shape[1] != 1:
            raise ValueError("y needs to have 1 label column")
        
        if self.scale:
            # Scale motif scores
            df_X = df_X.apply(lambda x: scale(x))
        
        idx = range(df_y.shape[0]) 

        y = df_y.iloc[idx]
        X = df_X.loc[y.index].values
        y = y.values.flatten()
            
        # Convert (putative) string labels 
        l = LabelEncoder()
        y = l.fit_transform(y)

        # Split data 
        X_train,X_test,y_train,y_test = train_test_split(X,y)

        sys.stderr.write("Setting parameters through cross-validation\n")
        # Determine best parameters based on CV
        self.clf.fit(X_train,y_train)
    
        sys.stdout.write("Average score ({} fold CV): {}\n".format(
                self.kfolds,
                self.clf.score(X_test, y_test)
                ))

        sys.stderr.write("Estimate coefficients using bootstrapping\n")

        # Estimate coefficients using bootstrappig
        #b = BaggingClassifier(self.clf.best_estimator_, 
        #        max_samples=0.75, n_jobs=-1, random_state=state)
        b = BaggingClassifier(self.clf.best_estimator_, 
                max_samples=0.75, n_jobs=-1)
        b.fit(X,y)
        
        # Get mean coefficients
        coeffs = np.array([e.coef_ for e in b.estimators_]).mean(axis=0)
        
        # Create dataframe of predicted coefficients 
        self.act_ = pd.DataFrame(coeffs.T)
        
        # Convert labels back to original names
        self.act_.columns = l.inverse_transform(range(len(l.classes_)))
        self.act_.index = df_X.columns
        
        # Permutations
        sys.stderr.write("Permutations\n")
        random_dfs = []
        for i in range(10):
            y_random = np.random.permutation(y)
            b.fit(X,y_random)
            coeffs = np.array([e.coef_ for e in b.estimators_]).mean(axis=0)
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

class MWMoap(object):
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
    
    def fit(self, df_X, df_y):
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
        fdr = multipletests(pvals.flatten(), 
                method="fdr_bh")[1].reshape(pvals.shape)
        
        # create output DataFrame
        self.act_ = pd.DataFrame(-np.log10(pvals.T), 
                columns=clusters, index=df_X.columns)

class KSMoap(object):
    def __init__(self):
        """Predict motif activities using Kolmogorov-Smirnov p-value
    
        This method compares the motif score distribution of each 
        cluster versus the motif score distribution of all other 
        clusters.
        
        Parameters
        ----------
       
        Attributes
        ----------
        act_ : DataFrame, shape (n_motifs, n_clusters)
            -log10 of the KS p-value, corrected for multiple
            testing using the Benjamini-Hochberg correction
        """
        self.act_ = None
        self.act_description = ("activity values: BH-corrected "
                               "-log10 KS p-value")
    
    def fit(self, df_X, df_y):
        if not df_y.shape[0] == df_X.shape[0]:
            raise ValueError("number of regions is not equal")
        if df_y.shape[1] != 1:
            raise ValueError("y needs to have 1 label column")
        
        # calculate Kolmogorov-Smirnov p-values
        pvals = []
        clusters  =  df_y[df_y.columns[0]].unique()
        for cluster in clusters:
            pos = df_X[df_y.iloc[:,0] == cluster]
            neg = df_X[df_y.iloc[:,0] != cluster]
            p = [ks_2samp(pos[m], neg[m])[1] for m in pos.columns]
            pvals.append(p)
        
        # correct for multipe testing
        pvals = np.array(pvals)
        fdr = multipletests(pvals.flatten(), 
                method="fdr_bh")[1].reshape(pvals.shape)
        
        # create output DataFrame
        self.act_ = pd.DataFrame(-np.log10(pvals.T), 
                columns=clusters, index=df_X.columns)

class ClassicMoap(object):
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
    
    def fit(self, df_X, df_y):
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
        fdr = multipletests(pvals.flatten(), 
                method="fdr_bh")[1].reshape(pvals.shape)
        
        # create output DataFrame
        self.act_ = pd.DataFrame(-np.log10(pvals.T), 
                columns=clusters, index=df_X.columns)

class RFMoap(object):
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
        self.act_description = ("activity values: feature importances "
                               "from fitted Random Forest model")

    def fit(self, df_X, df_y):
        if not df_y.shape[0] == df_X.shape[0]:
            raise ValueError("number of regions is not equal")
        if df_y.shape[1] != 1:
            raise ValueError("y needs to have 1 label column")

        le = LabelEncoder()
        y = le.fit_transform(df_y.iloc[:,0].values)

        clf = RandomForestClassifier(n_estimators=100)
        orc = OneVsRestClassifier(clf)
        orc.fit(df_X.values, y)

        test = np.array([c.feature_importances_ for c in orc.estimators_]).T

        # create output DataFrame
        self.act_ = pd.DataFrame(test,
                columns=le.inverse_transform(range(len(le.classes_))),
                index=df_X.columns)

class MaraMoap(object):
    def __init__(self, iterations=10000):
        """Predict motif activities using a MARA-like algorithm

        Bayesian regression model, fitted using MCMC

        Parameters
        ----------
        iterations : int, optional, default 100000
            number of iterations used in MCMC

        Attributes
        ----------
        act_ : DataFrame, shape (n_motifs, n_clusters)
            fitted motif activities
        zscore_ : DataFrame, shape (n_motifs, n_clusters)
            z-score (trace mean / trace std)
        act_trace_mean_ : DataFrame, shape (n_motifs, n_clusters)
            trace of the mean of the activities from the MCMC
        act_trace_std_ : DataFrame, shape (n_motifs, n_clusters)
            trace of the standard deviation of the activities from the MCMC
        """
        self.act_description = ("activity values: mean coefficient from "
                               "fitted model")
        
        self.iterations = iterations

        # initialize attributes
        self.ridge_ = None
        self.act_ = None 
        self.act_trace_mean_ = None
        self.act_trace_std_ = None
        self.zscore_ = None

    def fit(self, df_X, df_y):
        if not df_y.shape[0] == df_X.shape[0]:
            raise ValueError("number of regions is not equal")
        
        X = df_X.copy()
        y = df_y.copy()

        # subtract row- and column-wise mean
        y = y - y.mean(0)
        y = y.subtract(y.mean(1), 0)

        # subtract column-wise mean
        X = X.subtract(X.mean(1), 0)
        
        # ridge regression
        sys.stderr.write("ridge regression\n")
        clf = Ridge(alpha = 0.1, solver="svd")
        ridge_result = []
        for col in y.columns:
            clf.fit(X.values, y[col].values)
            ridge_result.append(clf.coef_.transpose())
        
        self.ridge_ = pd.DataFrame(ridge_result, 
                columns=df_X.columns, index=df_y.columns).T
        
        # standardize the coefficients
        self.ridge_ = self.ridge_.subtract(self.ridge_.mean(1), 0)
        
        # initialize empty DataFrames for result
        self.act_ = pd.DataFrame(index=df_X.columns)
        self.act_trace_mean_ = pd.DataFrame(index=df_X.columns)
        self.act_trace_std_ = pd.DataFrame(index=df_X.columns)
        self.zscore_ = pd.DataFrame(index=df_X.columns)
 
       
        # fit the activities
        for col in df_y.columns:
            self._fit_model(X, y, col)

    def _fit_model(self, df_X, df_y, col):
        sys.stderr.write("fitting {}\n".format(col))

        model = make_model(
                len(df_X.columns),
                self.ridge_[col].values,
                df_y[col].values,
                df_X.values)

        map_ = pm.MAP(model)
        map_.fit()
        mcmc = pm.MCMC(model)
        mcmc.sample( self.iterations , int(self.iterations * 0.1), 10)
        zscores = np.mean(mcmc.trace("a")[:,:], 0) / np.std(mcmc.trace("a")[:,:], 0)

        self.act_[col] = mcmc.a.value
        self.act_trace_mean_[col] = np.mean(mcmc.trace("a")[:,:], 0)
        self.act_trace_std_[col] = np.std(mcmc.trace("a")[:,:], 0)
        self.zscore_[col] = zscores
        
        mcmc.db.close()

class LassoMoap(object):
    def __init__(self, scale=True, kfolds=5, alpha_stepsize=1/3.0):
        """Predict motif activities using Lasso MultiTask regression

        Parameters
        ----------
        scale : boolean, optional, default True
            If ``True``, the motif scores will be scaled 
            before classification
 
        kfolds : integer, optional, default 5
            number of kfolds for parameter search
        
        alpha_stepsize : float, optional, default 0.333
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

        # initialize attributes
        self.act_ = None 
        self.sig_ = None 
    
        mtk = MultiTaskLasso()
        parameters = {
            "alpha": [np.exp(-x) for x in np.arange(0, 10, alpha_stepsize)],
        }
        self.clf = GridSearchCV(mtk, parameters, cv=kfolds, n_jobs=4)

    def fit(self, df_X, df_y):
        if not df_y.shape[0] == df_X.shape[0]:
            raise ValueError("number of regions is not equal")
        
        idx = range(df_y.shape[0])
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
        for i in range(10):
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

    
    def _get_coefs(self, X, y):
        n_samples = 0.75 * X.shape[0]
        max_samples = X.shape[0]
        m = self.clf.best_estimator_
        coefs = []
        for i in range(10):
            idx = np.random.randint(0, n_samples, max_samples)
            m.fit(X[idx], y[idx])
            coefs.append(m.coef_)
        coefs = np.array(coefs).mean(axis=0)
        return coefs

NCPUS=4

def fit_model(data):
    X,y,multi,alpha, C = data
    #print "fitting {} {}".format(X.shape, y.shape)
    # Set classifier options.
    clf = CDClassifier(penalty="l1/l2",
                       loss="squared_hinge",
                       multiclass=multi,
                       max_iter=20,
                       alpha=alpha,
                       C=C,
                       tol=1e-3)

    # Train the model.
    return clf.fit(X, y)

def eval_model(df, sets, motifs, alpha, nsample=1000, k=10, cutoff=0):
    ret = select_sets(df, sets)
    y = pd.DataFrame({"label":0}, index=df.index)
    for label, rows in enumerate(ret):
        y.loc[rows] = label + 1
    y = y[y["label"] > 0]
    y -= 1

    clf = CDClassifier(penalty="l1/l2",
                       loss="squared_hinge",
                       multiclass=len(sets) > 2,
                       max_iter=20,
                       alpha=alpha,
                       C=1.0 / motifs.shape[0],
                       tol=1e-3)

    accs = []
    fractions = []

    for i in np.arange(k):

        idx = np.random.choice(range(y.shape[0]), nsample, replace=True)

        y_pred = y.iloc[idx[:nsample * 0.8 + 1]]
        X_pred = motifs.loc[y_pred.index].values
        y_pred = y_pred.values.flatten()

        y_test = y.iloc[idx[nsample * 0.8 + 1:]]
        X_test = motifs.loc[y_test.index].values
        y_test = y_test.values.flatten()

        # train the model
        clf.fit(X_pred, y_pred)

        acc = clf.score(X_test, y_test)
        fraction = clf.n_nonzero(percentage=True)

        accs.append(acc)
        fractions.append(fraction)

    #print alpha, accs, fractions
    return alpha, np.median(accs), np.median(fractions)

def select_sets(df, sets, threshold=1):
    abs_diff = 1 
    ret = []
    for s in sets:
        other = [c for c in df.columns if not c in s]
        x = df[(df[s] >= threshold).any(1) & ((df[s].max(1) - df[other].max(1)) >= abs_diff)].index
        ret.append(x)
    return ret

class MoreMoap(object):
    def __init__(self, scale=True):
        """Predict motif activities using Lasso MultiTask regression
        In development, doesn't really work consistently.
        """
        self.scale = scale

    def fit(self, df_X, df_y):
        
        y = df_y.apply(scale, 0)
        dist = pairwise_distances(y.values.T)
        L = linkage(dist, method="ward")
        
        self.dfs = dict(
                [(exp, pd.DataFrame(index=df_X.columns)) for exp in df_y.columns])

        for nclus in range(3, len(df_y.columns) + 1):
            labels = fcluster(L, nclus, 'maxclust')
            print labels
            if max(labels) < nclus:
                sys.stderr.write("remaining clusters are too similar")
                break
            sets = []
            for i in range(1, nclus + 1):
                sets.append(list(df_y.columns[labels == i]))
            print sets
            nbootstrap = 100
            result = self._run_bootstrap_model(
                        y, sets, df_X, nsample=1000, nbootstrap=nbootstrap
                        )
            if result is None:
                sys.stderr.write("no result for {} clusters".format(nclus))
                continue 

            if result.shape[1] < (nbootstrap * nclus) * 0.5:
                sys.stderr.write("not enough for bootstrap\n")
                break
            
            for i in range(nclus):
                print '***'
                print sets
                print result.columns[:len(sets)]
                print '***'
                cols = result.columns[range(i,result.shape[1], nclus)]
                print cols
                act = result[cols].mean(1)
                for col in df_y.columns[labels == i + 1]:
                    print col
                    self.dfs[col][str(nclus)] = act.values

            print self.dfs["NK"].shape
            print self.dfs["NK"].columns 

        self.act_ = pd.DataFrame(index=df_X.columns)
        for col in self.dfs.keys():
            self.act_[col] = self.dfs[col].mean(1)


    def _select_alpha(self, df, sets, motifs, nsample=1000):
        f = partial(eval_model, df, sets, motifs,
            nsample=int(nsample * 1.25), k=5)
    
        pool = Pool(NCPUS)
        alphas = [np.exp(-x) for x in np.arange(0, 10, 1/3.0)]
        results = pool.map(f, alphas)
        pool.close()
        pool.join()
    
        results = np.array(results)
        alpha, acc, fraction = results[np.argmax(results[:,1])]
    
        sys.stderr.write("alpha {}  accuracy: {}  fraction: {}\n".format(alpha, acc, fraction))
        return alpha
    

    def _run_bootstrap_model(self, df, sets, motifs, nsample=1000, cutoff=0, nbootstrap=100):

        ret = select_sets(df, sets)
        y = pd.DataFrame({"label":0}, index=df.index)
        for label, rows in enumerate(ret):
            y.loc[rows] = label + 1
        y = y[y["label"] > 0]
        
        if y.shape[0] == 0:
            sys.stderr.write("no sets with these filters\n")
            return
        
        if nsample > (len(y) / 2):
            nsample = len(y) / 2
            sys.stderr.write("setting nsample to {}\n".format(nsample))
        else:
            sys.stderr.write("nsample = {}\n".format(nsample))
    
        coef = pd.DataFrame(index=motifs.columns)
    
        sys.stderr.write("Selecting alpha\n")
        alpha = self._select_alpha(df, sets, motifs, nsample=nsample)
    
        acc  = []
        fraction = []
    
        data = []
    
        sys.stderr.write("Running {} bootstraps\n".format(nbootstrap))
        for i in range(nbootstrap):
    
            # Sample with replacement
            idx = np.random.choice(range(len(y.index)), nsample)
    
            # Create test dataset
            y_small = y.iloc[idx]
            if len(np.unique(y_small)) != len(sets):
                continue
            
            X_small = motifs.loc[y_small.index]
    
            y_small -= 1
            y_small = y_small.values.flatten()
    
            X = X_small.values
            data.append((X, y_small, len(sets) > 2, alpha, 1.0 / motifs.shape[0]))
    
        pool = Pool(NCPUS)
        
        results = pool.map(fit_model, data)
    
        #close the pool and wait for the work to finish 
        pool.close()
        pool.join()
        for i,(clf,(X,y_small,_,_,_)) in enumerate(zip(results, data)):
            # Accuracy
            #print y_small.
            acc.append(clf.score(X, y_small))
    
            # Percentage of selected features
            fraction.append(clf.n_nonzero(percentage=True))
             
            c = clf.coef_
            names = ["_".join(s) for s in sets]
            d = dict(zip(names, c))
            c = pd.DataFrame(d, index=motifs.columns).fillna(0)[names]
            coef = coef.join(c, rsuffix=i, how="left")
            #c = c[(abs(c) > cutoff).any(1)]
            #c = c.join(m2f)
    
            #counts.loc[c.index.values] += 1
            #print counts.sort("count").tail(1)
    
        print "Average accuracy", np.mean(acc)
        print "Average fraction", np.mean(fraction)
        return coef

def moap(inputfile, method="classic", scoring="score", outfile=None, motiffile=None, pwmfile=None, genome=None, cutoff=0.95):
    """ Run a single motif activity prediction algorithm.
    
    Parameters
    ----------
    
    inputfile : str
        File with regions (chr:start-end) in first column and either cluster 
        name in second column or a table with values.
    
    method : str, optional
        Motif activity method to use. Any of 'classic', 'ks', 'lasso', 
        'lightning', 'mara', 'rf'. Default is 'classic'. 
    
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
    
    cutoff : float, optional
        Cutoff for motif scanning
    
    Returns
    -------
    
    pandas DataFrame with motif activity
    """

    if scoring not in ['score', 'count']:
        raise ValueError("valid values are 'score' and 'count'")
    
    config = MotifConfig()

    m2f = None
    
    # read data
    df = pd.read_table(inputfile, index_col=0)

    if method in CLUSTER_METHODS:
        if df.shape[1] != 1:
            raise ValueError("1 column expected for {}".format(method))
    else:
        if np.dtype('object') in set(df.dtypes):
            raise ValueError(
                    "columns should all be numeric for {}".format(method))
        if method not in VALUE_METHODS:
            raise ValueError("method {} not valid".format(method))

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
            for row in s.count(list(df.index), cutoff=cutoff):
                scores.append(row)
        else:
            for row in s.best_score(list(df.index)):
                scores.append(row)

        motifs = pd.DataFrame(scores, index=df.index, columns=motif_names)
    else:
        motifs = pd.read_table(motiffile, index_col=0)   

    clf = None
    if method == "ks":
        clf = KSMoap()
    if method == "mwu":
        clf = MWMoap()
    if method == "rf":
        clf = RFMoap()
    if method == "lasso":
        clf = LassoMoap()
    if method == "lightning":
        clf = LightningMoap()
    if method == "mara":
        clf = MaraMoap()
    if method == "more":
        clf = MoreMoap()
    if method == "classic":
        clf = ClassicMoap()

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
