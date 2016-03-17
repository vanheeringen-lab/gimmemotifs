#!/usr/bin/env python
import sys
import argparse

import pandas as pd 
import numpy as np
from scipy.stats import scoreatpercentile,ks_2samp
from statsmodels.sandbox.stats.multicomp import multipletests

# scikit-learn
from sklearn.grid_search import GridSearchCV
from sklearn.preprocessing import scale, LabelEncoder
from sklearn.cross_validation import train_test_split
from sklearn.ensemble import BaggingClassifier
from sklearn.linear_model import Ridge

from lightning.classification import CDClassifier

import pymc as pm

from gimmemotifs.motif import read_motifs
from gimmemotifs.scanner import Scanner
from gimmemotifs.mara import make_model

class LightningClassifier(object):
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
        
        #self.cdc = CDClassifier(random_state=args.seed)
        self.cdc = CDClassifier()
        
        self.parameters = {
            "penalty": ["l1/l2"],
            "loss": ["squared_hinge"],
            "multiclass":[True],
            "max_iter":[20],
            "alpha": [np.exp(-x) for x in np.arange(0, 10, 1/3.0)],
            "C":[1.0 / motifs.shape[0], 0.5, 1.0],
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

class KSClassifier(object):
    def __init__(self):
        """Predict motif activities using lighting CDClassifier

        Parameters
        ----------
       
        Attributes
        ----------
        act_ : DataFrame, shape (n_motifs, n_clusters)
            -log10 of the KS p-value, corrected for multiple
            testing using the Benjamini-Hochberg correction
        """
        self.act_ = None
    
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

class MaraClassifier(object):
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

if __name__ == "__main__":

    motif_file = "/home/simon/prj/cis-bp/cis-bp.vertebrate.clusters.v3.0.pwm"
    map_file = "/home/simon/prj/cis-bp/cis-bp.vertebrate.clusters.v3.0.motif2factors.txt"
    # Load mapping of motif to transcription factors
    m2f = pd.read_csv(map_file, sep="\t", names=["motif","factors"], index_col=0)
    label_file = "test_rpkm_table.mm10.txt"

    outfile = "vla.out"
    #nsample = args.nsample
    #state = np.random.RandomState(seed=args.seed)

    # Read labels
    sys.stderr.write("read labels\n")
    df = pd.read_table(label_file, index_col=0)

    # initialize scanner
    s = Scanner()
    s.set_motifs(motif_file)
    s.set_genome("mm10")
    
    # scan for motifs
    sys.stderr.write("scanning for motifs\n")
    motif_names = [m.id for m in read_motifs(open(motif_file))]
    scores = []
    for row in s.best_score(list(df.index)):
        scores.append(row)
    motifs = pd.DataFrame(scores, index=df.index, columns=motif_names)   
    
    #clf = LightningClassifier()
    #clf = KSClassifier()
    clf = MaraClassifier(iterations=50)

    clf.fit(motifs, df)
   
    print clf.ridge_.sort_values("NK")[["NK"]].tail()
    print clf.act_.sort_values("NK")[["NK"]].tail()
    #print clf.act_.sort_values("trophoblast")

    #if nsample > 0:
    #    idx = np.random.choice(range(df.shape[0]), nsample, replace=False)
    #else:
    
    
    # Write output
    #df.to_csv(outfile, sep="\t")
