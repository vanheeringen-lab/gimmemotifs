# Copyright (c) 2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
""" Module for motif activity prediction """


def warn(*args, **kwargs):
    pass


import warnings

warnings.warn = warn
warnings.filterwarnings("ignore", message="sklearn.externals.joblib is deprecated")

import os
import sys

try:
    from itertools import izip
except ImportError:
    izip = zip
import logging

import pandas as pd
import numpy as np
from scipy.stats import hypergeom, mannwhitneyu
from statsmodels.sandbox.stats.multicomp import multipletests
from tqdm.auto import tqdm

# scikit-learn
from sklearn.ensemble import RandomForestClassifier
from sklearn.multiclass import OneVsRestClassifier
from sklearn.linear_model import MultiTaskLassoCV, BayesianRidge
from sklearn.multioutput import MultiOutputRegressor
from sklearn.preprocessing import scale, LabelEncoder
from sklearn.svm import LinearSVR
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

import xgboost

from gimmemotifs import __version__
from gimmemotifs.motif import read_motifs
from gimmemotifs.scanner import scan_regionfile_to_table
from gimmemotifs.config import MotifConfig
from gimmemotifs.utils import pfmfile_location

import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

logger = logging.getLogger("gimme.maelstrom")


class Moap(object):
    """Moap base class.

    Motif activity prediction.
    """

    _predictors = {}
    name = None

    @classmethod
    def create(cls, name, ncpus=None):
        """Create a Moap instance based on the predictor name.

        Parameters
        ----------
        name : str
            Name of the predictor (eg. Xgboost, BayesianRidge, ...)

        ncpus : int, optional
            Number of threads. Default is the number specified in the config.

        Returns
        -------
        moap : Moap instance
            moap instance.
        """
        try:
            return cls._predictors[name.lower()](ncpus=ncpus)
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


@register_predictor("BayesianRidge")
class BayesianRidgeMoap(Moap):
    def __init__(self, scale=True, ncpus=None):
        """Predict motif activities using Bayesian Ridge Regression.

        Parameters
        ----------
        scale : boolean, optional, default True
            If ``True``, the motif scores will be scaled
            before classification.

        ncpus : int, optional
            Number of threads. Default is the number specified in the config.

        Attributes
        ----------
        act_ : DataFrame, shape (n_motifs, n_clusters)
            Coefficients of the regression model.
        """

        self.act_description = "activity values: coefficients of the" "regression model"

        if ncpus is None:
            ncpus = int(MotifConfig().get_default_params().get("ncpus", 2))
        self.ncpus = ncpus
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
            df_X[:] = scale(df_X, axis=0)

        # logger.debug("Scaling y")

        # Normalize across samples and features
        # y = df_y.apply(scale, 1).apply(scale, 0)
        y = df_y

        X = df_X.loc[y.index]

        model = BayesianRidge()
        logger.debug("Fitting model")
        coefs = []
        for col in tqdm(y.columns, total=len(y.columns)):
            model.fit(X, y[col])
            coefs.append(model.coef_)
        logger.info("Done")

        self.act_ = pd.DataFrame(coefs, columns=X.columns, index=y.columns).T


@register_predictor("Xgboost")
class XgboostRegressionMoap(Moap):
    def __init__(self, scale=True, ncpus=None):
        """Predict motif activities using XGBoost.

        Parameters
        ----------
        scale : boolean, optional, default True
            If ``True``, the motif scores will be scaled
            before classification

        ncpus : int, optional
            Number of threads. Default is the number specified in the config.

        Attributes
        ----------
        act_ : DataFrame, shape (n_motifs, n_clusters)
            Feature scores.
        """

        self.act_description = "activity values: feature scores from" "fitted model"

        if ncpus is None:
            ncpus = int(MotifConfig().get_default_params().get("ncpus", 2))
        self.ncpus = ncpus
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
            df_X[:] = scale(df_X, axis=0)

        # Normalize across samples and features
        # y = df_y.apply(scale, 1).apply(scale, 0)
        y = df_y
        X = df_X.loc[y.index]

        # Define model
        xgb = xgboost.XGBRegressor(
            n_estimators=500,
            learning_rate=0.01,
            nthread=self.ncpus,
            min_child_weight=2,
            max_depth=3,
            subsample=0.8,
            colsample_bytree=0.8,
            objective="reg:squarederror",
        )

        logger.debug("xgb: 0%")

        self.act_ = pd.DataFrame(index=X.columns)

        # Fit model
        for i, col in enumerate(tqdm(y.columns)):
            xgb.fit(X, y[col].values)
            d = xgb.get_booster().get_fscore()
            self.act_[col] = [d.get(m, 0) for m in X.columns]

            for motif in self.act_.index:
                if self.act_.loc[motif, col] != 0:
                    high = df_y.loc[
                        df_X[motif] >= df_X[motif].quantile(0.75), col
                    ].mean()
                    low = df_y.loc[
                        df_X[motif] <= df_X[motif].quantile(0.25), col
                    ].mean()
                    if low > high:
                        self.act_.loc[motif, col] *= -1

            logger.debug("..{}%".format(int(float(i + 1) / len(y.columns) * 100)))
        logger.info("Done")


@register_predictor("MWU")
class MWUMoap(Moap):
    def __init__(self, *args, **kwargs):
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
        self.act_description = (
            "activity values: BH-corrected " "-log10 Mann-Whitney U p-value"
        )
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
        clusters = df_y[df_y.columns[0]].unique()
        for cluster in clusters:
            pos = df_X[df_y.iloc[:, 0] == cluster]
            neg = df_X[df_y.iloc[:, 0] != cluster]
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
        fpr = multipletests(pvals.flatten(), method="fdr_bh")[1].reshape(pvals.shape)

        # create output DataFrame
        self.act_ = pd.DataFrame(-np.log10(fpr.T), columns=clusters, index=df_X.columns)

        logger.info("Done")


@register_predictor("Hypergeom")
class HypergeomMoap(Moap):
    def __init__(self, *args, **kwargs):
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
        self.act_description = (
            "activity values: -log10-transformed, BH-corrected "
            "hypergeometric p-values"
        )
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
            pos = df_X[df_y.iloc[:, 0] == cluster]
            neg = df_X[df_y.iloc[:, 0] != cluster]

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
        fpr = multipletests(pvals.flatten(), method="fdr_bh")[1].reshape(pvals.shape)

        # create output DataFrame
        self.act_ = pd.DataFrame(-np.log10(fpr.T), columns=clusters, index=df_X.columns)

        logger.info("Done")


@register_predictor("RF")
class RFMoap(Moap):
    def __init__(self, ncpus=None):
        """Predict motif activities using a random forest classifier

        Parameters
        ----------
        ncpus : int, optional
            Number of threads. Default is the number specified in the config.

        Attributes
        ----------
        act_ : DataFrame, shape (n_motifs, n_clusters)
            feature importances from the model

        """
        self.act_ = None
        if ncpus is None:
            ncpus = int(MotifConfig().get_default_params().get("ncpus", 2))
        self.ncpus = ncpus
        self.act_description = (
            "activity values: feature importances " "from fitted Random Forest model"
        )
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
        y = le.fit_transform(df_y.iloc[:, 0].values)

        clf = RandomForestClassifier(n_estimators=100, n_jobs=self.ncpus)

        # Multiclass
        if len(le.classes_) > 2:
            orc = OneVsRestClassifier(clf)
            orc.fit(df_X.values, y)

            importances = np.array([c.feature_importances_ for c in orc.estimators_]).T
        else:  # Only two classes
            clf.fit(df_X.values, y)
            importances = np.array(
                [clf.feature_importances_, clf.feature_importances_]
            ).T

        for i, _ in enumerate(le.classes_):
            diff = df_X.loc[y == i].quantile(q=0.75) - df_X.loc[y != i].quantile(q=0.75)
            sign = (diff >= 0) * 2 - 1
            importances[:, i] *= sign

        # create output DataFrame
        self.act_ = pd.DataFrame(
            importances,
            columns=le.inverse_transform(range(len(le.classes_))),
            index=df_X.columns,
        )
        logger.info("Done")


@register_predictor("MultiTaskLasso")
class MultiTaskLassoMoap(Moap):
    def __init__(self, scale=True, ncpus=None):
        """Predict motif activities using MultiTaskLasso.

        Parameters
        ----------
        scale : boolean, optional, default True
            If ``True``, the motif scores will be scaled
            before classification.

        ncpus : int, optional
            Number of threads. Default is the number specified in the config.

        Attributes
        ----------
        act_ : DataFrame, shape (n_motifs, n_clusters)
            Coefficients of the regression model.
        """

        self.act_description = "activity values: coefficients of the" "regression model"

        if ncpus is None:
            ncpus = int(MotifConfig().get_default_params().get("ncpus", 2))
        self.ncpus = ncpus
        self.scale = scale
        self.act_ = None
        self.pref_table = "score"
        self.supported_tables = ["score", "count"]
        self.ptype = "regression"

    def fit(self, df_X, df_y):
        logger.info("Fitting MultiTaskLasso")

        if not df_y.shape[0] == df_X.shape[0]:
            raise ValueError("number of regions is not equal")

        if self.scale:
            logger.debug("Scaling motif scores")
            # Scale motif scores
            df_X.loc[:, :] = scale(df_X, axis=0)

        # logger.debug("Scaling y")

        # Normalize across samples and features
        # y = df_y.apply(scale, 1).apply(scale, 0)
        y = df_y

        X = df_X.loc[y.index]

        model = Pipeline(
            [
                ("scale", StandardScaler()),
                (
                    "reg",
                    MultiTaskLassoCV(
                        fit_intercept=False, n_alphas=20, n_jobs=self.ncpus
                    ),
                ),
            ]
        )
        logger.debug("Fitting model")
        model.fit(df_X, df_y)
        logger.info("Done")

        self.act_ = pd.DataFrame(
            model.steps[1][1].coef_, index=y.columns, columns=X.columns
        ).T

    def predict(self, df_X):
        return df_X.dot(self.act_.loc[df_X.columns])


@register_predictor("SVR")
class SVRMoap(Moap):
    def __init__(self, scale=True, ncpus=None):
        """Predict motif activities using Support Vector Regression.

        Parameters
        ----------
        scale : boolean, optional, default True
            If ``True``, the motif scores will be scaled
            before classification.

        ncpus : int, optional
            Number of threads. Default is the number specified in the config.

        Attributes
        ----------
        act_ : DataFrame, shape (n_motifs, n_clusters)
            SVR weights.
        """

        self.act_description = "activity values: SVR weights"

        if ncpus is None:
            ncpus = int(MotifConfig().get_default_params().get("ncpus", 2))
        self.ncpus = ncpus
        self.scale = scale
        self.act_ = None
        self.pref_table = "score"
        self.supported_tables = ["score", "count"]
        self.ptype = "regression"

    def fit(self, df_X, df_y):
        logger.info("Fitting SVR")

        if not df_y.shape[0] == df_X.shape[0]:
            raise ValueError("number of regions is not equal")

        if self.scale:
            logger.debug("Scaling motif scores")
            # Scale motif scores
            df_X.loc[:, :] = scale(df_X, axis=0)

        # logger.debug("Scaling y")

        # Normalize across samples and features
        # y = df_y.apply(scale, 1).apply(scale, 0)
        y = df_y
        self.columns = df_y.columns
        X = df_X.loc[y.index]

        clf = LinearSVR()
        self.model = MultiOutputRegressor(clf, n_jobs=1)
        logger.debug("Fitting model")
        self.model.fit(df_X, df_y)
        logger.info("Done")

        self.act_ = pd.DataFrame(
            {c: e.coef_ for c, e in zip(df_y.columns, self.model.estimators_)},
            index=X.columns,
        )

    def predict(self, df_X):
        # print(self.model.predict(df_X) )

        return pd.DataFrame(
            self.model.predict(df_X), index=df_X.index, columns=self.columns
        )


def moap(
    inputfile,
    method="hypergeom",
    scoring=None,
    outfile=None,
    motiffile=None,
    pfmfile=None,
    genome=None,
    fpr=0.01,
    ncpus=None,
    subsample=None,
    zscore=True,
    gc=True,
):
    """Run a single motif activity prediction algorithm.

    Parameters
    ----------
    inputfile : str
        :1File with regions (chr:start-end) in first column and either cluster
        name in second column or a table with values.

    method : str, optional
        Motif activity method to use. Any of 'hypergeom', 'lasso',
        'bayesianridge',
        'rf', 'xgboost'. Default is 'hypergeom'.

    scoring:  str, optional
        Either 'score' or 'count'

    outfile : str, optional
        Name of outputfile to save the fitted activity values.

    motiffile : str, optional
        Table with motif scan results. First column should be exactly the same
        regions as in the inputfile.

    pfmfile : str, optional
        File with motifs in pwm format. Required when motiffile is not
        supplied.

    genome : str, optional
        Genome name, as indexed by gimme. Required when motiffile is not
        supplied

    fpr : float, optional
        FPR for motif scanning

    ncpus : int, optional
        Number of threads to use. Default is the number specified in the config.

    zscore : bool, optional
        Use z-score normalized motif scores.

    gc : bool optional
        Use GC% bins for z-score.

    Returns
    -------
    pandas DataFrame with motif activity
    """

    if scoring and scoring not in ["score", "count"]:
        raise ValueError("valid values are 'score' and 'count'")

    if inputfile.endswith("feather"):
        df = pd.read_feather(inputfile)
        df = df.set_index(df.columns[0])
    else:
        # read data
        df = pd.read_table(inputfile, index_col=0, comment="#")

    clf = Moap.create(method, ncpus=ncpus)

    if clf.ptype == "classification":
        if df.shape[1] != 1:
            raise ValueError("1 column expected for {}".format(method))
    else:
        if np.dtype("object") in set(df.dtypes):
            raise ValueError("columns should all be numeric for {}".format(method))

    if motiffile is None:
        if genome is None:
            raise ValueError("need a genome")

        pfmfile = pfmfile_location(pfmfile)
        try:
            motifs = read_motifs(pfmfile)
        except Exception:
            sys.stderr.write("can't read motifs from {}".format(pfmfile))
            raise

        # scan for motifs
        motif_names = [m.id for m in read_motifs(pfmfile)]
        scores = []
        if method == "classic" or scoring == "count":
            logger.info("motif scanning (scores)")
            scores = scan_regionfile_to_table(
                inputfile,
                genome,
                "count",
                pfmfile=pfmfile,
                ncpus=ncpus,
                zscore=zscore,
                gc=gc,
            )
        else:
            logger.info("motif scanning (scores)")
            scores = scan_regionfile_to_table(
                inputfile,
                genome,
                "score",
                pfmfile=pfmfile,
                ncpus=ncpus,
                zscore=zscore,
                gc=gc,
            )
        motifs = pd.DataFrame(scores, index=df.index, columns=motif_names)

    elif isinstance(motiffile, pd.DataFrame):
        motifs = motiffile
    else:
        motifs = pd.read_table(motiffile, index_col=0, comment="#")

    if outfile and os.path.exists(outfile):
        out = pd.read_table(outfile, index_col=0, comment="#")
        ncols = df.shape[1]
        if ncols == 1:
            ncols = len(df.iloc[:, 0].unique())

        if out.shape[0] == motifs.shape[1] and out.shape[1] == ncols:
            logger.warn("%s output already exists... skipping", method)
            return out

    if subsample is not None:
        n = int(subsample * df.shape[0])
        logger.debug("Subsampling %d regions", n)
        df = df.sample(n)

    motifs = motifs.loc[df.index]

    clf.fit(motifs, df)

    if outfile:
        with open(outfile, "w") as f:
            f.write("# maelstrom - GimmeMotifs version {}\n".format(__version__))
            f.write("# method: {} with motif {}\n".format(method, scoring))
            if genome:
                f.write("# genome: {}\n".format(genome))
            if isinstance(motiffile, str):
                f.write("# motif table: {}\n".format(motiffile))
            f.write("# {}\n".format(clf.act_description))

        with open(outfile, "a") as f:
            clf.act_.to_csv(f, sep="\t")

    return clf.act_
