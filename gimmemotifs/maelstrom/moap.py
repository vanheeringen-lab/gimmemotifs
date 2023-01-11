# Copyright (c) 2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
""" Module for motif activity prediction """
import logging
import os

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, mannwhitneyu
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import BayesianRidge, MultiTaskLassoCV
from sklearn.multiclass import OneVsRestClassifier
from sklearn.multioutput import MultiOutputRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, StandardScaler, scale
from sklearn.svm import LinearSVR
from statsmodels.stats.multitest import multipletests
from tqdm.auto import tqdm

from gimmemotifs import __version__
from gimmemotifs.config import MotifConfig
from gimmemotifs.motif import read_motifs
from gimmemotifs.scanner import scan_regionfile_to_table
from gimmemotifs.utils import pfmfile_location

try:
    import xgboost  # noqa: optional

    _has_xgboost = True
except ImportError:
    _has_xgboost = False

logger = logging.getLogger("gimme.maelstrom")


class Moap(object):
    """Moap base class.

    Motif activity prediction.
    """

    _predictors = {}
    name = None

    @classmethod
    def create(cls, name, ncpus=None, **kwargs):
        """Create a Moap instance based on the predictor name.

        Parameters
        ----------
        name : str
            Name of the predictor (eg. Xgboost, BayesianRidge, ...)

        ncpus : int, optional
            Number of threads. Default is the number specified in the config.

        kwargs : any, optional
            keyword arguments passed to predictor.
            Unused arguments are dropped.

        Returns
        -------
        moap : Moap instance
            moap instance.
        """
        try:
            obj = cls._predictors[name.lower()]
        except KeyError:
            raise Exception("Unknown class")

        # filter for kwargs used by predictors
        accepted_kwargs = obj.__init__.__code__.co_varnames
        for k in list(kwargs):
            if k not in accepted_kwargs:
                del kwargs[k]

        return obj(ncpus=ncpus, **kwargs)

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
    def list_predictors(cls):
        """List available predictors."""
        return list(cls._predictors.keys())

    @classmethod
    def list_classification_predictors(cls):
        """List available classification predictors."""
        preds = cls._predictors.values()
        return [x.name for x in preds if x.ptype == "classification"]

    @classmethod
    def list_regression_predictors(cls):
        """List available regression predictors."""
        preds = cls._predictors.values()
        return [x.name for x in preds if x.ptype == "regression"]


register_predictor = Moap.register_predictor


@register_predictor("BayesianRidge")
class BayesianRidgeMoap(Moap):
    act_ = None
    act_description = "activity values: coefficients of the regression model"
    pref_table = "score"
    supported_tables = ["score", "count"]
    ptype = "regression"

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
        if ncpus is None:
            ncpus = int(MotifConfig().get_default_params().get("ncpus", 2))
        self.ncpus = ncpus
        self.scale = scale

    def fit(self, df_X, df_y):
        logger.info("Fitting BayesianRidge")

        if not df_y.shape[0] == df_X.shape[0]:
            raise ValueError("number of regions is not equal")

        if self.scale:
            logger.debug("Scaling motif scores")
            # Scale motif scores
            df_X[:] = scale(df_X, axis=0)

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
    act_ = None
    act_description = "activity values: feature scores from fitted model"
    pref_table = "score"
    supported_tables = ["score", "count"]
    ptype = "regression"

    def __init__(self, scale=True, ncpus=None, random_state=None):
        """Predict motif activities using XGBoost.

        Parameters
        ----------
        scale : boolean, optional, default True
            If ``True``, the motif scores will be scaled
            before classification

        ncpus : int, optional
            Number of threads. Default is the number specified in the config.

        random_state : numpy.random.RandomState object, optional
            make predictions deterministic.

        Attributes
        ----------
        act_ : DataFrame, shape (n_motifs, n_clusters)
            Feature scores.
        """
        if _has_xgboost is False:
            raise ImportError("Optional dependency xgboost is required.")

        if ncpus is None:
            ncpus = int(MotifConfig().get_default_params().get("ncpus", 2))
        self.ncpus = ncpus
        self.scale = scale
        self.random_state = random_state

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
            # xgb docs: nondeterministic with booster=gblinear & shotgun updater
            random_state=self.random_state,
        )

        self.act_ = pd.DataFrame(index=X.columns)

        # Fit model
        for col in tqdm(y.columns):
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

        logger.info("Done")


@register_predictor("MWU")
class MWUMoap(Moap):
    act_ = None
    act_description = "activity values: BH-corrected -log10 Mann-Whitney U p-value"
    pref_table = "score"
    supported_tables = ["score"]
    ptype = "classification"

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
        pass

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
                    logger.error(str(e))
                    logger.error(f"motif {m} failed, setting to p = 1")
                    p.append(1)
            pvals.append(p)

        # correct for multiple testing
        pvals = np.array(pvals)
        fpr = multipletests(pvals.flatten(), method="fdr_bh")[1].reshape(pvals.shape)

        # create output DataFrame
        self.act_ = pd.DataFrame(-np.log10(fpr.T), columns=clusters, index=df_X.columns)

        logger.info("Done")


@register_predictor("Hypergeom")
class HypergeomMoap(Moap):
    act_ = None
    act_description = (
        "activity values: -log10-transformed, BH-corrected hypergeometric p-values"
    )
    pref_table = "count"
    supported_tables = ["count"]
    ptype = "classification"

    def __init__(self, random_state=None, *args, **kwargs):
        """Predict motif activities using hypergeometric p-value

        Parameters
        ----------

        Attributes
        ----------
        act_ : DataFrame, shape (n_motifs, n_clusters)
            -log10 of the hypergeometric p-value, corrected for multiple
            testing using the Benjamini-Hochberg correction
        """
        self.random_state = random_state

    def fit(self, df_X, df_y):
        logger.info("Fitting Hypergeom")
        if not df_y.shape[0] == df_X.shape[0]:
            raise ValueError("number of regions is not equal")
        if df_y.shape[1] != 1:
            raise ValueError("y needs to have 1 label column")

        if set(df_X.dtypes) != {np.dtype(int)}:
            raise ValueError("need motif counts, not scores")

        # calculate hypergeometric p-values
        pvals = []
        clusters = df_y[df_y.columns[0]].unique()
        M = df_X.shape[0]
        hypergeometric = hypergeom
        # not sure if these have an effect (but it doesn't hurt)
        hypergeometric.seed = self.random_state
        hypergeometric.random_state = self.random_state
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
                p.append(hypergeometric.sf(x, M, n, N))

            pvals.append(p)

        # correct for multiple testing
        pvals = np.array(pvals)
        fpr = multipletests(pvals.flatten(), method="fdr_bh")[1].reshape(pvals.shape)

        # create output DataFrame
        self.act_ = pd.DataFrame(-np.log10(fpr.T), columns=clusters, index=df_X.columns)

        logger.info("Done")


@register_predictor("RF")
class RFMoap(Moap):
    act_ = None
    act_description = (
        "activity values: feature importances from fitted Random Forest model"
    )
    pref_table = "score"
    supported_tables = ["score", "count"]
    ptype = "classification"

    def __init__(self, ncpus=None, random_state=None):
        """Predict motif activities using a random forest classifier

        Parameters
        ----------
        ncpus : int, optional
            Number of threads. Default is the number specified in the config.

        random_state : numpy.random.RandomState object, optional
            make predictions deterministic.

        Attributes
        ----------
        act_ : DataFrame, shape (n_motifs, n_clusters)
            feature importances from the model
        """
        if ncpus is None:
            ncpus = int(MotifConfig().get_default_params().get("ncpus", 2))
        self.ncpus = ncpus
        self.random_state = random_state

    def fit(self, df_X, df_y):
        logger.info("Fitting RF")
        if not df_y.shape[0] == df_X.shape[0]:
            raise ValueError("number of regions is not equal")
        if df_y.shape[1] != 1:
            raise ValueError("y needs to have 1 label column")

        le = LabelEncoder()
        y = le.fit_transform(df_y.iloc[:, 0].values)

        clf = RandomForestClassifier(
            n_estimators=100, n_jobs=self.ncpus, random_state=self.random_state
        )

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
    act_ = None
    act_description = "activity values: coefficients of the regression model"
    pref_table = "score"
    supported_tables = ["score", "count"]
    ptype = "regression"

    def __init__(self, scale=True, ncpus=None, random_state=None):
        """Predict motif activities using MultiTaskLasso.

        Parameters
        ----------
        scale : boolean, optional, default True
            If ``True``, the motif scores will be scaled
            before classification.

        ncpus : int, optional
            Number of threads. Default is the number specified in the config.

        random_state : numpy.random.RandomState object, optional
            make predictions deterministic.

        Attributes
        ----------
        act_ : DataFrame, shape (n_motifs, n_clusters)
            Coefficients of the regression model.
        """
        if ncpus is None:
            ncpus = int(MotifConfig().get_default_params().get("ncpus", 2))
        self.ncpus = ncpus
        self.scale = scale
        self.random_state = random_state

    def fit(self, df_X, df_y):
        logger.info("Fitting MultiTaskLasso")

        if not df_y.shape[0] == df_X.shape[0]:
            raise ValueError("number of regions is not equal")

        if self.scale:
            logger.debug("Scaling motif scores")
            # Scale motif scores
            df_X.loc[:, :] = scale(df_X, axis=0)

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
                        fit_intercept=False,
                        n_alphas=20,
                        n_jobs=self.ncpus,
                        random_state=self.random_state,
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
    act_ = None
    act_description = "activity values: SVR weights"
    pref_table = "score"
    supported_tables = ["score", "count"]
    ptype = "regression"

    def __init__(self, scale=True, ncpus=None, random_state=None):
        """Predict motif activities using Support Vector Regression.

        Parameters
        ----------
        scale : boolean, optional, default True
            If ``True``, the motif scores will be scaled
            before classification.

        ncpus : int, optional
            Number of threads. Default is the number specified in the config.

        random_state : numpy.random.RandomState object, optional
            make predictions deterministic.

        Attributes
        ----------
        act_ : DataFrame, shape (n_motifs, n_clusters)
            SVR weights.
        """
        if ncpus is None:
            ncpus = int(MotifConfig().get_default_params().get("ncpus", 2))
        self.ncpus = ncpus
        self.scale = scale
        self.columns = None
        self.model = None
        self.random_state = random_state

    def fit(self, df_X, df_y):
        logger.info("Fitting SVR")

        if not df_y.shape[0] == df_X.shape[0]:
            raise ValueError("number of regions is not equal")

        if self.scale:
            logger.debug("Scaling motif scores")
            # Scale motif scores
            df_X.loc[:, :] = scale(df_X, axis=0)

        # Normalize across samples and features
        # y = df_y.apply(scale, 1).apply(scale, 0)
        y = df_y
        self.columns = df_y.columns
        X = df_X.loc[y.index]

        clf = LinearSVR(random_state=self.random_state)
        self.model = MultiOutputRegressor(clf, n_jobs=1)
        logger.debug("Fitting model")
        self.model.fit(df_X, df_y)
        logger.info("Done")

        self.act_ = pd.DataFrame(
            {c: e.coef_ for c, e in zip(df_y.columns, self.model.estimators_)},
            index=X.columns,
        )

    def predict(self, df_X):
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
    zscore=True,
    gc=True,
    subsample=None,
    random_state=None,
    ncpus=None,
    progress=None,
):
    """Run a single motif activity prediction algorithm.

     Parameters
     ----------
     inputfile : str
         :1File with regions (chr:start-end) in first column and either cluster
         name in second column or a table with values.

     method : str, optional
         Motif activity method to use. Any of
         'bayesianridge', 'xgboost', 'mwu', 'hypergeom',
         'rf', 'multitasklasso', 'svr'. Default is 'hypergeom'.

     scoring:  str, optional
         Either 'score' or 'count'

     outfile : str, optional
         Name of outputfile to save the fitted activity values.

     motiffile : str, optional
         Table with motif scan results. First column should be exactly the same
         regions as in the inputfile.

     pfmfile : str, optional
         File with motifs in pfm format. Required when motiffile is not
         supplied.

     genome : str, optional
         Genome name, as indexed by gimme. Required when motiffile is not
         supplied.

     zscore : bool, optional
         Use z-score normalized motif scores.

     gc : bool, optional
         Equally distribute GC percentages in background sequences.

     subsample : float, optional
         Fraction of regions to use.

     random_state : numpy.random.RandomState object, optional
         make predictions deterministic (where possible).

     ncpus : int, optional
         Number of threads to use.
         Default is the number specified in the config.

    progress : bool or None, optional
         provide progress bars for long computations.

     Returns
     -------
     pandas DataFrame with motif activity
    """

    if scoring and scoring not in ["score", "count"]:
        raise ValueError("valid values are 'score' and 'count'")

    # read data
    if inputfile.endswith("feather"):
        df = pd.read_feather(inputfile)
        df = df.set_index(df.columns[0])
    else:
        df = pd.read_table(inputfile, index_col=0, comment="#")

    clf = Moap.create(method, ncpus=ncpus, random_state=random_state)

    if clf.ptype == "classification":
        if df.shape[1] != 1:
            raise ValueError(f"1 column expected for {method}")
    else:
        if np.dtype("object") in set(df.dtypes):
            raise ValueError(f"columns should all be numeric for {method}")

    if motiffile is None:
        if genome is None:
            raise ValueError("need a genome")

        pfmfile = pfmfile_location(pfmfile)
        try:
            _ = read_motifs(pfmfile)
        except Exception:
            logger.error(f"can't read motifs from {pfmfile}")
            raise

        # scan for motifs
        motif_names = [m.id for m in read_motifs(pfmfile)]
        if method == "classic" or scoring == "count":
            logger.info("motif scanning (counts)")
            scores = scan_regionfile_to_table(
                inputfile,
                genome,
                "count",
                pfmfile=pfmfile,
                ncpus=ncpus,
                zscore=zscore,
                gc=gc,
                random_state=random_state,
                progress=progress,
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
                random_state=random_state,
                progress=progress,
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
            logger.warning(f"{method} output already exists... skipping")
            return out

    if subsample is not None:
        n = int(subsample * df.shape[0])
        logger.debug(f"Subsampling {n} regions")
        df = df.sample(n, random_state=random_state)

    motifs = motifs.loc[df.index]

    clf.fit(motifs, df)

    if outfile:
        with open(outfile, "w") as f:
            f.write(f"# maelstrom - GimmeMotifs version {__version__}\n")
            f.write(f"# method: {method} with motif {scoring}\n")
            if genome:
                f.write(f"# genome: {genome}\n")
            if isinstance(motiffile, str):
                f.write(f"# motif table: {motiffile}\n")
            f.write(f"# {clf.act_description}\n")

        with open(outfile, "a") as f:
            clf.act_.to_csv(f, sep="\t")

    return clf.act_
