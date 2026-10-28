"""
Runtime test: overhead of SHAP and approximate SHAP on top of GRNBoost2.

For each of `n_targets` targets, a local copy of the GRNBoost2 regression
(sklearn GradientBoostingRegressor with signifikante's SGBM_KWARGS and the
out-of-bag early stopping monitor) is trained on one bootstrap of the data.
The trained model is then explained on each of the remaining bootstraps with
  - exact TreeSHAP          (shap.TreeExplainer(model).shap_values(X))
  - approximate SHAP        (Saabas, shap_values(X, approximate=True))
and the runtimes of training, SHAP and approximate SHAP are reported separately.

Bootstraps are either
  - kallisto quantification bootstraps (--kallisto_dir): bootstrap b is the
    samples x features TPM matrix built from bootstrap/bs{b} of every sample's
    abundance.h5. The model is trained on bs0 and explained on bs1..bs{B-1}.
  - resampled samples (--transcript_file or --synthetic): bootstrap b draws
    the samples with replacement.

Example:
    python -m alternet.runners.shap_overhead_test \
        --kallisto_dir ../GSE130036 --regulator_file data/regulator_db.tsv \
        --biomart_file data/biomart.txt --level gene \
        --results_path results/ --experiment_name shap_test

    # without real data
    python -m alternet.runners.shap_overhead_test --synthetic \
        --results_path results/ --experiment_name shap_test_synthetic
"""

import argparse
import glob
import os
import os.path as op
import time

import h5py
import numpy as np
import pandas as pd
import shap
from scipy.stats import spearmanr
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.preprocessing import StandardScaler
import xgboost as xgb
from xgboost import XGBRegressor

from alternet.data_preprocessing import standardize_dataframe


# Copied from signifikante.core so the model matches grnboost2 exactly
SGBM_KWARGS = {
    'learning_rate': 0.01,
    'n_estimators': 5000,  # can be arbitrarily large
    'max_features': 0.1,
    'subsample': 0.9
}
EARLY_STOP_WINDOW_LENGTH = 25
DEMON_SEED = 666


class EarlyStopMonitor:
    """Copy of signifikante.core.EarlyStopMonitor (stop when the mean OOB improvement over the window is < 0)."""

    def __init__(self, window_length=EARLY_STOP_WINDOW_LENGTH):
        self.window_length = window_length

    def __call__(self, current_round, regressor, _):
        if current_round >= self.window_length - 1:
            lo = max(0, current_round - self.window_length + 1)
            hi = current_round + 1
            return np.mean(regressor.oob_improvement_[lo:hi]) < 0
        return False


class EvalSetMonitor:
    """
    GRNBoost2's window rule on a dedicated evaluation set instead of the out-of-bag rows:
    stop when the mean improvement of the eval-set MSE over the last `window_length` trees is < 0.
    """

    def __init__(self, X_eval, y_eval, window_length=EARLY_STOP_WINDOW_LENGTH):
        self.X_eval = X_eval
        self.y_eval = y_eval
        self.window_length = window_length
        self.pred = None
        self.improvement = []

    def __call__(self, current_round, regressor, _):
        if self.pred is None:
            self.pred = regressor.init_.predict(self.X_eval).astype(float)
            self.loss = np.mean((self.pred - self.y_eval) ** 2)
        # add the newest tree to the running eval-set prediction
        self.pred += regressor.learning_rate * regressor.estimators_[current_round, 0].predict(self.X_eval)
        loss = np.mean((self.pred - self.y_eval) ** 2)
        self.improvement.append(self.loss - loss)
        self.loss = loss
        if current_round >= self.window_length - 1:
            return np.mean(self.improvement[-self.window_length:]) < 0
        return False


def fit_grnboost2(X, y, seed=DEMON_SEED, early_stop_window_length=EARLY_STOP_WINDOW_LENGTH, X_eval=None, y_eval=None):
    """
    Train one GRNBoost2 regression model (same estimator and early stopping as signifikante.grnboost2).
    With X_eval/y_eval, the early stopping window is evaluated on that set instead of the out-of-bag rows.
    """
    regressor = GradientBoostingRegressor(random_state=seed, **SGBM_KWARGS)
    if X_eval is None:
        monitor = EarlyStopMonitor(early_stop_window_length)
    else:
        monitor = EvalSetMonitor(X_eval, y_eval, early_stop_window_length)
    regressor.fit(X, y, monitor=monitor)
    return regressor


def grnboost2_importances(regressor):
    """GRNBoost2 importance: feature importances de-normalized by the number of trees."""
    return regressor.feature_importances_ * len(regressor.estimators_)


# XGBoost configured like GRNBoost2: depth-3 trees, lr 0.01, 90% row subsampling, 10% of features per split
# (sklearn's max_features is per split -> colsample_bynode) and no L2 penalty like sklearn. 'exact' does not
# support colsample_bynode; 'hist' with 256 bins splits exactly as long as a feature has <= 256 distinct values
XGB_GRNBOOST2_KWARGS = {
    'objective': 'reg:squarederror',
    'learning_rate': 0.01,
    'n_estimators': 5000,  # upper bound, early stopping on the test samples picks the number of trees
    'max_depth': 3,
    'subsample': 0.9,
    'colsample_bynode': 0.1,
    'reg_lambda': 0.0,
    'tree_method': 'hist',
    'max_bin': 256,
    'n_jobs': 1,
}


class GRNBoost2Model:
    """GRNBoost2 copy; early stopping on out-of-bag rows (or an eval set), SHAP via the shap package."""

    def __init__(self, early_stopping_rounds=EARLY_STOP_WINDOW_LENGTH):
        self.early_stopping_rounds = early_stopping_rounds

    def fit(self, X, y, X_eval=None, y_eval=None, seed=DEMON_SEED):
        if X_eval is not None and len(X_eval) == 0:
            X_eval = y_eval = None
        self.model = fit_grnboost2(X, y, seed=seed, early_stop_window_length=self.early_stopping_rounds,
                                   X_eval=X_eval, y_eval=y_eval)
        self.n_trees = len(self.model.estimators_)
        return self

    def init_explainer(self):
        self.explainer = shap.TreeExplainer(self.model)

    def shap_values(self, X):
        return self.explainer.shap_values(X, check_additivity=False)

    def approx_shap_values(self, X):
        return self.explainer.shap_values(X, approximate=True, check_additivity=False)

    def expected_value(self):
        return np.ravel(self.explainer.expected_value)[0]

    def predict(self, X):
        return self.model.predict(X)

    def importances(self):
        return grnboost2_importances(self.model)


class XGBoostModel:
    """XGBoost with early stopping on a dedicated eval set; SHAP via xgboost's built-in pred_contribs."""

    def __init__(self, early_stopping_rounds=EARLY_STOP_WINDOW_LENGTH):
        self.early_stopping_rounds = early_stopping_rounds

    def fit(self, X, y, X_eval=None, y_eval=None, seed=DEMON_SEED):
        if X_eval is None or len(X_eval) == 0:
            raise ValueError("XGBoost needs an eval set for early stopping (--n_eval_bootstraps or --test_fraction).")
        self.model = XGBRegressor(random_state=seed, early_stopping_rounds=self.early_stopping_rounds,
                                  eval_metric='rmse', **XGB_GRNBOOST2_KWARGS)
        self.model.fit(X, y, eval_set=[(X_eval, y_eval)], verbose=False)
        self.n_trees = self.model.best_iteration + 1
        return self

    def init_explainer(self):
        # the booster keeps early_stopping_rounds trees past the best iteration; drop them
        self.booster = self.model.get_booster()[:self.n_trees]

    def _contribs(self, X, approx):
        # last column is the bias (expected value)
        return self.booster.predict(xgb.DMatrix(X), pred_contribs=True, approx_contribs=approx)

    def shap_values(self, X):
        c = self._contribs(X, approx=False)
        self._bias = c[0, -1]
        return c[:, :-1]

    def approx_shap_values(self, X):
        return self._contribs(X, approx=True)[:, :-1]

    def expected_value(self):
        return self._bias

    def predict(self, X):
        return self.booster.predict(xgb.DMatrix(X))

    def importances(self):
        # total gain is the counterpart of GRNBoost2's de-normalized impurity importance
        scores = self.booster.get_score(importance_type='total_gain')
        imp = np.zeros(self.model.n_features_in_)
        for k, v in scores.items():
            imp[int(k[1:])] = v
        return imp


MODELS = {'grnboost2': GRNBoost2Model, 'xgboost': XGBoostModel}


def mse(y_pred, y):
    return np.mean((y_pred - y) ** 2) if len(y) else np.nan


def run_target(regulators, target, train, explain, seed, model_type='grnboost2', test_samples=(),
               early_stopping_rounds=EARLY_STOP_WINDOW_LENGTH, eval_set=None):
    """
    Train on `train`, explain with SHAP and approximate SHAP on each matrix in `explain`.

    Parameters:
        regulators (list): regulator columns.
        target (str): target column.
        train (pd.DataFrame): samples x features matrix of the training bootstrap(s); the index holds sample names.
        explain (list of pd.DataFrame): samples x features matrices of the remaining bootstraps.
        model_type (str): 'grnboost2' or 'xgboost'.
        test_samples (collection): samples held out of training, with all their bootstrap rows.
            Without eval_set, XGBoost uses them for early stopping; for both models they give the held-out MSE.
        eval_set (pd.DataFrame): stacked bootstraps used only for early stopping (both models);
            rows of test samples are left out of it.

    Returns:
        timing (dict): runtimes in seconds, model metadata and prediction errors.
        edges (pd.DataFrame): one row per regulator with importance and mean |SHAP| statistics.
        boot (pd.DataFrame): one row per used edge, mean |SHAP| in each explained bootstrap (explain0, explain1, ...).
    """
    # like arboreto/signifikante, the target is never its own regulator
    regs = [r for r in regulators if r != target]
    is_test = train.index.isin(list(test_samples))
    fit_rows, test_rows = train[~is_test], train[is_test]

    if eval_set is not None:
        eval_rows = eval_set[~eval_set.index.isin(list(test_samples))]
    elif model_type == 'xgboost':
        eval_rows = test_rows
    else:
        eval_rows = None  # GRNBoost2 default: out-of-bag early stopping

    model = MODELS[model_type](early_stopping_rounds)
    start = time.perf_counter()
    model.fit(fit_rows[regs].to_numpy(), fit_rows[target].to_numpy(),
              None if eval_rows is None else eval_rows[regs].to_numpy(),
              None if eval_rows is None else eval_rows[target].to_numpy(), seed=seed)
    t_train = time.perf_counter() - start

    start = time.perf_counter()
    model.init_explainer()
    t_init = time.perf_counter() - start

    t_shap = 0.0
    t_approx = 0.0
    shap_abs = []
    approx_abs = []
    for df in explain:
        Xb = df[regs].to_numpy()

        start = time.perf_counter()
        sv = model.shap_values(Xb)
        t_shap += time.perf_counter() - start

        start = time.perf_counter()
        sv_approx = model.approx_shap_values(Xb)
        t_approx += time.perf_counter() - start

        shap_abs.append(np.abs(sv).mean(axis=0))
        approx_abs.append(np.abs(sv_approx).mean(axis=0))

    # prediction error outside the timed region
    def error(df):
        if len(df) == 0:
            return np.nan
        return mse(model.predict(df[regs].to_numpy()), df[target].to_numpy())

    y_fit = fit_rows[target].to_numpy()
    mse_train = error(fit_rows)
    mse_test = error(test_rows)  # held-out samples, training bootstrap(s)
    mse_eval = error(eval_rows) if eval_set is not None else np.nan  # early-stopping bootstraps
    explain_is_test = [df.index.isin(list(test_samples)) for df in explain]
    mse_explain = np.array([error(df[~t]) for df, t in zip(explain, explain_is_test)])
    mse_explain_test = np.array([error(df[t]) for df, t in zip(explain, explain_is_test)])
    # references: predicting the training mean (no model), and how far the target itself moves between bootstraps
    mse_null = np.mean([mse(y_fit.mean(), df[target].to_numpy()) for df in explain])
    mse_null_test = mse(y_fit.mean(), test_rows[target].to_numpy())  # no-model reference for mse_test
    ref = train.iloc[:len(explain[0])]  # first training bootstrap
    target_shift = np.mean([mse(ref[target].to_numpy(), df[target].to_numpy()) for df in explain]) \
        if all(df.index.equals(ref.index) for df in explain) else np.nan

    # sanity check outside the timed region: exact SHAP values sum to the prediction
    additivity_err = np.max(np.abs(sv.sum(axis=1) + model.expected_value() - model.predict(Xb)))

    shap_abs = np.vstack(shap_abs)
    approx_abs = np.vstack(approx_abs)
    importance = model.importances()

    edges = pd.DataFrame({
        'source': regs,
        'target': target,
        'importance': importance,
        'shap_mean_abs': shap_abs.mean(axis=0),
        'shap_std_abs': shap_abs.std(axis=0),
        'approx_shap_mean_abs': approx_abs.mean(axis=0),
        'approx_shap_std_abs': approx_abs.std(axis=0),
    })

    # rank agreement only over regulators the model uses; the many unused ones are tied at 0 and inflate the correlation
    used = edges[edges['importance'] > 0]

    # per-bootstrap mean |SHAP| of the used edges (unused regulators are in no tree, so their SHAP is exactly 0)
    is_used = (importance > 0)
    boot = pd.DataFrame(shap_abs[:, is_used].T, columns=[f'explain{b}' for b in range(len(explain))])
    boot.insert(0, 'source', np.asarray(regs)[is_used])
    boot.insert(1, 'target', target)

    n_explain = len(explain)
    timing = {
        'target': target,
        'n_regulators': len(regs),
        'model': model_type,
        'n_train_rows': len(fit_rows),
        'n_test_rows': len(test_rows),
        'n_trees': model.n_trees,
        'n_explain_bootstraps': n_explain,
        't_train': t_train,
        't_explainer_init': t_init,
        't_shap': t_shap,
        't_approx_shap': t_approx,
        # overhead relative to training; the explainer is built once and is needed by both methods
        'overhead_shap': (t_init + t_shap) / t_train,
        'overhead_approx_shap': (t_init + t_approx) / t_train,
        't_shap_per_bootstrap': t_shap / n_explain,
        't_approx_shap_per_bootstrap': t_approx / n_explain,
        'n_used_regulators': len(used),
        'spearman_shap_vs_approx': spearmanr(used['shap_mean_abs'], used['approx_shap_mean_abs']).statistic,
        'spearman_shap_vs_importance': spearmanr(used['shap_mean_abs'], used['importance']).statistic,
        'max_additivity_error': additivity_err,
        'mse_train': mse_train,
        'mse_explain_mean': mse_explain.mean(),
        'mse_explain_min': mse_explain.min(),
        'mse_explain_max': mse_explain.max(),
        'mse_diff': mse_explain.mean() - mse_train,
        'mse_eval': mse_eval,
        'mse_test': mse_test,
        'mse_null_test': mse_null_test,
        'mse_explain_test_mean': mse_explain_test.mean() if len(test_rows) else np.nan,
        'mse_null_explain': mse_null,
        'mse_target_shift': target_shift,
    }
    return timing, edges, boot


def edge_agreement(edges, boot, max_log2_rank_ratio=1.0):
    """
    Per-edge agreement between GRNBoost2 importance and SHAP, for the edges the model uses.

    Within each target, edges are ranked (1 = strongest) by importance on the training bootstrap,
    by mean |SHAP| over the explained bootstraps and by mean |approximate SHAP|. The spread of the
    SHAP rank across the explained bootstraps measures stability under quantification noise.
    An edge has low agreement when |log2(rank_shap / rank_importance)| > max_log2_rank_ratio.
    """
    used = edges[edges['importance'] > 0].copy()
    by_target = used.groupby('target')
    used['n_used'] = by_target['source'].transform('size')
    used['rank_importance'] = by_target['importance'].rank(ascending=False, method='min')
    used['rank_shap'] = by_target['shap_mean_abs'].rank(ascending=False, method='min')
    used['rank_approx_shap'] = by_target['approx_shap_mean_abs'].rank(ascending=False, method='min')

    bs_cols = [c for c in boot.columns if c.startswith('explain')]
    ranks = boot.groupby('target')[bs_cols].rank(ascending=False, method='min')
    rank_stats = pd.DataFrame({
        'source': boot['source'],
        'target': boot['target'],
        'rank_shap_q05': ranks.quantile(0.05, axis=1),
        'rank_shap_median': ranks.median(axis=1),
        'rank_shap_q95': ranks.quantile(0.95, axis=1),
    })
    used = used.merge(rank_stats, on=['source', 'target'], how='left')

    used['log2_rank_ratio'] = np.log2(used['rank_shap'] / used['rank_importance'])
    used['low_agreement'] = used['log2_rank_ratio'].abs() > max_log2_rank_ratio
    return used.sort_values(['target', 'rank_importance']).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Bootstraps by resampling samples (TSV / synthetic input)
# ---------------------------------------------------------------------------

def resample_bootstraps(data, n_bootstraps, rng):
    """`n_bootstraps` copies of `data` with the samples drawn with replacement."""
    n = data.shape[0]
    return [data.iloc[rng.integers(0, n, size=n)] for _ in range(n_bootstraps)]


def load_transcript_data(transcript_file, regulator_file, level):
    """
    Standardized samples x features matrix plus regulator list, as in infer_networks_v3.

    level='gene':       TF/TF_SF genes -> genes          (canonical network)
    level='transcript': TF/SF/TF_SF transcripts -> transcripts (fully AS-aware network)
    """
    regulator_list = pd.read_csv(regulator_file, sep='\t')
    transcript_data = pd.read_csv(transcript_file, sep='\t', index_col=0)
    sample_cols = [c for c in transcript_data.columns if c not in ['transcript_id', 'gene_id']]

    if level == 'gene':
        matrix = transcript_data.groupby('gene_id')[sample_cols].sum().T
    else:
        matrix = transcript_data.set_index('transcript_id')[sample_cols].T

    # constant features cannot be standardized and carry no signal
    matrix = matrix.loc[:, matrix.std() > 0]
    data = standardize_dataframe(matrix)
    regulators = sorted(regulator_ids(regulator_list, level) & set(data.columns))
    return data, regulators


def synthetic_data(n_samples, n_regulators, n_genes, rng):
    """Random data where each gene depends on a few regulators, so the trees have structure."""
    R = rng.standard_normal((n_samples, n_regulators))
    genes = np.empty((n_samples, n_genes))
    for j in range(n_genes):
        idx = rng.choice(n_regulators, size=3, replace=False)
        w = rng.standard_normal(3)
        genes[:, j] = R[:, idx] @ w + 0.5 * np.tanh(R[:, idx[0]] * R[:, idx[1]]) + 0.5 * rng.standard_normal(n_samples)
    regulators = [f'REG{i}' for i in range(n_regulators)]
    data = pd.DataFrame(np.hstack([R, genes]), columns=regulators + [f'GENE{j}' for j in range(n_genes)])
    return standardize_dataframe(data), regulators


# ---------------------------------------------------------------------------
# Kallisto quantification bootstraps
# ---------------------------------------------------------------------------

def regulator_ids(regulator_list, level):
    """TF/TF_SF genes for the gene level, all regulator transcripts for the transcript level (as in infer_networks_v3)."""
    if level == 'gene':
        return set(regulator_list[regulator_list['Regulator_type'].isin(['TF', 'TF_SF'])]['Gene stable ID'])
    return set(regulator_list['Transcript stable ID'])


def kallisto_files(kallisto_dir):
    """Map sample name -> abundance.h5 for <kallisto_dir>/<sample>/**/abundance.h5 (or a single abundance.h5)."""
    if op.isfile(kallisto_dir):
        return {op.basename(op.dirname(kallisto_dir)): kallisto_dir}
    files = sorted(glob.glob(op.join(kallisto_dir, '*', '**', 'abundance.h5'), recursive=True))
    return {op.relpath(f, kallisto_dir).split(os.sep)[0]: f for f in files}


def counts_to_tpm(counts, eff_lengths):
    """TPM from estimated counts; works on a vector or on rows of a (bootstraps x transcripts) matrix."""
    rho = np.divide(counts, eff_lengths, out=np.zeros_like(counts, dtype=float), where=eff_lengths > 0)
    return rho / rho.sum(axis=-1, keepdims=True) * 1e6


def feature_index(transcript_ids, tx2gene, level):
    """Feature name for every kallisto transcript (version suffix stripped; gene ID at gene level)."""
    tx = pd.Index([t.split('.')[0] for t in transcript_ids])
    if level == 'gene':
        return pd.Index(tx.map(tx2gene))  # unmapped transcripts become NaN and are dropped later
    return tx


def read_kallisto_point_estimates(files):
    """Samples x transcripts TPM matrix from est_counts (version suffixes stripped), used to pick features."""
    rows = {}
    for sample, path in files.items():
        with h5py.File(path, 'r') as f:
            ids = [i.decode() for i in f['aux/ids'][:]]
            tpm = counts_to_tpm(f['est_counts'][:], f['aux/eff_lengths'][:])
        rows[sample] = pd.Series(tpm, index=feature_index(ids, None, 'transcript'))
    return pd.DataFrame(rows).T


def to_gene_level(transcript_matrix, tx2gene):
    """Sum transcript columns into gene columns; transcripts without a gene are dropped."""
    genes = transcript_matrix.columns.map(tx2gene)
    keep = genes.notna()
    return transcript_matrix.loc[:, keep].T.groupby(genes[keep]).sum().T


def supported_isoforms(annotation, point_tx, tx2gene, min_tpm, min_usage):
    """
    Supported isoforms per gene: protein-coding transcripts of protein-coding genes with
    TSL 1 or 2, in GENCODE basic, not a fragment, mean TPM >= min_tpm and mean usage >= min_usage.

    Returns:
        pd.DataFrame: one row per supported isoform (gene, transcript, names, TSL, mean TPM, usage).
    """
    ann = annotation[
        (annotation['Gene type'] == 'protein_coding')
        & (annotation['Transcript type'] == 'protein_coding')
        & annotation['Transcript support level (TSL)'].isin(['tsl1', 'tsl2'])
        & (annotation['GENCODE basic annotation'] == 'GENCODE basic')
        & (annotation['Fragment'].astype(str) != 'True')
    ]
    ann = ann[ann['Transcript stable ID'].isin(point_tx.columns)].copy()

    gene_mean = to_gene_level(point_tx, tx2gene).mean()
    ann['mean_tpm'] = point_tx[ann['Transcript stable ID']].mean().to_numpy()
    ann['usage'] = ann['mean_tpm'] / ann['Gene stable ID'].map(gene_mean).to_numpy()
    ann = ann[(ann['mean_tpm'] >= min_tpm) & (ann['usage'] >= min_usage)]
    return ann[['Gene stable ID', 'Gene name', 'Transcript stable ID', 'Transcript name',
                'Transcript support level (TSL)', 'mean_tpm', 'usage']]


def select_multi_isoform_genes(isoforms, exclude_genes, n_genes, rng):
    """Randomly pick `n_genes` genes with >= 2 supported isoforms, not in `exclude_genes`."""
    counts = isoforms.groupby('Gene stable ID').size()
    candidates = sorted(set(counts[counts >= 2].index) - set(exclude_genes))
    if len(candidates) < n_genes:
        raise ValueError(f"Only {len(candidates)} genes with >= 2 supported isoforms found.")
    genes = [str(g) for g in rng.choice(candidates, size=n_genes, replace=False)]
    return isoforms[isoforms['Gene stable ID'].isin(genes)].sort_values(['Gene stable ID', 'usage'],
                                                                        ascending=[True, False])


def read_kallisto_bootstraps(files, tx2gene, level, features, n_bootstraps):
    """
    One samples x features TPM matrix per kallisto bootstrap (bs0 .. bs{n_bootstraps-1}).

    Only `features` are kept, so memory stays small even with 200k transcripts.
    """
    features = pd.Index(features)
    data = np.empty((n_bootstraps, len(files), len(features)))
    for s, (sample, path) in enumerate(files.items()):
        with h5py.File(path, 'r') as f:
            ids = [i.decode() for i in f['aux/ids'][:]]
            eff_lengths = f['aux/eff_lengths'][:]
            available = len(f['bootstrap'])
            if available < n_bootstraps:
                raise ValueError(f"{path} has only {available} bootstraps, {n_bootstraps} requested.")
            boot = np.vstack([f[f'bootstrap/bs{b}'][:] for b in range(n_bootstraps)])
        tpm = counts_to_tpm(boot, eff_lengths)
        names = feature_index(ids, tx2gene, level)
        # sum transcripts into features (identity at transcript level) and keep only requested features
        col = features.get_indexer(names)
        keep = col >= 0
        out = np.zeros((n_bootstraps, len(features)))
        np.add.at(out.T, col[keep], tpm[:, keep].T)
        data[:, s, :] = out
    samples = list(files.keys())
    return [pd.DataFrame(data[b], index=samples, columns=features) for b in range(n_bootstraps)]


def standardize_bootstraps(train, explain):
    """Fit the scaler on the training bootstrap and apply the same transform to the others."""
    train = train.loc[:, train.std() > 0]
    scaler = StandardScaler().fit(train)

    def transform(df):
        return pd.DataFrame(scaler.transform(df[train.columns]), index=df.index, columns=train.columns)

    return transform(train), [transform(df) for df in explain]


def main():
    parser = argparse.ArgumentParser(description="Measure the runtime overhead of SHAP and approximate SHAP for GRNBoost2.")
    parser.add_argument("--kallisto_dir", type=str,
                        help="Folder with <sample>/**/abundance.h5 files; uses the kallisto bootstraps.")
    parser.add_argument("--biomart_file", type=str, default="data/biomart.txt", help="Transcript -> gene map (kallisto input).")
    parser.add_argument("--min_tpm", type=float, default=1.0,
                        help="Kallisto input: keep features with mean TPM >= min_tpm in the point estimates.")
    parser.add_argument("--multi_isoform_genes", type=int,
                        help="Kallisto input: pick this many non-regulator genes with >= 2 supported isoforms as targets "
                             "(the genes at --level gene, their supported isoforms at --level transcript).")
    parser.add_argument("--annotation_file", type=str, default="data/annotation_ensembl111.tsv",
                        help="Ensembl transcript annotation with TSL, GENCODE basic and fragment columns.")
    parser.add_argument("--min_usage", type=float, default=0.1,
                        help="Supported isoforms must have mean usage >= min_usage within their gene.")
    parser.add_argument("--max_log2_rank_ratio", type=float, default=1.0,
                        help="Edges whose SHAP rank differs from their importance rank by more than this factor "
                             "(log2; 1 = rank doubles or halves) are flagged as low agreement.")
    parser.add_argument("--transcript_file", type=str, help="Transcript TSV (transcript_id, gene_id, one column per sample).")
    parser.add_argument("--regulator_file", type=str, default="data/regulator_db.tsv")
    parser.add_argument("--level", choices=['gene', 'transcript'], default='gene',
                        help="gene: TF genes -> genes (canonical). transcript: regulator transcripts -> transcripts (fully_as).")
    parser.add_argument("--synthetic", action='store_true', help="Use synthetic data instead of --transcript_file.")
    parser.add_argument("--n_samples", type=int, default=300, help="Synthetic data only.")
    parser.add_argument("--n_regulators", type=int, default=1500, help="Synthetic data only.")
    parser.add_argument("--n_genes", type=int, default=50, help="Synthetic data only.")
    parser.add_argument("--n_targets", type=int, default=10)
    parser.add_argument("--targets", type=str, nargs='*', help="Explicit target IDs (overrides --n_targets).")
    parser.add_argument("--n_bootstraps", type=int, default=10,
                        help="Total bootstraps: n_train_bootstraps for training + the rest for SHAP.")
    parser.add_argument("--n_eval_bootstraps", type=int, default=0,
                        help="Bootstraps after the training ones, stacked into the early-stopping eval set "
                             "(both models). SHAP is computed on the bootstraps after those.")
    parser.add_argument("--n_train_bootstraps", type=int, default=1,
                        help="Bootstraps stacked as pseudoreplicates into the training matrix.")
    parser.add_argument("--model", choices=sorted(MODELS), default='grnboost2',
                        help="grnboost2: sklearn copy, early stopping on out-of-bag rows. "
                             "xgboost: same tree settings, early stopping on held-out test samples.")
    parser.add_argument("--test_fraction", type=float, default=0.0,
                        help="Fraction of samples held out of training with all their bootstrap rows "
                             "(required > 0 for --model xgboost; gives held-out MSE for both models).")
    parser.add_argument("--early_stopping_rounds", type=int, default=EARLY_STOP_WINDOW_LENGTH,
                        help="XGBoost: stop when the test RMSE has not improved for this many trees.")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--results_path", type=str, required=True)
    parser.add_argument("--experiment_name", type=str, required=True)
    args = parser.parse_args()

    k = args.n_train_bootstraps
    e = args.n_eval_bootstraps
    if k < 1 or e < 0 or args.n_bootstraps <= k + e:
        parser.error("Need --n_train_bootstraps >= 1 and --n_bootstraps > n_train + n_eval (at least one to explain).")
    n_explain = args.n_bootstraps - k - e
    if args.model == 'xgboost' and args.test_fraction <= 0 and e == 0:
        parser.error("--model xgboost needs an eval set: --n_eval_bootstraps > 0 or --test_fraction > 0.")

    rng = np.random.default_rng(args.seed)

    def pick_targets(columns, regulators):
        if args.targets:
            return args.targets
        candidates = [c for c in columns if c not in set(regulators)] or list(columns)
        return [str(t) for t in rng.choice(candidates, size=min(args.n_targets, len(candidates)), replace=False)]

    if args.kallisto_dir:
        files = kallisto_files(args.kallisto_dir)
        if not files:
            parser.error(f"No abundance.h5 found under {args.kallisto_dir}.")
        biomart = pd.read_csv(args.biomart_file, sep='\t')
        tx2gene = dict(zip(biomart['Transcript stable ID'], biomart['Gene stable ID']))
        regulator_list = pd.read_csv(args.regulator_file, sep='\t')

        point_tx = read_kallisto_point_estimates(files)
        point = to_gene_level(point_tx, tx2gene) if args.level == 'gene' else point_tx
        expressed = point.columns[point.mean() >= args.min_tpm]
        regulators = sorted(regulator_ids(regulator_list, args.level) & set(expressed))

        if args.multi_isoform_genes:
            # same seed -> same genes at both levels; targets are the genes or their supported isoforms
            annotation = pd.read_csv(args.annotation_file, sep='\t', low_memory=False)
            isoforms = supported_isoforms(annotation, point_tx, tx2gene, args.min_tpm, args.min_usage)
            selection = select_multi_isoform_genes(isoforms, set(regulator_list['Gene stable ID']),
                                                   args.multi_isoform_genes, rng)
            print("Selected genes and supported isoforms:")
            print(selection.to_string(index=False, float_format='%.3f'))
            os.makedirs(args.results_path, exist_ok=True)
            selection.to_csv(op.join(args.results_path, f"{args.experiment_name}.selection.tsv"), sep='\t', index=False)
            col = 'Gene stable ID' if args.level == 'gene' else 'Transcript stable ID'
            targets = list(dict.fromkeys(selection[col]))
        else:
            targets = pick_targets(expressed, regulators)

        boots = read_kallisto_bootstraps(files, tx2gene, args.level, sorted(set(regulators) | set(targets)),
                                         args.n_bootstraps)
        train, rest = standardize_bootstraps(pd.concat(boots[:k]), boots[k:])
        eval_set, explain = (pd.concat(rest[:e]) if e else None), rest[e:]
        # features that are constant in the training bootstrap are dropped by the scaler
        regulators = [r for r in regulators if r in train.columns]
        targets = [t for t in targets if t in train.columns]
        source = f"kallisto bootstraps ({len(files)} samples, {point.shape[1]} {args.level}s, {len(expressed)} expressed)"
    else:
        if args.synthetic:
            data, regulators = synthetic_data(args.n_samples, args.n_regulators, args.n_genes, rng)
        elif args.transcript_file:
            data, regulators = load_transcript_data(args.transcript_file, args.regulator_file, args.level)
        else:
            parser.error("Provide --kallisto_dir, --transcript_file or --synthetic.")
        targets = pick_targets(data.columns, regulators)
        boots = resample_bootstraps(data, args.n_bootstraps, rng)
        train, rest = pd.concat(boots[:k]), boots[k:]
        eval_set, explain = (pd.concat(rest[:e]) if e else None), rest[e:]
        source = "resampled samples"

    if len(regulators) == 0:
        raise ValueError("No regulators found in the data.")

    samples = train.index.unique()
    n_test = int(round(args.test_fraction * len(samples)))
    test_samples = [str(x) if isinstance(x, np.str_) else x for x in rng.choice(samples, size=n_test, replace=False)]
    if args.test_fraction > 0 and (n_test == 0 or n_test == len(samples)):
        parser.error(f"--test_fraction {args.test_fraction} leaves no train or no test samples ({len(samples)} samples).")

    print(f"{source}: {len(regulators)} regulators, {len(targets)} targets, "
          f"{args.n_bootstraps} bootstraps ({k} train -> {train.shape[0]} rows, {e} early-stopping eval, {n_explain} explain), "
          f"model={args.model}, {n_test} of {len(samples)} samples held out")

    timings = []
    edges = []
    boots_shap = []
    for i, target in enumerate(targets):
        timing, target_edges, target_boot = run_target(regulators, target, train, explain, seed=DEMON_SEED,
                                                       model_type=args.model, test_samples=test_samples,
                                                       early_stopping_rounds=args.early_stopping_rounds,
                                                       eval_set=eval_set)
        timings.append(timing)
        edges.append(target_edges)
        boots_shap.append(target_boot)
        print(f"[{i + 1}/{len(targets)}] {target}: trees={timing['n_trees']} "
              f"train={timing['t_train']:.2f}s shap={timing['t_shap']:.2f}s approx={timing['t_approx_shap']:.2f}s "
              f"(init={timing['t_explainer_init']:.3f}s)")

    timings = pd.DataFrame(timings)
    edges = pd.concat(edges, ignore_index=True)
    boots_shap = pd.concat(boots_shap, ignore_index=True)
    agreement = edge_agreement(edges, boots_shap, args.max_log2_rank_ratio)

    os.makedirs(args.results_path, exist_ok=True)
    timings.to_csv(op.join(args.results_path, f"{args.experiment_name}.timings.tsv"), sep='\t', index=False)
    edges.to_csv(op.join(args.results_path, f"{args.experiment_name}.edges.tsv"), sep='\t', index=False)
    boots_shap.to_csv(op.join(args.results_path, f"{args.experiment_name}.shap_bootstraps.tsv"), sep='\t', index=False)
    agreement.to_csv(op.join(args.results_path, f"{args.experiment_name}.agreement.tsv"), sep='\t', index=False)

    total_train = timings['t_train'].sum()
    total_init = timings['t_explainer_init'].sum()
    total_shap = timings['t_shap'].sum()
    total_approx = timings['t_approx_shap'].sum()
    print("\n=== Summary (all targets) ===")
    print(f"{f'training ({k} bootstraps)':<30}: {total_train:8.2f} s")
    print(f"explainer init                : {total_init:8.2f} s")
    print(f"{f'SHAP ({n_explain} bootstraps)':<30}: {total_shap:8.2f} s  "
          f"overhead {(total_init + total_shap) / total_train:6.2f}x training")
    print(f"{f'approx SHAP ({n_explain} bootstraps)':<30}: {total_approx:8.2f} s  "
          f"overhead {(total_init + total_approx) / total_train:6.2f}x training")
    print(f"SHAP / approx SHAP            : {total_shap / total_approx:8.2f}x")
    print(f"median Spearman(SHAP, approx) : {timings['spearman_shap_vs_approx'].median():8.3f}")
    print(f"median Spearman(SHAP, import.): {timings['spearman_shap_vs_importance'].median():8.3f}")
    print(f"max additivity error          : {timings['max_additivity_error'].max():.2e}")
    print(f"median MSE train / explain    : {timings['mse_train'].median():.3f} / {timings['mse_explain_mean'].median():.3f} "
          f"(no-model {timings['mse_null_explain'].median():.3f})")
    if e:
        print(f"median MSE eval bootstraps    : {timings['mse_eval'].median():.3f}")
    if n_test:
        print(f"median MSE held-out samples   : {timings['mse_test'].median():.3f} (training bootstraps) / "
              f"{timings['mse_explain_test_mean'].median():.3f} (explained bootstraps), "
              f"no-model {timings['mse_null_test'].median():.3f}")


if __name__ == '__main__':
    main()
