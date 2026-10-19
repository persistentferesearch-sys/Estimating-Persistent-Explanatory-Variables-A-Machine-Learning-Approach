from __future__ import annotations

import argparse
import hashlib
import json
import sys
import traceback
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from doubleml import DoubleMLData, DoubleMLPLR
from doubleml import __version__ as doubleml_version
from scipy import sparse
from scipy import __version__ as scipy_version
from scipy.stats import chi2
from sklearn.dummy import DummyRegressor
from sklearn.linear_model import LassoCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn import __version__ as sklearn_version

BACKEND_VERSION = "0.8.4-rc14"


def words(value: str | None) -> list[str]:
    return [x for x in (value or "").split() if x]


def stata_num(x: float) -> str:
    if x is None or not np.isfinite(x):
        return "."
    return f"{float(x):.17g}"


def stata_name(name: str) -> str:
    keep = [c if c.isalnum() or c == "_" else "_" for c in str(name)]
    out = "".join(keep).strip("_")
    if not out:
        out = "v"
    if out[0].isdigit():
        out = f"v_{out}"
    return out[:32]


def read_data(path: str, cols: list[str]) -> pd.DataFrame:
    df = pd.read_stata(path, convert_categoricals=False)
    missing = [c for c in cols if c and c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in input .dta: {missing}")
    return df.dropna(subset=[c for c in cols if c]).reset_index(drop=True)


def sparse_design(df: pd.DataFrame, absorb: list[str], controls: list[str]):
    blocks, names, groups = [], [], []
    n = len(df)
    for col in absorb:
        codes, levels = pd.factorize(df[col], sort=True)
        if (codes < 0).any():
            raise ValueError(f"Missing values remain in absorb variable {col}")
        block = sparse.csr_matrix((np.ones(n), (np.arange(n), codes)), shape=(n, len(levels)))
        blocks.append(block)
        names.extend([f"{col}_{stata_name(v)}" for v in levels])
        groups.extend([col] * len(levels))

    if controls:
        x = df[controls].to_numpy(dtype=float)
        blocks.append(sparse.csr_matrix(x))
        names.extend(controls)
        groups.extend(["numeric"] * len(controls))

    if not blocks:
        return sparse.csr_matrix((n, 0)), [], []
    return sparse.hstack(blocks, format="csr"), names, groups


def lasso_model(args, n_train: int):
    cv = min(max(2, int(args.lassocv_folds)), n_train)
    alphas = int(args.lassocv_alphas)
    version_bits = sklearn_version.split(".")
    major_minor = tuple(int(x) for x in version_bits[:2] if x.isdigit())
    alpha_args = {"alphas": alphas}
    if major_minor and major_minor < (1, 7):
        alpha_args = {"n_alphas": alphas, "alphas": None}
    return make_pipeline(
        StandardScaler(with_mean=False),
        LassoCV(
            cv=cv,
            eps=float(args.lassocv_eps),
            **alpha_args,
            fit_intercept=True,
            max_iter=int(args.lassocv_maxiter),
            tol=float(args.lassocv_tol),
            n_jobs=None,
            random_state=int(args.rseed),
        ),
    )


def fit_predict(X, y, train, pred, args, groups):
    train = np.asarray(train, dtype=int)
    pred = np.asarray(pred, dtype=int)
    if X.shape[1] == 0 or len(train) < 3 or np.nanstd(y[train]) < 1e-12:
        return np.repeat(float(np.nanmean(y[train])), len(pred)), Counter()

    model = lasso_model(args, len(train))
    model.fit(X[train], y[train])
    pred_y = model.predict(X[pred])
    coefs = model.named_steps["lassocv"].coef_
    sel = np.flatnonzero(np.abs(coefs) > float(args.coef_tol))
    return pred_y, Counter(groups[i] for i in sel)


def soften_binary(target, pred):
    vals = np.unique(np.asarray(target, dtype=float)[np.isfinite(target)])
    if len(vals) <= 2 and np.isin(vals, [0.0, 1.0]).all():
        return np.clip(pred, 1e-6, 1 - 1e-6)
    return pred


def soften_binary_matrix(targets, preds):
    out = np.asarray(preds, dtype=float).copy()
    for j in range(out.shape[1]):
        out[:, j] = soften_binary(targets[:, j], out[:, j])
    return out


def nuisance_design(X, d, j, treat_names, groups):
    """Match DoubleML's multi-treatment rule: other treatments enter X."""
    other = [h for h in range(d.shape[1]) if h != j]
    if not other:
        return X, groups
    extra = sparse.csr_matrix(d[:, other])
    extra_groups = [f"treatment:{treat_names[h]}" for h in other]
    return sparse.hstack((X, extra), format="csr"), groups + extra_groups


def nuisance_predictions(X, y, d, train, pred, args, groups, treat_names):
    n_pred, q = len(pred), d.shape[1]
    lhats, mhats, selected = np.empty((n_pred, q)), np.empty((n_pred, q)), []
    for j, name in enumerate(treat_names):
        design, design_groups = nuisance_design(X, d, j, treat_names, groups)
        lhats[:, j], lsel = fit_predict(design, y, train, pred, args, design_groups)
        mhat, msel = fit_predict(design, d[:, j], train, pred, args, design_groups)
        mhats[:, j] = soften_binary(d[:, j], mhat)
        selected.append({"treatment": name, "ml_l": dict(lsel), "ml_m": dict(msel)})
    return lhats, mhats, selected


def full_predictions(X, y, d, args, groups, treat_names):
    idx = np.arange(X.shape[0])
    return nuisance_predictions(X, y, d, idx, idx, args, groups, treat_names)


def fold_predictions(X, y, d, fold, args, groups, treat_names):
    n = X.shape[0]
    lhat = np.empty((n, d.shape[1]))
    mhats = np.empty((n, d.shape[1]))
    selected = []
    for f in sorted(pd.unique(fold)):
        test = np.where(fold == f)[0]
        train = np.where(fold != f)[0]
        lhat[test], mhats[test], fold_sel = nuisance_predictions(
            X, y, d, train, test, args, groups, treat_names
        )
        selected.append({"fold": int(f), "selected": fold_sel})
    return lhat, mhats, selected


def split_predictions(X, y, d, smpls, args, groups, treat_names):
    n = X.shape[0]
    lhat = np.empty((n, d.shape[1]))
    mhats = np.empty((n, d.shape[1]))
    coverage = np.zeros(n, dtype=int)
    selected = []
    for f, (train, test) in enumerate(smpls, start=1):
        train, test = np.asarray(train), np.asarray(test)
        if not len(train) or not len(test):
            raise ValueError(f"Empty train/test sample in score fold {f}")
        lhat[test], mhats[test], fold_sel = nuisance_predictions(
            X, y, d, train, test, args, groups, treat_names
        )
        coverage[test] += 1
        selected.append({"fold": f, "selected": fold_sel})
    if not np.all(coverage == 1):
        raise ValueError("DoubleML score folds do not cover every row exactly once")
    return lhat, mhats, selected


def row_folds_from_column(fold):
    idx = np.arange(len(fold))
    out = []
    for f in sorted(pd.unique(fold)):
        test = np.where(fold == f)[0]
        train = np.setdiff1d(idx, test, assume_unique=True)
        out.append((train, test))
    return out


def observation_split_audit(fold, k, seed):
    """Verify the observation-level XPO outer folds and record their geometry."""
    smpls = row_folds_from_column(fold)
    if len(smpls) != int(k):
        raise ValueError(f"Expected {k} observation score folds, found {len(smpls)}")

    coverage = np.zeros(len(fold), dtype=int)
    folds = []
    for f, (train, test) in enumerate(smpls, start=1):
        if not len(train) or not len(test):
            raise ValueError(f"Empty train/test sample in observation score fold {f}")
        if np.intersect1d(train, test).size:
            raise ValueError(f"Row overlap in observation score fold {f}")
        coverage[test] += 1
        folds.append({
            "fold": f,
            "n_train": int(len(train)),
            "n_test": int(len(test)),
            "train_index_sha256": index_hash(train),
            "test_index_sha256": index_hash(test),
        })
    if not np.all(coverage == 1):
        raise ValueError("Observation score folds do not cover every row exactly once")

    return {
        "seed": int(seed),
        "inference_contract": "observation_dml",
        "crossfit_unit": "observation",
        "nuisance_crossfit": True,
        "folds_per_dimension": len(smpls),
        "score_folds": len(smpls),
        "all_disjoint": True,
        "coverage_min": int(coverage.min()),
        "coverage_max": int(coverage.max()),
        "min_train": min(x["n_train"] for x in folds),
        "max_train": max(x["n_train"] for x in folds),
        "min_test": min(x["n_test"] for x in folds),
        "max_test": max(x["n_test"] for x in folds),
        "folds": folds,
    }


def make_dml(depvar, treat_names, y, d, clusters):
    data = {depvar: y, "_x_stub": np.zeros(len(y))}
    for j, t in enumerate(treat_names):
        data[t] = d[:, j]
    if clusters is not None:
        cluster_cols = []
        for j in range(clusters.shape[1]):
            name = f"_cluster{j + 1}"
            data[name] = clusters[:, j]
            cluster_cols.append(name)
        return DoubleMLData(
            pd.DataFrame(data),
            y_col=depvar,
            d_cols=treat_names,
            x_cols=["_x_stub"],
            cluster_cols=cluster_cols,
        )
    return DoubleMLData(pd.DataFrame(data), y_col=depvar, d_cols=treat_names, x_cols=["_x_stub"])


def external_predictions(treat_names, lhat, mhats):
    return {
        t: {"ml_l": lhat[:, j][:, None], "ml_m": mhats[:, j][:, None]}
        for j, t in enumerate(treat_names)
    }


def fit_doubleml(args, y, d, treat_names, lhat, mhats, fold):
    data = make_dml(args.depvar, treat_names, y, d, None)
    obj = DoubleMLPLR(
        data,
        ml_l=DummyRegressor(),
        ml_m=DummyRegressor(),
        n_folds=max(2, int(args.xfolds)),
        draw_sample_splitting=False,
    )

    if args.cmd == "sparse_xporeg":
        obj.set_sample_splitting(row_folds_from_column(fold))
    else:
        idx = np.arange(len(y))
        obj.set_sample_splitting((idx, idx))

    obj.fit(external_predictions=external_predictions(treat_names, lhat, mhats))
    return obj


def index_hash(values):
    raw = np.asarray(values, dtype=np.int64).tobytes()
    return hashlib.sha256(raw).hexdigest()


def cluster_split_audit(obj, clusters, cluster_names, k):
    smpls, cluster_smpls = obj.smpls[0], obj.smpls_cluster[0]
    m, n = len(cluster_names), len(clusters)
    if len(smpls) != k**m or len(cluster_smpls) != k**m:
        raise ValueError(f"Expected {k**m} cluster score folds")

    coverage = np.zeros(n, dtype=int)
    partitions = [[] for _ in range(m)]
    folds = []
    for f, ((train, test), (train_c, test_c)) in enumerate(zip(smpls, cluster_smpls), start=1):
        if not len(train) or not len(test):
            raise ValueError(f"Empty train/test sample in cluster score fold {f}")
        if np.intersect1d(train, test).size:
            raise ValueError(f"Row overlap in cluster score fold {f}")
        part_ids, disjoint = [], []
        train_mask = np.ones(n, dtype=bool)
        test_mask = np.ones(n, dtype=bool)
        for j in range(m):
            tr, te = np.asarray(train_c[j]), np.asarray(test_c[j])
            ok = not np.intersect1d(tr, te).size
            ok &= np.array_equal(np.union1d(tr, te), np.unique(clusters[:, j]))
            ok &= np.isin(clusters[train, j], tr).all() and np.isin(clusters[test, j], te).all()
            disjoint.append(bool(ok))
            train_mask &= np.isin(clusters[:, j], tr)
            test_mask &= np.isin(clusters[:, j], te)
            sig = tuple(te.tolist())
            known = [tuple(x) for x in partitions[j]]
            if sig not in known:
                partitions[j].append(te.tolist())
                known.append(sig)
            part_ids.append(known.index(sig) + 1)
        if not all(disjoint):
            raise ValueError(f"Cluster overlap in score fold {f}")
        rows_match = np.array_equal(np.flatnonzero(train_mask), train)
        rows_match &= np.array_equal(np.flatnonzero(test_mask), test)
        if not rows_match:
            raise ValueError(f"Row/cluster split mismatch in score fold {f}")
        coverage[test] += 1
        folds.append({
            "fold": f,
            "partition_ids": part_ids,
            "n_train": int(len(train)),
            "n_test": int(len(test)),
            "train_index_sha256": index_hash(train),
            "test_index_sha256": index_hash(test),
        })

    if not np.all(coverage == 1):
        raise ValueError("Cluster score folds do not cover every row exactly once")
    combinations = [tuple(item["partition_ids"]) for item in folds]
    expected = {
        tuple(x + 1 for x in combo)
        for combo in np.ndindex(*(k for _ in range(m)))
    }
    if len(set(combinations)) != k**m or set(combinations) != expected:
        raise ValueError("Cluster score folds are not the complete K^m Cartesian product")
    n_clusters = []
    repetitions = []
    for j, parts in enumerate(partitions):
        all_values = np.unique(clusters[:, j])
        joined = np.concatenate([np.asarray(x) for x in parts])
        if len(parts) != k or len(joined) != len(all_values) or not np.array_equal(np.unique(joined), all_values):
            raise ValueError(f"Cluster dimension {j + 1} is not a {k}-way partition")
        counts = Counter(combo[j] for combo in combinations)
        if set(counts.values()) != {k ** (m - 1)}:
            raise ValueError(f"Cluster dimension {j + 1} partitions have wrong repetition counts")
        repetitions.append(dict(sorted(counts.items())))
        n_clusters.append(int(len(all_values)))

    train_n = [x["n_train"] for x in folds]
    test_n = [x["n_test"] for x in folds]
    return {
        "cluster_names": cluster_names,
        "cluster_dims": m,
        "folds_per_dimension": k,
        "score_folds": len(folds),
        "n_clusters": n_clusters,
        "n_cluster_intersections": int(pd.DataFrame(clusters).drop_duplicates().shape[0]),
        "test_partitions": partitions,
        "partition_repetitions": repetitions,
        "cartesian_complete": True,
        "all_disjoint": True,
        "coverage_min": int(coverage.min()),
        "coverage_max": int(coverage.max()),
        "min_train": min(train_n),
        "max_train": max(train_n),
        "min_test": min(test_n),
        "max_test": max(test_n),
        "folds": folds,
    }


def fit_cluster_xpo(args, X, y, d, treat_names, clusters, cluster_names, groups):
    k = int(args.xfolds)
    for name, values in zip(cluster_names, clusters.T):
        if len(pd.unique(values)) < k:
            raise ValueError(f"xfolds({k}) exceeds clusters in {name}")
    np.random.seed(int(args.rseed))
    obj = DoubleMLPLR(
        make_dml(args.depvar, treat_names, y, d, clusters),
        ml_l=DummyRegressor(),
        ml_m=DummyRegressor(),
        n_folds=k,
    )
    lhat, mhats, selected = split_predictions(
        X, y, d, obj.smpls[0], args, groups, treat_names
    )
    if not np.isfinite(lhat).all() or not np.isfinite(mhats).all():
        raise ValueError("Nonfinite external prediction")
    obj.fit(external_predictions=external_predictions(treat_names, lhat, mhats))
    audit = cluster_split_audit(obj, clusters, cluster_names, k)
    audit.update({
        "seed": int(args.rseed),
        "inference_contract": "cluster_dml",
        "crossfit_unit": "cluster",
        "nuisance_crossfit": True,
    })
    return obj, selected, audit


def fit_cluster_po(args, y, d, treat_names, lhat, mhats, clusters, cluster_names):
    """Keep PO nuisance fits and coefficients observation weighted."""
    n_clusters = [len(pd.unique(clusters[:, j])) for j in range(clusters.shape[1])]
    if min(n_clusters) < 2:
        raise ValueError("Cluster-robust inference requires at least two clusters per dimension")
    obj = fit_doubleml(args, y, d, treat_names, lhat, mhats, None)
    audit = {
        "inference_contract": (
            "conventional_multiway_cluster" if len(cluster_names) == 2
            else "conventional_cluster_robust"
        ),
        "crossfit_unit": "none",
        "nuisance_crossfit": False,
        "cluster_dims": len(cluster_names),
        "n_clusters": n_clusters,
        "n_cluster_intersections": int(pd.DataFrame(clusters).drop_duplicates().shape[0]),
    }
    return obj, audit


def full_covariance(obj, clusters=None):
    """Extend DoubleML's scalar score variance to all treatment pairs."""
    psi_raw, deriv_raw = np.asarray(obj.psi), np.asarray(obj.psi_deriv)
    if obj.n_rep != 1 or psi_raw.ndim != 3 or psi_raw.shape[1] != 1:
        raise ValueError("Full covariance currently requires exactly one DoubleML repetition")
    if psi_raw.shape != deriv_raw.shape:
        raise ValueError("DoubleML score and derivative arrays have different shapes")
    psi, deriv = psi_raw[:, 0, :], deriv_raw[:, 0, :]
    if clusters is None:
        jacobian = deriv.mean(axis=0)
        meat = psi.T @ psi / len(psi)
        scale = len(psi)
    else:
        m = clusters.shape[1]
        if m not in {1, 2}:
            raise ValueError("Full cluster covariance requires one or two cluster dimensions")
        smpls, cluster_smpls = obj.smpls[0], obj.smpls_cluster[0]
        k = int(round(len(smpls) ** (1 / m)))
        if k**m != len(smpls):
            raise ValueError("Cluster score folds are not a K^m Cartesian product")
        meat = np.zeros((psi.shape[1], psi.shape[1]))
        jacobian = np.zeros(psi.shape[1])
        for (_, test), (_, test_clusters) in zip(smpls, cluster_smpls):
            if m == 1:
                first = np.asarray(test_clusters[0])
                const = 1 / len(first)
                for value in first:
                    score_sum = psi[clusters[:, 0] == value].sum(axis=0)
                    meat += const * np.outer(score_sum, score_sum)
                jacobian += deriv[test].sum(axis=0) / len(first)
            else:
                first, second = map(np.asarray, test_clusters)
                const = min(len(first), len(second)) / (len(first) * len(second)) ** 2
                for value in first:
                    rows = (clusters[:, 0] == value) & np.isin(clusters[:, 1], second)
                    score_sum = psi[rows].sum(axis=0)
                    meat += const * np.outer(score_sum, score_sum)
                for value in second:
                    rows = (clusters[:, 1] == value) & np.isin(clusters[:, 0], first)
                    score_sum = psi[rows].sum(axis=0)
                    meat += const * np.outer(score_sum, score_sum)
                jacobian += deriv[test].sum(axis=0) / (len(first) * len(second))
        meat /= k**m
        jacobian /= k**m
        scale = min(len(pd.unique(clusters[:, j])) for j in range(m))

    jacobian_tol = np.sqrt(np.finfo(float).eps) * max(1.0, np.abs(deriv).mean(axis=0).max())
    if np.any(np.abs(jacobian) <= jacobian_tol):
        raise ValueError(f"Weak or zero PLR score Jacobian (tolerance {jacobian_tol:.3g})")
    covariance = meat / (scale * np.outer(jacobian, jacobian))
    covariance = (covariance + covariance.T) / 2
    if not np.isfinite(covariance).all():
        raise ValueError("Nonfinite full score covariance")
    eigenvalues = np.linalg.eigvalsh(covariance)
    eigen_tol = 1e-10 * max(1.0, np.abs(eigenvalues).max())
    if eigenvalues.min() < -eigen_tol:
        raise ValueError("Full score covariance is materially indefinite")
    if not np.allclose(np.diag(covariance), np.square(obj.se), rtol=1e-9, atol=1e-12):
        raise ValueError("Full covariance diagonal disagrees with DoubleML standard errors")
    return covariance


def conventional_cluster_covariance(obj, clusters):
    """CR0/CGM covariance for observation-weighted PO scores."""
    psi_raw, deriv_raw = np.asarray(obj.psi), np.asarray(obj.psi_deriv)
    if obj.n_rep != 1 or psi_raw.ndim != 3 or psi_raw.shape[1] != 1:
        raise ValueError("Full covariance currently requires exactly one DoubleML repetition")
    if psi_raw.shape != deriv_raw.shape:
        raise ValueError("DoubleML score and derivative arrays have different shapes")
    psi, deriv = psi_raw[:, 0, :], deriv_raw[:, 0, :]
    clusters = np.asarray(clusters)
    if clusters.ndim == 1:
        clusters = clusters[:, None]
    if clusters.shape[1] not in {1, 2}:
        raise ValueError("Conventional cluster covariance requires one or two dimensions")

    def grouped_outer(columns):
        keys = pd.MultiIndex.from_arrays([clusters[:, j] for j in columns])
        codes, levels = pd.factorize(keys, sort=False)
        sums = np.zeros((len(levels), psi.shape[1]))
        np.add.at(sums, codes, psi)
        return sums.T @ sums

    meat = grouped_outer([0])
    if clusters.shape[1] == 2:
        meat += grouped_outer([1]) - grouped_outer([0, 1])
    jacobian = deriv.mean(axis=0)
    jacobian_tol = np.sqrt(np.finfo(float).eps) * max(1.0, np.abs(deriv).mean(axis=0).max())
    if np.any(np.abs(jacobian) <= jacobian_tol):
        raise ValueError(f"Weak or zero PLR score Jacobian (tolerance {jacobian_tol:.3g})")
    covariance = meat / (len(psi) ** 2 * np.outer(jacobian, jacobian))
    covariance = (covariance + covariance.T) / 2
    if not np.isfinite(covariance).all():
        raise ValueError("Nonfinite full score covariance")
    eigenvalues = np.linalg.eigvalsh(covariance)
    eigen_tol = 1e-10 * max(1.0, np.abs(eigenvalues).max())
    if eigenvalues.min() < -eigen_tol:
        raise ValueError("Conventional cluster covariance is materially indefinite")
    return covariance


def write_postdo(path, args, names, coef, covariance, n_obs, n_features, fold_audit=None):
    variances = np.diag(covariance)
    if not np.isfinite(covariance).all() or np.any(variances <= 0):
        raise ValueError("Final treatment covariance has nonpositive or nonfinite variance")
    colnames = " ".join(names)
    k = len(names)
    rank = int(np.linalg.matrix_rank(covariance))
    if rank != k:
        raise ValueError(f"Final treatment covariance is rank deficient ({rank} < {k})")
    wald = float(np.asarray(coef) @ np.linalg.pinv(covariance) @ np.asarray(coef))
    p_value = float(chi2.sf(wald, rank)) if rank else np.nan
    esample = f" esample({args.esample_var})" if getattr(args, "esample_var", "") else ""
    title = "Cross-fit " if args.cmd == "sparse_xporeg" else ""
    title += "partialling-out linear regression (Python DoubleML backend)"
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("ereturn clear\n")
        fh.write(f"matrix __sp_b = J(1, {k}, .)\n")
        for i, value in enumerate(coef, start=1):
            fh.write(f"matrix __sp_b[1, {i}] = {stata_num(value)}\n")
        fh.write(f"matrix colnames __sp_b = {colnames}\n")
        fh.write(f"matrix __sp_V = J({k}, {k}, 0)\n")
        for i in range(k):
            for j in range(k):
                fh.write(
                    f"matrix __sp_V[{i + 1}, {j + 1}] = {stata_num(covariance[i, j])}\n"
                )
        fh.write(f"matrix rownames __sp_V = {colnames}\n")
        fh.write(f"matrix colnames __sp_V = {colnames}\n")
        fh.write(f"ereturn post __sp_b __sp_V, obs({int(n_obs)}){esample}\n")
        fh.write("matrix __sp_cov = e(V)\n")
        fh.write("ereturn matrix sparse_covariance = __sp_cov\n")
        fh.write(f"ereturn local cmd \"{args.cmd}\"\n")
        fh.write(f'ereturn local title "{title}"\n')
        fh.write(f"ereturn local depvar \"{args.depvar}\"\n")
        fh.write(f"ereturn local varsofinterest \"{colnames}\"\n")
        fh.write("ereturn local sparse_backend \"python\"\n")
        fh.write(f'ereturn local sparse_version "{BACKEND_VERSION}"\n')
        fh.write('ereturn local python_model "DoubleMLPLR"\n')
        fh.write(f'ereturn local python_doubleml_version "{doubleml_version}"\n')
        fh.write(f'ereturn local python_sklearn_version "{sklearn_version}"\n')
        fh.write('ereturn local nuisance_inner_cv "row"\n')
        fh.write(f"ereturn local selection \"{args.selection}\"\n")
        fh.write(f"ereturn local learner \"{args.learner}\"\n")
        fh.write(f"ereturn local vce \"{args.vce}\"\n")
        if args.cluster:
            fh.write(f"ereturn local clustvar \"{words(args.cluster)[0]}\"\n")
            fh.write(f"ereturn local clustvars \"{args.cluster}\"\n")
            vcetype = "Multiway cluster-robust" if len(words(args.cluster)) == 2 else "Cluster-robust"
        else:
            vcetype = "Robust"
        fh.write(f'ereturn local vcetype "{vcetype}"\n')
        fh.write(f"ereturn scalar N = {int(n_obs)}\n")
        fh.write(f"ereturn scalar rank = {rank}\n")
        fh.write(f"ereturn scalar df = {rank}\n")
        fh.write(f"ereturn scalar df_m = {rank}\n")
        fh.write('ereturn local chi2type "Wald"\n')
        fh.write(f"ereturn scalar chi2 = {stata_num(wald)}\n")
        fh.write(f"ereturn scalar p = {stata_num(p_value)}\n")
        fh.write(f"ereturn scalar k_varsofinterest = {k}\n")
        fh.write(f"ereturn scalar sparse_features = {int(n_features)}\n")
        fh.write(f"ereturn scalar sparse_variance = {stata_num(variances[0])}\n")
        fh.write(f"ereturn scalar n_treatments = {k}\n")
        min_eigen = np.linalg.eigvalsh(covariance).min()
        fh.write(f"ereturn scalar sparse_vcov_min_eigen = {stata_num(min_eigen)}\n")
        fh.write("ereturn scalar finite = 1\n")
        fh.write(f"ereturn scalar post_rank_valid = {int(rank == k)}\n")
        fh.write("ereturn scalar sparse_covariance_full = 1\n")
        fh.write("ereturn scalar sparse_vcov_full = 1\n")
        if fold_audit:
            fh.write(f'ereturn local inference_contract "{fold_audit["inference_contract"]}"\n')
            fh.write(f'ereturn local crossfit_unit "{fold_audit["crossfit_unit"]}"\n')
            fh.write(f'ereturn scalar nuisance_crossfit = {int(fold_audit["nuisance_crossfit"])}\n')
            if "score_folds" in fold_audit:
                fh.write(f"ereturn scalar n_xfolds = {fold_audit['folds_per_dimension']}\n")
                fh.write(f"ereturn scalar n_score_folds = {fold_audit['score_folds']}\n")
                fh.write(f"ereturn scalar fold_coverage_min = {fold_audit['coverage_min']}\n")
                fh.write(f"ereturn scalar fold_coverage_max = {fold_audit['coverage_max']}\n")
                fh.write(f"ereturn scalar fold_min_train = {fold_audit['min_train']}\n")
                fh.write(f"ereturn scalar fold_min_test = {fold_audit['min_test']}\n")
                if fold_audit["crossfit_unit"] == "observation":
                    fh.write("ereturn scalar observation_disjoint = 1\n")
                    fh.write(f"ereturn scalar observation_fold_seed = {fold_audit['seed']}\n")
                else:
                    fh.write(f'ereturn local foldclustvars "{args.cluster}"\n')
                    fh.write("ereturn scalar cluster_disjoint = 1\n")
                    fh.write(f"ereturn scalar cluster_fold_seed = {fold_audit['seed']}\n")
                    fh.write(f"ereturn scalar n_xfolds_per_cluster = {fold_audit['folds_per_dimension']}\n")
                    fh.write(f"ereturn scalar cluster_coverage_min = {fold_audit['coverage_min']}\n")
                    fh.write(f"ereturn scalar cluster_coverage_max = {fold_audit['coverage_max']}\n")
                    fh.write(f"ereturn scalar cluster_min_train = {fold_audit['min_train']}\n")
                    fh.write(f"ereturn scalar cluster_min_test = {fold_audit['min_test']}\n")
            if "cluster_dims" in fold_audit:
                fh.write(f"ereturn scalar n_cluster_dims = {fold_audit['cluster_dims']}\n")
                fh.write(f"ereturn scalar N_clust = {min(fold_audit['n_clusters'])}\n")
                for j, value in enumerate(fold_audit["n_clusters"], start=1):
                    fh.write(f"ereturn scalar N_clust{j} = {value}\n")
                if fold_audit["cluster_dims"] == 2:
                    fh.write(f"ereturn scalar N_clust12 = {fold_audit['n_cluster_intersections']}\n")


def run(args):
    if args.cmd not in {"sparse_poreg", "sparse_xporeg"}:
        raise ValueError("--cmd must be sparse_poreg or sparse_xporeg")
    if int(args.xfolds) < 2:
        raise ValueError("--xfolds must be at least 2")
    if args.vce not in {"robust", "cluster"}:
        raise ValueError("--vce must be robust or cluster")
    esample_var = getattr(args, "esample_var", "")
    if esample_var and (
        len(esample_var) > 32
        or not (esample_var[0].isalpha() or esample_var[0] == "_")
        or any(not (c.isalnum() or c == "_") for c in esample_var)
    ):
        raise ValueError("--esample-var must be a valid Stata variable name")
    treat = words(args.treat)
    if not treat or len(set(treat)) != len(treat):
        raise ValueError("--treat requires one or more distinct columns")
    treat_names = [stata_name(x) for x in words(args.treat_labels) or treat]
    if len(treat_names) != len(treat):
        treat_names = [stata_name(x) for x in treat]
    if len(set(treat_names)) != len(treat_names):
        raise ValueError(f"Treatment names collide after Stata sanitization: {treat_names}")

    absorb = words(args.absorb)
    controls = words(args.controls)
    cluster_names = words(args.cluster)
    if len(cluster_names) > 2 or len(set(cluster_names)) != len(cluster_names):
        raise ValueError("--cluster accepts one or two distinct columns")
    if (args.vce == "cluster") != bool(cluster_names):
        raise ValueError("vce(cluster) requires --cluster, and robust forbids it")
    needed = [args.depvar] + treat + absorb + controls
    needed.extend(cluster_names)
    if args.fold_var:
        needed.append(args.fold_var)

    df = read_data(args.input, needed)
    y = df[args.depvar].to_numpy(dtype=float)
    d = df[treat].to_numpy(dtype=float)
    X, feature_names, feature_groups = sparse_design(df, absorb, controls)

    clusters = df[cluster_names].to_numpy() if cluster_names else None
    fold_audit = None
    if args.cmd == "sparse_xporeg" and clusters is not None:
        obj, selected, fold_audit = fit_cluster_xpo(
            args, X, y, d, treat_names, clusters, cluster_names, feature_groups
        )
    elif args.cmd == "sparse_xporeg":
        if not args.fold_var:
            raise ValueError("Robust sparse_xporeg requires --fold-var")
        fold = df[args.fold_var].to_numpy()
        fold_audit = observation_split_audit(fold, args.xfolds, args.rseed)
        lhat, mhats, selected = fold_predictions(
            X, y, d, fold, args, feature_groups, treat_names
        )
        mhats = soften_binary_matrix(d, mhats)
        obj = fit_doubleml(args, y, d, treat_names, lhat, mhats, fold)
    else:
        fold = None
        lhat, mhats, selected = full_predictions(
            X, y, d, args, feature_groups, treat_names
        )
        mhats = soften_binary_matrix(d, mhats)
        if clusters is None:
            obj = fit_doubleml(args, y, d, treat_names, lhat, mhats, fold)
        else:
            obj, fold_audit = fit_cluster_po(
                args, y, d, treat_names, lhat, mhats, clusters, cluster_names
            )
    if args.cmd == "sparse_poreg" and clusters is not None:
        covariance = conventional_cluster_covariance(obj, clusters)
    else:
        covariance = full_covariance(obj, clusters)
    write_postdo(
        args.postdo, args, treat_names, obj.coef, covariance, len(df), X.shape[1], fold_audit
    )
    posted_se = np.sqrt(np.diag(covariance))

    audit = {
        "cmd": args.cmd,
        "backend_version": BACKEND_VERSION,
        "N": int(len(df)),
        "features": int(X.shape[1]),
        "nnz": int(X.nnz),
        "absorb": absorb,
        "controls": controls,
        "treat": treat,
        "selection": args.selection,
        "learner": args.learner,
        "coef": obj.coef.tolist(),
        "se": posted_se.tolist(),
        "variance": np.diag(covariance).tolist(),
        "covariance": covariance.tolist(),
        "covariance_full": True,
        "covariance_max_diagonal_gap": float(
            np.max(np.abs(np.diag(covariance) - np.square(posted_se)))
        ),
        "doubleml_se_diagonal_gap": float(
            np.max(np.abs(np.diag(covariance) - np.square(obj.se)))
        ),
        "covariance_symmetric": bool(np.allclose(covariance, covariance.T)),
        "covariance_rank": int(np.linalg.matrix_rank(covariance)),
        "covariance_min_eigenvalue": float(np.linalg.eigvalsh(covariance).min()),
        "n_treatments": len(treat),
        "other_treatments_as_covariates": len(treat) > 1,
        "selected": selected,
        "folds": fold_audit if fold_audit and "score_folds" in fold_audit else None,
        "cluster_inference": fold_audit if fold_audit and "cluster_dims" in fold_audit else None,
        "cluster_folds": (
            fold_audit if fold_audit and "score_folds" in fold_audit and "cluster_dims" in fold_audit else None
        ),
        "versions": {
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy_version,
            "scikit_learn": sklearn_version,
            "doubleml": doubleml_version,
        },
    }
    Path(str(args.postdo) + ".audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cmd", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--postdo", required=True)
    parser.add_argument("--esample-var", default="")
    parser.add_argument("--depvar", required=True)
    parser.add_argument("--treat", required=True)
    parser.add_argument("--treat-labels", default="")
    parser.add_argument("--absorb", default="")
    parser.add_argument("--controls", default="")
    parser.add_argument("--selection", default="plugin")
    parser.add_argument("--vce", default="robust")
    parser.add_argument("--cluster", default="")
    parser.add_argument("--fold-var", default="")
    parser.add_argument("--xfolds", type=int, default=10)
    parser.add_argument("--rseed", type=int, default=1103)
    parser.add_argument("--max-iter", type=int, default=4000)
    parser.add_argument("--coef-tol", type=float, default=1e-3)
    parser.add_argument("--learner", default="lassocv")
    parser.add_argument("--lassocv-folds", type=int, default=5)
    parser.add_argument("--lassocv-alphas", type=int, default=100)
    parser.add_argument("--lassocv-eps", type=float, default=0.001)
    parser.add_argument("--lassocv-maxiter", type=int, default=5000)
    parser.add_argument("--lassocv-tol", type=float, default=1e-4)
    return parser.parse_args()


def main():
    args = parse_args()
    args.cmd = args.cmd.lower().strip()
    args.vce = args.vce.lower().strip()
    if args.learner.lower() != "lassocv":
        raise ValueError("Only learner(lassocv) is implemented in the Python backend V1")
    args.learner = args.learner.lower()
    args.selection = args.selection.lower()
    if args.selection != "lassocv":
        raise ValueError("Python backend currently implements selection(lassocv) only")
    args.cluster = args.cluster.strip()
    args.fold_var = args.fold_var.strip()
    run(args)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        postdo = None
        if "--postdo" in sys.argv:
            pos = sys.argv.index("--postdo") + 1
            if pos < len(sys.argv):
                postdo = sys.argv[pos]
        log = Path(str(postdo) + ".py.log") if postdo else Path("sparse_prl_xpo_backend.py.log")
        log.write_text(traceback.format_exc(), encoding="utf-8")
        raise
