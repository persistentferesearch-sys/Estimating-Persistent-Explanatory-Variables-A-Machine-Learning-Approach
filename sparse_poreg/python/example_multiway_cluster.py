"""Two-way cluster-disjoint PLR, following the DoubleML multiway example."""

import numpy as np
import pandas as pd
import doubleml
import scipy
import sklearn
from doubleml import DoubleMLData, DoubleMLPLR
from sklearn.base import clone
from sklearn.linear_model import LassoCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sparse_prl_xpo_backend import full_covariance


SEED = 20260717
K = 3
rng = np.random.default_rng(SEED)
g1_size = g2_size = 18
p = 10

# Complete two-way array with shocks along both cluster dimensions.
g1 = np.repeat(np.arange(g1_size), g2_size)
g2 = np.tile(np.arange(g2_size), g1_size)
n = len(g1)
a = (g1 + 2 * g2) % 3
X = (
    0.45 * rng.normal(size=(g1_size, p))[g1]
    + 0.45 * rng.normal(size=(g2_size, p))[g2]
    + 0.55 * rng.normal(size=(n, p))
)
beta = 0.5 ** np.arange(1, p + 1)
gamma = 0.4 ** np.arange(1, p + 1)
a_effect = np.array([-0.2, 0.0, 0.2])[a]
d1 = (
    X @ beta
    + 0.25 * a_effect
    + 0.30 * rng.normal(size=g1_size)[g1]
    + 0.30 * rng.normal(size=g2_size)[g2]
    + rng.normal(size=n)
)
d2 = (
    -0.35 * X @ beta
    + 0.20 * a_effect
    + 0.25 * rng.normal(size=g1_size)[g1]
    + 0.25 * rng.normal(size=g2_size)[g2]
    + rng.normal(size=n)
)
y = (
    d1 - 0.65 * d2
    + X @ gamma
    + 0.20 * a_effect
    + 0.30 * rng.normal(size=g1_size)[g1]
    + 0.30 * rng.normal(size=g2_size)[g2]
    + rng.normal(size=n)
)

x_cols = [f"x{j + 1}" for j in range(p)]
df = pd.DataFrame(X, columns=x_cols)
a_cols = []
for level in range(3):
    name = f"a_{level}"
    df[name] = (a == level).astype(float)
    a_cols.append(name)
df = df.assign(g1=g1, g2=g2, d1=d1, d2=d2, y=y)

# Stata analogue (same cluster split geometry, not numerical parity):
# sparse_xporeg y d1 d2, absorb(a) controls(x1-x10) selection(lassocv) \
#     vce(cluster g1 g2) foldcluster(g1 g2) xfolds(3) \
#     backend(python) learner(lassocv)
data = DoubleMLData(
    df,
    y_col="y",
    d_cols=["d1", "d2"],
    x_cols=x_cols + a_cols,
    cluster_cols=["g1", "g2"],
)
learner = make_pipeline(
    StandardScaler(),
    LassoCV(cv=3, alphas=50, max_iter=5000, tol=1e-5),
)

# DoubleML partitions each cluster dimension K ways, yielding K^2 score folds.
np.random.seed(SEED)
model = DoubleMLPLR(
    data,
    ml_l=clone(learner),
    ml_m=clone(learner),
    n_folds=K,
)
model.fit()

coef, se = model.coef, model.se
clusters = np.column_stack((g1, g2))
covariance = full_covariance(model, clusters)
assert np.isfinite(np.r_[coef, se, covariance.ravel()]).all() and np.all(se > 0)
assert covariance.shape == (2, 2)
assert np.allclose(covariance, covariance.T)
assert np.allclose(np.diag(covariance), se**2)
assert abs(covariance[0, 1]) > 0
assert model.n_folds == K**2
assert len(model.smpls[0]) == K**2
assert len(model.smpls_cluster[0]) == K**2

coverage = np.zeros(n, dtype=int)
partitions = [set(), set()]
for (train, test), (train_c, test_c) in zip(model.smpls[0], model.smpls_cluster[0]):
    assert len(train) and len(test) and not np.intersect1d(train, test).size
    train_mask = np.ones(n, dtype=bool)
    test_mask = np.ones(n, dtype=bool)
    for j in range(2):
        tr, te = np.asarray(train_c[j]), np.asarray(test_c[j])
        assert not np.intersect1d(tr, te).size
        assert set(np.union1d(tr, te)) == set(np.unique(clusters[:, j]))
        train_mask &= np.isin(clusters[:, j], tr)
        test_mask &= np.isin(clusters[:, j], te)
        partitions[j].add(tuple(te.tolist()))
    assert np.array_equal(np.flatnonzero(train_mask), train)
    assert np.array_equal(np.flatnonzero(test_mask), test)
    coverage[test] += 1
assert np.all(coverage == 1)
for j in range(2):
    values = [v for part in partitions[j] for v in part]
    assert len(partitions[j]) == K and set(values) == set(np.unique(clusters[:, j]))

print("coefficients=" + ",".join(f"{x:.8f}" for x in coef))
print("standard_errors=" + ",".join(f"{x:.8f}" for x in se))
print(f"covariance_offdiag={covariance[0, 1]:.8f}")
print(f"cluster_score_folds={model.n_folds}")
print(
    "versions="
    f"numpy:{np.__version__},pandas:{pd.__version__},scipy:{scipy.__version__},"
    f"scikit-learn:{sklearn.__version__},doubleml:{doubleml.__version__}"
)
