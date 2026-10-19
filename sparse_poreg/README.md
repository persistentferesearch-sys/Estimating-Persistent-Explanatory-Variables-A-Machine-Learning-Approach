# Sparse partialing-out for empirical finance

This Stata package helps estimate the association between a variable of interest—such as R&D spending—and an outcome such as patenting, while selecting from many possible controls and firm or year indicators.

The two main commands estimate **linear models**:

- `sparse_poreg` uses partialing-out lasso (PRL).
- `sparse_xporeg` uses cross-fitted partialing-out lasso (XPO/DML).

Both can call Stata’s official estimator or use this package’s sparse Mata implementation. A separate **Poisson command**, `sparse_popoisson`, is under development for count outcomes.

## Install and try it

Download and extract the package. In Stata, change the working directory to the extracted folder—the one containing `stata.toc`—then run:

```stata
net install sparse_poreg, from("`c(pwd)'") replace
do "tests/finance_quickstart.do"
```

The quick start creates artificial firm-year data. It runs the official and Mata routes on the same R&D–patent example and prints `PASS` if their coefficient and standard error agree. It does not use your research data or require Python.

The package requires Stata 16 or newer. The comparisons described below were run in Stata 19.

## Your first model

Suppose `lnnpatent` is the outcome, `RDAT` is the variable whose coefficient you want, and `PERMCO` and `fyear` identify firms and years:

```stata
sparse_poreg lnnpatent RDAT, ///
    absorb(PERMCO fyear) ///
    controls(lnME lnAge TobinQ ROA Leverage) ///
    vce(cluster PERMCO) ///
    selection(bic, gridminok) ///
    backend(auto)
```

Here, `RDAT` receives a reported coefficient and standard error. The variables in `controls()` are **candidates for selection**. `absorb(PERMCO fyear)` makes the observed firm and year indicator levels candidates too; it is not the same instruction as absorbing fixed effects in an ordinary within regression. `vce(cluster PERMCO)` clusters standard errors by firm.

Every variable placed between the outcome and the comma is a variable of interest and receives its own coefficient. If you want only the R&D coefficient, put the other financial variables in `controls()`, as above.

For cross-fitted XPO, add a fold count and seed:

```stata
sparse_xporeg lnnpatent RDAT, ///
    absorb(PERMCO fyear) ///
    controls(lnME lnAge TobinQ ROA Leverage) ///
    vce(cluster PERMCO) ///
    selection(bic, gridminok) ///
    xfolds(5) rseed(1103) ///
    backend(auto)
```

With `backend(auto)`, XPO uses official Stata when the expanded design fits. It will not silently switch to the experimental Mata XPO route.

## Which route should I choose?

| Goal | Option |
|---|---|
| Start a routine analysis | `backend(auto)` |
| Reproduce Stata’s official command through this package | `backend(official)` |
| Compare the sparse calculation with official Stata on a manageable design | `backend(mata) strict(stata)` |
| Run a very large linear PRL design with sparse safety checks | `backend(mata) strict(production)` |
| Deliberately run sparse XPO | `backend(mata) frontierok strict(stata)` for a parity comparison, or `strict(production)` for the scalable route |

`strict(stata)` is an option **within** the Mata backend. It is intended for comparisons with Stata’s selection conventions. `strict(production)` avoids some dense official screening and applies the package’s production checks; it should not be described as universally identical to official Stata.

Sparse XPO requires `frontierok` because that route remains experimental. For the official-style comparisons here, its five folds split **observations**, even when standard errors are clustered by firm. Researchers who need whole firms kept together across training and holdout samples should consult `help sparse_xporeg` for `foldcluster()`; that is a different design.

## A comparison using the paper’s public data

The package was also checked against the [paper’s public replication data](https://github.com/persistentferesearch-sys/persistent_fe_remendy). The published specification treats R&D and 12 other variables as variables of interest, selects among firm and year indicators, clusters by firm, and uses BIC. XPO uses five folds.

| Model | Official result for RDAT | Sparse Mata result for RDAT |
|---|---:|---:|
| Linear PRL coefficient | 0.1986367 | 0.198636660737147 |
| Linear PRL cluster SE | 0.0181304 | 0.0181303830386044 |
| Linear XPO coefficient | 0.2131053 | 0.213105298267939 |
| Linear XPO cluster SE | 0.0136054 | 0.0136053698613053 |

The Mata PRL fit also matched the **local official fit** across all 13 coefficients and the full covariance matrix to numerical precision. On that machine, Mata PRL took about 8 minutes and the local official fit about 6 hours. The Mata XPO fit matches all 13 coefficient and standard-error pairs **as printed in the published official log**. Because that log rounds the numbers, it does not establish full-precision XPO covariance equality, and it does not provide a comparable official XPO estimation timer.

The commands, data preparation, and evidence are in [RDPatentComparison.md](RDPatentComparison.md). These results concern this recorded specification; changing the sample, variable order, indicator coding, selection rule, or folds may change the comparison.

## Poisson outcomes and other options

For an outcome that is a count rather than a linear-model outcome, see `help sparse_popoisson`. This Poisson command is a **development feature** with separate validation. `sparse_poipoisson` is an alias for it; neither command is interchangeable with the linear commands above.

The package also includes an optional Python/DoubleML route using LassoCV. Because it uses a different selection method, its results are a sensitivity analysis rather than an official Stata plugin/BIC replication. Most readers can begin with the Stata and Mata routes without installing Python.

For command options and postestimation, type `help sparse_poreg`, `help sparse_xporeg`, or `help sparse_popoisson` in Stata.

This is an independent third-party package, not affiliated with or endorsed by StataCorp LLC. The package code is distributed under the MIT license.
