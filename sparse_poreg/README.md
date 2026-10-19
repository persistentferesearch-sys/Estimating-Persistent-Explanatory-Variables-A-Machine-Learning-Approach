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

## Quick start
```stata
version 19.0
clear all
set more off
args package_root
if trim(`"`package_root'"') != "" {
    local package_root = subinstr(`"`package_root'"', "\", "/", .)
    adopath ++ "`package_root'/stata"
}
capture which sparse_poreg
if _rc {
    display as error "sparse_poreg is not installed; follow START_HERE_FOR_COAUTHORS.txt"
    exit 111
}

* Artificial firm-year data; no private or real research data are used.
set seed 1103
set obs 1000
generate long observation = _n
generate int firm = ceil(_n/10)
generate int year = 2000 + mod(_n - 1, 10)
generate double size = rnormal()
generate double cash = rnormal()
generate double leverage = rnormal()
generate double rd = .40*size + .20*cash + rnormal()
generate double patents = .20*rd + .30*size - .10*leverage + ///
    .15*sin(firm) + .05*(year - 2000) + rnormal()

local model patents rd, absorb(firm year) controls(size cash leverage) ///
    selection(bic, gridminok) vce(cluster firm)

* On this small example, auto uses official Stata.
quietly sparse_poreg `model' backend(auto)
scalar b_auto = _b[rd]
scalar se_auto = _se[rd]

* The same model using the package's fast large-fixed-effect route.
sort observation, stable
quietly sparse_poreg `model' backend(mata) strict(stata)
scalar b_mata = _b[rd]
scalar se_mata = _se[rd]

display as text _newline "Simple R&D-patent example"
display as text "---------------------------------------------"
display as text "Small-data official route: coefficient = " ///
    as result %8.4f b_auto as text ", SE = " as result %8.4f se_auto
display as text "Fast Mata route:           coefficient = " ///
    as result %8.4f b_mata as text ", SE = " as result %8.4f se_mata

assert abs(b_auto - b_mata) < 1e-8
assert abs(se_auto - se_mata) < 1e-8
display as result _newline "PASS: the official and fast Mata results agree."

```

## R\&D-Patent relation

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



This is an independent third-party package, not affiliated with or endorsed by StataCorp LLC. The package code is distributed under the MIT license.
