# Sparse partialing-out for empirical finance

This Stata package offers two **linear** estimators for studies with many possible controls or firm and year indicators:

- `sparse_poreg`: partialing-out lasso (PRL).
- `sparse_xporeg`: cross-fitted partialing-out lasso (XPO/DML).

Each command can run Stata's official estimator or the package's sparse Mata calculation. The examples below use R&D spending and patenting, but you can replace those variables with your own outcome, variable of interest, controls, and cluster ID.

## 1. Download and install

1. On [this repository's main page](https://github.com/persistentferesearch-sys/Estimating-Persistent-Explanatory-Variables-A-Machine-Learning-Approach), select **Code > Download ZIP** and extract it.
2. Open Stata 19. Choose **File > Change Working Directory** and select the extracted `sparse_poreg` folder. You should see `stata.toc` and `sparse_poreg.pkg` there.
3. Copy these lines into Stata:

```stata
net install sparse_poreg, from("`c(pwd)'") replace
which sparse_poreg
which sparse_xporeg
```

The two `which` commands should show installed `.ado` files. Installation does not change the paper's dataset.

## 2. Try the finance quick start

Run the included example:

```stata
do "tests/finance_quickstart.do"
```

It creates artificial firm-year data, estimates a simple R&D–patent relationship, and compares the official and sparse Mata PRL routes. A successful run ends with **`PASS: the official and fast Mata results agree.`** You do not need to download research data for this step. The [quick-start source](tests/finance_quickstart.do) shows exactly how the artificial data and both estimates are created.

## 3. Understand the command

For a model with one main variable of interest:

```stata
sparse_poreg lnnpatent RDAT, ///
    absorb(PERMCO fyear) ///
    controls(lnME lnAge TobinQ ROA Leverage) ///
    selection(bic, gridminok) vce(cluster PERMCO) ///
    backend(auto)
```

`lnnpatent` is the outcome. `RDAT` receives a reported coefficient and standard error. The financial variables in `controls()` are **candidates for selection**. `absorb(PERMCO fyear)` makes the observed firm and year indicator levels candidates too; it does not perform an ordinary within fixed-effect transformation. `vce(cluster PERMCO)` clusters standard errors by firm.

**Every variable between the outcome and the comma gets its own coefficient.** If you want only the RDAT coefficient, put your other explanatory variables in `controls()`, as in this illustration. This one-target illustration is different from the paper's published multi-target specification below.

For cross-fitted XPO, use the same model and specify the folds:

```stata
sparse_xporeg lnnpatent RDAT, ///
    absorb(PERMCO fyear) ///
    controls(lnME lnAge TobinQ ROA Leverage) ///
    selection(bic, gridminok) vce(cluster PERMCO) ///
    xfolds(5) rseed(1103) backend(auto)
```

## 4. Choose a route

| Your goal | Use |
|---|---|
| Begin a routine model | `backend(auto)` |
| Run the official Stata estimator through this package | `backend(official)` |
| Compare Mata with official Stata on a manageable design | `backend(mata) strict(stata)` for PRL; add `frontierok` for XPO |
| Run a very large sparse PRL design with native safety checks | `backend(mata) strict(production)` |

`backend(auto)` uses official Stata for a manageable design. For XPO, it stops if official expansion is too large; choose sparse XPO explicitly with `backend(mata) frontierok`. `strict(stata)` is Mata's comparison mode. `strict(production)` uses a different, more scalable screening contract, so it is not a promise of identical estimates.

For the paper-style five-fold XPO comparison, folds split **observations**, even with firm-clustered standard errors. If the research design requires whole firms kept out of each training fold, read `help sparse_xporeg` about `foldcluster()`: that is a different analysis.

## 5. Run the paper's public RD patent example

The paper's [public data](../Replicate_data_MissingRDPatent_ML.dta) are in the parent folder of `sparse_poreg` in the downloaded repository. With the working directory still set to `sparse_poreg`, run:

```stata
use "../Replicate_data_MissingRDPatent_ML.dta", clear
sort PERMCO fyear
local y lnnpatent
local x RDAT lnME RD_missing lnAge lnK2L TobinQ ROA Leverage ///
    CASHAT KZidx InstOwn oms_HHidx oms_HHidx_square

* Linear PRL using sparse Mata.
sparse_poreg `y' `x', absorb(PERMCO fyear) ///
    vce(cluster PERMCO) selection(bic, gridminok) rseed(1103) ///
    backend(mata) strict(stata)

* Cross-fitted linear XPO using five folds.
sparse_xporeg `y' `x', absorb(PERMCO fyear) ///
    vce(cluster PERMCO) selection(bic, gridminok) ///
    xfolds(5) rseed(1103) backend(mata) frontierok strict(stata)
```

Here, **all 13 variables in `x` are targets**; each receives a coefficient. Firm and year indicators are candidate controls. These are substantial fits, especially XPO, so try the artificial quick start first. The data loaded with `use` are read from the repository; the commands above do not save changes to that file.

The published official and package Mata results for RDAT are:

| Linear model | Published official coefficient (SE) | Package Mata coefficient (SE) |
|---|---:|---:|
| PRL | 0.1986367 (0.0181304) | 0.198636660737147 (0.0181303830386044) |
| XPO | 0.2131053 (0.0136054) | 0.213105298267939 (0.0136053698613053) |

The PRL result also has a local full-precision official-versus-Mata check. The XPO comparison matches the published values **at their printed precision**. See [RDPatentComparison.md](RDPatentComparison.md) for the exact official commands, all comparison details, and the saved logs. Matching this example does not guarantee identical results after changing the data, design, variable order, or selection rule.

## More help

In Stata, type `help sparse_poreg` or `help sparse_xporeg`. The package is independent third-party software, not affiliated with or endorsed by StataCorp LLC. Its code is distributed under the [MIT license](LICENSE).
