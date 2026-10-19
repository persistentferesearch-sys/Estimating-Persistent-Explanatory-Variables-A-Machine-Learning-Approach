# RD patent example: official Stata and sparse Mata

This page gives the **exact linear PRL and XPO commands** used for the public-data comparison. Start with the shorter [README](README.md) if you are installing the package for the first time.

The dataset is [Replicate_data_MissingRDPatent_ML.dta](../Replicate_data_MissingRDPatent_ML.dta) in this repository. The recorded comparison used Stata 19, 86,341 observations, 11,544 firm clusters, BIC with `gridminok`, and seed 1103. The data were sorted by `PERMCO fyear` before estimation. XPO used five observation folds. No source data were changed.

## 1. Load the data and define the variables

From the `sparse_poreg` folder in the downloaded repository:

```stata
use "../Replicate_data_MissingRDPatent_ML.dta", clear
sort PERMCO fyear
local y lnnpatent
local x RDAT lnME RD_missing lnAge lnK2L TobinQ ROA Leverage ///
    CASHAT KZidx InstOwn oms_HHidx oms_HHidx_square
```

Every variable in `x` is a **target with its own reported coefficient**. Firm and year indicators are selection controls. This reproduces the paper's multi-target design; moving the other 12 variables into `controls()` would estimate a different specification.

## 2. Linear PRL: official and Mata commands

The paper's direct official Stata command is:

```stata
poregress `y' `x', controls(i.PERMCO i.fyear) ///
    vce(cluster PERMCO) selection(bic, gridminok) rseed(1103)
```

The same official estimator can be called through this package:

```stata
sparse_poreg `y' `x', controls(i.PERMCO i.fyear) ///
    vce(cluster PERMCO) selection(bic, gridminok) rseed(1103) ///
    backend(official)
```

For the sparse Mata calculation, run:

```stata
sparse_poreg `y' `x', absorb(PERMCO fyear) ///
    vce(cluster PERMCO) selection(bic, gridminok) rseed(1103) ///
    backend(mata) strict(stata)
```

Do not add `frontierok` to this clustered BIC PRL command; this release rejects it here. `strict(stata)` requests the Mata comparison mode.

## 3. Linear XPO: official and Mata commands

The paper's direct official Stata command is:

```stata
xporegress `y' `x', controls(i.PERMCO i.fyear) ///
    vce(cluster PERMCO) xfolds(5) rseed(1103) ///
    selection(bic, gridminok)
```

The package's official route calls the same estimator:

```stata
sparse_xporeg `y' `x', controls(i.PERMCO i.fyear) ///
    vce(cluster PERMCO) xfolds(5) rseed(1103) ///
    selection(bic, gridminok) backend(official)
```

The recorded sparse Mata command is:

```stata
sparse_xporeg `y' `x', absorb(PERMCO fyear) ///
    vce(cluster PERMCO) xfolds(5) rseed(1103) ///
    selection(bic, gridminok) backend(mata) frontierok strict(stata)
```

`frontierok` is required for the experimental Mata XPO route; it does not disable safety checks. These five folds split observations. Firm-clustered standard errors do **not** mean that entire firms are held out from the training folds.

## 4. Compare the results

| Route | RDAT coefficient | Cluster standard error | Timing evidence |
|---|---:|---:|---:|
| [Published official PRL](../ReplicateLog_Results/ReplicateLog_PRL_Linear.smcl) | 0.1986367 | 0.0181304 | Published log window; not an estimator timer |
| Package official PRL | 0.198636660737012 | 0.0181303830385999 | 21,739.120 seconds, local Stata timer |
| Package Mata PRL | 0.198636660737147 | 0.0181303830386044 | 477.897 seconds, local Stata timer |
| [Published official XPO](../ReplicateLog_Results/ReplicateLog_DML_Linear.smcl) | 0.2131053 | 0.0136054 | Published log window; not an estimator timer |
| Package Mata XPO | 0.213105298267939 | 0.0136053698613053 | 2,107.911 seconds, local Stata timer |

The local **PRL** official and Mata saved estimates have the same coefficient names and order. Across all 13 coefficients, their maximum absolute difference is `1.35197408823728e-13`; across every entry of the 13-by-13 covariance matrix, it is `2.08904074594507e-15`. The Mata PRL fit was approximately 45.5 times faster than the local official fit on that machine. See the [saved-matrix audit](results/PRL_saved_matrix_parity_audit.log).

The completed **Mata XPO** fit matches all 13 coefficient and standard-error pairs **as printed** in the paper's official log. See the [pair-by-pair audit](results/XPO_published_parity_audit.json). The published log rounds the numbers, so this comparison does **not** establish equality of full-precision XPO coefficients or the complete covariance matrix. It also does not give a reliable official-versus-Mata XPO speed ratio.

The official results are also summarized in the paper's [PRL table](../ReplicateLog_Results/PRL_Linear_estimates.csv) and [XPO table](../ReplicateLog_Results/DML_Linear_estimates.csv). Completed package logs are [official PRL](results/PRL_package_official_complete.log) and [Mata PRL/XPO](results/PO_XPO_package_mata_complete.log). The reference logs saved with this package are [PRL](results/PRL_published_official_GitHub.smcl) and [XPO](results/XPO_published_official_GitHub.smcl).

## What this comparison establishes

The official commands use base-omitted `i.PERMCO i.fyear`; the Mata commands build all observed levels using `absorb(PERMCO fyear)`. That distinction can matter in another dataset. For a meaningful comparison, hold the sample, sort order, target list, selection rule, VCE, seed, and fold count fixed. The numerical statements above apply to this recorded public-data example, not every possible design.
