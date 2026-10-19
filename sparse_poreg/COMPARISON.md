# Public-data linear PRL/XPO comparison

Source: the [paper's public repository](https://github.com/persistentferesearch-sys/persistent_fe_remendy), its `REPLIC~1.DO` and `REPLIC~1.DTA`, and the published linear PRL/XPO logs. Local package fits used a byte-identical data copy, sorted by `PERMCO fyear`, with Stata 19, 86,341 observations, 11,544 firm clusters, BIC with `gridminok`, and seed 1103. XPO used five observation folds. No source data were changed.

After `use "C:/path/to/REPLIC~1.DTA", clear` and `sort PERMCO fyear`, define:

```stata
local y lnnpatent
local x RDAT lnME RD_missing lnAge lnK2L TobinQ ROA Leverage ///
    CASHAT KZidx InstOwn oms_HHidx oms_HHidx_square
```

| Route | Stata command | RDAT coefficient | Cluster SE | Timing evidence |
|---|---|---:|---:|---:|
| [Published official PRL](../ReplicateLog_Results/REPLIC~4.SMC) ([table](../ReplicateLog_Results/PRL_LI~1.CSV)) | ``poregress `y' `x', controls(i.PERMCO i.fyear) vce(cluster PERMCO) selection(bic, gridminok) rseed(1103)`` | 0.1986367 | 0.0181304 | 4,201 s (1 h 10 m 1 s; log window) |
| Package official PRL | ``sparse_poreg `y' `x', controls(i.PERMCO i.fyear) vce(cluster PERMCO) selection(bic, gridminok) rseed(1103) backend(official)`` | 0.198636660737012 | 0.0181303830385999 | 21,739.120 s (6.04 h; local Stata timer) |
| Package Mata PRL | ``sparse_poreg `y' `x', absorb(PERMCO fyear) vce(cluster PERMCO) selection(bic, gridminok) rseed(1103) backend(mata) strict(stata)`` | 0.198636660737147 | 0.0181303830386044 | 477.897 s (7.97 min; local Stata timer) |
| [Published official DML/XPO](../ReplicateLog_Results/REPLIC~1.SMC) ([table](../ReplicateLog_Results/DML_LI~1.CSV)) | ``xporegress `y' `x', controls(i.PERMCO i.fyear) vce(cluster PERMCO) xfolds(5) rseed(1103) selection(bic, gridminok)`` | 0.2131053 | 0.0136054 | 12,986 s (3 h 36 m 26 s; log window) |
| Package Mata DML/XPO | ``sparse_xporeg `y' `x', absorb(PERMCO fyear) vce(cluster PERMCO) xfolds(5) rseed(1103) selection(bic, gridminok) backend(mata) frontierok strict(stata)`` | 0.213105298267939 | 0.0136053698613053 | 2,107.911 s (35.13 min; local Stata timer) |

The published PRL log opened at 09:46:49 and closed at 10:56:50 on 7 February 2025; the published XPO log opened at 10:56:50 and closed at 14:33:16. The resulting 4,201-second and 12,986-second intervals are **log windows**, not estimation timers. They include any commands around estimation and possible idle time. The local package figures are Stata timer outputs from a different machine, so these columns should not be treated as a cross-machine speed benchmark.

With the data loaded and sorted and the y/x locals defined above, each table command is complete and copyable. All x variables are variables of interest; firm/year indicators are selection controls. The official backend delegates to Stata's poregress/xporegress. In clustered BIC PRL, frontierok is rejected by rc14; the completed PRL run did not use it. XPO requires that flag.

**PRL parity:** the local official and Mata saved estimates have identical coefficient names/order. Across all 13 coefficients, maximum absolute difference is `1.35197408823728e-13`; across all entries of the 13-by-13 covariance matrix, `2.08904074594507e-15`. Both pass a `1e-10` threshold. The published official PRL result matches the displayed values. Mata PRL was 45.5 times faster on this machine.

**XPO parity:** the complete Mata run matches all 13 coefficient/SE pairs in the paper's official DML log exactly **as printed**. N, cluster count, 11,570 controls, 1,737 selected controls, and five folds also match. `results/XPO_published_parity_audit.json` records every pair. The published log rounds results, so this does not establish equality of full-precision estimates or the covariance matrix. No official-versus-Mata XPO speed ratio is available from these logs. `frontierok` remains required for Mata XPO.

`absorb(PERMCO fyear)` builds all observed indicator levels, whereas the paper specifies base-omitted `i.PERMCO i.fyear`. That distinction can matter on other samples; `strict(stata)` alone does not guarantee parity for arbitrary data. The claims above are confined to this dataset and the recorded checks.

## Included evidence and source binding

- `results/PRL_package_official_complete.log`: completed anonymized package official-backend PRL log.
- `results/PRL_saved_matrix_parity_audit.log`: saved-estimate coefficient-name, coefficient, and full covariance parity check for local official versus Mata PRL.
- `results/PO_XPO_package_mata_complete.log`: complete anonymized package Mata PRL and XPO log.
- `results/PRL_published_official_GitHub.smcl`, `results/XPO_published_official_GitHub.smcl`: the paper's official Stata logs, with only their original log-directory string redacted. These are **published reference** results. The adjacent CSV files are the paper's rounded tables.
- `results/XPO_published_parity_audit.json`: each published XPO coefficient/SE versus the complete Mata run, at displayed precision.

Current GitHub Git blob IDs verified at packaging: do-file `9376e2a83ab1dfd80c906308e0d54a1e910d747a`; data `c5bac67e3d494622d8220b0d376649026b27dbbe`; official PRL log `2040b8e9e7aec649b43c96fe5b8395d9d50b287d` (`REPLIC~4.SMC`); official XPO log `305714f5a460ac3d84b332647b0030be01b3588b` (`REPLIC~1.SMC`). The local source do-file/data copies matched those blobs. This comparison concerns **linear** models only; the paper's Poisson commands and the separate sparse Poisson development command are outside this ZIP.
