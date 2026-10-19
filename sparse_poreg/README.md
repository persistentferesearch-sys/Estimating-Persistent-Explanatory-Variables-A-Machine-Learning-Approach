# sparse_poreg 0.8.4-rc14: linear PRL and XPO

This independent Stata package provides `sparse_poreg` (linear post-regularization/PRL) and `sparse_xporeg` (linear cross-fit partialing-out/DML). Both have official Stata and sparse Mata backends. This ZIP contains the unchanged 22 July 2026 `0.8.4-rc14` program files plus a public-data example, anonymized logs, and comparison notes. It does **not** include research data, saved estimates, or the separate developmental Poisson command. It is not a StataCorp product.

## Install

Extract the ZIP. In Stata, change directory to the extracted folder containing `sparse_poreg.pkg` and run:

```stata
net install sparse_poreg, from("`c(pwd)'") replace
which sparse_poreg
which sparse_xporeg
do "tests/finance_quickstart.do"
```

The commands require Stata 16 or newer; the empirical logs here used Stata 19. The quick check uses generated data and Stata 19. The optional Python backend has separate requirements in `python/requirements-experimental.txt`.

## Paper's public linear example: sparse Mata first

Download `REPLIC~1.DTA` from the [paper replication repository](https://github.com/persistentferesearch-sys/persistent_fe_remendy), then run:

```stata
use "C:/path/to/REPLIC~1.DTA", clear
sort PERMCO fyear
local y lnnpatent
local x RDAT lnME RD_missing lnAge lnK2L TobinQ ROA Leverage ///
    CASHAT KZidx InstOwn oms_HHidx oms_HHidx_square

* Sparse Mata PRL/PO, parity-oriented mode.
sparse_poreg `y' `x', absorb(PERMCO fyear) ///
    vce(cluster PERMCO) selection(bic, gridminok) rseed(1103) ///
    backend(mata) strict(stata)

* Sparse Mata DML/XPO, experimental frontier.
sparse_xporeg `y' `x', absorb(PERMCO fyear) ///
    vce(cluster PERMCO) selection(bic, gridminok) rseed(1103) ///
    xfolds(5) backend(mata) frontierok strict(stata)
```

The `x` variables are variables of interest. Firm/year indicators are selection controls. The copy-ready [Mata do-file](examples/github_rdpatent_mata_prl_xpo.do) runs both commands. For clustered BIC PRL, this release rejects `frontierok`; the successful PRL run used `backend(mata) strict(stata)`. The XPO route requires `frontierok`. That flag is an experimental-route acknowledgement, not an accuracy setting.

To run the exact official Stata routes through this package, use `controls(i.PERMCO i.fyear)` with `backend(official)` in place of `absorb(...)` and `backend(mata) ...`. The earlier [official comparison do-file](examples/github_rdpatent_linear_po_xpo.do) includes those commands. Direct official Stata commands are `poregress` and `xporegress` with the same model options. Full official fits can take hours.

[COMPARISON.md](COMPARISON.md) gives results, timing, command-by-command evidence, and parity limits. `results/` contains complete package official-PRL and Mata-PRL/XPO logs, the paper's published official PRL/XPO logs and CSVs, and the XPO printed-precision audit. The included official XPO log is the **paper's published reference**, not a package-produced run. The Mata XPO route still requires `frontierok`.

See `docs/README.md` for backend boundaries and `UPLOAD.md` for self-upload instructions. `SHA256SUMS` covers the ZIP contents; the adjacent `.zip.sha256` checks the whole archive. The package source is MIT licensed.
