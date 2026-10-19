*! Public-data Mata PRL and five-fold XPO; Stata 19; no data writes
version 19.0
clear all
set more off
set maxvar 120000

args datafile
if `"`datafile'"' == "" {
    display as error "Supply the path to the paper's REPLIC~1.DTA."
    exit 198
}
confirm file `"`datafile'"'
use `"`datafile'"', clear
sort PERMCO fyear

local y lnnpatent
local x RDAT lnME RD_missing lnAge lnK2L TobinQ ROA Leverage ///
    CASHAT KZidx InstOwn oms_HHidx oms_HHidx_square

* Clustered BIC PRL accepts strict(stata), but not frontierok in rc14.
sparse_poreg `y' `x', absorb(PERMCO fyear) ///
    vce(cluster PERMCO) selection(bic, gridminok) rseed(1103) ///
    backend(mata) strict(stata)
estimates store mata_prl
display as result "MATA_PRL RDAT=" %21.15g _b[RDAT] ///
    " SE=" %21.15g _se[RDAT]

* XPO's experimental Mata route requires frontierok.
sparse_xporeg `y' `x', absorb(PERMCO fyear) ///
    vce(cluster PERMCO) selection(bic, gridminok) rseed(1103) ///
    xfolds(5) backend(mata) frontierok strict(stata)
estimates store mata_xpo
display as result "MATA_XPO RDAT=" %21.15g _b[RDAT] ///
    " SE=" %21.15g _se[RDAT]
display as result "GITHUB_MATA_PRL_XPO_COMPLETE"
