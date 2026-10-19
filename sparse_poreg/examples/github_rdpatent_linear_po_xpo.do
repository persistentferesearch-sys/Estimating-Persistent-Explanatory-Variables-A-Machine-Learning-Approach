*! Linear PO/XPO example for persistent_fe_remendy public replication data
version 19.0
clear all
set more off
set maxvar 120000

* Run from the extracted package folder after net install. The data are read
* from the path supplied by the user and are never changed on disk.
args datafile compare
if `"`datafile'"' == "" {
    display as error "Supply the path to the repository's REPLIC~1.DTA."
    exit 198
}
if "`compare'" == "" local compare 0
if !inlist("`compare'", "0", "1") {
    display as error "Second argument must be 0 or 1."
    exit 198
}
confirm file `"`datafile'"'
use `"`datafile'"', clear
sort PERMCO fyear

* The exact PO and DML runs match the paper's linear specifications. All
* variables after RDAT are variables of interest, not selection controls.
local outcome lnnpatent
local voi RDAT lnME RD_missing lnAge lnK2L TobinQ ROA Leverage ///
    CASHAT KZidx InstOwn oms_HHidx oms_HHidx_square
local controls i.PERMCO i.fyear
local common controls(`controls') vce(cluster PERMCO) ///
    selection(bic, gridminok) rseed(1103)

* Official poregress, run separately only when comparison is requested.
if `compare' {
    quietly poregress `outcome' `voi', `common'
    matrix official_po_b = e(b)
    matrix official_po_V = e(V)
    estimates store official_po
}

* 1. Exact PO: package wrapper for the official estimator and paper design.
sparse_poreg `outcome' `voi', `common' backend(official)
estimates store exact_po
display as result "PO: RDAT coefficient=" %12.8f _b[RDAT] ///
    "  clustered SE=" %12.8f _se[RDAT]
if `compare' {
    mata: assert(max(abs(st_matrix("official_po_b")-st_matrix("e(b)")))==0)
    mata: assert(max(abs(st_matrix("official_po_V")-st_matrix("e(V)")))==0)
    display as result "PO: package and direct official b/V match exactly."
}

* 2. Sparse Mata PO. absorb() uses all observed firm/year levels, unlike
* the paper's base-omitted i.PERMCO i.fyear; coefficients need not match.
* strict(stata) selects the package's official-parity audit behavior.
sparse_poreg `outcome' `voi', absorb(PERMCO fyear) ///
    vce(cluster PERMCO) selection(bic, gridminok) rseed(1103) ///
    backend(mata) strict(stata)
estimates store mata_po
display as result "Mata PO: RDAT coefficient=" %12.8f _b[RDAT] ///
    "  clustered SE=" %12.8f _se[RDAT]

* Official xporegress uses five observation folds in the paper.
if `compare' {
    quietly xporegress `outcome' `voi', `common' xfolds(5)
    matrix official_xpo_b = e(b)
    matrix official_xpo_V = e(V)
    estimates store official_xpo
}

* 3. DML/XPO: exact official five-fold paper specification.
sparse_xporeg `outcome' `voi', `common' xfolds(5) backend(official)
estimates store dml_xpo
display as result "XPO: RDAT coefficient=" %12.8f _b[RDAT] ///
    "  clustered SE=" %12.8f _se[RDAT]
if `compare' {
    mata: assert(max(abs(st_matrix("official_xpo_b")-st_matrix("e(b)")))==0)
    mata: assert(max(abs(st_matrix("official_xpo_V")-st_matrix("e(V)")))==0)
    display as result "XPO: package and direct official b/V match exactly."
}
display as result "GITHUB_LINEAR_PO_XPO_EXAMPLE_COMPLETE"
