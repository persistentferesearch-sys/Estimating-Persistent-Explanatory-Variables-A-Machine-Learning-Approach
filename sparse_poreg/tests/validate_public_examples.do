version 19.0
clear all
set more off
args pythonexe pythonrequired pkgarg

* Copy-ready public examples on a fold-safe crossed design.
* Usage: do validate_public_examples.do [pythonexe] [required] [pkgroot]
* Optional third argument is the source-tree root.  With no root argument,
* discover it from the working directory or the installed ado path.
local pkg = subinstr(trim(`"`pkgarg'"'), "\", "/", .)
local here = subinstr(`"`c(pwd)'"', "\", "/", .)
local adodir
if trim(`"`pkg'"') == "" {
    capture confirm file "`here'/stata/sparse_poreg.ado"
    if !_rc local pkg "`here'"
}
if trim(`"`pkg'"') == "" {
    capture confirm file "`here'/../stata/sparse_poreg.ado"
    if !_rc local pkg "`here'/.."
}
if trim(`"`pkg'"') == "" {
    capture quietly findfile sparse_poreg.ado
    if !_rc {
        local adofile = subinstr(`"`r(fn)'"', "\", "/", .)
        mata: st_local("adodir", pathgetparent(st_local("adofile")))
        if regexm(`"`adofile'"', "/stata/sparse_poreg[.]ado$") {
            local pkg = regexr(`"`adofile'"', "/stata/sparse_poreg[.]ado$", "")
        }
    }
}
if trim(`"`adodir'"') == "" & trim(`"`pkg'"') != "" local adodir "`pkg'/stata"
if trim(`"`adodir'"') == "" {
    di as err "cannot locate sparse_poreg.ado; pass the source root as the third argument"
    exit 601
}
adopath ++ "`adodir'"
local K 3
local seed 17072026
local paritytol 1e-10
scalar __pub_bdiff_max = 0
scalar __pub_vdiff_max = 0

capture program drop _pub_finite
program _pub_finite
    tempname b V
    matrix `b' = e(b)
    matrix `V' = e(V)
    mata: st_numscalar("__pub_fin", all(abs(st_matrix("`b'")) :< .) & all(abs(st_matrix("`V'")) :< .))
    mata: st_numscalar("__pub_sym", max(abs(st_matrix("`V'") :- st_matrix("`V'")')))
    mata: st_numscalar("__pub_vmin", min(symeigenvalues(st_matrix("`V'"))))
    assert scalar(__pub_fin) == 1
    assert scalar(__pub_sym) < 1e-10
    assert scalar(__pub_vmin) > -1e-12
end

capture program drop _pub_native
program _pub_native
    syntax, CMD(string) [PLUGIN XPO]
    assert "`e(cmd)'" == "`cmd'"
    assert "`e(sparse_backend)'" == "mata"
    assert "`e(sparse_strict)'" == "production"
    assert e(converged) == 1
    assert e(kkt_pass) == 1
    assert e(post_rank_valid) == 1
    assert e(finite) == 1
    assert e(sparse_route_n) == e(N)
    assert e(sparse_route_p) == e(sparse_features_raw)
    assert reldif(e(sparse_route_pn), ///
        e(sparse_route_p)/e(sparse_route_n)) < 1e-12
    assert "`e(sparse_backend_requested)'" == "mata"
    assert "`e(sparse_route_reason)'" == "explicit_mata"
    if "`plugin'" != "" assert e(plugin_pass) == 1
    if "`xpo'" != "" {
        assert e(training_support_safe) == 1
        assert e(training_support_relaxed) == 0
        assert e(min_train_support) > 1
    }
    _pub_finite
end

capture program drop _pub_python
program _pub_python
    syntax, CMD(string) [CLUSTERDIMS(integer 0)]
    tempname b V Vcopy offdiag
    local bnames
    local vrows
    local vcols
    matrix `b' = e(b)
    matrix `V' = e(V)
    matrix `Vcopy' = e(sparse_covariance)
    local bnames : colnames `b'
    local vrows : rownames `V'
    local vcols : colnames `V'
    scalar `offdiag' = abs(`V'[1,2])
    assert "`e(cmd)'" == "`cmd'"
    assert "`e(depvar)'" == "y"
    assert "`e(selection)'" == "lassocv"
    assert "`e(learner)'" == "lassocv"
    assert "`e(varsofinterest)'" == "`bnames'"
    assert "`vrows'" == "`bnames'" & "`vcols'" == "`bnames'"
    assert "`e(sparse_backend)'" == "python"
    assert "`e(sparse_backend_requested)'" == "python"
    assert "`e(sparse_route_reason)'" == "explicit_python"
    assert e(sparse_route_p) == e(sparse_features)
    assert e(sparse_route_n) == e(N)
    assert reldif(e(sparse_route_pn), ///
        e(sparse_route_p)/e(sparse_route_n)) < 1e-12
    assert trim(`"`e(python_exec)'"') != ""
    assert trim(`"`e(python_script)'"') != ""
    assert `"`e(python_executable)'"' == `"`e(python_exec)'"'
    capture confirm file `"`e(python_exec)'"'
    assert _rc == 0
    capture confirm file `"`e(python_script)'"'
    assert _rc == 0
    assert colsof(`b') == 2
    assert rowsof(`V') == 2 & colsof(`V') == 2
    assert `offdiag' > 1e-14
    assert mreldif(`V', `Vcopy') == 0
    assert e(sparse_covariance_full) == 1
    assert e(sparse_vcov_full) == 1
    assert e(n_treatments) == 2
    assert e(k_varsofinterest) == 2
    assert e(rank) == 2 & e(df_m) == 2
    assert e(post_rank_valid) == 1 & e(finite) == 1
    assert e(chi2) >= 0 & e(chi2) < .
    assert inrange(e(p), 0, 1)
    quietly count if e(sample)
    assert r(N) == e(N)
    quietly count if e(sample) & missing(x2)
    assert r(N) == 0
    assert e(N) < _N
    if `clusterdims' {
        assert "`e(inference_contract)'" == "cluster_dml"
        assert "`e(crossfit_unit)'" == "cluster"
        assert "`e(foldclustvars)'" == "`e(clustvars)'"
        assert e(cluster_disjoint) == 1
        assert e(n_cluster_dims) == `clusterdims'
        assert e(n_xfolds) == e(n_xfolds_per_cluster)
        assert e(n_score_folds) == e(n_xfolds)^`clusterdims'
        assert e(cluster_coverage_min) == 1
        assert e(cluster_coverage_max) == 1
        assert e(cluster_min_train) > 0
        assert e(cluster_min_test) > 0
    }
    _pub_finite
    assert abs(e(sparse_vcov_min_eigen) - scalar(__pub_vmin)) < ///
        1e-10*max(1, abs(scalar(__pub_vmin)))
end

* Every absorbed level occurs in every (g1,g2) cell, so all XPO training
* samples retain safe support under observation- and cluster-level folds.
local G1 12
local G2 9
local R 9
set obs `=`G1'*`G2'*`R''
gen long g1 = floor((_n - 1)/(`G2'*`R')) + 1
gen long g2 = floor(mod(_n - 1, `G2'*`R')/`R') + 1
gen byte rep = mod(_n - 1, `R') + 1
gen long a = rep
gen long b = 1 + mod(rep - 1 + g1 + 2*g2, `R')

set seed `seed'
gen double x1 = rnormal()
gen double x2 = rnormal()
gen double u = .30*sin(g1) + .20*cos(g2) + rnormal()
gen double d1 = .60*x1 - .20*x2 + .12*sin(a) + .08*cos(b) + rnormal()
gen double d2 = -.25*x1 + .50*x2 + .10*cos(a) - .07*sin(b) + rnormal()
gen double y = 1.10*d1 - .70*d2 + .40*x1 - .15*x2 + u

* Official and native PO/XPO: plugin/BIC x robust/one-way cluster.
foreach model in po xpo {
    local official `model'regress
    local sparse sparse_`model'reg
    local xopts
    local frontier
    local xflag
    if "`model'" == "xpo" {
        local xopts xfolds(`K') rseed(`seed')
        local frontier frontierok
        local xflag xpo
    }

    foreach sel in plugin bic {
        local selspec `sel'
        local pflag
        if "`sel'" == "plugin" local pflag plugin
        if "`sel'" == "bic" local selspec "bic, gridminok"

        foreach vc in robust cluster {
            local vcespec robust
            if "`vc'" == "cluster" local vcespec "cluster g1"

            tempname b_off V_off b_sp V_sp
            quietly `official' y d1 d2, controls(ibn.a ibn.b x1 x2) ///
                selection(`selspec') vce(`vcespec') `xopts'
            assert "`e(cmd)'" == "`official'"
            assert e(N) == _N
            matrix `b_off' = e(b)
            matrix `V_off' = e(V)
            _pub_finite

            quietly `sparse' y d1 d2, absorb(a b) controls(x1 x2) ///
                selection(`selspec') vce(`vcespec') backend(mata) ///
                strict(production) kernel(indexed) `frontier' `xopts'
            _pub_native, cmd(`sparse') `pflag' `xflag'
            matrix `b_sp' = e(b)
            matrix `V_sp' = e(V)
            mata: st_numscalar("__pub_bdiff", max(abs(st_matrix("`b_sp'") :- st_matrix("`b_off'"))))
            mata: st_numscalar("__pub_vdiff", max(abs(st_matrix("`V_sp'") :- st_matrix("`V_off'"))))
            assert scalar(__pub_bdiff) < `paritytol'
            assert scalar(__pub_vdiff) < `paritytol'
            scalar __pub_bdiff_max = max(scalar(__pub_bdiff_max), scalar(__pub_bdiff))
            scalar __pub_vdiff_max = max(scalar(__pub_vdiff_max), scalar(__pub_vdiff))
        }
    }
}
di as result "PASS: official/native 2x2 grids; max b/V difference = " ///
    %10.3e scalar(__pub_bdiff_max) " / " %10.3e scalar(__pub_vdiff_max)

* Native multiway examples: plugin and BIC use the same inference contracts.
foreach sel in plugin bic {
    local selspec `sel'
    local pflag
    if "`sel'" == "plugin" local pflag plugin
    if "`sel'" == "bic" local selspec "bic, gridminok"

    * Conventional two-way clustered PO inference.
    quietly sparse_poreg y d1 d2, absorb(a b) controls(x1 x2) ///
        selection(`selspec') vce(cluster g1 g2) ///
        backend(mata) frontierok strict(production) kernel(indexed)
    _pub_native, cmd(sparse_poreg) `pflag'
    assert "`e(inference_contract)'" == "multiway_cluster"
    assert e(n_cluster_dims) == 2
    assert e(N_clust) == min(`G1', `G2')

    * One-way cluster-disjoint XPO.
    quietly sparse_xporeg y d1 d2, absorb(a b) controls(x1 x2) ///
        selection(`selspec') vce(cluster g1) foldcluster(g1) ///
        xfolds(`K') rseed(`seed') backend(mata) frontierok ///
        strict(production) kernel(indexed)
    _pub_native, cmd(sparse_xporeg) `pflag' xpo
    assert "`e(inference_contract)'" == "cluster_dml"
    assert e(cluster_disjoint) == 1 & e(n_cluster_dims) == 1
    assert e(n_score_folds) == `K'

    * Two-way cluster-DML: K partitions per dimension give K^2 score folds.
    quietly sparse_xporeg y d1 d2, absorb(a b) controls(x1 x2) ///
        selection(`selspec') vce(cluster g1 g2) ///
        foldcluster(g1 g2) xfolds(`K') rseed(`seed') ///
        backend(mata) frontierok strict(production) kernel(indexed)
    _pub_native, cmd(sparse_xporeg) `pflag' xpo
    assert "`e(inference_contract)'" == "cluster_dml"
    assert e(cluster_disjoint) == 1 & e(n_cluster_dims) == 2
    assert e(n_score_folds) == `K'^2
}

* Optional Python/DoubleML examples.  Two treatments make the full covariance
* contract observable; missing controls verify that e(sample) is restored.
local py `"`pythonexe'"'
if trim(`"`py'"') == "" local py `"`c(python_exec)'"'
local pyscript
capture quietly findfile sparse_prl_xpo_backend.py
if !_rc local pyscript `"`r(fn)'"'
else if trim(`"`pkg'"') != "" local pyscript "`pkg'/python/sparse_prl_xpo_backend.py"
local pyready 1
capture confirm file `"`py'"'
if _rc local pyready 0
capture confirm file "`pyscript'"
if _rc local pyready 0

if !`pyready' & "`pythonrequired'" == "required" {
    di as err "required Python executable/backend script is unavailable"
    exit 601
}
if `pyready' {
    preserve
    gen double d2p = d2 + .35*d1
    replace x2 = . if mod(_n, 101) == 0
    tempname py_po_b
    local pyopts selection(lassocv) learner(lassocv) ///
        lassocvfolds(2) lassocvalphas(12) lassocvmaxiter(2000) ///
        backend(python) python(`"`py'"')

    quietly sparse_poreg y d1 d2p, absorb(a b) controls(x1 x2) ///
        vce(robust) `pyopts'
    _pub_python, cmd(sparse_poreg)
    matrix `py_po_b' = e(b)

    quietly sparse_poreg y d1 d2p, absorb(a b) controls(x1 x2) ///
        vce(cluster g1) rseed(`seed') `pyopts'
    _pub_python, cmd(sparse_poreg)
    assert mreldif(e(b), `py_po_b') == 0
    assert "`e(clustvars)'" == "g1"
    assert "`e(vcetype)'" == "Cluster-robust"
    assert "`e(inference_contract)'" == "conventional_cluster_robust"
    assert "`e(crossfit_unit)'" == "none"
    assert e(nuisance_crossfit) == 0
    local esc : e(scalars)
    assert !strpos(" `esc' ", " n_score_folds ")

    quietly sparse_poreg y d1 d2p, absorb(a b) controls(x1 x2) ///
        vce(cluster g1 g2) rseed(`seed') `pyopts'
    _pub_python, cmd(sparse_poreg)
    assert mreldif(e(b), `py_po_b') == 0
    assert "`e(clustvar)'" == "g1" & "`e(clustvars)'" == "g1 g2"
    assert "`e(vcetype)'" == "Multiway cluster-robust"
    assert "`e(inference_contract)'" == "conventional_multiway_cluster"
    assert "`e(crossfit_unit)'" == "none"
    assert e(nuisance_crossfit) == 0
    assert e(n_cluster_dims) == 2
    local esc : e(scalars)
    assert !strpos(" `esc' ", " n_score_folds ")
    assert e(N_clust1) == `G1' & e(N_clust2) == `G2'

    quietly sparse_xporeg y d1 d2p, absorb(a b) controls(x1 x2) ///
        vce(robust) xfolds(`K') rseed(`seed') `pyopts'
    _pub_python, cmd(sparse_xporeg)
    assert "`e(inference_contract)'" == "observation_dml"
    assert "`e(crossfit_unit)'" == "observation"
    assert e(nuisance_crossfit) == 1 & e(observation_disjoint) == 1
    assert e(n_xfolds) == `K' & e(n_score_folds) == `K'
    assert e(fold_coverage_min) == 1 & e(fold_coverage_max) == 1
    local esc : e(scalars)
    assert !strpos(" `esc' ", " cluster_disjoint ")

    * Matching foldcluster() is mandatory for clustered Python XPO.
    quietly sparse_xporeg y d1 d2p, absorb(a b) controls(x1 x2) ///
        vce(cluster g1) foldcluster(g1) ///
        xfolds(`K') rseed(`seed') `pyopts'
    _pub_python, cmd(sparse_xporeg) clusterdims(1)
    assert "`e(clustvars)'" == "g1"
    assert e(n_xfolds) == `K'
    assert e(N_clust) == `G1' & e(N_clust1) == `G1'

    quietly sparse_xporeg y d1 d2p, absorb(a b) controls(x1 x2) ///
        vce(cluster g1 g2) foldcluster(g1 g2) ///
        xfolds(`K') rseed(`seed') `pyopts'
    _pub_python, cmd(sparse_xporeg) clusterdims(2)
    assert "`e(clustvars)'" == "g1 g2"
    assert "`e(vcetype)'" == "Multiway cluster-robust"
    assert e(n_xfolds) == `K'
    assert e(N_clust) == min(`G1', `G2')
    assert e(N_clust1) == `G1' & e(N_clust2) == `G2'
    assert e(N_clust12) == `G1'*`G2'
    restore
    di as result "PASS: Python public examples"
}
else di as txt "SKIP: Python examples (configured executable/backend script unavailable)"

* Interface failures stop before Python is launched.
capture noisily sparse_xporeg y d1 d2, absorb(a b) controls(x1 x2) ///
    selection(lassocv) backend(python) vce(cluster g1) ///
    foldcluster(g2) xfolds(`K') python(`"`py'"')
assert _rc == 198
capture noisily sparse_xporeg y d1 d2, absorb(a b) controls(x1 x2) ///
    selection(lassocv) backend(python) vce(cluster g1 g1) ///
    foldcluster(g1) xfolds(`K') python(`"`py'"')
assert _rc == 198
capture noisily sparse_xporeg y d1 d2, absorb(a b) controls(x1 x2) ///
    selection(lassocv) backend(python) vce(cluster g1 g2 a) ///
    foldcluster(g1 g2) xfolds(`K') python(`"`py'"')
assert _rc == 198
capture noisily sparse_xporeg y d1 d2, absorb(a b) controls(x1 x2) ///
    selection(lassocv) backend(python) vce(cluster g1 g2) ///
    xfolds(`K') python(`"`py'"')
assert _rc == 198
capture noisily sparse_poreg y i.a, absorb(b) controls(x1 x2) ///
    selection(lassocv) backend(python) python(`"`py'"')
assert _rc == 198
capture noisily sparse_poreg y d1, absorb(a b) controls(c.x1#c.x2) ///
    selection(lassocv) backend(python) python(`"`py'"')
assert _rc == 198
capture noisily sparse_poreg y d1, absorb(a b) controls((x1) x2) ///
    selection(lassocv) backend(python) python(`"`py'"')
assert _rc == 198

di as result "PASS: public official, Mata, multiway, and optional Python examples"
