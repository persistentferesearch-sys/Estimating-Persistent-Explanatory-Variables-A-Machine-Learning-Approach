version 19.0
clear all
set more off
args pkgarg

* Copy-ready public examples on a fold-safe crossed design.
* Usage: do validate_public_examples.do [pkgroot]
* With no root argument,
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

di as result "PASS: public official and Mata examples"
