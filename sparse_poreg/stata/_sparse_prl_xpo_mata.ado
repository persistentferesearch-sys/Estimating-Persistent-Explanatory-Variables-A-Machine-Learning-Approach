*! _sparse_prl_xpo_mata 0.8.4-rc14
program _sparse_prl_xpo_mata, eclass
    version 16.0
    syntax varlist(min=2 numeric) [if] [in], CMD(string) ABSorb(varlist numeric) ///
        [CONTROLS(varlist numeric) SELection(string) VCE(string asis) ///
        STRICT(string) BICMode(string) FOLDVAR(varname numeric) ///
        FOLDVAR2(varname numeric) ///
        XFOLDS(integer 10) ///
        RSEED(integer 1103) MAXIter(integer 4000) COEFTol(real 1e-3) ///
        ACCuracy(string) KERNEL(string) XCHECK XCHECKPrefix(name) ///
        ALLOWKKTFAIL CLUSTERFOLDS]

    marksample touse
    gettoken depvar treat : varlist
    local controls_all
    foreach a of local absorb {
        local controls_all `controls_all' ibn.`a'
    }
    local controls_all `controls_all' `controls'
    local k_varsofinterest : word count `treat'

    local cmd = lower(trim("`cmd'"))
    if !inlist("`cmd'", "sparse_poreg", "sparse_xporeg") {
        di as err "cmd() must be sparse_poreg or sparse_xporeg"
        exit 198
    }
    local selection = lower(trim("`selection'"))
    if "`selection'" == "" local selection plugin
    if !inlist("`selection'", "plugin", "bic") {
        di as err "selection() must be plugin or bic"
        exit 198
    }
    local bicmode = lower(trim("`bicmode'"))
    if "`bicmode'" == "" local bicmode = cond("`selection'" == "bic", "runningmin", "none")
    if ("`selection'" == "bic" & !inlist("`bicmode'", "runningmin", "gridminok")) | ///
        ("`selection'" == "plugin" & "`bicmode'" != "none") {
        di as err "bicmode() must be runningmin or gridminok with selection(bic)"
        exit 198
    }
    local strict = lower(trim("`strict'"))
    if "`strict'" == "" local strict stata
    if !inlist("`strict'", "stata", "production") {
        di as err "strict() must be stata or production"
        exit 198
    }
    local vcespec = trim("`vce'")
    if "`vcespec'" == "" local vcespec robust
    gettoken vceword vcerest : vcespec
    local vceword = lower(trim("`vceword'"))
    local vcerest = trim("`vcerest'")
    if "`vceword'" == "robust" & "`vcerest'" == "" {
        local vcetype robust
        local clustvar
        local clustvars
        local nclust 0
    }
    else if "`vceword'" == "cluster" {
        local clustvars `vcerest'
        local nclust : word count `clustvars'
        if !inrange(`nclust', 1, 2) {
            di as err "vce() must be robust or cluster clustvar [clustvar2]"
            exit 198
        }
        gettoken clustvar clustrest : clustvars
        gettoken clustvar2 : clustrest
        if `nclust' == 2 & "`clustvar'" == "`clustvar2'" {
            di as err "cluster variables must be distinct"
            exit 198
        }
        local vcetype cluster
    }
    else {
        di as err "vce() must be robust or cluster clustvar [clustvar2]"
        exit 198
    }
    if "`allowkktfail'" != "" & !("`cmd'" == "sparse_xporeg" & ///
        "`selection'" == "plugin" & "`vcetype'" == "cluster" & ///
        "`strict'" != "production") {
        di as err "allowkktfail is restricted to experimental XPO plugin/cluster fits"
        exit 198
    }
    if "`cmd'" == "sparse_xporeg" & "`foldvar'" == "" {
        di as err "foldvar() is required for cmd(sparse_xporeg)"
        exit 198
    }
    if "`foldvar2'" != "" & "`clusterfolds'" == "" {
        di as err "foldvar2() requires clusterfolds"
        exit 198
    }
    if "`clusterfolds'" != "" {
        if "`cmd'" != "sparse_xporeg" | "`vcetype'" != "cluster" {
            di as err "clusterfolds requires clustered sparse_xporeg"
            exit 198
        }
        if (`nclust' == 1 & "`foldvar2'" != "") | ///
            (`nclust' == 2 & "`foldvar2'" == "") {
            di as err "cluster and fold dimensions do not match"
            exit 198
        }
    }
    if "`cmd'" == "sparse_xporeg" & "`vcetype'" == "cluster" & ///
        `nclust' == 2 & "`clusterfolds'" == "" {
        di as err "two-way cluster inference requires clusterfolds"
        exit 198
    }
    if "`xcheckprefix'" != "" & ("`cmd'" != "sparse_xporeg" | ///
        "`foldvar2'" != "") {
        di as err "xcheckprefix() is restricted to ordinary sparse_xporeg"
        exit 198
    }
    local accuracy = lower(trim("`accuracy'"))
    if "`accuracy'" == "" local accuracy high
    if !inlist("`accuracy'", "high", "fast") {
        di as err "accuracy() must be high or fast"
        exit 198
    }
    local kernel = lower(trim("`kernel'"))
    if "`kernel'" == "" local kernel auto
    if !inlist("`kernel'", "auto", "legacy", "indexed") {
        di as err "kernel() must be auto, legacy, or indexed"
        exit 198
    }
    if `maxiter' < 1 {
        di as err "maxiter() must be positive"
        exit 198
    }

    local markvars `absorb' `controls'
    if "`clustvars'" != "" local markvars `markvars' `clustvars'
    if "`foldvar'" != "" local markvars `markvars' `foldvar'
    if "`foldvar2'" != "" local markvars `markvars' `foldvar2'
    markout `touse' `markvars'
    quietly count if `touse'
    if r(N) == 0 error 2000

    mata: _spmata_run("`cmd'", "`depvar'", "`treat'", "`absorb'", ///
        "`controls'", "`touse'", "`selection'", "`vcetype'", "`clustvars'", ///
        "`foldvar'", "`foldvar2'", `xfolds', `rseed', `maxiter', `coeftol', ///
        "`accuracy'", "`strict'", "`kernel'", "`clusterfolds'" != "", ///
        "`bicmode'", "`xcheck'" != "", "`xcheckprefix'", ///
        `c(matsize)', `c(macrolen)')

    if scalar(__spmata_controls_sel_complete) != 1 {
        di as txt "warning: selected-control union exceeds c(macrolen); e(controls_sel) omitted"
    }

    if scalar(__spmata_converged) != 1 {
        di as err "Mata lasso solver did not converge within maxiter(`maxiter')"
        di as err "maximum normalized KKT residual = " %10.3e scalar(__spmata_kkt_rel_max)
        exit 430
    }
    if scalar(__spmata_kkt_pass) != 1 {
        if "`allowkktfail'" == "" {
            di as err "Mata lasso solution failed the KKT certificate"
            di as err "maximum normalized KKT residual = " %10.3e scalar(__spmata_kkt_rel_max)
            exit 430
        }
        di as txt "warning: experimental XPO plugin/cluster fit is not KKT-certified"
        di as txt "maximum normalized KKT residual = " %10.3e scalar(__spmata_kkt_rel_max)
    }
    if "`selection'" == "plugin" & scalar(__spmata_plugin_pass) != 1 {
        if "`strict'" == "production" {
            di as err "plugin loading updates did not converge or stabilize support"
            di as err "capped nuisance calls = " %9.0g scalar(__spmata_plugin_capped_calls)
            exit 430
        }
        di as txt "warning: plugin loading updates did not converge or stabilize support"
        di as txt "capped nuisance calls = " %9.0g scalar(__spmata_plugin_capped_calls)
    }
    if "`cmd'" == "sparse_poreg" & scalar(__spmata_k_treat_omitted) > 0 {
        di as txt "note: `__spmata_treat_omitted' omitted because of collinearity"
    }
    if scalar(__spmata_rank_valid) != 1 {
        di as err "rank-deficient final treatment regression"
        exit 498
    }
    if scalar(__spmata_finite) != 1 {
        di as err "Mata backend produced nonfinite estimates"
        exit 498
    }
    if "`cmd'" == "sparse_poreg" & scalar(__spmata_vce_psd_adjusted) == 1 {
        di as txt "warning: two-way cluster VCE was projected to positive semidefinite"
    }

    matrix colnames __spmata_b = `treat'
    matrix rownames __spmata_V = `treat'
    matrix colnames __spmata_V = `treat'
    local post_dof
    if "`cmd'" == "sparse_poreg" & "`vcetype'" == "cluster" & `nclust' == 2 {
        local post_dof dof(`=scalar(__spmata_N_clust) - 1')
    }
    ereturn post __spmata_b __spmata_V, obs(`=scalar(__spmata_N)') `post_dof'
    local wald_df = scalar(__spmata_vce_rank)
    ereturn scalar rank = `wald_df'
    if "`cmd'" == "sparse_poreg" & "`vcetype'" == "cluster" & `nclust' == 2 {
        ereturn scalar df_m = `wald_df'
        if `wald_df' > 0 {
            ereturn scalar F = scalar(__spmata_wald) / `wald_df'
            ereturn scalar p = Ftail(`wald_df', e(df_r), e(F))
        }
        else {
            ereturn scalar F = .
            ereturn scalar p = .
        }
    }
    else {
        ereturn scalar df = `wald_df'
        ereturn local chi2type "Wald"
        if `wald_df' > 0 {
            ereturn scalar chi2 = scalar(__spmata_wald)
            ereturn scalar p = chi2tail(`wald_df', e(chi2))
        }
        else {
            ereturn scalar chi2 = .
            ereturn scalar p = .
        }
    }
    ereturn local cmd "`cmd'"
    ereturn hidden local estat_cmd "_sparse_estat"
    if "`cmd'" == "sparse_poreg" ///
        ereturn local title "Partialing-out linear model"
    else ereturn local title "Cross-fit partialing-out linear model"
    ereturn local depvar "`depvar'"
    ereturn local varsofinterest "`treat'"
    ereturn local varsofinterest_omitted "`__spmata_treat_omitted'"
    ereturn scalar k_varsofinterest_omitted = scalar(__spmata_k_treat_omitted)
    ereturn local selection "`selection'"
    if "`selection'" == "bic" ereturn local bic_mode "`bicmode'"
    ereturn local vce "`vcetype'"
    if "`vcetype'" == "robust" ereturn local vcetype "Robust"
    else ereturn local vcetype "Cluster-robust"
    ereturn local sparse_backend "mata"
    ereturn local sparse_strict "`strict'"
    ereturn local sparse_absorb "`absorb'"
    ereturn local sparse_accuracy "`accuracy'"
    if scalar(__spmata_indexed) ereturn local sparse_kernel "indexed"
    else ereturn local sparse_kernel "legacy"
    ereturn scalar N = scalar(__spmata_N)
    ereturn scalar sparse_features = scalar(__spmata_p)
    ereturn scalar sparse_features_raw = scalar(__spmata_raw_p)
    ereturn scalar xcheck_features_kept = scalar(__spmata_xcheck_kept)
    ereturn scalar xcheck_features_dropped = scalar(__spmata_xcheck_dropped)
    ereturn scalar xcheck_features_kept_min = scalar(__spmata_xcheck_kept_min)
    ereturn scalar xcheck_features_kept_max = scalar(__spmata_xcheck_kept_max)
    ereturn scalar xcheck_features_dropped_min = scalar(__spmata_xcheck_dropped_min)
    ereturn scalar xcheck_features_dropped_max = scalar(__spmata_xcheck_dropped_max)
    if "`xcheckprefix'" != "" {
        ereturn local sparse_xcheck "stata_fold"
        ereturn local sparse_xcheck_order "fold_input"
        ereturn local kkt_scope "stata_fold_screened_design"
    }
    else if "`xcheck'" != "" {
        ereturn local sparse_xcheck "stata"
        ereturn local sparse_xcheck_order "input"
        ereturn local kkt_scope "stata_screened_design"
    }
    else {
        ereturn local sparse_xcheck "none"
        ereturn local kkt_scope "full_design"
    }
    ereturn scalar converged = scalar(__spmata_converged)
    ereturn scalar lasso_calls = scalar(__spmata_lasso_calls)
    ereturn scalar lasso_iterations = scalar(__spmata_lasso_iterations)
    ereturn scalar lasso_sweeps = scalar(__spmata_lasso_sweeps)
    ereturn scalar cd_full_sweeps = scalar(__spmata_cd_full_sweeps)
    ereturn scalar cd_active_sweeps = scalar(__spmata_cd_active_sweeps)
    ereturn scalar cd_kkt_scans = scalar(__spmata_cd_kkt_scans)
    ereturn scalar cd_kkt_fallbacks = scalar(__spmata_cd_kkt_fallbacks)
    ereturn scalar cd_active_max = scalar(__spmata_cd_active_max)
    ereturn scalar cd_polish_calls = scalar(__spmata_cd_polish_calls)
    ereturn scalar cd_polish_success = scalar(__spmata_cd_polish_success)
    ereturn scalar bic_kkt_refine_calls = scalar(__spmata_bic_kkt_refine_calls)
    ereturn scalar bic_kkt_refine_success = scalar(__spmata_bic_kkt_refine_success)
    ereturn scalar kkt_max = scalar(__spmata_kkt_max)
    ereturn scalar kkt_rel_max = scalar(__spmata_kkt_rel_max)
    ereturn scalar kkt_pass = scalar(__spmata_kkt_pass)
    ereturn scalar plugin_converged = scalar(__spmata_plugin_converged)
    ereturn scalar plugin_pass = scalar(__spmata_plugin_pass)
    ereturn scalar plugin_support_stable = scalar(__spmata_plugin_support_stable)
    ereturn scalar plugin_capped_calls = scalar(__spmata_plugin_capped_calls)
    ereturn scalar plugin_load_diff_max = scalar(__spmata_plugin_load_diff_max)
    ereturn scalar plugin_iterations = scalar(__spmata_plugin_iterations)
    ereturn scalar selected_max = scalar(__spmata_selected_max)
    ereturn scalar selected_final = scalar(__spmata_selected_final)
    ereturn scalar k_controls = scalar(__spmata_raw_p)
    ereturn scalar k_controls_sel = scalar(__spmata_k_controls_sel)
    ereturn local controls_sel `"`__spmata_controls_sel'"'
    ereturn scalar controls_sel_complete = scalar(__spmata_controls_sel_complete)
    ereturn local controls "`controls_all'"
    local n_lasso_fits = scalar(__spmata_lasso_fits)
    local selection_counts_complete = (`n_lasso_fits' <= c(matsize))
    if `selection_counts_complete' {
        tempname selection_counts
        matrix `selection_counts' = J(`n_lasso_fits', 6, .)
    }
    local lasso_depvars
    local lasso_selected_complete 1
    forvalues i = 1/`n_lasso_fits' {
        foreach x in d s k f g t c {
            local key __sl`x'`i'
            local fit_`x' `"``key''"'
        }
        local lasso_depvars `lasso_depvars' `fit_d'
        ereturn hidden local lasso_depvar_`i' `"`fit_d'"'
        ereturn hidden local lasso_selected_`i' `"`fit_s'"'
        ereturn hidden scalar lasso_sel_p_`i' = `fit_k'
        ereturn hidden scalar lasso_fold_`i' = `fit_f'
        ereturn hidden scalar lasso_fold2_`i' = `fit_g'
        ereturn hidden scalar lasso_selected_complete_`i' = `fit_c'
        if `selection_counts_complete' {
            matrix `selection_counts'[`i', 1] = `i'
            matrix `selection_counts'[`i', 2] = `fit_f'
            matrix `selection_counts'[`i', 3] = `fit_g'
            matrix `selection_counts'[`i', 4] = `fit_t'
            matrix `selection_counts'[`i', 5] = `fit_k'
            matrix `selection_counts'[`i', 6] = `fit_c'
        }
        if !`fit_c' local lasso_selected_complete 0
    }
    if `selection_counts_complete' {
        matrix colnames `selection_counts' = fit fold fold2 target selected complete
        ereturn matrix selection_counts = `selection_counts'
    }
    ereturn scalar selection_counts_complete = `selection_counts_complete'
    ereturn scalar k_lasso_fits = `n_lasso_fits'
    ereturn scalar lasso_selected_complete = `lasso_selected_complete'
    ereturn local lasso_depvars "`lasso_depvars'"
    ereturn local lasso_depvars_uniq "`depvar' `treat'"
    ereturn scalar k_varsofinterest = `k_varsofinterest'
    if scalar(__spmata_support_matrix) {
        matrix colnames __spmata_support = feature absorb_index raw_level
        ereturn matrix sparse_support = __spmata_support
    }
    ereturn scalar sparse_support_matrix = scalar(__spmata_support_matrix)
    ereturn scalar bic_calls = scalar(__spmata_bic_calls)
    ereturn scalar bic_identified = scalar(__spmata_bic_identified)
    ereturn scalar bic_stop = scalar(__spmata_bic_stop)
    ereturn scalar bic_grid = scalar(__spmata_bic_grid)
    ereturn scalar bic_saturation_rescues = scalar(__spmata_bic_saturation_rescues)
    ereturn scalar bic_knots_total = scalar(__spmata_bic_knots_total)
    ereturn scalar bic_knots_max = scalar(__spmata_bic_knots_max)
    ereturn scalar post_rank_valid = scalar(__spmata_rank_valid)
    ereturn scalar nuisance_rank_deficient = scalar(__spmata_nuisance_rank_deficient)
    ereturn scalar post_condition = scalar(__spmata_post_condition)
    ereturn scalar final_rank = scalar(__spmata_final_rank)
    ereturn scalar final_condition = scalar(__spmata_final_condition)
    ereturn scalar finite = scalar(__spmata_finite)
    ereturn scalar vce_psd_adjusted = scalar(__spmata_vce_psd_adjusted)
    ereturn scalar vce_raw_min_eigen = scalar(__spmata_vce_raw_min_eigen)
    ereturn scalar vce_rank = scalar(__spmata_vce_rank)
    ereturn scalar n_xfolds = `xfolds'
    ereturn scalar n_score_folds = scalar(__spmata_score_folds)
    ereturn scalar n_cluster_dims = scalar(__spmata_cluster_ways)
    if "`clusterfolds'" != "" {
        ereturn local crossfit_unit "cluster"
        ereturn local foldclustvars "`clustvars'"
    }
    else if "`cmd'" == "sparse_xporeg" ereturn local crossfit_unit "observation"
    if "`clustvars'" != "" {
        ereturn local clustvar "`clustvar'"
        ereturn local clustvars "`clustvars'"
        ereturn scalar N_clust = scalar(__spmata_N_clust)
        ereturn scalar N_clust1 = scalar(__spmata_N_clust1)
        if scalar(__spmata_cluster_ways) == 2 {
            ereturn scalar N_clust2 = scalar(__spmata_N_clust2)
            ereturn scalar N_clust12 = scalar(__spmata_N_clust12)
        }
    }
end

capture mata: mata drop _spmata_group()
capture mata: mata drop _spmata_same_partition()
capture mata: mata drop _spmata_diag_init()
capture mata: mata drop _spmata_support_add()
capture mata: mata drop _spmata_selection_record()
capture mata: mata drop _spmata_support_names()
capture mata: mata drop _spmata_support_table()
capture mata: mata drop _spmata_apply_xcheck()
capture mata: mata drop _spmata_build()
capture mata: mata drop _spmata_build_kernel()
capture mata: mata drop _spmata_build_kernel_opts()
capture mata: mata drop _spmata_weights()
capture mata: mata drop _spmata_mean()
capture mata: mata drop _spmata_scale()
capture mata: mata drop _spmata_xb()
capture mata: mata drop _spmata_grad()
capture mata: mata drop _spmata_kkt()
capture mata: mata drop _spmata_lasso_status()
capture mata: mata drop _spmata_kkt_status()
capture mata: mata drop _spmata_lip()
capture mata: mata drop _spmata_soft()
capture mata: mata drop _spmata_keep()
capture mata: mata drop _spmata_keep_strict()
capture mata: mata drop _spmata_col()
capture mata: mata drop _spmata_localpos()
capture mata: mata drop _spmata_cd_cache()
capture mata: mata drop _spmata_lasso_cd()
capture mata: mata drop _spmata_lasso_cd_cached()
capture mata: mata drop _spmata_lasso_cd_cyclic_cached()
capture mata: mata drop _spmata_lasso_cd_screened_cached()
capture mata: mata drop _spmata_cd_active_block()
capture mata: mata drop _spmata_cd_polish()
capture mata: mata drop _spmata_cd_sweep()
capture mata: mata drop _spmata_cd_sweep_indexed()
capture mata: mata drop _spmata_cd_sweep_indexed_cached()
capture mata: mata drop _spmata_cd_sweep_factor_cached()
capture mata: mata drop _spmata_cd_sweep_absblock_cached()
capture mata: mata drop _spmata_lasso()
capture mata: mata drop _spmata_lasso_diag()
capture mata: mata drop _spmata_plugin_polish_diag()
capture mata: mata drop _spmata_loadings()
capture mata: mata drop _spmata_loadings_indexed()
capture mata: mata drop _spmata_plugin_lambda()
capture mata: mata drop _spmata_dense()
capture mata: mata drop _spmata_fitcheck()
capture mata: mata drop _spmata_independent_columns()
capture mata: mata drop _spmata_resid_ols()
capture mata: mata drop _spmata_resid_wols()
capture mata: mata drop _spmata_resid_wols_diag()
capture mata: mata drop _spmata_plugin_resid_diag()
capture mata: mata drop _spmata_predict_ols()
capture mata: mata drop _spmata_topcorr()
capture mata: mata drop _spmata_plugin_select_core()
capture mata: mata drop _spmata_plugin_select()
capture mata: mata drop _spmata_bic_kkt_refine()
capture mata: mata drop _spmata_bic_select()
capture mata: mata drop _spmata_bic_dump_path()
capture mata: mata drop _spmata_probe_run()
capture mata: mata drop _spmata_prl()
capture mata: mata drop _spmata_xpo()
capture mata: mata drop _spmata_xpo2()
capture mata: mata drop _spmata_sandwich()
capture mata: mata drop _spmata_cluster_vce()
capture mata: mata drop _spmata_sandwich2()
capture mata: mata drop _spmata_dml2()
capture mata: mata drop _spmata_dml_cluster1()
capture mata: mata drop _spmata_dml_cluster2()
capture mata: mata drop _spmata_run()

mata:
struct spgroup {
    real colvector code
    real colvector count
    real colvector level
    real scalar G
}

struct spdesign {
    real matrix acode
    real rowvector nlev
    real rowvector offset
    real rowvector avarpos
    real colvector astart
    real colvector acount
    real colvector aobs
    real colvector alevel
    real matrix C
    string rowvector avars
    string rowvector cvars
    string rowvector xalias
    real rowvector xalias_col
    real scalar p
    real scalar raw_p
    real scalar abs_p
    real scalar indexed
    real rowvector xkeep
}

struct spdiag {
    real scalar converged, lasso_calls, lasso_iterations, lasso_sweeps
    real scalar cd_full_sweeps, cd_active_sweeps, cd_kkt_scans
    real scalar cd_kkt_fallbacks, cd_active_max, cd_polish_calls
    real scalar cd_polish_success, bic_kkt_refine_calls, bic_kkt_refine_success
    real scalar kkt_max, kkt_rel_max, kkt_pass
    real scalar plugin_calls, plugin_converged, plugin_pass
    real scalar plugin_support_stable, plugin_iterations, plugin_capped_calls
    real scalar plugin_load_diff_max
    real scalar selected_max, selected_final
    real scalar bic_gridminok, bic_calls, bic_identified, bic_stop, bic_grid
    real scalar bic_knots_total, bic_knots_max, bic_saturation_rescues
    real scalar rank_valid, post_fits, post_condition, nuisance_rank_deficient
    real scalar final_rank, final_condition, finite
    real scalar vce_psd_adjusted, vce_raw_min_eigen
    real colvector selected_mask, omitted_treat
}

struct spcdcache {
    real colvector apos, densepos, wx
    real rowvector xval, xmean, hdiag, ordv, densemean, denseh
    real matrix densex
}

struct spdiag scalar _spmata_diag_init()
{
    struct spdiag scalar diag

    diag.converged = 1
    diag.plugin_converged = 1
    diag.plugin_pass = 1
    diag.plugin_support_stable = 1
    diag.rank_valid = 1
    diag.kkt_pass = 1
    diag.lasso_calls = 0
    diag.lasso_iterations = 0
    diag.lasso_sweeps = 0
    diag.cd_full_sweeps = 0
    diag.cd_active_sweeps = 0
    diag.cd_kkt_scans = 0
    diag.cd_kkt_fallbacks = 0
    diag.cd_active_max = 0
    diag.cd_polish_calls = 0
    diag.cd_polish_success = 0
    diag.bic_kkt_refine_calls = 0
    diag.bic_kkt_refine_success = 0
    diag.kkt_max = 0
    diag.kkt_rel_max = 0
    diag.plugin_calls = 0
    diag.plugin_iterations = 0
    diag.plugin_capped_calls = 0
    diag.plugin_load_diff_max = 0
    diag.selected_max = 0
    diag.selected_final = 0
    diag.bic_gridminok = 0
    diag.bic_calls = 0
    diag.bic_identified = 0
    diag.bic_stop = 0
    diag.bic_grid = 0
    diag.bic_knots_total = 0
    diag.bic_knots_max = 0
    diag.bic_saturation_rescues = 0
    diag.post_fits = 0
    diag.post_condition = 0
    diag.nuisance_rank_deficient = 0
    diag.final_rank = .
    diag.final_condition = .
    diag.finite = .
    diag.vce_psd_adjusted = 0
    diag.vce_raw_min_eigen = .
    diag.selected_mask = J(0, 1, .)
    diag.omitted_treat = J(0, 1, .)
    return(diag)
}

struct spgroup scalar _spmata_group(real colvector v)
{
    real scalar n, i, g, pos
    real colvector ord, code, count, level
    struct spgroup scalar out

    n = rows(v)
    ord = order(v, 1)
    code = J(n, 1, .)
    count = J(n, 1, 0)
    level = J(n, 1, .)
    g = 0
    for (i = 1; i <= n; i++) {
        pos = ord[i]
        if (i == 1) {
            g++
            level[g] = v[pos]
        }
        else if (v[pos] != v[ord[i - 1]]) {
            g++
            level[g] = v[pos]
        }
        code[pos] = g
        count[g] = count[g] + 1
    }
    out.code = code
    out.count = count[1::g]
    out.level = level[1::g]
    out.G = g
    return(out)
}

real scalar _spmata_same_partition(real colvector x, real colvector y,
    real scalar G)
{
    real scalar i
    real colvector map

    map = J(G, 1, 0)
    for (i = 1; i <= rows(x); i++) {
        if (!map[x[i]]) map[x[i]] = y[i]
        else if (map[x[i]] != y[i]) return(0)
    }
    return(1)
}

struct spdesign scalar _spmata_build_kernel_opts(real colvector idx,
    string scalar absorb, string scalar controls, string scalar kernel,
    real scalar aliases)
{
    real scalar a, c, i, j, k, n, p0, raw_abs_p, base, pos, col, duplicate
    real scalar alias_base
    string rowvector avars, cvars
    real matrix A
    real colvector cursor, alias_code
    struct spgroup scalar g
    struct spdesign scalar S

    avars = tokens(absorb)
    cvars = tokens(controls)
    n = rows(idx)
    a = cols(avars)
    c = cols(cvars)
    S.acode = J(n, 0, .)
    S.nlev = J(1, 0, .)
    S.avarpos = J(1, 0, .)
    S.alevel = J(0, 1, .)
    S.avars = J(1, 0, "")
    S.cvars = cvars
    S.xalias = J(1, 0, "")
    S.xalias_col = J(1, 0, .)
    raw_abs_p = 0
    if (a) {
        A = st_data(idx, avars)
        for (j = 1; j <= a; j++) {
            g = _spmata_group(A[, j])
            raw_abs_p = raw_abs_p + g.G
            duplicate = 0
            for (k = 1; k <= cols(S.acode); k++) {
                if (g.G == S.nlev[k]) {
                    duplicate = _spmata_same_partition(g.code, S.acode[, k], g.G)
                }
                if (duplicate) {
                    break
                }
            }
            if (!duplicate) {
                S.acode = S.acode, g.code
                S.nlev = S.nlev, g.G
                S.avarpos = S.avarpos, j
                S.alevel = S.alevel \ g.level
                S.avars = S.avars, avars[j]
            }
            else if (aliases) {
                // Stata may retain a collinear block through either variable
                // name.  Keep aliases for the discarded partition so its
                // e(allvars) mask maps onto the compact retained block.
                alias_code = J(g.G, 1, .)
                for (i = 1; i <= n; i++)
                    alias_code[g.code[i]] = S.acode[i, k]
                alias_base = 0
                for (i = 1; i < k; i++) alias_base = alias_base + S.nlev[i]
                for (i = 1; i <= g.G; i++) {
                    S.xalias = S.xalias,
                        strtrim(strofreal(g.level[i], "%21.0g")) +
                        "bn." + avars[j]
                    S.xalias_col = S.xalias_col,
                        alias_base + alias_code[i]
                }
            }
        }
        a = cols(S.acode)
        S.offset = J(1, a, 0)
        p0 = 0
        for (j = 1; j <= a; j++) {
            S.offset[j] = p0
            p0 = p0 + S.nlev[j]
        }
        S.abs_p = p0
        // One user rule: Mata uses indexed unless legacy is forced for tests.
        S.indexed = (kernel != "legacy")
        if (S.indexed) {
            // Inverted supports for absorbed indicators.  aobs stores the
            // original-row indices for every level, ordered by row within
            // level; astart/acount locate a feature's slice.  Storage is
            // O(n * #absorb), and CD visits only nonzero observations.
            S.astart = J(S.abs_p, 1, 0)
            S.acount = J(S.abs_p, 1, 0)
            S.aobs = J(n * a, 1, 0)
            base = 0
            for (j = 1; j <= a; j++) {
                for (i = 1; i <= n; i++) {
                    col = S.offset[j] + S.acode[i, j]
                    S.acount[col] = S.acount[col] + 1
                }
                pos = base + 1
                for (k = 1; k <= S.nlev[j]; k++) {
                    col = S.offset[j] + k
                    S.astart[col] = pos
                    pos = pos + S.acount[col]
                }
                cursor = S.astart[(S.offset[j] + 1)::
                    (S.offset[j] + S.nlev[j])]
                for (i = 1; i <= n; i++) {
                    k = S.acode[i, j]
                    S.aobs[cursor[k]] = i
                    cursor[k] = cursor[k] + 1
                }
                base = base + n
            }
        }
        else {
            S.astart = J(0, 1, .)
            S.acount = J(0, 1, .)
            S.aobs = J(0, 1, .)
        }
    }
    else {
        S.abs_p = 0
        S.indexed = 0
        S.offset = J(1, 0, .)
        S.astart = J(0, 1, .)
        S.acount = J(0, 1, .)
        S.aobs = J(0, 1, .)
    }
    S.C = c ? st_data(idx, cvars) : J(n, 0, .)
    S.p = S.abs_p + c
    S.raw_p = raw_abs_p + c
    S.xkeep = J(1, S.p, 1)
    return(S)
}

struct spdesign scalar _spmata_build_kernel(real colvector idx,
    string scalar absorb, string scalar controls, string scalar kernel)
{
    return(_spmata_build_kernel_opts(idx, absorb, controls, kernel, 0))
}

// Backward-compatible constructor used by diagnostic/probe do-files.
struct spdesign scalar _spmata_build(real colvector idx, string scalar absorb,
    string scalar controls)
{
    return(_spmata_build_kernel(idx, absorb, controls, "legacy"))
}

void _spmata_support_add(struct spdiag scalar diag, real colvector selected)
{
    if (rows(selected)) diag.selected_mask[selected] = J(rows(selected), 1, 1)
}

string rowvector _spmata_support_names(struct spdesign scalar S,
    real colvector selected)
{
    real scalar i, col, a
    string rowvector names

    names = J(1, rows(selected), "")
    for (i = 1; i <= rows(selected); i++) {
        col = selected[i]
        if (col > S.abs_p) names[i] = S.cvars[col - S.abs_p]
        else {
            a = 1
            while (col > S.offset[a] + S.nlev[a]) a++
            names[i] = strtrim(strofreal(S.alevel[col], "%21.0g")) +
                "bn." + S.avars[a]
        }
    }
    return(names)
}

void _spmata_selection_record(struct spdesign scalar S,
    real colvector selected, string scalar depvar, real scalar fit,
    real scalar fold, real scalar fold2, real scalar target,
    real scalar macrolen)
{
    real scalar complete
    string scalar i, support

    i = strtrim(strofreal(fit, "%21.0f"))
    support = invtokens(_spmata_support_names(S, selected))
    complete = strlen(support) <= macrolen
    st_local("__sld" + i, depvar)
    st_local("__sls" + i, complete ? support : "")
    st_local("__slk" + i, strtrim(strofreal(rows(selected), "%21.0f")))
    st_local("__slf" + i, strtrim(strofreal(fold, "%21.0g")))
    st_local("__slg" + i, strtrim(strofreal(fold2, "%21.0g")))
    st_local("__slt" + i, strtrim(strofreal(target, "%21.0f")))
    st_local("__slc" + i, strofreal(complete, "%1.0f"))
    st_numscalar("__spmata_lasso_fits", fit)
}

void _spmata_apply_xcheck(struct spdesign scalar S, string scalar spec)
{
    transmorphic scalar map
    real scalar i, col, kept, dropped, old
    string rowvector names, keep

    S.xkeep = J(1, S.p, 1)
    if (strtrim(spec) != "") {
        names = _spmata_support_names(S, 1::S.p)
        keep = tokens(spec)
        map = asarray_create("string", 1)
        for (i = 1; i <= cols(names); i++) asarray(map, names[i], i)
        for (i = 1; i <= cols(S.xalias); i++)
            asarray(map, S.xalias[i], S.xalias_col[i])
        S.xkeep = J(1, S.p, 0)
        for (i = 1; i <= cols(keep); i++) {
            if (!asarray_contains(map, keep[i]))
                _error(3498, "unmapped Stata xcheck feature: " + keep[i])
            col = asarray(map, keep[i])
            S.xkeep[col] = 1
        }
        if (sum(S.xkeep) != cols(keep))
            _error(3498, "duplicate Stata xcheck features")
    }
    kept = sum(S.xkeep)
    dropped = S.p - kept
    st_numscalar("__spmata_xcheck_kept", kept)
    st_numscalar("__spmata_xcheck_dropped", dropped)
    old = st_numscalar("__spmata_xcheck_kept_min")
    st_numscalar("__spmata_xcheck_kept_min", missing(old) ? kept : min((old, kept)))
    old = st_numscalar("__spmata_xcheck_kept_max")
    st_numscalar("__spmata_xcheck_kept_max", missing(old) ? kept : max((old, kept)))
    old = st_numscalar("__spmata_xcheck_dropped_min")
    st_numscalar("__spmata_xcheck_dropped_min", missing(old) ? dropped : min((old, dropped)))
    old = st_numscalar("__spmata_xcheck_dropped_max")
    st_numscalar("__spmata_xcheck_dropped_max", missing(old) ? dropped : max((old, dropped)))
}

real matrix _spmata_support_table(struct spdesign scalar S,
    real colvector selected)
{
    real scalar i, col, a
    real matrix out

    out = J(rows(selected), 3, .)
    for (i = 1; i <= rows(selected); i++) {
        col = selected[i]
        out[i, 1] = col
        if (col <= S.abs_p) {
            a = 1
            while (col > S.offset[a] + S.nlev[a]) a++
            out[i, 2] = S.avarpos[a]
            out[i, 3] = S.alevel[col]
        }
        else out[i, 2] = 0
    }
    return(out)
}

real colvector _spmata_weights(struct spgroup scalar cg)
{
    real scalar i
    real colvector w

    w = J(rows(cg.code), 1, .)
    for (i = 1; i <= rows(w); i++) w[i] = 1 / (cg.G * cg.count[cg.code[i]])
    return(w :/ sum(w))
}

real rowvector _spmata_mean(struct spdesign scalar S, real colvector rowsel, real colvector w)
{
    real scalar i, j, k, col, r
    real rowvector mu

    mu = J(1, S.p, 0)
    if (S.indexed) {
        for (i = 1; i <= rows(rowsel); i++) {
            r = rowsel[i]
            for (j = 1; j <= cols(S.acode); j++) {
                col = S.offset[j] + S.acode[r, j]
                mu[col] = mu[col] + w[i]
            }
        }
    }
    else {
        for (j = 1; j <= cols(S.acode); j++) {
            for (k = 1; k <= S.nlev[j]; k++) {
                col = S.offset[j] + k
                mu[col] = sum(w :* (S.acode[rowsel, j] :== k))
            }
        }
    }
    for (j = 1; j <= cols(S.C); j++) {
        mu[S.abs_p + j] = sum(w :* S.C[rowsel, j])
    }
    return(mu)
}

real rowvector _spmata_scale(struct spdesign scalar S, real colvector rowsel, real colvector w)
{
    real scalar i, j, k, c0, col, r, xmean, x2
    real rowvector scale, xsum

    scale = J(1, S.p, 1)
    if (S.indexed) {
        xsum = J(1, S.abs_p, 0)
        for (i = 1; i <= rows(rowsel); i++) {
            r = rowsel[i]
            for (j = 1; j <= cols(S.acode); j++) {
                col = S.offset[j] + S.acode[r, j]
                xsum[col] = xsum[col] + w[i]
            }
        }
        for (col = 1; col <= S.abs_p; col++) {
            xmean = xsum[col]
            scale[col] = sqrt(max((xmean - xmean^2, 1e-12)))
        }
    }
    else {
        for (j = 1; j <= cols(S.acode); j++) {
            for (k = 1; k <= S.nlev[j]; k++) {
                col = S.offset[j] + k
                xmean = sum(w :* (S.acode[rowsel, j] :== k))
                scale[col] = sqrt(max((xmean - xmean^2, 1e-12)))
            }
        }
    }
    c0 = S.abs_p
    for (j = 1; j <= cols(S.C); j++) {
        xmean = sum(w :* S.C[rowsel, j])
        x2 = sum(w :* (S.C[rowsel, j] :^ 2))
        scale[c0 + j] = sqrt(max((x2 - xmean^2, 1e-12)))
    }
    return(scale)
}

real colvector _spmata_xb(struct spdesign scalar S, real colvector rowsel,
    real colvector beta, real rowvector scale, real rowvector load)
{
    real scalar i, j, r, col
    real colvector out, ba, bc, cidx

    out = J(rows(rowsel), 1, 0)
    ba = S.abs_p ? (beta[1::S.abs_p] :/
        (scale[1::S.abs_p] :* load[1::S.abs_p])') : J(0, 1, .)
    for (j = 1; j <= cols(S.acode); j++) {
        for (i = 1; i <= rows(rowsel); i++) {
            r = rowsel[i]
            col = S.offset[j] + S.acode[r, j]
            out[i] = out[i] + ba[col]
        }
    }
    if (cols(S.C)) {
        cidx = (S.abs_p + 1)::S.p
        bc = beta[cidx] :/ scale[cidx]' :/ load[cidx]'
        out = out + S.C[rowsel, .] * bc
    }
    return(out)
}

real colvector _spmata_grad(struct spdesign scalar S, real colvector rowsel,
    real colvector v, real colvector w, real rowvector scale, real rowvector load)
{
    real scalar i, j, r, col
    real colvector g, cidx

    g = J(S.p, 1, 0)
    for (j = 1; j <= cols(S.acode); j++) {
        for (i = 1; i <= rows(rowsel); i++) {
            r = rowsel[i]
            col = S.offset[j] + S.acode[r, j]
            g[col] = g[col] + w[i] * v[i] / (scale[col] * load[col])
        }
    }
    if (cols(S.C)) {
        cidx = (S.abs_p + 1)::S.p
        g[cidx] = (S.C[rowsel, .]' * (w :* v)) :/ scale[cidx]' :/ load[cidx]'
    }
    return(g)
}

real scalar _spmata_kkt(struct spdesign scalar S, real colvector rowsel,
    real colvector y, real colvector w, real rowvector scale, real rowvector load,
    real scalar lambda, real colvector b)
{
    real scalar a
    real colvector g, pred, resid, violation

    // Match the load floor used by every lasso solver.
    load = load :+ (load :< 1e-8) :* (1e-8 :- load)
    pred = _spmata_xb(S, rowsel, b, scale, load)
    a = sum(w :* (y :- pred))
    resid = y :- a :- pred
    g = _spmata_grad(S, rowsel, resid, w, scale, load)
    violation = (b :> 0) :* abs(g :- lambda) +
        (b :< 0) :* abs(g :+ lambda) +
        (b :== 0) :* (abs(g) :- lambda) :* (abs(g) :> lambda)
    return(max(violation))
}

void _spmata_lasso_status(struct spdiag scalar diag, real scalar iterations,
    real scalar sweeps, real scalar converged)
{
    diag.lasso_calls = diag.lasso_calls + 1
    diag.lasso_iterations = diag.lasso_iterations + iterations
    diag.lasso_sweeps = diag.lasso_sweeps + sweeps
    diag.converged = diag.converged & converged
}

void _spmata_kkt_status(struct spdiag scalar diag, real scalar kkt,
    real scalar lambda)
{
    real scalar rel

    // Only the selected nuisance solution is scanned.  Scanning every BIC
    // knot or plugin reloading step is costly and not decision-relevant.
    rel = kkt / max((lambda, 1e-12))
    diag.kkt_max = max((diag.kkt_max, kkt))
    diag.kkt_rel_max = max((diag.kkt_rel_max, rel))
    diag.kkt_pass = diag.kkt_pass & (rel <= 1e-4)
}

real scalar _spmata_lip(struct spdesign scalar S, real colvector rowsel,
    real colvector w, real rowvector scale, real rowvector load)
{
    real scalar i, L
    real colvector b, xb, g

    b = J(S.p, 1, 1 / sqrt(S.p))
    for (i = 1; i <= 12; i++) {
        xb = _spmata_xb(S, rowsel, b, scale, load)
        g = _spmata_grad(S, rowsel, xb, w, scale, load)
        L = sqrt(sum(g :^ 2))
        if (L <= 1e-12) return(1)
        b = g / L
    }
    xb = _spmata_xb(S, rowsel, b, scale, load)
    L = sum(w :* (xb :^ 2))
    return(max((L, 1e-6)))
}

real colvector _spmata_soft(real colvector z, real scalar t)
{
    real colvector a
    a = abs(z) :- t
    return(((z :> 0) - (z :< 0)) :* a :* (a :> 0))
}

real colvector _spmata_keep(real colvector beta, real rowvector load,
    real scalar tol, real scalar maxk)
{
    real colvector score, ord, keep
    tol = min((tol, 1e-8))
    score = abs(beta :/ load')
    keep = selectindex(score :> tol)
    if (rows(keep) > maxk) {
        ord = order(score[keep], -1)
        keep = keep[ord[1::maxk]]
    }
    return(keep)
}

real colvector _spmata_keep_strict(real colvector beta, real rowvector load,
    real scalar tol, real scalar maxk)
{
    real colvector score, ord, keep
    tol = max((tol, 1e-8))
    score = abs(beta :/ load')
    keep = selectindex(score :> tol)
    if (rows(keep) > maxk) {
        ord = order(score[keep], -1)
        keep = keep[ord[1::maxk]]
    }
    return(keep)
}

real colvector _spmata_topcorr(struct spdesign scalar S, real colvector rowsel,
    real colvector y, real colvector w, real rowvector scale, real scalar ncols)
{
    real colvector yc, score, ord
    real scalar take

    if (S.p == 0 | ncols <= 0) return(J(0, 1, .))
    yc = y :- sum(w :* y)
    score = abs(_spmata_grad(S, rowsel, yc, w, scale, J(1, S.p, 1)))
    take = min((ncols, rows(score)))
    ord = order(score, -1)
    return(ord[1::take])
}

real colvector _spmata_col(struct spdesign scalar S, real colvector rowsel,
    real scalar col, real rowvector scale, real rowvector load)
{
    real scalar a, lev

    if (col <= S.abs_p) {
        for (a = 1; a <= cols(S.acode); a++) {
            if (col > S.offset[a] & col <= S.offset[a] + S.nlev[a]) {
                lev = col - S.offset[a]
                return((S.acode[rowsel, a] :== lev) :/ (scale[col] * load[col]))
            }
        }
    }
    return(S.C[rowsel, col - S.abs_p] :/ (scale[col] * load[col]))
}

real colvector _spmata_localpos(struct spdesign scalar S, real colvector rowsel)
{
    real colvector out

    out = J(rows(S.acode), 1, 0)
    out[rowsel] = 1::rows(rowsel)
    return(out)
}

struct spcdcache scalar _spmata_cd_cache(struct spdesign scalar S,
    real colvector rowsel, real colvector w, real rowvector scale,
    real rowvector load)
{
    struct spcdcache scalar cache
    real scalar j, c, pos, i, lo, hi, nt, xval, x2, q, xm
    real colvector localpos, xj

    localpos = _spmata_localpos(S, rowsel)
    cache.apos = localpos[S.aobs]
    cache.densepos = 1::rows(rowsel)
    cache.wx = J(rows(S.aobs), 1, 0)
    cache.xval = J(1, S.abs_p, .)
    cache.xmean = J(1, S.abs_p, 0)
    cache.hdiag = J(1, S.abs_p, 0)
    cache.densemean = J(1, cols(S.C), 0)
    cache.denseh = J(1, cols(S.C), 0)
    cache.densex = J(0, 0, .)
    if (S.abs_p == 0 | S.abs_p == S.p) cache.ordv = (1::S.p)'
    else cache.ordv = ((S.abs_p + 1)::S.p)', (1::S.abs_p)'
    for (j = 1; j <= S.abs_p; j++) {
        lo = S.astart[j]
        hi = lo + S.acount[j] - 1
        xval = 1 / (scale[j] * load[j])
        x2 = xval^2
        q = 0
        nt = 0
        for (pos = lo; pos <= hi; pos++) {
            i = cache.apos[pos]
            if (i == 0) continue
            nt++
            cache.wx[pos] = w[i] * xval
            q = q + w[i]
        }
        xm = q * xval
        cache.xval[j] = xval
        cache.xmean[j] = xm
        if (nt > 0 & nt < rows(rowsel)) {
            cache.hdiag[j] = max((q * (1 - q) * x2, 0))
        }
    }
    // Dense controls are revisited in every active sweep. Cache their
    // prescaled columns when the workspace uses at most about 64 MiB; for
    // wider designs retain the same arithmetic without the extra storage.
    if (cols(S.C) & rows(rowsel) * cols(S.C) <= 8e6) {
        cache.densex = J(rows(rowsel), cols(S.C), .)
    }
    for (c = 1; c <= cols(S.C); c++) {
        j = S.abs_p + c
        xj = S.C[rowsel, c] :/ (scale[j] * load[j])
        xm = sum(w :* xj)
        cache.densemean[c] = xm
        cache.denseh[c] = sum(w :* ((xj :- xm) :^ 2))
        if (rows(cache.densex)) cache.densex[, c] = xj
    }
    return(cache)
}

// One Gauss-Seidel sweep over the coordinates in ordv (activeonly=1
// restricts to currently nonzero coefficients). Returns the largest
// relative coefficient change; b and pred are updated in place.
real scalar _spmata_cd_sweep(struct spdesign scalar S, real colvector rowsel,
    real colvector y, real colvector w, real rowvector scale, real rowvector load,
    real scalar lambda, real colvector b, real colvector pred,
    real rowvector ordv, real scalar activeonly)
{
    real scalar jj, j, nt, a, q, xval, xm, z, h, bnew, diff, maxdiff
    real colvector xj, r

    a = sum(w :* (y :- pred))
    maxdiff = 0
    for (jj = 1; jj <= cols(ordv); jj++) {
        j = ordv[jj]
        if (activeonly & b[j] == 0) continue
        xj = _spmata_col(S, rowsel, j, scale, load)
        if (j <= S.abs_p) {
            nt = sum(xj :!= 0)
            q = sum(w :* (xj :!= 0))
            xval = 1 / (scale[j] * load[j])
            xm = q * xval
            h = (nt > 0 & nt < rows(xj)) ?
                max((q * (1 - q) * xval^2, 0)) : 0
        }
        else {
            xm = sum(w :* xj)
            h = sum(w :* ((xj :- xm) :^ 2))
        }
        r = y :- a :- pred + (xj :- xm) * b[j]
        z = sum(w :* xj :* r)
        if (h <= 1e-12) bnew = 0
        else if (z > lambda) bnew = (z - lambda) / h
        else if (z < -lambda) bnew = (z + lambda) / h
        else bnew = 0
        diff = bnew - b[j]
        if (diff != 0) {
            // [LASSO] lasso p.166 step 2b: RELATIVE coefficient change,
            // reldif convention |a-b|/(1+|b|), tolerance 1e-7.
            maxdiff = max((maxdiff, abs(diff) / (1 + abs(b[j]))))
            pred = pred + xj * diff
            a = a - xm * diff
            b[j] = bnew
        }
    }
    return(maxdiff)
}

// Support-indexed equivalent of _spmata_cd_sweep().  Absorbed indicators
// touch only the observations in their inverted-support slice.  Continuous
// controls retain the legacy dense arithmetic.  Support indices are ordered
// exactly as rowsel, so nonzero terms enter each sum in the legacy order.
real scalar _spmata_cd_sweep_indexed(struct spdesign scalar S,
    real colvector rowsel, real colvector localpos, real colvector y,
    real colvector w, real rowvector scale, real rowvector load,
    real scalar lambda, real colvector b, real colvector pred,
    real rowvector ordv, real scalar activeonly)
{
    real scalar jj, j, i, pos, nt, a, q, xm, z, h, bnew, diff, maxdiff
    real scalar lo, hi, xval, x2, rval
    real colvector xj, r, loc

    a = sum(w :* (y :- pred))
    maxdiff = 0
    for (jj = 1; jj <= cols(ordv); jj++) {
        j = ordv[jj]
        if (activeonly & b[j] == 0) continue
        if (j <= S.abs_p) {
            lo = S.astart[j]
            hi = lo + S.acount[j] - 1
            xval = 1 / (scale[j] * load[j])
            x2 = xval^2
            q = 0
            nt = 0
            for (pos = lo; pos <= hi; pos++) {
                i = localpos[S.aobs[pos]]
                if (i == 0) continue
                nt++
                q = q + w[i]
            }
            xm = q * xval
            h = (nt > 0 & nt < rows(rowsel)) ?
                max((q * (1 - q) * x2, 0)) : 0
            z = 0
            if (nt > 0 & h > 1e-12) {
                for (pos = lo; pos <= hi; pos++) {
                    i = localpos[S.aobs[pos]]
                    if (i == 0) continue
                    rval = y[i] - a - pred[i] + (xval - xm) * b[j]
                    z = z + w[i] * xval * rval
                }
            }
        }
        else {
            loc = 1::rows(rowsel)
            xj = _spmata_col(S, rowsel, j, scale, load)
            xm = sum(w :* xj)
            h = sum(w :* ((xj :- xm) :^ 2))
            r = y :- a :- pred + (xj :- xm) * b[j]
            z = sum(w :* xj :* r)
        }
        if (h <= 1e-12) bnew = 0
        else if (z > lambda) bnew = (z - lambda) / h
        else if (z < -lambda) bnew = (z + lambda) / h
        else bnew = 0
        diff = bnew - b[j]
        if (diff != 0) {
            maxdiff = max((maxdiff, abs(diff) / (1 + abs(b[j]))))
            if (j <= S.abs_p) {
                for (pos = lo; pos <= hi; pos++) {
                    i = localpos[S.aobs[pos]]
                    if (i > 0) pred[i] = pred[i] + xval * diff
                }
            }
            else pred[loc] = pred[loc] + xj * diff
            a = a - xm * diff
            b[j] = bnew
        }
    }
    return(maxdiff)
}

real scalar _spmata_cd_sweep_indexed_cached(struct spdesign scalar S,
    real colvector rowsel, real colvector y, real colvector w,
    real rowvector scale, real rowvector load, real scalar lambda,
    real colvector b, real colvector pred, real rowvector ordv,
    real scalar activeonly, struct spcdcache scalar cache)
{
    real scalar jj, j, i, pos, a, xm, z, h, bnew, diff, maxdiff
    real scalar lo, hi, xval, rval
    real colvector xj, r, loc

    a = sum(w :* (y :- pred))
    maxdiff = 0
    for (jj = 1; jj <= cols(ordv); jj++) {
        j = ordv[jj]
        if (activeonly & b[j] == 0) continue
        if (j <= S.abs_p) {
            lo = S.astart[j]
            hi = lo + S.acount[j] - 1
            xval = cache.xval[j]
            xm = cache.xmean[j]
            h = cache.hdiag[j]
            z = 0
            if (h > 1e-12) {
                for (pos = lo; pos <= hi; pos++) {
                    i = cache.apos[pos]
                    if (i == 0) continue
                    rval = y[i] - a - pred[i] + (xval - xm) * b[j]
                    z = z + cache.wx[pos] * rval
                }
            }
        }
        else {
            loc = cache.densepos
            if (rows(cache.densex)) xj = cache.densex[, j - S.abs_p]
            else xj = _spmata_col(S, rowsel, j, scale, load)
            xm = cache.densemean[j - S.abs_p]
            r = y :- a :- pred + (xj :- xm) * b[j]
            z = sum(w :* xj :* r)
            h = cache.denseh[j - S.abs_p]
        }
        if (h <= 1e-12) bnew = 0
        else if (z > lambda) bnew = (z - lambda) / h
        else if (z < -lambda) bnew = (z + lambda) / h
        else bnew = 0
        diff = bnew - b[j]
        if (diff != 0) {
            maxdiff = max((maxdiff, abs(diff) / (1 + abs(b[j]))))
            if (j <= S.abs_p) {
                for (pos = lo; pos <= hi; pos++) {
                    i = cache.apos[pos]
                    if (i > 0) pred[i] = pred[i] + xval * diff
                }
            }
            else pred[loc] = pred[loc] + xj * diff
            a = a - xm * diff
            b[j] = bnew
        }
    }
    return(maxdiff)
}

// Coordinate-equivalent sweep for pure absorbed-factor designs.  Levels of
// one factor have disjoint supports, so their residual sums are accumulated
// once while the scalar intercept and coefficients retain Gauss-Seidel order.
real scalar _spmata_cd_sweep_factor_cached(struct spdesign scalar S,
    real colvector rowsel, real colvector y, real colvector w,
    real scalar lambda, real colvector b, real colvector pred,
    real rowvector ordv, struct spcdcache scalar cache)
{
    real scalar a, f, lev, j, i, lo, hi, xm, h, z, bnew, diff, maxdiff
    real colvector allowed, code, sumr, db, xval, r

    allowed = J(S.abs_p, 1, 0)
    allowed[ordv'] = J(cols(ordv), 1, 1)
    a = sum(w :* (y :- pred))
    maxdiff = 0
    for (f = 1; f <= cols(S.nlev); f++) {
        lo = S.offset[f] + 1
        hi = lo + S.nlev[f] - 1
        code = S.acode[rowsel, f]
        sumr = J(S.nlev[f], 1, 0)
        r = y :- pred
        for (i = 1; i <= rows(code); i++)
            sumr[code[i]] = sumr[code[i]] + w[i] * r[i]

        db = J(S.nlev[f], 1, 0)
        xval = cache.xval[lo::hi]'
        for (lev = 1; lev <= S.nlev[f]; lev++) {
            j = lo + lev - 1
            if (!allowed[j]) continue
            xm = cache.xmean[j]
            h = cache.hdiag[j]
            z = xval[lev] * sumr[lev] - a * xm +
                xm * (xval[lev] - xm) * b[j]
            if (h <= 1e-12) bnew = 0
            else if (z > lambda) bnew = (z - lambda) / h
            else if (z < -lambda) bnew = (z + lambda) / h
            else bnew = 0
            diff = bnew - b[j]
            if (diff != 0) {
                maxdiff = max((maxdiff, abs(diff) / (1 + abs(b[j]))))
                db[lev] = diff
                a = a - xm * diff
                b[j] = bnew
            }
        }
        pred = pred + xval[code] :* db[code]
    }
    return(maxdiff)
}

// Jointly update every active level of each absorbed factor.  Conditional on
// the other factors, the intercept is a scalar monotone root and each level is
// a soft-threshold update.  This avoids slow one-level cycling on dense FE
// supports; ordinary and small working sets keep the established sweep.
real scalar _spmata_cd_sweep_absblock_cached(struct spdesign scalar S,
    real colvector rowsel, real colvector y, real colvector w,
    real rowvector scale, real rowvector load, real scalar lambda,
    real colvector b, real colvector pred, real rowvector activeord,
    struct spcdcache scalar cache)
{
    real scalar a, h, j, i, pos, lo, hi, first, last, xval, q, s
    real scalar left, right, mid, F, z, bnew, diff, maxdiff
    real colvector allowed, basepred, meanj, qj, newb, present

    allowed = J(S.abs_p, 1, 0)
    allowed[activeord'] = J(cols(activeord), 1, 1)
    maxdiff = 0
    for (a = 1; a <= cols(S.nlev); a++) {
        lo = S.offset[a] + 1
        hi = lo + S.nlev[a] - 1
        basepred = pred
        meanj = qj = newb = J(S.nlev[a], 1, 0)

        for (j = lo; j <= hi; j++) {
            xval = cache.xval[j]
            first = S.astart[j]
            last = first + S.acount[j] - 1
            if (b[j] != 0) {
                for (pos = first; pos <= last; pos++) {
                    i = cache.apos[pos]
                    if (i > 0) basepred[i] = basepred[i] - xval * b[j]
                }
            }
            q = xval > 0 ? cache.xmean[j] / xval : 0
            qj[j - lo + 1] = q
            if (q > 0) {
                s = 0
                for (pos = first; pos <= last; pos++) {
                    i = cache.apos[pos]
                    if (i > 0) s = s + w[i] * (y[i] - basepred[i])
                }
                meanj[j - lo + 1] = s / q
            }
        }
        present = selectindex(qj :> 0)
        if (!rows(present)) continue
        left = min(meanj[present])
        right = max(meanj[present])
        for (h = 1; h <= 50; h++) {
            mid = (left + right) / 2
            F = 0
            for (j = lo; j <= hi; j++) {
                q = qj[j - lo + 1]
                if (q <= 0) continue
                xval = cache.xval[j]
                bnew = 0
                if (allowed[j]) {
                    z = q * xval * (meanj[j - lo + 1] - mid)
                    if (z > lambda) bnew = (z - lambda) / (q * xval^2)
                    else if (z < -lambda) bnew = (z + lambda) / (q * xval^2)
                }
                F = F + q * (meanj[j - lo + 1] - mid - xval * bnew)
            }
            if (F > 0) left = mid
            else right = mid
        }

        mid = (left + right) / 2
        pred = basepred
        for (j = lo; j <= hi; j++) {
            q = qj[j - lo + 1]
            xval = cache.xval[j]
            bnew = 0
            if (q > 0 & allowed[j]) {
                z = q * xval * (meanj[j - lo + 1] - mid)
                if (z > lambda) bnew = (z - lambda) / (q * xval^2)
                else if (z < -lambda) bnew = (z + lambda) / (q * xval^2)
            }
            newb[j - lo + 1] = bnew
            diff = bnew - b[j]
            maxdiff = max((maxdiff, abs(diff) / (1 + abs(b[j]))))
            if (bnew != 0) {
                first = S.astart[j]
                last = first + S.acount[j] - 1
                for (pos = first; pos <= last; pos++) {
                    i = cache.apos[pos]
                    if (i > 0) pred[i] = pred[i] + xval * bnew
                }
            }
        }
        b[lo::hi] = newb
    }
    return(maxdiff)
}

real matrix _spmata_cd_active_block(struct spdesign scalar S,
    real colvector rowsel, real rowvector activeord, real scalar lo,
    real scalar hi, real rowvector scale, real rowvector load,
    struct spcdcache scalar cache)
{
    real scalar h, j, pos, i, c, first, last
    real matrix X
    real colvector blockrows

    X = J(hi - lo + 1, cols(activeord), 0)
    blockrows = rowsel[lo::hi]
    for (h = 1; h <= cols(activeord); h++) {
        j = activeord[h]
        if (j <= S.abs_p) {
            first = S.astart[j]
            last = first + S.acount[j] - 1
            for (pos = first; pos <= last; pos++) {
                i = cache.apos[pos]
                if (i >= lo & i <= hi) X[i - lo + 1, h] = cache.xval[j]
            }
        }
        else {
            c = j - S.abs_p
            if (rows(cache.densex)) X[, h] = cache.densex[lo::hi, c]
            else X[, h] = S.C[blockrows, c] / (scale[j] * load[j])
        }
    }
    return(X)
}

// A sign-stable working set satisfies one weighted KKT linear system.
// Accumulate it with a bounded design block; ambiguous cases fall back to CD.
real scalar _spmata_cd_polish(struct spdesign scalar S,
    real colvector rowsel, real colvector y, real colvector w,
    real rowvector scale, real rowvector load, real scalar lambda,
    real colvector b, real colvector pred, real rowvector activeord,
    struct spcdcache scalar cache)
{
    real scalar lo, hi, k, drops, wsum, ysum, bs
    real matrix X, H, G
    real colvector rhs, q, coef, newb, keep, rr, prednew
    real rowvector work, xsum

    k = cols(activeord)
    // Ordinary solver polishing is still called only for k<=64.  The larger
    // ceiling lets the explicit plugin KKT rescue polish empirical supports
    // without changing the selected path; its caller rejects support changes.
    if (k == 0 | k > 512) return(0)
    bs = max((1, floor(8e6 / k)))
    H = J(k, k, 0)
    q = J(k, 1, 0)
    xsum = J(1, k, 0)
    for (lo = 1; lo <= rows(y); lo = lo + bs) {
        hi = min((lo + bs - 1, rows(y)))
        rr = lo::hi
        X = _spmata_cd_active_block(S, rowsel, activeord, lo, hi,
            scale, load, cache)
        H = H + quadcross(X, w[rr], X)
        q = q + quadcross(X, w[rr], y[rr])
        xsum = xsum + quadcross(w[rr], X)
    }
    X = J(0, 0, .)
    wsum = sum(w)
    ysum = quadcross(w, y)
    H = H - xsum' * xsum / wsum
    q = q - xsum' * ysum / wsum

    work = (1::k)'
    drops = 0
    while (cols(work)) {
        G = H[work', work']
        if (rank(G) < cols(work) | cond(G) >= 1e12) return(0)
        rhs = q[work'] - lambda * sign(b[activeord[work]'])
        coef = qrsolve(G, rhs)
        if (!all(abs(coef) :< .)) return(0)
        newb = coef
        keep = newb :* b[activeord[work]'] :> 0
        if (all(keep)) {
            prednew = J(rows(y), 1, 0)
            bs = max((1, floor(8e6 / cols(work))))
            for (lo = 1; lo <= rows(y); lo = lo + bs) {
                hi = min((lo + bs - 1, rows(y)))
                rr = lo::hi
                X = _spmata_cd_active_block(S, rowsel,
                    activeord[work], lo, hi, scale, load, cache)
                prednew[rr] = X * newb
            }
            b[activeord'] = J(cols(activeord), 1, 0)
            b[activeord[work]'] = newb
            pred = prednew
            return(drops + 1)
        }
        drops = drops + sum(!keep)
        work = select(work', keep)'
    }
    b[activeord'] = J(cols(activeord), 1, 0)
    pred = J(rows(y), 1, 0)
    return(drops + 1)
}

real colvector _spmata_lasso_cd_cached(struct spdesign scalar S,
    real colvector rowsel, real colvector y, real colvector w,
    real rowvector scale, real rowvector load, real scalar lambda,
    real scalar maxiter, real colvector binit, real scalar polish,
    real scalar absblock,
    struct spcdcache scalar cache, struct spdiag scalar diag)
{
    real scalar iter, inner, maxdiff, iterations, sweeps, converged, hasviol
    real scalar polished, useblock, blockkkt
    real scalar a
    real colvector b, pred, resid, g
    real rowvector activeord, blockord

    if (rows(binit) == S.p) {
        b = binit
        pred = _spmata_xb(S, rowsel, b, scale, load)
    }
    else {
        b = J(S.p, 1, 0)
        pred = J(rows(y), 1, 0)
    }
    iterations = 0
    sweeps = 0
    converged = 0
    for (iter = 1; iter <= maxiter; iter++) {
        iterations++
        if (S.abs_p == S.p & st_global("SPMATA_ABSFACTOR") != "0")
            maxdiff = _spmata_cd_sweep_factor_cached(S, rowsel, y, w,
                lambda, b, pred, cache.ordv, cache)
        else maxdiff = _spmata_cd_sweep_indexed_cached(S, rowsel, y, w,
                scale, load, lambda, b, pred, cache.ordv, 0, cache)
        sweeps++
        diag.cd_full_sweeps = diag.cd_full_sweeps + 1
        if (maxdiff < 1e-7) {
            converged = 1
            break
        }
        activeord = select(cache.ordv', b[cache.ordv'] :!= 0)'
        if (absblock & S.abs_p == S.p & st_global("SPMATA_ABSBLOCK") != "0" &
                cols(activeord) > floor(.9 * S.p)) {
            blockord = select(cache.ordv',
                abs(b[cache.ordv'] :/ load[cache.ordv']') :> 1e-8)'
            if (cols(blockord) > 256 & cols(blockord) <= floor(.9 * S.p)) {
                b = b :* (abs(b :/ load') :> 1e-8)
                pred = _spmata_xb(S, rowsel, b, scale, load)
                activeord = blockord
            }
        }
        diag.cd_active_max = max((diag.cd_active_max, cols(activeord)))
        polished = 0
        if ((polish == 2 | (polish == 1 & S.p >= 50000)) &
                cols(activeord) <= 64) {
            diag.cd_polish_calls = diag.cd_polish_calls + 1
            polished = _spmata_cd_polish(S, rowsel, y, w, scale, load,
                lambda, b, pred, activeord, cache)
            diag.cd_polish_success = diag.cd_polish_success + (polished > 0)
        }
        useblock = absblock & S.abs_p == S.p & cols(activeord) > 256 &
            st_global("SPMATA_ABSBLOCK") != "0"
        maxdiff = 1
        if (useblock) {
            // Coefficients are not unique when all levels of an absorbed
            // factor are active.  Stop the joint updates on the unique KKT
            // fit, then retain the ordinary coordinate path as a fallback.
            for (inner = 1; inner <= min((maxiter, 50)); inner++) {
                maxdiff = _spmata_cd_sweep_absblock_cached(S, rowsel, y, w,
                    scale, load, lambda, b, pred, activeord, cache)
                sweeps++
                diag.cd_active_sweeps = diag.cd_active_sweeps + 1
                blockkkt = _spmata_kkt(S, rowsel, y, w, scale, load,
                    lambda, b)
                diag.cd_kkt_scans = diag.cd_kkt_scans + 1
                if (blockkkt <= max((1e-10, lambda * 1e-6))) {
                    maxdiff = 0
                    break
                }
                if (maxdiff < 1e-7) break
            }
        }
        if (!useblock | maxdiff >= 1e-7) {
            if (useblock) diag.cd_kkt_fallbacks = diag.cd_kkt_fallbacks + 1
            maxdiff = 1
            for (inner = 1; inner <= maxiter & maxdiff >= 1e-7; inner++) {
                if (S.abs_p == S.p & st_global("SPMATA_ABSFACTOR") != "0")
                    maxdiff = _spmata_cd_sweep_factor_cached(S, rowsel,
                        y, w, lambda, b, pred, activeord, cache)
                else maxdiff = _spmata_cd_sweep_indexed_cached(S, rowsel,
                        y, w, scale, load, lambda, b, pred, activeord, 1,
                        cache)
                sweeps++
                diag.cd_active_sweeps = diag.cd_active_sweeps + 1
            }
        }
        if (maxdiff >= 1e-7) continue

        // Match the active-coordinate updates of the old terminal full
        // sweep, then replace its expensive inactive-coordinate work with
        // one vectorized KKT scan. Violations return to the established
        // ordered full sweep on the next outer iteration.
        if (S.abs_p == S.p & st_global("SPMATA_ABSFACTOR") != "0")
            maxdiff = _spmata_cd_sweep_factor_cached(S, rowsel, y, w,
                lambda, b, pred, activeord, cache)
        else maxdiff = _spmata_cd_sweep_indexed_cached(S, rowsel, y, w,
                scale, load, lambda, b, pred, activeord, 1, cache)
        sweeps++
        diag.cd_active_sweeps = diag.cd_active_sweeps + 1
        if (maxdiff >= 1e-7) continue
        a = sum(w :* (y :- pred))
        resid = y :- a :- pred
        g = _spmata_grad(S, rowsel, resid, w, scale, load)
        diag.cd_kkt_scans = diag.cd_kkt_scans + 1
        hasviol = !all(abs(g) :< .) |
            any((b :== 0) :& (abs(g) :> lambda))
        if (!hasviol) {
            converged = 1
            break
        }
        diag.cd_kkt_fallbacks = diag.cd_kkt_fallbacks + 1
    }
    _spmata_lasso_status(diag, iterations, sweeps, converged)
    return(b)
}

real colvector _spmata_lasso_cd_cyclic_cached(struct spdesign scalar S,
    real colvector rowsel, real colvector y, real colvector w,
    real rowvector scale, real rowvector load, real scalar lambda,
    real scalar maxiter, real colvector binit, struct spcdcache scalar cache,
    struct spdiag scalar diag)
{
    real scalar iter, maxdiff, converged
    real colvector b, pred, b0, pred0

    if (rows(binit) == S.p) {
        b = binit
        pred = _spmata_xb(S, rowsel, b, scale, load)
    }
    else {
        b = J(S.p, 1, 0)
        pred = J(rows(y), 1, 0)
    }
    converged = 0
    for (iter = 1; iter <= maxiter; iter++) {
        maxdiff = _spmata_cd_sweep_indexed_cached(S, rowsel, y, w,
            scale, load, lambda, b, pred, cache.ordv, 0, cache)
        if (maxdiff < 1e-7) {
            converged = 1
            break
        }
    }
    _spmata_lasso_status(diag, iter, iter, converged)
    return(b)
}

real colvector _spmata_lasso_cd_screened_cached(struct spdesign scalar S,
    real colvector rowsel, real colvector y, real colvector w,
    real rowvector scale, real rowvector load, real scalar lambda,
    real scalar maxiter, real colvector binit, real rowvector workord,
    struct spcdcache scalar cache, struct spdiag scalar diag)
{
    real scalar iter, maxdiff, converged
    real colvector b, pred

    b = rows(binit) == S.p ? binit : J(S.p, 1, 0)
    pred = _spmata_xb(S, rowsel, b, scale, load)
    converged = 0
    for (iter = 1; iter <= maxiter; iter++) {
        maxdiff = _spmata_cd_sweep_indexed_cached(S, rowsel, y, w,
            scale, load, lambda, b, pred, workord, 0, cache)
        if (maxdiff < 1e-7) {
            converged = 1
            break
        }
    }
    _spmata_lasso_status(diag, iter, iter, converged)
    return(b)
}

real colvector _spmata_lasso_cd(struct spdesign scalar S, real colvector rowsel,
    real colvector y, real colvector w, real rowvector scale, real rowvector load,
    real scalar lambda, real scalar maxiter, real colvector binit,
    struct spdiag scalar diag)
{
    struct spcdcache scalar cache
    real scalar iter, inner, maxdiff, iterations, sweeps, converged
    real colvector b, pred, localpos
    real rowvector activeord, ordv

    if (S.indexed) {
        cache = _spmata_cd_cache(S, rowsel, w, scale, load)
        return(_spmata_lasso_cd_cached(S, rowsel, y, w, scale, load,
            lambda, maxiter, binit, 0, 0, cache, diag))
    }

    // [LASSO] lasso fitting p.225: warm starts — model(lambda_i) is the
    // starting point for model(lambda_i+1). binit carries the previous
    // knot's solution; empty binit = cold start (plugin call sites).
    if (rows(binit) == S.p) {
        b = binit
        pred = _spmata_xb(S, rowsel, b, scale, load)
    }
    else {
        b = J(S.p, 1, 0)
        pred = J(rows(y), 1, 0)
    }
    // Sweep CONTROLS first, then absorb-factor levels — the order of the
    // official design matrix (controls($cc) ibn.($fc)); decides which
    // collinear representative enters at near-ties.
    if (S.abs_p == 0)          ordv = (1::S.p)'
    else if (S.abs_p == S.p)   ordv = (1::S.p)'
    else ordv = ((S.abs_p + 1)::S.p)', (1::S.abs_p)'
    localpos = S.indexed ? _spmata_localpos(S, rowsel) : J(0, 1, .)
    iterations = 0
    sweeps = 0
    converged = 0
    // glmnet-style active-set strategy (Friedman-Hastie-Tibshirani 2010,
    // sec. 2.6): after each full sweep, iterate ONLY the active set to
    // convergence, then re-run a full sweep; converged when a full sweep
    // itself meets the tolerance. Resolves near-tie entry events in both
    // directions (breathe react folds, v0.8.x Finding 3 residual).
    for (iter = 1; iter <= maxiter; iter++) {
        iterations++
        if (S.indexed) {
            maxdiff = _spmata_cd_sweep_indexed(S, rowsel, localpos, y, w,
                scale, load, lambda, b, pred, ordv, 0)
        }
        else {
            maxdiff = _spmata_cd_sweep(S, rowsel, y, w, scale, load, lambda,
                b, pred, ordv, 0)
        }
        sweeps++
        if (maxdiff < 1e-7) {
            converged = 1
            break
        }
        activeord = select(ordv', b[ordv'] :!= 0)'
        for (inner = 1; inner <= maxiter; inner++) {
            if (S.indexed) {
                maxdiff = _spmata_cd_sweep_indexed(S, rowsel, localpos, y, w,
                    scale, load, lambda, b, pred, activeord, 1)
            }
            else {
                maxdiff = _spmata_cd_sweep(S, rowsel, y, w, scale, load, lambda,
                    b, pred, activeord, 1)
            }
            sweeps++
            if (maxdiff < 1e-7) break
        }
    }
    _spmata_lasso_status(diag, iterations, sweeps, converged)
    return(b)
}

real colvector _spmata_lasso_diag(struct spdesign scalar S, real colvector rowsel,
    real colvector y, real colvector w, real rowvector scale, real rowvector load,
    real scalar lambda, real scalar maxiter, real scalar high, real colvector binit,
    struct spdiag scalar diag)
{
    real scalar i, L, a, t, tnew, diff, itermax, iterations
    real colvector b, bold, theta, pred, resid, grad

    load = load :+ (load :< 1e-8) :* (1e-8 :- load)
    if (high) return(_spmata_lasso_cd(S, rowsel, y, w, scale, load,
        lambda, maxiter, binit, diag))
    L = 2 * _spmata_lip(S, rowsel, w, scale, load)
    b = (rows(binit) == S.p ? binit : J(S.p, 1, 0))
    theta = b
    t = 1
    itermax = high ? maxiter : min((maxiter, 150))
    iterations = 0
    for (i = 1; i <= itermax; i++) {
        iterations++
        pred = _spmata_xb(S, rowsel, theta, scale, load)
        a = sum(w :* (y :- pred))
        resid = y :- a :- pred
        grad = -_spmata_grad(S, rowsel, resid, w, scale, load)
        bold = b
        b = _spmata_soft(theta :- grad / L, lambda / L)
        diff = max(abs(b :- bold))
        if (diff < 1e-7) break
        tnew = (1 + sqrt(1 + 4 * t^2)) / 2
        theta = b + ((t - 1) / tnew) * (b :- bold)
        t = tnew
    }
    _spmata_lasso_status(diag, iterations, iterations, diff < 1e-7)
    return(b)
}

// Rescue a failing certificate without changing its selected support.
// Exact active-set polish is followed by the ordinary full scan.
real colvector _spmata_plugin_polish_diag(struct spdesign scalar S,
    real colvector rowsel, real colvector y, real colvector w,
    real rowvector scale, real rowvector load, real scalar lambda,
    real scalar maxiter, real colvector beta, struct spdiag scalar diag)
{
    real scalar polished
    real colvector pred
    real rowvector activeord
    struct spcdcache scalar cache

    load = load :+ (load :< 1e-8) :* (1e-8 :- load)
    cache = _spmata_cd_cache(S, rowsel, w, scale, load)
    pred = _spmata_xb(S, rowsel, beta, scale, load)
    activeord = select(cache.ordv', beta[cache.ordv'] :!= 0)'
    diag.cd_polish_calls++
    polished = _spmata_cd_polish(S, rowsel, y, w, scale, load,
        lambda, beta, pred, activeord, cache)
    diag.cd_polish_success = diag.cd_polish_success + (polished > 0)
    if (!polished) return(beta)
    return(_spmata_lasso_cd_cached(S, rowsel, y, w, scale, load,
        lambda, maxiter, beta, 0, 0, cache, diag))
}

// Stable probe API; production calls _spmata_lasso_diag() to collect status.
real colvector _spmata_lasso(struct spdesign scalar S, real colvector rowsel,
    real colvector y, real colvector w, real rowvector scale, real rowvector load,
    real scalar lambda, real scalar maxiter, real scalar high, real colvector binit)
{
    struct spdiag scalar diag
    diag = _spmata_diag_init()
    return(_spmata_lasso_diag(S, rowsel, y, w, scale, load, lambda,
        maxiter, high, binit, diag))
}

real rowvector _spmata_loadings_indexed(struct spdesign scalar S,
    real colvector rowsel, real colvector resid, struct spgroup scalar cg,
    real rowvector scale, real rowvector mu, real scalar centered)
{
    real scalar i, j, pos, col, g, lo, hi, nt, cent, base, ss, value, rbar2
    real colvector localpos, rbar, raw, stamp, touched
    real rowvector load

    localpos = _spmata_localpos(S, rowsel)
    rbar = J(cg.G, 1, 0)
    for (i = 1; i <= rows(rowsel); i++) {
        g = cg.code[i]
        rbar[g] = rbar[g] + resid[i] / cg.count[g]
    }
    rbar2 = sum(rbar :^ 2)
    load = J(1, S.p, 0)
    raw = J(cg.G, 1, 0)
    stamp = J(cg.G, 1, 0)
    touched = J(cg.G, 1, 0)

    // Absorbed indicators: start from the all-cluster centering term and
    // replace it only for clusters touched by the indicator.  No G-by-p
    // matrix is formed; total work follows nonzeros plus touched clusters.
    for (col = 1; col <= S.abs_p; col++) {
        cent = centered ? mu[col] / scale[col] : 0
        ss = centered ? cent^2 * rbar2 : 0
        lo = S.astart[col]
        hi = lo + S.acount[col] - 1
        nt = 0
        for (pos = lo; pos <= hi; pos++) {
            i = localpos[S.aobs[pos]]
            if (i == 0) continue
            g = cg.code[i]
            if (stamp[g] != col) {
                nt++
                touched[nt] = g
                stamp[g] = col
                raw[g] = 0
            }
            raw[g] = raw[g] + resid[i] / (scale[col] * cg.count[g])
        }
        for (j = 1; j <= nt; j++) {
            g = touched[j]
            base = -cent * rbar[g]
            value = raw[g] + base
            ss = ss + value^2 - base^2
        }
        load[col] = sqrt(max((ss / cg.G, 0)))
    }

    // Ordinary controls are few and dense; retain the legacy cluster-score
    // accumulation for them using one reusable G-vector.
    for (j = 1; j <= cols(S.C); j++) {
        col = S.abs_p + j
        raw = J(cg.G, 1, 0)
        for (i = 1; i <= rows(rowsel); i++) {
            g = cg.code[i]
            raw[g] = raw[g] + S.C[rowsel[i], j] * resid[i] /
                (scale[col] * cg.count[g])
        }
        if (centered) raw = raw :- rbar * (mu[col] / scale[col])
        load[col] = sqrt(sum(raw :^ 2) / cg.G)
    }
    return(load :+ (load :< 1e-10) :* (1e-10 :- load))
}

real rowvector _spmata_loadings(struct spdesign scalar S, real colvector rowsel,
    real colvector resid, struct spgroup scalar cg, real rowvector scale,
    real rowvector mu, real scalar centered)
{
    real scalar i, j, r, col
    real matrix scores
    real rowvector load
    real colvector rbar

    if (S.indexed) {
        return(_spmata_loadings_indexed(S, rowsel, resid, cg, scale, mu,
            centered))
    }

    // Equal-cluster scores of standardized x times the supplied residual.
    scores = J(cg.G, S.p, 0)
    rbar = J(cg.G, 1, 0)
    for (i = 1; i <= rows(rowsel); i++) {
        r = rowsel[i]
        rbar[cg.code[i]] = rbar[cg.code[i]] + resid[i] / cg.count[cg.code[i]]
        for (j = 1; j <= cols(S.acode); j++) {
            col = S.offset[j] + S.acode[r, j]
            scores[cg.code[i], col] = scores[cg.code[i], col] + resid[i] / (scale[col] * cg.count[cg.code[i]])
        }
        for (j = 1; j <= cols(S.C); j++) {
            col = S.abs_p + j
            scores[cg.code[i], col] = scores[cg.code[i], col] +
                S.C[r, j] * resid[i] / (scale[col] * cg.count[cg.code[i]])
        }
    }
    if (centered) scores = scores :- rbar * (mu :/ scale)
    load = sqrt(colsum(scores :^ 2) / cg.G)
    return(load :+ (load :< 1e-10) :* (1e-10 :- load))
}

real scalar _spmata_plugin_lambda(real scalar n_eff, real scalar p)
{
    real scalar gamma
    gamma = 0.1 / log(max((n_eff, p, 3)))
    return(1.1 / sqrt(n_eff) * invnormal(1 - gamma / (2 * p)))
}

real matrix _spmata_dense(struct spdesign scalar S, real colvector rowsel,
    real colvector selected, real rowvector scale)
{
    real scalar h, j, a, lev, col
    real matrix Z

    if (rows(selected) == 0) return(J(rows(rowsel), 0, .))
    Z = J(rows(rowsel), rows(selected), 0)
    for (h = 1; h <= rows(selected); h++) {
        col = selected[h]
        if (col <= S.abs_p) {
            for (a = 1; a <= cols(S.acode); a++) {
                if (col > S.offset[a] & col <= S.offset[a] + S.nlev[a]) {
                    lev = col - S.offset[a]
                    Z[, h] = (S.acode[rowsel, a] :== lev) :/ scale[col]
                    break
                }
            }
        }
        else {
            j = col - S.abs_p
            Z[, h] = S.C[rowsel, j] :/ scale[col]
        }
    }
    return(Z)
}

void _spmata_fitcheck(real matrix Z, struct spdiag scalar diag,
    real scalar final)
{
    real scalar r, c

    r = rank(Z)
    c = cond(Z)
    diag.post_fits = diag.post_fits + 1
    if (final) {
        diag.rank_valid = diag.rank_valid & (r == cols(Z)) & (c < .)
        diag.final_rank = r
        diag.final_condition = c
    }
    else {
        diag.nuisance_rank_deficient = diag.nuisance_rank_deficient + (r < cols(Z))
        if (c < .) diag.post_condition = max((diag.post_condition, c))
    }
}

real colvector _spmata_independent_columns(real matrix D,
    real scalar intercept)
{
    real scalar j, r, rnew, cnew
    real colvector kept
    real matrix Z, candidate

    Z = intercept ? J(rows(D), 1, 1) : J(rows(D), 0, .)
    kept = J(0, 1, .)
    r = rank(Z)
    for (j = 1; j <= cols(D); j++) {
        candidate = (Z, D[, j])
        rnew = rank(candidate)
        cnew = cond(candidate)
        if (rnew > r & cnew < .) {
            kept = kept \ j
            Z = candidate
            r = rnew
        }
    }
    return(kept)
}

real colvector _spmata_resid_ols(struct spdesign scalar S, real colvector rowsel,
    real colvector y, real colvector selected, real rowvector scale,
    struct spdiag scalar diag)
{
    real matrix Z
    real colvector b

    if (rows(selected) == 0) {
        _spmata_fitcheck(J(rows(y), 1, 1), diag, 0)
        return(y :- mean(y))
    }
    Z = J(rows(y), 1, 1), _spmata_dense(S, rowsel, selected, scale)
    _spmata_fitcheck(Z, diag, 0)
    b = qrsolve(Z, y)
    return(y :- Z * b)
}

real colvector _spmata_resid_wols_diag(struct spdesign scalar S, real colvector rowsel,
    real colvector y, real colvector selected, real rowvector scale,
    real colvector w, struct spdiag scalar diag, real scalar check)
{
    real matrix Z
    real colvector b, sw

    if (rows(selected) == 0) {
        if (check) _spmata_fitcheck(sqrt(w), diag, 0)
        return(y :- sum(w :* y))
    }
    Z = J(rows(y), 1, 1), _spmata_dense(S, rowsel, selected, scale)
    sw = sqrt(w)
    if (check) _spmata_fitcheck(Z :* (sw * J(1, cols(Z), 1)), diag, 0)
    b = qrsolve(Z :* (sw * J(1, cols(Z), 1)), y :* sw)
    return(y :- Z * b)
}

// Stable probe API; initialization fits are intentionally not rank-gated.
real colvector _spmata_resid_wols(struct spdesign scalar S, real colvector rowsel,
    real colvector y, real colvector selected, real rowvector scale,
    real colvector w)
{
    struct spdiag scalar diag
    diag = _spmata_diag_init()
    return(_spmata_resid_wols_diag(S, rowsel, y, selected, scale, w, diag, 0))
}

real colvector _spmata_plugin_resid_diag(struct spdesign scalar S,
    real colvector rowsel, real colvector y, real colvector selected,
    real rowvector scale, real colvector w, real scalar clustered,
    struct spdiag scalar diag, real scalar check)
{
    real colvector ow, resid

    if (!clustered) {
        return(_spmata_resid_wols_diag(S, rowsel, y, selected, scale,
            w, diag, check))
    }
    ow = J(rows(rowsel), 1, 1 / rows(rowsel))
    resid = _spmata_resid_wols_diag(S, rowsel, y, selected, scale,
        ow, diag, check)
    // Stata's unbalanced cluster plugin fits ordinary post-lasso OLS,
    // then recenters that residual for the equal-cluster objective.
    if (clustered) resid = resid :- sum(w :* resid)
    return(resid)
}

real colvector _spmata_predict_ols(struct spdesign scalar S, real colvector train,
    real colvector test, real colvector y, real colvector selected, real rowvector scale,
    struct spdiag scalar diag)
{
    real matrix Z0, Z1
    real colvector b

    if (rows(selected) == 0) {
        _spmata_fitcheck(J(rows(train), 1, 1), diag, 0)
        return(J(rows(test), 1, mean(y)))
    }
    Z0 = J(rows(train), 1, 1), _spmata_dense(S, train, selected, scale)
    Z1 = J(rows(test), 1, 1), _spmata_dense(S, test, selected, scale)
    _spmata_fitcheck(Z0, diag, 0)
    b = qrsolve(Z0, y)
    return(Z1 * b)
}

real colvector _spmata_plugin_select_core(struct spdesign scalar S, real colvector rowsel,
    real colvector y, struct spgroup scalar cg, real colvector w,
    real rowvector scale, real scalar maxiter, real scalar tol,
    real scalar high, real scalar strict_tol, struct spdiag scalar diag)
{
    real scalar step, diff, steps, lambda, clustered, used, converged, kkt, check_kkt
    real scalar support_stable, certified, path_support_same
    real rowvector load, beta_load, path_load, newload, mu
    real colvector beta, path_beta, polished_beta, selected, path_selected
    real colvector check_selected
    real colvector prev_selected, resid

    lambda = _spmata_plugin_lambda(cg.G, sum(S.xkeep))
    clustered = (cg.G < rows(rowsel))
    mu = _spmata_mean(S, rowsel, w)
    selected = _spmata_topcorr(S, rowsel, y, w, scale, S.p)
    selected = select(selected, (S.xkeep[1, selected])')
    if (rows(selected) > min((5, rows(y) - 2)))
        selected = selected[1::min((5, rows(y) - 2))]
    resid = _spmata_plugin_resid_diag(S, rowsel, y, selected, scale, w,
        clustered, diag, 0)
    load = _spmata_loadings(S, rowsel, resid, cg, scale, mu, 1) :+
        (1 :- S.xkeep) :* 1e100
    beta_load = load
    steps = high ? 15 : 3
    // [LASSO] lasso, Methods and formulas, Algorithm 2 step 4: after the
    // first iteration the penalty loadings are normalized by 1/(N - s_hat),
    // s_hat = nonzero count of the current lasso. Our loadings divide by
    // cg.G, so the robust case (cg.G == N) needs a sqrt(N/(N-s)) rescale.
    // Clustered loadings instead use ordinary post-lasso OLS residuals
    // recentered under equal-cluster weights; the helper handles that path.
    used = 0
    converged = 0
    certified = 0
    support_stable = 0
    prev_selected = J(0, 1, .)
    for (step = 1; step <= steps; step++) {
        used++
        beta_load = load
        beta = _spmata_lasso_diag(S, rowsel, y, w, scale, load, lambda,
            maxiter, high, J(0, 1, .), diag)
        if (strict_tol > 0) selected = _spmata_keep_strict(beta, load, strict_tol, rows(y) - 2)
        else selected = _spmata_keep(beta, load, tol, rows(y) - 2)
        support_stable = rows(selected) == rows(prev_selected)
        if (support_stable & rows(selected))
            support_stable = all(selected :== prev_selected)
        prev_selected = selected
        resid = _spmata_plugin_resid_diag(S, rowsel, y, selected, scale, w,
            clustered, diag, 1)
        newload = _spmata_loadings(S, rowsel, resid, cg, scale, mu, 1) :+
            (1 :- S.xkeep) :* 1e100
        if (!clustered) {
            newload = newload * sqrt(cg.G / max((cg.G - rows(selected), 1)))
        }
        diff = max(abs(newload :- load))
        load = newload
        if (diff < 1e-8) {
            converged = 1
            break
        }
    }
    // At the update cap, certify one solve on the final loadings.
    if (!converged) {
        path_selected = selected
        path_beta = beta
        path_load = beta_load
        prev_selected = path_selected
        beta_load = load
        beta = _spmata_lasso_diag(S, rowsel, y, w, scale, beta_load,
            lambda, maxiter, high, J(0, 1, .), diag)
        if (strict_tol > 0)
            check_selected = _spmata_keep_strict(beta, beta_load, strict_tol, rows(y) - 2)
        else check_selected = _spmata_keep(beta, beta_load, tol, rows(y) - 2)
        check_kkt = _spmata_kkt(S, rowsel, y, w, scale, beta_load, lambda, beta)
        if (check_kkt / max((lambda, 1e-12)) > 1e-4 & high & S.indexed) {
            beta = _spmata_plugin_polish_diag(S, rowsel, y, w, scale,
                beta_load, lambda, maxiter, beta, diag)
            if (strict_tol > 0)
                check_selected = _spmata_keep_strict(beta, beta_load,
                    strict_tol, rows(y) - 2)
            else check_selected = _spmata_keep(beta, beta_load, tol,
                rows(y) - 2)
            check_kkt = _spmata_kkt(S, rowsel, y, w, scale, beta_load,
                lambda, beta)
        }
        support_stable = rows(check_selected) == rows(prev_selected)
        if (support_stable & rows(check_selected))
            support_stable = all(check_selected :== prev_selected)
        certified = support_stable &
            (check_kkt / max((lambda, 1e-12)) <= 1e-4)
        // The cold solve is a safety certificate, not an estimator update.
        // Preserve Stata's step-15 path; production rejects if the check differs.
        selected = path_selected
        beta = path_beta
        beta_load = path_load
        kkt = _spmata_kkt(S, rowsel, y, w, scale, beta_load, lambda, beta)
    }
    else kkt = _spmata_kkt(S, rowsel, y, w, scale, beta_load, lambda, beta)
    // Post-lasso predictions use only selected support.  If the finite Stata
    // path has the right support but a loose beta, polish beta solely for its
    // final KKT certificate and refuse any support change.
    if (kkt / max((lambda, 1e-12)) > 1e-4 & high & S.indexed) {
        polished_beta = _spmata_plugin_polish_diag(S, rowsel, y, w, scale,
            beta_load, lambda, maxiter, beta, diag)
        if (strict_tol > 0)
            check_selected = _spmata_keep_strict(polished_beta, beta_load,
                strict_tol, rows(y) - 2)
        else check_selected = _spmata_keep(polished_beta, beta_load, tol,
            rows(y) - 2)
        path_support_same = rows(check_selected) == rows(selected)
        if (path_support_same & rows(selected))
            path_support_same = all(check_selected :== selected)
        if (path_support_same) {
            beta = polished_beta
            kkt = _spmata_kkt(S, rowsel, y, w, scale, beta_load, lambda, beta)
        }
    }
    diag.plugin_calls = diag.plugin_calls + 1
    diag.plugin_iterations = diag.plugin_iterations + used
    diag.plugin_converged = diag.plugin_converged & converged
    diag.plugin_support_stable = diag.plugin_support_stable & support_stable
    diag.plugin_pass = diag.plugin_pass & (converged | certified)
    diag.plugin_capped_calls = diag.plugin_capped_calls + !converged
    diag.plugin_load_diff_max = max((diag.plugin_load_diff_max, diff))
    _spmata_kkt_status(diag, kkt, lambda)
    return(selected)
}

real colvector _spmata_plugin_select(struct spdesign scalar S, real colvector rowsel,
    real colvector y, struct spgroup scalar cg, real colvector w,
    real rowvector scale, real scalar maxiter, real scalar tol,
    real scalar high, real scalar strict_tol, struct spdiag scalar diag)
{
    return(_spmata_plugin_select_core(S, rowsel, y, cg, w, scale,
        maxiter, tol, high, strict_tol, diag))
}

// Refine a selected BIC knot on the full lasso vector.  A sub-threshold dummy
// coefficient can be needed to satisfy KKT even though it is not reported in
// the post-lasso support.  The caller retains this vector only if that reported
// support is unchanged and the full KKT scan passes.
real colvector _spmata_bic_kkt_refine(struct spdesign scalar S,
    real colvector rowsel, real colvector y, real colvector w,
    real rowvector scale, real rowvector load, real scalar lambda,
    real scalar maxiter, real colvector beta, struct spdiag scalar diag)
{
    struct spcdcache scalar cache
    real scalar iter, maxdiff, kkt, cap, polished
    real colvector b, pred, b0, pred0
    real rowvector activeord

    load = load :+ (load :< 1e-8) :* (1e-8 :- load)
    cache = _spmata_cd_cache(S, rowsel, w, scale, load)
    b = beta
    pred = _spmata_xb(S, rowsel, b, scale, load)
    cap = min((maxiter, 1000))
    diag.bic_kkt_refine_calls = diag.bic_kkt_refine_calls + 1

    // Solve the chosen knot's complete nonzero working set first.  Unlike a
    // reported-support polish, this retains sub-threshold coefficients needed
    // to resolve dummy collinearity.  Ambiguous or large systems fall through
    // to the established full coordinate sweeps.
    activeord = select(cache.ordv', b[cache.ordv'] :!= 0)'
    // Dense cross-products pay only for a genuinely small BIC working set;
    // near-saturated XPO paths are faster on the indexed coordinate fallback.
    if (cols(activeord) & cols(activeord) <= 256) {
        b0 = b
        pred0 = pred
        diag.cd_polish_calls = diag.cd_polish_calls + 1
        polished = _spmata_cd_polish(S, rowsel, y, w, scale, load,
            lambda, b, pred, activeord, cache)
        diag.cd_polish_success = diag.cd_polish_success + (polished > 0)
        if (polished) {
            kkt = _spmata_kkt(S, rowsel, y, w, scale, load, lambda, b)
            if (kkt / max((lambda, 1e-12)) <= 1e-4) {
                diag.bic_kkt_refine_success =
                    diag.bic_kkt_refine_success + 1
                return(b)
            }
        }
        b = b0
        pred = pred0
    }
    for (iter = 1; iter <= cap; iter++) {
        maxdiff = _spmata_cd_sweep_indexed_cached(S, rowsel, y, w,
            scale, load, lambda, b, pred, cache.ordv, 0, cache)
        diag.lasso_iterations = diag.lasso_iterations + 1
        diag.lasso_sweeps = diag.lasso_sweeps + 1
        diag.cd_full_sweeps = diag.cd_full_sweeps + 1
        if (iter == 1 | mod(iter, 10) == 0 | maxdiff < 1e-10) {
            kkt = _spmata_kkt(S, rowsel, y, w, scale, load, lambda, b)
            if (kkt / max((lambda, 1e-12)) <= 1e-4) {
                diag.bic_kkt_refine_success = diag.bic_kkt_refine_success + 1
                return(b)
            }
        }
        if (maxdiff < 1e-13) break
    }
    return(b)
}

real colvector _spmata_bic_select(struct spdesign scalar S, real colvector rowsel,
    real colvector y, real colvector w, real rowvector scale, real scalar n_eff,
    real scalar maxiter, real scalar tol, real scalar high, real scalar is_xpo,
    real scalar production, struct spdiag scalar diag)
{
    struct spcdcache scalar cache
    real scalar i, k, bic, bestbic, consec, rel, devchange, identified, prevrss
    real scalar lambda, lambda0, rss, tss, bestkkt, bestkkt_before, candidate_kkt
    real scalar rescue, support_same
    real scalar maxk, grid, ratio, clustered, bestlambda, bestk, besti
    real scalar profile, f0, a0, s0, p0, ps0, evaluated, stop_reached
    real rowvector load
    real colvector yc, score, beta, selected, best, bestbeta, pred, resid, bprev
    real colvector polished_beta, check_selected
    string scalar source

    // v0.7.0: Stata's DOCUMENTED stopping rule, one rule for every context
    // (po/xpo, robust/cluster). Source: [LASSO] lasso, option bictolerance(#)
    // (default 1e-2), and [LASSO] lasso fitting, "Selection using the BIC
    // function":
    //   "A minimum is identified when the BIC function rises above a nominal
    //    minimum for at least two smaller lambdas with a relative difference
    //    in the BIC function greater than #."
    // The rule operates on Stata's ABSOLUTE-scale BIC
    //   BIC = -2 lnL + (k+1) ln(n_eff)
    //       = n_eff (1 + ln 2pi) + n_eff ln(RSS_w) + (k+1) ln(n_eff)
    // (Gaussian lnL; k+1 counts the intercept; n_eff = G under cluster, n
    // otherwise; RSS_w is the weighted residual sum of squares with weights
    // summing to 1).  The Stata 19 empirical tolerance microscope confirms
    // the relative rule and also shows that gridminok keeps the running BIC
    // minimum when the deviance stop is reached.  See
    // tests/probe_official_bictol_empirical.do.
    //
    // Secondary rule: stop when the Gaussian deviance-ratio change
    // |RSS_i - RSS_(i-1)| / TSS is below stop(1e-5).
    //
    // The selection rule is common to PO and XPO.  `is_xpo` only enables the
    // absorbed-block acceleration.
    load = S.xkeep :+ (1 :- S.xkeep) :* 1e100
    if (high & S.indexed) cache = _spmata_cd_cache(S, rowsel, w, scale, load)
    yc = y :- sum(w :* y)
    score = abs(_spmata_grad(S, rowsel, yc, w, scale, load))
    lambda0 = max((max(score), 1e-8))
    bestbic = .
    best = J(0, 1, .)
    bestbeta = J(S.p, 1, 0)
    bestlambda = .
    bestk = .
    besti = .
    maxk = rows(y) - 2
    tss = sum(w :* (yc :^ 2))
    clustered = (n_eff < rows(y))
    n_eff = max((n_eff, 3))
    grid = 100
    ratio = (sum(S.xkeep) < rows(y)) ? 1e-4 : 1e-2
    consec = 0
    identified = 0
    evaluated = 0
    stop_reached = 0
    prevrss = .
    bprev = J(0, 1, .)
    profile = (st_global("SPMATA_PROFILE") == "1")
    for (i = 1; i <= grid; i++) {
        lambda = exp(log(lambda0) + (i - 1) * (log(lambda0 * ratio) - log(lambda0)) / (grid - 1))
        f0 = diag.cd_full_sweeps
        a0 = diag.cd_active_sweeps
        s0 = diag.cd_kkt_scans
        p0 = diag.cd_polish_calls
        ps0 = diag.cd_polish_success
        if (high & S.indexed) {
            // Stata's finite cyclic path determines near-tie BIC supports in
            // manageable PO designs. Keep the active path for large PO/XPO.
            if (!is_xpo & S.p > 256 & S.p <= 20000)
                beta = _spmata_lasso_cd_cyclic_cached(S, rowsel, y, w,
                    scale, load, lambda, maxiter, bprev, cache, diag)
            else beta = _spmata_lasso_cd_cached(S, rowsel, y, w, scale,
                load, lambda, maxiter, bprev, clustered, is_xpo, cache, diag)
        }
        else beta = _spmata_lasso_diag(S, rowsel, y, w, scale, load,
            lambda, maxiter, high, bprev, diag)
        bprev = beta
        selected = _spmata_keep(beta, load, tol, maxk)
        pred = _spmata_xb(S, rowsel, beta, scale, load)
        resid = y :- sum(w :* (y :- pred)) :- pred
        rss = sum(w :* (resid :^ 2))
        devchange = i > 1 ? abs(rss - prevrss) /
            max((abs(tss), 1e-300)) : .
        k = rows(selected)
        bic = n_eff * (1 + log(2 * pi())) + n_eff * log(max((rss, 1e-300))) +
            (k + 1) * log(n_eff)
        evaluated = i
        if (profile) printf("  bic_knot: i=%g lambda=%9.6f k=%g full=%g active=%g kkt=%g polish=%g/%g\n",
            i, lambda, k, diag.cd_full_sweeps - f0,
            diag.cd_active_sweeps - a0, diag.cd_kkt_scans - s0,
            diag.cd_polish_success - ps0, diag.cd_polish_calls - p0)
        if (i == 1 | bic < bestbic) {
            bestbic = bic
            best = selected
            bestbeta = beta
            bestlambda = lambda
            bestk = k
            besti = i
            consec = 0
        }
        else {
            rel = (bic - bestbic) / (1 + abs(bestbic))
            consec = (rel > 0.01) ? consec + 1 : 0
            if (consec >= 2) {
                identified = 1
                break
            }
        }
        if (i > 1 & devchange < 1e-5) {
            stop_reached = 1
            break
        }
        prevrss = rss
    }
    diag.bic_calls++
    diag.bic_knots_total = diag.bic_knots_total + evaluated
    diag.bic_knots_max = max((diag.bic_knots_max, evaluated))
    if (identified) {
        diag.bic_identified++
        source = "minimum"
    }
    else {
        if (stop_reached) {
            diag.bic_stop++
            source = "running_min"
        }
        else {
            diag.bic_grid++
            source = diag.bic_gridminok ? "lambda_gmin" : "running_min"
            if (diag.bic_gridminok) {
                // A nearly saturated last knot is numerically unstable in a
                // non-unique absorbed design.  Keep the genuine BIC minimum.
                if (is_xpo & k >= floor(.95 * S.p) &
                        bestk < floor(.95 * S.p)) {
                    source = "running_min_saturated"
                    diag.bic_saturation_rescues++
                }
                else {
                    best = selected
                    bestbeta = beta
                    bestlambda = lambda
                    bestk = k
                    besti = evaluated
                    bestbic = bic
                }
            }
        }
    }
    printf("  bic_select: source=%s knots=%g lambda0=%9.6f n_eff=%g clust=%g lam*=%9.6f i*=%g k*=%g bic*=%g\n",
        source, evaluated, lambda0, n_eff, clustered, bestlambda, besti,
        bestk, bestbic)
    // The path's coefficient-difference stop can be too loose at a very small
    // selected lambda.  Repair only the chosen solution, and accept the repair
    // only when its selected support is unchanged.
    bestkkt = _spmata_kkt(S, rowsel, y, w, scale, load, bestlambda, bestbeta)
    bestkkt_before = bestkkt
    candidate_kkt = .
    rescue = 0
    support_same = 1
    if (high & S.indexed & bestkkt / max((bestlambda, 1e-12)) > 1e-4) {
        rescue = 1
        polished_beta = _spmata_bic_kkt_refine(S, rowsel, y, w, scale, load,
            bestlambda, maxiter, bestbeta, diag)
        check_selected = _spmata_keep(polished_beta, load, tol, maxk)
        candidate_kkt = _spmata_kkt(S, rowsel, y, w, scale, load,
            bestlambda, polished_beta)
        support_same = (rows(check_selected) == rows(best))
        if (support_same & rows(best)) support_same = all(check_selected :== best)
        if (support_same) {
            bestbeta = polished_beta
            bestkkt = _spmata_kkt(S, rowsel, y, w, scale, load,
                bestlambda, bestbeta)
        }
    }
    if (st_global("SPMATA_BIC_KKT_PROFILE") == "1")
        printf("  bic_selected_kkt: before=%g candidate=%g after=%g lambda=%g rescue=%g support_same=%g\n",
            bestkkt_before / max((bestlambda, 1e-12)),
            candidate_kkt / max((bestlambda, 1e-12)),
            bestkkt / max((bestlambda, 1e-12)), bestlambda, rescue,
            support_same)
    _spmata_kkt_status(diag, bestkkt, bestlambda)
    return(best)
}

// Debug: fit the full 100-knot BIC path with no stopping, dump each knot to
// an outfile. Used by _sparse_bic_probe.ado to reverse-engineer Stata's rule.
// Format: one row per knot, columns i,lambda,k,r2,bic,tss,rss.
void _spmata_bic_dump_path(struct spdesign scalar S, real colvector rowsel,
    real colvector y, real colvector w, real rowvector scale, real scalar n_eff,
    real scalar maxiter, real scalar tol, real scalar high, string scalar outpath,
    string scalar supportpath, real scalar omitcol, string scalar omitcols,
    real scalar solver,
    struct spdiag scalar diag)
{
    struct spcdcache scalar cache
    real scalar i, k, bic, bicabs, lambda, lambda0, rss, tss, r2
    real scalar maxk, grid, ratio, clustered, loss, l1, objective, kkt
    real scalar minactive, prevlambda, threshold
    real rowvector load
    real colvector yc, score, beta, selected, pred, resid, bprev, prevgrad
    real colvector omitted
    real colvector work
    real rowvector workord
    real scalar fh, sfh, h, bestbic, besti
    string rowvector names

    load = J(1, S.p, 1)
    if (omitcol >= 1 & omitcol <= S.p) load[omitcol] = 1e100
    if (strtrim(omitcols) != "") {
        omitted = strtoreal(tokens(omitcols))'
        omitted = select(omitted, omitted :>= 1 :& omitted :<= S.p)
        if (rows(omitted)) load[omitted'] = J(1, rows(omitted), 1e100)
    }
    if (high & S.indexed) cache = _spmata_cd_cache(S, rowsel, w, scale, load)
    yc = y :- sum(w :* y)
    score = abs(_spmata_grad(S, rowsel, yc, w, scale, load))
    lambda0 = max((max(score), 1e-8))
    maxk = rows(y) - 2
    tss = sum(w :* (yc :^ 2))
    clustered = (n_eff < rows(y))
    n_eff = max((n_eff, 3))
    grid = 100
    ratio = (S.p < rows(y)) ? 1e-4 : 1e-2

    fh = fopen(outpath, "w")
    fput(fh, "i,lambda,k,r2,bic,bic_abs,tss,rss,loss,l1,objective,kkt_abs,kkt_rel,min_abs_active")
    if (supportpath != "") {
        sfh = fopen(supportpath, "w")
        fput(sfh, "i,col,name,beta_standardized,coefficient")
    }
    bestbic = .
    besti = .
    bprev = J(0, 1, .)
    prevlambda = .
    prevgrad = J(0, 1, .)
    for (i = 1; i <= grid; i++) {
        lambda = exp(log(lambda0) + (i - 1) * (log(lambda0 * ratio) - log(lambda0)) / (grid - 1))
        if (high & S.indexed) {
            if (solver == 1) beta = _spmata_lasso_cd_cyclic_cached(S, rowsel,
                y, w, scale, load, lambda, maxiter, bprev, cache, diag)
            else if (solver == 2) {
                if (i == 1) workord = cache.ordv
                else {
                    threshold = 2 * lambda - prevlambda
                    work = (abs(prevgrad) :>= threshold) :|
                        (bprev :!= 0)
                    workord = select(cache.ordv', work[cache.ordv'] :== 1)'
                }
                beta = _spmata_lasso_cd_screened_cached(S, rowsel, y, w,
                    scale, load, lambda, maxiter, bprev, workord, cache, diag)
            }
            else beta = _spmata_lasso_cd_cached(S, rowsel, y, w, scale,
                load, lambda, maxiter, bprev, clustered, 1, cache, diag)
        }
        else beta = _spmata_lasso_diag(S, rowsel, y, w, scale, load,
            lambda, maxiter, high, bprev, diag)
        bprev = beta
        selected = _spmata_keep(beta, load, tol, maxk)
        pred = _spmata_xb(S, rowsel, beta, scale, load)
        resid = y :- sum(w :* (y :- pred)) :- pred
        prevgrad = _spmata_grad(S, rowsel, resid, w, scale, load)
        prevlambda = lambda
        rss = sum(w :* (resid :^ 2))
        r2 = (tss > 0 ? 1 - rss / tss : 0)
        k = rows(selected)
        loss = .5 * rss
        l1 = sum(abs(beta))
        objective = loss + lambda * l1
        kkt = _spmata_kkt(S, rowsel, y, w, scale, load, lambda, beta)
        minactive = k ? min(abs(beta[selected])) : .
        bicabs = n_eff * (1 + log(2 * pi())) +
            n_eff * log(max((rss, 1e-300))) + (k + 1) * log(n_eff)
        if (supportpath != "" & k) {
            names = _spmata_support_names(S, selected)
            for (h = 1; h <= k; h++) {
                fput(sfh, sprintf("%g,%g,%s,%21.15g,%21.15g", i,
                    selected[h], names[h], beta[selected[h]],
                    beta[selected[h]] / scale[selected[h]]))
            }
        }
        bic = n_eff * log(max((1 - r2, 1e-300))) + k * log(n_eff)
        fput(fh, sprintf("%g,%21.15g,%g,%21.15g,%21.15g,%21.15g,%21.15g,%21.15g,%21.15g,%21.15g,%21.15g,%21.15g,%21.15g,%21.15g",
            i, lambda, k, r2, bic, bicabs, tss, rss, loss, l1,
            objective, kkt, kkt / max((lambda, 1e-12)), minactive))
        if (i == 1 | bic < bestbic) {
            bestbic = bic
            besti = i
        }
        // Cost cap only: 15 knots past the running min is far beyond any
        // stopping point Stata has ever shown; deep-lambda knots are the
        // expensive ones and irrelevant to rule reverse-engineering.
        if (i - besti >= 15) break
    }
    fclose(fh)
    if (supportpath != "") fclose(sfh)
}

real matrix _spmata_prl(struct spdesign scalar S, real colvector y, real matrix D,
    struct spgroup scalar cg, real colvector w, string scalar selection,
    real scalar maxiter, real scalar tol, real scalar high, real scalar production,
    struct spdiag scalar diag, string scalar depvar, string rowvector treat,
    real scalar macrolen)
{
    real scalar j
    real colvector rowsel, selected, target, ry
    real matrix rd
    real rowvector scale

    rowsel = 1::rows(y)
    ry = J(rows(y), 1, .)
    rd = J(rows(y), cols(D), .)
    scale = _spmata_scale(S, rowsel, w)
    for (j = 0; j <= cols(D); j++) {
        target = (j == 0 ? y : D[, j])
        if (selection == "plugin") selected = _spmata_plugin_select(S, rowsel,
            target, cg, w, scale, maxiter, tol, high, 0, diag)
        else selected = _spmata_bic_select(S, rowsel, target, w, scale, cg.G,
            maxiter, tol, high, 0, production, diag)
        _spmata_support_add(diag, selected)
        _spmata_selection_record(S, selected,
            (j == 0 ? depvar : treat[j]), j + 1, ., ., j, macrolen)
        diag.selected_final = rows(selected)
        diag.selected_max = max((diag.selected_max, diag.selected_final))
        printf("%s target=%s selected=%g\n", selection, (j == 0 ? "y" : "d" + strofreal(j)), rows(selected))
        if (j == 0) ry = _spmata_resid_ols(S, rowsel, target, selected, scale, diag)
        else rd[, j] = _spmata_resid_ols(S, rowsel, target, selected, scale, diag)
    }
    return((ry, rd))
}

real matrix _spmata_xpo(struct spdesign scalar S, real colvector y, real matrix D,
    real colvector fold, struct spgroup scalar fullcg, string scalar selection,
    real scalar maxiter, real scalar tol, real scalar high, real scalar production,
    string scalar xcheckprefix, struct spdiag scalar diag,
    string scalar depvar, string rowvector treat, real scalar macrolen)
{
    real scalar f, j
    real colvector folds, train, test, target, yhat
    real matrix dhat
    real rowvector scale
    struct spgroup scalar cg
    real colvector w, selected
    string scalar xglobal, xspec

    folds = uniqrows(sort(fold, 1))
    yhat = J(rows(y), 1, 0)
    dhat = J(rows(y), cols(D), 0)
    for (f = 1; f <= rows(folds); f++) {
        if (xcheckprefix != "") {
            xglobal = xcheckprefix + strtrim(strofreal(f))
            xspec = st_global(xglobal)
            if (strtrim(xspec) == "")
                _error(3498, "missing or empty XPO xcheck global: " + xglobal)
            _spmata_apply_xcheck(S, xspec)
        }
        train = selectindex(fold :!= folds[f])
        test = selectindex(fold :== folds[f])
        cg = _spmata_group(fullcg.code[train])
        w = (fullcg.G > 1 ? _spmata_weights(cg) : J(rows(train), 1, 1 / rows(train)))
        scale = _spmata_scale(S, train, w)
        for (j = 0; j <= cols(D); j++) {
            target = (j == 0 ? y : D[, j])
            if (selection == "plugin") selected = _spmata_plugin_select(S,
                train, target[train], cg, w, scale, maxiter, tol, high, 0, diag)
            else selected = _spmata_bic_select(S, train, target[train], w,
                scale, cg.G, maxiter, tol, high, 1, production, diag)
            _spmata_support_add(diag, selected)
            _spmata_selection_record(S, selected,
                (j == 0 ? depvar : treat[j]),
                (f - 1) * (cols(D) + 1) + j + 1, folds[f], ., j,
                macrolen)
            diag.selected_final = rows(selected)
            diag.selected_max = max((diag.selected_max, diag.selected_final))
            printf("fold=%g %s target=%s selected=%g\n", f, selection, (j == 0 ? "y" : "d" + strofreal(j)), rows(selected))
            if (j == 0) yhat[test] = _spmata_predict_ols(S, train, test,
                target[train], selected, scale, diag)
            else dhat[test, j] = _spmata_predict_ols(S, train, test,
                target[train], selected, scale, diag)
        }
    }
    return((y :- yhat, D :- dhat))
}

real matrix _spmata_xpo2(struct spdesign scalar S, real colvector y, real matrix D,
    real colvector fold1, real colvector fold2, struct spgroup scalar fullcg,
    string scalar selection, real scalar maxiter, real scalar tol,
    real scalar high, real scalar production, struct spdiag scalar diag,
    string scalar depvar, string rowvector treat, real scalar macrolen)
{
    real scalar f1, f2, j
    real colvector folds1, folds2, train, test, target, yhat
    real matrix dhat
    real rowvector scale
    struct spgroup scalar cg
    real colvector w, selected

    folds1 = uniqrows(sort(fold1, 1))
    folds2 = uniqrows(sort(fold2, 1))
    yhat = J(rows(y), 1, 0)
    dhat = J(rows(y), cols(D), 0)
    for (f1 = 1; f1 <= rows(folds1); f1++) {
        for (f2 = 1; f2 <= rows(folds2); f2++) {
            train = selectindex((fold1 :!= folds1[f1]) :&
                (fold2 :!= folds2[f2]))
            test = selectindex((fold1 :== folds1[f1]) :&
                (fold2 :== folds2[f2]))
            if (!rows(train) | !rows(test))
                _error(3498, "empty two-way cluster-DML train/test cell")
            cg = _spmata_group(fullcg.code[train])
            w = (fullcg.G > 1 ? _spmata_weights(cg) :
                J(rows(train), 1, 1 / rows(train)))
            scale = _spmata_scale(S, train, w)
            for (j = 0; j <= cols(D); j++) {
                target = (j == 0 ? y : D[, j])
                if (selection == "plugin") selected = _spmata_plugin_select(S,
                    train, target[train], cg, w, scale, maxiter, tol, high,
                    0, diag)
                else selected = _spmata_bic_select(S, train, target[train], w,
                    scale, cg.G, maxiter, tol, high, 1, production, diag)
                _spmata_support_add(diag, selected)
                _spmata_selection_record(S, selected,
                    (j == 0 ? depvar : treat[j]),
                    ((f1 - 1) * rows(folds2) + f2 - 1) *
                        (cols(D) + 1) + j + 1,
                    folds1[f1], folds2[f2], j, macrolen)
                diag.selected_final = rows(selected)
                diag.selected_max = max((diag.selected_max,
                    diag.selected_final))
                printf("fold=(%g,%g) %s target=%s selected=%g\n", f1, f2,
                    selection, (j == 0 ? "y" : "d" + strofreal(j)),
                    rows(selected))
                if (j == 0) yhat[test] = _spmata_predict_ols(S, train, test,
                    target[train], selected, scale, diag)
                else dhat[test, j] = _spmata_predict_ols(S, train, test,
                    target[train], selected, scale, diag)
            }
        }
    }
    return((y :- yhat, D :- dhat))
}

real matrix _spmata_sandwich(real matrix D, real colvector y, string scalar vce,
    real colvector clcode, real scalar G, real scalar intercept,
    struct spdiag scalar diag)
{
    real scalar i, n, q, qfit
    real matrix Z, Dr, bread, meat, scores, V, Vfull
    real rowvector beta, bfull
    real colvector resid, kept, omitted, coefidx, mask

    n = rows(y)
    q = cols(D)
    kept = _spmata_independent_columns(D, intercept)
    mask = J(q, 1, 1)
    if (rows(kept)) mask[kept] = J(rows(kept), 1, 0)
    omitted = selectindex(mask)
    diag.omitted_treat = omitted
    qfit = rows(kept)
    Dr = qfit ? D[, kept] : J(n, 0, .)
    Z = intercept ? (J(n, 1, 1), Dr) : Dr
    _spmata_fitcheck(Z, diag, 1)
    beta = qrsolve(Z, y)'
    resid = y - Z * beta'
    bread = invsym(Z' * Z)
    if (vce == "cluster") {
        scores = J(G, cols(Z), 0)
        for (i = 1; i <= n; i++) scores[clcode[i], .] = scores[clcode[i], .] + resid[i] * Z[i, .]
        meat = scores' * scores
        V = bread * meat * bread
        if (G > 1 & n > qfit) V = V * (G / (G - 1)) * ((n - 1) / (n - qfit))
    }
    else {
        meat = J(cols(Z), cols(Z), 0)
        for (i = 1; i <= n; i++) meat = meat + (resid[i]^2) * (Z[i, .]' * Z[i, .])
        V = bread * meat * bread
        if (n > qfit) V = V * n / (n - qfit)
    }
    if (intercept) {
        coefidx = qfit ? 2::(qfit + 1) : J(0, 1, .)
        beta = qfit ? beta[1, coefidx'] : J(1, 0, .)
        V = qfit ? V[coefidx, coefidx] : J(0, 0, .)
    }
    bfull = J(1, q, 0)
    Vfull = J(q, q, 0)
    if (qfit) {
        bfull[1, kept'] = beta
        Vfull[kept, kept] = V
    }
    return((bfull \ Vfull))
}

real matrix _spmata_cluster_vce(real matrix Z, real colvector resid,
    real colvector clcode, real scalar G, real scalar n, real scalar q)
{
    real scalar i
    real matrix scores, bread, V

    if (G <= 1) _error(498, "clustered VCE requires at least two clusters")
    scores = J(G, cols(Z), 0)
    for (i = 1; i <= n; i++) {
        scores[clcode[i], .] = scores[clcode[i], .] + resid[i] * Z[i, .]
    }
    bread = invsym(Z' * Z)
    V = bread * (scores' * scores) * bread
    if (n > q) V = V * (G / (G - 1)) * ((n - 1) / (n - q))
    return(V)
}

real matrix _spmata_sandwich2(real matrix D, real colvector y,
    real colvector cluster1, real colvector cluster2, real scalar intercept,
    struct spdiag scalar diag)
{
    real scalar n, q, qfit
    real rowvector beta, bfull, L
    real colvector resid, kept, omitted, coefidx, mask
    real matrix Z, Dr, V, Vfull, Vwork, Q
    struct spgroup scalar cg1, cg2, cg12

    n = rows(y)
    q = cols(D)
    kept = _spmata_independent_columns(D, intercept)
    mask = J(q, 1, 1)
    if (rows(kept)) mask[kept] = J(rows(kept), 1, 0)
    omitted = selectindex(mask)
    diag.omitted_treat = omitted
    qfit = rows(kept)
    Dr = qfit ? D[, kept] : J(n, 0, .)
    Z = intercept ? (J(n, 1, 1), Dr) : Dr
    _spmata_fitcheck(Z, diag, 1)
    beta = qrsolve(Z, y)'
    resid = y - Z * beta'
    cg1 = _spmata_group(cluster1)
    cg2 = _spmata_group(cluster2)
    cg12 = _spmata_group((cg1.code :- 1) :* cg2.G :+ cg2.code)

    V = _spmata_cluster_vce(Z, resid, cg1.code, cg1.G, n, qfit) +
        _spmata_cluster_vce(Z, resid, cg2.code, cg2.G, n, qfit) -
        _spmata_cluster_vce(Z, resid, cg12.code, cg12.G, n, qfit)
    V = (V + V') / 2
    L = symeigenvalues(V)
    diag.vce_raw_min_eigen = min(L)
    if (diag.vce_raw_min_eigen < 0) {
        Vwork = V
        symeigensystem(Vwork, Q, L)
        L = L :* (L :> 0)
        V = Q * diag(L) * Q'
        V = (V + V') / 2
        diag.vce_psd_adjusted = 1
    }
    if (intercept) {
        coefidx = qfit ? 2::(qfit + 1) : J(0, 1, .)
        beta = qfit ? beta[1, coefidx'] : J(1, 0, .)
        V = qfit ? V[coefidx, coefidx] : J(0, 0, .)
    }
    bfull = J(1, q, 0)
    Vfull = J(q, q, 0)
    if (qfit) {
        bfull[1, kept'] = beta
        Vfull[kept, kept] = V
    }
    return((bfull \ Vfull))
}

real matrix _spmata_dml2(real matrix D, real colvector y, string scalar vce,
    real colvector clcode, real scalar G, real colvector fold,
    struct spdiag scalar diag)
{
    real scalar i, g, n, q, balanced
    real rowvector beta
    real colvector resid, foldweight, pairweight, sqrtweight
    real matrix bread, meat, scores, V, Dw
    struct spgroup scalar fg

    n = rows(y)
    q = cols(D)
    _spmata_fitcheck(D, diag, 1)
    beta = qrsolve(D, y)'
    resid = y - D * beta'

    // xporegress averages the Jacobian and score second moments over
    // folds, giving each fold weight 1/K even when n_k differs.  The old
    // pooled sandwich is identical only when every fold has the same size.
    // Keep that old path for balanced folds to preserve the certified
    // equal-fold baseline exactly, including floating-point accumulation.
    fg = _spmata_group(fold)
    balanced = (min(fg.count) == max(fg.count))
    if (balanced) {
        bread = invsym(D' * D)
        if (vce == "cluster") {
            scores = J(G, q, 0)
            for (i = 1; i <= n; i++) scores[clcode[i], .] = scores[clcode[i], .] + resid[i] * D[i, .]
            meat = scores' * scores
        }
        else {
            meat = J(q, q, 0)
            for (i = 1; i <= n; i++) meat = meat + (resid[i]^2) * (D[i, .]' * D[i, .])
        }
    }
    else {
        // n/(K*n_k) converts the fold-average moments into unnormalised
        // bread/meat matrices whose sandwich has the usual O(1/n) scale.
        foldweight = J(n, 1, .)
        for (i = 1; i <= n; i++) {
            foldweight[i] = n / (fg.G * fg.count[fg.code[i]])
        }
        sqrtweight = sqrt(foldweight)
        Dw = D :* (sqrtweight * J(1, q, 1))
        bread = invsym(Dw' * Dw)
        if (vce == "cluster") {
            // clcode indexes cluster-by-fold pairs.  Scores are first
            // summed within each pair, then weighted by that fold's n_k.
            scores = J(G, q, 0)
            pairweight = J(G, 1, .)
            for (i = 1; i <= n; i++) {
                scores[clcode[i], .] = scores[clcode[i], .] + resid[i] * D[i, .]
                pairweight[clcode[i]] = foldweight[i]
            }
            meat = J(q, q, 0)
            for (g = 1; g <= G; g++) {
                meat = meat + pairweight[g] * (scores[g, .]' * scores[g, .])
            }
        }
        else {
            meat = J(q, q, 0)
            for (i = 1; i <= n; i++) {
                meat = meat + foldweight[i] * (resid[i]^2) * (D[i, .]' * D[i, .])
            }
        }
    }
    V = bread * meat * bread
    return((beta \ V))
}

real matrix _spmata_dml_cluster1(real matrix D, real colvector y,
    real colvector cluster, real colvector fold, struct spdiag scalar diag)
{
    real scalar i, g, n, q, K
    real rowvector beta
    real colvector cfold, Gk, weight, sqrtweight, resid
    real matrix Dw, bread, meat, scores, V
    struct spgroup scalar cg, fg, cfg

    n = rows(y)
    q = cols(D)
    cg = _spmata_group(cluster)
    fg = _spmata_group(fold)
    cfold = J(cg.G, 1, .)
    for (i = 1; i <= n; i++) {
        g = cg.code[i]
        if (missing(cfold[g])) cfold[g] = fold[i]
        else if (cfold[g] != fold[i])
            _error(3498, "cluster leaks across outer folds")
    }
    cfg = _spmata_group(cfold)
    K = cfg.G
    if (K != fg.G) _error(3498, "cluster-fold partition is incomplete")
    Gk = J(K, 1, 0)
    for (g = 1; g <= cg.G; g++)
        Gk[cfg.code[g]] = Gk[cfg.code[g]] + 1

    weight = J(n, 1, .)
    for (i = 1; i <= n; i++)
        weight[i] = 1 / (K * Gk[cfg.code[cg.code[i]]])
    sqrtweight = sqrt(weight)
    Dw = D :* (sqrtweight * J(1, q, 1))
    _spmata_fitcheck(Dw, diag, 1)
    beta = qrsolve(Dw, sqrtweight :* y)'
    resid = y - D * beta'
    bread = invsym(Dw' * Dw)

    scores = J(cg.G, q, 0)
    for (i = 1; i <= n; i++)
        scores[cg.code[i], .] = scores[cg.code[i], .] +
            resid[i] * D[i, .]
    meat = J(q, q, 0)
    for (g = 1; g <= cg.G; g++)
        meat = meat + (scores[g, .]' * scores[g, .]) /
            (K * Gk[cfg.code[g]])
    V = bread * meat * bread / cg.G
    return((beta \ V))
}

real matrix _spmata_dml_cluster2(real matrix D, real colvector y,
    real colvector cluster1, real colvector cluster2,
    real colvector fold1, real colvector fold2, struct spdiag scalar diag)
{
    real scalar i, g, h, k, l, n, q, K1, K2, pos, c
    real rowvector beta
    real colvector cfold1, cfold2, G1k, G2l, weight, sqrtweight, resid
    real matrix Dw, bread, meat, scores1, scores2, V
    struct spgroup scalar cg1, cg2, fg1, fg2, cfg1, cfg2

    n = rows(y)
    q = cols(D)
    cg1 = _spmata_group(cluster1)
    cg2 = _spmata_group(cluster2)
    fg1 = _spmata_group(fold1)
    fg2 = _spmata_group(fold2)
    cfold1 = J(cg1.G, 1, .)
    cfold2 = J(cg2.G, 1, .)
    for (i = 1; i <= n; i++) {
        g = cg1.code[i]
        h = cg2.code[i]
        if (missing(cfold1[g])) cfold1[g] = fold1[i]
        else if (cfold1[g] != fold1[i])
            _error(3498, "first cluster dimension leaks across folds")
        if (missing(cfold2[h])) cfold2[h] = fold2[i]
        else if (cfold2[h] != fold2[i])
            _error(3498, "second cluster dimension leaks across folds")
    }
    cfg1 = _spmata_group(cfold1)
    cfg2 = _spmata_group(cfold2)
    K1 = cfg1.G
    K2 = cfg2.G
    if (K1 != fg1.G | K2 != fg2.G)
        _error(3498, "two-way cluster-fold partition is incomplete")
    G1k = J(K1, 1, 0)
    G2l = J(K2, 1, 0)
    for (g = 1; g <= cg1.G; g++)
        G1k[cfg1.code[g]] = G1k[cfg1.code[g]] + 1
    for (h = 1; h <= cg2.G; h++)
        G2l[cfg2.code[h]] = G2l[cfg2.code[h]] + 1

    weight = J(n, 1, .)
    for (i = 1; i <= n; i++)
        weight[i] = 1 / (K1 * K2 * G1k[fg1.code[i]] * G2l[fg2.code[i]])
    sqrtweight = sqrt(weight)
    Dw = D :* (sqrtweight * J(1, q, 1))
    _spmata_fitcheck(Dw, diag, 1)
    beta = qrsolve(Dw, sqrtweight :* y)'
    resid = y - D * beta'
    bread = invsym(Dw' * Dw)

    scores1 = J(cg1.G * K2, q, 0)
    scores2 = J(cg2.G * K1, q, 0)
    for (i = 1; i <= n; i++) {
        pos = (cg1.code[i] - 1) * K2 + fg2.code[i]
        scores1[pos, .] = scores1[pos, .] + resid[i] * D[i, .]
        pos = (cg2.code[i] - 1) * K1 + fg1.code[i]
        scores2[pos, .] = scores2[pos, .] + resid[i] * D[i, .]
    }
    meat = J(q, q, 0)
    for (g = 1; g <= cg1.G; g++) {
        k = cfg1.code[g]
        for (l = 1; l <= K2; l++) {
            c = min((G1k[k], G2l[l])) /
                ((G1k[k] * G2l[l])^2 * K1 * K2)
            pos = (g - 1) * K2 + l
            meat = meat + c * (scores1[pos, .]' * scores1[pos, .])
        }
    }
    for (h = 1; h <= cg2.G; h++) {
        l = cfg2.code[h]
        for (k = 1; k <= K1; k++) {
            c = min((G1k[k], G2l[l])) /
                ((G1k[k] * G2l[l])^2 * K1 * K2)
            pos = (h - 1) * K1 + k
            meat = meat + c * (scores2[pos, .]' * scores2[pos, .])
        }
    }
    V = bread * meat * bread / min((cg1.G, cg2.G))
    return((beta \ V))
}

// Used by _sparse_bic_probe.ado. Builds sparse design, runs full-path dump.
void _spmata_probe_run(string scalar target, string scalar absorb,
    string scalar controls, string scalar cluster, string scalar touse,
    real scalar maxiter, real scalar tol, real scalar high,
    string scalar outpath, string scalar kernel, string scalar supportpath,
    string scalar weighting, real scalar omitcol, string scalar omitcols,
    real scalar solver)
{
    real colvector idx, y, clraw, w
    real scalar n_eff
    real rowvector scale
    struct spdesign scalar S
    struct spdiag scalar diag
    struct spgroup scalar cg

    idx = selectindex(st_data(., touse) :== 1)
    y = st_data(idx, target)
    S = _spmata_build_kernel(idx, absorb, controls, kernel)
    if (cluster == "") {
        clraw = 1::rows(y)
        cg = _spmata_group(clraw)
        n_eff = rows(y)
        w = J(rows(y), 1, 1 / rows(y))
    }
    else {
        clraw = st_data(idx, cluster)
        cg = _spmata_group(clraw)
        n_eff = cg.G
        if (weighting == "observation") w = J(rows(y), 1, 1 / rows(y))
        else w = _spmata_weights(cg)
    }
    scale = _spmata_scale(S, 1::rows(y), w)
    diag = _spmata_diag_init()
    _spmata_bic_dump_path(S, 1::rows(y), y, w, scale, n_eff, maxiter,
        tol, high, outpath, supportpath, omitcol, omitcols, solver, diag)

    st_numscalar("__spprobe_n", rows(y))
    st_numscalar("__spprobe_p", S.p)
    st_numscalar("__spprobe_g", cg.G)
}

void _spmata_run(string scalar cmd, string scalar depvar, string scalar treats,
    string scalar absorb, string scalar controls, string scalar touse,
    string scalar selection, string scalar vce, string scalar clusters,
    string scalar foldvar, string scalar foldvar2, real scalar xfolds,
    real scalar rseed,
    real scalar maxiter, real scalar tol, string scalar accuracy,
    string scalar strict, string scalar kernel, real scalar clusterfolds,
    string scalar bicmode, real scalar xcheck, string scalar xcheckprefix,
    real scalar matsize,
    real scalar macrolen)
{
    real colvector idx, y, clraw, fold, fold2, w, ry, support
    real matrix C, D, rd, V, R, out
    real rowvector b
    real scalar high, production, nclust, G1, G2, G12
    string rowvector tvars, cvars, support_names, omitted_names
    string scalar support_text
    struct spdesign scalar S
    struct spdiag scalar diag
    struct spgroup scalar cg, cg1, cg2, cg12, vcg

    // strict = "stata" (default) -> production=0, preserves v0.3.0/v0.4.0 behavior
    //   under `strict(stata)`; matches baseline_20260704.csv to the digit.
    // strict = "production"      -> production=1, uses documented BIC path rule
    //   in _spmata_bic_select (see docs/design_notes.md).
    production = (strict == "production")
    high = (accuracy == "high")
    diag = _spmata_diag_init()
    idx = selectindex(st_data(., touse) :== 1)
    tvars = tokens(treats)
    y = st_data(idx, depvar)
    D = st_data(idx, tvars)
    st_numscalar("__spmata_xcheck_kept_min", .)
    st_numscalar("__spmata_xcheck_kept_max", .)
    st_numscalar("__spmata_xcheck_dropped_min", .)
    st_numscalar("__spmata_xcheck_dropped_max", .)
    S = _spmata_build_kernel_opts(idx, absorb, controls, kernel,
        xcheck | xcheckprefix != "")
    if (xcheckprefix == "")
        _spmata_apply_xcheck(S, xcheck ? st_global("e(allvars)") : "")
    diag.bic_gridminok = (bicmode == "gridminok")
    diag.selected_mask = J(S.p, 1, 0)
    st_numscalar("__spmata_lasso_fits", 0)
    cvars = tokens(clusters)
    nclust = cols(cvars)
    G1 = G2 = G12 = .
    if (!nclust) {
        clraw = 1::rows(y)
        cg = _spmata_group(clraw)
    }
    else if (nclust == 1) {
        clraw = st_data(idx, cvars[1])
        cg = _spmata_group(clraw)
        cg1 = cg
        G1 = cg.G
    }
    else {
        C = st_data(idx, cvars)
        cg1 = _spmata_group(C[, 1])
        cg2 = _spmata_group(C[, 2])
        cg12 = _spmata_group((cg1.code :- 1) :* cg2.G :+ cg2.code)
        G1 = cg1.G
        G2 = cg2.G
        G12 = cg12.G
        // DoubleML learners are observation-weighted; clustering enters the
        // disjoint split and orthogonal-score covariance below.
        clraw = 1::rows(y)
        cg = _spmata_group(clraw)
    }
    w = (vce == "cluster" & nclust == 1 ? _spmata_weights(cg) :
        J(rows(y), 1, 1 / rows(y)))
    if (cmd == "sparse_poreg") {
        R = _spmata_prl(S, y, D, cg, w, selection, maxiter, tol, high,
            production, diag, depvar, tvars, macrolen)
        ry = R[, 1]
        rd = R[, (2::(cols(D) + 1))]
        if (vce == "cluster" & nclust == 2)
            out = _spmata_sandwich2(rd, ry, C[, 1], C[, 2], 1, diag)
        else out = _spmata_sandwich(rd, ry, vce, cg.code, cg.G, 1, diag)
    }
    else {
        fold = st_data(idx, foldvar)
        if (clusterfolds & nclust == 2) {
            fold2 = st_data(idx, foldvar2)
            R = _spmata_xpo2(S, y, D, fold, fold2, cg, selection, maxiter,
                tol, high, production, diag, depvar, tvars, macrolen)
        }
        else R = _spmata_xpo(S, y, D, fold, cg, selection, maxiter, tol,
            high, production, xcheckprefix, diag, depvar, tvars, macrolen)
        ry = R[, 1]
        rd = R[, (2::(cols(D) + 1))]
        if (clusterfolds & nclust == 2)
            out = _spmata_dml_cluster2(rd, ry, C[, 1], C[, 2], fold,
                fold2, diag)
        else if (clusterfolds & nclust == 1)
            out = _spmata_dml_cluster1(rd, ry, cg.code, fold, diag)
        else if (vce == "cluster") {
            vcg = _spmata_group(cg.code :* (max(fold) + 1) :+ fold)
            out = _spmata_dml2(rd, ry, vce, vcg.code, vcg.G, fold, diag)
        }
        else out = _spmata_dml2(rd, ry, vce, cg.code, cg.G, fold, diag)
    }
    b = out[1, .]
    V = out[(2::(cols(D) + 1)), .]
    diag.finite = all(abs(b) :< .) & all(abs(V) :< .) &
        all(abs(ry) :< .) & all(abs(rd) :< .)
    support = selectindex(diag.selected_mask)
    support_names = _spmata_support_names(S, support)
    omitted_names = rows(diag.omitted_treat) ?
        tvars[diag.omitted_treat'] : J(1, 0, "")
    support_text = invtokens(support_names)
    if (xcheckprefix != "") {
        st_numscalar("__spmata_xcheck_kept",
            st_numscalar("__spmata_xcheck_kept_min"))
        st_numscalar("__spmata_xcheck_dropped",
            st_numscalar("__spmata_xcheck_dropped_max"))
    }

    st_matrix("__spmata_b", b)
    st_matrix("__spmata_V", V)
    if (rows(support) & rows(support) <= matsize) st_matrix("__spmata_support",
        _spmata_support_table(S, support))
    st_local("__spmata_controls_sel",
        strlen(support_text) <= macrolen ? support_text : "")
    st_local("__spmata_treat_omitted", invtokens(omitted_names))
    st_numscalar("__spmata_N", rows(y))
    st_numscalar("__spmata_p", S.p)
    st_numscalar("__spmata_raw_p", S.raw_p)
    st_numscalar("__spmata_indexed", S.indexed)
    st_numscalar("__spmata_cluster_ways", nclust)
    st_numscalar("__spmata_N_clust1", G1)
    st_numscalar("__spmata_N_clust2", G2)
    st_numscalar("__spmata_N_clust12", G12)
    st_numscalar("__spmata_N_clust", nclust == 2 ? min((G1, G2)) : G1)
    st_numscalar("__spmata_vce_psd_adjusted", diag.vce_psd_adjusted)
    st_numscalar("__spmata_vce_raw_min_eigen", diag.vce_raw_min_eigen)
    st_numscalar("__spmata_vce_rank", rank(V))
    st_numscalar("__spmata_wald", rank(V) ?
        max((0, b * invsym(V) * b')) : .)
    st_numscalar("__spmata_score_folds", cmd == "sparse_xporeg" ?
        (clusterfolds & nclust == 2 ? xfolds^2 : xfolds) : .)
    st_numscalar("__spmata_converged", diag.converged)
    st_numscalar("__spmata_lasso_calls", diag.lasso_calls)
    st_numscalar("__spmata_lasso_iterations", diag.lasso_iterations)
    st_numscalar("__spmata_lasso_sweeps", diag.lasso_sweeps)
    st_numscalar("__spmata_cd_full_sweeps", diag.cd_full_sweeps)
    st_numscalar("__spmata_cd_active_sweeps", diag.cd_active_sweeps)
    st_numscalar("__spmata_cd_kkt_scans", diag.cd_kkt_scans)
    st_numscalar("__spmata_cd_kkt_fallbacks", diag.cd_kkt_fallbacks)
    st_numscalar("__spmata_cd_active_max", diag.cd_active_max)
    st_numscalar("__spmata_cd_polish_calls", diag.cd_polish_calls)
    st_numscalar("__spmata_cd_polish_success", diag.cd_polish_success)
    st_numscalar("__spmata_bic_kkt_refine_calls", diag.bic_kkt_refine_calls)
    st_numscalar("__spmata_bic_kkt_refine_success", diag.bic_kkt_refine_success)
    st_numscalar("__spmata_kkt_max", diag.kkt_max)
    st_numscalar("__spmata_kkt_rel_max", diag.kkt_rel_max)
    st_numscalar("__spmata_kkt_pass", diag.kkt_pass)
    st_numscalar("__spmata_plugin_converged",
        diag.plugin_calls ? diag.plugin_converged : .)
    st_numscalar("__spmata_plugin_pass",
        diag.plugin_calls ? diag.plugin_pass : .)
    st_numscalar("__spmata_plugin_support_stable",
        diag.plugin_calls ? diag.plugin_support_stable : .)
    st_numscalar("__spmata_plugin_capped_calls",
        diag.plugin_calls ? diag.plugin_capped_calls : .)
    st_numscalar("__spmata_plugin_load_diff_max",
        diag.plugin_calls ? diag.plugin_load_diff_max : .)
    st_numscalar("__spmata_plugin_iterations", diag.plugin_iterations)
    st_numscalar("__spmata_selected_max", diag.selected_max)
    st_numscalar("__spmata_selected_final", diag.selected_final)
    st_numscalar("__spmata_k_controls_sel", rows(support))
    st_numscalar("__spmata_k_treat_omitted", rows(diag.omitted_treat))
    st_numscalar("__spmata_controls_sel_complete",
        strlen(support_text) <= macrolen)
    st_numscalar("__spmata_support_matrix",
        rows(support) > 0 & rows(support) <= matsize)
    st_numscalar("__spmata_bic_calls", diag.bic_calls)
    st_numscalar("__spmata_bic_identified", diag.bic_identified)
    st_numscalar("__spmata_bic_stop", diag.bic_stop)
    st_numscalar("__spmata_bic_grid", diag.bic_grid)
    st_numscalar("__spmata_bic_saturation_rescues",
        diag.bic_saturation_rescues)
    st_numscalar("__spmata_bic_knots_total", diag.bic_knots_total)
    st_numscalar("__spmata_bic_knots_max", diag.bic_knots_max)
    st_numscalar("__spmata_rank_valid", diag.rank_valid)
    st_numscalar("__spmata_nuisance_rank_deficient", diag.nuisance_rank_deficient)
    st_numscalar("__spmata_post_condition", diag.post_condition)
    st_numscalar("__spmata_final_rank", diag.final_rank)
    st_numscalar("__spmata_final_condition", diag.final_condition)
    st_numscalar("__spmata_finite", diag.finite)
}
end
