*! sparse_xporeg 0.8.4-rc14
program sparse_xporeg, eclass
    version 16.0
    local sparse_cmdline `"sparse_xporeg `0'"'
    syntax varlist(fv min=2 numeric) [if] [in], ///
        [ABSorb(varlist) CONTROLS(string asis) SELection(string) VCE(string asis) ///
        XFOLDS(integer 10) RSEED(integer 1103) BACKend(string) STRICT(string) ///
        DENSEMax(integer 20000) PYthon(string) MAXIter(integer 4000) ///
        COEFTol(real 1e-3) ACCuracy(string) KERNEL(string) SQRTlasso ///
        SEMI MISSINGOK FRONTIEROK FOLDCLUSTER(varlist numeric min=1 max=2) ///
        LEARNER(string) LASSOCVFolds(integer 5) ///
        LASSOCVAlphas(integer 100) LASSOCVEps(real 0.001) ///
        LASSOCVMaxiter(integer 5000) LASSOCVTol(real 1e-4) *]

    marksample touse
    gettoken depvar treat : varlist

    local requested_backend = lower(trim("`backend'"))
    if "`requested_backend'" == "" local requested_backend auto
    if "`requested_backend'" == "stata" local requested_backend official
    local run_backend `requested_backend'
    if !inlist("`run_backend'", "auto", "official", "mata", "python") {
        di as err "backend() must be auto, official, mata, or python"
        exit 198
    }
    local selspec = lower(trim(`"`selection'"'))
    if `"`selspec'"' == "" local selspec plugin
    local comma = strpos(`"`selspec'"', ",")
    local selopts
    if `comma' {
        local selection = trim(substr(`"`selspec'"', 1, `comma' - 1))
        local selopts = trim(substr(`"`selspec'"', `comma' + 1, strlen(`"`selspec'"')))
    }
    else local selection `"`selspec'"'
    if !inlist("`selection'", "plugin", "bic", "cv", "adaptive", "lassocv") {
        di as err "selection() must be plugin, bic, cv, adaptive, or lassocv"
        exit 198
    }
    local bicmode = cond("`selection'" == "bic", "runningmin", "none")
    if `"`selopts'"' != "" {
        if "`selection'" != "bic" | "`selopts'" != "gridminok" {
            di as err "only selection(bic, gridminok) has a native suboption"
            exit 198
        }
        local bicmode gridminok
    }
    else if `comma' {
        di as err "selection() suboption is missing"
        exit 198
    }
    if "`selection'" == "lassocv" & "`run_backend'" != "python" {
        di as err "selection(lassocv) requires backend(python)"
        exit 198
    }
    if "`run_backend'" == "python" & "`selection'" != "lassocv" {
        di as err "backend(python) currently implements selection(lassocv) only"
        exit 198
    }
    if "`sqrtlasso'" != "" & "`selection'" == "adaptive" {
        di as err "selection(adaptive) cannot be combined with sqrtlasso"
        exit 198
    }

    local strict = lower(trim("`strict'"))
    if "`strict'" == "" local strict stata
    if !inlist("`strict'", "stata", "production") {
        di as err "strict() must be stata or production"
        exit 198
    }
    local accuracy = lower(trim("`accuracy'"))
    if "`accuracy'" == "" local accuracy high
    if !inlist("`accuracy'", "high", "fast") {
        di as err "accuracy() must be high or fast"
        exit 198
    }
    if `xfolds' < 2 {
        di as err "xfolds() must be at least 2"
        exit 198
    }
    local kernel = lower(trim("`kernel'"))
    if "`kernel'" == "" local kernel auto
    if !inlist("`kernel'", "auto", "legacy", "indexed") {
        di as err "kernel() must be auto, legacy, or indexed"
        exit 198
    }

    if `rseed' < 0 | `rseed' > 2147483647 {
        di as err "rseed() must be between 0 and 2147483647"
        exit 198
    }
    if `densemax' < 1 {
        di as err "densemax() must be positive"
        exit 198
    }
    if `maxiter' < 1 | `coeftol' <= 0 {
        di as err "maxiter() and coeftol() must be positive"
        exit 198
    }
    if `lassocvfolds' < 2 | `lassocvalphas' < 2 | ///
        `lassocveps' <= 0 | `lassocveps' > 1 | ///
        `lassocvmaxiter' < 1 | `lassocvtol' <= 0 {
        di as err "invalid Python lassocv tuning option"
        exit 198
    }

    local vce = trim(`"`vce'"')
    if "`vce'" == "" local vce robust
    gettoken vcetype vcerest : vce
    local vcetype = lower("`vcetype'")
    if "`vcetype'" == "robust" & trim("`vcerest'") == "" {
        local vcetype robust
        local clustvar
        local clustvar2
        local clustvars
        local nclust 0
        local vce robust
    }
    else if "`vcetype'" == "cluster" {
        local clustvars = trim("`vcerest'")
        local nclust : word count `clustvars'
        if !inrange(`nclust', 1, 2) {
            di as err "vce() must be robust or cluster clustvar [clustvar2]"
            exit 198
        }
        capture unab clustvars : `clustvars'
        local clust_rc = _rc
        if `clust_rc' {
            di as err "cluster variable not found"
            exit 111
        }
        local nclust : word count `clustvars'
        if !inrange(`nclust', 1, 2) {
            di as err "vce() accepts at most two cluster variables"
            exit 198
        }
        local unique_clustvars : list uniq clustvars
        local n_unique : word count `unique_clustvars'
        if `n_unique' != `nclust' {
            di as err "cluster variables must be distinct"
            exit 198
        }
        if `nclust' == 2 local clustvars : list sort clustvars
        gettoken clustvar clustrest : clustvars
        gettoken clustvar2 : clustrest
        local vcetype cluster
        local vce cluster `clustvars'
    }
    else {
        di as err "vce() must be robust or cluster clustvar [clustvar2]"
        exit 198
    }

    local clusterfolds
    if "`foldcluster'" != "" {
        if "`vcetype'" != "cluster" {
            di as err "foldcluster() requires vce(cluster ...)"
            exit 198
        }
        local nfoldcluster : word count `foldcluster'
        local fold_only : list foldcluster - clustvars
        local vce_only : list clustvars - foldcluster
        if `nfoldcluster' != `nclust' | "`fold_only'" != "" | "`vce_only'" != "" {
            di as err "foldcluster() and vce(cluster ...) must contain the same variables"
            exit 198
        }
        local foldcluster `clustvars'
        local clusterfolds clusterfolds
    }
    if `nclust' == 2 & "`clusterfolds'" == "" {
        di as err "two-way vce(cluster ...) requires foldcluster()"
        exit 198
    }
    if "`clusterfolds'" != "" & "`run_backend'" == "mata" & ///
        "`frontierok'" == "" {
        di as err "Mata foldcluster() requires frontierok"
        exit 198
    }
    if "`clusterfolds'" != "" & !inlist("`run_backend'", "mata", "python") {
        di as err "foldcluster() requires explicit backend(mata) or backend(python)"
        exit 198
    }
    if "`run_backend'" == "python" & "`vcetype'" == "cluster" & ///
        "`clusterfolds'" == "" {
        di as err "clustered backend(python) requires matching foldcluster()"
        exit 198
    }
    if "`frontierok'" != "" & "`run_backend'" != "mata" {
        di as err "frontierok requires explicit backend(mata)"
        exit 198
    }

    local fold_markvars `depvar' `treat' `absorb' `controls'
    if "`clustvars'" != "" local fold_markvars `fold_markvars' `clustvars'

    local route_absorb
    local route_controls
    if "`absorb'" != "" local route_absorb absorb(`absorb')
    if trim(`"`controls'"') != "" local route_controls controls(`controls')
    if inlist("`run_backend'", "auto", "official") {
        _sparse_official_route, touse(`touse') densemax(`densemax') ///
            `route_absorb' `route_controls'
        local dense_controls `"`s(dense_controls)'"'
        local numeric_absorb `s(numeric_absorb)'
        local auto_official `s(auto_official)'
        local route_base `s(route_reason)'
        local route_diag_p `s(dense_n)'
        local route_diag_n `s(route_n)'
    }
    else _sparse_official_route, touse(`touse') densemax(`densemax') limitsonly
    local route_cap `s(official_cap)'
    local route_maxvar `s(maxvar_headroom)'
    local route_reason explicit_`requested_backend'
    if "`run_backend'" == "auto" {
        if !`auto_official' {
            di as err "auto cannot route large XPO to the experimental Mata engine"
            di as err "expanded P=`route_diag_p', N=`route_diag_n', official cap=`route_cap'"
            di as err "use backend(mata) frontierok, or backend(python) selection(lassocv)"
            exit 198
        }
        local run_backend official
        local route_reason auto_official_`route_base'
    }
    if "`run_backend'" == "official" {
        if !`numeric_absorb' {
            di as err "backend(official) requires numeric absorb() variables"
            exit 198
        }
        xporegress `depvar' `treat' `if' `in', controls(`dense_controls') ///
            selection(`selspec') vce(`vce') xfolds(`xfolds') ///
            rseed(`rseed') `sqrtlasso' `semi' `missingok' `options'
        _sparse_route_post, requested(`requested_backend') ///
            reason(`route_reason') p(`=e(k_controls)') n(`=e(N)') ///
            densemax(`densemax') officialcap(`route_cap') ///
            maxvarheadroom(`route_maxvar')
        ereturn local sparse_backend "official"
        ereturn local sparse_backend_alias "stata"
        ereturn local sparse_absorb "`absorb'"
        ereturn local sparse_version "0.8.4-rc14"
        if "`requested_backend'" == "auto" di as txt ///
            "auto route: official; P=" %12.0gc e(sparse_route_p) ///
            " N=" %12.0gc e(sparse_route_n) " P/N=" %9.4f e(sparse_route_pn)
        exit
    }
    if "`run_backend'" == "mata" {
        if !inlist("`selection'", "plugin", "bic") {
            di as err "backend(mata) supports selection(plugin|bic) only"
            exit 198
        }
        if "`absorb'" == "" {
            di as err "absorb() is required with backend(mata)"
            exit 198
        }
        foreach a of local absorb {
            capture confirm numeric variable `a'
            if _rc {
                di as err "backend(mata) requires numeric absorb() variables"
                exit 109
            }
        }
        if "`vcetype'" == "cluster" {
            foreach c of local clustvars {
                capture confirm numeric variable `c'
                if _rc {
                    di as err "backend(mata) requires numeric cluster variables"
                    exit 109
                }
            }
        }
        if trim(`"`controls'"') != "" {
            if strpos(`"`controls'"', "(") | strpos(`"`controls'"', ")") | ///
                strpos(`"`controls'"', ".") | strpos(`"`controls'"', "#") {
                di as err "backend(mata) controls() must be a plain numeric varlist"
                exit 198
            }
            capture unab mata_controls : `controls'
            if _rc {
                di as err "backend(mata) controls() must be a plain numeric varlist"
                exit 198
            }
            foreach c of local mata_controls {
                capture confirm numeric variable `c'
                if _rc {
                    di as err "backend(mata) controls() must be a plain numeric varlist"
                    exit 109
                }
            }
            local controls `mata_controls'
        }
        if "`frontierok'" == "" {
            di as err "native XPO is experimental; specify frontierok to opt in"
            di as err "for certified inference use backend(official) or a feasible backend(auto) design"
            exit 198
        }
        if trim(`"`sqrtlasso' `semi' `missingok' `options'"') != "" {
            di as err "sqrtlasso, semi, missingok, resample(), technique(), and official passthrough options require backend(official)"
            exit 198
        }
        local guard_controls
        local guard_cluster
        if "`controls'" != "" local guard_controls controls(`controls')
        if "`clustvars'" != "" local guard_cluster cluster(`clustvars')
        local accopt
        if "`accuracy'" != "" local accopt accuracy(`accuracy')
        local xcheckopt
        local xcheck_seconds .
        local xcheck_prefix
        tempvar fold
        if "`clusterfolds'" != "" & `nclust' == 2 tempvar fold2
        preserve
            keep if `touse'
            keep `depvar' `treat' `absorb' `controls' `clustvars'
            tempvar fold_touse
            gen byte `fold_touse' = 1
            markout `fold_touse' `fold_markvars'
            keep if `fold_touse'
            quietly count
            local n_complete = r(N)
            if `n_complete' == 0 {
                restore
                di as err "no complete observations"
                exit 2000
            }
            if `xfolds' > `n_complete' {
                restore
                di as err "xfolds() cannot exceed the complete-case sample size"
                exit 198
            }

            local guard_folds `fold'
            local fold2opt
            if "`clusterfolds'" != "" {
                if `nclust' == 1 {
                    _sparse_cluster_folds `fold', clusters(`clustvars') ///
                        xfolds(`xfolds') rseed(`rseed')
                }
                else {
                    _sparse_cluster_folds `fold' `fold2', ///
                        clusters(`clustvars') xfolds(`xfolds') rseed(`rseed')
                    local guard_folds `fold' `fold2'
                    local fold2opt foldvar2(`fold2')
                }

                forvalues k1 = 1/`xfolds' {
                    if `nclust' == 1 {
                        quietly count if `fold' == `k1'
                        if !r(N) {
                            restore
                            di as err "cluster score fold `k1' is empty"
                            exit 498
                        }
                        quietly count if `fold' != `k1'
                        if !r(N) {
                            restore
                            di as err "cluster training fold `k1' is empty"
                            exit 498
                        }
                    }
                    else {
                        forvalues k2 = 1/`xfolds' {
                            quietly count if `fold' == `k1' & `fold2' == `k2'
                            if !r(N) {
                                restore
                                di as err "two-way cluster score fold (`k1',`k2') is empty"
                                exit 498
                            }
                            quietly count if `fold' != `k1' & `fold2' != `k2'
                            if !r(N) {
                                restore
                                di as err "two-way cluster training fold (`k1',`k2') is empty"
                                exit 498
                            }
                        }
                    }
                }
            }
            else {
                set seed `rseed'
                _dml_part_data, gr(`fold') xfolds(`xfolds')
            }
            _sparse_native_guard, touse(`fold_touse') absorb(`absorb') ///
                `guard_controls' `guard_cluster' foldvar(`guard_folds')
            local min_abs = r(min_abs_count)
            local min_train = r(min_train_count)
            local support_relaxed = (`min_train' <= 1)
            if `support_relaxed' & "`strict'" == "production" {
                restore
                di as err "strict(production) requires at least two observations per absorbed level in every training fold"
                di as err "merge/drop rare levels or use strict(stata) for official-compatible zero/singleton handling"
                exit 498
            }
            if `support_relaxed' {
                di as txt "warning: strict(stata) is using official-compatible absent/singleton training-level handling"
                di as txt "use strict(production) to require minimum training support of two"
            }
            if "`strict'" == "stata" & "`clusterfolds'" == "" {
                local xvars
                foreach a of local absorb {
                    local xvars `xvars' i.`a'
                }
                local xvars `xvars' `controls'
                local xcluster
                if `nclust' == 1 local xcluster cluster(`clustvar')
                tempname xprefix
                local xsuffix = substr("`xprefix'", 3, .)
                local xcheck_prefix "SPXCHK`xsuffix'_"
                local xstart = clock("`c(current_date)' `c(current_time)'", "DMY hms")
                forvalues k = 1/`xfolds' {
                    capture quietly lasso linear `depvar' `xvars' if `fold' != `k', ///
                        selection(none) lambdas(1000000) `xcluster' nolog
                    local xrc = _rc
                    if `xrc' {
                        forvalues h = 1/`k' {
                            capture macro drop `xcheck_prefix'`h'
                        }
                        restore
                        di as err "Stata fold `k' xcheck prepass failed (return code `xrc')"
                        di as err "increase maxvar or use strict(production)"
                        exit `xrc'
                    }
                    local xspec `"`e(allvars)'"'
                    mata: st_global("`xcheck_prefix'`k'", st_local("xspec"))
                }
                local xstop = clock("`c(current_date)' `c(current_time)'", "DMY hms")
                local xcheck_seconds = (`xstop' - `xstart') / 1000
                local xcheckopt xcheck xcheckprefix(`xcheck_prefix')
            }
            local kktopt
            if "`selection'" == "plugin" & "`vcetype'" == "cluster" & ///
                "`frontierok'" != "" & "`strict'" != "production" ///
                local kktopt allowkktfail
            capture noisily _sparse_prl_xpo_mata `depvar' `treat', ///
                cmd(sparse_xporeg) ///
                absorb(`absorb') controls(`controls') selection(`selection') ///
                bicmode(`bicmode') ///
                vce(`vce') strict(`strict') foldvar(`fold') `fold2opt' ///
                xfolds(`xfolds') rseed(`rseed') `clusterfolds' ///
                maxiter(`maxiter') coeftol(`coeftol') kernel(`kernel') ///
                `xcheckopt' `accopt' `kktopt'
            local native_rc = _rc
            if "`xcheck_prefix'" != "" {
                forvalues k = 1/`xfolds' {
                    capture macro drop `xcheck_prefix'`k'
                }
            }
            if `native_rc' {
                restore
                exit `native_rc'
            }
        restore
        _sparse_route_post, requested(`requested_backend') ///
            reason(`route_reason') p(`=e(sparse_features_raw)') n(`=e(N)') ///
            densemax(`densemax') officialcap(`route_cap') ///
            maxvarheadroom(`route_maxvar')
        if "`clusterfolds'" != "" {
            ereturn local inference_contract "cluster_dml"
            ereturn local crossfit_unit "cluster"
            ereturn local foldclustvars "`clustvars'"
            ereturn scalar cluster_disjoint = 1
            ereturn scalar n_xfolds_per_cluster = `xfolds'
            ereturn scalar n_cluster_dims = `nclust'
            ereturn scalar n_score_folds = `xfolds'^`nclust'
        }
        ereturn local sparse_backend "mata"
        ereturn local sparse_strict "`strict'"
        ereturn local sparse_absorb "`absorb'"
        ereturn local sparse_version "0.8.4-rc14"
        ereturn scalar min_abs_support = `min_abs'
        ereturn scalar min_train_support = `min_train'
        ereturn scalar training_support_safe = (`min_train' > 1)
        ereturn scalar training_support_relaxed = `support_relaxed'
        ereturn scalar xcheck_seconds = `xcheck_seconds'
        ereturn display
        exit
    }

    if trim(`"`sqrtlasso' `semi' `missingok' `frontierok' `options'"') != "" {
        di as err "backend(python) is experimental and accepts lassocv options only"
        exit 198
    }
    local learner = lower(trim("`learner'"))
    if "`learner'" == "" local learner lassocv
    if "`learner'" != "lassocv" {
        di as err "backend(python) currently implements learner(lassocv) only"
        exit 198
    }

    if trim(`"`python'"') == "" _sparse_python_paths
    else _sparse_python_paths, python("`python'")
    local pyexec `"`s(python_exec)'"'
    local pyscript `"`s(python_script)'"'
    tempfile input postdo
    tempvar fold

    if strpos(`"`treat'"', ".") | strpos(`"`treat'"', "#") {
        di as err "backend(python) requires plain numeric treatment variables"
        exit 198
    }
    if strpos(`"`controls'"', "(") | strpos(`"`controls'"', ")") | ///
            strpos(`"`controls'"', "#") {
        di as err "backend(python) does not support always-in or interaction controls"
        exit 198
    }
    fvrevar `treat' if `touse'
    local py_treat `r(varlist)'
    local py_absorb `absorb'
    local py_controls
    foreach c of local controls {
        if substr("`c'", 1, 2) == "c." {
            local c = substr("`c'", 3, .)
        }
        local raw = subinstr("`c'", "ibn.", "", .)
        local raw = subinstr("`raw'", "i.", "", .)
        if "`raw'" != "`c'" {
            local py_absorb `py_absorb' `raw'
        }
        else {
            if strpos("`c'", ".") {
                di as err "unsupported Python control term: `c'"
                exit 198
            }
            capture unab cvars : `c'
            if _rc {
                di as err "Python control variable not found: `c'"
                exit 111
            }
            local py_controls `py_controls' `cvars'
        }
    }
    local py_absorb : list uniq py_absorb
    local py_controls : list uniq py_controls
    foreach v of local py_treat {
        capture confirm numeric variable `v'
        if _rc {
            di as err "backend(python) requires numeric treatment variables"
            exit 109
        }
    }
    foreach v of local py_absorb {
        capture confirm numeric variable `v'
        if _rc {
            di as err "backend(python) requires numeric absorb variables"
            exit 109
        }
    }
    foreach v of local py_controls {
        capture confirm numeric variable `v'
        if _rc {
            di as err "backend(python) requires numeric controls"
            exit 109
        }
    }
    local savevars `depvar' `py_treat' `py_absorb' `py_controls' `clustvars'
    local savevars : list uniq savevars
    tempvar py_sample
    quietly gen byte `py_sample' = `touse'
    quietly markout `py_sample' `savevars', strok
    quietly count if `py_sample'
    if r(N) == 0 {
        di as err "no complete observations"
        exit 2000
    }
    local py_foldvar
    preserve
        keep if `py_sample'
        keep `savevars'
        if "`vcetype'" == "robust" {
            set seed `rseed'
            _dml_part_data, gr(`fold') xfolds(`xfolds')
            local py_foldvar `fold'
        }
        save "`input'", replace
    restore

    local shellcmd `""`pyexec'" "`pyscript'" --cmd sparse_xporeg --input "`input'" --postdo "`postdo'" --esample-var "`py_sample'" --depvar "`depvar'" --treat "`py_treat'" --treat-labels "`py_treat'" --absorb "`py_absorb'" --controls "`py_controls'" --selection "`selection'" --vce "`vcetype'" --cluster "`clustvars'" --fold-var "`py_foldvar'" --xfolds `xfolds' --rseed `rseed' --max-iter `maxiter' --coef-tol `coeftol' --learner "`learner'" --lassocv-folds `lassocvfolds' --lassocv-alphas `lassocvalphas' --lassocv-eps `lassocveps' --lassocv-maxiter `lassocvmaxiter' --lassocv-tol `lassocvtol'"'
    tempfile launchbat launchps1 rcfile
    local launchbat "`launchbat'.bat"
    local launchps1 "`launchps1'.ps1"
    local rcfile "`rcfile'.rc"
    capture erase "`rcfile'"
    capture erase "`postdo'"
    file open __bat using "`launchbat'", write text replace
    file write __bat "@echo off" _n
    file write __bat `"`shellcmd'"' _n
    file write __bat `"echo %ERRORLEVEL% > "`rcfile'""' _n
    file close __bat
    file open __ps1 using "`launchps1'", write text replace
    file write __ps1 `"& '`launchbat''"' _n
    file close __ps1
    capture winexec powershell.exe -NoProfile -ExecutionPolicy Bypass -WindowStyle Hidden -File "`launchps1'"
    if _rc capture quietly shell powershell.exe -NoProfile -ExecutionPolicy Bypass -WindowStyle Hidden -File "`launchps1'"
    local py_rc = _rc
    if !`py_rc' {
        local spins 0
        while `spins' < 172800 {
            capture confirm file "`rcfile'"
            if !_rc continue, break
            sleep 500
            local ++spins
        }
        capture confirm file "`rcfile'"
        if _rc local py_rc = 498
        else {
            file open __rc using "`rcfile'", read text
            file read __rc line
            file close __rc
            local py_rc = real(trim("`line'"))
            if missing(`py_rc') local py_rc = 498
        }
    }
    capture confirm file "`postdo'"
    local post_missing = _rc
    if `py_rc' | `post_missing' {
        di as err "sparse_xporeg backend failed"
        if `py_rc' di as err "python rc=`py_rc'"
        capture noisily type "`postdo'.py.log"
        exit 498
    }

    include "`postdo'"
    _sparse_route_post, requested(`requested_backend') ///
        reason(`route_reason') p(`=e(sparse_features)') n(`=e(N)') ///
        densemax(`densemax') officialcap(`route_cap') ///
        maxvarheadroom(`route_maxvar')
    ereturn local sparse_backend "python"
    ereturn local sparse_absorb "`absorb'"
    ereturn local sparse_version "0.8.4-rc14"
    ereturn local python_exec `"`pyexec'"'
    ereturn local python_executable `"`pyexec'"'
    ereturn local python_script `"`pyscript'"'
    ereturn local cmdline `"`sparse_cmdline'"'
    ereturn display
end
