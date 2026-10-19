*! sparse_poreg 0.8.4-rc14
program sparse_poreg, eclass
    version 16.0
    local sparse_cmdline `"sparse_poreg `0'"'
    syntax varlist(fv min=2 numeric) [if] [in], ///
        [ABSorb(varlist) CONTROLS(string asis) SELection(string) VCE(string asis) ///
        BACKend(string) STRICT(string) DENSEMax(integer 20000) ///
        MAXIter(integer 4000) COEFTol(real 1e-3) ACCuracy(string) KERNEL(string) ///
        SQRTlasso SEMI MISSINGOK FRONTIEROK RSEED(integer -1) *]

    marksample touse
    gettoken depvar treat : varlist

    local requested_backend = lower(trim("`backend'"))
    if "`requested_backend'" == "" local requested_backend auto
    if "`requested_backend'" == "stata" local requested_backend official
    local run_backend `requested_backend'
    if !inlist("`run_backend'", "auto", "official", "mata") {
        di as err "backend() must be auto, official, or mata"
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
    if !inlist("`selection'", "plugin", "bic", "cv", "adaptive") {
        di as err "selection() must be plugin, bic, cv, or adaptive"
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
    local kernel = lower(trim("`kernel'"))
    if "`kernel'" == "" local kernel auto
    if !inlist("`kernel'", "auto", "legacy", "indexed") {
        di as err "kernel() must be auto, legacy, or indexed"
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
    if `rseed' < -1 | `rseed' > 2147483647 {
        di as err "rseed() must be between 0 and 2147483647"
        exit 198
    }
    local seed_only
    if `rseed' >= 0 local seed_only rseed()

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
        capture unab clustvars : `vcerest'
        if _rc {
            di as err "cluster variable not found"
            exit 111
        }
        local nclust : word count `clustvars'
        if !inrange(`nclust', 1, 2) {
            di as err "vce() must be robust or cluster clustvar [clustvar2]"
            exit 198
        }
        local unique : list uniq clustvars
        local nunique : word count `unique'
        if `nunique' != `nclust' {
            di as err "cluster variables must be distinct"
            exit 198
        }
        local clustvars : list sort clustvars
        gettoken clustvar clustrest : clustvars
        gettoken clustvar2 : clustrest
        local vcetype cluster
        local vce cluster `clustvars'
    }
    else {
        di as err "vce() must be robust or cluster clustvar [clustvar2]"
        exit 198
    }
    local po_multiway = "`vcetype'" == "cluster" & `nclust' == 2
    if `po_multiway' & !("`run_backend'" == "mata" & "`frontierok'" != "") {
        di as err "two-way clustering requires backend(mata) frontierok"
        exit 198
    }
    if "`frontierok'" != "" & !("`run_backend'" == "mata" & ///
        "`vcetype'" == "cluster" & ///
        ("`selection'" == "plugin" | `po_multiway')) {
        di as err "frontierok requires a supported backend(mata) clustered fit"
        exit 198
    }

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
    }
    else _sparse_official_route, touse(`touse') densemax(`densemax') limitsonly
    local route_cap `s(official_cap)'
    local route_maxvar `s(maxvar_headroom)'
    local route_reason explicit_`requested_backend'
    if "`run_backend'" == "auto" {
        if `auto_official' local run_backend official
        else local run_backend mata
        local route_reason auto_`run_backend'_`route_base'
    }
    if "`run_backend'" == "official" {
        if !`numeric_absorb' {
            di as err "backend(official) requires numeric absorb() variables"
            exit 198
        }
        local seedopt
        if `rseed' >= 0 local seedopt rseed(`rseed')
        poregress `depvar' `treat' `if' `in', controls(`dense_controls') ///
            selection(`selspec') vce(`vce') `sqrtlasso' `semi' ///
            `missingok' `seedopt' `options'
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
        if trim(`"`sqrtlasso' `semi' `missingok' `options'"') != "" {
            di as err "sqrtlasso, semi, missingok, and official passthrough options require backend(official)"
            exit 198
        }
        local native_markvars `absorb' `controls'
        if "`clustvars'" != "" local native_markvars `native_markvars' `clustvars'
        quietly markout `touse' `native_markvars'
        quietly count if `touse'
        if r(N) == 0 {
            di as err "no complete observations"
            exit 2000
        }
        local xcheckopt
        local xcheck_seconds .
        if "`strict'" == "stata" & "`selection'" == "bic" {
            if `nclust' > 1 {
                di as err "strict(stata) xcheck is unavailable for two-way clustering"
                di as err "use strict(production) for native multiway inference"
                exit 198
            }
            local xvars
            foreach a of local absorb {
                local xvars `xvars' i.`a'
            }
            local xvars `xvars' `controls'
            local xcluster
            if `nclust' == 1 local xcluster cluster(`clustvar')
            local xstart = clock("`c(current_date)' `c(current_time)'", "DMY hms")
            capture quietly lasso linear `depvar' `xvars' if `touse', ///
                selection(none) lambdas(1000000) `xcluster' nolog
            local xrc = _rc
            if `xrc' {
                di as err "Stata xcheck prepass failed (return code `xrc')"
                di as err "increase maxvar or use strict(production)"
                exit `xrc'
            }
            local xcheckopt xcheck
            local xstop = clock("`c(current_date)' `c(current_time)'", "DMY hms")
            local xcheck_seconds = (`xstop' - `xstart') / 1000
        }
        if "`selection'" == "plugin" & "`vcetype'" == "cluster" {
            local guard_controls
            if "`controls'" != "" local guard_controls controls(`controls')
            _sparse_native_guard, touse(`touse') absorb(`absorb') ///
                `guard_controls' cluster(`clustvars')
            if r(min_abs_count) <= 1 {
                di as err "native plugin/cluster does not support singleton absorbed levels"
                di as err "use selection(bic), merge/drop singleton levels, or backend(official)"
                exit 498
            }
        }
        local accopt
        if "`accuracy'" != "" local accopt accuracy(`accuracy')
        capture noisily _sparse_prl_xpo_mata `depvar' `treat' if `touse', ///
            cmd(sparse_poreg) absorb(`absorb') controls(`controls') ///
            selection(`selection') bicmode(`bicmode') vce(`vce') ///
            strict(`strict') maxiter(`maxiter') ///
            coeftol(`coeftol') kernel(`kernel') `xcheckopt' `accopt'
        local native_rc = _rc
        if `native_rc' exit `native_rc'
        _sparse_route_post, requested(`requested_backend') ///
            reason(`route_reason') p(`=e(sparse_features_raw)') n(`=e(N)') ///
            densemax(`densemax') officialcap(`route_cap') ///
            maxvarheadroom(`route_maxvar')
        ereturn local sparse_backend "mata"
        ereturn local sparse_strict "`strict'"
        ereturn local sparse_absorb "`absorb'"
        ereturn local sparse_version "0.8.4-rc14"
        ereturn scalar xcheck_seconds = `xcheck_seconds'
        if `rseed' >= 0 ereturn scalar rseed = `rseed'
        if `po_multiway' {
            ereturn local inference_contract "multiway_cluster"
            ereturn local nuisance_weighting "observation"
            ereturn local vcetype "Multiway cluster-robust"
            ereturn local clustvar "`clustvar'"
            ereturn local clustvars "`clustvars'"
            ereturn scalar multiway_cluster = 1
        }
        if "`requested_backend'" == "auto" di as txt ///
            "auto route: Mata; P=" %12.0gc e(sparse_route_p) ///
            " N=" %12.0gc e(sparse_route_n) " P/N=" %9.4f e(sparse_route_pn)
        ereturn display
        exit
    }

end
