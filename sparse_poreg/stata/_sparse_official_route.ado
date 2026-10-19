*! _sparse_official_route 0.8.4-rc14  22jul2026
program _sparse_official_route, sclass
    version 16.0
    syntax, TOUSE(varname) DENSEMAX(integer) ///
        [ABSORB(varlist) CONTROLS(string asis) LIMITSONLY]

    quietly count if `touse'
    local route_n = r(N)
    quietly describe
    local maxvar_headroom = max(0, c(maxvar) - r(k) - 50)
    local cap = min(`densemax', `maxvar_headroom')

    if "`limitsonly'" != "" {
        sreturn local dense_controls
        sreturn local numeric_absorb .
        sreturn local expand_ok 0
        sreturn local dense_n .
        sreturn local expanded_p .
        sreturn local route_pn .
        sreturn local route_n `route_n'
        sreturn local densemax `densemax'
        sreturn local maxvar_headroom `maxvar_headroom'
        sreturn local official_cap `cap'
        sreturn local official_headroom .
        sreturn local route_reason limits_only
        sreturn local auto_official 0
        exit
    }

    local absorb_fv
    local numeric_absorb 1
    foreach a of local absorb {
        capture confirm numeric variable `a'
        if _rc local numeric_absorb 0
        local absorb_fv `absorb_fv' ibn.`a'
    }

    local combined `"`controls' `absorb_fv'"'
    if trim(`"`combined'"') == "" {
        di as err "absorb() or controls() is required"
        exit 198
    }

    _lasso_parse_controls `combined'
    local always `s(ainclude)'
    local other `s(othervars)'
    local indepvars `s(indepvars)'
    if trim(`"`always'"') != "" local dense_controls `"(`always') `other'"'
    else local dense_controls `"`other'"'

    local dense_n = `cap' + 1
    local expand_ok 0
    if `numeric_absorb' {
        capture quietly fvexpand `indepvars' if `touse'
        if !_rc {
            local expanded `r(varlist)'
            local expanded : list uniq expanded
            local dense_n : word count `expanded'
            local expand_ok 1
            if `dense_n' <= `cap' {
                preserve
                    quietly keep if `touse'
                    tempvar complete
                    quietly gen byte `complete' = 1
                    capture quietly fvrevar `indepvars'
                    if _rc local expand_ok 0
                    else {
                        quietly markout `complete' `r(varlist)'
                        capture quietly fvexpand `indepvars' if `complete'
                        if _rc local expand_ok 0
                        else {
                            local expanded `r(varlist)'
                            local expanded : list uniq expanded
                            local dense_n : word count `expanded'
                        }
                    }
                restore
            }
        }
    }

    local auto_official = (`numeric_absorb' & `expand_ok' & `dense_n' <= `cap')
    if !`numeric_absorb' local route_reason non_numeric_absorb
    else if !`expand_ok' local route_reason expansion_failed
    else if `dense_n' > `cap' local route_reason expanded_p_exceeds_cap
    else local route_reason within_cap
    local expanded_p = cond(`expand_ok', `dense_n', .)
    local route_pn = cond(`expand_ok' & `route_n' > 0, `dense_n'/`route_n', .)
    local official_headroom = `cap' - `dense_n'

    sreturn local dense_controls `"`dense_controls'"'
    sreturn local numeric_absorb `numeric_absorb'
    sreturn local expand_ok `expand_ok'
    sreturn local dense_n `dense_n'
    sreturn local expanded_p `expanded_p'
    sreturn local route_pn `route_pn'
    sreturn local route_n `route_n'
    sreturn local densemax `densemax'
    sreturn local maxvar_headroom `maxvar_headroom'
    sreturn local official_cap `cap'
    sreturn local official_headroom `official_headroom'
    sreturn local route_reason `route_reason'
    sreturn local auto_official `auto_official'
end
