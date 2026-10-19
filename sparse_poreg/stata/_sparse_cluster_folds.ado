*! _sparse_cluster_folds 0.8.4-rc7  16jul2026
program _sparse_cluster_folds, rclass
    version 16.0
    syntax newvarlist(min=1 max=2) [if] [in], ///
        CLUSTERS(varlist numeric min=1 max=2) ///
        XFOLDS(numlist integer min=1 max=2) RSEED(integer)

    local nc : word count `clusters'
    local nv : word count `varlist'
    local nk : word count `xfolds'
    if `nv' != `nc' | !inlist(`nk', 1, `nc') {
        di as err "specify one new fold variable per cluster variable"
        exit 198
    }
    if `nc' == 2 & `nk' == 1 local xfolds `xfolds' `xfolds'
    if `rseed' < 0 | `rseed' > 2147483647 {
        di as err "rseed() must be between 0 and 2147483647"
        exit 198
    }

    marksample touse, novarlist
    quietly markout `touse' `clusters'
    forvalues j = 1/`nc' {
        local c : word `j' of `clusters'
        local k : word `j' of `xfolds'
        if `k' < 2 {
            di as err "xfolds() must contain integers greater than one"
            exit 198
        }
        tempvar check
        quietly egen byte `check' = tag(`c') if `touse'
        quietly count if `check' == 1
        local g`j' = r(N)
        quietly drop `check'
        if `g`j'' < `k' {
            di as err "cluster variable `c' has fewer groups than folds"
            exit 498
        }
        local k`j' = `k'
    }

    tempvar row tag u rank
    quietly gen long `row' = _n
    set seed `rseed'
    forvalues j = 1/`nc' {
        local c : word `j' of `clusters'
        local f : word `j' of `varlist'
        quietly sort `c' `row'
        quietly egen byte `tag' = tag(`c') if `touse'
        quietly gen double `u' = runiform() if `tag' == 1
        quietly sort `u' `c'
        quietly gen long `rank' = sum(`tag' == 1)
        quietly gen long `f' = mod(`rank' - 1, `k`j'') + 1 if `tag'
        quietly bysort `c' (`f'): replace `f' = `f'[1] if `touse'
        quietly drop `tag' `u' `rank'
    }
    quietly sort `row'

    return scalar N_dimensions = `nc'
    forvalues j = 1/`nc' {
        return scalar N_cluster`j' = `g`j''
        return scalar N_folds`j' = `k`j''
    }
    return local clustervars "`clusters'"
    return local foldvars "`varlist'"
end
