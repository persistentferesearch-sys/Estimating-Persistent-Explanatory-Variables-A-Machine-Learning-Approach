*! _sparse_native_guard 0.8.4-rc7  16jul2026
program _sparse_native_guard, rclass
    version 16.0
    syntax, TOUSE(varname) ABSORB(varlist) ///
        [CONTROLS(varlist) CLUSTER(varlist numeric max=2) ///
        FOLDVAR(varlist numeric min=1 max=2)]

    preserve
        quietly keep if `touse'
        tempvar complete count within1 excluded train
        quietly gen byte `complete' = 1
        quietly markout `complete' `absorb' `controls' `cluster' `foldvar'
        quietly keep if `complete'
        quietly count
        local n_complete = r(N)
        local mincount .
        local mintrain .
        local nf : word count `foldvar'
        if `n_complete' {
            foreach a of local absorb {
                quietly bysort `a': gen long `count' = _N
                quietly summarize `count', meanonly
                if missing(`mincount') | r(min) < `mincount' local mincount = r(min)
                if `nf' == 1 {
                    quietly bysort `a' `foldvar': gen long `within1' = _N
                    quietly gen long `train' = `count' - `within1'
                    quietly summarize `train', meanonly
                    if missing(`mintrain') | r(min) < `mintrain' local mintrain = r(min)
                    quietly drop `within1' `train'
                }
                else if `nf' == 2 {
                    gettoken fold1 fold2 : foldvar
                    quietly levelsof `fold1', local(folds1)
                    quietly levelsof `fold2', local(folds2)
                    foreach f1 of local folds1 {
                        foreach f2 of local folds2 {
                            quietly gen byte `excluded' = ///
                                `fold1' == `f1' | `fold2' == `f2'
                            quietly bysort `a': egen long `within1' = total(`excluded')
                            quietly gen long `train' = `count' - `within1'
                            quietly summarize `train', meanonly
                            if missing(`mintrain') | r(min) < `mintrain' ///
                                local mintrain = r(min)
                            quietly drop `excluded' `within1' `train'
                        }
                    }
                }
                quietly drop `count'
            }
        }
    restore
    return scalar N_complete = `n_complete'
    return scalar min_abs_count = `mincount'
    return scalar min_train_count = `mintrain'
end
