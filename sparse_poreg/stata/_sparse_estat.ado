*! _sparse_estat 0.8.4-rc14
program _sparse_estat, rclass
    version 16.0
    return clear
    if !inlist("`e(cmd)'", "sparse_poreg", "sparse_xporeg") exit 301

    gettoken sub rest : 0, parse(" ,")
    local n = length("`sub'")
    if "`sub'" == substr("selection", 1, max(3, `n')) {
        _sparse_estat_selection `rest'
        return add
        exit
    }
    estat_default `0'
end

program _sparse_estat_selection, rclass
    version 16.0
    syntax [, SUPPORT]

    capture confirm scalar e(k_lasso_fits)
    if _rc {
        di as err "per-fit selection details are unavailable"
        exit 321
    }

    local fits = e(k_lasso_fits)
    di as txt _newline "Nuisance-selection details"
    di as txt "Fit" _col(8) "Target" _col(42) "Fold 1" ///
        _col(52) "Fold 2" _col(62) "Selected"
    forvalues i = 1/`fits' {
        local dep `"`e(lasso_depvar_`i')'"'
        local f1 = e(lasso_fold_`i')
        local f2 = e(lasso_fold2_`i')
        local f1txt -
        local f2txt -
        if `f1' < . {
            local f1txt = string(`f1', "%9.0g")
            local f1txt = strtrim("`f1txt'")
        }
        if `f2' < . {
            local f2txt = string(`f2', "%9.0g")
            local f2txt = strtrim("`f2txt'")
        }
        di as res %4.0f `i' _col(8) "`dep'" _col(42) "`f1txt'" ///
            _col(52) "`f2txt'" _col(62) %9.0g e(lasso_sel_p_`i')
        if "`support'" != "" {
            if e(lasso_selected_complete_`i') {
                local selected `"`e(lasso_selected_`i')'"'
                if `"`selected'"' == "" local selected (none)
                di as txt "     support: " as res `"`selected'"'
            }
            else di as txt "     support: omitted; exceeds c(macrolen)"
        }
    }
    di as txt "Union selected across fits = " as res %9.0g e(k_controls_sel)

    local counts_complete = e(selection_counts_complete)
    if `counts_complete' {
        tempname counts
        matrix `counts' = e(selection_counts)
        return matrix counts = `counts'
    }
    else di as txt "Selection-count matrix omitted; fits exceed c(matsize)"
    return scalar selection_counts_complete = `counts_complete'
    return scalar k_lasso_fits = `fits'
    return scalar k_controls_sel = e(k_controls_sel)
end
