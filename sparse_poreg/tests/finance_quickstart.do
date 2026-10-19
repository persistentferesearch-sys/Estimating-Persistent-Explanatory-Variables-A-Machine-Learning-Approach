version 19.0
clear all
set more off
args package_root
if trim(`"`package_root'"') != "" {
    local package_root = subinstr(`"`package_root'"', "\", "/", .)
    adopath ++ "`package_root'/stata"
}
capture which sparse_poreg
if _rc {
    display as error "sparse_poreg is not installed; see README.md for installation"
    exit 111
}

* Artificial firm-year data; no private or real research data are used.
set seed 1103
set obs 1000
generate long observation = _n
generate int firm = ceil(_n/10)
generate int year = 2000 + mod(_n - 1, 10)
generate double size = rnormal()
generate double cash = rnormal()
generate double leverage = rnormal()
generate double rd = .40*size + .20*cash + rnormal()
generate double patents = .20*rd + .30*size - .10*leverage + ///
    .15*sin(firm) + .05*(year - 2000) + rnormal()

local model patents rd, absorb(firm year) controls(size cash leverage) ///
    selection(bic, gridminok) vce(cluster firm)

* On this small example, auto uses official Stata.
quietly sparse_poreg `model' backend(auto)
scalar b_auto = _b[rd]
scalar se_auto = _se[rd]

* The same model using the package's fast large-fixed-effect route.
sort observation, stable
quietly sparse_poreg `model' backend(mata) strict(stata)
scalar b_mata = _b[rd]
scalar se_mata = _se[rd]

display as text _newline "Simple R&D-patent example"
display as text "---------------------------------------------"
display as text "Small-data official route: coefficient = " ///
    as result %8.4f b_auto as text ", SE = " as result %8.4f se_auto
display as text "Fast Mata route:           coefficient = " ///
    as result %8.4f b_mata as text ", SE = " as result %8.4f se_mata

assert abs(b_auto - b_mata) < 1e-8
assert abs(se_auto - se_mata) < 1e-8
display as result _newline "PASS: the official and fast Mata results agree."
