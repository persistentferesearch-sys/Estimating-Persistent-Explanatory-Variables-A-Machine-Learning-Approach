*! _sparse_python_paths 0.8.4-rc14  22jul2026
program _sparse_python_paths, sclass
    version 16.0
    syntax [, PYthon(string)]

    local pyexec = trim(`"`python'"')
    if `"`pyexec'"' == "" {
        capture local configured `"`c(python_exec)'"'
        if !_rc local pyexec = trim(`"`configured'"')
    }
    foreach code in 34 39 96 {
        local pyexec = subinstr(`"`pyexec'"', char(`code'), "", .)
    }
    if `"`pyexec'"' == "" {
        di as err `"configure Python or specify python("path/to/python.exe")"'
        exit 601
    }
    capture confirm file `"`pyexec'"'
    if _rc {
        di as err "Python executable not found: `pyexec'"
        exit 601
    }

    capture findfile sparse_prl_xpo_backend.py
    if !_rc local pyscript `"`r(fn)'"'
    else {
        capture findfile sparse_poreg.ado
        if _rc {
            di as err "cannot locate sparse_poreg.ado to resolve the Python backend"
            exit 601
        }
        local caller `"`r(fn)'"'
        mata: st_local("pyscript", pathresolve(pathgetparent(st_local("caller")), "../python/sparse_prl_xpo_backend.py"))
    }
    capture confirm file `"`pyscript'"'
    if _rc {
        di as err "Python backend script not found: sparse_prl_xpo_backend.py"
        exit 601
    }

    sreturn local python_exec `"`pyexec'"'
    sreturn local python_script `"`pyscript'"'
    sreturn local executable `"`pyexec'"'
    sreturn local script `"`pyscript'"'
end
