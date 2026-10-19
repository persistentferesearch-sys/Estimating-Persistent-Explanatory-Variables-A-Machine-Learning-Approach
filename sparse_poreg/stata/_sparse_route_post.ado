*! _sparse_route_post 0.8.4-rc14  22jul2026
program _sparse_route_post, eclass
    version 16.0
    syntax, REQUESTED(string) REASON(string) P(real) N(real) ///
        DENSEMAX(real) OFFICIALCAP(real) MAXVARHEADROOM(real)

    ereturn local sparse_backend_requested "`requested'"
    ereturn local sparse_route_reason "`reason'"
    ereturn scalar sparse_route_p = `p'
    ereturn scalar sparse_route_n = `n'
    ereturn scalar sparse_route_pn = cond(`n' > 0, `p'/`n', .)
    ereturn scalar sparse_route_densemax = `densemax'
    ereturn scalar sparse_route_official_cap = `officialcap'
    ereturn scalar sparse_route_maxvar_headroom = `maxvarheadroom'
    ereturn scalar sparse_route_official_headroom = `officialcap' - `p'
end
