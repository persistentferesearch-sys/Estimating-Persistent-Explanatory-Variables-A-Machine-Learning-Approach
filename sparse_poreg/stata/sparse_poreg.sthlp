{smcl}
{* *! version 0.8.4-rc14  22jul2026}{...}
{vieweralsosee "sparse_xporeg" "help sparse_xporeg"}{...}
{vieweralsosee "poregress" "help poregress"}{...}
{vieweralsosee "Stata lasso inference examples" "browse https://www.stata.com/manuals/lassoinferenceexamples.pdf"}{...}
{viewerjumpto "Quick start" "sparse_poreg##quick"}{...}
{viewerjumpto "Syntax" "sparse_poreg##syntax"}{...}
{viewerjumpto "Backends" "sparse_poreg##backends"}{...}
{viewerjumpto "Options" "sparse_poreg##options"}{...}
{viewerjumpto "Examples" "sparse_poreg##examples"}{...}
{viewerjumpto "Stored results" "sparse_poreg##results"}{...}

{title:Title}

{p2colset 5 22 24 2}{...}
{p2col:{bf:sparse_poreg}}Partialing-out lasso with official, sparse Mata, and
experimental Python backends{p_end}
{p2colreset}{...}


{marker quick}{...}
{title:Quick start}

{pstd}
Use {cmd:backend(auto)} for ordinary work.  It runs official
{helpb poregress} when the expanded controls fit safely in Stata and otherwise
uses the sparse Mata engine when the requested options are supported.

{phang2}{cmd:. use https://www.stata-press.com/data/r19/lassoex, clear}{p_end}
{phang2}{cmd:. sparse_poreg y d1 d2, absorb(f1-f30) controls(x1-x100) ///}{p_end}
{phang2}{cmd:>     selection(plugin) vce(robust) backend(auto)}{p_end}

{pstd}
For a genuinely large dummy set, request Mata production checks explicitly.

{phang2}{cmd:. sparse_poreg y d1 d2, absorb(firm year) controls(x1-x20) ///}{p_end}
{phang2}{cmd:>     selection(bic, gridminok) vce(cluster firm) ///}{p_end}
{phang2}{cmd:>     backend(mata) strict(production)}{p_end}

{phang2}{cmd:. estat selection}{p_end}
{phang2}{cmd:. estat selection, support}{p_end}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:sparse_poreg} {depvar} {it:varsofinterest} {ifin}{cmd:,}
[{opth abs:orb(varlist)}
{cmd:controls(}[{cmd:(}{it:alwaysvars}{cmd:)}] {it:othervars}{cmd:)}
{cmd:selection(}{cmd:plugin}|{cmd:bic}[{cmd:, gridminok}]|
{cmd:cv}|{cmd:adaptive}|{cmd:lassocv}{cmd:)}
{cmd:vce(robust}|{cmd:cluster} {it:clustvar}
[{it:clustvar2}]{cmd:)}
{cmd:backend(auto}|{cmd:official}|{cmd:mata}|{cmd:python)}
{opt densemax(#)} {opt strict(stata|production)} {opt frontierok}
{opt accuracy(high|fast)}
{opt rseed(#)} {opt maxiter(#)} {opt coeftol(#)}
{opt python(path)} {opt learner(lassocv)}
{opt lassocvfolds(#)} {opt lassocvalphas(#)}
{opt lassocveps(#)} {opt lassocvmaxiter(#)}
{opt lassocvtol(#)}]

{pstd}
At least one of {opt absorb()} or {opt controls()} is required.  Defaults are
{cmd:selection(plugin)}, {cmd:vce(robust)}, and {cmd:backend(auto)}.


{title:Description}

{pstd}
{cmd:sparse_poreg} estimates a partialing-out linear model.  The variables of
interest are not penalized.  Candidate controls are selected for the outcome
and for each variable of interest, and their union is used for orthogonalized
estimation.

{pstd}
Every observed level of each variable in {opt absorb()} is a penalized
candidate control.  Thus {cmd:absorb(f1-f30)} matches
{cmd:ibn.(f1-f30)}, not the base-omitting {cmd:i.(f1-f30)} used in the
Stata manual.  The examples below use the manual's data and variables, but
use full-level indicators so the official and Mata designs are identical.


{marker backends}{...}
{title:Backends and expected results}

{synoptset 20 tabbed}{...}
{synopthdr:backend}
{synoptline}
{synopt:{cmd:auto}}Recommended default.  Official when feasible; otherwise
Mata for supported PO specifications.{p_end}
{synopt:{cmd:official}}Official {cmd:poregress}: plugin/BIC/CV/adaptive,
robust or one-way cluster, factor and always-in controls, and official
postestimation.{p_end}
{synopt:{cmd:mata}}Sparse plugin/BIC for large numeric controls and absorbed
IDs.  Robust and one-way cluster are production targets; two-way clustering
is explicit frontier functionality.{p_end}
{synopt:{cmd:python}}Experimental DoubleML/LassoCV PLR with robust, one-way,
or two-way covariance.  It is a sensitivity estimator, not plugin/BIC
replication.{p_end}
{synoptline}

{pstd}
With the same full-level design, selection rule, VCE, sample, and observation
order, supported official and Mata PO examples are expected to agree to
numerical precision.  {cmd:strict(stata)} is the parity-audit mode.
{cmd:strict(production)} is the scalable safety mode and may use a different
full-design contract when Stata's dense screen is infeasible.

{pstd}
Python selects by LassoCV, so its estimates need not equal plugin or BIC.
Python adds a direct sensitivity route and posts the full covariance across
multiple variables of interest.

{pstd}
Official {cmd:poregress} has no two-way cluster option.  Two-way examples
therefore show Mata and Python only.


{title:Which backend should I use?}

{pstd}
P is the number of candidate controls after Stata expands firm, person, year,
and other categorical variables into dummies.  It can be much larger than the
number of variables typed in the command.  After a successful fit, see it with

{phang2}{cmd:. display e(sparse_route_p)}{p_end}

{pstd}
The copy-ready fast Mata options are

{phang2}{cmd:backend(mata) strict(production)}{p_end}

{pstd}
In plain English, {cmd:backend(mata)} uses the package's sparse engine, and
{cmd:strict(production)} stops if a safety check fails.  When {cmd:absorb()} is
used, Mata chooses its fixed-effect computation method automatically.  There is
nothing else for routine users to choose.

{pstd}{bf:Simple PO rule of thumb}{p_end}

{phang}{bf:Plugin:} If P is below 400, use {cmd:backend(auto)}.  If P is 400
or more, try the fast Mata options shown above.  At P=412, Mata was 1.76 to
2.09 times faster in this release's tests.{p_end}

{phang}{bf:BIC:} If P is below 800, use {cmd:backend(auto)}.  If P is 800 or
more, try the fast Mata options shown above.  Mata was about 1.57 times faster
at P=812, 2.39 to 2.74 times faster at P=1,612, and 5.07 to 5.12 times faster
at P=3,212 in the tested designs.{p_end}

{phang}{bf:Two-way clustering:} Official {cmd:poregress} does not provide it.
Use {cmd:backend(mata) frontierok strict(production)}, or use Python LassoCV as
a sensitivity analysis.  Here {cmd:frontierok} confirms that you intentionally
chose a package extension not available in official Stata; it does not turn
off the production checks.{p_end}

{phang}{bf:Official expansion does not fit:} Use
the fast Mata options shown above regardless of P.  The official limit is the
smaller of {cmd:densemax()} (20,000 by default) and the available
{cmd:set maxvar} headroom.{p_end}

{pstd}
These are practical starting points, not universal cutoffs.  Runtime also
depends on sample size, the number of variables of interest, and the BIC path.
Near a cutoff, time both backends on the exact model.  Below a cutoff,
{cmd:backend(auto)} favors official-Stata reproducibility; it does not promise
that official Stata is faster on every design.


{marker options}{...}
{title:Options}

{dlgtab:Model and selection}

{phang}
{opth absorb(varlist)} specifies numeric categorical IDs.  Their full observed
indicator sets are penalized candidates.

{phang}
{cmd:controls(}[{cmd:(}{it:alwaysvars}{cmd:)}] {it:othervars}{cmd:)}
specifies controls.  Parenthesized controls are unpenalized on the official
path.  Mata accepts a plain numeric varlist only.

{phang}
{cmd:selection(plugin)} is the default.  {cmd:selection(bic, gridminok)}
selects an identified BIC minimum and enables the documented grid-minimum
rule.  CV and adaptive are official-only; LassoCV is Python-only.

{dlgtab:Inference}

{phang}
{cmd:vce(robust)} is the default.  Official and Mata support one-way
{cmd:vce(cluster g)}.  Mata two-way {cmd:vce(cluster g1 g2)} requires
{cmd:frontierok}.  Python supports robust, one-way, and two-way covariance.

{dlgtab:Routing and safety}

{phang}
{opt densemax(#)} sets the user's official-expansion ceiling; default 20,000.
It cannot override {cmd:set maxvar}.

{phang}
{cmd:strict(stata)} is the default compatibility mode.  Native BIC uses one
official design screen and preserves input order.  Sort explicitly before
official/Mata parity comparisons.

{phang}
{cmd:strict(production)} skips dense screening and requires full KKT, plugin,
rank, and finite-result safety checks.  Uncertified native fits return
{cmd:r(430)} rather than estimates.

{phang}
{cmd:accuracy(high|fast)}, {opt maxiter(#)}, and {opt coeftol(#)} control the
Mata solver.  Defaults are high, 4,000, and .001.

{dlgtab:Official and Python options}

{pstd}
The official route forwards supported {cmd:poregress} options, including
{cmd:sqrtlasso}, {cmd:semi}, {cmd:missingok}, lasso overrides, and reporting
or optimization options.  Official-only options are never silently ignored.

{pstd}
Python requires {cmd:selection(lassocv)} and {cmd:learner(lassocv)}.
Supply the executable with {cmd:python()} or configure Stata's Python.
The {cmd:lassocv*()} options control inner LassoCV tuning.


{marker examples}{...}
{title:Examples}

{pstd}{bf:Start here: run the simple finance example}{p_end}

{phang2}{cmd:. findfile finance_quickstart.do}{p_end}
{phang2}{cmd:. do "`r(fn)'"}{p_end}

{pstd}
It creates artificial firm-year data and compares the small-data official
route with the fast Mata route.  A green {cmd:PASS} means both results agree.

{pstd}{bf:1. Replay the Stata-manual data and variable lists}{p_end}

{phang2}{cmd:. use https://www.stata-press.com/data/r19/lassoex, clear}{p_end}
{phang2}{cmd:. local x x1-x100}{p_end}
{phang2}{cmd:. local f f1-f30}{p_end}
{phang2}{cmd:. generate int g1 = mod(_n - 1, 12) + 1}{p_end}
{phang2}{cmd:. generate int g2 = mod(floor((_n - 1)/12), 9) + 1}{p_end}

{pstd}
The artificial {cmd:g1} and {cmd:g2} variables make the one- and two-way
syntax reproducible.  Replace them with scientifically meaningful clusters
in real work.

{pstd}{bf:2. Official backend: complete plugin/BIC by robust/one-way grid}{p_end}

{phang2}{cmd:. sparse_poreg y d1 d2, absorb(`f') controls(`x') selection(plugin) vce(robust) backend(official)}{p_end}
{phang2}{cmd:. sparse_poreg y d1 d2, absorb(`f') controls(`x') selection(plugin) vce(cluster g1) backend(official)}{p_end}
{phang2}{cmd:. sparse_poreg y d1 d2, absorb(`f') controls(`x') selection(bic, gridminok) vce(robust) backend(official)}{p_end}
{phang2}{cmd:. sparse_poreg y d1 d2, absorb(`f') controls(`x') selection(bic, gridminok) vce(cluster g1) backend(official)}{p_end}
{phang2}{cmd:. lassoinfo, each}{p_end}

{pstd}{bf:3. Mata backend: the same four estimators}{p_end}

{phang2}{cmd:. sparse_poreg y d1 d2, absorb(`f') controls(`x') selection(plugin) vce(robust) backend(mata) strict(stata)}{p_end}
{phang2}{cmd:. sparse_poreg y d1 d2, absorb(`f') controls(`x') selection(plugin) vce(cluster g1) backend(mata) strict(stata)}{p_end}
{phang2}{cmd:. sparse_poreg y d1 d2, absorb(`f') controls(`x') selection(bic, gridminok) vce(robust) backend(mata) strict(stata)}{p_end}
{phang2}{cmd:. sparse_poreg y d1 d2, absorb(`f') controls(`x') selection(bic, gridminok) vce(cluster g1) backend(mata) strict(stata)}{p_end}
{phang2}{cmd:. estat selection}{p_end}

{pstd}
Use {cmd:strict(stata)} above to audit parity on this manageable design.
For a large sparse design, change it to {cmd:strict(production)}.

{pstd}{bf:4. Verify official-versus-Mata numerical agreement}{p_end}

{phang2}{cmd:. sparse_poreg y d1 d2, absorb(`f') controls(`x') selection(bic, gridminok) vce(cluster g1) backend(official)}{p_end}
{phang2}{cmd:. matrix b_off = e(b)}{p_end}
{phang2}{cmd:. matrix V_off = e(V)}{p_end}
{phang2}{cmd:. sparse_poreg y d1 d2, absorb(`f') controls(`x') selection(bic, gridminok) vce(cluster g1) backend(mata) strict(stata)}{p_end}
{phang2}{cmd:. matrix b_mata = e(b)}{p_end}
{phang2}{cmd:. matrix V_mata = e(V)}{p_end}
{phang2}{cmd:. mata: max(abs(st_matrix("b_off") :- st_matrix("b_mata")))}{p_end}
{phang2}{cmd:. mata: max(abs(st_matrix("V_off") :- st_matrix("V_mata")))}{p_end}

{pstd}
The last two values should be near machine precision.  The comparison is
meaningful because both commands use the same full-level design, selection
rule, VCE, sample, and order.

{pstd}{bf:5. Mata two-way clustering: plugin and BIC}{p_end}

{phang2}{cmd:. sparse_poreg y d1 d2, absorb(`f') controls(`x') selection(plugin) vce(cluster g1 g2) backend(mata) frontierok strict(production)}{p_end}
{phang2}{cmd:. sparse_poreg y d1 d2, absorb(`f') controls(`x') selection(bic, gridminok) vce(cluster g1 g2) backend(mata) frontierok strict(production)}{p_end}

{pstd}
Official {cmd:poregress} has no corresponding two-way command.

{pstd}{bf:6. Python: robust, one-way, and two-way LassoCV sensitivity fits}{p_end}

{phang2}{cmd:. sparse_poreg y d1 d2, absorb(`f') controls(`x') selection(lassocv) learner(lassocv) vce(robust) backend(python) python("C:/path/to/python.exe")}{p_end}
{phang2}{cmd:. sparse_poreg y d1 d2, absorb(`f') controls(`x') selection(lassocv) learner(lassocv) vce(cluster g1) backend(python) python("C:/path/to/python.exe")}{p_end}
{phang2}{cmd:. sparse_poreg y d1 d2, absorb(`f') controls(`x') selection(lassocv) learner(lassocv) vce(cluster g1 g2) backend(python) python("C:/path/to/python.exe")}{p_end}
{phang2}{cmd:. matrix list e(V)}{p_end}

{pstd}
The Python coefficients are not expected to equal plugin/BIC because the
selection method is LassoCV.  Python does post a complete covariance matrix
for {cmd:d1} and {cmd:d2}, including their off-diagonal covariance.

{pstd}{bf:7. Run the complete installed validation fixture}{p_end}

{phang2}{cmd:. findfile validate_public_examples.do}{p_end}
{phang2}{cmd:. local examples "`r(fn)'"}{p_end}
{phang2}{cmd:. do "`examples'"}{p_end}


{title:Postestimation and safety}

{pstd}
After an official fit, use {cmd:lassoinfo, each}.  After a Mata fit, use
{cmd:estat selection} or {cmd:estat selection, support}.  A selected union is
the set chosen by at least one nuisance model, not necessarily every model.

{pstd}
Mata returns {cmd:r(430)} when the final lasso is not KKT/plugin certified.
Clustered plugin selection also rejects singleton absorbed levels.  Combine
or remove singleton categories, use BIC, or use the official backend.

{pstd}
If a variable of interest is collinear after partialing out, Mata follows
official {cmd:poregress}: it reports that coefficient as omitted, stores zero
in the corresponding row and column of {cmd:e(b)} and {cmd:e(V)}, and retains
the estimable coefficients.  The omitted names are saved in
{cmd:e(varsofinterest_omitted)}.


{marker results}{...}
{title:Stored results}

{synoptset 31 tabbed}{...}
{synopt:{cmd:e(b)}}coefficients for variables of interest{p_end}
{synopt:{cmd:e(V)}}full covariance matrix{p_end}
{synopt:{cmd:e(sparse_backend)}}backend actually used{p_end}
{synopt:{cmd:e(sparse_route_reason)}}routing decision{p_end}
{synopt:{cmd:e(sparse_route_p)}}expanded candidate-control count P{p_end}
{synopt:{cmd:e(sparse_route_n)}}estimation-sample size N{p_end}
{synopt:{cmd:e(sparse_route_pn)}}P/N{p_end}
{synopt:{cmd:e(k_controls_sel)}}selected-control union size{p_end}
{synopt:{cmd:e(controls_sel)}}selected union, when it fits a macro{p_end}
{synopt:{cmd:e(k_varsofinterest_omitted)}}number of collinear variables of interest omitted{p_end}
{synopt:{cmd:e(varsofinterest_omitted)}}names of omitted variables of interest{p_end}
{synopt:{cmd:e(selection_counts)}}target-specific selection audit{p_end}
{synopt:{cmd:e(kkt_pass)}}native KKT certificate{p_end}
{synopt:{cmd:e(plugin_pass)}}native plugin certificate{p_end}
{synopt:{cmd:e(post_rank_valid)}}final rank check{p_end}
{synopt:{cmd:e(finite)}}finite-result check{p_end}


{title:License and trademark}

{pstd}
This independent third-party package is distributed under the MIT license.
It is not affiliated with or endorsed by StataCorp LLC.  Stata is a registered
trademark of StataCorp LLC.  The license covers this package only, not Stata
or other third-party software.


{title:Also see}

{psee}
{helpb sparse_xporeg}, {helpb poregress}, and the package
{cmd:docs/README.md}
