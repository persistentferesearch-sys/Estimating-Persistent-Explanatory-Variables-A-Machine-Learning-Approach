{smcl}
{* *! version 0.8.4-rc14  22jul2026}{...}
{vieweralsosee "sparse_poreg" "help sparse_poreg"}{...}
{vieweralsosee "xporegress" "help xporegress"}{...}
{vieweralsosee "Stata lasso inference examples" "browse https://www.stata.com/manuals/lassoinferenceexamples.pdf"}{...}
{viewerjumpto "Quick start" "sparse_xporeg##quick"}{...}
{viewerjumpto "Syntax" "sparse_xporeg##syntax"}{...}
{viewerjumpto "Backends" "sparse_xporeg##backends"}{...}
{viewerjumpto "Folds" "sparse_xporeg##folds"}{...}
{viewerjumpto "Options" "sparse_xporeg##options"}{...}
{viewerjumpto "Examples" "sparse_xporeg##examples"}{...}
{viewerjumpto "Stored results" "sparse_xporeg##results"}{...}

{title:Title}

{p2colset 5 23 25 2}{...}
{p2col:{bf:sparse_xporeg}}Cross-fit partialing-out lasso with official,
sparse Mata, and experimental Python backends{p_end}
{p2colreset}{...}


{marker quick}{...}
{title:Quick start}

{pstd}
Use {cmd:backend(auto)} for small and medium designs.  It delegates the entire
fit to official {helpb xporegress}.

{phang2}{cmd:. use https://www.stata-press.com/data/r19/lassoex, clear}{p_end}
{phang2}{cmd:. sparse_xporeg y d1 d2, absorb(f1-f30) controls(x1-x100) ///}{p_end}
{phang2}{cmd:>     selection(plugin) vce(robust) xfolds(10) ///}{p_end}
{phang2}{cmd:>     rseed(12345) backend(auto)}{p_end}

{pstd}
The package never switches silently from official XPO to experimental Mata
XPO.  You must choose it explicitly with {cmd:backend(mata) frontierok}.

{phang2}{cmd:. sparse_xporeg y d1 d2, absorb(firm year) controls(x1-x20) ///}{p_end}
{phang2}{cmd:>     selection(bic, gridminok) vce(cluster firm) ///}{p_end}
{phang2}{cmd:>     xfolds(5) rseed(1103) backend(mata) frontierok ///}{p_end}
{phang2}{cmd:>     strict(production)}{p_end}

{phang2}{cmd:. estat selection}{p_end}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:sparse_xporeg} {depvar} {it:varsofinterest} {ifin}{cmd:,}
[{opth abs:orb(varlist)}
{cmd:controls(}[{cmd:(}{it:alwaysvars}{cmd:)}] {it:othervars}{cmd:)}
{cmd:selection(}{cmd:plugin}|{cmd:bic}[{cmd:, gridminok}]|
{cmd:cv}|{cmd:adaptive}|{cmd:lassocv}{cmd:)}
{cmd:vce(robust}|{cmd:cluster} {it:clustvar}
[{it:clustvar2}]{cmd:)}
{opt xfolds(#)} {opt rseed(#)} {opth foldcluster(varlist)}
{cmd:backend(auto}|{cmd:official}|{cmd:mata}|{cmd:python)}
{opt densemax(#)} {opt strict(stata|production)} {opt frontierok}
{opt accuracy(high|fast)}
{opt maxiter(#)} {opt coeftol(#)}
{opt python(path)} {opt learner(lassocv)}
{opt lassocvfolds(#)} {opt lassocvalphas(#)}
{opt lassocveps(#)} {opt lassocvmaxiter(#)}
{opt lassocvtol(#)}]

{pstd}
At least one of {opt absorb()} or {opt controls()} is required.  Defaults are
{cmd:selection(plugin)}, {cmd:vce(robust)}, {cmd:xfolds(10)},
{cmd:rseed(1103)}, and {cmd:backend(auto)}.


{title:Description}

{pstd}
{cmd:sparse_xporeg} estimates a cross-fit partialing-out linear model.
Nuisance models are learned outside each held-out fold, orthogonal scores are
formed in the holdout sample, and the scores are combined using a DML2-style
estimate.

{pstd}
Every observed level of each variable in {opt absorb()} is a penalized
candidate control.  Thus {cmd:absorb(f1-f30)} matches
{cmd:ibn.(f1-f30)}, not the base-omitting {cmd:i.(f1-f30)} in the Stata
manual.  The examples below replay the manual's data and variables with
full-level indicators so official and Mata use the same design.


{marker backends}{...}
{title:Backends and expected results}

{synoptset 20 tabbed}{...}
{synopthdr:backend}
{synoptline}
{synopt:{cmd:auto}}Recommended for manageable designs.  Runs official Stata.
If official expansion is infeasible, auto stops and prints how to opt in to
the package's experimental Mata XPO.{p_end}
{synopt:{cmd:official}}Official {cmd:xporegress}: plugin/BIC/CV/adaptive,
robust or one-way cluster, observation folds, resampling, DML1/DML2, and
official postestimation.{p_end}
{synopt:{cmd:mata}}Experimental sparse plugin/BIC XPO.  Supports robust,
one-way cluster, and one-/two-way cluster-disjoint folds.  Always requires
{cmd:frontierok}, which confirms that the user intentionally chose this
experimental route.{p_end}
{synopt:{cmd:python}}Experimental DoubleML/LassoCV XPO.  Supports robust
observation folds and one-/two-way cluster-disjoint folds.  It is not
plugin/BIC replication.{p_end}
{synoptline}

{pstd}
For ordinary observation folds, the same full-level design, selection rule,
VCE, seed, fold count, sample, and order should give official and Mata results
that agree to numerical precision in supported parity cases.
{cmd:strict(stata)} is the parity-audit mode.

{pstd}
Python selects by LassoCV, so its results need not equal plugin or BIC.
It adds a sensitivity estimator with full cross-treatment covariance and
cluster-disjoint one-/two-way score folds.

{pstd}
Official {cmd:xporegress} has no two-way cluster or multiway cluster-disjoint
route.  Those examples therefore show Mata and Python only.


{marker folds}{...}
{title:Observation and cluster-disjoint folds}

{phang}
Without {cmd:foldcluster()}, official and Mata XPO split observations.
With {cmd:vce(cluster g)}, a cluster may occur in both training and holdout
samples.  This is the official-parity contract.

{phang}
{cmd:foldcluster(g)} assigns whole clusters to K outer folds.  It must match
{cmd:vce(cluster g)}.

{phang}
{cmd:foldcluster(g1 g2)} partitions each cluster dimension K ways and creates
K-squared score cells.  For each cell, training excludes observations sharing
either held-out partition.  With {cmd:xfolds(3)}, there are nine score cells.

{pstd}
Two-way cluster-disjoint XPO can be expensive: with q variables of interest,
it fits K-squared times q+1 nuisance models, and BIC searches a lambda path
inside each fit.  Choose K deliberately and benchmark locally.

{pstd}
Native {cmd:strict(production)} audits every training sample and returns
{cmd:r(498)} if an absorbed level has zero or one training observation.
This is a safety refusal, not a convergence failure.  Use fewer folds, repair
rare levels, or use official XPO when its observation-fold contract is
scientifically appropriate.


{title:Which backend should I use?}

{pstd}
P is the number of candidate controls after Stata expands firm, person, year,
and other categorical variables into dummies.  After a successful fit, see it
with

{phang2}{cmd:. display e(sparse_route_p)}{p_end}

{pstd}{bf:Simple XPO rule}{p_end}

{phang}{bf:Routine robust or one-way clustered XPO:} Use
{cmd:backend(auto)}.  It uses official {cmd:xporegress} whenever the expanded
design fits.{p_end}

{phang}{bf:Auto says the design is too large:} Auto stops and prints P.  It
does not switch silently to Mata XPO.  Simplify the design, or deliberately
choose the experimental Mata route described next.{p_end}

{phang}{bf:Experimental Mata XPO:} Use
{cmd:backend(mata) frontierok strict(production)}.  Every training fold must
pass the support and numerical safety checks.{p_end}

{pstd}
The three options mean:{p_end}

{phang}{cmd:backend(mata)} uses the package's experimental sparse XPO.{p_end}
{phang}{cmd:frontierok} confirms that you intentionally chose it; it does not
turn off a check.{p_end}
{phang}{cmd:strict(production)} stops instead of reporting estimates when a
safety check fails.{p_end}
{pstd}
When {cmd:absorb()} is used, Mata XPO chooses its fixed-effect computation
method automatically.  There is nothing else for routine users to choose.

{pstd}
Use this Mata command only when all three statements are true:{p_end}

{phang}1. Official XPO is too large for Stata, or you need folds that keep all
observations from the same firm or cluster together.{p_end}
{phang}2. In every training sample, each fixed-effect group (for example, each
firm, person, or year) has at least two observations.{p_end}
{phang}3. {cmd:strict(production)} finishes without an error.  If a support or
numerical check fails, the command stops instead of reporting estimates.{p_end}

{pstd}
Passing the checks means the calculation succeeded safely for that sample; it
does not guarantee that Mata XPO will match official {cmd:xporegress} in every
possible data set.

{phang}{bf:Folds that keep whole firms or clusters together:} Add matching
{cmd:vce(cluster ...)} and {cmd:foldcluster(...)} options to the Mata command.
One cluster variable gives one-way folds; two variables give two-way folds.
Official {cmd:xporegress} has no equivalent two-way estimator.{p_end}

{pstd}
There is no release-approved P-only speed cutoff for Mata XPO.  Cross-fitting
repeats the control-selection regressions in every fold; two-way folds repeat
them across every pair of folds.  Time the exact model locally, and use Mata
XPO only under the three conditions above.


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

{dlgtab:Folds and inference}

{phang}
{opt xfolds(#)} sets K; default 10 and minimum 2.  In two-way
cluster-disjoint XPO it creates K-squared score cells.

{phang}
{opt rseed(#)} sets the outer split seed; default 1103.

{phang}
{opth foldcluster(varlist)} requests one- or two-way cluster-disjoint outer
folds.  Its variables must exactly match {cmd:vce(cluster ...)}.

{phang}
Official {cmd:resample}[{cmd:(}{it:#}{cmd:)}] and
{cmd:technique(dml1|dml2)} are forwarded only to {cmd:xporegress}.
Mata is DML2-style and does not implement official resampling.

{dlgtab:Routing and safety}

{phang}
{opt densemax(#)} sets the official-expansion ceiling; default 20,000.
It cannot override {cmd:set maxvar}.

{phang}
{cmd:strict(stata)} is the default parity mode.  With observation folds,
native XPO runs an official design screen inside each training fold.
This is dense and is not the scalable large-P route.

{phang}
{cmd:strict(production)} skips official screens and requires safe training
support plus optimizer, selection, rank, and finite-result checks.

{phang}
Every Mata XPO fit requires {cmd:backend(mata) frontierok}.
Here {cmd:frontierok} confirms an intentional experimental choice; it does not
relax any safety check.  It is invalid with auto, official, or Python backends.

{cmd:accuracy()}, {cmd:maxiter()}, and {cmd:coeftol()} control the native
solver.

{dlgtab:Python}

{pstd}
Python requires {cmd:selection(lassocv)} and {cmd:learner(lassocv)}.
Supply the executable with {cmd:python()} or configure Stata's Python.
Clustered Python XPO requires matching {cmd:vce(cluster ...)} and
{cmd:foldcluster(...)}.  LassoCV inner tuning folds remain row-wise.


{marker examples}{...}
{title:Examples}

{pstd}{bf:1. Replay the Stata-manual data and variable lists}{p_end}

{phang2}{cmd:. use https://www.stata-press.com/data/r19/lassoex, clear}{p_end}
{phang2}{cmd:. local x x1-x100}{p_end}
{phang2}{cmd:. local f f1-f30}{p_end}
{phang2}{cmd:. generate int g1 = mod(_n - 1, 12) + 1}{p_end}
{phang2}{cmd:. generate int g2 = mod(floor((_n - 1)/12), 9) + 1}{p_end}

{pstd}
The generated clusters make the syntax reproducible; replace them with real
cluster IDs in empirical work.  All comparisons below use
{cmd:xfolds(3) rseed(1103)} so official and Mata replay the same folds.

{pstd}{bf:2. Official backend: plugin/BIC by robust/one-way cluster}{p_end}

{phang2}{cmd:. sparse_xporeg y d1 d2, absorb(`f') controls(`x') selection(plugin) vce(robust) xfolds(3) rseed(1103) backend(official)}{p_end}
{phang2}{cmd:. sparse_xporeg y d1 d2, absorb(`f') controls(`x') selection(plugin) vce(cluster g1) xfolds(3) rseed(1103) backend(official)}{p_end}
{phang2}{cmd:. sparse_xporeg y d1 d2, absorb(`f') controls(`x') selection(bic, gridminok) vce(robust) xfolds(3) rseed(1103) backend(official)}{p_end}
{phang2}{cmd:. sparse_xporeg y d1 d2, absorb(`f') controls(`x') selection(bic, gridminok) vce(cluster g1) xfolds(3) rseed(1103) backend(official)}{p_end}
{phang2}{cmd:. lassoinfo, each}{p_end}

{pstd}{bf:3. Mata observation-fold backend: the same four estimators}{p_end}

{phang2}{cmd:. sparse_xporeg y d1 d2, absorb(`f') controls(`x') selection(plugin) vce(robust) xfolds(3) rseed(1103) backend(mata) frontierok strict(stata)}{p_end}
{phang2}{cmd:. sparse_xporeg y d1 d2, absorb(`f') controls(`x') selection(plugin) vce(cluster g1) xfolds(3) rseed(1103) backend(mata) frontierok strict(stata)}{p_end}
{phang2}{cmd:. sparse_xporeg y d1 d2, absorb(`f') controls(`x') selection(bic, gridminok) vce(robust) xfolds(3) rseed(1103) backend(mata) frontierok strict(stata)}{p_end}
{phang2}{cmd:. sparse_xporeg y d1 d2, absorb(`f') controls(`x') selection(bic, gridminok) vce(cluster g1) xfolds(3) rseed(1103) backend(mata) frontierok strict(stata)}{p_end}
{phang2}{cmd:. estat selection}{p_end}

{pstd}
These clustered commands still use observation folds; clustering affects
selection and covariance.  That is the official-parity contract.

{pstd}{bf:4. Verify official-versus-Mata numerical agreement}{p_end}

{phang2}{cmd:. sparse_xporeg y d1 d2, absorb(`f') controls(`x') selection(bic, gridminok) vce(cluster g1) xfolds(3) rseed(1103) backend(official)}{p_end}
{phang2}{cmd:. matrix b_off = e(b)}{p_end}
{phang2}{cmd:. matrix V_off = e(V)}{p_end}
{phang2}{cmd:. sparse_xporeg y d1 d2, absorb(`f') controls(`x') selection(bic, gridminok) vce(cluster g1) xfolds(3) rseed(1103) backend(mata) frontierok strict(stata)}{p_end}
{phang2}{cmd:. matrix b_mata = e(b)}{p_end}
{phang2}{cmd:. matrix V_mata = e(V)}{p_end}
{phang2}{cmd:. mata: max(abs(st_matrix("b_off") :- st_matrix("b_mata")))}{p_end}
{phang2}{cmd:. mata: max(abs(st_matrix("V_off") :- st_matrix("V_mata")))}{p_end}

{pstd}
The last two values should be near machine precision in a supported parity
case.  If not, first check sample, sort order, factor expansion, seed, fold
count, selection rule, VCE, and safety diagnostics.

{pstd}{bf:5. Mata cluster-disjoint XPO on lassoex}{p_end}

{phang2}{cmd:. sparse_xporeg y d1 d2, absorb(`f') controls(`x') selection(plugin) vce(cluster g1) foldcluster(g1) xfolds(3) rseed(1103) backend(mata) frontierok strict(production)}{p_end}
{phang2}{cmd:. sparse_xporeg y d1 d2, absorb(`f') controls(`x') selection(bic, gridminok) vce(cluster g1) foldcluster(g1) xfolds(3) rseed(1103) backend(mata) frontierok strict(production)}{p_end}
{phang2}{cmd:. sparse_xporeg y d1 d2, absorb(`f') controls(`x') selection(bic, gridminok) vce(cluster g1 g2) foldcluster(g1 g2) xfolds(3) rseed(1103) backend(mata) frontierok strict(production)}{p_end}
{phang2}{cmd:. display e(n_score_folds) e(min_train_support)}{p_end}

{pstd}
The two-way BIC command creates nine score folds.  The row-number clusters are
only a syntax demonstration.  On this artificial design, two-way plugin
correctly fails its production loading certificate with {cmd:r(430)}.
Official {cmd:xporegress} has no corresponding two-way estimator.

{pstd}{bf:6. Python: robust and one-/two-way cluster-disjoint LassoCV}{p_end}

{phang2}{cmd:. sparse_xporeg y d1 d2, absorb(`f') controls(`x') selection(lassocv) learner(lassocv) vce(robust) xfolds(3) rseed(1103) backend(python) python("C:/path/to/python.exe")}{p_end}
{phang2}{cmd:. sparse_xporeg y d1 d2, absorb(`f') controls(`x') selection(lassocv) learner(lassocv) vce(cluster g1) foldcluster(g1) xfolds(3) rseed(1103) backend(python) python("C:/path/to/python.exe")}{p_end}
{phang2}{cmd:. sparse_xporeg y d1 d2, absorb(`f') controls(`x') selection(lassocv) learner(lassocv) vce(cluster g1 g2) foldcluster(g1 g2) xfolds(3) rseed(1103) backend(python) python("C:/path/to/python.exe")}{p_end}
{phang2}{cmd:. matrix list e(V)}{p_end}

{pstd}
Python results are not expected to equal plugin/BIC.  They provide a LassoCV
sensitivity estimate and the full covariance between {cmd:d1} and {cmd:d2}.

{pstd}{bf:7. Fold-safe Mata two-way plugin and BIC}{p_end}

{phang2}{cmd:. findfile validate_public_examples.do}{p_end}
{phang2}{cmd:. local examples "`r(fn)'"}{p_end}
{phang2}{cmd:. do "`examples'"}{p_end}
{phang2}{cmd:. sparse_xporeg y d1 d2, absorb(a b) controls(x1 x2) selection(plugin) vce(cluster g1 g2) foldcluster(g1 g2) xfolds(3) rseed(1103) backend(mata) frontierok strict(production)}{p_end}
{phang2}{cmd:. sparse_xporeg y d1 d2, absorb(a b) controls(x1 x2) selection(bic, gridminok) vce(cluster g1 g2) foldcluster(g1 g2) xfolds(3) rseed(1103) backend(mata) frontierok strict(production)}{p_end}

{pstd}
The installed fixture gives every absorbed level safe support in every
two-way training sample, so both production commands return certified
results.


{title:Postestimation and safety}

{pstd}
After official XPO, use {cmd:lassoinfo, each}.  After Mata XPO, use
{cmd:estat selection} or {cmd:estat selection, support}.  Mata also reports
minimum full-sample and training-fold absorbed-level support.

{pstd}
Mata returns {cmd:r(430)} for uncertified KKT/plugin paths and
{cmd:r(498)} for unsafe training support.  Under explicit
{cmd:strict(stata)}, clustered plugin frontier fits may warn and post
{cmd:e(kkt_pass)=0} or {cmd:e(plugin_pass)=0}; production refuses them.


{marker results}{...}
{title:Stored results}

{synoptset 32 tabbed}{...}
{synopt:{cmd:e(b)}}coefficients for variables of interest{p_end}
{synopt:{cmd:e(V)}}full covariance matrix{p_end}
{synopt:{cmd:e(sparse_backend)}}backend actually used{p_end}
{synopt:{cmd:e(sparse_route_reason)}}routing decision{p_end}
{synopt:{cmd:e(sparse_route_p)}}expanded candidate-control count P{p_end}
{synopt:{cmd:e(sparse_route_n)}}estimation-sample size N{p_end}
{synopt:{cmd:e(sparse_route_pn)}}P/N{p_end}
{synopt:{cmd:e(k_controls_sel)}}selected-control union size{p_end}
{synopt:{cmd:e(selection_counts)}}fold/target selection audit{p_end}
{synopt:{cmd:e(n_score_folds)}}number of outer score folds{p_end}
{synopt:{cmd:e(min_train_support)}}minimum absorbed-level training support{p_end}
{synopt:{cmd:e(training_support_safe)}}training-support certificate{p_end}
{synopt:{cmd:e(cluster_disjoint)}}cluster-disjoint-fold indicator{p_end}
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
{helpb sparse_poreg}, {helpb xporegress}, and the package
{cmd:docs/README.md}
