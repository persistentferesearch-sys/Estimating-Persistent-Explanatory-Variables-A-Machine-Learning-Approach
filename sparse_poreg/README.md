# sparse_poreg_pkg

`sparse_poreg` and `sparse_xporeg` use two deliberately different contracts:
official Stata for manageable designs and a narrower sparse Mata engine for
large absorbed-dummy sets.

Current release: `0.8.4-rc14`. It retains the rc13 routing, Python/DoubleML,
and certification baseline, and adds official-compatible omission of
collinear PO variables of interest in the Mata backend.  Omitted coefficients
and their covariance rows/columns are zero-filled, estimable coefficients are
retained, and the names are saved in `e(varsofinterest_omitted)`.  The
focused robust, one-way, two-way, public-example, and empirical ownership
validations pass.  See `CERTIFICATION_20260715.md` for the rc13 manifest and
documented frontiers.

This is an independent third-party package, not affiliated with or endorsed
by StataCorp LLC. Stata is a registered trademark of StataCorp LLC. References
to official commands describe compatibility testing; the MIT license covers
this package only, not Stata or other third-party software.

## Install

The repository root is a Stata net-install site:

```stata
net install sparse_poreg, from("C:/path/to/extracted/package") replace
help sparse_poreg
help sparse_xporeg
```

For a nontechnical coauthor check, open `START_HERE_FOR_COAUTHORS.txt`. After
selecting the extracted folder with Stata's **File > Change Working
Directory** menu, the complete check is only:

```stata
net install sparse_poreg, from("`c(pwd)'") replace
do "tests/finance_quickstart.do"
```

The example uses artificial firm-year data and prints a green `PASS`; it needs
no Python, PowerShell, or research data.

`stata.toc` and `sparse_poreg.pkg` install both commands, all runtime helpers,
the experimental Python bridge, and the copy-ready public fixture. The
isolated PLUS-tree smoke test is `tests/validate_net_install.do`.

## Which backend should you use?

`P` means the number of candidate controls after Stata expands firm, person,
year, and other categorical variables into dummies. It may be much larger than
the number of variables typed in the command. After any successful fit, type:

```stata
display "requested/chosen: `e(sparse_backend_requested)' / `e(sparse_backend)'"
display "P=" e(sparse_route_p) " N=" e(sparse_route_n) ///
    " P/N=" e(sparse_route_pn)
display "official cap=" e(sparse_route_official_cap) ///
    " headroom=" e(sparse_route_official_headroom)
```

### Practical PO speed rule

| model | P | command to use first | tested timing evidence |
|---|---:|---|---|
| PO, plugin | below 400 | `backend(auto)` | simplest official-first workflow |
| PO, plugin | 400 or more | `backend(mata) strict(production)` | Mata was 1.76-2.09x faster at P=412 |
| PO, BIC | below 800 | `backend(auto)` | simplest official-first workflow |
| PO, BIC | 800 or more | `backend(mata) strict(production)` | about 1.57x at P=812, 2.39-2.74x at P=1,612, and 5.07-5.12x at P=3,212 |
| PO, two-way clustering | any P | `backend(mata) frontierok strict(production)` | official `poregress` has no two-way estimator |
| PO, official expansion does not fit | any P | `backend(mata) strict(production)` | required sparse route |

These are conservative starting points from this release's synthetic tests,
not universal laws. Near a cutoff, time the exact model on both backends.
Below a cutoff, `backend(auto)` favors official-Stata reproducibility; it does
not promise that official Stata is always faster. `densemax()` (20,000 by
default) is a feasibility ceiling, not a speed cutoff. The timing comparisons
used Mata's automatic method. Routine users have nothing else to choose. The
diagnostic [timing matrix](../logs/benchmark_kernel_choice.csv) confirms exact
results; the alternative saved less than one second in its only faster case.

### Practical XPO rule

Use `backend(auto)` for routine robust or one-way clustered XPO whenever the
official expansion fits. Auto deliberately stops instead of silently changing
to experimental Mata XPO. Use explicit
`backend(mata) frontierok strict(production)` only when official
XPO is too large or you need folds that keep each firm or cluster together,
every fixed-effect group has at least two observations in every training
sample, and all production checks pass. `frontierok` simply confirms that you
intentionally selected experimental Mata XPO; it does not disable a safety
check. There is no release-approved P-only speed cutoff for Mata XPO. Two-way
folds that keep clusters together require explicit Mata or Python syntax;
official `xporegress` has no equivalent estimator.

### Why the cutoffs are approximate

N alone is not a useful switch. P controls official feasibility, P/N describes
statistical dimension, and N x P approximates one dense design pass. Actual
work also multiplies by the nuisance fits: with q treatments, PO fits q+1
targets, K-fold XPO fits K(q+1), and two-way K-by-K XPO fits K^2(q+1); BIC
also searches a lambda path.

`rc13_route_kernelauto_20260719` is current-source routing evidence, not a
universal cutoff. The detailed tested anchors are retained below for audit and
local benchmarking.

| profile | documented low-density anchor | documented high-density anchor | interpretation |
|---|---|---|---|
| PO plugin (robust and cluster) | `N=1000, P=112, P/N=.112, N*P=112000` | `N=1000, P=412, P/N=.412, N*P=412000` | observed sparse-performance anchor |
| PO BIC, robust | `N=8000, P=812, P/N=.1015, N*P=6496000` | `N=4000, P=1612, P/N=.403, N*P=6448000` | conservative cross-run sparse-performance anchor |
| PO BIC, cluster | `N=8000, P=812, P/N=.1015, N*P=6496000` | `N=1000, P=412, P/N=.412, N*P=412000` | observed sparse-performance anchor |
| XPO plugin (robust and cluster) | `N=1000, P=112, P/N=.112, N*P=112000` | `N=1000, P=412, P/N=.412, N*P=412000` | official default; sparse remains frontier |
| XPO BIC (robust and cluster) | `N=8000, P=812, P/N=.1015, N*P=6496000` | `N=1000, P=412, P/N=.412, N*P=412000` | official default; sparse remains frontier |

The current run directly places the high-density PO-BIC/robust anchor at
N=4,000/P=1,612 and the low-density PO-BIC anchors at N=8,000/P=812. All
anchors are advisory and conditional.

The official route preserves official results and postestimation while
forwarding factor and always-in controls, official selection methods,
`sqrtlasso`, `semi`, `missingok`, lasso overrides, reporting options, and XPO
`resample()` / `technique()` options. Official `poregress` and `xporegress`
support robust or at most one cluster variable; two-way clustering requires
explicit Mata frontier or Python syntax.

```stata
adopath ++ "C:/path/to/extracted/package/stata"
use https://www.stata-press.com/data/r19/breathe, clear
local cc no2_home age age0 sev_home green_home noise_school sev_school precip siblings_old siblings_young
local fc sex grade overweight lbweight breastfeed msmoke meducation feducation

* Default plugin/robust; manageable designs route to official Stata.
sparse_poreg react no2_class, absorb(`fc') controls(`cc') backend(auto)
sparse_xporeg react no2_class, absorb(`fc') controls(`cc') ///
    xfolds(10) rseed(12345) backend(auto)
```

## Large dummy sets: native Mata

The native large-design certification target uses numeric IDs in `absorb()`, a
plain numeric varlist in `controls()`, multiple treatments with full
covariance, and these PO profiles:

| command | selection | VCE | status |
|---|---|---|---|
| `sparse_poreg` | plugin | robust | native certification target |
| `sparse_poreg` | BIC | robust | native certification target |
| `sparse_poreg` | BIC | cluster | native certification target |
| `sparse_poreg` | plugin | cluster | native certification target; singleton absorbed levels rejected |
| `sparse_poreg` | plugin or BIC | two-way cluster | research frontier; requires `backend(mata) frontierok` |
| `sparse_xporeg` | any native profile | robust or cluster | research frontier; requires `backend(mata) frontierok` |

The installed public fixture creates a fold-safe data set and executes the
complete official/Mata and multiway grids:

```stata
findfile validate_public_examples.do
local examples "`r(fn)'"
do "`examples'"

* Plugin/BIC by robust and one-way cluster, PO and XPO.
foreach sel in plugin bic {
    local s "`sel'"
    if "`sel'" == "bic" local s "bic, gridminok"
    sparse_poreg y d1 d2, absorb(a b) controls(x1 x2) ///
        selection(`s') vce(robust) backend(mata) strict(production)
    sparse_poreg y d1 d2, absorb(a b) controls(x1 x2) ///
        selection(`s') vce(cluster g1) backend(mata) strict(production)
    sparse_xporeg y d1 d2, absorb(a b) controls(x1 x2) ///
        selection(`s') vce(robust) xfolds(3) rseed(1103) ///
        backend(mata) frontierok strict(production)
    sparse_xporeg y d1 d2, absorb(a b) controls(x1 x2) ///
        selection(`s') vce(cluster g1) xfolds(3) rseed(1103) ///
        backend(mata) frontierok strict(production)
}

* Plugin/BIC multiway: conventional PO and K^2 cluster-DML XPO.
foreach sel in plugin bic {
    local s "`sel'"
    if "`sel'" == "bic" local s "bic, gridminok"
    sparse_poreg y d1 d2, absorb(a b) controls(x1 x2) ///
        selection(`s') vce(cluster g1 g2) ///
        backend(mata) frontierok strict(production)
    sparse_xporeg y d1 d2, absorb(a b) controls(x1 x2) ///
        selection(`s') vce(cluster g1 g2) foldcluster(g1 g2) ///
        xfolds(3) rseed(1103) backend(mata) frontierok strict(production)
}
```

For native PO BIC, rc10 `strict(stata)` mirrors Stata's hidden design screen
once, then reuses that retained-column mask for every sparse nuisance fit. It
preserves input order and posts the kept/dropped counts, screening time, order
contract, and KKT scope. In a rank-deficient dummy design, Stata's screen can
depend on row order, so reproduce the official data preparation explicitly:

```stata
sort PERMCO fyear
sparse_poreg lnnpatent RDAT lnME RD_missing lnAge lnK2L TobinQ ROA ///
    Leverage CASHAT KZidx InstOwn oms_HHidx oms_HHidx_square, ///
    absorb(PERMCO fyear) vce(cluster PERMCO) ///
    selection(bic, gridminok) rseed(1103) ///
    backend(mata) strict(stata)
```

For ordinary observation-fold native XPO, rc11 `strict(stata)` runs the same
official screen separately inside each training fold and applies that fold's
mask only to that fold's nuisance fits. It warns, rather than aborts, when an
absorbed level is absent or singleton in a training sample because official
Stata handles that fold-local design. This mode requires each expanded fold
screen to fit under `set maxvar`.

`strict(production)` skips official screening, certifies the complete native
design, and requires every absorbed level to have at least two observations in
every XPO training sample; otherwise it returns `r(498)`. The commands never
sort observations internally.

`backend(auto)` never silently opts into native XPO. A large native XPO run
therefore stops and asks for explicit `backend(mata) frontierok`. Native PO
plugin/cluster retains its full-sample singleton guard and otherwise supports
balanced and unbalanced clusters. Legacy one-way PO calls may still specify
`frontierok`; it is a deprecated no-op there and a required opt-in for
two-way PO.

Historical rc6 work closed the then-measured XPO selection-path differences.
Duplicate absorbed
partitions are compacted before selection, and coordinate descent profiles the
intercept with centered Hessians. The K=3/5/7/10 plugin/BIC by robust/cluster
matrix now matches official Stata to `1.776e-15` in coefficients and
`8.674e-18` in covariance; `breathe` XPO/BIC/robust matches to `7.550e-15`
and `1.388e-16`. These are retained historical tests, not a blanket current
claim for every XPO design.

The root causes, implementation changes, and remaining safety boundary are
summarized in `XPO_PARITY_FIX_20260715.md`.

Without `foldcluster()`, native Mata XPO retains official-compatible
observation-level outer folds even under `vce(cluster ...)`. rc11 adds the
fold-local official screen described above. Historical rc7 work added explicit
`foldcluster()` research modes: seeded one-way cluster-disjoint folds and
two-way `K^2` folds with Chiang--DoubleML moment and score-covariance
normalization. In a two-way score cell, training excludes rows sharing either
held-out cluster partition. See `CLUSTER_DML_XPO_20260716.md` and the
[DoubleML example](https://docs.doubleml.org/stable/examples/py_double_ml_multiway_cluster.html).
`foldcluster()` has no official fold-screen oracle and remains a separate
frontier contract; `strict(production)` is preferred. Native XPO remains
frontier because broad nominal-coverage evidence is incomplete, two-way
`K^2` BIC can be expensive, and plugin fits can still fail their loading-update
certificate. The current 30k scale cells pass; the two-way BIC cell took
3,058--3,086 estimator seconds per repeat.

The executable Python reference is `python/example_multiway_cluster.py`. It
uses `DoubleMLData(cluster_cols=["g1", "g2"])` and
two treatments in `DoubleMLPLR(n_folds=3)`, verifies nine score folds and a
full 2-by-2 score covariance, and maps directly to the
Stata `foldcluster(g1 g2) xfolds(3)` geometry. This is an alignment of split and
score contracts, not a numerical parity claim: the reference uses LassoCV,
while native Stata examples use plugin or BIC.

## Experimental Python backend

Install `python/requirements-experimental.txt`, which pins DoubleML 0.11.3 and
the tested scientific-Python stack. Then either pass a portable executable
with `python()` or configure Stata's `c(python_exec)`. The wrapper finds the
backend script on the ado path or, in a source checkout, in the sibling
`python/` directory.

After `net install`, locate the installed requirements and copy-ready public
fixture without assuming an installation directory:

```stata
findfile requirements-experimental.txt
display "`r(fn)'"
findfile validate_public_examples.do
local examples "`r(fn)'"
do "`examples'" "`c(python_exec)'" required
```

This is a direct DoubleML PLR contract with scikit-learn LassoCV nuisance
learners. It accepts multiple treatments and returns the familiar Stata
results: all treatment coefficients are in `e(b)`, and the complete symmetric
cross-treatment covariance is in both `e(V)` and `e(sparse_covariance)`.
Off-diagonal treatment covariances are estimated, not filled with zero;
`e(sparse_covariance_full)=1` records that contract.

The Python source has 18 tests covering robust, one-way, and two-way
independent score-sum covariance oracles; symmetric positive-semidefinite
two-treatment covariance with nonzero off-diagonal entries; PO and XPO paths;
cluster-disjoint K and K^2 score folds; ragged-design reproducibility; explicit
fold-count/empty-cell refusals; unbalanced clustered-PO coefficient invariance;
rank-deficient covariance refusal; and zero-Jacobian rejection. The
standalone Stata public preflight covers all PO/XPO robust, one-way, and two-way
examples, including invalid factor-treatment, interaction, always-in, and
stale-`e()` cases. The standalone suite and the final hash-bound full
certification both pass.

```stata
local py "C:/path/to/python.exe"

sparse_poreg y d1 d2, absorb(a b) controls(x1 x2) ///
    selection(lassocv) learner(lassocv) vce(robust) ///
    lassocvfolds(3) lassocvalphas(30) ///
    backend(python) python("`py'")

sparse_poreg y d1 d2, absorb(a b) controls(x1 x2) ///
    selection(lassocv) learner(lassocv) vce(cluster g1) ///
    lassocvfolds(3) backend(python) python("`py'")

sparse_poreg y d1 d2, absorb(a b) controls(x1 x2) ///
    selection(lassocv) learner(lassocv) vce(cluster g1 g2) ///
    lassocvfolds(3) backend(python) python("`py'")

sparse_xporeg y d1 d2, absorb(a b) controls(x1 x2) ///
    selection(lassocv) learner(lassocv) vce(robust) ///
    lassocvfolds(3) lassocvalphas(30) xfolds(3) rseed(1103) ///
    backend(python) python("`py'")

sparse_xporeg y d1 d2, absorb(a b) controls(x1 x2) ///
    selection(lassocv) learner(lassocv) vce(cluster g1) ///
    foldcluster(g1) lassocvfolds(3) xfolds(3) rseed(1103) ///
    backend(python) python("`py'")

sparse_xporeg y d1 d2, absorb(a b) controls(x1 x2) ///
    selection(lassocv) learner(lassocv) vce(cluster g1 g2) ///
    foldcluster(g1 g2) lassocvfolds(3) xfolds(3) rseed(1103) ///
    backend(python) python("`py'")

matrix list e(V)
```

Python PO supports robust, one-way, and two-way clustered inference. Its
full-sample nuisance fits and coefficient are observation-weighted and do not
change with `vce()`, cluster variables, or `rseed()`; clustering changes only
the covariance and no K-way cluster partition is used. One-way PO uses a CR0
score sandwich and two-way PO uses
the conventional CGM marginal-plus-marginal-minus-intersection sandwich,
without Stata small-cluster corrections. The fits are not cross-fit
(`e(nuisance_crossfit)=0`, `e(crossfit_unit)="none"`).
By contrast, Python XPO uses observation-level outer folds under `vce(robust)`.
For that route, `e(crossfit_unit)="observation"`,
`e(nuisance_crossfit)=1`, `e(observation_disjoint)=1`, and
`e(n_xfolds)=e(n_score_folds)` expose the audited split. Generic
`e(fold_coverage_*)` and `e(fold_min_*)` scalars report its coverage and minimum
train/test sizes. `e(cluster_disjoint)` is not posted because no cluster
dimension exists. Clustered XPO
requires `foldcluster()` to match the one or two variables in `vce(cluster
...)`. DoubleML partitions each cluster dimension K ways, where K is
`xfolds()`: one-way clustering therefore has K score folds and two-way
clustering has the K^2 Cartesian-product score cells. For every two-way cell,
training excludes observations sharing either held-out cluster partition.
The bridge audits row coverage and the cluster-disjoint training contract.
For clustered XPO, `e(cluster_disjoint)=1` and
`e(n_score_folds)=xfolds()^e(n_cluster_dims)` make that geometry visible.
Here K is not silently reduced: every cluster dimension needs at least K unique
values, and every realized Cartesian score/training cell must be nonempty.
Reduce `xfolds()` when sparse cluster intersections violate that gate.

The cluster-disjoint guarantee concerns the outer score folds. LassoCV's inner
hyperparameter CV remains row-wise, so it can place rows from one cluster in
different inner folds. This limitation is explicit and is one reason the
backend remains experimental.

Clustered Python XPO follows DoubleML 0.11.3 semantics. Its two-way covariance
is the Chiang--DoubleML positive sum of the marginal cluster-score components
with fold normalization, not PO's conventional CGM sandwich. Neither path adds
Stata's small-cluster finite-sample corrections. Rank-deficient treatment
covariance is refused before Stata results are posted.

Only `selection(lassocv) learner(lassocv)` is implemented on this path. Python
is not an implementation of Stata's plugin or BIC penalty path, so matching
syntax, folds, `e(b)`, and full `e(V)` output does not imply official
`poregress`/`xporegress` numerical parity.

The bridge accepts plain or `c.` numeric controls and maps `i.` or `ibn.`
categorical controls to absorbed indicators. Treatments must be plain numeric
variables. Factor treatments, `#`/`##` interactions, parenthesized always-in
controls, and nonnumeric treatments, controls, or absorbed IDs are rejected.
Unlike the Mata backend, Python does not post exact selected absorbed IDs or
`e(controls_sel)`. LassoCV warnings and iteration metadata are also not a
production convergence gate; these are explicit experimental limitations.

rc8 adds conventional two-way clustered PO inference. It uses marginal-plus-
marginal-minus-intersection covariance, Stata-style PSD repair, and
`min(G1,G2)-1` degrees of freedom. Nuisance selection is observation-weighted,
so this is a frontier inference extension rather than official `poregress`
parity. See `PO_MULTIWAY_20260716.md`.

## Reliability diagnostics

Successful Mata fits post:

- convergence, lasso calls, iterations, and sweeps;
- coordinate-descent full/active sweeps, inactive KKT scans/fallbacks, maximum
  active size, and block-Gram polish calls/successes;
- final-solution KKT maximum, relative maximum, and pass status;
- plugin convergence and iterations, plus `e(plugin_pass)`,
  `e(plugin_support_stable)`, `e(plugin_capped_calls)`, and
  `e(plugin_load_diff_max)`;
- XPO support minima and mode flags: `e(min_abs_support)`,
  `e(min_train_support)`, `e(training_support_safe)`, and
  `e(training_support_relaxed)`;
- Stata-screen metadata: `e(sparse_xcheck)`, `e(kkt_scope)`,
  `e(xcheck_seconds)`, and kept/dropped counts or foldwise ranges;
- selected-set sizes, `e(controls_sel)`, `e(k_controls_sel)`, and raw absorbed
  IDs in `e(sparse_support)` when the union fits the current `matsize`;
- native per-fit targets, folds, counts, and supports in
  `e(selection_counts)`, `e(lasso_depvar_#)`, `e(lasso_sel_p_#)`, and
  `e(lasso_selected_#)`; `e(selection_counts_complete)=0` omits only the
  matrix when the fit count exceeds `c(matsize)`;
- `e(controls_sel_complete)` so macro-length overflow is explicit rather than
  fatal;
- BIC path-ending counts (`e(bic_identified)`, `e(bic_stop)`, and
  `e(bic_grid)`);
- nuisance and final-rank diagnostics plus condition numbers;
- finite-result, backend, raw/effective feature-count, and kernel metadata;
- multiway cluster counts, raw minimum VCE eigenvalue, PSD-repair status, and
  covariance rank when applicable.

Iteration exhaustion returns `r(430)` before estimates are posted. A
rank-deficient final treatment regression or nonfinite result also fails.
Nuisance rank deficiency is recorded because fitted values can remain
identified despite redundant controls.

KKT is a numerical optimality diagnostic, not a proof of official parity or
valid inference. The test suite checks those properties separately.

`e(plugin_iterations)` counts loading-update iterations across nuisance fits;
the default per-fit cap is 15. It is not a selected-control count.

`e(controls_sel)` is the union: a firm/person level appears when at least one
outcome or treatment nuisance lasso selects it (and, for XPO, in any fold).
It is not the number selected by every regression.  After a native Mata fit,
`estat selection` displays the individual nuisance counts and
`estat selection, support` adds their selected terms.  XPO rows are
fold-major; two-way cluster-disjoint fits report both fold dimensions.
For an official `poregress` or `xporegress` result, the corresponding display
is `lassoinfo, each`.

PO profiles have a hard KKT gate. Every native
`selection(plugin) strict(production)` fit also requires `e(plugin_pass)=1`;
a capped plugin call
that neither converges nor passes its cold support/KKT certificate returns
`r(430)`. XPO additionally applies its training-support guard. Only the
experimental clustered-XPO plugin route under `strict(stata)` may return with
a warning and `e(kkt_pass)=0` or `e(plugin_pass)=0`; production refuses it.

## Current rc13 evidence

Full run `rc13_cert_full_mit4_20260719` passes with process exit 0 and stable
source hashes. Its 96 cells contain 94 production-exact matches, one
intentional refusal, and one diagnostic frontier; maximum production `b` and
`V` differences are `9.104e-15` and `3.903e-18`. Quick run
`rc13_cert_quick_coauthor4_20260719` binds the one-rule Mata computation
default, coauthor quick start, manual-data examples, and final status text to
the unchanged release source.

The current routing run `rc13_route_kernelauto_20260719` completed in 1353.481
seconds: 144/144 isolated engine cells, 48/48 paired comparisons, 16/16
guidance rows, 21/21 stable source rows, and 144/144 support files. It has
`overall_pass=1`, no timeouts, hard failures, parity failures, or safety
failures, and all comparisons are exact. The maximum coefficient and covariance
differences are `5.99520433297585e-15` and `1.00000709752843e-17`. Evidence is
in the [routing manifest](../logs/benchmark_routing_manifest_rc13_route_kernelauto_20260719.json)
and matching `raw`, `resources`, `pairs`, `guidance`, and `sources` CSVs.

The Python 18-test suite and core Stata public preflight pass in `python.log`,
including the factor-treatment, interaction, always-in, and stale-`e()` gates.
The fixed-point state unit test passes exact fixed/cycle/budget classifications.
The source-stable trajectory gate `rc13_fixedpoint_pobicycle_20260718` also
passes. Its two exact-repeat 10k runs classify nine cycles with zero solver,
KKT, verification, or exhaustion failures in 24.388 and 25.077 seconds. The
30k run verifies ten fixed points and five primitive period-2 cycles in 77.589
seconds. All seven source rows are stable. Evidence is in the
[trajectory runner](../logs/probe_xpo_plugin_fixedpoint_rc13_fixedpoint_pobicycle_20260718_runner.csv)
and matching source/per-run CSVs.

The current-source public `strict(production)` refusal tests now pass at both
scales. The 10k case returns `r(430)` with full/training support `3/2`,
solver/KKT/rank `1/1/1`, and plugin pass/capped calls/iterations `0/9/135`.
The 30k case returns `r(430)` with support `5/4`, the same `1/1/1` safety
diagnostics, and `0/5/99` in 58.503 seconds. See the
[10k refusal log](../logs/validate_xpo_plugin_cycle_refusal.log) and
[30k refusal log](../logs/validate_xpo_plugin_cycle_refusal_30k.log). The final
full certification binds these tests to its unchanged source snapshot.

The gate never invents a fixed point. `strict(production)` refuses a detected
cycle; opt-in `strict(stata)` frontier use may preserve the finite step-15 path
with a warning. A standalone engineering probe shows that the nonmutating cold
KKT rescue leaves support unchanged while reducing relative KKT error from
`5.365761e-4` to `2.096615e-14`; release status relies on the trajectory and
public-refusal gates. The final full-certification RunId is
`rc13_cert_full_mit4_20260719`; its manifest and summary are complete.

The corrected full release matrix `rc13_scale_full_pobicycle_20260718` passes:
40/40 shards, 80/80 cells, 40/40 exact-repeat groups, 35 stable source rows,
76 ordinary passes, and four expected production refusals. Its 16 official
parity cells have maximum `b` and `V` differences `4.663e-15` and `1.016e-18`.

| current-source scale profile | repeat time | boundary |
|---|---:|---|
| PO/plugin robust, one million controls | 117.394--122.798 s | production pass |
| PO/plugin cluster, one million controls | 118.236--118.857 s | production pass |
| PO/BIC robust, one million controls | 779.426--805.182 s | production pass |
| PO/BIC cluster, one million controls | 613.192--617.038 s | production pass |
| XPO/BIC one-way cluster-disjoint, 30k controls | 238.910--240.547 s | production pass |
| XPO/BIC two-way `K^2`, 30k controls | 3,057.789--3,086.429 s | feasible, not a fast default |

Peak RSS is 1.012 GiB for the production family and 0.354 GiB for the native
family. A separate diagnostic refusal reaches 26.922 GiB and is excluded from
production recommendations. See the
[release manifest](../logs/benchmark_indexed_scale_release_rc13_scale_full_pobicycle_20260718_manifest.csv)
and [validated union](../logs/benchmark_indexed_scale_release_rc13_scale_full_pobicycle_20260718_union.json).
The old unsafe-order 10k/30k artifacts remain historical only.

### Current-source all-firm XPO-BIC/robust diagnostic

The paired run `rc13_devstop_xpobic_pair_r5_20260718` uses 11,830 observations,
455 firms, five canonical folds, and `rseed(1103)`. Correcting the deviance
stop to `|RSS_i-RSS_(i-1)|/TSS < 1e-5` gives 69/70 exact fold-target supports
and an exact 481-control final union. The only difference is
`20757.PERMCO` in the fold-2 `lnppent_emp` nuisance fit; dropping it fails the
native KKT gate, so the safe result retains it.

| sparse s | official s | speedup | max `b` diff | max SE diff | max `V` diff |
|---:|---:|---:|---:|---:|---:|
| 442.235 | 569.694 | 1.29x | `4.048e-5` | `1.699e-6` | `1.515e-7` |

Both engine cells and all artifacts are hash-bound and valid. The deliberately
exact pair gate rejects this diagnostic because one fold support differs; it
is close-parity evidence, not exact release evidence.

### Current-source full empirical XPO evidence

The final-source full-sample XPO-BIC/cluster comparisons are exact:

| sample | sparse s | official s | selected | max `b` / `V` diff | result |
|---|---:|---:|---:|---:|---|
| patenting | 264.018 | 2,052.546 | 274 | `6.700e-14` / `8.601e-15` | exact, live 7.77x |
| no-missing-R&D | 380.778 | 3,862.323 | 429 | `1.961e-13` / `5.201e-15` | exact, live 10.14x |
| all-firm | 1,137.910 | 15,125.485 | 1,737 | `1.402e-14` / `3.500e-16` | exact; validated official anchor |

The all-firm live official cell exceeded its two-hour cap. Its completed
official-only artifact and 4h12m timing were cryptographically revalidated
against the unchanged Stata executable, official ado files, data, manifests,
and artifact hashes. This proves current sparse numerical parity, but its
indicative 13.29x timing ratio is not labeled a synchronized speedup. The
three run IDs are `rc13_emp_xpo_pat_full_biccl2_20260719`,
`rc13_emp_xpo_nomiss_full_biccl_20260719`, and
`rc13_emp_xpo_full_biccl_20260719`; every manifest has all release gates set.

Minimum training-fold support is zero in all three full samples. These are
therefore `strict(stata)` compatibility results, not production-safe XPO
claims. `strict(production)` is expected to refuse these folds with `r(498)`.

### Retained rc12 empirical and inference evidence

The hash-bound empirical run `rc12_emp_final_20260717_051000` uses BIC,
one-way clustered XPO, five folds, and `rseed(1103)`. All six safe and full
pairs have exact support and numerical parity:

| scope/sample | sparse s | official s | speedup | max `b` / `V` diff | selected | status |
|---|---:|---:|---:|---:|---:|---|
| safe/allfirm | 21.145 | 73.567 | 3.48x | `2.753e-14` / `4.413e-15` | 50 | exact support |
| safe/patenting | 17.726 | 50.592 | 2.85x | `1.754e-14` / `5.884e-15` | 37 | exact support |
| safe/nominssrd | 7.638 | 25.482 | 3.34x | `1.110e-14` / `9.659e-15` | 19 | exact support |
| full/patenting | 189.060 | 1,781.973 | 9.43x | `6.756e-14` / `8.507e-15` | 274 | exact support |
| full/nominssrd | 248.440 | 3,410.070 | 13.73x | `1.965e-13` / `5.156e-15` | 429 | exact support |
| full/allfirm | 683.101 | 15,125.485 | 22.14x | `1.443e-14` / `3.563e-16` | 1,737 | exact support |

The three full sparse fits pass convergence, finite, rank, and KKT gates, but
their minimum training support is zero. They therefore use relaxed
`strict(stata)` compatibility handling and are not `strict(production)`
evidence. See `benchmark_empirical_xpo_rc12_emp_final_20260717_051000.csv` and
its metadata, runner, and artifact index under `logs/`. The final master and
runner each contain 18 clean rows; all 75 indexed artifact hashes validate.
The full/allfirm support comparison is exact, with zero support-only tokens on
either side.

The replacement current-source harness is
[`run_benchmark_empirical_matrix.ps1`](../tests/run_benchmark_empirical_matrix.ps1).
Its concise pair output records sparse and official seconds, speedup, maximum
`b`/`V` differences, full support-set overlap, canonical fold/safety contracts,
and source/artifact hashes. A one-row final manifest closes each run. The
contract and non-estimation preflight are documented in
[`EMPIRICAL_RELEASE_MATRIX.md`](EMPIRICAL_RELEASE_MATRIX.md). The three
requested full-sample XPO-BIC/cluster runs now have `release_evidence=1`,
`overall_run_gate=1`, and `overall_release_gate=1`; the broader 90-cell and
144-cell empirical campaigns remain optional.

The current-source 500-repetition inference diagnostic completes 5,999 of
6,000 attempted fits. The sole failure is the intentional XPO/plugin/two-way
production `r(430)` for an unstable capped loading update. Coverage ranges
from 0.234 to 0.998 across profiles; this one `N=480` heteroskedastic DGP is a
diagnostic, not nominal-coverage certification. Targeted PO/XPO BIC cluster1
comparisons are numerically exact against official Stata and share the low
coverage. Evidence is in
`coverage_inference_rc12_full_20260717_043000.csv`, its manifest/runner/artifact
files, `probe_coverage_bic_cluster1_summary.csv`, and
`probe_coverage_xpo_plugin_cluster2_failure.log`.

At 10k features, the rc12 non-mutating PO probe remains historical evidence.
Its XPO discrepancy is not a current parity boundary: the corrected safe-order
probe above shows machine-precision `b,V` parity. The remaining XPO
plugin/cluster boundary is the failed loading-cycle certificate and required
production refusal, not coefficient or covariance drift.

## Earlier rc11 evidence - historical

The following artifacts were produced before the planned rc11 version/hash
freeze. They preserve implementation lineage but are not current rc13 release
evidence.

| artifact | result | boundary |
|---|---|---|
| `benchmark_empirical_xpo_20260716_230600.csv` | all three safe empirical samples have exact support; maximum `b`/`V` differences are `2.753e-14`/`9.659e-15`; sparse runs are 2.9--3.6x faster | ordinary-fold `strict(stata)` only |
| `benchmark_empirical_xpo_20260717_004000.csv` | all three full sparse fits pass convergence, finite, rank, and KKT gates | no paired full official results yet; relaxed training support used |
| `benchmark_indexed_scale_20260717_001402.csv` | unsafe fold-order artifact; its approximate XPO plugin differences are superseded by the corrected safe-order exact results | pre-freeze harness lineage only |
| `benchmark_indexed_scale_20260717_000825.csv` | all six native one-/two-way cluster scale cells pass truth, finite, rank, KKT, and production plugin gates | pre-freeze 10k native run |
| `coverage_inference_20260717_010202.csv` | 50-rep pilot completed; coverage varies materially by profile | diagnostic pilot, not certification |
| `probe_coverage_bic_cluster1_summary.csv` | sparse and official PO/XPO BIC cluster1 results are numerically exact; both share the pilot's poor coverage | identifies DGP/estimator behavior, not a sparse parity defect |
| `probe_coverage_xpo_plugin_cluster2_failure.log` | sequential isolation reproduces the single rep-39 capped plugin failure as intentional `r(430)` under production safety | safety gate worked as designed |

Older rc6--rc8 evidence below is retained as historical implementation
lineage.

- rc6 centered-CD and duplicate-partition gates pass. The official parity
  matrix is in `validate_indexed_full_16_Jul_2026_000913.log`; the final
  `breathe` XPO/BIC/robust probe is
  `breathe_xpo_bic_robust_16_Jul_2026_000412.log`; and the 30,002-feature
  benchmark passes in `benchmark_indexed_large_16_Jul_2026_001007.log`.
- The historical rc6 quick tier passes in `run_certification.log`. NewB-newE add 32
  supported official/native cells with maximum coefficient and covariance
  differences `4.885e-15` and `2.277e-18`; newA records the intended XPO
  singleton-fold refusals.
- The rc7 quick tier also passes the independent one-/two-way cluster-DML
  covariance oracle, seeded public `foldcluster()` fits, cluster-order and
  row-order checks, and algebraically valid/invalid training-support fixtures.
  Historical logs are `validate_cluster_disjoint_xpo_16_Jul_2026_114027.log` and
  `validate_xpo_support_guard_16_Jul_2026_113322.log`.
- The rc8 PO two-way oracle passes plugin/BIC, both cluster orders, full
  two-treatment covariance, equivalent partitions, metadata, and refusal
  cases in `validate_po_multiway_vce_16_Jul_2026_121329.log`. The 30,002-feature
  robust and two-way PO smoke passes in
  `benchmark_indexed_large_16_Jul_2026_121450.log`.
- The clustered-plugin residual/loading rule is identified independently by
  `probe_p7_continuous_oracle.do`: all seven active continuous-column checks
  match official loadings within `1.22e-14`. Unbalanced newA/newB/newD PO and
  supported XPO comparisons now match official coefficients and covariance at
  machine precision. See `CLUSTER_PLUGIN_FIX_20260715.md`.
- The historical rc4 solver-focus gate remains in `validate_rc4_focus.log`.
- `benchmark_indexed_scale_20260715_203958.csv` passes frozen 10k
  plugin/robust and BIC/cluster cells. Plugin matches official Stata to
  `1.665e-15` in coefficients and `8.403e-19` in covariance.
- rc4 replaces thousands of slow clustered active sweeps with guarded
  block-Gram polish and a vectorized inactive KKT scan.

| PO BIC/cluster | estimator | wall | peak RSS | KKT fallback |
|---|---:|---:|---:|---:|
| 500k features | 275.279 s | 276.664 s | 594.5 MiB | 0 |
| 1m features | 579.990 s | 582.568 s | 947.5 MiB | 0 |

Both rows pass convergence, finite, rank, covariance, KKT, and truth `<0.05`
gates, and exactly reproduce prior `b`, `V`, and solver counters. The old
500k clustered-BIC timeout is resolved. Plugin/robust remains the recommended
very-large profile based on the current scale evidence; the frozen rc4 rows
remain historical supporting evidence.

The scale driver isolates cells in fresh processes, validates result identity,
supports per-cell timeouts, measures RSS, and checks exact repeatability. The
corrected 10k/30k and general scale cells pass the hash-bound release harness.
The full certification and final documentation binding also pass; the release
builder creates and independently verifies the archive checksum.

## Reproduce

```powershell
$trajectoryRunId = 'rc13_fixedpoint_pobicycle_20260718'
.\tests\run_certification.ps1 -Tier quick -RunId rc13_cert_quick_YYYYMMDD_HHMMSS -RoutingRunId ROUTING_RUN_ID
.\tests\run_probe_xpo_plugin_fixedpoint_trajectory.ps1 -RunId $trajectoryRunId
.\tests\run_benchmark_routing.ps1 -Mode full -RunId rc13_route_YYYYMMDD_HHMMSS
.\tests\run_certification.ps1 -Tier full -RunId rc13_cert_full_YYYYMMDD_HHMMSS -RoutingRunId ROUTING_RUN_ID -TrajectoryRunId $trajectoryRunId
.\tests\validate_benchmark_indexed_scale_release_plan.ps1
.\tests\run_benchmark_indexed_scale_release.ps1 -PlanOnly
.\tests\run_benchmark_indexed_scale_release.ps1 `
    -RunId rc13_scale_release_YYYYMMDD `
    -TrajectoryRunId $trajectoryRunId `
    -CellTimeoutSeconds 7200
# After an interruption, repeat the preceding command with -Resume.
powershell -NoProfile -ExecutionPolicy Bypass -File .\tests\validate_empirical_matrix_harness.ps1
.\tests\run_benchmark_empirical_matrix.ps1 -Scope safe,full `
    -Engine sparse,production,official -AllowOfficialFull `
    -RequireParity -PlanOnly
# Build and independently verify the MIT-licensed release archive:
.\scripts\build_release.ps1
.\scripts\verify_release.ps1 -Archive .\dist\sparse_poreg-0.8.4-rc14.zip
```

The scale release launcher expands the checked-in nonrectangular specification
to exactly 40 serial shards and 80 cells. It freezes the plan and executable
source hashes, classifies parity/native/refusal evidence, and validates the
completed shard union. Use `run_benchmark_indexed_scale.ps1` directly only for
diagnostic subsets.

See `CERTIFICATION_20260715.md` and `ROADMAP_20260715.md` for current
boundaries. `CANDIDATE_RC2_20260715.md`, `CANDIDATE_RC3_20260715.md`,
`CANDIDATE_RC4_20260715.md`, and `LEGACY_VALIDATION_20260716.md` are frozen
historical snapshots. The Python LassoCV backend remains experimental: it
posts full treatment covariance but is not an official plugin/BIC substitute.

The package is distributed under the MIT license in the root `LICENSE`. The
release builder uses an explicit allowlist and excludes the three empirical
`.dta` files under `tests/`; they are test evidence, not public payload.



# sparse_poreg 0.8.4-rc14: linear PRL and XPO

This independent Stata package provides `sparse_poreg` (linear post-regularization/PRL) and `sparse_xporeg` (linear cross-fit partialing-out/DML). Both have official Stata and sparse Mata backends. This ZIP contains the unchanged 22 July 2026 `0.8.4-rc14` program files plus a public-data example, anonymized logs, and comparison notes. It does **not** include research data, saved estimates, or the separate developmental Poisson command. It is not a StataCorp product.

## Install

Extract the ZIP. In Stata, change directory to the extracted folder containing `sparse_poreg.pkg` and run:

```stata
net install sparse_poreg, from("`c(pwd)'") replace
which sparse_poreg
which sparse_xporeg
do "tests/finance_quickstart.do"
```

The commands require Stata 16 or newer; the empirical logs here used Stata 19. The quick check uses generated data and Stata 19. The optional Python backend has separate requirements in `python/requirements-experimental.txt`.

## Paper's public linear example: sparse Mata first

Download `REPLIC~1.DTA` from the [paper replication repository](https://github.com/persistentferesearch-sys/persistent_fe_remendy), then run:

```stata
use "C:/path/to/REPLIC~1.DTA", clear
sort PERMCO fyear
local y lnnpatent
local x RDAT lnME RD_missing lnAge lnK2L TobinQ ROA Leverage ///
    CASHAT KZidx InstOwn oms_HHidx oms_HHidx_square

* Sparse Mata PRL/PO, parity-oriented mode.
sparse_poreg `y' `x', absorb(PERMCO fyear) ///
    vce(cluster PERMCO) selection(bic, gridminok) rseed(1103) ///
    backend(mata) strict(stata)

* Sparse Mata DML/XPO, experimental frontier.
sparse_xporeg `y' `x', absorb(PERMCO fyear) ///
    vce(cluster PERMCO) selection(bic, gridminok) rseed(1103) ///
    xfolds(5) backend(mata) frontierok strict(stata)
```

The `x` variables are variables of interest. Firm/year indicators are selection controls. The copy-ready [Mata do-file](examples/github_rdpatent_mata_prl_xpo.do) runs both commands. For clustered BIC PRL, this release rejects `frontierok`; the successful PRL run used `backend(mata) strict(stata)`. The XPO route requires `frontierok`. That flag is an experimental-route acknowledgement, not an accuracy setting.

To run the exact official Stata routes through this package, use `controls(i.PERMCO i.fyear)` with `backend(official)` in place of `absorb(...)` and `backend(mata) ...`. The earlier [official comparison do-file](examples/github_rdpatent_linear_po_xpo.do) includes those commands. Direct official Stata commands are `poregress` and `xporegress` with the same model options. Full official fits can take hours.

[COMPARISON.md](COMPARISON.md) gives results, timing, command-by-command evidence, and parity limits. `results/` contains complete package official-PRL and Mata-PRL/XPO logs, the paper's published official PRL/XPO logs and CSVs, and the XPO printed-precision audit. The included official XPO log is the **paper's published reference**, not a package-produced run. The Mata XPO route still requires `frontierok`.

See `docs/README.md` for backend boundaries and `UPLOAD.md` for self-upload instructions. `SHA256SUMS` covers the ZIP contents; the adjacent `.zip.sha256` checks the whole archive. The package source is MIT licensed.
