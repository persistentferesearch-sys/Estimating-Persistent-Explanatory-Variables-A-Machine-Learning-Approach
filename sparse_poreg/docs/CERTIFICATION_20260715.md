# Certification record - updated 2026-07-22

Status: `0.8.4-rc14` release.  The rc13 certification baseline is retained.
rc14 adds official-compatible omission of collinear PO variables of interest
in the native Mata final regression.

## RC14 PO omitted-target fix

Official `poregress` omits an unidentified variable of interest and continues
with the estimable coefficients.  Native PO previously returned `r(498)` for
the entire model.  rc14 now keeps a full-rank target subset, zero-fills the
omitted coefficient and covariance row/column in their original order, and
posts `e(k_varsofinterest_omitted)` and
`e(varsofinterest_omitted)`.

| rc14 gate | evidence | status |
|---|---|---|
| synthetic robust/one-way/two-way omission | `validate_po_omitted_target.log` | **PASS**; retained robust and one-way `b,V` match official to `1e-10`; two-way output is finite and rank-valid |
| indexed kernel | `validate_indexed_kernel_22_Jul_2026_095202.log` | **PASS** |
| PO multiway covariance | `validate_po_multiway_vce_22_Jul_2026_095205.log` | **PASS** |
| current breathe matrix | `validate_breathe_current_22_Jul_2026_095318.log` | **PASS** |
| historical 56-cell matrix | `validate_matrix_22Jul2026_095319.log` | **PASS** |
| public official/Mata/Python examples | `_run_public_examples_current.log` | **PASS** |
| four formerly failed ownership jobs | launcher tag `po_omitted_fix_20260722` | **PASS 4/4**; `sigdum_hhp` omitted as official Stata does |

Three of those four ownership coefficients match official Stata to numerical
precision.  The all-named firm+person-FE BIC coefficient remains a separate
selection-path frontier (`0.2156234` native versus `0.2250702` official); the
rank/omission failure itself is fixed.

## RC13 release evidence

Release evidence must share the same runtime hashes. Each harness manifest
also freezes its runner, driver, Stata executable, official ado files, and
inputs; a single undifferentiated "source hash" is not sufficient.

| rc13 gate | run ID / evidence | status |
|---|---|---|
| final source-binding quick certification | `rc13_cert_quick_coauthor4_20260719` | **PASS**; one-rule Mata computation default, coauthor quick start, manual-data help examples, and final status text bound to unchanged release source |
| prior synthetic routing profile | `rc13_route_full_20260717_1405` | complete supporting evidence; predates current runtime |
| current-runtime routing profile | `rc13_route_kernelauto_20260719` | **PASS**; 144 cells, 48 pairs, 16 guidance rows, all exact, 21 stable sources |
| Mata computation-method check | `benchmark_kernel_choice.csv` | **PASS**; one automatic indexed method for users, exact support/numerical parity, and less than one second lost in the only legacy-faster case |
| Python/Stata public preflight | `python.log`; Python 18-test suite | **PASS** in final hash-bound full tier |
| fixed-point state unit test | `validate_plugin_fixedpoint_state.log` | PASS |
| public production cycle refusal | `validate_xpo_plugin_cycle_refusal.log`; `validate_xpo_plugin_cycle_refusal_30k.log` | **PASS** standalone at 10k/30k and bound by final full tier |
| 10k/30k fixed-point trajectory | `rc13_fixedpoint_pobicycle_20260718` | **PASS**; exact 10k repeat, verified 30k fixed points/cycles, zero failures, seven stable sources |
| full release-scale matrix | `rc13_scale_full_pobicycle_20260718` | **PASS**; 40/40 shards, 80/80 cells, 40 repeat groups, 16 exact official-parity cells, 76 passes, four expected production refusals, stable sources |
| final full certification | `rc13_cert_full_mit4_20260719` | **PASS**; exit 0, source stable, 96 cells, 94 production exact, one intentional refusal, one diagnostic frontier |
| corrected 10k safe-order pilot | `benchmark_indexed_scale_safeorder_probe_20260717.csv` | superseded supporting evidence; exact robust/cluster `b,V` is now included in the release matrix |
| old 30k cycle scale shard | `rc13_scale_cycle30_final_20260717_1721` | **SUPERSEDED** unsafe-order artifact |
| retained empirical XPO parity/speed | `rc12_emp_final_20260717_051000` | prior-runtime evidence; six exact-support matches, 2.85x--22.14x faster |
| current all-firm XPO-BIC/robust diagnostic | `rc13_devstop_xpobic_pair_r5_20260718` | valid paired cells; 1.29x faster, exact final union, 69/70 fold supports; one KKT-required extra dummy, so exact pair gate correctly rejects |
| current all-firm PO-BIC/robust parity | `rc13_emp_pobic_safe_final_20260719` | **PASS**; 481/481 exact support, max `b`/`V` differences `3.297e-14`/`7.008e-16`; official remains faster on this manageable design |
| current full XPO-BIC/cluster parity | `rc13_emp_xpo_*_full_biccl*_20260719` | **PASS**; all three empirical pairs have exact support and numerical parity; two live timing pairs are 7.77x and 10.14x faster; all-firm uses a hash-validated official oracle after the live two-hour timeout |
| retained inference diagnostic | `rc12_full_20260717_043000` | prior-runtime evidence; 5,999/6,000 fits and one classified production refusal |
| install metadata and placement | `stata.toc`, `sparse_poreg.pkg`; `validate_net_install_20260717_215350.log` | **PASS**; isolated PLUS tree finds both commands, `_sparse_estat`, Python files, and public fixture; bridge resolver uses installed `PLUS/py`, not the checkout |
| source/archive manifest | `certification_rc13_manifest_rc13_cert_full_mit4_20260719.csv`; release `SHA256SUMS` | **PASS**; source manifest complete; archive builder verifies every payload hash |

Current rc13 runtime hashes for the final certification are:

| runtime component | SHA-256 |
|---|---|
| `stata/_sparse_prl_xpo_mata.ado` | `E8F3CD1B5C66F5A5CF40F20AC4108C9A356D552505BB306E64F04940A6E6F43C` |
| `stata/sparse_poreg.ado` | `8A5496F1BD64092AEEBA2B9EF4B606CC87CA3024D78320E4441C2CA52ACD7699` |
| `stata/sparse_xporeg.ado` | `5BEB1A38108E24776A655C6C90EF07441C3AECB3C748750870A54313F4C013E6` |
| `python/sparse_prl_xpo_backend.py` | `DBF947208A0D3B3CD20E3885B335C2ED702308C18250BD5FA2C20D86F2B0F765` |

## Current-source routing evidence

The source-stable run `rc13_route_kernelauto_20260719` completed in 1353.481
seconds with 144/144 engine rows, 48/48 paired comparisons, 16/16 guidance
rows, 21/21 source rows, and 144/144 support files. It has `overall_pass=1`, no
timeout, hard, source-drift, parity, or safety failure, and every comparison is
exact. Maximum coefficient and covariance differences are
`5.99520433297585e-15` and `1.00000709752843e-17`.

| profile | documented low-density anchor | documented high-density anchor | interpretation |
|---|---|---|---|
| PO plugin (robust and cluster) | `N=1000, P=112, P/N=.112, N*P=112000` | `N=1000, P=412, P/N=.412, N*P=412000` | observed sparse-performance anchor |
| PO BIC, robust | `N=8000, P=812, P/N=.1015, N*P=6496000` | `N=4000, P=1612, P/N=.403, N*P=6448000` | conservative cross-run sparse-performance anchor |
| PO BIC, cluster | `N=8000, P=812, P/N=.1015, N*P=6496000` | `N=1000, P=412, P/N=.412, N*P=412000` | observed sparse-performance anchor |
| XPO plugin (robust and cluster) | `N=1000, P=112, P/N=.112, N*P=112000` | `N=1000, P=412, P/N=.412, N*P=412000` | official default; sparse remains frontier |
| XPO BIC (robust and cluster) | `N=8000, P=812, P/N=.1015, N*P=6496000` | `N=1000, P=412, P/N=.412, N*P=412000` | official default; sparse remains frontier |

An anchor is the earliest tested point whose complete tested suffix had exact
parity, safety, and at least 1.5x speed. The current run directly places
high-density PO-BIC/robust at N=4,000/P=1,612 and low-density PO-BIC at
N=8,000/P=812. The manifest labels the evidence
`advisory_conditional_not_auto_route`: no interpolation, extrapolation, or
automatic routing is justified. Evidence is in
the [routing manifest](../logs/benchmark_routing_manifest_rc13_route_kernelauto_20260719.json)
and its matching raw, resource, pair, guidance, and source CSVs.

## Current rc13 Python preflight and fixed-point gate

The Python backend accepts multiple treatments and posts the full symmetric
cross-treatment covariance in `e(V)` and `e(sparse_covariance)`. It supports PO
robust/one-way/two-way inference and XPO robust or cluster-disjoint K/K^2
folds. Its 18-test suite covers independent robust, one-way, and two-way
score-sum covariance oracles, PSD/full-covariance output, nonzero off-diagonal
covariance, PO/XPO paths, balanced and ragged K^2 geometry, deterministic fold
replay, explicit invalid-fold/empty-cell refusals, and zero-Jacobian rejection.

The standalone core Stata public preflight passes PO and XPO robust, one-way,
and two-way examples in `python.log`, including full covariance, `e(sample)`,
routing/provenance metadata, cluster geometry, factor-treatment, interaction,
always-in refusal, and stale-`e()` cleanup gates. These tests are included in
the passing final hash-bound full tier. The fixed-point state test passes fixed, cycle,
and budget-refusal classifications in `validate_plugin_fixedpoint_state.log`.

Full certification `rc13_cert_full_mit4_20260719` passes with exit 0 and stable
source hashes. Its summary records 96 total cells, 95 numerical matches, 94
production-exact cells, one intentional refusal, one diagnostic frontier, and
maximum production `b`/`V` differences `9.10382880193e-15`/
`3.90312782095e-18`. Quick certification
`rc13_cert_quick_coauthor4_20260719` binds the one-rule Mata computation
default, coauthor quick start, manual-data examples, and final status text to
the unchanged release source.

The source-stable trajectory classification/safety gate
`rc13_fixedpoint_pobicycle_20260718`
passes. The two 10k runs are byte-for-byte repeats: each finds nine cycles, has
zero solver/KKT/verification/exhaustion failures, and takes 24.388 and 25.077
seconds. Fold-target `1|y` has primitive period 6; the other eight have period
2. The 30k run verifies ten main fixed points plus five primitive period-2
cycles with zero failures in 77.589 seconds. Its exact cycle keys are `1|d1`,
`2|y`, `2|d2`, `4|y`, and `5|y`. All seven required source rows, including the
Stata executable, have identical before/after hashes. The runner and runtime
hashes are
`0985F72EBFFD63E24A8642A000D904725196A0F070063159F5809E73E421D370`
and `970F1AFA3673ACB5116D68B52C8FB6196832E30EB41E8FD21A8C785179568FD9`.
The hardened runner requires each child process to exit 0 and emit its final
PASS marker, checks exact fold-target classifications and cycle lengths, and
compares exact hashes for repeated `anchors`, `states`, `summary`, and `support`
artifacts. Evidence is in the
[trajectory runner](../logs/probe_xpo_plugin_fixedpoint_rc13_fixedpoint_pobicycle_20260718_runner.csv),
matching `sources.csv`, and per-run CSVs. The earlier
`rc13_fixedpoint_final2_20260717_1523` and
`rc13_fixedpoint_gate_20260717_1455` runs are superseded engineering evidence.

The standalone public `strict(production)` refusal tests now pass at both
scales. At 10k, production returns `r(430)` with full/training support `3/2`,
solver/KKT/rank `1/1/1`, and plugin pass/capped calls/iterations `0/9/135`. At
30k it returns `r(430)` with support `5/4`, the same `1/1/1` diagnostics, and
`0/5/99` in 58.503 seconds. Evidence is in the
[10k refusal log](../logs/validate_xpo_plugin_cycle_refusal.log) and
[30k refusal log](../logs/validate_xpo_plugin_cycle_refusal_30k.log). The final
full certification binds these logs/tests to its unchanged source snapshot.
The policy does not invent a fixed point: production refuses cycles,
while opt-in `strict(stata)` frontier use may preserve the finite step-15 path
with a warning. A standalone engineering probe shows that the nonmutating cold
KKT rescue leaves support unchanged while reducing relative KKT error from
`5.365761e-4` to `2.096615e-14`; the release decision relies on the trajectory
and public-refusal gates.

## Current hash-bound release-scale matrix

Run `rc13_scale_full_pobicycle_20260718` passes the corrected release launcher:
40/40 serial shards, 80/80 cells, and 40/40 exact-repeat groups. All 35 source
rows are stable. The 16 official-parity cells are exact to numerical precision;
maximum coefficient and covariance differences are `4.663e-15` and
`1.016e-18`. The union contains 76 ordinary passes and four expected
production refusals, with no timeout or unexpected classification.

| current-source profile | repeat result | measured boundary |
|---|---|---|
| PO/plugin, robust, one million controls | 117.394--122.798 s | `N=1,499,994`; truth error `0.01278`; production pass |
| PO/plugin, cluster, one million controls | 118.236--118.857 s | `N=1,499,994`; truth error `0.01277`; production pass |
| PO/BIC, robust, one million controls | 779.426--805.182 s | `N=1,499,994`; production pass |
| PO/BIC, cluster, one million controls | 613.192--617.038 s | `N=1,499,994`; production pass |
| XPO/BIC, one-way cluster-disjoint, 30k controls | 238.910--240.547 s | `N=149,980`; 656 nuisance fits; production pass |
| XPO/BIC, two-way `K^2`, 30k controls | 3,057.789--3,086.429 s | `N=374,950`; 3,872 nuisance fits; feasible, but not a fast default |

Peak RSS is 1.012 GiB for the production family and 0.354 GiB for the native
family. The separate 30k clustered-plugin trajectory diagnostic reaches
26.922 GiB; it is a deliberate refusal/trajectory cell and is excluded from
production recommendations. The authoritative
[release manifest](../logs/benchmark_indexed_scale_release_rc13_scale_full_pobicycle_20260718_manifest.csv)
and [validated union](../logs/benchmark_indexed_scale_release_rc13_scale_full_pobicycle_20260718_union.json)
record the plan/specification hashes and `overall_release_gate=1`.

The earlier unsafe-order 10k/30k shards remain historical only. The corrected
driver restores and asserts `obsid` order, records the fitted
`e(min_train_support)`, and accepts a cycle refusal only when the frozen
trajectory evidence and public diagnostics agree.

Python controls may be plain or `c.` numeric variables; `i.` and `ibn.`
categorical controls are mapped to absorbed indicators. Treatments must be
plain numeric variables. Factor treatments, interactions (`#`/`##`),
parenthesized always-in controls, and nonnumeric inputs are explicit refusals.

## Release-target and frontier contracts

| native profile | status |
|---|---|
| PO plugin robust | release certification target |
| PO BIC robust | release certification target |
| PO BIC cluster | release certification target |
| PO plugin cluster | release certification target; singleton absorbed levels rejected |
| PO plugin/BIC two-way cluster | frontier; explicit `backend(mata) frontierok` |
| all XPO profiles | frontier; explicit `backend(mata) frontierok` |

Official `poregress` / `xporegress` remains the target for safely expandable
small and medium designs. For ordinary-fold native XPO, `strict(stata)` runs
an official retained-column screen inside each training fold and warns on
absent/singleton training levels. `strict(production)` skips that screen and
requires minimum training support of two. `foldcluster()` has no official
screen and remains a separate frontier. PO plugin/cluster still refuses
full-sample singleton absorbed levels but supports balanced and unbalanced
nonsingleton clusters.

## Current-source full empirical XPO parity

The final-source full-sample BIC, one-way-cluster XPO comparisons use five
canonical seeded folds and `rseed(1103)`. All three exact pair gates pass:

| sample | sparse s | official s | speed evidence | selected | max `b` diff | max `V` diff |
|---|---:|---:|---|---:|---:|---:|
| patenting | 264.018 | 2,052.546 | live, 7.77x | 274 | `6.700e-14` | `8.601e-15` |
| no-missing-R&D | 380.778 | 3,862.323 | live, 10.14x | 429 | `1.961e-13` | `5.201e-15` |
| all-firm | 1,137.910 | 15,125.485 | validated anchor; indicative 13.29x | 1,737 | `1.402e-14` | `3.500e-16` |

The all-firm live official cell reached fold 3 but exceeded its two-hour cell
cap. `-UseFrozenOfficialRc12` then validated the completed official-only
artifact against the current Stata executable and official ado hashes, all
input-data hashes, the frozen runner/manifests, and every official artifact.
The pair therefore supports current-source sparse numerical parity, while its
pair timing remains correctly marked non-comparable.

All three manifests report `release_evidence=1`, `overall_run_gate=1`, and
`overall_release_gate=1`. Every support is exact (274, 429, and 1,737), and
the patenting/no-missing-R&D fold and timing comparisons are fully live. All
three designs have minimum training-fold support zero, so these are
`strict(stata)` compatibility claims. They are not `strict(production)`
claims; the expected production response is `r(498)`.

## Retained rc12 empirical XPO evidence

The hash-bound runner uses BIC selection with `gridminok`, one-way cluster VCE,
five observation folds, and `rseed(1103)`; the exact commands are in
`tests/benchmark_empirical_xpo_cell.do`.

| scope/sample | sparse s | official s | speedup | selected | max `b` diff | max `V` diff | support |
|---|---:|---:|---:|---:|---:|---:|---|
| safe/allfirm | 21.145 | 73.567 | 3.48x | 50 | `2.753e-14` | `4.413e-15` | exact |
| safe/patenting | 17.726 | 50.592 | 2.85x | 37 | `1.754e-14` | `5.884e-15` | exact |
| safe/nominssrd | 7.638 | 25.482 | 3.34x | 19 | `1.110e-14` | `9.659e-15` | exact |
| full/patenting | 189.060 | 1,781.973 | 9.43x | 274 | `6.756e-14` | `8.507e-15` | exact |
| full/nominssrd | 248.440 | 3,410.070 | 13.73x | 429 | `1.965e-13` | `5.156e-15` | exact |
| full/allfirm | 683.101 | 15,125.485 | 22.14x | 1,737 | `1.44328993201e-14` | `3.56268833879e-16` | exact |

All artifact, input, and source hashes and all support-set comparisons pass.
The final master and runner each contain 18 clean rows; the artifact index has
75 rows and every indexed hash validates. The sparse support files contain
exactly `k_selected` unique tokens. Official and sparse files order those
tokens differently, so their raw support-file hashes differ even though the
sets are exact. Full/allfirm has zero support-only tokens on either side.
Evidence:

- `logs/benchmark_empirical_xpo_rc12_emp_final_20260717_051000.csv`
- `logs/benchmark_empirical_xpo_runner_rc12_emp_final_20260717_051000.csv`
- `logs/benchmark_empirical_xpo_artifacts_rc12_emp_final_20260717_051000.csv`
- `logs/benchmark_empirical_xpo_rc12_emp_final_20260717_051000.csv.meta.json`

The three full sparse fits pass convergence, finite, final-rank, and KKT gates,
but each has minimum training support zero and uses relaxed `strict(stata)`.
They are compatibility/frontier evidence, not `strict(production)` evidence.
The allfirm sparse fit peaked at 11.736 GiB RSS and took 683.101 seconds; its
official counterpart took 15,125.485 seconds, a 22.14x sparse speedup.

Current-source replacement evidence is produced by
`tests/run_benchmark_empirical_matrix.ps1`. Its pair CSV includes timing
comparability, `b`/`V` and support parity, the canonical seeded-fold and safety
contract, and source/artifact hashes. The three current full-sample manifests
listed above satisfy their complete pair, stable-source, and release gates.
`tests/validate_empirical_matrix_harness.ps1` statically checks this schema and
the canonical-order anchors, then exercises the real `-PlanOnly` path without
launching Stata. See `EMPIRICAL_RELEASE_MATRIX.md`. The broader 90-cell
baseline and 144-cell expanded empirical matrices remain separate optional
campaigns; the three requested full XPO-BIC/cluster pairs are complete.

The non-estimation baseline plan contains 90 cells and 30 sparse/official
pairs. It covers every plugin/BIC by robust/cluster profile on all three safe
samples plus BIC/cluster on each full sample. The ultimate expanded full-scope
plan adds `-AllowFullMatrix` and contains 144 cells and 48 pairs; it remains
separate pending work and must not be inferred from a passing 90-cell manifest.

## Retained rc12 inference diagnostic

`coverage_inference_rc12_full_20260717_043000.csv` runs 500 repetitions for
each of 12 profiles on one `N=480` heteroskedastic DGP:

| command/selection | robust | cluster1 | cluster2 |
|---|---:|---:|---:|
| PO plugin | 500/500, 0.918 | 500/500, 0.902 | 500/500, 0.898 |
| PO BIC | 500/500, 0.920 | 500/500, 0.234 | 500/500, 0.906 |
| XPO plugin | 500/500, 0.940 | 500/500, 0.874 | 499/500, 0.996 |
| XPO BIC | 500/500, 0.926 | 500/500, 0.480 | 500/500, 0.998 |

The entries are successful fits and empirical 95% coverage. They diagnose this
specific DGP; they are not nominal-coverage certification. The only refusal is
the intentional XPO/plugin/two-way `r(430)` for an unstable capped loading
update. Targeted sparse and official PO/XPO BIC cluster1 fits are numerically
exact and share the low coverage, so that result is not a sparse-parity defect.
See the matching runner, manifest, and artifact CSVs plus
`probe_coverage_bic_cluster1_summary.csv` and
`probe_coverage_xpo_plugin_cluster2_failure.log` under `logs/`.

At 10k features, the rc12 PO plugin/cluster probe remains historical evidence.
The old XPO discrepancy is superseded by the corrected safe-order pilot, which
is exact for both robust and clustered `b,V`. Clustered XPO remains a
loading-cycle safety frontier: production refusal is required even though the
reported estimate and covariance match official Stata.

## Pre-freeze rc11 evidence - historical/superseded

These artifacts validate the implementation direction but were created before
the planned rc11 version/hash freeze. They must not be cited as current rc13
release evidence.

### Ordinary-fold empirical XPO

`benchmark_empirical_xpo_20260716_230600.csv` pairs official and sparse BIC,
one-way-cluster XPO on three support-safe empirical samples:

| sample | official seconds | sparse seconds | max `b` diff | max `V` diff | support |
|---|---:|---:|---:|---:|---|
| allfirm | 78.727 | 21.745 | `2.753e-14` | `4.413e-15` | exact |
| patenting | 52.057 | 17.921 | `1.754e-14` | `5.884e-15` | exact |
| nominssrd | 25.240 | 7.711 | `1.110e-14` | `9.659e-15` | exact |

`benchmark_empirical_xpo_20260717_004000.csv` records the three full sparse
fits before the freeze. All pass convergence, finite, final-rank, and KKT
gates, but all use relaxed training support and have no paired full official
result yet:

| sample | N | absorbed groups | seconds | selected | min train | xcheck seconds |
|---|---:|---:|---:|---:|---:|---:|
| allfirm | 86,341 | 11,544 | 658.024 | 1,737 | 0 | 198 |
| patenting | 45,913 | 4,312 | 184.510 | 274 | 0 | 38 |
| nominssrd | 49,217 | 7,129 | 245.367 | 429 | 0 | 67 |

### Scale and inference diagnostics

- `benchmark_indexed_scale_20260717_001402.csv` is retained only as pre-freeze
  harness lineage. Its approximate XPO plugin differences came from the unsafe
  second fold permutation and are superseded by the corrected safe-order exact
  results.
- `benchmark_indexed_scale_20260717_000825.csv` is the pre-freeze native
  cluster run. All six one-/two-way PO/XPO cells pass truth, convergence,
  finite, final-rank, KKT, and production plugin gates.
- `coverage_inference_20260717_010202.csv` is a 50-rep diagnostic pilot, not a
  certification study. Its 95% coverage rates are:

  | command/selection | robust | cluster1 | cluster2 |
  |---|---:|---:|---:|
  | PO plugin | 0.88 | 0.94 | 0.92 |
  | PO BIC | 0.90 | 0.24 | 0.94 |
  | XPO plugin | 0.92 | 0.94 | 1.00 (49/50 successful) |
  | XPO BIC | 0.90 | 0.42 | 1.00 |

  The small pilot is useful for routing diagnostics, not for a release-level
  coverage claim.
- `probe_coverage_bic_cluster1_summary.csv` shows sparse and official PO/XPO
  BIC cluster1 estimates, covariance, and support agree to numerical precision;
  both share the pilot's poor coverage. This is not a sparse-parity defect.
- `probe_coverage_xpo_plugin_cluster2_failure.log` sequentially isolates rep
  39 and confirms the sole plugin/two-way pilot failure is the intentional
  production `r(430)` gate: one capped loading update neither converges nor
  stabilizes support.

## Historical rc8 verification snapshot

The 2026-07-16 quick tier passes after the PO multiway implementation. It
includes every prior official-parity and cluster-DML gate, the independent
two-treatment inclusion-exclusion oracle, plugin/BIC public paths, reversed and
equivalent cluster partitions, ragged intersections, PSD repair, interface
refusals, and the 30,002-feature benchmark.

| rc8 artifact | SHA-256 |
|---|---|
| `stata/_sparse_prl_xpo_mata.ado` | `444E328425B4C6DC718F3FB209DE4D96CA092A90829595C7B2EED772D90D182A` |
| `stata/sparse_poreg.ado` | `CE9BA2EC1E26296D5DED63DC8288CB798CBDDB81E3C71858D183BDE76D4A3F88` |
| `stata/sparse_poreg.sthlp` | `E6719CAF13423B9CAFC8350C63650FB994B04086FF8404C2CBA092FC9D1E1241` |
| `tests/validate_po_multiway_vce.do` | `31DBB7AB4E250A8DF1B3D6D987F840F5CB8EF13BD54DAF96EB8EABF885EE03E3` |
| `tests/benchmark_indexed_large.do` | `33854A302A29C54B636474B920E94B146283A96F4971B7C3233679E7735FA41B` |
| `tests/run_certification.do` | `375E25091B327D822B7C7F141CE8A29FE9050D3EC48B6A2DB6B1D52FCBF2294F` |
| `docs/PO_MULTIWAY_20260716.md` | `FE8A5B1C9140F922E817DB81ABB7537C4859B3D6305A277848D1331519B776AB` |
| `logs/validate_po_multiway_vce_16_Jul_2026_121329.log` | `A86DEA01DDBF855BB9B4CB4FD35BBF9501123F51AE1ECBF4DB3496D20523D2E1` |
| `logs/benchmark_indexed_large_16_Jul_2026_121450.log` | `62BD8D4018CC9A065D2CFB683880269F38EE5FB71281FEC7A6EA4725BDE5E73E` |
| `run_certification.log` | `82BA4E974555BF2428B4FCD35B043C4954BD466462AAC16D7A16C4A6EFC2B6C1` |

## Frozen rc7 verification snapshot

The 2026-07-16 quick tier passes after the cluster-disjoint implementation. It
includes the unchanged official-parity suite, independent one- and two-way
cluster-DML covariance oracles, public BIC and strict-plugin tests, isolated
support failures, and the 30,002-feature indexed-kernel benchmark. A later
focused rerun also passes the final row-order and cluster-order checks.

| rc7 artifact | SHA-256 |
|---|---|
| `stata/_sparse_prl_xpo_mata.ado` | `7FFEB23B10936E81C52147D7D96720D4D702E96C7A426D8F2E80D0086232A55F` |
| `stata/_sparse_cluster_folds.ado` | `1DFBB6F92EF0B7AAD6646BA399754CC902AE448382A3E50D4C5AC727C27564EE` |
| `stata/_sparse_native_guard.ado` | `A88273749CF6CB0BD383A91C3E358A92DEFDEFE7C22F00F1F92F3859BE9EFE7E` |
| `stata/sparse_xporeg.ado` | `D507250A62EB07FC2AAED214A5B5A553C12D090A75E3754B94B5839C90772098` |
| `stata/sparse_xporeg.sthlp` | `4BEF5FEAB842A37174779BD435ACE629C0706B6417A5023E668FF33421E28EDC` |
| `tests/validate_cluster_dml_vce.do` | `87AC99C9F580C7515A2C86BFE73337BEF1EDB87A1C173A7B9B71493C1E123939` |
| `tests/validate_cluster_disjoint_xpo.do` | `8067B5140B7910F1F3B2A0054A2FD00DBCBBAB44E8FE451BB9FC22137835B2D5` |
| `tests/validate_xpo_support_guard.do` | `183852D4A81A6E4D91A19D7F01704908FC0B5DC12E61511EFCEF0062A55A3F4E` |
| `tests/run_certification.do` | `18A0A87CD1AFA10EAEE57F8D25856ECB391E3FD09345544D88F542756EC42C50` |
| `run_certification.log` | `F5CB35239DFED27200F5DE530F11BC587430511BE7433D1C43CAB669ACD81123` |

## rc6 parity evidence

The deterministic continuous-column oracle recovers all seven official
loadings with maximum relative error `1.22e-14`. Unbalanced newA, newB, and
newD supported PO/XPO plugin-cluster comparisons are at machine precision,
and indexed/legacy kernel parity passes. See
`CLUSTER_PLUGIN_FIX_20260715.md` for the exact rule, logs, and remaining XPO
boundary.

The current K=3/5/7/10 plugin/BIC by robust/cluster matrix matches official
Stata to `1.776e-15` in coefficients and `8.674e-18` in covariance. Breathe
XPO/BIC/robust matches to `7.550e-15` and `1.388e-16`; all 20 fold selections
match. Key final-source logs are
`validate_indexed_full_16_Jul_2026_000913.log`,
`validate_xfold_vce_16_Jul_2026_000849.log`,
`breathe_xpo_bic_robust_16_Jul_2026_000412.log`, and
`benchmark_indexed_large_16_Jul_2026_001007.log`.

| frozen rc6 candidate artifact | SHA-256 |
|---|---|
| `stata/_sparse_prl_xpo_mata.ado` | `C307877B09D61827A2D2360C84EA2943AAED0DD46DB4FE335709E424BFDFBED9` |
| `stata/_sparse_plugin_probe.ado` | `A31327ABEF3B5F3F720A9F31CB0222A9C4BF9509A243CA2C39F8D110D79131B1` |
| `stata/sparse_poreg.ado` | `B4BE3FE4C7AD85586C428E57E6F82EAE1E3EE9AB9E36BEA0F8657DE28A5C2C28` |
| `stata/sparse_xporeg.ado` | `44D31E60F9EDF54C2C17944CE71A63C68C2825E0449BBBCD7DB22ADFF0333DA8` |
| `stata/sparse_poreg.sthlp` | `03804C5AF3563A5ADE57398EA03E76448D554A57EE3280C6B2B70241804AC127` |
| `stata/sparse_xporeg.sthlp` | `DB915C56ED22E82F2DF9CC8D610F104EC25CDDE06CE44888F63A7F1F1EB0F8A9` |
| `tests/validate_indexed_full.do` | `5566474CEFC79775182752A38D390AA8EF7995F5914B952163E060D2F2803E54` |
| `tests/validate_new5.do` | `3A19BC2DF35F8D984BA93F1214C16A91DD3850443330D840619A9DA9AA634EA8` |
| `run_certification.log` | `FD6C5D16886E83F2A62139F803A76F97C88AFB625C89D85589CB17B82E3FF424` |

## Frozen rc4 source snapshot (historical)

| artifact | SHA-256 |
|---|---|
| `stata/_sparse_prl_xpo_mata.ado` | `F44F4A7D53766F0452A58FCB57A17ADD87B8379B091E0A9BF5178CEF087DFBF4` |
| `stata/sparse_poreg.ado` | `B8BD76172D60C630FD64DACA34E47F0C0E4089A91ED61BD5A79242A5F6C1C03D` |
| `stata/sparse_xporeg.ado` | `DC96A19E27B3C2A8BEEA7B9737FF0565FA1F9EBEE0C8B14E70D4FDD2EC4E5AB2` |
| `stata/_sparse_official_route.ado` | `429E38DACC5DD94040EA29D77F9269DBCB0B96EEBF944E289FE3222BE545CE1B` |
| `stata/_sparse_native_guard.ado` | `7C356AC9627B06E9B12220920B23CAAB3B5C8403171221FA5FE93275E84C3B75` |
| `python/sparse_prl_xpo_backend.py` | `41ECA045E1A5501A6E75A2A6318FFE919A2BF0962BD0047E960D2550DFD5ED63` |

## Frozen rc2 baseline (historical)

- Parent Git base: `1a4c73366796980bb226f455ee8d1ffb9a4e8386`
- Candidate version: `0.8.4-rc2`
- Stata: 19 MP for Windows

| artifact | SHA-256 |
|---|---|
| `stata/_sparse_prl_xpo_mata.ado` | `1ED702DE858B8FFE355CF15AA3B8BFE6C1BE92241391F6FE0620AD57542A6552` |
| `stata/sparse_poreg.ado` | `DBEA691D2DA8623F13D08824ED8B4990E22B404AFD0E4F7B1F56F0FD077B7F4F` |
| `stata/sparse_xporeg.ado` | `0891FDDCE8650914D28CC4B3D172E68121257450629BC1CA5619AF02F104CEEB` |
| `python/sparse_prl_xpo_backend.py` | `3DA6CAD1CD215B97E0CC82B36AA6BB2A670BB3C3F8F65D3F1CE3B0D848A19436` |
| `logs/validate_indexed_kernel_15_Jul_2026_120224.log` | `CCE1E91A0507AC75F6BC7E781D97798D78D152606F51EB4D673F32C0EB955CED` |
| `logs/benchmark_indexed_large_15_Jul_2026_120300.log` | `FB85745E8737F08619B6C38A880F213A451F0AFA30CAB59620A9F26B1869ABFC` |

The deleted transfer manifest was checked before deletion: all 19 canonical
artifacts matched their recorded hashes.

## Historical applied verification evidence

The full historical-data rerun completed on 2026-07-16.  It records 39
machine-precision matches and 17 intended singleton-support guards across the
56-cell N grid, plus fresh newA-newE and breathe results.  See
`LEGACY_VALIDATION_20260716.md` for the exact classification and log names.

- The final quick tier passes in `run_certification.log`; its principal
  component logs are `validate_official_interface_16_Jul_2026_121326.log`,
  `validate_diagnostics_16_Jul_2026_121328.log`,
  `validate_po_multiway_vce_16_Jul_2026_121329.log`,
  `validate_xfold_vce_16_Jul_2026_121331.log`,
  `validate_indexed_kernel_16_Jul_2026_121343.log`,
  `validate_indexed_full_16_Jul_2026_121356.log`, and
  `benchmark_indexed_large_16_Jul_2026_121450.log`.
- `validate_rc4_focus.log` passes exact full/subset, robust/cluster, lasso,
  path, BIC, block-Gram, sign-drop, and guarded-fallback comparisons with
  maximum difference zero.
- Duplicate-partition compaction reduces the adversarial expanded design from
  71 raw to 53 effective controls, matching Stata's selection dimension.
- Centered coordinate descent closes the prior safe-support and breathe XPO
  selection gaps at machine precision. Native XPO still requires frontier
  opt-in because broad coverage remains incomplete, `K^2` BIC is expensive,
  and plugin/cluster can fail its KKT certificate.
- Final newB-newE grids complete 32 supported cells with maximum coefficient
  and covariance differences `4.885e-15` and `2.277e-18`. NewA completes its
  four PO comparisons and records the required native singleton-fold refusal
  for all four XPO profiles.
- `benchmark_indexed_scale_20260715_203958.csv` and its RSS companion pass the
  frozen 10k plugin/robust and BIC/cluster cells. The plugin official oracle
  differs by `1.665e-15` in coefficients and `8.403e-19` in covariance.
- The 30,002-feature distinct-support smoke passes with finite symmetric
  covariance and convergence/KKT/rank diagnostics.

## Frozen clustered-BIC scale evidence

| target | estimator | wall | peak RSS | truth error |
|---|---:|---:|---:|---:|
| 500k | 275.279 s | 276.664 s | 623,407,104 bytes (594.5 MiB) | 0.012562 |
| 1m | 579.990 s | 582.568 s | 993,562,624 bytes (947.5 MiB) | 0.012772 |

The rows are in `benchmark_indexed_scale_20260715_204518.csv` and
`benchmark_indexed_scale_20260715_205029.csv`, with their RSS companions.
Both have `rc=0`, `gate_rc=0`, exact feature counts, indexed kernel,
convergence, finite, rank, covariance, KKT, and truth `<0.05` gates.

- 500k: 107 full and 208 active sweeps, 104 inactive KKT scans, zero
  fallbacks, active maximum 21, polish 104/104, relative KKT `1.0011e-12`.
- 1m: 113 full and 220 active sweeps, 110 inactive KKT scans, zero fallbacks,
  active maximum 36, polish 110/110, relative KKT `5.1214e-13`.

The 500k `b`, `V`, and counters exactly match the prior two-run evidence in
`benchmark_indexed_scale_20260715_193332.csv`. The 1m row exactly matches the
two prior runs in `benchmark_indexed_scale_20260715_192323.csv` and
`benchmark_indexed_scale_20260715_194249.csv`.

Earlier-source plugin/robust and BIC/robust runs reached 1m, but they are not
substitutes for the remaining frozen rc4 repeated ladder.

## Large empirical PO BIC parity (rc10, pre-freeze)

`strict(stata)` now harvests Stata's retained design once and applies the same
mask to every native PO-BIC nuisance path. On the sorted 86,341-observation,
11,570-candidate all-firms replication design, all 14 official nuisance
subspaces expose the same 11,263 retained columns; 307 firm indicators are
dropped and all 26 year indicators remain.

The full native result is `b(RDAT)=0.198636661`,
`SE(RDAT)=0.018130383`, and union support 1,241. Against the saved Stata 19
oracle, maximum absolute differences are `1.352e-13` for `b` and `2.089e-15`
for `V`. Total sparse compatibility runtime was 224.947 seconds versus
7,124.101 seconds for the previously saved official run. The raw local data,
after the exact rename/clean and `sort PERMCO fyear`, independently returned
the same coefficient, SE, mask counts, and union support in 215.489 seconds.

Evidence:

- `logs/benchmark_github_sparse_po_xcheck.log`
- `logs/benchmark_allfirm_raw_sorted_sparse_po_xcheck.log`
- `tests/verify_esrf_common_mask.do`
- `logs/validate_one_N1000_16Jul2026_203106.csv` (all eight PO/XPO cells pass;
  maximum `b` and `V` differences remain at numerical precision)

This is an input-order compatibility contract, not a claim that Stata's hidden
screen is permutation-invariant. To reproduce the historical 0.1986367 fit,
sort exactly as the GitHub runner before both estimators. `strict(production)`
continues to use the full native design without the official screening pass.

## PO multiway boundary

rc8 adds explicit `vce(cluster g1 g2)` behind `backend(mata) frontierok`.
The final VCE is the sum of both marginal cluster sandwiches minus the observed
intersection sandwich, with per-component native one-way corrections,
Stata-style PSD repair, and `min(G1,G2)-1` degrees of freedom. Plugin and BIC
nuisance selection are observation-weighted, so this is not an official
`poregress` or multiway-plugin parity claim. The independent full-V, public
order, equivalent-partition, ragged-intersection, PSD, and refusal gates pass.
See `PO_MULTIWAY_20260716.md`.

## XPO clustering boundary

The official-compatible default still uses observation-level outer folds under
one-way clustered VCE. For rc11 `strict(stata)`, each ordinary training fold
gets its own official retained-column screen; the resulting fold-local mask is
used only for that fold. This path posts
`e(sparse_xcheck)="stata_fold"`,
`e(kkt_scope)="stata_fold_screened_design"`, screening time, and foldwise
kept/dropped ranges.

Historical rc7 work added explicit `foldcluster()`: seeded one-way
cluster partitions and true two-way `K^2` train/test cells with Chiang--DoubleML
score covariance. The public BIC and strict-plugin paths, independent full-V
formula oracle, row/order invariance, valid/invalid support fixtures, and the
current 30k release-scale cells pass. See `CLUSTER_DML_XPO_20260716.md` and
`validate_cluster_disjoint_xpo_16_Jul_2026_114027.log`. This mode is not an
official-Stata parity claim and has no fold-local official screen. It remains a
frontier because broad nominal-coverage evidence is incomplete and the 30k
two-way BIC cell takes 3,058--3,086 seconds. `strict(production)` is preferred.

## Reliability controls now present

- max-iteration exhaustion fails before `ereturn post`;
- nonfinite output and rank-deficient final treatment regressions fail;
- final-solution KKT, plugin convergence, selection, rank, condition, and
  finite diagnostics are posted;
- plugin support stability, capped loading-call count, maximum loading change,
  and combined pass status are posted;
- XPO minimum full/training support, safe/relaxed mode, official-screen scope,
  foldwise kept/dropped ranges, and screen time are posted;
- coordinate-descent full/active sweeps, inactive KKT scans/fallbacks, maximum
  active size, and block-Gram polish calls/successes are posted;
- nuisance rank deficiency is reported separately;
- PO multiway cluster counts, raw VCE minimum eigenvalue, PSD-repair status,
  covariance rank, and denominator degrees of freedom are posted;
- XPO `strict(production)` and PO plugin/cluster guards prevent known unsafe
  singleton-support structures; XPO `strict(stata)` records relaxed handling.

PO profiles have a hard KKT gate. Every native `selection(plugin)
strict(production)` fit also rejects a loading path that neither converges nor
passes its cold support/KKT certificate. XPO additionally enforces its
training-support guard. Explicit `strict(stata)` research fits may warn and
post failed plugin diagnostics.

The scale driver isolates every cell in a fresh process, validates row
identity, supports per-cell timeouts, measures RSS, checks exact repeats, and
rejects truth error above `0.05`.

## Retained release boundaries

1. Full certification `rc13_cert_full_mit4_20260719` passed on its recorded
   source snapshot. Current routing `rc13_route_kernelauto_20260719` and quick
   certification `rc13_cert_quick_coauthor4_20260719` bind the automatic
   indexed default; the full suite already exercised that engine explicitly.
2. Native XPO and PO multiway inference remains documented as a frontier. The completed 500-repetition
   diagnostic does not justify a blanket nominal-coverage claim.
3. Keep Python experimental: it is LassoCV rather than plugin/BIC, inner CV is
   row-wise, its two-way DoubleML positive-sum VCE is not the CGM/Stata
   minus-intersection sandwich or small-cluster correction, it does not post
   exact selected absorbed IDs, and LassoCV warnings are not a production
   convergence gate.
4. Keep synthetic routing anchors advisory; they are not universal N/P cutoffs.
5. The release builder uses an explicit allowlist, writes `SHA256SUMS`, and
   independently verifies every archived payload hash.
6. MIT is selected in the root `LICENSE` and must remain in every archive.

The release builder's explicit allowlist excludes the three empirical `.dta`
files under `tests/`; they remain private test evidence rather than public
payload.
