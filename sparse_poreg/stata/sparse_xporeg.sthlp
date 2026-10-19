{smcl}
{* *! sparse_xporeg 0.8.4-rc14}
{vieweralsosee "sparse_poreg" "help sparse_poreg"}
{vieweralsosee "xporegress" "help xporegress"}

{title:Title}
{pstd}{bf:sparse_xporeg} -- Cross-fitted linear partialing-out lasso with official and sparse Mata routes{p_end}

{title:Syntax}
{p 8 16 2}{cmd:sparse_xporeg} {it:outcome targets} {ifin}{cmd:,}
[{cmd:absorb(}{it:numeric_ids}{cmd:)}
{cmd:controls(}{it:candidate_controls}{cmd:)}
{cmd:selection(plugin|bic[, gridminok]|cv|adaptive)}
{cmd:vce(robust|cluster id [id2])}
{cmd:xfolds(#)} {cmd:rseed(#)} {cmd:foldcluster(ids)}
{cmd:backend(auto|official|mata)}
{cmd:strict(stata|production)} {cmd:frontierok}]{p_end}

{title:Description}
{pstd}This command estimates a cross-fitted linear model. It learns nuisance models outside each held-out fold and combines the holdout scores. Every variable after the outcome and before the comma is a reported target. Variables in {cmd:controls()} and observed levels of IDs in {cmd:absorb()} are candidate controls, not conventional absorbed fixed effects.{p_end}

{title:Choosing a route}
{phang}{cmd:backend(auto)} uses official {helpb xporegress} when the expanded design fits. If it does not fit, auto stops; it never silently switches to Mata XPO.{p_end}
{phang}{cmd:backend(official)} calls official {helpb xporegress}. Official CV, adaptive selection, resampling, and DML technique options are available here.{p_end}
{phang}{cmd:backend(mata) frontierok} deliberately chooses the experimental sparse route. Mata supports plugin and BIC selection. {cmd:frontierok} does not disable safety checks.{p_end}
{pstd}{cmd:strict(stata)} uses fold-local Stata-compatible screens for comparisons. {cmd:strict(production)} uses native checks and requires sufficient training support. Matching one tested example does not guarantee parity on another sample.{p_end}

{title:Folds and inference}
{pstd}{cmd:xfolds(#)} sets the number of outer folds; default 10. {cmd:rseed(#)} makes the split reproducible; default 1103. Without {cmd:foldcluster()}, official and Mata XPO split observations, even with clustered standard errors. Therefore observations from the same firm can appear in both training and holdout samples.{p_end}
{pstd}To keep entire clusters together, use matching {cmd:vce(cluster firm)} and {cmd:foldcluster(firm)} with explicit Mata syntax. This is a different design from the official-parity example. Two cluster variables create K-squared score cells and can be expensive.{p_end}

{title:Empirical finance example}
{phang2}{cmd:. sort PERMCO fyear}{p_end}
{phang2}{cmd:. sparse_xporeg lnnpatent RDAT, absorb(PERMCO fyear) ///}{p_end}
{phang2}{cmd:>     controls(lnME lnAge TobinQ ROA Leverage) ///}{p_end}
{phang2}{cmd:>     selection(bic, gridminok) vce(cluster PERMCO) ///}{p_end}
{phang2}{cmd:>     xfolds(5) rseed(1103) backend(mata) frontierok strict(stata)}{p_end}
{pstd}This reports a coefficient for RDAT. See the package README and RDPatentComparison.md for the published public-data specification, which has multiple targets.{p_end}

{title:After estimation}
{phang2}{cmd:. estat selection}{p_end}
{phang2}{cmd:. matrix list e(b)}{p_end}
{phang2}{cmd:. matrix list e(V)}{p_end}
{pstd}Mata XPO can refuse an unsafe fold or an uncertified numerical result rather than report estimates. Official results retain official postestimation.{p_end}
