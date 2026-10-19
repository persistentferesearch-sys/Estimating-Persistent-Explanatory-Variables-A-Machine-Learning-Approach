{smcl}
{* *! sparse_poreg 0.8.4-rc14}
{vieweralsosee "sparse_xporeg" "help sparse_xporeg"}
{vieweralsosee "poregress" "help poregress"}

{title:Title}
{pstd}{bf:sparse_poreg} -- Linear partialing-out lasso with official and sparse Mata routes{p_end}

{title:Quick start}
{pstd}For an artificial firm-year example, type:{p_end}
{phang2}{cmd:. findfile finance_quickstart.do}{p_end}
{phang2}{cmd:. do "`r(fn)'"}{p_end}
{pstd}The example compares the coefficient and standard error from the official and Mata routes.{p_end}

{title:Syntax}
{p 8 16 2}{cmd:sparse_poreg} {it:outcome targets} {ifin}{cmd:,}
[{cmd:absorb(}{it:numeric_ids}{cmd:)}
{cmd:controls(}{it:candidate_controls}{cmd:)}
{cmd:selection(plugin|bic[, gridminok]|cv|adaptive)}
{cmd:vce(robust|cluster id [id2])}
{cmd:backend(auto|official|mata)}
{cmd:strict(stata|production)} {cmd:frontierok}
{cmd:rseed(#)} {cmd:densemax(#)}]{p_end}

{title:What the variables mean}
{pstd}The first variable is the outcome. Every variable after it and before the comma is a target with a reported coefficient. Variables in {cmd:controls()} are candidates for selection. Each observed level of an ID in {cmd:absorb()} is also a candidate indicator. This is not conventional within-transformation fixed-effect absorption.{p_end}

{title:Choosing a route}
{phang}{cmd:backend(auto)} runs official {helpb poregress} when the expanded design fits; otherwise it uses the supported sparse Mata route. This is the default.{p_end}
{phang}{cmd:backend(official)} calls official {helpb poregress}. It accepts official selection and postestimation options, subject to Stata's own limits.{p_end}
{phang}{cmd:backend(mata)} uses sparse indicators. It supports plugin and BIC selection with robust or one-way clustered inference. Two-way clustering requires {cmd:frontierok} and is a separate package extension.{p_end}
{pstd}{cmd:strict(stata)} is the Mata comparison mode. Use the same sample, observation order, indicator design, selection rule, and VCE when comparing with official Stata. {cmd:strict(production)} avoids the dense official screen and applies native safety checks; it is not an exact-parity promise.{p_end}

{title:Options most readers need}
{phang}{cmd:selection(plugin)} is the default. {cmd:selection(bic, gridminok)} uses BIC. CV and adaptive selection are available through the official route only.{p_end}
{phang}{cmd:vce(robust)} is the default. {cmd:vce(cluster firm)} clusters by firm. Two-way clustering is available only with explicit Mata frontier syntax.{p_end}
{phang}{cmd:densemax(#)} sets the maximum official expansion considered by routing; the default is 20,000. Stata's {cmd:maxvar} limit still applies.{p_end}
{phang}{cmd:rseed(#)} forwards a seed where relevant. For numerical comparisons, sort the data and use the same seed on both routes.{p_end}

{title:Empirical finance example}
{phang2}{cmd:. sort PERMCO fyear}{p_end}
{phang2}{cmd:. sparse_poreg lnnpatent RDAT, absorb(PERMCO fyear) ///}{p_end}
{phang2}{cmd:>     controls(lnME lnAge TobinQ ROA Leverage) ///}{p_end}
{phang2}{cmd:>     selection(bic, gridminok) vce(cluster PERMCO) ///}{p_end}
{phang2}{cmd:>     backend(mata) strict(stata)}{p_end}
{pstd}This reports a coefficient for RDAT. The financial variables and firm/year indicators are candidate controls. See the package README for the complete public-data replication with multiple targets.{p_end}

{title:After estimation}
{phang2}{cmd:. estat selection}{p_end}
{phang2}{cmd:. estat selection, support}{p_end}
{phang2}{cmd:. matrix list e(b)}{p_end}
{phang2}{cmd:. matrix list e(V)}{p_end}
{pstd}Mata results post the selected-control count and numerical checks. An uncertified fit can stop with an error rather than report estimates. Official results retain official postestimation.{p_end}
