# Estimating Persistent Explanatory Variables: A Machine Learning Approach

This repository contains the public data, Stata code, and saved results for *Estimating Persistent Explanatory Variables: A Machine Learning Approach*.

The paper's replication code covers four approaches:

1. Linear and Poisson regressions with and without firm fixed effects
2. Within, between, and overall R-squared
3. The adjusted Hausman–Taylor approach
4. Post-regularization LASSO (PRL) and double machine learning (DML)

The [sparse_poreg folder](sparse_poreg/) is a separate Stata package for **linear** PRL and cross-fitted DML. Its [step-by-step guide](sparse_poreg/README.md) starts with artificial finance data and then shows the public R&D–patent example.

## Data and code

- [Replicate_code_MissingRDPatent_ML.do](Replicate_code_MissingRDPatent_ML.do) is the paper's replication script.
- [Replicate_data_MissingRDPatent_ML.dta](Replicate_data_MissingRDPatent_ML.dta) is the public replication dataset.
- [panelr2.ado](panelr2.ado) calculates overall, between, and within R-squared.
- [ReplicateLog_Results](ReplicateLog_Results/) contains saved logs and results.
- [sparse_poreg](sparse_poreg/) contains the separate linear Stata package.

## Reproduce the paper's results

1. Select **Code → Download ZIP** above and extract the repository.
2. Create a folder named `my_replication_output` inside the extracted repository. This keeps newly generated results separate from the saved results.
3. Open `Replicate_code_MissingRDPatent_ML.do` in Stata's Do-file Editor. Near the top, replace its computer-specific `cd` command and data-loading command with the lines below. Replace the example path with the actual location of your extracted repository.

   ```stata
   cd "C:/path/to/repository/my_replication_output"
   clear all
   adopath ++ ".."
   set maxvar 120000
   use "../Replicate_data_MissingRDPatent_ML.dta", clear
   ```

4. Install the commands used by the script if they are missing: `reghdfe`, `ppmlhdfe`, and `estout` (which provides `estadd` and `esttab`).
5. Run the edited do-file in Stata. It writes new logs and CSV results to `my_replication_output`.

The supplied replication script specifies Stata 17 or later. Some models can take substantial time; the saved outputs in [ReplicateLog_Results](ReplicateLog_Results/) let you inspect the paper's results without rerunning every model.

## Sample data construction

`PERMCO` is a masked firm identifier. The following variable descriptions summarize the public data; see the paper and original data sources for further details.
  - `patent`: Number of patents (USPTO PatentsView)
  - `lnnpatent`: Logarithm of one plus the number of patents (USPTO PatentsView)
  - `RDAT`: R&D expenses divided by total assets (Compustat)
  - `lnME`: Logarithm of market equity (CRSP)
  - `RD_missing`: Indicator that R&D expenses are missing (Compustat)
  - `lnAge`: Logarithm of years a firm exists in Compustat
  - `lnK2L`: Logarithm of net property, plant, and equipment divided by employees (Compustat)
  - `TobinQ`: Firm market value divided by replacement cost (Compustat and CRSP)
  - `ROA`: Return on assets (Compustat)
  - `Leverage`: Long-term plus current debt divided by total assets (Compustat)
  - `CASHAT`: Cash divided by total assets (Compustat)
  - `KZidx`: Financial constraint index (Compustat)
  - `InstOwn`: Institutional ownership divided by shares outstanding (CRSP and Thomson/Refinitiv)
  - `oms_HHidx`: One minus the Herfindahl–Hirschman index based on three-digit industry classification (Compustat)
  - `oms_HHidx_square`: Square of `oms_HHidx` (Compustat)


## Saved results

The [ReplicateLog_Results folder](ReplicateLog_Results/) contains the supplied Stata logs and CSV estimates for the linear, Poisson, PRL, and DML analyses. GitHub displays several uploaded files under abbreviated names. The replication script shows the descriptive output names it generates when run locally.

## Try the linear sparse package

You can try the package without running the full paper replication. Follow the [sparse_poreg installation guide](sparse_poreg/README.md), which includes a small artificial-data check and the public R&D–patent example. The detailed comparison with Stata's official linear commands is in [RDPatentComparison.md](sparse_poreg/RDPatentComparison.md).

## Software

The supplied paper replication script specifies Stata 17 or later and uses Stata's official `poregress`, `xporegress`, `popoisson`, and `xpopoisson` commands. It also calls `reghdfe`, `ppmlhdfe`, `estadd`, and `esttab`; install the corresponding community packages if they are not available. The relative `adopath` line above makes the included `panelr2.ado` available to that script.

The separate [sparse_poreg package](sparse_poreg/) documents its own installation and requirements.

## Citation

Please cite *Estimating Persistent Explanatory Variables: A Machine Learning Approach* when using these data.

