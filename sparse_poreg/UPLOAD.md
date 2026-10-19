# Upload using your other GitHub account

This prepared sharing edition contains no `.dta`, local `.ster`, or unredacted local run log. The paper data stay in the [public replication repository](https://github.com/persistentferesearch-sys/persistent_fe_remendy).

1. Sign in to the GitHub account you want to use and create a repository for `sparse_poreg`. To let visitors browse and install the source, upload the **contents of the extracted package folder** so `README.md`, `sparse_poreg.pkg`, `stata/`, `docs/`, and `results/` are at the repository root. In GitHub's web UI choose **Add file -> Upload files**.
2. Attach the prepared ZIP and adjacent `.zip.sha256` to a GitHub release: **Releases -> Draft a new release**, choose a tag such as `v0.8.4-rc14-sharing`, then attach both files. Mark it as a pre-release if you want the experimental Mata XPO boundary prominent. A ZIP uploaded only as one repository file will not display its contents as browsable source.
3. Review the rendered `README.md`, `COMPARISON.md`, and `results/` under that account before publishing. The package has no automatic posting action.

To check the archive before upload, from the extracted folder run:

```powershell
powershell -File scripts/verify_release.ps1 -Archive "C:/path/to/sparse_poreg-0.8.4-rc14-share-20261019.zip"
```

To install from a downloaded ZIP, extract it and run in Stata:

```stata
cd "C:/path/to/extracted/sparse_poreg-0.8.4-rc14-share-20261019"
net install sparse_poreg, from("`c(pwd)'") replace
do "tests/finance_quickstart.do"
```

The result logs are evidence, not installation inputs. Download the public `.dta` separately for the paper example.
