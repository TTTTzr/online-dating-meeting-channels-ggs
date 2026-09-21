# Replication of Figure 1

This package reproduces the main figure in the accepted manuscript from the original GGS-II Wave 1 files. It covers the weighted descriptive analysis and its confidence intervals.

## Data

Obtain the harmonized Stata releases through the [GGP Data Portal](https://www.ggp-i.org/data-catalog/), following its registration and data-access agreement. The analysis uses the file inventory available on 23 August 2026. Use these versions, not newer releases:

| Country or region | Required filename |
| --- | --- |
| Argentina (Buenos Aires) | GGSII_Wave1_AR_V_1_2.dta |
| Austria | GGSII_Wave1_AT_V_1_1.dta |
| Czechia | GGSII_Wave1_CZ_V_2_1.dta |
| Estonia | GGSII_Wave1_EE_V_2_2.dta |
| France | GGSII_Wave1_FR_V_1_0.dta |
| Hong Kong | GGSII_Wave1_HK_V_1_0.dta |
| Iceland | GGSII_Wave1_ISL_V_1_0.dta |
| Netherlands | GGSII_Wave1_NL_V_1_0.dta |
| Taiwan | GGSII_Wave1_TW_V_1_0.dta |
| United Kingdom | GGSII_Wave1_UK_V_1_2.dta |
| Uruguay | GGSII_Wave1_UY_V_1_4.dta |

Place the files in one directory, either together or in their original unzipped subfolders. The script requires exactly one copy of each named file. The supplement documents the context inclusion and exclusion criteria; this script uses that fixed set.

Restricted microdata are not included. The script does not save respondent-level or union-level data. Only figures, aggregate results, sample counts, and file checksums are exported. Dataset citations are in the article; see also the [GGP citation and acknowledgment requirements](https://www.ggp-i.org/data/citation-guidelines/).

## Run

Use Python 3.12. From this folder:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python replicate_main.py --data-dir /path/to/ggs_data
```

On Windows, replace the activation command with `.venv\Scripts\activate` (Command Prompt) or `.venv\Scripts\Activate.ps1` (PowerShell). Replace `/path/to/ggs_data` with your data directory; quote paths containing spaces. Outputs are written to `output/` beside the script. An optional `--output-dir` selects another destination. Running again replaces the generated files in that output directory.

Only NumPy, pandas, and Matplotlib are direct dependencies. There is no dependency on the original project scripts or a prebuilt analytic dataset. Arial is used if installed, with Helvetica and DejaVu Sans as fallbacks; font availability can cause small cosmetic differences across computers.

## Analysis

- Stack current co-residential unions (`dem30a == 1`, `dem30by`, `dem22a`) and up to 20 prior unions (`lhi04_y1`–`lhi04_y20`, `lhi04a_1`–`lhi04a_20`). Retain union starts in 2002–2021, age at start 15–69 (union year minus birth year), a valid meeting channel, and positive respondent weight. GGS special missing values are excluded.
- Meeting codes: online dating = 4; other online = 5; friends/family = 10–11; school/work = 1–2; public/leisure = 3 and 6–9; residual other = 12. All six channels remain in the denominator.
- Periods are 2002–04, 2005–07, 2008–10, 2011–13, 2014–16, 2017–19, and 2020–21. Period refers to co-residence start.
- Each union receives its respondent's GGS weight. Shares are weighted channel totals divided by the total weight in the context-period. No compositional controls or regression model enter the main figure.
- Use 500 Poisson(1) bootstrap replicates, clustered on respondent within context. Each respondent's multiplier applies to all their unions in all periods. Confidence limits are the 2.5th and 97.5th percentiles. Panel B contrasts are calculated within each replicate before taking percentiles, not by subtracting separate confidence limits. These intervals do not incorporate additional PSU or stratum information.
- Panel A displays online-dating shares; Panel B displays 2017–19 minus 2002–04 changes in five channels, omitting residual other from the display only. Estimates are descriptive, not causal displacement effects.

The base seed is 8262026; context seeds add 10007 times the zero-based alphabetical context index, matching the analysis.

## Outputs

- `figure_1.pdf` and `figure_1.png`: the main figure.
- `main_figure_shares.csv`: all six channel shares and confidence limits for every context-period (proportions, not percentages).
- `main_figure_changes.csv`: all six pre-pandemic changes and confidence limits (proportion differences; multiply by 100 for percentage points).
- `sample_counts.csv`: unions and distinct respondents within each context-period. Respondent counts must not be summed across periods.
- `findings.csv`: online-dating changes and the largest offline decline in each context, in percentage points.
- `run_summary.txt` and `input_files.sha256`: sample checks, software versions, random seed, and input-file identifiers.

The expected sample is 38,772 unions from 31,775 respondents, including 3,374 unions in 2020–21.

Tested on 16 September 2026 with Python 3.12.14 and the versions in `requirements.txt`. All shares, pre-pandemic changes, confidence limits, and context-period counts reproduce the published analysis to numerical precision. A separate run in a fresh environment produced identical numerical outputs. Runtime was about 35 seconds on the test machine, excluding installation.
