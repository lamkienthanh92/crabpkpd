# CRAB Combination PK/PD Simulator

Browser-based Monte Carlo simulator for colistin, meropenem and ampicillin/sulbactam against carbapenem-resistant *Acinetobacter baumannii* (CRAB). For each formula (colistin + meropenem, colistin + ampicillin/sulbactam, meropenem + ampicillin/sulbactam, and the triple combination) it estimates the probability that every drug reaches its own target concentration in the same virtual patient, in plasma, lung epithelial lining fluid, cerebrospinal fluid and soft tissue, with safety limits on meropenem trough and colistin steady-state concentrations.

**Research use only. Not for choosing doses in individual patients.**

- Live app: https://crabpkpd.netlify.app/
- Source: https://github.com/lamkienthanh92/crabpkpd
- Archive (all versions): https://doi.org/10.5281/zenodo.24196427

## Run

Open `index.html` in a browser. No build step and no network access are needed (React, ReactDOM and htm are in `vendor/`).

## Deploy on Netlify

Connect the GitHub repository (or drag the folder to Netlify Drop). `netlify.toml` publishes the repository root; the build command is empty.

## Contents

| Path | Content |
| --- | --- |
| `index.html` | the complete application (unminified source inside) |
| `vendor/` | React 18.3.1, ReactDOM 18.3.1, htm 3.1.1 |
| `data/isolates_39.csv` | checkerboard concentrations of 39 isolates (isolates 15, 26, 27, 28: meropenem in the triple combination corrected to 2 µg/mL after a data-entry error; original entries are included) |
| `engine/` | simulation engine (`core_en.js`) and the Node.js script that reproduces the manuscript's results: `cd engine && node run_formulas.js` (writes `results_formulas.json`) |
| `docs/Supplementary_Tables.xlsx` | model parameters and sources, verification against published values, all results |

Options of `run_formulas.js` (environment variables; defaults reproduce the base case): `SCEN` colistin scenarios as `name:free fraction:target multiplier,...` (default `S0:0.49:1,S3:0.085:4`; the four-scenario analysis uses `S0:0.49:1,S1:0.085:1,S2:0.49:4,S3:0.085:4`), `MERO_T` meropenem time target (default 0.4; 0.6 and 1.0 in the sensitivity analyses), `ELF_K` and `CSF_K` multipliers of the meropenem lung and CSF ratios (0.5556/1.4815 and 0.5556/1.7778 give 0.15/0.05 and 0.40/0.16), `IND_CMS=1` colistimethate dosed by each patient's own clearance, `ISO` isolate file (`iso35.json` omits isolates 15, 26, 27, 28), `OUT` output file. Colistin lung and soft-tissue concentrations are 0.26 and 0.12 times the free plasma concentration (pig microdialysis, Rottbøll & Friis 2016); CSF is 0.09 times total plasma (Ziaka 2013). Regimens are selected by renal stratum (nominal clearance 155, 110, 70, 40, 20 mL/min); pharmacokinetics use each virtual patient's own clearance.

## Model sources

Meropenem: Crandon et al., Intensive Care Med 2011. Sulbactam: Jaruratanasirikul et al., Antimicrob Agents Chemother 2016. Colistin: Nation et al., Clin Infect Dis 2017. All parameters are listed in `docs/Supplementary_Tables.xlsx` (sheet S1).

## Citation

Cite the software with DOI 10.5281/zenodo.24196427 (all versions; see `CITATION.cff`). Please also cite: Huynh T, Lam T, Pham N, Luong L, Jusselme MD. J Med Microbiol 2025;74:002086, https://doi.org/10.1099/jmm.0.002086 (source of the isolates and checkerboard method).

## Licence

Code: MIT (`LICENSE`). Data in `data/`: CC BY 4.0. Third-party components: `THIRD_PARTY_LICENSES.md`.

## Before publishing

- Replace the author ORCID placeholder in `CITATION.cff`.
- On Zenodo, enable the GitHub integration for the repository and create a release (for example `v1.0.0`); Zenodo reads `.zenodo.json`. Add the Zenodo DOI to the manuscript's data availability statement.
- Add the manuscript DOI to `.zenodo.json` (`related_identifiers`) once available.


## Sensitivity analysis in the app

The *Sensitivity analysis* tab re-runs the selected regimen under the colistin free fraction and target scenarios (A to D), meropenem targets of 40, 60 and 100%, low, base and high meropenem lung and CSF penetration, and two sulbactam lung ratios, and exports the table as CSV.
