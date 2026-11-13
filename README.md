# CRAB Combination PK/PD Simulator

Browser-based Monte Carlo simulator for colistin, meropenem and ampicillin/sulbactam against carbapenem-resistant *Acinetobacter baumannii* (CRAB). For each formula (colistin + meropenem, colistin + ampicillin/sulbactam, meropenem + ampicillin/sulbactam, and the triple combination) it estimates the probability that every drug reaches its own target concentration in the same virtual patient, in plasma, lung epithelial lining fluid, cerebrospinal fluid and soft tissue, with safety limits on meropenem trough and colistin steady-state concentrations.

**Research use only. Not for choosing doses in individual patients.**

## Run

Open `index.html` in a browser. No build step and no network access are needed (React, ReactDOM and htm are in `vendor/`).

## Deploy on Netlify

Connect the GitHub repository (or drag the folder to Netlify Drop). `netlify.toml` publishes the repository root; the build command is empty.

## Contents

| Path | Content |
| --- | --- |
| `index.html` | the complete application (unminified source inside) |
| `vendor/` | React 18.3.1, ReactDOM 18.3.1, htm 3.1.1 |
| `data/isolates_39.csv` | checkerboard concentrations of 39 isolates (isolates 15, 26, 27, 28: meropenem in the triple combination set to 2 µg/mL; original records are included) |
| `engine/` | simulation engine (`core_en.js`) and the Node.js script that reproduces the manuscript's results: `cd engine && node run_formulas.js` (writes `results_formulas.json`) |
| `docs/Supplementary_Tables.xlsx` | model parameters and sources, verification against published values, all results |

## Model sources

Meropenem: Crandon et al., Intensive Care Med 2011. Sulbactam: Jaruratanasirikul et al., Antimicrob Agents Chemother 2016. Colistin: Nation et al., Clin Infect Dis 2017. All parameters are listed in `docs/Supplementary_Tables.xlsx` (sheet S1).

## Citation

See `CITATION.cff`. Please also cite: Huynh T, Lam T, Pham N, Luong L, Jusselme MD. J Med Microbiol 2025;74:002086, https://doi.org/10.1099/jmm.0.002086 (source of the isolates and checkerboard method).

## Licence

Code: MIT (`LICENSE`). Data in `data/`: CC BY 4.0. Third-party components: `THIRD_PARTY_LICENSES.md`.

## Before publishing

- Replace the author ORCID placeholder in `CITATION.cff`.
- On Zenodo, enable the GitHub integration for the repository and create a release (for example `v1.0.0`); Zenodo reads `.zenodo.json`. Add the Zenodo DOI to the manuscript's data availability statement.
- Add the manuscript DOI to `.zenodo.json` (`related_identifiers`) once available.
