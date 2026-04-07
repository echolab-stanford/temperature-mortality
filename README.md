# temperature-mortality

Code repository for "Understanding and addressing temperature impacts on mortality". This repo contains the main scripts to replicate all statistical results (Figures 1-6) in the paper. The underlying EU and Mexico mortality data are public and provided below with merged temperature data; the US data are restricted use and available by application to the CDC and are not provided. This project is licensed under the terms of the MIT license.

---

## Repository Structure

```
temperature-mortality/
|-- config/
|   +-- countries.R                        # Per-country model configuration (single source of truth)
|-- data/
|-- docs/
|   +-- CODEBOOK.md                        # Column reference for all cleaned input files
|-- processed/                             # Model outputs (.rds) -- not tracked by git
|-- fig/                                   # Generated figures -- not tracked by git
|-- run_all.R                              # Master orchestration script (start here)
|-- us_mex_eu_response_attribution.R       # Step 1: estimate models, produce .rds outputs
|-- us_mex_eu_response_attribution_plot.R  # Step 2: generate Figs 1-3
|-- age_lag_response_us_mexico.R           # Step 3: Fig 4
|-- by_cause_us_mexico.R                   # Step 3: Fig 5
|-- us_mex_eu_maps.R                       # Step 3: Fig 6
|-- temperature-mortality.Rproj            # Open this in RStudio/Positron
+-- renv.lock                              # Exact package versions for reproducibility
```

> `data/` and `processed/` are listed in `.gitignore` and are never committed.
> Always open the project via `temperature-mortality.Rproj` so `here::here()` resolves paths correctly.

---

## System Requirements

**Software:**
- R >= 4.5.1
- `renv` (for package management -- installed automatically if missing)

**Operating System:** Tested on macOS Sequoia 15.6.1. Should work on any OS supporting R.

**Hardware:**
- >= 24 GB RAM recommended (bootstrap steps load large panels into memory)
- >= 15 GB free disk space for data files

**Approximate runtimes** (modern workstation, 24 GB RAM):

| Task | Approximate time |
|---|---|
| Package installation via `renv::restore()` | 5-10 min |
| Analysis script per country (`si_folder = "winsor"`) | 2-4 hours |
| Plot script (Figs 1-3) | 5-10 min |
| Supplementary analysis per country per `si_folder` | 30-60 min |

---

## Installation

```bash
git clone https://github.com/echolab-stanford/temperature-mortality.git
cd temperature-mortality
```

Open `temperature-mortality.Rproj` in RStudio or Positron, then restore exact package versions:

```r
# Install renv if needed
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

renv::restore()
```

---

## Replication: Quickest Path (Figs 1-3 only)

If you want to skip re-running the bootstrap (2-4 hours per country), download pre-computed `.rds` files directly.

### Step 0 -- Download pre-computed processed outputs

Create the directory structure, then download each country's `.rds` files into `processed/<COUNTRY>/`:

```bash
mkdir -p processed/US processed/EU processed/MEX
```

| Country | Dropbox link |
|---|---|
| US | [LINK TBD] |
| EU | [LINK TBD] |
| MEX | [LINK TBD] |

Place all `.rds` files directly inside `processed/<COUNTRY>/` (e.g., `processed/US/US_Pooled_results.rds`).

### Step 1 -- Generate Figs 1-3

```r
source(here::here("us_mex_eu_response_attribution_plot.R"))
```

Outputs saved to `fig/combined/`:
- `fig_1.pdf` -- Pooled dose-response curves + % attributable deaths
- `fig_2.pdf` -- Heterogeneity by age, climate, income, decade
- `fig_3.pdf` -- NOTE: needs rewrite (see TODO comment in script)

---

## Replication: Full Pipeline (all figures, from raw data)

### Step 0 -- Download raw data

#### Figures 1-3 data

Place files in `data/<COUNTRY>/cleaned/`. See [docs/CODEBOOK.md](docs/CODEBOOK.md) for expected column schema.

**Mexico** (`data/MEX/cleaned/`):

1. `MEXFullPanel_complete.pq` [~1.57 GB](https://www.dropbox.com/scl/fi/iqw0eet98ivok94fbita4/MEXFullPanel_complete.pq?rlkey=v9c4ufv8a13ctutdnbug3v5uc&st=hqzn0l4t&dl=0)
2. `MEXDaily_temp_covariates.pq` [~933 MB](https://www.dropbox.com/scl/fi/fuqk3iq4rcmh51kc07avh/MEXDaily_temp_covariates.pq?rlkey=bk6us32lwds3x1f3jlgdzrh2d&st=8vtvwjof&dl=0)
3. `MEXTempBins_monthly_complete.pq` [~6.5 MB](https://www.dropbox.com/scl/fi/mwf7h6v1l940ttfy3ddbw/MEXTempBins_monthly_complete.pq?rlkey=w6c27kvnttswuhq98wh88tbfj&st=xelsvetb&dl=0)

**EU** (`data/EU/cleaned/`):

1. `EUFullPanel_complete.pq` [~2.7 GB](https://www.dropbox.com/scl/fi/108e8ohy269thyhjsp90j/EUFullPanel_complete.pq?rlkey=r552u8o05x5k8eedpt8d7e2io&st=3qosof02&dl=0)
2. `EUDaily_temp_covariates.pq` [~402 MB](https://www.dropbox.com/scl/fi/28mzqrpn8d225440ytigp/EUDaily_temp_covariates.pq?rlkey=10460u0qlbehc27kdnjq2yu4j&st=765v3639&dl=0)
3. `EUTempBins_monthly_complete.pq` [~12.1 MB](https://www.dropbox.com/scl/fi/zr35tn5pczh1vyljyq7o2/EUTempBins_monthly_complete.pq?rlkey=46683aso0umbizbfb7zhixc1z&st=oh5pnx0v&dl=0)

**US** (`data/US/cleaned/`):

1. `USFullPanel_complete.pq` -- *Restricted use; apply to CDC. Not available for public download.*
2. `USDaily_temp_covariates.pq` [~1.77 GB](https://www.dropbox.com/scl/fi/ww8h0u2yuv4gb6oq1kyfl/USDaily_temp_covariates.pq?rlkey=vc8izr0i7jjmh5ywo19tz6t2x&e=1&dl=0)
3. `USTempBins_monthly_complete.pq` [~27 MB](https://www.dropbox.com/scl/fi/q5qaridb8zpea11ki861y/USTempBins_monthly_complete.pq?rlkey=om98hd2j117u8ckbrjybd88u5&dl=0)

#### Figures 4-6 data

**US** (`data/US/`):

1. `singleage_adjustment.csv` [~8 KB](https://www.dropbox.com/scl/fi/74053cyzuoha47j3dfa0k/singleage_adjustment.csv?rlkey=266eo3eugc0chixzfurmk7a3n&dl=0)
2. `pop_panel_singleage.pq` [~28 MB](https://www.dropbox.com/scl/fi/zdes5va6xlmpdbr6dk9n6/pop_panel_singleage.pq?rlkey=f4len1nozbshac0q3bzg168ps&dl=0)
3. `us_temp.pq` [~1.96 GB](https://www.dropbox.com/scl/fi/4zjz4gos7sd8bbul1mw2y/us_temp.pq?rlkey=k2xv5hzwa6xslnmcyfgnnw0wn&dl=0)
4. `us_precip.pq` [~823 MB](https://www.dropbox.com/scl/fi/ppmm5er5k86c1l7mq14tt/us_precip.pq?rlkey=3wndfp4w202mksbjlousik3c5&dl=0)
5. `fully_monthly.pq` -- *Restricted use; not available for public download.*
6. `CONUS_consistent_counties.gpkg` [~16 MB](https://www.dropbox.com/scl/fi/ls3dqkf8xfpxqjbvp4ksc/CONUS_consistent_counties.gpkg?rlkey=7oftm4mxve65vqe5jlrxok1aj&dl=0)

**Mexico** (`data/MEX/`):

1. `daily_death_pop_income.pq` [~367 MB](https://www.dropbox.com/scl/fi/szfowpm04zjmizvw866yq/daily_death_pop_income.pq?rlkey=h5bjwtmf7kew0zj0f5by38ic8&dl=0)
2. `mexico_temp.pq` [~1.52 GB](https://www.dropbox.com/scl/fi/hv37en54tbyhgbaogfzuc/mexico_temp.pq?rlkey=7nul80l3auua579j3tmc4hq2v&dl=0)
3. `mexico_precip.pq` [~653 MB](https://www.dropbox.com/scl/fi/g5ikw804xf4l2mtxcgzr9/mexico_precip.pq?rlkey=j2z2y136k5b2d0j0j33ballmw&dl=0)
4. `adm2_to_order.csv` [~29 KB](https://www.dropbox.com/scl/fi/06my16kddud0pg9pb1jag/adm2_to_order.csv?rlkey=xmxfo2ymhp45tw2r0j3gt0vca&dl=0)
5. `monthly_death_pop_income.pq` [~302 MB](https://www.dropbox.com/scl/fi/kdmxy0kyblzdpjc8rudtz/monthly_death_pop_income.pq?rlkey=487n9d7gclc1hcv65ufwkfsbl&dl=0)
6. `monthly_weather.pq` [~2.33 GB](https://www.dropbox.com/scl/fi/a25gcczxtx7jhu7vv1g0p/monthly_weather.pq?rlkey=5vpl3qdxu0jijzl3d6gni2dlc&dl=0)
7. `full_monthly.pq` [~1.69 GB](https://www.dropbox.com/scl/fi/emswcm2vhz26cpmgao44r/full_monthly.pq?rlkey=qt0jlm01yarkogb2bn5q05342&dl=0)
8. `full_monthly_small.pq` [~54 MB](https://www.dropbox.com/scl/fi/siwphhybw5ubuttfw7r9e/full_monthly_small.pq?rlkey=tlbdwnf0yfwlu211ubfd6nlo0&dl=0)
9. `mexico.gpkg` [~6.84 MB](https://www.dropbox.com/scl/fi/ksbdkrjd76ayhqw3o64ov/mexico.gpkg?rlkey=3xq7rnbl6c1g5vfa78zwgvw2i&dl=0)

**EU** (`data/EU/`):

1. `full_weekly.pq` [~146 MB](https://www.dropbox.com/scl/fi/txd2ajouz26hjo4u2x7yt/full_weekly.pq?rlkey=kre5ziioq7c9g5uw38yo2kubm&dl=0)
2. `eu_nuts3_2024.gpkg` [~29 MB](https://www.dropbox.com/scl/fi/kr61s2lurcsia3jmldjo7/eu_nuts3_2024.gpkg?rlkey=au9iwm8gj0av6xvx4qdeevf4k&dl=0)
3. `eu_avg_temps.pq` [~582 KB](https://www.dropbox.com/scl/fi/e4d0h3kam7el7nb7syrzx/eu_avg_temps.pq?rlkey=pusolhe8r17g39fnsjw4vtlz4&dl=0)
4. `eu_temp_limits.pq` [~24 KB](https://www.dropbox.com/scl/fi/73pgc9dj4c4xoqbse7fhy/eu_temp_limits.pq?rlkey=bsbpt29bvbd67tbhy632w2btb&dl=0)

**World map** (`data/`):

1. `geoBoundariesCGAZ_ADM0.gpkg` [~155 MB](https://www.dropbox.com/scl/fi/8geqx7i7q4m90snvdgmsh/geoBoundariesCGAZ_ADM0.gpkg?rlkey=81mz0k6nbnhbp1o0uhb8xtdhj&dl=0)

### Step 1 -- Run analysis script (one country at a time)

Open `us_mex_eu_response_attribution.R`. Set `country_name` and `si_folder` at the top, then source.

For Figs 1-3, run each country with `si_folder = "winsor"`:

```r
# In us_mex_eu_response_attribution.R:
country_name <- "US"    # then repeat for "EU" and "MEX"
si_folder    <- "winsor"
```

Outputs are saved to `processed/<COUNTRY>/`.

### Step 2 -- Generate Figs 1-3

```r
source(here::here("us_mex_eu_response_attribution_plot.R"))
```

### Step 3 -- Generate Figs 4-6

```r
source(here::here("age_lag_response_us_mexico.R"))  # Fig 4
source(here::here("by_cause_us_mexico.R"))           # Fig 5
source(here::here("us_mex_eu_maps.R"))               # Fig 6
```

---

## Figure Index

| Figure | Script | Function / section | si_folder |
|---|---|---|---|
| Fig 1 | `us_mex_eu_response_attribution_plot.R` | `plot_pooled_fig()` | winsor (default) |
| Fig 2 | `us_mex_eu_response_attribution_plot.R` | `plot_response()` + `plot_annual_deaths()` | winsor (default) |
| Fig 3 | `us_mex_eu_response_attribution_plot.R` | `plot_deaths_at_bins()` -- needs rewrite | winsor (default) |
| Fig 4 | `age_lag_response_us_mexico.R` | -- | -- |
| Fig 5 | `by_cause_us_mexico.R` | -- | -- |
| Fig 6 | `us_mex_eu_maps.R` | -- | -- |
| SI Fig S1 | `us_mex_eu_response_attribution_plot.R` | `plot_pooled_fig()` + `plot_response()` | log_rate |
| SI Fig S2 | `us_mex_eu_response_attribution_plot.R` | `plot_response()` + `plot_annual_deaths()` | age_standardized |
| SI Fig S9 | `us_mex_eu_response_attribution_plot.R` | decade_fixed1970 section | winsor, US only |

---

## Adding a New Country

1. Uncomment the country entry in `config/countries.R` and fill in the parameters.
2. Place processed `.rds` files in `processed/<COUNTRY>/` (one file per model in `models_to_run`).
3. The plot script will automatically include the new country in any figure listed in its `figs` field.
4. See [docs/CODEBOOK.md](docs/CODEBOOK.md) for the expected input schema if re-running from raw data.
