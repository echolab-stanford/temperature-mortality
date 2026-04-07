# Data Codebook

This document describes all columns expected in the cleaned input `.pq` files for `us_mex_eu_response_attribution.R`.
Files are placed in `data/<COUNTRY>/cleaned/` and are available via Dropbox (see README.md for links).

---

## Files per Country

| File | Description | Temporal resolution |
|---|---|---|
| `<COUNTRY>FullPanel_complete.pq` | Main panel: mortality rates, population, income, temperature polynomials | Monthly (US, MEX) / Weekly (EU) |
| `<COUNTRY>TempBins_monthly_complete.pq` | Pre-computed population-weighted person-days per 1°C temperature bin | Monthly (US, MEX) / Weekly (EU) |
| `<COUNTRY>Daily_temp_covariates.pq` | Daily temperature for MMT calculation and hot/cold attribution | Daily |

---

## `<COUNTRY>FullPanel_complete.pq` — Column Reference

| Column | Type | Unit | Description |
|---|---|---|---|
| `adm_id` | character | — | Unique administrative unit ID (county for US/MEX, NUTS-3 for EU) |
| `year` | integer | — | Calendar year |
| `month` | integer | 1–12 | Month of observation (US, MEX) |
| `week` | integer | 1–53 | ISO week number (EU only) |
| `rate_all` | numeric | deaths / person | All-age mortality rate |
| `rate_young` | numeric | deaths / person | Mortality rate, age < 5 (or country-specific definition) |
| `rate_adult` | numeric | deaths / person | Mortality rate, ages 5–64 (or country-specific) |
| `rate_elderly` | numeric | deaths / person | Mortality rate, age ≥ 65 (or country-specific) |
| `n_deaths_all` | numeric | count | All-age death count |
| `n_deaths_young` | numeric | count | Young age death count |
| `n_deaths_adult` | numeric | count | Adult age death count |
| `n_deaths_elderly` | numeric | count | Elderly age death count |
| `pop_all` | numeric | persons | Total population |
| `pop_young` | numeric | persons | Young age population |
| `pop_adult` | numeric | persons | Adult population |
| `pop_elderly` | numeric | persons | Elderly population |
| `temp_d1_l0` | numeric | °C | Temperature (1st-degree polynomial, lag 0) |
| `temp_d2_l0` | numeric | °C² | Temperature squared (lag 0) |
| `temp_d3_l0` | numeric | °C³ | Temperature cubed (lag 0) |
| `temp_d4_l0` | numeric | °C⁴ | Temperature to 4th power (lag 0) |
| `temp_d1_t1` – `temp_d1_t4` | numeric | °C | Lagged temperature polynomial terms (EU only; t1=lag 1 week, etc.) |
| `precip_d1_l0` | numeric | mm | Precipitation (1st-degree polynomial, lag 0) |
| `precip_d2_l0` | numeric | mm² | Precipitation squared (lag 0) |
| `mean_tm` | numeric | °C | Long-run mean annual temperature (for climate heterogeneity model) |
| `mean_income` | numeric | local currency | Mean household or per-capita income (for income heterogeneity model) |
| `baseline_mean` | numeric | deaths / person | Pre-sample-period mean mortality rate (for baseline heterogeneity model) |
| `temp_percentile` | integer | 0–100 | Percentile bin of `mean_tm` for the admin unit |
| `income_percentile` | integer | 0–100 | Percentile bin of `mean_income` |
| `baseline_percentile` | integer | 0–100 | Percentile bin of `baseline_mean` |
| `decade` | integer | e.g. 1970 | Decade indicator (floor of year to nearest 10) |
| `age_std_rate` | numeric | deaths / person | Age-standardized mortality rate (US, MEX only; SI analysis) |

---

## `<COUNTRY>TempBins_monthly_complete.pq` — Column Reference

| Column | Type | Description |
|---|---|---|
| `adm_id` | character | Administrative unit ID |
| `year` | integer | Calendar year |
| `month` / `week` | integer | Temporal period (matches resolution of FullPanel) |
| `-49.5`, `-48.5`, … `49.5` | numeric | Population × days spent in each 1°C bin centered at that midpoint |

Temperature bin columns are named by their midpoint value (e.g., `-49.5` covers −50 to −49°C).

---

## `<COUNTRY>Daily_temp_covariates.pq` — Column Reference

| Column | Type | Description |
|---|---|---|
| `adm_id` | character | Administrative unit ID |
| `year` | integer | Calendar year |
| `month` | integer | Month |
| `day` | integer | Day of month |
| `temp1` | numeric | Daily mean temperature (°C) — used for MMT identification |
| `temp2` | numeric | Daily temperature squared |
| `temp3` | numeric | Daily temperature cubed |
| `temp4` | numeric | Daily temperature to 4th power |
| `temp_percentile` | integer | Percentile bin of unit's `mean_tm` (matches FullPanel) |
| `income_percentile` | integer | Percentile bin of unit's `mean_income` |
| `baseline_percentile` | integer | Percentile bin of unit's `baseline_mean` |
| `decade` | integer | Decade indicator |

---

## Age Group Definitions

Age cutoffs may vary slightly by country depending on data availability.
Confirm exact cutoffs in the data preparation scripts for each country.

| Group | Approximate ages |
|---|---|
| Young | < 5 years |
| Adult | 5–64 years |
| Elderly | ≥ 65 years |

---

## Adding a New Country

A new country's cleaned data must conform to this schema before the analysis script can run.
Key checklist:

- [ ] All rate and count columns present for each age group available
- [ ] `adm_id` is a stable unique identifier
- [ ] Temperature polynomial columns (`temp_d1_l0` through `temp_d4_l0`) pre-computed
- [ ] `decade` column computed as `floor(year / 10) * 10`
- [ ] Percentile columns computed within the country's own distribution
- [ ] Temporal resolution matches `time_agg_var` in `config/countries.R`
