# temperature-mortality

Code repository for "Understanding and addressing temperature impacts on mortality". This repo contains main scripts to replicate main statistical results (Figure 1, 2, 3, S1, S2, and S9) in the paper. The underlying EU and Mexico mortality data are public and provided below with merged temperature data; the US data are restricted use and available by application to the CDC and are not provided.

### System Requirements

- All software dependencies and operating systems (including version numbers). 
- Versions the software has been tested on. 
- Any required non-standard hardware.

### Installation Guide

Instructions and typical install time on a normal desktop computer.

1. Clone this repository:

   ```bash
   git clone https://github.com/yourusername/temperature-mortality.git
   cd temperature-mortality
   ```

2. Set working directory in R:
- Recommended: Open temperature-mortality.Rproj in RStudio (create via File → New Project → Existing Directory if needed)
- Alternative: Use `setwd()` to point to your cloned directory: `setwd("~/path/to/temperature-mortality")`

### Replication Instructions

### Step 1: Download and Organize Data

Download each country's data and place files in `data/[country_name]/cleaned`, where `country_name` is MEX, EU, or US. For each country, the data include three parquet files:

- `[country_name]Daily_temp_covariates.pq` contains daily observed temperatures
- `[country_name]FullPanel_complete.pq` contains aggregated and merged monthly panel (weekly for EU)
- `[country_name]TempBins_monthly_complete.pq` contains pre-calculated person days population exposure

**Mexico Data:**

1. `MEXFullPanel_complete.pq` [~1.57GB](https://www.dropbox.com/scl/fi/iqw0eet98ivok94fbita4/MEXFullPanel_complete.pq?rlkey=v9c4ufv8a13ctutdnbug3v5uc&st=hqzn0l4t&dl=0)
2. `MEXDaily_temp_covariates.pq` [~933.4MB](https://www.dropbox.com/scl/fi/fuqk3iq4rcmh51kc07avh/MEXDaily_temp_covariates.pq?rlkey=bk6us32lwds3x1f3jlgdzrh2d&st=8vtvwjof&dl=0)
3. `MEXTempBins_monthly_complete.pq` [~6.5MB](https://www.dropbox.com/scl/fi/mwf7h6v1l940ttfy3ddbw/MEXTempBins_monthly_complete.pq?rlkey=w6c27kvnttswuhq98wh88tbfj&st=xelsvetb&dl=0)

**EU Data:**

1. `EUFullPanel_complete.pq` [~2.7GB](https://www.dropbox.com/scl/fi/108e8ohy269thyhjsp90j/EUFullPanel_complete.pq?rlkey=r552u8o05x5k8eedpt8d7e2io&st=3qosof02&dl=0)
2. `EUDaily_temp_covariates.pq` [~402MB](https://www.dropbox.com/scl/fi/28mzqrpn8d225440ytigp/EUDaily_temp_covariates.pq?rlkey=10460u0qlbehc27kdnjq2yu4j&st=765v3639&dl=0)
3. `EUTempBins_monthly_complete.pq` [~12.1MB](https://www.dropbox.com/scl/fi/zr35tn5pczh1vyljyq7o2/EUTempBins_monthly_complete.pq?rlkey=46683aso0umbizbfb7zhixc1z&st=oh5pnx0v&dl=0)

**US Data:**

The US data are restricted use and available by application to the CDC and are not provided.

### Step 2: Generate Analysis Dataframes

Run `MortalityTemp_combined_analysis_v2.R` to generate dataframes for plotting figures in the paper. The resulting outputs including bootstrap and plot dataframes will be saved in `processed/[country_name]`.

The `MortalityTemp_combined_analysis_v2.R` script requires several configurations. After each configuration change, the script needs to be run again:

**For Figure 1, Figure 2, Figure 3, Figure S9:**

- Set `country_name` = ["US", "EU", or "MEX"]
- Set `si_folder` = ["winsor"]
- Dataframes will be saved in `processed/[country_name]`

**For Figure S1:**

- Set `country_name` = ["US", "EU", or "MEX"]
- Set `si_folder` = ["log_rate"]
- Dataframes will be saved in `processed/[country_name]/log_rate` (log mortality rate as outcome)

**For Figure S2:**

- Set `country_name` = ["US", "EU", or "MEX"]
- Set `si_folder` = ["age_standardized"]
- Dataframes will be saved in `processed/[country_name]/age_standardized` (age standardized rate as outcome)

### Step 3: Generate Figures

After completing Step 2, run `MortalityTemp_combined_plots_v2.R` to generate figures. The resulting figures will be saved in `fig/combined/`:

- Figure 1: `fig_1.pdf`
- Figure 2: `fig_2.pdf`
- Figure 3: `fig_3.pdf`
- Figure S1: `fig_2_log_rate.pdf`
- Figure S2: `fig_2_age_std_us_mex.pdf`
- Figure S9: `decade_fixed1970.pdf`

### Expected Output and Runtime

Expected output and expected run time for execution on a normal desktop computer: ...
