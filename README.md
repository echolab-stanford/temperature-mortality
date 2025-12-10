# temperature-mortality

Code repository for "Understanding and addressing temperature impacts on mortality". This repo contains main scripts to replicate main statistical results (Figures 1-6) in the paper. The underlying EU and Mexico mortality data are public and provided below with merged temperature data; the US data are restricted use and available by application to the CDC and are not provided. This project is licensed under the terms of the MIT license. 

### System Requirements

**Software Dependencies:**
- R version 4.5.1 (or later)
- Required R packages:

  - `arrow` - for reading parquet files
  - `fastverse`
  - `tidyverse`
  - `fixest`
  - `matrixStats` 
  - `ISOweek`
  - `pammtools`
  - `tictoc`
  - `sf`
  - `ggnewscale`

**Operating System:**

- Tested on: macOS Sequoia 15.6.1
- Should work on any OS that supports R

**Hardware:**

- Minimum 24GB RAM recommended (for large datasets)
- At least 15GB free disk space for data files

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

3. Typical install time ~ 5-10 mins
- Cloning the repo itself is very fast (<1 minute) since it's just code
- The main time is installing R packages (5-10 minutes typical for common packages)
- If users already have some packages installed, it will be faster

### Replication Instructions

### Figures 1, 2, and 3

#### Step 1: Download and Organize Data

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

1. `USFullPanel_complete.pq` - The US data are restricted use and available by application to the CDC and are not provided.
2. `USDaily_temp_covariates.pq` [~1.77GB](https://www.dropbox.com/scl/fi/ww8h0u2yuv4gb6oq1kyfl/USDaily_temp_covariates.pq?rlkey=vc8izr0i7jjmh5ywo19tz6t2x&e=1&dl=0)
3. `USTempBins_monthly_complete.pq` [~27.03MB](https://www.dropbox.com/scl/fi/q5qaridb8zpea11ki861y/USTempBins_monthly_complete.pq?rlkey=om98hd2j117u8ckbrjybd88u5&dl=0)



#### Step 2: Generate Analysis Dataframes

Run `us_mex_eu_response_attribution.R` to generate dataframes for plotting figures. The resulting outputs including bootstrap and plot dataframes will be saved in `processed/[country_name]`.

The `us_mex_eu_response_attribution_plots.R` script requires several configurations. After each configuration change, the script needs to be run again:

**For Figure 1, 2, and 3:**

- Set `country_name` = ["US", "EU", or "MEX"]
- Set `si_folder` = ["winsor"]
- Dataframes will be saved in `processed/[country_name]`


#### Step 3: Generate Figures

After completing Step 2, run `MortalityTemp_combined_plots_v2.R` to generate figures. The resulting figures will be saved in `fig/combined/`:

- Figure 1: `fig_1.pdf`
- Figure 2: `fig_2.pdf`
- Figure 3: `fig_3.pdf`

### Figures 4, 5, and 6

#### Step 1: Download and Organize Data

Download each country's data and place files in `data/[country_name]`, where `country_name` is MEX or US.

US data: 

1. `singleage_adjustment.csv`[~8.35KB](https://www.dropbox.com/scl/fi/74053cyzuoha47j3dfa0k/singleage_adjustment.csv?rlkey=266eo3eugc0chixzfurmk7a3n&dl=0)
2. `pop_panel_singleage.pq`[~27.96MB](https://www.dropbox.com/scl/fi/zdes5va6xlmpdbr6dk9n6/pop_panel_singleage.pq?rlkey=f4len1nozbshac0q3bzg168ps&dl=0)
3. `us_temp.pq` [~1.96GB](https://www.dropbox.com/scl/fi/4zjz4gos7sd8bbul1mw2y/us_temp.pq?rlkey=k2xv5hzwa6xslnmcyfgnnw0wn&dl=0)
4. `us_precip.pq` [~823.33MB](https://www.dropbox.com/scl/fi/ppmm5er5k86c1l7mq14tt/us_precip.pq?rlkey=3wndfp4w202mksbjlousik3c5&dl=0)
5. `fully_monthly.pq` - data is not public so cannot be shared for replication
6. `CONUS_consistent_counties` [~15.98MB](https://www.dropbox.com/scl/fi/ls3dqkf8xfpxqjbvp4ksc/CONUS_consistent_counties.gpkg?rlkey=7oftm4mxve65vqe5jlrxok1aj&dl=0)

Mexico data:

1. `daily_death_pop_income.pq` [~366.5MB](https://www.dropbox.com/scl/fi/szfowpm04zjmizvw866yq/daily_death_pop_income.pq?rlkey=h5bjwtmf7kew0zj0f5by38ic8&dl=0)
2. `mexico_temp.pq` [~1.52GB](https://www.dropbox.com/scl/fi/hv37en54tbyhgbaogfzuc/mexico_temp.pq?rlkey=7nul80l3auua579j3tmc4hq2v&dl=0)
3. `mexico_precip.pq` [~653.19MB](https://www.dropbox.com/scl/fi/g5ikw804xf4l2mtxcgzr9/mexico_precip.pq?rlkey=j2z2y136k5b2d0j0j33ballmw&dl=0)
4. `adm2_to_order.csv` [~29.42KB](https://www.dropbox.com/scl/fi/06my16kddud0pg9pb1jag/adm2_to_order.csv?rlkey=xmxfo2ymhp45tw2r0j3gt0vca&dl=0)
5. `monthly_death_pop_income.pq` [~302.46MB](https://www.dropbox.com/scl/fi/kdmxy0kyblzdpjc8rudtz/monthly_death_pop_income.pq?rlkey=487n9d7gclc1hcv65ufwkfsbl&dl=0)
6. `monthly_weather.pq` [~2.33GB](https://www.dropbox.com/scl/fi/a25gcczxtx7jhu7vv1g0p/monthly_weather.pq?rlkey=5vpl3qdxu0jijzl3d6gni2dlc&dl=0)
7. `full_monthly.pq` [~1.69GB](https://www.dropbox.com/scl/fi/emswcm2vhz26cpmgao44r/full_monthly.pq?rlkey=qt0jlm01yarkogb2bn5q05342&dl=0)
8. `full_monthly_small.pq` [~54.42MB](https://www.dropbox.com/scl/fi/siwphhybw5ubuttfw7r9e/full_monthly_small.pq?rlkey=tlbdwnf0yfwlu211ubfd6nlo0&dl=0)
9. `mexico.gpkg` [~6.84MB](https://www.dropbox.com/scl/fi/ksbdkrjd76ayhqw3o64ov/mexico.gpkg?rlkey=3xq7rnbl6c1g5vfa78zwgvw2i&dl=0)

Europe data: 

1. `full_weekly.pq` [~146.MB](https://www.dropbox.com/scl/fi/txd2ajouz26hjo4u2x7yt/full_weekly.pq?rlkey=kre5ziioq7c9g5uw38yo2kubm&dl=0)
2. `eu_nuts3_2024.gpkg` [~29.46MB](https://www.dropbox.com/scl/fi/kr61s2lurcsia3jmldjo7/eu_nuts3_2024.gpkg?rlkey=au9iwm8gj0av6xvx4qdeevf4k&dl=0)
3. `eu_avg_temps.pq` [~582.34KB](https://www.dropbox.com/scl/fi/e4d0h3kam7el7nb7syrzx/eu_avg_temps.pq?rlkey=pusolhe8r17g39fnsjw4vtlz4&dl=0)
4. `eu_temp_limits.pq` [~23.98KB](https://www.dropbox.com/scl/fi/73pgc9dj4c4xoqbse7fhy/eu_temp_limits.pq?rlkey=bsbpt29bvbd67tbhy632w2btb&dl=0)

World map:

1. `geoBoundariesCGAZ_ADM0` [~154.63](https://www.dropbox.com/scl/fi/8geqx7i7q4m90snvdgmsh/geoBoundariesCGAZ_ADM0.gpkg?rlkey=81mz0k6nbnhbp1o0uhb8xtdhj&dl=0)

#### Step 2: Generate Figures

Run the following scripts to generate figures. The resulting figures will be saved in `fig/combined/`:

- For Figure 4: Run `age_lag_response_us_mexico.R` -> outputs `lag_age_response.pdf`
- For Figure 5: Run `by_cause_us_mexico.R` -> outputs `cause.pdf`
- For Figure 6: Run `us_mex_eu_maps.R` -> outputs `county_maps_*`
