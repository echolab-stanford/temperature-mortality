################################################################################
# run_all.R
#
# Master orchestration script for the temperature-mortality analysis.
# Documents the full execution order for reproducing all main figures.
#
# BEFORE RUNNING:
#   1. Open the project via temperature-mortality.Rproj (ensures here() works)
#   2. Run renv::restore() to install exact package versions
#   3. Download processed .rds files from Dropbox into processed/<COUNTRY>/
#      (see README.md for links), OR download raw data and run Step 1 below
#
# APPROXIMATE RUNTIMES (on a modern workstation, 24GB RAM):
#   Step 1 (analysis, per country):  ~2–4 hours  (bootstrap-heavy)
#   Step 2 (plots):                  ~5–10 min
#   Step 3 (supplementary figures):  ~30–60 min per country
#
# All scripts source config/countries.R internally.
################################################################################

################################################################################
# STEP 1: Generate processed .rds files from raw data
# Run us_mex_eu_response_attribution.R once per country.
# Skip this step if you downloaded pre-computed .rds files from Dropbox.
#
# Set country_name and si_folder inside the script before sourcing.
#
# Required si_folder runs for full replication:
#   - "winsor"           -> main figures (Figs 1-3)
#   - "log_rate"         -> SI Fig S1
#   - "age_standardized" -> SI Fig S2 (US and MEX only)
################################################################################

# for (country in c("US", "EU", "MEX")) {
#   message("=== Running analysis for: ", country, " ===")
#   # Open us_mex_eu_response_attribution.R, set country_name = country,
#   # set si_folder = "winsor", then source:
#   # source(here::here("us_mex_eu_response_attribution.R"))
# }

################################################################################
# STEP 2: Generate main figures (Figs 1, 2, 3)
# Reads from processed/<COUNTRY>/*_results.rds
# Outputs to fig/combined/
################################################################################

# source(here::here("us_mex_eu_response_attribution_plot.R"))

# Expected outputs:
#   fig/combined/fig_1.pdf  -- Pooled dose-response + % attributable deaths
#   fig/combined/fig_2.pdf  -- Heterogeneity by age, climate, income, decade
#   fig/combined/fig_3.pdf  -- TODO: needs rewrite (see TODO in script)

################################################################################
# STEP 3: Generate supplementary figures
# Run the main analysis script with each si_folder before running plot script.
################################################################################

# SI Fig S1 (log rate):
#   Set si_folder = "log_rate" in us_mex_eu_response_attribution.R for each country
#   Then run the log_rate section at the bottom of us_mex_eu_response_attribution_plot.R

# SI Fig S2 (age standardized, US + MEX only):
#   Set si_folder = "age_standardized" for US and MEX
#   Then run the age_standardized section at the bottom of us_mex_eu_response_attribution_plot.R

################################################################################
# STEP 4: Additional figures (Figs 4, 5, 6)
# These scripts use different raw data inputs (see README.md § Figures 4-6).
################################################################################

# Fig 4 -- Distributed lag and age-specific response curves (US + MEX)
# source(here::here("age_lag_response_us_mexico.R"))

# Fig 5 -- Cause-specific mortality (US + MEX)
# source(here::here("by_cause_us_mexico.R"))

# Fig 6 -- Geographic maps of temperature-mortality sensitivity (US + MEX + EU)
# source(here::here("us_mex_eu_maps.R"))
