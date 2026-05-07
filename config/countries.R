################################################################################
# config/countries.R
#
# Central configuration for all countries in the analysis.
# Each entry in `country_config` defines:
#
#   time_agg_var   - temporal aggregation unit used in the panel ("week" or "month")
#   fixed_effects  - fixest FE formula string
#   rhs_vars       - character vector of RHS temperature/precip polynomial terms
#   models_to_run  - which models to estimate for this country
#   figs           - which main figures this country appears in
#
# To add a new country:
#   1. Add an entry below (use the TEMPLATE block as a guide)
#   2. Ensure processed/<COUNTRY>/ contains the required .rds files
#      (see data/CODEBOOK.md for expected input format)
#   3. Set models_to_run = c("Pooled","Young","Adult","Elderly") until
#      heterogeneity analysis data is available
################################################################################

country_config <- list(

  US = list(
    time_agg_var  = "month",
    fixed_effects = "adm_id^year + adm_id^month",
    rhs_vars      = c(
      "temp_d1_l0", "temp_d2_l0", "temp_d3_l0", "temp_d4_l0",
      "precip_d1_l0", "precip_d2_l0"
    ),
    models_to_run = c("Pooled", "Young", "Adult", "Elderly",
                      "temp", "income", "decade", "baseline"),
    figs          = c(1, 2, 3),
    year_range    = NULL   # NULL = no filter; or c(start, end)
  ),

  EU = list(
    time_agg_var  = "week",
    fixed_effects = "adm_id^year + adm_id^week",
    rhs_vars      = c(
      "temp_d1_l0", "temp_d1_t1", "temp_d1_t2", "temp_d1_t3",
      "temp_d2_l0", "temp_d2_t1", "temp_d2_t2", "temp_d2_t3",
      "temp_d3_l0", "temp_d3_t1", "temp_d3_t2", "temp_d3_t3",
      "temp_d4_l0", "temp_d4_t1", "temp_d4_t2", "temp_d4_t3",
      "precip_d1_l0", "precip_d1_t1", "precip_d1_t2", "precip_d1_t3",
      "precip_d2_l0", "precip_d2_t1", "precip_d2_t2", "precip_d2_t3"
    ),
    models_to_run = c("Pooled", "Young", "Adult", "Elderly",
                      "temp", "income", "decade", "baseline"),
    figs          = c(1, 2, 3),
    year_range    = NULL
  ),

  MEX = list(
    time_agg_var  = "month",
    fixed_effects = "adm_id^year + adm_id^month",
    rhs_vars      = c(
      "temp_d1_l0", "temp_d2_l0", "temp_d3_l0", "temp_d4_l0",
      "precip_d1_l0", "precip_d2_l0"
    ),
    models_to_run = c("Pooled", "Young", "Adult", "Elderly",
                      "temp", "income", "decade", "baseline"),
    figs          = c(1, 2, 3),
    year_range    = NULL
  ),

  # ── NEW COUNTRIES ────────
  #
  # Each new country only needs Figs 1–3 and runs the four core age models.
  # Steps to activate:
  #   1. Uncomment the block
  #   2. Confirm time resolution and available age groups
  #   3. Place processed .rds files in processed/<COUNTRY>/
  #
  # TEMPLATE:
  # ,
  BRA = list(
    time_agg_var  = "month",                        # TODO: confirm
    fixed_effects = "adm_id^year + adm_id^month",  # TODO: confirm
    rhs_vars      = c(                              # TODO: confirm
      "temp_d1_l0", "temp_d2_l0", "temp_d3_l0", "temp_d4_l0",
      "precip_d1_l0", "precip_d2_l0"
    ),
    models_to_run = c("Pooled", "Young", "Adult", "Elderly",
                      "temp", "income", "decade", "baseline"),
    figs          = c(1, 2, 3),
    year_range    = c(1998, 2019)  # subset to available income data years
  ),

  # ── REGIONAL POOLS ────────────────────────────────────────────────────────────

  LA = list(
    time_agg_var  = "month",
    fixed_effects = "adm_id^year + adm_id^month",
    rhs_vars      = c(
      "temp_d1_l0", "temp_d2_l0", "temp_d3_l0", "temp_d4_l0",
      "precip_d1_l0", "precip_d2_l0"
    ),
    models_to_run = c("Pooled"),      # pooled region: Pooled only for now
    figs          = c(1, 3),          # appears in Fig 1 and Fig 3
    year_range    = NULL
  ),

  # ── MULTI-DATASET COUNTRIES ──────────────────────────────────────────────────
  # Entries below use dataset_tag as the config key (not country_name).
  # Each entry carries country_name = "THL" so attribution.R knows the actual
  # country code for file naming and path construction.
  # figs = c() excludes them from auto-loading in attribution_plot.R until
  # the preferred dataset is chosen.

  THL_1997_2025 = list(
    country_name  = "THL",                          # actual country code (for file names / paths)
    time_agg_var  = "month",
    fixed_effects = "adm_id^year + adm_id^month",
    rhs_vars      = c(
      "temp_d1_l0", "temp_d2_l0", "temp_d3_l0", "temp_d4_l0",
      "precip_d1_l0", "precip_d2_l0"
    ),
    models_to_run = c("Pooled"),
    figs          = c(),                             # excluded from auto-loading; activate when ready
    year_range    = NULL
  ),

  THL_2005_2025 = list(
    country_name  = "THL",
    time_agg_var  = "month",
    # year_str_break included as FE; update formula once FE structure is confirmed
    fixed_effects = "adm_id^year + adm_id^month",
    rhs_vars      = c(
      "temp_d1_l0", "temp_d2_l0", "temp_d3_l0", "temp_d4_l0",
      "precip_d1_l0", "precip_d2_l0"
    ),
    models_to_run = c("Pooled"),
    figs          = c(),                             # excluded from auto-loading; activate when ready
    year_range    = NULL
  )

  #
  # CHL = list(
  #   time_agg_var  = "month",
  #   fixed_effects = "adm_id^year + adm_id^month",
  #   rhs_vars      = c(
  #     "temp_d1_l0", "temp_d2_l0", "temp_d3_l0", "temp_d4_l0",
  #     "precip_d1_l0", "precip_d2_l0"
  #   ),
  #   models_to_run = c("Pooled", "Young", "Adult", "Elderly"),
  #   figs          = c(1, 2, 3)
  # ),
  #
  # MNG = list(
  #   time_agg_var  = "month",
  #   fixed_effects = "adm_id^year + adm_id^month",
  #   rhs_vars      = c(
  #     "temp_d1_l0", "temp_d2_l0", "temp_d3_l0", "temp_d4_l0",
  #     "precip_d1_l0", "precip_d2_l0"
  #   ),
  #   models_to_run = c("Pooled", "Young", "Adult", "Elderly"),
  #   figs          = c(1, 2, 3)
  # ),
  #
  # PHL = list(
  #   time_agg_var  = "month",
  #   fixed_effects = "adm_id^year + adm_id^month",
  #   rhs_vars      = c(
  #     "temp_d1_l0", "temp_d2_l0", "temp_d3_l0", "temp_d4_l0",
  #     "precip_d1_l0", "precip_d2_l0"
  #   ),
  #   models_to_run = c("Pooled", "Young", "Adult", "Elderly"),
  #   figs          = c(1, 2, 3)
  # ),
  #
  # THL = list(
  #   time_agg_var  = "month",
  #   fixed_effects = "adm_id^year + adm_id^month",
  #   rhs_vars      = c(
  #     "temp_d1_l0", "temp_d2_l0", "temp_d3_l0", "temp_d4_l0",
  #     "precip_d1_l0", "precip_d2_l0"
  #   ),
  #   models_to_run = c("Pooled", "Young", "Adult", "Elderly"),
  #   figs          = c(1, 2, 3)
  # ),
  #
  # ZAF = list(
  #   time_agg_var  = "month",
  #   fixed_effects = "adm_id^year + adm_id^month",
  #   rhs_vars      = c(
  #     "temp_d1_l0", "temp_d2_l0", "temp_d3_l0", "temp_d4_l0",
  #     "precip_d1_l0", "precip_d2_l0"
  #   ),
  #   models_to_run = c("Pooled", "Young", "Adult", "Elderly"),
  #   figs          = c(1, 2, 3)
  # )
)
