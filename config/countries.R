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
    figs          = c(1, 2, 3)
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
    figs          = c(1, 2, 3)
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
    figs          = c(1, 2, 3)
  )

  # ── NEW COUNTRIES (placeholder — fill in once processed data is ready) ────────
  #
  # Each new country only needs Figs 1–3 and runs the four core age models.
  # Steps to activate:
  #   1. Uncomment the block
  #   2. Confirm time resolution and available age groups
  #   3. Place processed .rds files in processed/<COUNTRY>/
  #
  # TEMPLATE:
  # ,
  # BRA = list(
  #   time_agg_var  = "month",                        # TODO: confirm
  #   fixed_effects = "adm_id^year + adm_id^month",  # TODO: confirm
  #   rhs_vars      = c(                              # TODO: confirm
  #     "temp_d1_l0", "temp_d2_l0", "temp_d3_l0", "temp_d4_l0",
  #     "precip_d1_l0", "precip_d2_l0"
  #   ),
  #   models_to_run = c("Pooled", "Young", "Adult", "Elderly"),
  #   figs          = c(1, 2, 3)
  # ),
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
