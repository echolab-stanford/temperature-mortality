# responses by cause

# Load packages (fastverse loads many data/table-related packages)
pacman::p_load(fastverse, tidyverse, arrow, fixest)

# Set working directory
setwd('/Users/aw/Dropbox/Research/temp_mortality')

# Set seed for reproducibility
set.seed(42)

# Use multithreading for fixest regressions
fixest::setFixest_nthreads(8)

# Define "not in" operator
`%notin%` <- Negate(`%in%`)

# Tensor product operator for spline-like model matrices
"%c%" <- function(a,b) {
  mgcv::tensor.prod.model.matrix(list(as.matrix(a),as.matrix(b)))
}

# Arrow option to pull vectors by default
options(arrow.pull_as_vector = TRUE)

# Clean theme for plotting
theme_clean <- function(){
  theme_classic() +
    theme(axis.text=element_text(size=10),axis.title=element_text(size=10),
          axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)),
          axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)),
          plot.title=element_text(color="black",size=12,hjust=0.5,face="bold"),
          axis.line = element_line(color = 'black',linewidth  = 0.35),
          axis.ticks = element_line(colour = "black",linewidth = 0.35),
          axis.text.x=element_text(color="black"),
          axis.text.y=element_text(color="black"))
}

# Okabe–Ito colorblind-safe palette
okabe_colors <- c('#D55E00','#0072B2','#009E73','#CC79A7','#000000','#E69F00','#F0E442','#4DA3E5','#555555','#D9D9D9','#93C373','#C42500','#DDAF98','#CBB8D8')

# --- Mexico mortality by cause ---

# Load monthly mortality file, construct adm1, restrict to ≤2019
mexico_mort <- read_parquet('data/mexico/full_monthly.pq') %>% 
  mutate(adm1 = adm2 %/% 1e3) %>% 
  filter(year <= 2019)
# {https://www.dropbox.com/scl/fi/emswcm2vhz26cpmgao44r/full_monthly.pq?rlkey=qt0jlm01yarkogb2bn5q05342&dl=0}

# Population (average per year)
mexico_pops <- mexico_mort %>% 
  group_by(adm2, year) %>% 
  summarize(pop = mean(total_pop)) %>% 
  as.data.table()

# Weather dataset (ERA5 temperature exposure)
mexico_weather <- open_dataset('data/mexico/met/mexico_temp.pq') %>% 
  filter(year >= 1990, year <= 2019) %>% 
  select(temp = order_1, year, month, day, poly_id) %>% 
  mutate(year = as.integer(year), poly_id = as.integer(poly_id)) %>% 
  collect()

# Mapping between grid polygons and administrative units
adm2_mapping <- open_dataset('data/mexico/geometry/adm2_to_order.csv', format = 'csv')

# Merge spatial mapping into weather
mexico_weather %<>% 
  merge(adm2_mapping %>% rename(poly_id = order), by = 'poly_id') %>% 
  select(-poly_id)

# Merge population into weather dataset
mexico_weather %<>% merge(mexico_pops, by = c('adm2','year'), all.x = TRUE)

# Aggregate daily exposures into temperature × year person-days
mexico_exposure <- mexico_weather %>% 
  mutate(temp = round(temp)) %>% 
  group_by(temp, year) %>% 
  summarize(person_days = sum(pop, na.rm = TRUE))

# Prediction temperature vector
pred_temps_vector <- -50:50

# Range for selecting minimum mortality temperature (MMT)
mmt_options <- which(pred_temps_vector <= 32 & pred_temps_vector >= 5)

# Option: use a fixed MMT or estimate it from model
fixed_mmt <- FALSE
mmt_if_fixed <- 0

# Polynomial basis (degree 4)
pred_temps_mat <- poly(pred_temps_vector, degree = 4, raw = TRUE)

# Cause lists and labels for plotting
causes <- c('cardio','cancer','resp','brain','accident','inf_immune','nutrition','misc_birth_blood_bone_gastro_ext','other','kidney','poison','suicide','homicide','sudden')
cause_labels <- c('Cardiovascular','Cancer','Respiratory','Neurodegenerative','Accident','Infections/\nimmune conditions','Nutrition-related','Blood/Bone/\nBirth-related','Other','Renal','Poisoning','Suicide','Homicide','Sudden/unexplained')

# Lookup for facet labels
cause_lookup <- as_labeller(function(x) cause_labels[match(x, causes)])

# Compute regression weights (share of population in each year)
mexico_mort[, weight := total_pop/sum(total_pop), by = 'year']

# Initialize result containers
cause_plot_df_mexico <- list()
deaths_by_cause_year_mexico <- list()

# Loop through causes and estimate temperature–mortality curve
for (this_cause in causes) {
  cat(this_cause, '\n')
  
  # Fixed-effects regression: cause-specific mortality vs temp/precip splines
  cause_naive <- feols(
    as.formula(paste0(this_cause, '/total_pop ~ temp1_l0 + temp2_l0 + temp3_l0 + temp4_l0 + precip1_l0 + precip2_l0 | adm2^year + adm2^month')),
    data = mexico_mort, 
    weights = ~ weight, 
    cluster = 'adm1'
  )
  
  # Extract temperature-related coefficients
  which_temp <- grepl('temp', names(coef(cause_naive)))
  temp_coefs <- coef(cause_naive)[which_temp]
  temp_vcov <- vcov(cause_naive)[which_temp, which_temp]
  
  # Compute MMT (minimum predicted mortality temperature)
  cause_mmt <- pred_temps_vector[mmt_options][which.min((pred_temps_mat %*% temp_coefs)[mmt_options])]
  if (fixed_mmt) {cause_mmt <- mmt_if_fixed}
  
  # Center polynomial at the MMT
  pred_temps_centered <- scale(pred_temps_mat, center = poly(cause_mmt, degree = 4, raw = TRUE), scale = FALSE)
  
  # Store curve for plotting
  cause_plot_df_mexico[[length(cause_plot_df_mexico)+1]] <- 
    data.table(
      cause = this_cause,
      temp = pred_temps_vector,
      mort_rate = c(pred_temps_centered %*% temp_coefs),
      se = sqrt(pmax(0, rowSums(pred_temps_centered * (pred_temps_centered %*% temp_vcov))))
    ) %>% mutate(side = fifelse(temp < cause_mmt, 'cold', 'hot'))
  
  # Compute estimated deaths by year × temperature
  deaths_by_cause_year_mexico[[length(deaths_by_cause_year_mexico)+1]] <- 
    data.table(
      cause = this_cause,
      temp = rep(pred_temps_vector, length(unique(mexico_exposure$year))),
      year = rep(unique(mexico_exposure$year), each = length(pred_temps_vector)),
      mort_rate = rep(c(pred_temps_centered %*% temp_coefs), length(unique(mexico_exposure$year)))
    ) %>%
    merge(mexico_exposure, by = c('temp','year')) %>%
    mutate(deaths = person_days * mort_rate) %>%
    mutate(side = fifelse(temp < cause_mmt, 'cold', 'hot'))
}

# Bind results
cause_plot_df_mexico <- rbindlist(cause_plot_df_mexico)
deaths_by_cause_year_mexico <- rbindlist(deaths_by_cause_year_mexico)

# --- U.S. mortality by cause ---

# Load U.S. mortality data --- data is not public so cannot be shared for replication
us_mort <- read_parquet('data/usa/full_monthly.pq') %>% 
  mutate(state = fips %/% 1e3) %>% 
  as.data.table()

# U.S. population per FIPS-year
us_pops <- us_mort %>% group_by(fips, year) %>% summarize(pop = mean(pop_total))

# U.S. temperature data (ERA5)
us_weather <- open_dataset('data/usa/met/era5/us_temp.pq') %>% 
  filter(year >= 1969, year <= 2019) %>% 
  select(temp = order_1, year, month, day, fips = poly_id) %>% 
  mutate(year = as.integer(year), fips = as.integer(fips)) %>% 
  collect()

# Merge weather with population
us_weather %<>% merge(us_pops, by = c('fips','year'), all.x = TRUE)

# Construct exposure matrix: person-days by temperature × year
us_exposure <- us_weather %>% 
  mutate(temp = round(temp)) %>% 
  group_by(temp, year) %>% 
  summarize(person_days = sum(pop, na.rm = TRUE))

# Initialize lists
cause_plot_df_us <- list()
deaths_by_cause_year_us <- list()

# Loop over causes (same approach as Mexico)
for (this_cause in causes) {
  cat(this_cause, '\n')
  
  # FE regression for the U.S.
  cause_naive <- feols(
    as.formula(paste0(this_cause, '/pop_total ~ temp_d1_l0 + temp_d2_l0 + temp_d3_l0 + temp_d4_l0 + precip_d1_l0 + precip_d2_l0 | fips^year + fips^month')),
    data = us_mort %>% filter(year < 2020),
    weights = ~ pop_total,
    cluster = 'state'
  )
  
  # Extract coefficients
  which_temp <- grepl('temp', names(coef(cause_naive)))
  temp_coefs <- coef(cause_naive)[which_temp]
  temp_vcov <- vcov(cause_naive)[which_temp, which_temp]
  
  # Compute MMT
  cause_mmt <- pred_temps_vector[mmt_options][which.min((pred_temps_mat %*% temp_coefs)[mmt_options])]
  if (fixed_mmt) {cause_mmt <- mmt_if_fixed}
  
  # Center polynomial
  pred_temps_centered <- scale(pred_temps_mat, center = poly(cause_mmt, degree = 4, raw = TRUE), scale = FALSE)
  
  # Store prediction curves
  cause_plot_df_us[[length(cause_plot_df_us)+1]] <- 
    data.table(
      cause = this_cause,
      temp = pred_temps_vector,
      mort_rate = c(pred_temps_centered %*% temp_coefs),
      se = sqrt(pmax(0, rowSums(pred_temps_centered * (pred_temps_centered %*% temp_vcov))))
    ) %>% mutate(side = fifelse(temp < cause_mmt, 'cold', 'hot'))
  
  # Compute deaths by year × temperature
  deaths_by_cause_year_us[[length(deaths_by_cause_year_us)+1]] <- 
    data.table(
      cause = this_cause,
      temp = rep(pred_temps_vector, length(unique(us_exposure$year))),
      year = rep(unique(us_exposure$year), each = length(pred_temps_vector)),
      mort_rate = rep(c(pred_temps_centered %*% temp_coefs), length(unique(us_exposure$year)))
    ) %>%
    merge(us_exposure, by = c('temp','year')) %>%
    mutate(deaths = person_days * mort_rate) %>%
    mutate(side = fifelse(temp < cause_mmt, 'cold', 'hot'))
}

# Bind U.S. results
cause_plot_df_us <- rbindlist(cause_plot_df_us)
deaths_by_cause_year_us <- rbindlist(deaths_by_cause_year_us)

# Combine Mexico and U.S. data for plotting
cause_plot_df <- cause_plot_df_mexico %>% 
  mutate(region = 'Mexico') %>% 
  rbind(cause_plot_df_us %>% mutate(region = 'USA'))

# Multiplier for presentation scale (per 100k)
multiplier <- 1e5

# Plot all causes (Mexico + USA) with confidence bands
cause_multiples <- ggplot(
  cause_plot_df %>% 
    filter((region == 'Mexico' & temp > 5 & temp < 40) | 
             (region == 'USA' & temp > -10 & temp < 40)) %>% 
    mutate(cause = factor(cause, levels = c('cardio','cancer','resp','brain','accident','inf_immune','nutrition','misc_birth_blood_bone_gastro_ext','other','kidney','poison','suicide','homicide','sudden'))),
  aes(
    x = temp, 
    y = mort_rate * multiplier, 
    color = region, 
    fill = region,
    ymin = (mort_rate - 1.96 * se) * multiplier,
    ymax = (mort_rate + 1.96 * se) * multiplier
  )
) +
  geom_hline(yintercept = 0, color = 'black', alpha = 0.25) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_line() +
  facet_wrap(~cause, nrow = 3, scales = 'free', labeller = labeller(cause = cause_lookup)) +
  theme_clean() +
  ylab('Change in mortality per 100,000') +
  xlab('Temperature (ºC)') +
  scale_color_manual(values = okabe_colors, name = 'Region') +
  scale_fill_manual(values = okabe_colors, name = 'Region') +
  coord_cartesian(ylim = c(0, NA)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black", face = "bold", size = 8))

# Print figure
cause_multiples

# Save plot
ggsave(cause_multiples, file = 'output/figures/presentation/cause.pdf', height = 6, width = 10, units = 'in')