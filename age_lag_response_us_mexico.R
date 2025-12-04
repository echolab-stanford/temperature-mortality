# age and lag response plot (U.S. and Mexico)

# Load packages (fastverse pulls in data.table, etc.; tidyverse for data wrangling; arrow for parquet; pammtools/dlnm for lag structures;
# fixest for regressions; tictoc for timing; sf for spatial; ggnewscale for multiple scales in ggplot)
pacman::p_load(fastverse, tidyverse, arrow, pammtools, fixest, tictoc, sf, ggnewscale)

# Set working directory
setwd('/Users/aw/Dropbox/Research/temp_mortality')

# Set seed for reproducibility
set.seed(42)

# Use multiple threads for fixest
fixest::setFixest_nthreads(8)

# Negated %in% operator
`%notin%` <- Negate(`%in%`)

# Arrow: pull vectors as R vectors instead of Arrow arrays
options(arrow.pull_as_vector = TRUE)

# Custom ggplot theme for clean plots
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

# Color scales for lag and age plots
lag_color_scale <- c('#D55E00','#0072B2','#009E73','#CC79A7','#000000','#E69F00','#F0E442','#4DA3E5','#555555','#D9D9D9','#93C373','#C42500','#DDAF98','#CBB8D8')
age_color_scale1 <- viridis::viridis(7)
age_color_scale2 <- c(age_color_scale1[1:2],'red',age_color_scale1[3:7])

## us lag response ------------------------------------------------------------------

# Read single-age adjustment to convert single-year ages into weights
age_adj <- fread('data/usa/pop/singleage_adjustment.csv', skip = 1)
# {https://www.dropbox.com/scl/fi/74053cyzuoha47j3dfa0k/singleage_adjustment.csv?rlkey=266eo3eugc0chixzfurmk7a3n&dl=0}

# Drop header/footer rows
age_adj <- age_adj[2:(nrow(age_adj)-4)]

# Keep only age and count columns
age_adj %<>% select(age = V1, count = V2)

# Parse counts as integer (remove commas)
age_adj %<>% mutate(count = as.integer(gsub(',','',count)))

# Explicitly set ages 0–90
age_adj$age <- 0:90

# Create age weights based on population share
age_adj$weight <- age_adj$count/sum(age_adj$count)

# Keep only age and weight
age_adj %<>% select(-count)

# Read county-year population totals (single-age panel aggregated)
pop_panel <- read_parquet('data/usa/pop/pop_panel_singleage.pq')
# {https://www.dropbox.com/scl/fi/zdes5va6xlmpdbr6dk9n6/pop_panel_singleage.pq?rlkey=f4len1nozbshac0q3bzg168ps&dl=0}
pop_panel <- pop_panel[, .(pop_total = sum(pop)), by = list(fips, year)]

# Read ERA5 temperature and precipitation datasets
temp <- open_dataset('data/usa/met/era5/us_temp.pq', format = 'parquet')
# {https://www.dropbox.com/scl/fi/4zjz4gos7sd8bbul1mw2y/us_temp.pq?rlkey=k2xv5hzwa6xslnmcyfgnnw0wn&dl=0}
precip <- open_dataset('data/usa/met/era5/us_precip.pq', format = 'parquet')
# {https://www.dropbox.com/scl/fi/ppmm5er5k86c1l7mq14tt/us_precip.pq?rlkey=3wndfp4w202mksbjlousik3c5&dl=0}

# Build daily temperature data with 4 polynomial orders
temp %<>% mutate(date = ymd(paste0(year, '-', month, '-', day))) %>% select(fips = poly_id, date, temp1 = order_1, temp2 = order_2, temp3 = order_3, temp4 = order_4)

# Build daily precipitation data with 2 polynomial orders
precip %<>% mutate(date = ymd(paste0(year, '-', month, '-', day))) %>% select(fips = poly_id, date, precip1 = order_1, precip2 = order_2)

# Merge temp and precip (older path; superseded by full.pq below but kept)
weather <- merge(temp %>% mutate(fips = as.integer(fips)), precip %>% mutate(fips = as.integer(fips)), by = c('fips','date')) %>%
  mutate(year = year(date)) %>% 
  filter(year > 1970, year < 1989) %>% 
  rename(temp_d1 = temp1, temp_d2 = temp2, temp_d3 = temp3, temp_d4 = temp4, precip_d1 = precip1, precip_d2 = precip2) %>% 
  group_by(fips, date) %>% 
  summarize(across(temp_d1:precip_d2, ~ mean(.x))) %>%
  as.data.table()

# Read U.S. individual-level mortality (no multiple causes), apply geographic harmonization,
# and restrict to contiguous states and valid FIPS --- data is private and thus not available for replication
mort <- arrow::open_dataset('-----', format = 'parquet') %>%
  filter(year > 1971, year < 1989) %>%
  mutate(fips_res =
           case_when(fips_res %in% c(4027, 4910) ~ 4012, #boundary changes (Census), SEER restriction
                     fips_res %in% c(8005, 8013, 8014, 8031, 8059, 8123, 8911, 8912, 8913, 8914) ~ 8001, #boundary changes (Census), SEER restriction
                     fips_res %in% c(12086) ~ 12025, #name change (Census)
                     fips_res %in% c(24033) ~ 24031, #boundary changes (Census)
                     fips_res %in% c(29193) ~ 29186, #SEER restriction
                     fips_res %in% c(30067, 30113) ~ 30031, #boundary changes (Census), SEER restriction
                     fips_res %in% c(35061, 35910) ~ 35006, #boundary changes (Census), SEER restriction
                     fips_res %in% c(36047, 36061, 36081, 36085, 36910) ~ 36005, #SEER restriction
                     fips_res %in% c(46131) ~ 46071, #boundary changes (Census)
                     fips_res %in% c(46113) ~ 46102, #name change (Census), SEER restriction
                     fips_res %in% c(51540) ~ 51003, #boundary changes (Census)
                     fips_res %in% c(51560, 51580, 51916) ~ 51005, #boundary changes (Census), SEER restriction
                     fips_res %in% c(51059, 51510, 51600, 51610, 51918) ~ 51013, #boundary changes (Census), SEER restriction
                     fips_res %in% c(51790, 51820) ~ 51015, #boundary changes (Census)
                     fips_res %in% c(51031, 51515, 51680, 51917) ~ 51019, #boundary changes (Census), SEER restriction
                     fips_res %in% c(51149, 51570, 51670, 51730) ~ 51053, #boundary changes (Census)
                     fips_res %in% c(51840) ~ 51069, #boundary changes (Census)
                     fips_res %in% c(51595) ~ 51081, #boundary changes (Census)
                     fips_res %in% c(51780) ~ 51083, #boundary changes (Census), SEER restriction
                     fips_res %in% c(51690, 51913) ~ 51089, #SEER restriction
                     fips_res %in% c(51199, 51700, 51735, 51830, 51911) ~ 51095, #boundary changes (Census), SEER restriction
                     fips_res %in% c(51750) ~ 51121, #boundary changes (Census)
                     fips_res %in% c(51800) ~ 51123, #boundary changes (Census)
                     fips_res %in% c(51590) ~ 51143, #boundary changes (Census)
                     fips_res %in% c(51683, 51685, 51910) ~ 51153, #boundary changes (Census), SEER restriction
                     fips_res %in% c(51770, 51775) ~ 51161, #boundary changes (Census)
                     fips_res %in% c(51530, 51678) ~ 51163, #boundary changes (Census)
                     fips_res %in% c(51660, 51914) ~ 51165, #boundary changes (Census), SEER restriction
                     fips_res %in% c(51620) ~ 51175, #boundary changes (Census)
                     fips_res %in% c(51179, 51630, 51915) ~ 51177, #boundary changes (Census), SEER restriction
                     fips_res %in% c(51520) ~ 51191, #boundary changes (Census)
                     fips_res %in% c(55115) ~ 55078, #SEER restriction
                     fips_res %in% c(51640) ~ 51035, #BLS restriction
                     fips_res %in% c(51720) ~ 51195, #BLS restriction
                     fips_res %in% c(55115) ~ 55078, #BLS restriction
                     fips_res == 32025 ~ 32510, #NCHS to FIPS mapping issue
                     fips_res == 16089 ~ 16021, #NCHS to FIPS mapping issue
                     fips_res == 56047 ~ 56029, #NCHS to FIPS mapping issue
                     TRUE ~ fips_res),
         fips_occur =
           case_when(fips_occur %in% c(4027, 4910) ~ 4012, #boundary changes (Census), SEER restriction
                     fips_occur %in% c(8005, 8013, 8014, 8031, 8059, 8123, 8911, 8912, 8913, 8914) ~ 8001, #boundary changes (Census), SEER restriction
                     fips_occur %in% c(12086) ~ 12025, #name change (Census)
                     fips_occur %in% c(24033) ~ 24031, #boundary changes (Census)
                     fips_occur %in% c(29193) ~ 29186, #SEER restriction
                     fips_occur %in% c(30067, 30113) ~ 30031, #boundary changes (Census), SEER restriction
                     fips_occur %in% c(35061, 35910) ~ 35006, #boundary changes (Census), SEER restriction
                     fips_occur %in% c(36047, 36061, 36081, 36085, 36910) ~ 36005, #SEER restriction
                     fips_occur %in% c(46131) ~ 46071, #boundary changes (Census)
                     fips_occur %in% c(46113) ~ 46102, #name change (Census), SEER restriction
                     fips_occur %in% c(51540) ~ 51003, #boundary changes (Census)
                     fips_occur %in% c(51560, 51580, 51916) ~ 51005, #boundary changes (Census), SEER restriction
                     fips_occur %in% c(51059, 51510, 51600, 51610, 51918) ~ 51013, #boundary changes (Census), SEER restriction
                     fips_occur %in% c(51790, 51820) ~ 51015, #boundary changes (Census)
                     fips_occur %in% c(51031, 51515, 51680, 51917) ~ 51019, #boundary changes (Census), SEER restriction
                     fips_occur %in% c(51149, 51570, 51670, 51730) ~ 51053, #boundary changes (Census)
                     fips_occur %in% c(51840) ~ 51069, #boundary changes (Census)
                     fips_occur %in% c(51595) ~ 51081, #boundary changes (Census)
                     fips_occur %in% c(51780) ~ 51083, #boundary changes (Census), SEER restriction
                     fips_occur %in% c(51690, 51913) ~ 51089, #SEER restriction
                     fips_occur %in% c(51199, 51700, 51735, 51830, 51911) ~ 51095, #boundary changes (Census), SEER restriction
                     fips_occur %in% c(51750) ~ 51121, #boundary changes (Census)
                     fips_occur %in% c(51800) ~ 51123, #boundary changes (Census)
                     fips_occur %in% c(51590) ~ 51143, #boundary changes (Census)
                     fips_occur %in% c(51683, 51685, 51910) ~ 51153, #boundary changes (Census), SEER restriction
                     fips_occur %in% c(51770, 51775) ~ 51161, #boundary changes (Census)
                     fips_occur %in% c(51530, 51678) ~ 51163, #boundary changes (Census)
                     fips_occur %in% c(51660, 51914) ~ 51165, #boundary changes (Census), SEER restriction
                     fips_occur %in% c(51620) ~ 51175, #boundary changes (Census)
                     fips_occur %in% c(51179, 51630, 51915) ~ 51177, #boundary changes (Census), SEER restriction
                     fips_occur %in% c(51520) ~ 51191, #boundary changes (Census)
                     fips_occur %in% c(55115) ~ 55078, #SEER restriction
                     fips_occur %in% c(51640) ~ 51035, #BLS restriction
                     fips_occur %in% c(51720) ~ 51195, #BLS restriction
                     fips_occur %in% c(55115) ~ 55078, #BLS restriction
                     fips_occur == 32025 ~ 32510, #NCHS to FIPS mapping issue
                     fips_occur == 16089 ~ 16021, #NCHS to FIPS mapping issue
                     fips_occur == 56047 ~ 56029, #NCHS to FIPS mapping issue
                     TRUE ~ fips_occur)) %>%
  rename(fips = fips_res) %>% 
  filter(fips != 0,                            # drop unknown FIPS
         fips %/% 1e3 != 2,                    # drop AK
         fips %/% 1e3 != 15,                   # drop HI
         fips != 13999)                        # drop special federal region code in Georgia

# Aggregate to daily county-level total deaths
mort_total <- mort %>% group_by(fips, year, month, day) %>% summarize(deaths_total = n()) %>% collect() %>% as.data.table()

# Attach a date field and keep minimal columns
mort_total %<>% mutate(date = ymd(paste0(year, '-', month, '-', day))) %>% select(fips, date, deaths_total)
mort_total <- mort_total[!is.na(date)]

# Build a full county-day panel for study period
us_df <- CJ(fips = unique(mort_total$fips), date = seq.Date(from = as.Date('1972-01-01'), to = as.Date('1988-12-31'), by = 'day'))

# Merge deaths (fill missing with zero)
us_df %<>% merge(mort_total, by = c('fips','date'), all.x = T)
us_df[, deaths_total := ifelse(is.na(deaths_total), 0, deaths_total)]

# Extract year from date
us_df[, year := year(date)]

# Merge county-year population and daily weather into panel
us_df %<>% merge(pop_panel, by = c('fips','year'), all.x = T)
us_df %<>% merge(weather, by = c('fips','date'), all.x = T)

# Recompute year and month just to ensure consistency
us_df[, `:=`(year = year(date), month = month(date))]

# Compute monthly mean temp_d1 by county-year-month
fips_monthyear_means <- us_df[, .(monthmean = mean(temp_d1)), by = list(fips, month, year)]

# Create a time index (month index) across entire sample
fips_monthyear_means[, index := year*12 + month]
fips_monthyear_means[, index := index - min(index) + 1]

# Estimate linear temperature trend over time within each county-month combination
per_county <- fips_monthyear_means %>%
  group_by(fips, month) %>%
  do({
    fit <- lm(monthmean ~ 1 + index, data = .)
    tidied <- broom::tidy(fit)
    slope_row <- tidied[tidied$term == "index",]
    data.frame(
      beta_hat   = slope_row$estimate,
      sigma2_hat = slope_row$std.error^2
    )
  }) %>%
  ungroup()

# Maximum lag for U.S. lag-response
max_lag <- 60

# Build lag basis for distributed lag non-linear model
predlag <- 0:max_lag
lag_df <- 5
basislag <- do.call(dlnm::onebasis,c(list(x = predlag), list(fun = "bs", knots = dlnm::equalknots(0:max_lag, nk = lag_df), intercept = TRUE)))

# Polynomial order for temperature exposure
temp_poly_order <- 4

# Allocate matrix for temperature-lag basis predictors
temp_basis <- matrix(0, nrow=nrow(us_df), ncol=temp_poly_order*ncol(basislag))

# Build cross-basis: polynomial in temp × spline in lag, by county
for (v in seq(length=temp_poly_order)) {
  print(paste0('poly ', v))
  # Build lagged temp matrix (0:max_lag) grouped by fips
  mat <- as.matrix(flag(us_df[[paste0('temp_d',v)]], 0:max_lag, g=us_df$fips))
  for (l in seq(length=ncol(basislag))) {
    print(paste0('lag ', l))
    # For each lag-basis column, multiply lagged temp by that basis
    temp_basis[,ncol(basislag)*(v-1)+l] <- mat %*% (basislag[,l])
  }
}
rm(mat)

# Convert to data.table and name columns
temp_basis <- as.data.table(temp_basis)
temp_basis_colnames <- c(t(outer(paste0('temp',1:temp_poly_order, '_l'), 1:ncol(basislag), FUN = paste0)))
colnames(temp_basis) <- temp_basis_colnames

# Attach identifiers
temp_basis[, `:=`(fips = us_df$fips, date = us_df$date)]

# Polynomial order for precipitation exposure
precip_poly_order <- 2

# Allocate matrix for precipitation-lag basis predictors
precip_basis <- matrix(0, nrow=nrow(us_df), ncol=precip_poly_order*ncol(basislag))

# Build precipitation cross-basis
for (v in seq(length=precip_poly_order)) {
  print(paste0('poly ', v))
  mat <- as.matrix(flag(us_df[[paste0('precip_d',v)]], 0:max_lag, g=us_df$fips))
  for (l in seq(length=ncol(basislag))) {
    print(paste0('lag ', l))
    precip_basis[,ncol(basislag)*(v-1)+l] <- mat %*% (basislag[,l])
  }
}
rm(mat)

# Convert to data.table and name columns
precip_basis <- as.data.table(precip_basis)
precip_basis_colnames <- c(t(outer(paste0('precip',1:precip_poly_order, '_l'), 1:ncol(basislag), FUN = paste0)))
colnames(precip_basis) <- precip_basis_colnames

# Merge temp and precip lag-basis into a single dataset
precip_basis[, `:=`(fips = us_df$fips, date = us_df$date)]
weather_basis <- temp_basis %>% merge(precip_basis, by = c('fips','date'))

# Free up memory
rm(temp_basis)
rm(precip_basis)
gc()

# Drop raw weather vars (temp_d*, precip_d*) and merge basis terms into main panel
us_df %<>% select(-starts_with('temp'), -starts_with('precip')) %>% merge(weather_basis, by = c('fips','date'))

# Add state FIPS (first two digits)
us_df[, state := substr(fips, 1, 2)]

# More cleanup
rm(weather_basis)
gc()

# Create week and month identifiers
us_df[, week := isoweek(date)]
us_df[, month := month(date)]

# Estimate U.S. lag-response model with high-dimensional fixed effects and population weights
tictoc::tic()
model <- fixest::feols(
  as.formula(paste0(
    "deaths_total/pop_total ~ ",
    paste(temp_basis_colnames, collapse = " + "),
    " + ", paste(precip_basis_colnames, collapse = " + "),
    " | fips^year + fips^week + date"
  )),
  weights = ~ pop_total,
  data = us_df,
  se = "standard"
)

## Select temperature coefficients and clustered variance-covariance matrix
indices <- grepl("temp", names(coef(model)))
coefs <- coef(model)[indices]
vcovmat <- vcov(model, cluster = "state")[indices,indices] #add clustering at state level
tictoc::toc()

# Prediction grid over temperature and lag
predvar <- -50:50
predlag <- 0:max_lag

# Build polynomial temperature basis for predictions
predtemps <- poly(rep(predvar, length(predlag)), degree = temp_poly_order, raw = TRUE)

# Build lag basis at prediction grid (logknots for smoother tails)
basispredlag <- do.call(dlnm::onebasis,c(list(x = rep(predlag, each = length(predvar))), list(fun = "bs", knots = dlnm::logknots(0:max_lag, nk = lag_df), intercept = TRUE)))

# Tensor product (temp poly × lag spline)
tensor <- predtemps%c%basispredlag

# Initialize cumulative effect and standard error matrices
xpred <- 0
cumfit <- cumse <- indfit <- matrix(0,length(predvar),length(predlag))

# Compute forward cumulative effect over lags (0 → max_lag)
for (i in seq(length(predlag))) {
  ind <- seq(length(predvar)) + length(predvar)*(i-1)
  xpred <- xpred + tensor[ind,,drop=FALSE]
  cumfit[,i] <- xpred%*%coefs
  cumse[,i] <- sqrt(pmax(0,rowSums(xpred*(xpred%*%vcovmat))))
}
cumfit <- c(cumfit)
cumse <- c(cumse)

# Store prediction surface for U.S. lag-response
pred_df <- data.table(pred = cumfit,
                      lag = rep(0:max_lag, each = length(predvar)),
                      temp = rep(predvar, length(0:max_lag)),
                      se = cumse)

# Identify MMT (minimum mortality temperature) over plausible positive range
mmt_options <- which(predvar < 30 & predvar > 0)
mmt <- pred_df[lag == max_lag]$temp[mmt_options][which.min(pred_df[lag == max_lag]$pred[mmt_options])]

# Center polynomial basis at MMT
predtemps_center <- scale(predtemps, center = poly(mmt, degree = temp_poly_order, raw = T), scale = F)
tensor <- predtemps_center%c%basispredlag

# Recompute cumulative effects after centering at MMT
xpred <- 0
cumfit <- cumse <- matrix(0,length(predvar),length(predlag))
for (i in seq(length(predlag))) {
  ind <- seq(length(predvar)) + length(predvar)*(i-1)
  xpred <- xpred + tensor[ind,,drop=FALSE]
  cumfit[,i] <- xpred%*%coefs
  cumse[,i] <- sqrt(pmax(0,rowSums(xpred*(xpred%*%vcovmat))))
}
cumfit <- c(cumfit)
cumse <- c(cumse)

# Final prediction surface (centered)
pred_df <- data.table(pred = cumfit,
                      lag = rep(0:max_lag, each = length(predvar)),
                      temp = rep(predvar, length(0:max_lag)),
                      se = cumse)

# Save U.S. lag-response predictions
write_parquet(pred_df, sink = 'output/usa/lag_response.pq')

## mexico lag response -------------------------------------------------------------

# Read Mexico daily deaths/population and construct total deaths/pop, plus weights
df <- arrow::open_dataset("data/mexico/daily_death_pop_income.pq", format = "parquet") %>% 
  rename(date = date.x) %>% 
  select(-date.y) %>% 
  mutate(deaths = young + adult + elderly, 
         pop = young_pop + adult_pop + elderly_pop, 
         week = isoweek(date), 
         adm1 = adm2 %/% 1e3) %>% 
  collect()
# {https://www.dropbox.com/scl/fi/szfowpm04zjmizvw866yq/daily_death_pop_income.pq?rlkey=h5bjwtmf7kew0zj0f5by38ic8&dl=0}

# Sort by geography and date
df <- df[order(adm2,date)]

# Build daily weights by share of national population
df[, weight := pop / sum(pop), by = date]

# Maximum lag for Mexico lag-response
max_lag <- 28

# Read daily Mexican temperature and precipitation grids
temp <- open_dataset('data/mexico/met/mexico_temp.pq', format = 'parquet')
# {https://www.dropbox.com/scl/fi/hv37en54tbyhgbaogfzuc/mexico_temp.pq?rlkey=7nul80l3auua579j3tmc4hq2v&dl=0}
precip <- open_dataset('data/mexico/met/mexico_precip.pq', format = 'parquet')
# {https://www.dropbox.com/scl/fi/g5ikw804xf4l2mtxcgzr9/mexico_precip.pq?rlkey=j2z2y136k5b2d0j0j33ballmw&dl=0}

# Mapping from polygon order to ADM2 codes
adm2_mapping <- open_dataset('data/mexico/geometry/adm2_to_order.csv', format = 'csv')
# {https://www.dropbox.com/scl/fi/06my16kddud0pg9pb1jag/adm2_to_order.csv?rlkey=xmxfo2ymhp45tw2r0j3gt0vca&dl=0}

# Merge ADM2 mapping and construct daily temperature series (poly 1–4)
temp %<>% merge(adm2_mapping, by.x = 'poly_id', by.y = 'order') %>% 
  mutate(date = ymd(paste0(year, '-', month, '-', day))) %>% 
  dplyr::select(adm2, date, temp1 = order_1, temp2 = order_2, temp3 = order_3, temp4 = order_4)

# Merge ADM2 mapping and construct daily precipitation series (poly 1–2)
precip %<>% merge(adm2_mapping, by.x = 'poly_id', by.y = 'order') %>% 
  mutate(date = ymd(paste0(year, '-', month, '-', day))) %>% 
  dplyr::select(adm2, date, precip1 = order_1, precip2 = order_2)

# Merge daily temperature and precipitation for Mexico
weather <- merge(temp, precip, by = c('adm2','date'))

# Restrict to modern period for Mexico
weather %<>% filter(date >= '1997-01-01', date <= '2023-12-31')

# Standardize ADM2 codes as 6-character strings
weather[, adm2 := str_pad(as.character(adm2), 6, 'left', '0')]

# Merge weather with deaths/population (right join on weather: include days with weather but possibly missing deaths)
df %<>% merge(weather, by = c('adm2','date'), all.y = T)

# Build lag basis for Mexico (logknots)
predlag <- 0:max_lag
lag_df <- 4
basislag <- do.call(dlnm::onebasis,c(list(x = predlag), list(fun = "bs", knots = dlnm::logknots(0:max_lag, nk = lag_df), intercept = TRUE)))

# Cross-basis for temperature in Mexico
temp_poly_order <- 4
temp_basis <- matrix(0, nrow=nrow(df), ncol=temp_poly_order*ncol(basislag))
for (v in seq(length=temp_poly_order)) {
  print(paste0('poly ', v))
  mat <- as.matrix(flag(df[[paste0('temp',v)]], 0:max_lag, g=df$adm2))
  for (l in seq(length=ncol(basislag))) {
    print(paste0('lag ', l))
    temp_basis[,ncol(basislag)*(v-1)+l] <- mat %*% (basislag[,l])
  }
}
rm(mat)
temp_basis <- as.data.table(temp_basis)
temp_basis_colnames <- c(t(outer(paste0('temp',1:temp_poly_order, '_l'), 1:ncol(basislag), FUN = paste0)))
colnames(temp_basis) <- temp_basis_colnames
temp_basis[, `:=`(adm2 = df$adm2, date = df$date)]

# Cross-basis for precipitation in Mexico
precip_poly_order <- 2
precip_basis <- matrix(0, nrow=nrow(df), ncol=precip_poly_order*ncol(basislag))
for (v in seq(length=precip_poly_order)) {
  print(paste0('poly ', v))
  mat <- as.matrix(flag(df[[paste0('precip',v)]], 0:max_lag, g=df$adm2))
  for (l in seq(length=ncol(basislag))) {
    print(paste0('lag ', l))
    precip_basis[,ncol(basislag)*(v-1)+l] <- mat %*% (basislag[,l])
  }
}
rm(mat)
precip_basis <- as.data.table(precip_basis)
precip_basis_colnames <- c(t(outer(paste0('precip',1:precip_poly_order, '_l'), 1:ncol(basislag), FUN = paste0)))
colnames(precip_basis) <- precip_basis_colnames
precip_basis[, `:=`(adm2 = df$adm2, date = df$date)]

# For Mexico: combine temp and precip basis (dropping geographic identifiers from temp_basis and re-attaching later)
weather_basis <- temp_basis %>% select(-adm2, -date) %>% cbind(precip_basis)

# Free memory
rm(temp_basis)
rm(precip_basis)
gc()

# Construct analysis dataset for Mexico and merge in lag-basis terms
df %<>% select(adm1, adm2, date, week, month, year, weight, deaths, pop) %>% merge(weather_basis, by = c('adm2','date'))

# More cleanup
rm(weather_basis)
gc()

# Estimate Mexico lag-response model with fixed effects and population weights
tictoc::tic()
model <- fixest::feols(
  as.formula(paste0(
    "deaths/pop ~ ",
    paste(temp_basis_colnames, collapse = " + "),
    " + ", paste(precip_basis_colnames, collapse = " + "),
    " | adm2^year + adm1^week + date"
  )),
  weights = ~ weight,
  data = df,
  se = "standard"
)

## Extract temperature coefficients and cluster-robust vcov at ADM1 (state) level
indices <- grepl("temp", names(coef(model)))
coefs <- coef(model)[indices]
vcovmat <- vcov(model, cluster = "adm1")[indices,indices]
tictoc::toc()

# Prediction grid for Mexico
predvar <- -20:40
predlag <- 0:max_lag

# Build polynomial temp basis
predtemps <- poly(rep(predvar, length(predlag)), degree = temp_poly_order, raw = TRUE)

# Lag spline basis with logknots
basispredlag <- do.call(dlnm::onebasis,c(list(x = rep(predlag, each = length(predvar))), list(fun = "bs", knots = dlnm::logknots(0:max_lag, nk = lag_df), intercept = TRUE)))

# Tensor product for Mexico
tensor <- predtemps%c%basispredlag

# Compute cumulative effects over lags
xpred <- 0
cumfit <- cumse <- matrix(0,length(predvar),length(predlag))
for (i in seq(length(predlag))) {
  ind <- seq(length(predvar)) + length(predvar)*(i-1)
  xpred <- xpred + tensor[ind,,drop=FALSE]
  cumfit[,i] <- xpred%*%coefs
  cumse[,i] <- sqrt(pmax(0,rowSums(xpred*(xpred%*%vcovmat))))
}
cumfit <- c(cumfit)
cumse <- c(cumse)

# Store prediction surface for Mexico
pred_df <- data.table(pred = cumfit,
                      lag = rep(0:max_lag, each = length(predvar)),
                      temp = rep(predvar, length(0:max_lag)),
                      se = cumse)

# Identify MMT for Mexico
mmt_options <- which(predvar < 30 & predvar > 0)
mmt <- pred_df[lag == max_lag]$temp[mmt_options][which.min(pred_df[lag == max_lag]$pred[mmt_options])]

# Center temp polynomials at MMT
predtemps_center <- scale(predtemps, center = poly(mmt, degree = temp_poly_order, raw = T), scale = F)
tensor <- predtemps_center%c%basispredlag

# Recompute cumulative response centered at MMT
xpred <- 0
cumfit <- cumse <- matrix(0,length(predvar),length(predlag))
for (i in seq(length(predlag))) {
  ind <- seq(length(predvar)) + length(predvar)*(i-1)
  xpred <- xpred + tensor[ind,,drop=FALSE]
  cumfit[,i] <- xpred%*%coefs
  cumse[,i] <- sqrt(pmax(0,rowSums(xpred*(xpred%*%vcovmat))))
}
cumfit <- c(cumfit)
cumse <- c(cumse)

# Final prediction surface (centered) for Mexico
pred_df <- data.table(pred = cumfit,
                      lag = rep(0:max_lag, each = length(predvar)),
                      temp = rep(predvar, length(0:max_lag)),
                      se = cumse)

# Save Mexico lag-response predictions
write_parquet(pred_df, sink = 'output/mexico/lag_response.pq')

## us response by age --------------------------------------------------------------

# Read U.S. population panel by age group
pop_panel <- read_parquet('data/usa/pop/pop_panel_singleage.pq')

# Total county-year population
total_pop <- pop_panel[, .(pop_total = sum(pop)), by = list(fips, year)]

# Reshape population from long (age_group) to wide (pop_<age_group>)
pop_panel <- pop_panel %>% 
  select(-state) %>% 
  mutate(age_group = paste0('pop_',age_group)) %>% 
  pivot_wider(id_cols = c('fips','year'), names_from = 'age_group', values_from = 'pop')

# Attach total population
pop_panel %<>% merge(total_pop, by = c('fips','year'))

# Read full monthly panel (deaths and population by age group, plus weather) --- mortality
# data is private and thus not available for replication
month_panel <- read_parquet('----') %>% 
  mutate(state = fips %/% 1e3) %>% 
  mutate(
    # Construct broader age groups for some analyses
    deaths_5_24 = deaths_5_14 + deaths_15_24, pop_5_24 = pop_5_14 + pop_15_24,
    deaths_25_44 = deaths_25_34 + deaths_35_44, pop_25_44 = pop_25_34 + pop_35_44,
    deaths_5_44 = deaths_5_24 + deaths_25_44, pop_5_44 = pop_5_24 + pop_25_44,
    deaths_5_54 = deaths_5_44 + deaths_45_54, pop_5_54 = pop_5_44 + pop_45_54,
    deaths_5_24_noviolence = deaths_5_14_noviolence + deaths_15_24_noviolence,
    deaths_25_44_noviolence = deaths_25_34_noviolence + deaths_35_44_noviolence,
    deaths_5_44_noviolence = deaths_5_24_noviolence + deaths_25_44_noviolence,
    deaths_5_54_noviolence = deaths_5_44_noviolence + deaths_45_54_noviolence
  )

# Fine age groups and labels for U.S. age-response
age_groups <- c('lt1', '1_4', '5_14', '15_24', '25_34', '35_44', '5_44', '45_54', '55_64', '65_74', '75_84', 'gt85')
age_group_labels <- c('Under 1', '1 to 4', '5 to 14', '15 to 24', '25 to 34', '35 to 44', '5 to 44', '45 to 54', '55 to 64', '65 to 74', '75 to 84', 'Above 84')

# Temperature prediction grid for age-response
pred_temps_vector <- -50:50
pred_temps <- poly(pred_temps_vector, degree = 4, raw = T)

# List to collect age-specific prediction datasets
pred_df <- list()

# Loop over age groups, including all causes
for (this_age_group in age_groups) {
  
  # Build variable names for deaths and population
  death_var <- paste0('deaths_',this_age_group)
  pop_var <- paste0('pop_',this_age_group)
  
  # Estimate temperature response with contemporaneous (lag 0) basis terms and FE
  model <- feols(as.formula(paste0(death_var, '/', pop_var, ' ~ temp_d1_l0 + temp_d2_l0 + temp_d3_l0 + temp_d4_l0 + precip_d1_l0 + precip_d2_l0 | fips^year + fips^month')), data = month_panel, weights = month_panel[[pop_var]], cluster = ~ state)
  
  # Subset non-missing cells with positive population
  nonmissing <- month_panel %>% filter(!is.na(!!sym(death_var)), !is.na(!!sym(pop_var)), !!sym(pop_var) != 0)
  
  # Typical mortality rate (baseline scaling)
  typical_rate <- fmean(nonmissing[[death_var]]/nonmissing[[pop_var]], na.rm = T)
  
  # Pick out temperature coefficients and vcov
  which_temp <- grepl('temp', names(coef(model)))
  
  model_coef <- coef(model)[which_temp]
  model_vcov <- vcov(model)[which_temp,which_temp]
  
  # Predict on raw polynomial basis
  preds <-  pred_temps%*%model_coef
  
  # Find MMT within plausible range
  mmt_options <- which(pred_temps_vector < 30 & pred_temps_vector > 0)
  mmt <- pred_temps_vector[mmt_options][which.min(preds[mmt_options])]
  
  # Center temperature polynomial at MMT
  pred_temps_scaled <- scale(pred_temps, center = poly(mmt, degree = 4, raw = T), scale = F)
  
  # Classify temperatures as cold/hot relative to MMT (not used later but kept)
  sides <- c(rep('cold',(min(mmt_options)-1) + which.min(preds[mmt_options])),rep('hot',length(pred_temps_vector) - ((min(mmt_options)-1) + which.min(preds[mmt_options]))))
  
  # Recompute predictions after centering
  preds <- pred_temps_scaled%*%model_coef
  
  # Delta-method standard errors
  ses <- sqrt(pmax(0,rowSums(pred_temps_scaled*(pred_temps_scaled%*%model_vcov))))
  
  # Label for this age group
  age_group_label <- age_group_labels[match(this_age_group, age_groups)]
  
  # Store results; violence = 1 indicates all-cause
  pred_df[[length(pred_df)+1]] <- data.table(temp = pred_temps_vector, pred = c(preds), age_group = age_group_label, se = ses, typical_rate = typical_rate, mmt = mmt, violence = 1)
  
}

# Loop over age groups, excluding violence-related deaths
for (this_age_group in age_groups) {
  
  # For non-violence we use deaths_<age>_noviolence
  death_var <- paste0('deaths_',this_age_group,'_noviolence')
  pop_var <- paste0('pop_',this_age_group)
  
  # Estimate temperature response for non-violent causes
  model <- feols(as.formula(paste0(death_var, '/', pop_var, ' ~ temp_d1_l0 + temp_d2_l0 + temp_d3_l0 + temp_d4_l0 + precip_d1_l0 + precip_d2_l0 | fips^year + fips^month')), data = month_panel, weights = month_panel[[pop_var]], cluster = ~ state)
  
  # Subset non-missing cells with positive population
  nonmissing <- month_panel %>% filter(!is.na(!!sym(death_var)), !is.na(!!sym(pop_var)), !!sym(pop_var) != 0)
  
  # Typical mortality rate (non-violent)
  typical_rate <- fmean(nonmissing[[death_var]]/nonmissing[[pop_var]], na.rm = T)
  
  # Temperature coefficients and vcov
  which_temp <- grepl('temp', names(coef(model)))
  
  model_coef <- coef(model)[which_temp]
  model_vcov <- vcov(model)[which_temp,which_temp]
  
  # Predict on raw polynomial basis
  preds <-  pred_temps%*%model_coef
  
  # Find MMT
  mmt_options <- which(pred_temps_vector < 30 & pred_temps_vector > 0)
  mmt <- pred_temps_vector[mmt_options][which.min(preds[mmt_options])]
  
  # Center temperature polynomial at MMT
  pred_temps_scaled <- scale(pred_temps, center = poly(mmt, degree = 4, raw = T), scale = F)
  
  # Cold/hot classification (not used downstream)
  sides <- c(rep('cold',(min(mmt_options)-1) + which.min(preds[mmt_options])),rep('hot',length(pred_temps_vector) - ((min(mmt_options)-1) + which.min(preds[mmt_options]))))
  
  # Recompute predictions after centering
  preds <- pred_temps_scaled%*%model_coef
  
  # Standard errors via delta method
  ses <- sqrt(pmax(0,rowSums(pred_temps_scaled*(pred_temps_scaled%*%model_vcov))))
  
  # Age group label
  age_group_label <- age_group_labels[match(this_age_group, age_groups)]
  
  # Store non-violence predictions; violence = 0
  pred_df[[length(pred_df)+1]] <- data.table(temp = pred_temps_vector, pred = c(preds), age_group = age_group_label, se = ses, typical_rate = typical_rate, mmt = mmt, violence = 0)
  
}

# Combine all-age-group predictions and order age factor
pred_df <- rbindlist(pred_df) %>% mutate(age_group = factor(age_group, levels = age_group_labels))

# Convert from log-risk scale to percent change in daily mortality:
# multiply by 30.437 (approx. average days per month) and divide by typical rate
pred_df %<>% mutate(pred = 30.437*(pred/typical_rate), se = 30.437*(se/typical_rate))

# Save U.S. age-response predictions
write_parquet(pred_df, sink = 'output/usa/age_response.pq')

## mexico response by age ----------------------------------------------------------

# Read Mexico monthly deaths and population, define age groups and merge with monthly weather
mexico_mort <- read_parquet('data/mexico/monthly_death_pop_income.pq') %>% 
  mutate(adm1 = substr(adm2, 1, 3), adm2 = as.integer(adm2)) %>% 
  mutate(
    deaths_lt1 = alt1, pop_lt1 = alt1_pop,
    deaths_1_4 = alt5,  pop_1_4 = alt5_pop,
    deaths_5_44  = alt15 + alt25 + alt35 + alt45,
    pop_5_44     = alt15_pop + alt25_pop + alt35_pop + alt45_pop,
    deaths_45_54 = alt55, pop_45_54 = alt55_pop,
    deaths_55_64 = alt65, pop_55_64 = alt65_pop,
    deaths_65_74 = alt75, pop_65_74 = alt75_pop,
    deaths_75_84 = alt85, pop_75_84 = alt85_pop,
    deaths_gt85  = altinf, pop_gt85 = altinf_pop,
    adm1 = substr(adm2, 1, 3)
  ) %>% 
  select(adm2, adm1, year, month, deaths_lt1, pop_lt1, deaths_1_4, pop_1_4, deaths_5_44, pop_5_44, deaths_45_54, pop_45_54, deaths_55_64, pop_55_64, deaths_65_74, pop_65_74, deaths_75_84, pop_75_84, deaths_gt85, pop_gt85)
# {https://www.dropbox.com/scl/fi/kdmxy0kyblzdpjc8rudtz/monthly_death_pop_income.pq?rlkey=487n9d7gclc1hcv65ufwkfsbl&dl=0}

# Monthly Mexican weather with lag 0 metrics
monthly_weather <- read_parquet('data/mexico/met/monthly_weather.pq') %>% select(adm2, year, month, ends_with('_l0'))
# {https://www.dropbox.com/scl/fi/a25gcczxtx7jhu7vv1g0p/monthly_weather.pq?rlkey=5vpl3qdxu0jijzl3d6gni2dlc&dl=0}

# Combine monthly deaths/pop and weather
mexico_month_panel <- merge(mexico_mort, monthly_weather, by = c('adm2','year','month'))

# Mexican age groups and labels (broader than U.S. in working ages)
mex_age_groups <- c('lt1', '1_4', '5_44', '45_54', '55_64', '65_74', '75_84', 'gt85')
mex_age_group_labels <- c('Under 1', '1 to 4', '5 to 44', '45 to 54', '55 to 64', '65 to 74', '75 to 84', 'Above 84')

# Temperature prediction grid for Mexico age-response
pred_temps_vector <- -50:50
pred_temps <- poly(pred_temps_vector, degree = 4, raw = T)

# List to store Mexican age-response predictions
mex_pred_df <- list()

# Loop across Mexican age groups
for (this_age_group in mex_age_groups) {
  
  # Build deaths and population variable names
  death_var <- paste0('deaths_',this_age_group)
  pop_var <- paste0('pop_',this_age_group)
  
  # Temperature response model with lag 0 for Mexico, FE on adm2×year and adm2×month
  model <- feols(as.formula(paste0(death_var, '/', pop_var, ' ~ temp1_l0 + temp2_l0 + temp3_l0 + temp4_l0 + precip1_l0 + precip2_l0 | adm2^year + adm2^month')), data = mexico_month_panel, weights = mexico_month_panel[[pop_var]], cluster = ~ adm1)
  
  # Subset non-missing cells
  nonmissing <- mexico_month_panel %>% filter(!is.na(!!sym(death_var)), !is.na(!!sym(pop_var)), !!sym(pop_var) != 0)
  
  # Typical monthly mortality rate
  typical_rate <- fmean(nonmissing[[death_var]]/nonmissing[[pop_var]], na.rm = T)
  
  # Extract temperature coefficients and vcov
  which_temp <- grepl('temp', names(coef(model)))
  
  model_coef <- coef(model)[which_temp]
  model_vcov <- vcov(model)[which_temp,which_temp]
  
  # Predictions on raw temp polynomial
  preds <-  pred_temps%*%model_coef
  
  # Identify MMT
  mmt_options <- which(pred_temps_vector < 30 & pred_temps_vector > 0)
  mmt <- pred_temps_vector[mmt_options][which.min(preds[mmt_options])]
  
  # Center temp polynomial at MMT
  pred_temps_scaled <- scale(pred_temps, center = poly(mmt, degree = 4, raw = T), scale = F)
  
  # Cold/hot classification (not used later)
  sides <- c(rep('cold',(min(mmt_options)-1) + which.min(preds[mmt_options])),rep('hot',length(pred_temps_vector) - ((min(mmt_options)-1) + which.min(preds[mmt_options]))))
  
  # Recompute predictions after centering
  preds <- pred_temps_scaled%*%model_coef
  
  # Standard errors
  ses <- sqrt(pmax(0,rowSums(pred_temps_scaled*(pred_temps_scaled%*%model_vcov))))
  
  # Age group label
  age_group_label <- mex_age_group_labels[match(this_age_group, mex_age_groups)]
  
  # Store Mexican age-response; violence = 1 (all causes)
  mex_pred_df[[length(mex_pred_df)+1]] <- data.table(temp = pred_temps_vector, pred = c(preds), age_group = age_group_label, se = ses, typical_rate = typical_rate, mmt = mmt, violence = 1)
  
}

# Combine Mexican age-response predictions and order factor
mex_pred_df <- rbindlist(mex_pred_df) %>% mutate(age_group = factor(age_group, levels = mex_age_group_labels))

# Convert to percent change in daily mortality using 30.437 scaling
mex_pred_df %<>% mutate(pred = 30.437*(pred/typical_rate), se = 30.437*(se/typical_rate))

# Save Mexico age-response predictions
write_parquet(mex_pred_df, sink = 'output/mexico/age_response.pq')

## full figure ---------------------------------------------------------------------

# Read age- and lag-response prediction outputs for both countries
mex_age <- read_parquet('output/mexico/age_response.pq')
us_age <- read_parquet('output/usa/age_response.pq')
mex_lag <- read_parquet('output/mexico/lag_response.pq')
us_lag <- read_parquet('output/usa/lag_response.pq')

# Stack U.S. and Mexico age-response results, filter to violence == 1 and label type
age_plot_df <- rbind(mex_age %>% mutate(geo = 'Mexico'), us_age %>% mutate(geo = 'U.S.')) %>% 
  filter(violence == 1) %>% 
  mutate(type = 'Age') %>% 
  select(pred, age_group, temp, se, geo)

# Stack U.S. and Mexico lag-response results
lag_plot_df <- rbind(mex_lag %>% mutate(geo = 'Mexico'), us_lag %>% mutate(geo = 'U.S.')) %>% 
  mutate(type = 'Lag') %>% 
  select(pred, lag, temp, se, geo)

# Age-response panel (U.S. vs Mexico), restricted to key age groups and relevant temp ranges
age_plots <- ggplot(data = age_plot_df %>% filter((geo == 'U.S.' & temp >= 0 & temp <= 35) | (geo == 'Mexico' & temp >= 5 & temp <= 35),
                                                  age_group %notin% c('5 to 14', '15 to 24', '25 to 34', '35 to 44'))) +
  geom_hline(yintercept = 0, color = 'gray90') +
  geom_ribbon(aes(x = temp, y = pred, ymin = pred-1.96*se, ymax = pred+1.96*se, color = age_group, fill = age_group), alpha = 0.1, color = NA) +
  geom_line(aes(x = temp, y = pred, color = age_group)) +
  scale_color_manual(values = age_color_scale2, name = 'Age') + 
  scale_fill_manual(values = age_color_scale2, name = 'Age') + 
  theme_clean() + 
  ylab('Percent change in daily mortality') + 
  scale_y_continuous(labels = scales::percent_format()) + 
  xlab('Temperature (ºC)') +
  facet_wrap(facets = vars(geo), nrow = 1, scales = 'free') +
  theme(strip.background = element_blank(), strip.text = element_text(color = "black", face = "bold", size = 10))

# Lag-response panel (U.S. vs Mexico), selected lags by 7-day intervals and scaled per 100,000
lag_plots <- ggplot(data = lag_plot_df %>% filter((geo == 'U.S.' & temp >= 0 & temp <= 35) | (geo == 'Mexico' & temp >= 5 & temp <= 35), lag %in% seq(0,28,7)),
                    aes(x = temp, y = pred*1e5, ymin = pred*1e5 - 1.96*se*1e5, ymax = pred*1e5 + 1.96*se*1e5, color = as.factor(lag), fill = as.factor(lag))) +
  geom_ribbon(alpha = 0.1, color = NA) +
  geom_hline(yintercept = 0, color = 'gray90') +
  geom_line() +
  scale_fill_manual(values = lag_color_scale, name = 'Lag') + 
  scale_color_manual(values = lag_color_scale, name = 'Lag') + 
  theme_clean() +
  ylab('Deaths per 100,000') + 
  xlab('Temperature (ºC)') +
  facet_wrap(facets = vars(geo), nrow = 1, , scales = 'free') +
  theme(strip.background = element_blank(), strip.text = element_text(color = "black", face = "bold", size = 10))

# Combine age and lag panels into a single figure
full_plot <- cowplot::plot_grid(age_plots, lag_plots, nrow = 2, axis  = 'tlbr',align= 'hv')

# Save final figure as PDF
ggsave(full_plot, file = 'output/figures/lag_age_response.pdf', height = 8, width = 10, units = 'in')