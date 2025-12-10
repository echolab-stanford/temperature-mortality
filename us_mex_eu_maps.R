# maps of single-day effects in the U.S., Mexico, and EU
pacman::p_load(fastverse, tidyverse, arrow, pammtools, fixest, tictoc, sf, ggnewscale)

# Set working directory
setwd("~/path/to/temperature-mortality")

set.seed(42)
fixest::setFixest_nthreads(8)
`%notin%` <- Negate(`%in%`)
options(arrow.pull_as_vector = TRUE)
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

pred_temps_vector <- -50:50
pred_temps <- poly(pred_temps_vector, degree = 4, raw = T)

#### MEX
mex_df <- read_parquet('data/MEX/full_monthly_small.pq')
mex_df %<>% mutate(rate = deaths/pop) #%>% mutate(rate = ifelse(rate > fnth(rate, 0.99), fnth(rate, 0.99), rate)) #incomes are already in PPP- and inflation-adjusted international dollars

mex_reg <- feols(rate ~
                   temp1_l0 + temp2_l0 + temp3_l0 + temp4_l0 +
                   (temp1_l0 + temp2_l0 + temp3_l0 + temp4_l0):avg_temp +
                   (temp1_l0 + temp2_l0 + temp3_l0 + temp4_l0):log_pc_income +
                   (temp1_l0 + temp2_l0 + temp3_l0 + temp4_l0):log_bmr +
                   precip1_l0 + precip2_l0
                 | adm_id^year + adm_id^month, data = mex_df %>% filter(year < 2020), se = 'standard')

summary(mex_reg)
mex_reg <- coef(mex_reg)[grepl('temp',names(coef(mex_reg)))]

adm_list <- unique(mex_df$adm_id)
mex_adm_values <- list()
for (this_adm in adm_list) {
  
  cat(this_adm, '\n')
  
  for_2023 <- mex_df %>% filter(adm_id == this_adm, year == 2019, month == 1) %>% select(avg_temp, sd_temp, log_pc_income, log_bmr)
  for_2023$typical_rate <- mex_df %>% filter(adm_id == this_adm, year > 2009) %>% group_by(year) %>% summarize(rate_total = sum(rate)/365.25, .groups = 'keep') %>% ungroup() %>% summarize(rate_total = mean(rate_total), .groups = 'keep') %>% pull(rate_total)
  preds <-
    pred_temps%*%mex_reg[1:4] +
    for_2023$avg_temp*pred_temps%*%mex_reg[5:8] +
    for_2023$log_pc_income*pred_temps%*%mex_reg[9:12] +
    for_2023$log_bmr*pred_temps%*%mex_reg[13:16]
  mmt_options <- which(pred_temps_vector > for_2023$avg_temp - for_2023$sd_temp & pred_temps_vector < for_2023$avg_temp + for_2023$sd_temp & pred_temps_vector < 30 & pred_temps_vector > -10)
  mmt <- pred_temps_vector[mmt_options][which.min(preds[mmt_options])]
  sides <- c(rep('cold',(min(mmt_options)-1) + which.min(preds[mmt_options])),rep('hot',length(pred_temps_vector) - ((min(mmt_options)-1) + which.min(preds[mmt_options]))))
  preds <- preds - min(preds[mmt_options])
  pred_df <- data.table(temp = pred_temps_vector, pred = c(preds))
  
  mex_adm_values[[length(mex_adm_values)+1]] <- data.table(adm2 = this_adm,
                                                           temp = pred_temps_vector,
                                                           avg_temp = for_2023$avg_temp,
                                                           pc_income = for_2023$pc_income,
                                                           bmr = exp(for_2023$log_bmr),
                                                           percent_pred = c(preds)/for_2023$typical_rate,
                                                           mmt = mmt,
                                                           side = sides)
}
mex_adm_values <- rbindlist(mex_adm_values)

mex_weather <- open_dataset('data/MEX/mexico_temp.pq', format = 'parquet') %>% filter(year > 1979) %>% select(temp = order_1, year, order = poly_id) %>% collect()
mex_adm2_mapping <- open_dataset('data/MEX/adm2_to_order.csv', format = 'csv')
mex_weather %<>% merge(mex_adm2_mapping, by = 'order') %>% select(-order, temp, year, poly_id = adm2)

min_max_temps <- mex_weather %>% group_by(poly_id) %>% summarize(min_temp = min(temp, na.rm = T), max_temp = max(temp, na.rm = T)) %>% collect()

mex_adm_values %<>% merge(min_max_temps %>% mutate(adm2 = as.integer(poly_id)) %>% select(-poly_id), by = 'adm2') %>% mutate(exclude = ifelse(temp < min_temp - 1 | temp > max_temp + 1, 1, 0))

mex_polygons <- st_read('data/MEX/mexico.gpkg') %>% mutate(adm2 = as.integer(adm2))

mex_polygons_multiples <- rbind(mex_polygons %>% merge(mex_adm_values %>% filter(temp == 35), by = 'adm2', all.x = T),
                                mex_polygons %>% merge(mex_adm_values %>% filter(temp == 30), by = 'adm2', all.x = T),
                                mex_polygons %>% merge(mex_adm_values %>% filter(temp == 25), by = 'adm2', all.x = T),
                                mex_polygons %>% merge(mex_adm_values %>% filter(temp == 20), by = 'adm2', all.x = T),
                                mex_polygons %>% merge(mex_adm_values %>% filter(temp == 15), by = 'adm2', all.x = T),
                                mex_polygons %>% merge(mex_adm_values %>% filter(temp == 10), by = 'adm2', all.x = T),
                                mex_polygons %>% merge(mex_adm_values %>% filter(temp == 5), by = 'adm2', all.x = T))

mex_polygons_multiples %<>% mutate(temp = paste0(as.character(temp), 'ºC day')) %>% mutate(temp = factor(temp, levels = rev(unique(temp))))

temp_plot_mex_1030 <- ggplot() +
  geom_sf(data = mex_polygons_multiples %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(side == 'cold', temp %in% c('10ºC day','30ºC day')), aes(geometry = geom, fill = percent_pred), color = 'gray80', linewidth = 0.05) +
  scale_fill_gradientn(labels = scales::percent_format(), name = '% change, cold', limits = c(0,.5), oob = scales::squish, na.value = 'white', breaks = c(0, 0.25, 0.5), values = scales::rescale(seq(0, 0.5, length.out = 6)), colors = c('#eaeccc','#c2e4ef','#98cae1','#6ea6cd','#4a7bb7','#364b9a'), guide = guide_colorbar(order = 1)) +
  new_scale_fill() +
  geom_sf(data = mex_polygons_multiples %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(side == 'hot', temp %in% c('10ºC day','30ºC day')), aes(geometry = geom, fill = percent_pred), color = 'gray80', linewidth = 0.05) +
  scale_fill_gradientn(labels = scales::percent_format(), name = '% change, heat', limits = c(0,.3), oob = scales::squish, na.value = 'white', breaks = c(0, 0.15, 0.3), values = scales::rescale(seq(0, 0.3, length.out = 6)), colors = c('#eaeccc','#FEDA8b','#fdb366','#f67e4b','#dd3d2d','#a50026'), guide = guide_colorbar(order = 2)) +
  facet_wrap(facets = vars(as.factor(temp)), ncol = 2) + theme_clean() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line = element_blank()) + ggtitle('Change in daily mortality rate relative to temperature of minimum risk')
ggsave(temp_plot_mex_1030, file = 'fig/combined/county_maps_mex_1030.png', width = 10, height = 5, units = 'in')

#### EU
eu_panel <- read_parquet('data/EU/full_weekly.pq')
eu_panel %<>% mutate(rate = deaths/pop) %>% mutate(rate = ifelse(rate > fnth(rate, 0.99), fnth(rate, 0.99), rate))

eu_avg_temps <- read_parquet('data/EU/eu_avg_temps.pq')
eu_temp_limits <- read_parquet('data/EU/eu_temp_limits.pq')

eur_conversion_factor <- 1.11 # to convert to dollars (avg EUR to USD exchange rate in 2015, which is the index year)

eu_panel$log_pc_income <- log(eu_panel$pc_income*eur_conversion_factor)
eu_panel$income <- exp(eu_panel$log_pc_income)
eu_panel$bmr <- exp(eu_panel$log_bmr)

eu_reg <- feols(rate ~
                  temp_d1_l0 + temp_d2_l0 + temp_d3_l0 + temp_d4_l0 +
                  (temp_d1_l0 + temp_d2_l0 + temp_d3_l0 + temp_d4_l0):avg_temp +
                  (temp_d1_l0 + temp_d2_l0 + temp_d3_l0 + temp_d4_l0):log_pc_income +
                  (temp_d1_l0 + temp_d2_l0 + temp_d3_l0 + temp_d4_l0):log_bmr +
                  temp_d1_l1 + temp_d2_l1 + temp_d3_l1 + temp_d4_l1 +
                  (temp_d1_l1 + temp_d2_l1 + temp_d3_l1 + temp_d4_l1):avg_temp +
                  (temp_d1_l1 + temp_d2_l1 + temp_d3_l1 + temp_d4_l1):log_pc_income +
                  (temp_d1_l1 + temp_d2_l1 + temp_d3_l1 + temp_d4_l1):log_bmr +
                  temp_d1_l2 + temp_d2_l2 + temp_d3_l2 + temp_d4_l2 +
                  (temp_d1_l2 + temp_d2_l2 + temp_d3_l2 + temp_d4_l2):avg_temp +
                  (temp_d1_l2 + temp_d2_l2 + temp_d3_l2 + temp_d4_l2):log_pc_income +
                  (temp_d1_l2 + temp_d2_l2 + temp_d3_l2 + temp_d4_l2):log_bmr +
                  temp_d1_l3 + temp_d2_l3 + temp_d3_l3 + temp_d4_l3 +
                  (temp_d1_l3 + temp_d2_l3 + temp_d3_l3 + temp_d4_l3):avg_temp +
                  (temp_d1_l3 + temp_d2_l3 + temp_d3_l3 + temp_d4_l3):log_pc_income +
                  (temp_d1_l3 + temp_d2_l3 + temp_d3_l3 + temp_d4_l3):log_bmr +
                  precip_d1_l0 + precip_d2_l0 +
                  precip_d1_l1 + precip_d2_l1 +
                  precip_d1_l2 + precip_d2_l2 +
                  precip_d1_l3 + precip_d2_l3
                | adm_id^year + adm_id^week, data = eu_panel, se = 'standard')
summary(eu_reg)
eu_reg <- coef(eu_reg)[grepl('temp',names(coef(eu_reg)))]

nuts_list <- unique(eu_panel %>% filter(!is.na(rate), !is.na(pc_income)) %>% pull(adm_id))
eu_nuts_values <- list()
for (this_nuts in nuts_list) {
  
  cat(this_nuts, '\n')
  
  nuts_temp_sd <- eu_avg_temps %>% filter(adm_id == this_nuts, year == 2023) %>% pull(sd_temp)
  nuts_temp_avg <- eu_avg_temps %>% filter(adm_id == this_nuts, year == 2023) %>% pull(avg_temp)
  
  recent_year <- eu_panel %>% filter(adm_id == this_nuts, year > 2012) %>% select(income, bmr) %>% summarize(log_pc_income = log(mean(income, na.rm = T)), log_bmr = log(mean(bmr, na.rm = T)))
  recent_year$typical_rate <- eu_panel %>% filter(adm_id == this_nuts, year > 2012) %>% group_by(year) %>% summarize(rate_total = sum(rate)/365.25, .groups = 'keep') %>% ungroup() %>% summarize(rate_total = mean(rate_total), .groups = 'keep') %>% pull(rate_total)
  recent_year$avg_temp <- nuts_temp_avg
  recent_year$sd_temp <- nuts_temp_sd
  preds <-
    pred_temps%*%eu_reg[1:4] +
    pred_temps%*%eu_reg[5:8] +
    pred_temps%*%eu_reg[9:12] +
    pred_temps%*%eu_reg[13:16] +
    recent_year$avg_temp*pred_temps%*%eu_reg[17:20] +
    recent_year$log_pc_income*pred_temps%*%eu_reg[21:24] +
    recent_year$log_bmr*pred_temps%*%eu_reg[25:28] +
    recent_year$avg_temp*pred_temps%*%eu_reg[29:32] +
    recent_year$log_pc_income*pred_temps%*%eu_reg[33:36] +
    recent_year$log_bmr*pred_temps%*%eu_reg[37:40] +
    recent_year$avg_temp*pred_temps%*%eu_reg[41:44] +
    recent_year$log_pc_income*pred_temps%*%eu_reg[45:48] +
    recent_year$log_bmr*pred_temps%*%eu_reg[49:52] +
    recent_year$avg_temp*pred_temps%*%eu_reg[53:56] +
    recent_year$log_pc_income*pred_temps%*%eu_reg[57:60] +
    recent_year$log_bmr*pred_temps%*%eu_reg[61:64]
  mmt_options <- which(pred_temps_vector > nuts_temp_avg - nuts_temp_sd & pred_temps_vector < nuts_temp_avg + nuts_temp_sd & pred_temps_vector < 30 & pred_temps_vector > -10)
  mmt <- pred_temps_vector[mmt_options][which.min(preds[mmt_options])]
  sides <- c(rep('cold',(min(mmt_options)-1) + max(round(length(preds[mmt_options])/2),which.min(preds[mmt_options]), na.rm = T)),rep('hot',length(pred_temps_vector) - ((min(mmt_options)-1) + max(round(length(preds[mmt_options])/2), which.min(preds[mmt_options]), na.rm = T))))
  preds <- preds - min(preds[mmt_options])
  pred_df <- data.table(temp = pred_temps_vector, pred = c(preds))
  
  eu_nuts_values[[length(eu_nuts_values)+1]] <- data.table(nuts = this_nuts,
                                                           temp = pred_temps_vector,
                                                           avg_temp = recent_year$avg_temp,
                                                           pc_income = recent_year$pc_income,
                                                           bmr = recent_year$bmr,
                                                           percent_pred = c(preds)/recent_year$typical_rate,
                                                           mmt = mmt,
                                                           side = sides)
}
eu_nuts_values <- rbindlist(eu_nuts_values)

eu_nuts_values %<>% merge(eu_temp_limits, by.x = 'nuts', by.y = 'adm_id', all.x = T) %>% mutate(exclude = ifelse(temp < min_temp - 1 | temp > max_temp + 1, 1, 0))

eu_polygons <- st_read('data/EU/eu_nuts3_2024.gpkg') %>% mutate(nuts = NUTS_ID)
eu_polygons_multiples <- rbind(eu_polygons %>% merge(eu_nuts_values %>% filter(temp == 35), by = 'nuts', all.x = T),
                               eu_polygons %>% merge(eu_nuts_values %>% filter(temp == 30), by = 'nuts', all.x = T),
                               eu_polygons %>% merge(eu_nuts_values %>% filter(temp == 25), by = 'nuts', all.x = T),
                               eu_polygons %>% merge(eu_nuts_values %>% filter(temp == 20), by = 'nuts', all.x = T),
                               eu_polygons %>% merge(eu_nuts_values %>% filter(temp == 15), by = 'nuts', all.x = T),
                               eu_polygons %>% merge(eu_nuts_values %>% filter(temp == 10), by = 'nuts', all.x = T),
                               eu_polygons %>% merge(eu_nuts_values %>% filter(temp == 5), by = 'nuts', all.x = T),
                               eu_polygons %>% merge(eu_nuts_values %>% filter(temp == 0), by = 'nuts', all.x = T),
                               eu_polygons %>% merge(eu_nuts_values %>% filter(temp == -10), by = 'nuts', all.x = T))
eu_polygons_multiples %<>% filter(substr(nuts, 1, 2) %notin% c('TR'))

eu_polygons_multiples %<>% mutate(temp = paste0(as.character(temp), 'ºC day')) %>% mutate(temp = factor(temp, levels = rev(unique(temp))))

world <- st_read('data/WORLD/geoBoundariesCGAZ_ADM0.gpkg')

temp_plot_eu_030 <- ggplot() +
  geom_sf(data = world, aes(geometry = geom), fill = 'gray95', color = NA) +
  geom_sf(data = eu_polygons_multiples %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(side == 'cold', temp %in% c('0ºC day','25ºC day')), aes(geometry = geom, fill = percent_pred), color = 'gray85', linewidth = 0.05) +
  scale_fill_gradientn(labels = scales::percent_format(), name = '% change, cold', limits = c(0,.6), oob = scales::squish, na.value = 'white', breaks = c(0, 0.3, 0.6), values = scales::rescale(seq(0, 0.6, length.out = 6)), colors = c('#eaeccc','#c2e4ef','#98cae1','#6ea6cd','#4a7bb7','#364b9a'), guide = guide_colorbar(order = 1)) +
  new_scale_fill() +
  geom_sf(data = eu_polygons_multiples %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(side == 'hot', temp %in% c('0ºC day','25ºC day')), aes(geometry = geom, fill = percent_pred), color = 'gray85', linewidth = 0.05) +
  scale_fill_gradientn(labels = scales::percent_format(), name = '% change, heat', limits = c(0,0.2), oob = scales::squish, na.value = 'white', breaks = c(0, 0.1, 0.2), values = scales::rescale(seq(0, 0.2, length.out = 6)), colors = c('#eaeccc','#FEDA8b','#fdb366','#f67e4b','#dd3d2d','#a50026'), guide = guide_colorbar(order = 2)) +
  facet_wrap(facets = vars(as.factor(temp)), ncol = 2) + theme_clean() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line = element_blank()) + ggtitle('Change in daily mortality rate relative to temperature of minimum risk') + coord_sf(xlim = c(-14,31), ylim = c(33,75)) + theme(aspect.ratio = 1)
ggsave(temp_plot_eu_030, file = 'fig/combined/county_maps_eu_030.png', width = 10, height = 5, units = 'in')

#### U.S. --- mortality data for the U.S. is not publicly available, so it cannot be provided for this replication
month_panel <- read_parquet('---') %>% select(fips:precip_d2_l0)

county_polygons <- st_read('data/US/CONUS_consistent_counties.gpkg') %>% rename(fips = GEOID) %>% mutate(fips = as.integer(fips))

baseline_mort_rates <- list()
for (this_year in 1968:2023) {
  cat('now working on', this_year, '\n')
  baseline_mort_rates[[length(baseline_mort_rates)+1]] <- month_panel %>% filter(year < this_year + 4, year > this_year - 4) %>% group_by(fips) %>% summarize(bmr_adj = mean(age_adj_rate, na.rm = T)) %>% as.data.table() %>% mutate(year = this_year)
}
baseline_mort_rates <- rbindlist(baseline_mort_rates)

month_panel %<>% merge(baseline_mort_rates, by = c('year','fips'))

month_panel %<>% mutate(mort_rate = deaths_total/pop_total, homicide_rate = homicide/pop_total, log_pc_income = log(pc_income), log_bmr = log(bmr), log_bmr_adj = log(bmr_adj))

# n_bins <- 5
reg_df <- month_panel %>%
  mutate(mort_rate = ifelse(mort_rate > fnth(mort_rate, 0.99), fnth(mort_rate, 0.99), mort_rate)) %>%
  mutate(age_adj_rate_w = ifelse(age_adj_rate > fnth(age_adj_rate, 0.99), fnth(age_adj_rate, 0.99), age_adj_rate))

reg_df %<>% group_by(year) %>% mutate(weight = pop_total/sum(pop_total)) %>% as.data.table()

int_reg <- feols(age_adj_rate ~
                   temp_d1_l0 + temp_d2_l0 + temp_d3_l0 + temp_d4_l0 +
                   (temp_d1_l0 + temp_d2_l0 + temp_d3_l0 + temp_d4_l0):avg_temp +
                   (temp_d1_l0 + temp_d2_l0 + temp_d3_l0 + temp_d4_l0):log_pc_income +
                   (temp_d1_l0 + temp_d2_l0 + temp_d3_l0 + temp_d4_l0):bmr_adj +
                   precip_d1_l0 + precip_d2_l0
                 | fips^year + fips^month, data = reg_df, se = 'standard')
summary(int_reg)
int_reg <- coef(int_reg)[grepl('temp_',names(coef(int_reg)))]

pred_temps_vector <- seq(-50,50)
pred_temps <- poly(pred_temps_vector, degree = 4, raw = T)
fips_list <- unique(month_panel$fips)
county_values <- list()
for (this_fips in fips_list) {
  
  cat(this_fips, '\n')
  
  for_2023 <- reg_df %>% filter(fips == this_fips, year == 2023, month == 1) %>% select(pc_income, avg_temp, sd_temp, bmr, avg_age, log_pc_income, log_bmr, bmr_adj)
  for_2023$typical_rate <- reg_df %>% filter(fips == this_fips, year > 2013) %>% group_by(year) %>% summarize(rate_total = sum(age_adj_rate)/365.25, .groups = 'keep') %>% ungroup() %>% summarize(rate_total = mean(rate_total), .groups = 'keep') %>% pull(rate_total)
  preds <-
    pred_temps%*%int_reg[1:4] +
    for_2023$avg_temp*pred_temps%*%int_reg[5:8] +
    for_2023$log_pc_income*pred_temps%*%int_reg[9:12] +
    for_2023$bmr_adj*pred_temps%*%int_reg[13:16]
  mmt_options <- which(pred_temps_vector > for_2023$avg_temp - for_2023$sd_temp & pred_temps_vector < for_2023$avg_temp + for_2023$sd_temp & pred_temps_vector < 30 & pred_temps_vector > -10)
  mmt <- pred_temps_vector[mmt_options][which.min(preds[mmt_options])]
  # mmt <- ifelse(pred_temps_vector[mmt_options][which.min(preds[mmt_options])] > pred_temps_vector[mmt_options[1]] & pred_temps_vector[mmt_options][which.min(preds[mmt_options])] < pred_temps_vector[mmt_options[length(mmt_options)]],  pred_temps_vector[mmt_options][which.min(preds[mmt_options])], pred_temps_vector[mmt_options[1]])
  sides <- c(rep('cold',(min(mmt_options)-1) + which.min(preds[mmt_options])),rep('hot',length(pred_temps_vector) - ((min(mmt_options)-1) + which.min(preds[mmt_options]))))
  preds <- preds - preds[pred_temps_vector == mmt]
  pred_df <- data.table(temp = pred_temps_vector, pred = c(preds))
  
  county_values[[length(county_values)+1]] <- data.table(fips = this_fips,
                                                         temp = pred_temps_vector,
                                                         pred = c(preds),
                                                         avg_temp = for_2023$avg_temp,
                                                         # sd_temp = for_2023$sd_temp,
                                                         pc_income = exp(for_2023$log_pc_income),
                                                         bmr = for_2023$bmr_adj,
                                                         # avg_age = for_2023$avg_age,
                                                         percent_pred = c(preds)/for_2023$typical_rate,
                                                         side = sides
  )
}
county_values <- rbindlist(county_values)

min_max_temps <- open_dataset('data/US/us_temp.pq', format = 'parquet') %>% group_by(poly_id) %>% summarize(min_temp = min(order_1, na.rm = T), max_temp = max(order_1, na.rm = T)) %>% collect()

county_values %<>% merge(min_max_temps %>% mutate(fips = as.integer(poly_id)) %>% select(-poly_id), by = 'fips') %>% mutate(exclude = ifelse(temp < min_temp - 1 | temp > max_temp + 1, 1, 0))

county_polygons_multiples <- rbind(county_polygons %>% merge(county_values %>% filter(temp == 35), by = 'fips', all.x = T),
                                   county_polygons %>% merge(county_values %>% filter(temp == 30), by = 'fips', all.x = T),
                                   county_polygons %>% merge(county_values %>% filter(temp == 25), by = 'fips', all.x = T),
                                   county_polygons %>% merge(county_values %>% filter(temp == 20), by = 'fips', all.x = T),
                                   county_polygons %>% merge(county_values %>% filter(temp == 15), by = 'fips', all.x = T),
                                   county_polygons %>% merge(county_values %>% filter(temp == 10), by = 'fips', all.x = T),
                                   county_polygons %>% merge(county_values %>% filter(temp == 5), by = 'fips', all.x = T),
                                   county_polygons %>% merge(county_values %>% filter(temp == 0), by = 'fips', all.x = T),
                                   county_polygons %>% merge(county_values %>% filter(temp == -10), by = 'fips', all.x = T))

county_polygons_multiples %<>% mutate(temp = paste0(as.character(temp), 'ºC day')) %>% mutate(temp = factor(temp, levels = rev(unique(temp))))

temp_plot_bottomcode_030 <- ggplot() +
  geom_sf(data = county_polygons_multiples %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(side == 'cold', temp %in% c('0ºC day','30ºC day')), aes(geometry = geom, fill = percent_pred), color = 'gray80', linewidth = 0.1) +
  scale_fill_gradientn(labels = scales::percent_format(), name = '% change, cold', limits = c(0,.3), oob = scales::squish, na.value = 'white', breaks = c(0, 0.15, 0.3), values = scales::rescale(seq(0, 0.3, length.out = 6)), colors = c('#eaeccc','#c2e4ef','#98cae1','#6ea6cd','#4a7bb7','#364b9a'), guide = guide_colorbar(order = 1)) +
  new_scale_fill() +
  geom_sf(data = county_polygons_multiples %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(side == 'hot', temp %in% c('0ºC day','30ºC day')), aes(geometry = geom, fill = percent_pred), color = 'gray80', linewidth = 0.1) +
  scale_fill_gradientn(labels = scales::percent_format(), name = '% change, heat', limits = c(0,.1), oob = scales::squish, na.value = 'white', breaks = c(0, 0.05, 0.1), values = scales::rescale(seq(0, 0.1, length.out = 6)), colors = c('#eaeccc','#FEDA8b','#fdb366','#f67e4b','#dd3d2d','#a50026'), guide = guide_colorbar(order = 2)) +
  facet_wrap(facets = vars(as.factor(temp)), ncol = 2) + theme_clean() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line = element_blank()) + ggtitle('Change in daily mortality rate relative to temperature of minimum risk')
ggsave(temp_plot_bottomcode_030, file = 'fig/combined/county_maps4_030.png', width = 10, height = 5, units = 'in')