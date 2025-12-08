# maps of single-day effects in the U.S., Mexico, and EU
pacman::p_load(fastverse, tidyverse, arrow, pammtools, fixest, tictoc, sf, ggnewscale)
# Set working directory
# setwd("~/path/to/temperature-mortality")
# setwd('/Users/aw/Dropbox/Research/temp_mortality')
setwd("~/BurkeLab Dropbox/projects/temperature-mortality")

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
okabe_colors <- c('#D55E00','#0072B2','#009E73','#CC79A7','#000000','#E69F00','#F0E442','#4DA3E5','#555555','#D9D9D9')
map_to_nuts2024 <- function(x) {case_when(x %in% c('AT224') ~ 'AT225',
                                          x %in% c('BE221','BE222','BE224') ~ 'BE225',
                                          x %in% c('BE324','BE327') ~ 'BE328',
                                          x %in% c('BE321','BE322','BE325','BE326') ~ 'BE329',
                                          x %in% c('DEG04','DEG0F','DEG0S') ~ 'DEG0T',
                                          x %in% c('DEG0B','DEG0N','DEG0P','DEG0Q') ~ 'DEG0R',
                                          x %in% c('DEG0H','DEG0I','DEG0V') ~ 'DEG0U',
                                          x %in% c('EE006') ~ 'EE009',
                                          x %in% c('EE007') ~ 'EE00A',
                                          x %in% c('FI193','FI197','FI198') ~ 'FI19B',
                                          x %in% c('FI1C3','FI1C4','FI1C6') ~ 'FI1C7',
                                          x %in% c('FI194','FI195','FI199') ~ 'FI19A',
                                          x %in% c('FI19D1','FI1D2','FI1DA','FI1DB','FI1D3') ~ 'FI1DC',
                                          x %in% c('HR041') ~ 'HR050',
                                          x %in% c('HR046') ~ 'HR061',
                                          x %in% c('HR044') ~ 'HR062',
                                          x %in% c('HR045') ~ 'HR063',
                                          x %in% c('HR043') ~ 'HR064',
                                          x %in% c('HR042') ~ 'HR065',
                                          x %in% c('HR047') ~ 'HR021',
                                          x %in% c('HR048') ~ 'HR022',
                                          x %in% c('HR049') ~ 'HR023',
                                          x %in% c('HR04A') ~ 'HR024',
                                          x %in% c('HR04B') ~ 'HR025',
                                          x %in% c('HR04C') ~ 'HR026',
                                          x %in% c('HR04D') ~ 'HR027',
                                          x %in% c('HR04E') ~ 'HR028',
                                          x %in% c('ITG25','ITG29') ~ 'ITG2D',
                                          x %in% c('ITG2A','ITG2B','ITG2C','ITG2E','ITG2F','ITG2G','ITG26','ITG27','ITG28') ~ 'ITG2H',
                                          x %in% c('ITH59') ~ 'ITI31',
                                          x %in% c('LV003','LV006','LV007','LV008','LV00B','LV00A') ~ 'LV00C',
                                          x %in% c('NL111','NL113','NL114') ~ 'NL115',
                                          x %in% c('NL124','NL125','NL127') ~ 'NL128',
                                          x %in% c('NL310','NL33A','NL350') ~ 'NL364',
                                          x %in% c('NL332') ~ 'NL361',
                                          x %in% c('NL333') ~ 'NL362',
                                          x %in% c('NL337') ~ 'NL363',
                                          x %in% c('NL33B') ~ 'NL365',
                                          x %in% c('NL33C') ~ 'NL366',
                                          x %in% c('NL324','NL329','NL32A') ~ 'NL32B',
                                          x %in% c('NL412','NL413','NL415') ~ 'NL416',
                                          x %in% c('NO011') ~ 'NO081',
                                          x %in% c('NO012','NO031','NO032','NO083','NO084','NO085') ~ 'NO082',
                                          x %in% c('NO021','NO022') ~ 'NO020',
                                          x %in% c('NO033','NO034','NO093','NO094') ~ 'NO091',
                                          x %in% c('NO041','NO042') ~ 'NO092',
                                          x %in% c('NO043') ~ 'NO0A1',
                                          x %in% c('NO051','NO052') ~ 'NO0A2',
                                          x %in% c('NO053') ~ 'NO0A3',
                                          x %in% c('NO061','NO062') ~ 'NO060',
                                          x %in% c('NO072','NO073') ~ 'NO074',
                                          x %in% c('PT16H','PT16I','PT195') ~ 'PT1D2',
                                          x %in% c('PT170','PT1A0') ~ 'PT1B0',
                                          x %in% c('PT16B') ~ 'PT1D1',
                                          x %in% c('PT16D') ~ 'PT191',
                                          x %in% c('PT16E') ~ 'PT192',
                                          x %in% c('PT16F') ~ 'PT193',
                                          x %in% c('PT16G') ~ 'PT194',
                                          x %in% c('PT181') ~ 'PT1C1',
                                          x %in% c('PT184') ~ 'PT1C2',
                                          x %in% c('PT185') ~ 'PT1D3',
                                          x %in% c('PT186') ~ 'PT1C3',
                                          x %in% c('PT187') ~ 'PT1C4',
                                          x %in% c('PT16J') ~ 'PT196',
                                          x %in% c('UKK21','UKK22','UKK24') ~ 'UKK25',
                                          x %in% c('UKN10') ~ 'UKN0A',
                                          x %in% c('UKN11') ~ 'UKN0B',
                                          x %in% c('UKN12') ~ 'UKN0C',
                                          x %in% c('UKN13') ~ 'UKN0D',
                                          x %in% c('UKN14') ~ 'UKN0E',
                                          x %in% c('UKN15') ~ 'UKN0F',
                                          x %in% c('UKN16') ~ 'UKN0G',
                                          TRUE ~ x)}

weights <- data.table(age_group = c('lt1', '1_4', '5_14', '15_24', '25_34', '35_44', '45_54', '55_64', '65_74', '75_84', 'gt85'),
                      group_size = c(3794901, 15191619, 39976619, 38076743, 37233437, 44659185, 37030152, 23961506, 18135514, 12314793, 4259173))
weights %<>% mutate(weight = group_size/sum(group_size))
"%c%" <- function(a,b) {
  mgcv::tensor.prod.model.matrix(list(as.matrix(a),as.matrix(b)))
}
pdens <- function(x) {plot(density(x, na.rm = T))}

pred_temps_vector <- -50:50
pred_temps <- poly(pred_temps_vector, degree = 4, raw = T)

#### MEX
mex_df <- read_parquet('data/MEX/full_monthly.pq')
mex_df %<>% mutate(rate = deaths/pop) #%>% mutate(rate = ifelse(rate > fnth(rate, 0.99), fnth(rate, 0.99), rate)) #incomes are already in PPP- and inflation-adjusted international dollars

mex_reg <- feols(rate ~
                   temp1_l0 + temp2_l0 + temp3_l0 + temp4_l0 +
                   (temp1_l0 + temp2_l0 + temp3_l0 + temp4_l0):avg_temp +
                   (temp1_l0 + temp2_l0 + temp3_l0 + temp4_l0):log_pc_income +
                   (temp1_l0 + temp2_l0 + temp3_l0 + temp4_l0):log_bmr +
                   precip1_l0 + precip2_l0
                 | adm_id^year + adm_id^month, data = mex_df, se = 'standard')

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
                                                           bmr = for_2023$bmr,
                                                           percent_pred = c(preds),
                                                           mmt = mmt,
                                                           side = sides)
}
mex_adm_values <- rbindlist(mex_adm_values)

mex_weather <- open_dataset('data/MEX/mexico_temp.pq', format = 'parquet') %>% filter(year > 1979) %>% select(temp = order_1, year, order = poly_id) %>% collect()
mex_adm2_mapping <- open_dataset('data/MEX/adm2_to_order.csv', format = 'csv')
mex_weather %<>% merge(mex_adm2_mapping, by = 'order') %>% select(-order, temp, year, poly_id = adm2)

min_max_temps <- mex_weather %>% group_by(poly_id) %>% summarize(min_temp = min(temp, na.rm = T), max_temp = max(temp, na.rm = T)) %>% collect()

mex_adm_values %<>% merge(min_max_temps %>% mutate(adm2 = as.integer(poly_id)) %>% select(-poly_id), by = 'adm2') %>% mutate(exclude = ifelse(temp < min_temp - 1 | temp > max_temp + 1, 1, 0))

mex_polygons <- st_read('/Users/aw/Dropbox/Research/wbt_mortality/data/shapefiles/harmonized_1990.shp') %>% mutate(adm2 = as.integer(adm2))
# {https://www.dropbox.com/scl/fi/ksbdkrjd76ayhqw3o64ov/mexico.gpkg?rlkey=3xq7rnbl6c1g5vfa78zwgvw2i&dl=0}

mex_polygons_multiples <- rbind(mex_polygons %>% merge(mex_adm_values %>% filter(temp == 35), by = 'adm2', all.x = T),
                                mex_polygons %>% merge(mex_adm_values %>% filter(temp == 30), by = 'adm2', all.x = T),
                                mex_polygons %>% merge(mex_adm_values %>% filter(temp == 25), by = 'adm2', all.x = T),
                                mex_polygons %>% merge(mex_adm_values %>% filter(temp == 20), by = 'adm2', all.x = T),
                                mex_polygons %>% merge(mex_adm_values %>% filter(temp == 15), by = 'adm2', all.x = T),
                                mex_polygons %>% merge(mex_adm_values %>% filter(temp == 10), by = 'adm2', all.x = T),
                                mex_polygons %>% merge(mex_adm_values %>% filter(temp == 5), by = 'adm2', all.x = T))

mex_polygons_multiples %<>% mutate(temp = paste0(as.character(temp), 'ºC day')) %>% mutate(temp = factor(temp, levels = rev(unique(temp))))

temp_plot_mex <- ggplot() +
  geom_sf(data = mex_polygons_multiples %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(side == 'hot'), aes(geometry = geometry, fill = percent_pred*30.437), color = 'gray80', linewidth = 0.1) +
  scale_fill_gradientn(labels = scales::percent_format(), name = '% change, heat', limits = c(0,.5), oob = scales::squish, na.value = 'white', breaks = c(0, 0.25, 0.5), values = scales::rescale(seq(0, 0.5, length.out = 6)), colors = c('#eaeccc','#FEDA8b','#fdb366','#f67e4b','#dd3d2d','#a50026')) +
  new_scale_fill() +
  geom_sf(data = mex_polygons_multiples %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(side == 'cold'), aes(geometry = geometry, fill = percent_pred*30.437), color = 'gray80', linewidth = 0.1) +
  scale_fill_gradientn(labels = scales::percent_format(), name = '% change, cold', limits = c(0,.4), oob = scales::squish, na.value = 'white', breaks = c(0, 0.2, 0.4), values = scales::rescale(seq(0, 0.4, length.out = 6)), colors = c('#eaeccc','#c2e4ef','#98cae1','#6ea6cd','#4a7bb7','#364b9a')) +
  facet_wrap(facets = vars(as.factor(temp)), nrow = 3) + theme_clean() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line = element_blank()) + ggtitle('Change in daily mortality rate relative to temperature of minimum risk')
ggsave(temp_plot_mex, file = 'fig/combined/county_maps_mex.png', width = 10, height = 8, units = 'in')

temp_plot_mex_1025 <- ggplot() +
  geom_sf(data = mex_polygons_multiples %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(side == 'cold', temp %in% c('10ºC day','25ºC day')), aes(geometry = geometry, fill = percent_pred*30.437), color = 'gray80', linewidth = 0.05) +
  scale_fill_gradientn(labels = scales::percent_format(), name = '% change, cold', limits = c(0,.4), oob = scales::squish, na.value = 'white', breaks = c(0, 0.2, 0.4), values = scales::rescale(seq(0, 0.4, length.out = 6)), colors = c('#eaeccc','#c2e4ef','#98cae1','#6ea6cd','#4a7bb7','#364b9a'), guide = guide_colorbar(order = 1)) +
  new_scale_fill() +
  geom_sf(data = mex_polygons_multiples %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(side == 'hot', temp %in% c('10ºC day','25ºC day')), aes(geometry = geometry, fill = percent_pred*30.437), color = 'gray80', linewidth = 0.05) +
  scale_fill_gradientn(labels = scales::percent_format(), name = '% change, heat', limits = c(0,.3), oob = scales::squish, na.value = 'white', breaks = c(0, 0.15, 0.3), values = scales::rescale(seq(0, 0.3, length.out = 6)), colors = c('#eaeccc','#FEDA8b','#fdb366','#f67e4b','#dd3d2d','#a50026'), guide = guide_colorbar(order = 2)) +
  facet_wrap(facets = vars(as.factor(temp)), ncol = 2) + theme_clean() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line = element_blank()) + ggtitle('Change in daily mortality rate relative to temperature of minimum risk')
ggsave(temp_plot_mex_1025, file = 'fig/combined/county_maps_mex_1025.png', width = 10, height = 5, units = 'in')

temp_plot_mex_3025 <- ggplot() +
  geom_sf(data = mex_polygons_multiples %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(side == 'cold', temp %in% c('10ºC day','25ºC day')), aes(geometry = geometry, fill = percent_pred*30.437), color = 'gray80', linewidth = 0.05) +
  scale_fill_gradientn(labels = scales::percent_format(), name = '% change, cold', limits = c(0,.4), oob = scales::squish, na.value = 'white', breaks = c(0, 0.2, 0.4), values = scales::rescale(seq(0, 0.4, length.out = 6)), colors = c('#eaeccc','#c2e4ef','#98cae1','#6ea6cd','#4a7bb7','#364b9a'), guide = guide_colorbar(order = 1)) +
  new_scale_fill() +
  geom_sf(data = mex_polygons_multiples %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(side == 'hot', temp %in% c('10ºC day','25ºC day')), aes(geometry = geometry, fill = percent_pred*30.437), color = 'gray80', linewidth = 0.05) +
  scale_fill_gradientn(labels = scales::percent_format(), name = '% change, heat', limits = c(0,.3), oob = scales::squish, na.value = 'white', breaks = c(0, 0.15, 0.3), values = scales::rescale(seq(0, 0.3, length.out = 6)), colors = c('#eaeccc','#FEDA8b','#fdb366','#f67e4b','#dd3d2d','#a50026'), guide = guide_colorbar(order = 2)) +
  facet_wrap(facets = vars(as.factor(temp)), ncol = 2) + theme_clean() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line = element_blank()) + ggtitle('Change in daily mortality rate relative to temperature of minimum risk')
ggsave(temp_plot_mex_1025, file = 'fig/combined/county_maps_mex_1025.png', width = 10, height = 5, units = 'in')

#### EU
eu_panel <- read_parquet('data/EU/full_weekly.pq')
# {https://www.dropbox.com/scl/fi/txd2ajouz26hjo4u2x7yt/full_weekly.pq?rlkey=kre5ziioq7c9g5uw38yo2kubm&dl=0}
eu_panel %<>% mutate(rate = deaths/pop) %>% mutate(rate = ifelse(rate > fnth(rate, 0.99), fnth(rate, 0.99), rate))

eur_conversion_factor <- 1.11 # to convert to dollars (avg EUR to USD exchange rate in 2015, which is the index year)

eu_panel$log_pc_income <- log(eu_panel$pc_income*eur_conversion_factor)

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
  
  recent_year <- eu_panel %>% filter(adm_id == this_nuts, year > 2012) %>% select(log_pc_income, log_bmr) %>% summarize(log_pc_income = mean(log_pc_income, na.rm = T), log_bmr = mean(log_bmr, na.rm = T))
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
                                                           percent_pred = c(preds),
                                                           mmt = mmt,
                                                           side = sides)
}
eu_nuts_values <- rbindlist(eu_nuts_values)

eu_nuts_values %<>% merge(eu_temp_limits, by.x = 'nuts', by.y = 'adm_id', all.x = T) %>% mutate(exclude = ifelse(temp < min_temp - 1 | temp > max_temp + 1, 1, 0))

eu_polygons <- st_read('data/EU/eu_nuts3_2024.gpkg') %>% mutate(nuts = NUTS_ID)
# {https://www.dropbox.com/scl/fi/kr61s2lurcsia3jmldjo7/eu_nuts3_2024.gpkg?rlkey=au9iwm8gj0av6xvx4qdeevf4k&dl=0}
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
# {https://www.dropbox.com/scl/fi/8geqx7i7q4m90snvdgmsh/geoBoundariesCGAZ_ADM0.gpkg?rlkey=81mz0k6nbnhbp1o0uhb8xtdhj&dl=0}

temp_plot_eu <- ggplot() +
  geom_sf(data = world, aes(geometry = geom), fill = 'gray95', color = NA) +
  geom_sf(data = eu_polygons_multiples %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(side == 'hot'), aes(geometry = geom, fill = percent_pred*7), color = 'gray85', linewidth = 0.05) +
  scale_fill_gradientn(labels = scales::percent_format(), name = '% change, heat', limits = c(0, 0.5), oob = scales::squish, na.value = 'white', breaks = c(0, 0.25, 0.5), values = scales::rescale(seq(0, 0.5, length.out = 6)), colors = c('#eaeccc','#FEDA8b','#fdb366','#f67e4b','#dd3d2d','#a50026')) +
  new_scale_fill() +
  geom_sf(data = eu_polygons_multiples %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(side == 'cold'), aes(geometry = geom, fill = percent_pred*7), color = 'gray85', linewidth = 0.05) +
  scale_fill_gradientn(labels = scales::percent_format(), name = '% change, cold', limits = c(0, 0.5), oob = scales::squish, na.value = 'white', breaks = c(0, 0.25, 0.5), values = scales::rescale(seq(0, 0.5, length.out = 6)), colors = c('#eaeccc','#c2e4ef','#98cae1','#6ea6cd','#4a7bb7','#364b9a')) +
  facet_wrap(facets = vars(as.factor(temp)), nrow = 3) + theme_clean() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line = element_blank()) + ggtitle('Change in daily mortality rate relative to temperature of minimum risk') + coord_sf(xlim = c(-14,31), ylim = c(33,75)) + theme(aspect.ratio = 1)
ggsave(temp_plot_eu, file = 'fig/combined/county_maps_eu.png', width = 10, height = 8, units = 'in')

temp_plot_eu_030 <- ggplot() +
  geom_sf(data = world, aes(geometry = geom), fill = 'gray95', color = NA) +
  geom_sf(data = eu_polygons_multiples %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(side == 'cold', temp %in% c('0ºC day','25ºC day')), aes(geometry = geom, fill = percent_pred*7), color = 'gray85', linewidth = 0.05) +
  scale_fill_gradientn(labels = scales::percent_format(), name = '% change, cold', limits = c(0,.6), oob = scales::squish, na.value = 'white', breaks = c(0, 0.3, 0.6), values = scales::rescale(seq(0, 0.6, length.out = 6)), colors = c('#eaeccc','#c2e4ef','#98cae1','#6ea6cd','#4a7bb7','#364b9a'), guide = guide_colorbar(order = 1)) +
  new_scale_fill() +
  geom_sf(data = eu_polygons_multiples %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(side == 'hot', temp %in% c('0ºC day','25ºC day')), aes(geometry = geom, fill = percent_pred*7), color = 'gray85', linewidth = 0.05) +
  scale_fill_gradientn(labels = scales::percent_format(), name = '% change, heat', limits = c(0,0.2), oob = scales::squish, na.value = 'white', breaks = c(0, 0.1, 0.2), values = scales::rescale(seq(0, 0.2, length.out = 6)), colors = c('#eaeccc','#FEDA8b','#fdb366','#f67e4b','#dd3d2d','#a50026'), guide = guide_colorbar(order = 2)) +
  facet_wrap(facets = vars(as.factor(temp)), ncol = 2) + theme_clean() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line = element_blank()) + ggtitle('Change in daily mortality rate relative to temperature of minimum risk') + coord_sf(xlim = c(-14,31), ylim = c(33,75)) + theme(aspect.ratio = 1)
ggsave(temp_plot_eu_030, file = 'fig/combined/county_maps_eu_030.png', width = 10, height = 5, units = 'in')

eu_totals_panel <- CJ(fips = unique(county_totals$fips), side = c('hot','cold'))
county_totals_panel %<>% merge(county_totals, by = c('fips','side'), all.x = T) %>% mutate(deaths = fifelse(is.na(deaths), 0, deaths))

#### U.S. --- mortality data for the U.S. is not publicly available, so it cannot be provided for this replication
month_panel <- read_parquet('----') %>% select(fips:precip_d2_l0)

month_panel %<>% mutate(mort_rate = deaths_total/pop_total, homicide_rate = homicide/pop_total, log_pc_income = log(pc_income), log_bmr = log(bmr))

n_bins <- 5
reg_df <- month_panel %>%
  mutate(mort_rate = ifelse(mort_rate > fnth(mort_rate, 0.99), fnth(mort_rate, 0.99), mort_rate)) %>%
  mutate(avg_temp_bins = factor(cut(avg_temp, c(fmin(avg_temp)-1, lapply(1:(n_bins-1), function(x) {fnth(avg_temp, x/n_bins, na.rm = T)}), fmax(avg_temp)+1), include.lowest = T, labels = F)), pc_income_bins = factor(cut(log_pc_income, c(fmin(log_pc_income)-1, lapply(1:(n_bins-1), function(x) {fnth(log_pc_income, x/n_bins, na.rm = T)}), fmax(log_pc_income)+1), include.lowest = T, labels = F))) %>% as.data.table()

reg_df %<>% group_by(year) %>% mutate(weight = pop_total/sum(pop_total)) %>% as.data.table()

int_reg <- feols(homicide_rate ~
                   temp_d1_l0 + temp_d2_l0 + temp_d3_l0 + temp_d4_l0 +
                   (temp_d1_l0 + temp_d2_l0 + temp_d3_l0 + temp_d4_l0):avg_temp +
                   (temp_d1_l0 + temp_d2_l0 + temp_d3_l0 + temp_d4_l0):log_pc_income +
                   (temp_d1_l0 + temp_d2_l0 + temp_d3_l0 + temp_d4_l0):log_bmr +
                   precip_d1_l0 + precip_d2_l0
                 | fips^year + fips^month, data = reg_df, se = 'standard')
summary(int_reg)

int_reg <- coef(int_reg)[grepl('temp_',names(coef(int_reg)))]

avg_temps <- seq(5,20,5)
log_pc_incomes <- seq(9,11.5,0.5)
log_bmrs <- seq(-5.5,-4,0.5)

plot_df <- list()
for (this_avg_temp in avg_temps) {
  for (this_log_pc_income in log_pc_incomes) {
    for (this_log_bmr in log_bmrs) {
      preds <-
        pred_temps%*%int_reg[1:4] +
        this_avg_temp*pred_temps%*%int_reg[5:8] +
        this_log_pc_income*pred_temps%*%int_reg[9:12] +
        this_log_bmr*pred_temps%*%int_reg[13:16]
      # se <- pred_temps%*%int_reg[1:4] +
      #   this_avg_temp*pred_temps%*%int_reg[5:8] +
      #   this_log_pc_income*pred_temps%*%int_reg[9:12] +
      #   this_log_bmr*pred_temps%*%int_reg[13:16]
      which_mmt <- ifelse(pred_temps_vector[pred_temps_vector >= 0 & pred_temps_vector <= 30][which.min(preds[pred_temps_vector >= 0 & pred_temps_vector <= 30])] > 0 & pred_temps_vector[pred_temps_vector >= 0 & pred_temps_vector <= 30][which.min(preds[pred_temps_vector >= 0 & pred_temps_vector <= 30])] < 30,  pred_temps_vector[pred_temps_vector >= 0 & pred_temps_vector <= 30][which.min(preds[pred_temps_vector >= 0 & pred_temps_vector <= 30])], 0)
      preds <- preds - preds[pred_temps_vector == which_mmt]
      plot_df[[length(plot_df)+1]] <- data.table(avg_temp = this_avg_temp, income = round(exp(this_log_pc_income)), bmr = round(exp(this_log_bmr), digits = 4)*1e5, pred = c(preds), temp = pred_temps_vector)
    }
  }
}
plot_df <- rbindlist(plot_df)

ggplot(plot_df %>% filter(temp > 0, temp < 35), aes(x = temp, y = pred, color = as.factor(income))) + geom_line() + facet_grid(rows = vars(avg_temp), cols = vars(bmr)) + theme_minimal()

pred_temps_vector <- seq(-50,50,0.1)
pred_temps <- poly(pred_temps_vector, degree = 4, raw = T)
fips_list <- unique(month_panel$fips)
county_values <- list()
county_totals <- list()
for (this_fips in fips_list) {
  
  cat(this_fips, '\n')
  
  this_county_weather <- weather[fips == this_fips]
  
  for_2023 <- reg_df %>% filter(fips == this_fips, year == 2023, month == 1) %>% select(pc_income, avg_temp, sd_temp, bmr, avg_age, log_pc_income, log_bmr)
  for_2023$typical_rate <- reg_df %>% filter(fips == this_fips, year > 2013) %>% group_by(year) %>% summarize(rate_total = sum(mort_rate)/365.25, .groups = 'keep') %>% ungroup() %>% summarize(rate_total = mean(rate_total), .groups = 'keep') %>% pull(rate_total)
  preds <-
    pred_temps%*%int_reg[1:4] +
    for_2023$avg_temp*pred_temps%*%int_reg[5:8] +
    for_2023$log_pc_income*pred_temps%*%int_reg[9:12] +
    for_2023$log_bmr*pred_temps%*%int_reg[13:16]
  mmt_options <- which(pred_temps_vector > for_2023$avg_temp - for_2023$sd_temp & pred_temps_vector < for_2023$avg_temp + for_2023$sd_temp & pred_temps_vector < 30 & pred_temps_vector > -10)
  mmt <- ifelse(pred_temps_vector[mmt_options][which.min(preds[mmt_options])] > pred_temps_vector[mmt_options[1]] & pred_temps_vector[mmt_options][which.min(preds[mmt_options])] < pred_temps_vector[mmt_options[length(mmt_options)]],  pred_temps_vector[mmt_options][which.min(preds[mmt_options])], pred_temps_vector[mmt_options[1]])
  # sides <- c(rep('cold',(min(mmt_options)-1) + which.min(preds[mmt_options])),rep('hot',length(pred_temps_vector) - ((min(mmt_options)-1) + which.min(preds[mmt_options]))))
  preds <- preds - preds[pred_temps_vector == mmt]
  pred_df <- data.table(temp = pred_temps_vector, pred = c(preds))
  
  county_values[[length(county_values)+1]] <- data.table(fips = this_fips,
                                                         temp = pred_temps_vector,
                                                         pred = c(preds),
                                                         avg_temp = for_2023$avg_temp,
                                                         # sd_temp = for_2023$sd_temp,
                                                         pc_income = for_2023$pc_income,
                                                         # bmr = for_2023$bmr,
                                                         # avg_age = for_2023$avg_age,
                                                         percent_pred = c(preds)
                                                         # side = sides
  )
}
county_values <- rbindlist(county_values)

min_max_temps <- open_dataset('data/US/us_temp.pq', format = 'parquet') %>% group_by(poly_id) %>% summarize(min_temp = min(order_1, na.rm = T), max_temp = max(order_1, na.rm = T)) %>% collect()

county_values %<>% merge(min_max_temps %>% mutate(fips = as.integer(poly_id)) %>% select(-poly_id), by = 'fips') %>% mutate(exclude = ifelse(temp < min_temp - 1 | temp > max_temp + 1, 1, 0))

county_polygons <- st_read('data/US/CONUS_consistent_counties.gpkg') %>% rename(fips = GEOID) %>% mutate(fips = as.integer(fips))
# {https://www.dropbox.com/scl/fi/ls3dqkf8xfpxqjbvp4ksc/CONUS_consistent_counties.gpkg?rlkey=7oftm4mxve65vqe5jlrxok1aj&dl=0}

county_polygons_multiples <- rbind(county_polygons %>% merge(county_values %>% filter(temp == 35), by = 'fips', all.x = T),
                                   county_polygons %>% merge(county_values %>% filter(temp == 30), by = 'fips', all.x = T),
                                   county_polygons %>% merge(county_values %>% filter(temp == 25), by = 'fips', all.x = T),
                                   county_polygons %>% merge(county_values %>% filter(temp == 20), by = 'fips', all.x = T),
                                   county_polygons %>% merge(county_values %>% filter(temp == 15), by = 'fips', all.x = T),
                                   county_polygons %>% merge(county_values %>% filter(temp == 10), by = 'fips', all.x = T),
                                   county_polygons %>% merge(county_values %>% filter(temp == 5), by = 'fips', all.x = T),
                                   county_polygons %>% merge(county_values %>% filter(temp == 0), by = 'fips', all.x = T),
                                   county_polygons %>% merge(county_values %>% filter(temp == -10), by = 'fips', all.x = T))

county_polygons_multiples %<>% mutate(temp_lab = paste0(as.character(temp), 'ºC day')) %>% mutate(temp = factor(temp, levels = rev(unique(temp))))

temp_plot_bottomcode <- ggplot() +
  geom_sf(data = county_polygons_multiples %>% filter(temp %in% c(35,30,25,20,10,0)) %>% mutate(pred = ifelse(exclude == 1, NA, pred)), aes(geometry = geom, fill = pred*1e5/0.01576786), color = 'gray80', linewidth = 0.1) +
  scale_fill_gradientn(labels = scales::percent_format(), 'Change in homicide rate (%)', limits = c(-.75,.75), oob = scales::squish, na.value = 'white', breaks = c(-.75,-0.5,0,0.5,.75), values = scales::rescale(seq(-.75, .75, length.out = 5)), colors = c('#364b9a','#98cae1','#eaeccc','#fdb366','#a50026')) +
  facet_wrap(facets = vars(as.factor(temp_lab)), nrow = 3) + theme_clean() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line = element_blank()) + ggtitle('Change in daily mortality rate relative to temperature of minimum risk')

temp_plot_bottomcode <- ggplot() +
  geom_sf(data = county_polygons_multiples %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(side == 'hot'), aes(geometry = geom, fill = percent_pred*30.437), color = 'gray80', linewidth = 0.1) +
  scale_fill_gradientn(labels = scales::percent_format(), name = '% change, heat', limits = c(0,.3), oob = scales::squish, na.value = 'white', breaks = c(0, 0.15, 0.3), values = scales::rescale(seq(0, 0.3, length.out = 6)), colors = c('#eaeccc','#FEDA8b','#fdb366','#f67e4b','#dd3d2d','#a50026')) +
  new_scale_fill() +
  geom_sf(data = county_polygons_multiples %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(side == 'cold'), aes(geometry = geom, fill = percent_pred*30.437), color = 'gray80', linewidth = 0.1) +
  scale_fill_gradientn(labels = scales::percent_format(), name = '% change, cold', limits = c(0,.3), oob = scales::squish, na.value = 'white', breaks = c(0, 0.15, 0.3), values = scales::rescale(seq(0, 0.3, length.out = 6)), colors = c('#eaeccc','#c2e4ef','#98cae1','#6ea6cd','#4a7bb7','#364b9a')) +
  facet_wrap(facets = vars(as.factor(temp)), nrow = 3) + theme_clean() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line = element_blank()) + ggtitle('Change in daily mortality rate relative to temperature of minimum risk')
ggsave(temp_plot_bottomcode, file = 'fig/combined/county_maps4.png', width = 10, height = 8, units = 'in')

temp_plot_bottomcode_030 <- ggplot() +
  geom_sf(data = county_polygons_multiples %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(side == 'cold', temp %in% c('0ºC day','30ºC day')), aes(geometry = geom, fill = percent_pred*30.437), color = 'gray80', linewidth = 0.1) +
  scale_fill_gradientn(labels = scales::percent_format(), name = '% change, cold', limits = c(0,.25), oob = scales::squish, na.value = 'white', breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25), values = scales::rescale(seq(0, 0.25, length.out = 6)), colors = c('#eaeccc','#c2e4ef','#98cae1','#6ea6cd','#4a7bb7','#364b9a'), guide = guide_colorbar(order = 1)) +
  new_scale_fill() +
  geom_sf(data = county_polygons_multiples %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(side == 'hot', temp %in% c('0ºC day','30ºC day')), aes(geometry = geom, fill = percent_pred*30.437), color = 'gray80', linewidth = 0.1) +
  scale_fill_gradientn(labels = scales::percent_format(), name = '% change, heat', limits = c(0,.1), oob = scales::squish, na.value = 'white', breaks = c(0, 0.05, 0.1), values = scales::rescale(seq(0, 0.1, length.out = 6)), colors = c('#eaeccc','#FEDA8b','#fdb366','#f67e4b','#dd3d2d','#a50026'), guide = guide_colorbar(order = 2)) +
  facet_wrap(facets = vars(as.factor(temp)), ncol = 2) + theme_clean() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line = element_blank()) + ggtitle('Change in daily mortality rate relative to temperature of minimum risk')
ggsave(temp_plot_bottomcode_030, file = 'fig/combined/county_maps4_030.png', width = 10, height = 5, units = 'in')

hot_temp <- 27
cold_temp <- 5
mex_us_df <- rbind(county_polygons %>% merge(county_values %>% filter(temp == cold_temp), by = 'fips', all.x = T) %>%
                     mutate(temp = paste0(as.character(temp), 'ºC day')) %>% mutate(temp = factor(temp, levels = rev(unique(temp)))) %>%
                     select(adm_id = fips, temp, percent_pred, exclude, geom, side),
                   mex_polygons %>% merge(mex_adm_values %>% filter(temp == cold_temp), by = 'adm2', all.x = T) %>%
                     mutate(temp = paste0(as.character(temp), 'ºC day')) %>% mutate(temp = factor(temp, levels = rev(unique(temp)))) %>%
                     select(adm_id = adm2, temp, percent_pred, exclude, geom = geometry, side),
                   county_polygons %>% merge(county_values %>% filter(temp == hot_temp), by = 'fips', all.x = T) %>%
                     mutate(temp = paste0(as.character(temp), 'ºC day')) %>% mutate(temp = factor(temp, levels = rev(unique(temp)))) %>%
                     select(adm_id = fips, temp, percent_pred, exclude, geom, side),
                   mex_polygons %>% merge(mex_adm_values %>% filter(temp == hot_temp), by = 'adm2', all.x = T) %>%
                     mutate(temp = paste0(as.character(temp), 'ºC day')) %>% mutate(temp = factor(temp, levels = rev(unique(temp)))) %>%
                     select(adm_id = adm2, temp, percent_pred, exclude, geom = geometry, side))

us_mex_plot <- ggplot() +
  geom_sf(data = mex_us_df %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(side == 'cold', temp %in% c(paste0(cold_temp,'ºC day'),paste0(hot_temp,'ºC day'))), aes(geometry = geom, fill = percent_pred*30.437), color = 'gray80', linewidth = 0.1) +
  scale_fill_gradientn(labels = scales::percent_format(), name = '% change, cold', limits = c(0,.25), oob = scales::squish, na.value = 'white', breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25), values = scales::rescale(seq(0, 0.25, length.out = 6)), colors = c('#eaeccc','#c2e4ef','#98cae1','#6ea6cd','#4a7bb7','#364b9a'), guide = guide_colorbar(order = 1)) +
  new_scale_fill() +
  geom_sf(data = mex_us_df %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(side == 'hot', temp %in% c(paste0(cold_temp,'ºC day'),paste0(hot_temp,'ºC day'))), aes(geometry = geom, fill = percent_pred*30.437), color = 'gray80', linewidth = 0.1) +
  scale_fill_gradientn(labels = scales::percent_format(), name = '% change, heat', limits = c(0,.1), oob = scales::squish, na.value = 'white', breaks = c(0, 0.05, 0.1), values = scales::rescale(seq(0, 0.1, length.out = 6)), colors = c('#eaeccc','#FEDA8b','#fdb366','#f67e4b','#dd3d2d','#a50026'), guide = guide_colorbar(order = 2)) +
  facet_wrap(facets = vars(as.factor(temp)), ncol = 2) + theme_clean() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line = element_blank()) + ggtitle('Change in daily mortality rate relative to temperature of minimum risk')
us_mex_plot

scale_breaks <- c(0, 0.1, 0.2)
us_mex_plot_heat <- ggplot() +
  geom_sf(data = mex_us_df %>% mutate(percent_pred = ifelse(exclude == 1, NA, percent_pred)) %>% filter(temp == paste0(hot_temp,'ºC day')), aes(geometry = geom, fill = percent_pred*30.437), color = 'gray80', linewidth = 0.05) +
  scale_fill_gradientn(labels = c(paste0(100*scale_breaks[1:(length(scale_breaks)-1)], '%'), paste0('>',100*scale_breaks[length(scale_breaks)],'%')), name = '% change', limits = c(0,scale_breaks[length(scale_breaks)]), oob = scales::squish, na.value = 'white', breaks = scale_breaks, values = scales::rescale(seq(0, scale_breaks[length(scale_breaks)], length.out = 3)), colors = c('gray95','#dd3d2d','#520000')) +
  theme_clean() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line = element_blank()) + ggtitle(paste0('Change in daily mortality rate at ', hot_temp, 'ºC,\nrelative to temperature of minimum risk'))
ggsave(us_mex_plot_heat, file = 'fig/combined/county_maps_heat.png', width = 10, height = 5, units = 'in')