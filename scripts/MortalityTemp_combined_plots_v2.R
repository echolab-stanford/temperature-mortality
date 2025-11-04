################################################################################
# Last Updated By: Khusel Avirmed 
# Date: 2/9/2025
################################################################################

rm(list = ls())
gc()

pacman::p_load(fastverse, tidyverse, arrow, MetBrewer, cowplot, patchwork)
setwd("~/BurkeLab Dropbox/projects/temperature-mortality-descriptive")
script_dir <- 'script'

# load(paste(processed_country_dir, paste0(country_name, "AnnualBinned_person_day_exposure.RData"), sep = "/"))

# # bring in helper functions to perform analysis
# source(paste(script_dir, "helper_functions.R", sep = "/"))

read_country_results <- function(country_name, si_folder = NULL){
  
  # country_name <- "US"
  # si_folder <- NULL # log_rate, same_year, age_standardized, or NULL
  # country_name <- country_name
  
  if (!is.null(si_folder)) {
    processed_country_dir <- paste("processed", country_name, sep = "/")
    processed_country_dir <- file.path(processed_country_dir, si_folder)
  } else {
    processed_country_dir <- paste("processed", country_name, sep = "/")
  }
  
  print(processed_country_dir)
  pooled_result <- readRDS(  paste(processed_country_dir, paste0(country_name, "_Pooled_results.rds"), sep = "/"))
  young_result <- readRDS(  paste(processed_country_dir, paste0(country_name, "_Young_results.rds"), sep = "/"))
  adult_result <- readRDS(  paste(processed_country_dir, paste0(country_name, "_Adult_results.rds"), sep = "/"))
  elderly_result <- readRDS(  paste(processed_country_dir, paste0(country_name, "_Elderly_results.rds"), sep = "/"))
  temp_result <- readRDS(  paste(processed_country_dir, paste0(country_name, "_temp_results.rds"), sep = "/"))
  income_result <- readRDS(  paste(processed_country_dir, paste0(country_name, "_income_results.rds"), sep = "/"))
  decade_result <- readRDS(  paste(processed_country_dir, paste0(country_name, "_decade_results.rds"), sep = "/"))
  baseline_result <- readRDS(  paste(processed_country_dir, paste0(country_name, "_baseline_results.rds"), sep = "/"))
  
  all_response <- rbind(
    pooled_result$response, 
    young_result$response,
    adult_result$response,
    elderly_result$response, 
    temp_result$response, 
    income_result$response, 
    decade_result$response, 
    baseline_result$response
  )
  all_response[, country := country_name]
  
  all_deaths <- rbind(
    pooled_result$deaths, 
    young_result$deaths,
    adult_result$deaths,
    elderly_result$deaths,
    temp_result$deaths, 
    income_result$deaths, 
    decade_result$deaths, 
    baseline_result$deaths
  )
  all_deaths[, country := country_name]
  
  deaths_by_year_side <- rbind(
    pooled_result$deaths_by_year_side, 
    young_result$deaths_by_year_side,
    adult_result$deaths_by_year_side,
    elderly_result$deaths_by_year_side,
    temp_result$deaths_by_year_side, 
    income_result$deaths_by_year_side, 
    decade_result$deaths_by_year_side, 
    baseline_result$deaths_by_year_side
  )
  deaths_by_year_side[, country := country_name]
  
  return(list(response = all_response, deaths = all_deaths, deaths_by_year_side = deaths_by_year_side))
}

us <- read_country_results(country_name = "US", si_folder = NULL)
mex <- read_country_results(country_name = "MEX", si_folder = NULL)
eu <- read_country_results(country_name = "EU", si_folder = NULL)

combined_response <- rbind(us$response, mex$response, eu$response)
combined_deaths <- rbind(us$deaths, mex$deaths, eu$deaths)
combined_deaths_by_year_side <- rbind(us$deaths_by_year_side, mex$deaths_by_year_side, eu$deaths_by_year_side)

combined_deaths_by_year_side %>%
  filter(model == "Pooled") %>%
  group_by(country) %>%
  filter(year %in% tail(sort(unique(year)), 10)) %>%
  group_by(country, side) %>%
  summarise(
    temp_deaths_10yrs = mean(deaths),
    total_deaths_10yrs = mean(total_deaths),
    .groups = "drop"
  ) %>% as.data.table() 

combined_deaths_by_year_side %>%
  filter(model == "Pooled") %>%
  group_by(country) %>%
  filter(year == max(year)) %>% 
  ungroup() %>% 
  select(year, side,temp_deaths = deaths, total_deaths, country)

###################################################################################
#   MAIN FIGURE 1: Pooled Response & Pooled Percentage Attributable Deaths by Side
###################################################################################
plot_pooled_fig <- function(pooled_dataframe, deaths_dataframe, si_folder = NULL) {
  pooled_dataframe <- combined_response_log_rate
  deaths_dataframe <- combined_death_log_rate
  # si_folder = "log_rate"
  
  # pooled_dataframe <- combined_response
  # deaths_dataframe <- combined_deaths
  # si_folder = NULL
  
  pooled_response_df <- pooled_dataframe %>% filter(model == "Pooled")
  pooled_deaths_df <- deaths_dataframe %>% filter(model == "Pooled")
  
  pooled_response_df <- pooled_response_df %>% 
    mutate(country = factor(country, levels = c("US", "EU", "MEX")))
  # mutate(country = paste0(country, " (MMT = ", mmt, "°C)"),
  #        country = factor(country, levels = c("US (MMT = 26°C)", "EU (MMT = 20°C)", "MEX (MMT = 24°C)"))) %>%
  
  if (!is.null(si_folder) && si_folder == "log_rate") {
    pooled_response_df <- pooled_response_df %>% group_by(country) %>% 
      # filter(bins >= t1 & bins <= t99) %>% 
      mutate(
        pred_q025 = expm1(pred_q025),
        pred = expm1(pred),
        pred_q975 = expm1(pred_q975)
      ) %>% 
      ungroup() %>% 
      mutate(country = factor(country, level = c("US", "EU", "MEX")))
    
    y_lab_text <- "Percent change in deaths"
    
  } else if (!is.null(si_folder) && si_folder == "same_year") {
    pooled_response_df <- pooled_response_df %>% 
      mutate(
        share_pop = round(pop_deg_days / sum(pop_deg_days), 4)
      ) %>% 
      filter(share_pop != 0) %>% 
      mutate(
        pred_q025 = pred_q025 * 1e+5,
        pred = pred * 1e+5,
        pred_q975 = pred_q975 * 1e+5
      )
    
    y_lab_text <- "Change in deaths \n per 100,000 people"
    
  } else {  # Handles both NULL and any other unexpected value
    pooled_response_df <- pooled_response_df %>% 
      group_by(country) %>% 
      # filter(bins >= t1 & bins <= t99) %>% 
      mutate(
        pred_q025 = pred_q025 * 1e+5,
        pred = pred * 1e+5,
        pred_q975 = pred_q975 * 1e+5
      ) %>% 
      ungroup() %>% 
      mutate(country = factor(country, level = c("US", "EU", "MEX")))
    
    
    y_lab_text <- "Change in deaths \n per 100,000 people"
  }
  
  
  # Define the clean theme
  theme_clean <- function() {
    theme_classic() +
      theme(
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        plot.title = element_text(color = "black", size = 12, hjust = 0.5, face = "bold"),
        axis.line = element_line(color = 'black', linewidth = 0.35),
        axis.ticks = element_line(colour = "black", linewidth = 0.35),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        strip.text = element_text(size = 10, face = "bold")
      )
  }
  
  unique_levels <- sort(unique(pooled_response_df$country))
  color_palette <- c("#009988", "#EE7733", "#0072B2")
  # color_palette <- rev(met.brewer("Homer1", length(unique_levels)))
  # color_palette <- setNames(c("#1f77b4", "#ff7f0e", "#2ca02c"), unique_levels)
  
  # Create the response plot
  response <- pooled_response_df %>% 
    filter(pop_deg_days != 0 & !is.na(pop_deg_days)) %>% filter(bins >= t1) %>% 
    ggplot(aes(x = bins)) +
    geom_line(aes(y = pred, color = country), linewidth = 1) +
    geom_ribbon(aes(ymin = pred_q025, ymax = pred_q975, fill = country), alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    # geom_col(aes(y = -(pop_deg_days*100)*0.2, fill = country), position = "identity",color = "black", linetype = "solid", alpha = 0.3, width = 1) + 
    labs(x = "Temperature (°C)", y = y_lab_text) +
    scale_color_manual(values = color_palette) +
    scale_fill_manual(values = color_palette) +
    theme_clean() +
    # facet_grid(rows = vars(model_name), cols = vars(country), scales = "free_y") +
    theme(
      strip.text.x = element_text(size = 10, face = "bold"),
      strip.text.y = element_blank(),
      axis.text.x = element_blank(), 
      axis.title.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_line(color = "gray80", linetype = "dotted", linewidth = 0.5), 
      panel.grid.minor.x = element_blank(), 
      panel.grid.major.y = element_line(color = "gray90", linetype = "solid", linewidth = 0.5), 
      panel.grid.minor.y = element_blank(),
      legend.background = element_rect(fill = "transparent", colour = NA), 
      legend.position = "inside", 
      legend.position.inside = c(0.2, 0.7),
      legend.title = element_blank()
    )
  
  exposure <- pooled_response_df %>% 
    filter(pop_deg_days != 0 & !is.na(pop_deg_days)) %>% filter(bins >= t1) %>% 
    ggplot(aes(x = bins, weight = pop_deg_days, fill = country, color = country)) +
    geom_density(alpha = 0.3, bw = 1) + 
    scale_color_manual(values = setNames(color_palette, unique_levels)) +
    scale_fill_manual(values = setNames(color_palette, unique_levels)) +
    labs(x = "Temperature (°C)", y = " \n ") +
    theme_minimal() + 
    theme(
      legend.position = "none",
      axis.title.y = element_blank(), 
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(t = 0, r = 5, b = 0, l = 20)
    )
  
  plot_grid(response, exposure, ncol = 1)
  final_plot <- plot_grid(response, exposure, ncol = 1, rel_heights = c(2, 0.5), align = 'vh')
  
  breaks <- c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14)
  deaths <- pooled_deaths_df %>% 
    mutate(annual_share_total_mortality = deaths / total_deaths ) %>% 
    mutate(country = factor(country, levels = c("US", "EU", "MEX"))) %>%
    ggplot(aes(x = country, y = annual_share_total_mortality, fill = side)) + 
    geom_bar(stat = 'identity', position = 'stack', width = 0.9, alpha = 1) + 
    scale_fill_manual(values = c("cold" = "#D3D3D3", "hot" = "#FFCC99"), labels = c("cold" = "Cold", "hot" = "Hot"), name = 'Mortality due to ...') + 
    theme_clean() +
    scale_y_continuous(
      breaks = breaks,
      labels = scales::percent_format()
    ) + 
    ylab('Percent of total annual deaths') +
    theme(
      strip.text.x = element_text(size = 10, face = "bold"),
      strip.text.y = element_blank(),
      axis.title.x = element_blank()
    )
  
  plot_grid(final_plot, deaths, ncol = 2)
}

fig_1 <- plot_pooled_fig(combined_response, combined_deaths, si_folder = NULL)
fig_1_dir <- file.path("fig/combined/new_approach/", paste0("fig_1.pdf"))
ggsave(plot = fig_1, filename = fig_1_dir, width = 10, height = 5, units = "in", dpi = 300)

plot_response <- function(response_df, si_folder = NULL) {
  
  # response_df <- combined_response %>% filter(country == "US") %>% filter(model != "Pooled")
  # si_folder <- NULL
  # si_folder = "log_rate"
  # response_df <- combined_response_age_standardized %>% filter(country == "US") %>% filter(model != "Pooled" & model != "Adult" & model != "Elderly" & model != "Young")
  response_df <- response_df %>% filter(!is.na(pop_deg_days))
  # unique(response_df$model_name)
  # unique(response_df$model_level)
  # response_df <- all_response %>% filter(country == "MEX" & model != "Pooled")
  # response_df <- response_df %>% mutate(t99 = t99+2)
  
  response_df <- response_df %>% 
    # mutate(pop_deg_days = round(pop_deg_days, 5)) %>%
    # filter(pop_deg_days != 0 ) %>%
    # mutate(max_pop = max(pop_deg_days) + 0.02,
    #        pop_deg_days = pop_deg_days * max_pop) %>% 
    filter(key %in% c("pooled", "1970", "1980", "1990", "2000", "2010", '2020', '10%', '50%', '90%')) %>%
    mutate(group_level = case_when(
      model == 'Pooled' ~ "Pooled",
      model == 'temp' ~ "Climate",
      model == 'income' ~ "Income",
      model == 'baseline' ~ "Baseline Mortality Rate",
      model == 'Young' ~ "Age", 
      model == 'Adult' ~ "Age", 
      model == 'Elderly' ~ "Age", 
      model == 'decade' ~ "Decade",
      .default = as.character(model)
    )) %>% 
    mutate(key = paste(model, key))
  
  country_name <- unique(response_df$country)
  print(country_name)
  if (is.null(si_folder)){
    response_df <- response_df %>% 
      mutate(pred_q025 = pred_q025*1e+5,
             pred = pred*1e+5,
             pred_q975 = pred_q975*1e+5)
    y_lab_text <- "Change in deaths \n per 100,000 people"
  } else if (si_folder == "log_rate") {
    response_df <- response_df %>%
      mutate(pred_q025 = expm1(pred_q025),
             pred = expm1(pred),
             pred_q975 = expm1(pred_q975))
    if (country_name == "EU") {
      response_df <- response_df %>%
        mutate(pred_q025 = pred_q025/4.2,
               pred = pred/4.2,
               pred_q975 = pred_q975/4.2)
    }
    y_lab_text <- "Percent change in\nmonthly mortality rate"
  } else if (si_folder == "age_standardized") {
    response_df <- response_df %>% 
      mutate(pred_q025 = pred_q025*1e+5,
             pred = pred*1e+5,
             pred_q975 = pred_q975*1e+5)
    y_lab_text <- "Change in deaths \n per 100,000 people"
  }
  
  # if (!is.null(si_folder)) {
  #   response_df <- response_df %>%
  #     mutate(pred_q025 = expm1(pred_q025),
  #          pred = expm1(pred),
  #          pred_q975 = expm1(pred_q975))
  #   y_lab_text <- "Percent change in deaths"
  # } else {
  #   response_df <- response_df %>% 
  #     mutate(pred_q025 = pred_q025*1e+5,
  #            pred = pred*1e+5,
  #            pred_q975 = pred_q975*1e+5)
  #   y_lab_text <- "Change in deaths \n per 100,000 people"
  # }
  # 
  # Split data by level_new and create plots
  response_df_list <- response_df %>%
    mutate(key = case_when(
      key == 'Pooled pooled' ~ "Pooled",
      key == 'Pooled (Lag) pooled' ~ "Pooled lagged",
      key == 'Pooled (lag) pooled' ~ "Pooled lagged",
      key == 'temp 10%' ~ "Temp 10th",
      key == 'temp 50%' ~ "Temp 50th",
      key == 'temp 90%' ~ "Temp 90th",
      key == 'income 10%' ~ "Income 10th",
      key == 'income 50%' ~ "Income 50th",
      key == 'income 90%' ~ "Income 90th",
      key == 'baseline 10%' ~ "Baseline 10th",
      key == 'baseline 50%' ~ "Baseline 50th",
      key == 'baseline 90%' ~ "Baseline 90th",
      key == "Adult pooled" ~ "Adult", 
      key == "Elderly pooled" ~ "Elderly", 
      key == "Young pooled" ~ "Young",
      key == 'decade 1970' ~ "1970s",
      key == 'decade 1980' ~ "1980s",
      key == 'decade 1990' ~ "1990s",
      key == 'decade 2000' ~ "2000s",
      key == 'decade 2010' ~ "2010s",
      key == 'decade 2020' ~ "2020s",
      .default = as.character(key)
    )) %>%
    mutate(group_level = factor(group_level, levels = c("Pooled", "Age", "Climate", "Income", "Baseline Mortality Rate", "Decade"))) %>% 
    group_split(group_level)
  
  plot_list <- lapply(response_df_list, function(df){
    # df <- response_df_list[[5]]
    # df <- response_df_list[[2]] %>% filter(model_name == "Young")
    # df <- response_df_list[[1]]
    
    group_level <- c(unique(df$group_level))
    country <- unique(df$country)
    
    if (group_level == "Age"){
      df <- df %>% mutate(key = factor(key,  levels = c("Young", "Adult", "Elderly"))) %>% 
        arrange(bins, key)
    } else if (group_level == "Decade") {
      df <- df %>% mutate(key = factor(key,  levels = c("1970s", "1980s", "1990s", "2000s", "2010s", "2020s"))) %>% 
        arrange(key)
    }
    
    
    if (group_level == "Age") {
      axis_title_y_element <- element_text(size = 8, family = 'Helvetica')
    } else {
      axis_title_y_element <- element_blank()
    }
    
    if (si_folder == "age_standardized") {
      if (group_level == "Climate") {
        axis_title_y_element <- element_text(size = 8, family = 'Helvetica')
      } else {
        axis_title_y_element <- element_blank()
      }
    }
    
    if (country == "US"){
      strip_text_element <- element_text(size = 10, face = "bold", family = 'Helvetica')
    } else {
      strip_text_element <- element_blank()
    }
    
    
    unique_levels <- unique(df$key)
    n_length <- length(unique_levels)
    color_palette <- c("#a62f00", "#df7700", "#32b2da", "#6ad5e8", "#fff179")
    color_palette <- rev(color_palette[1:n_length])
    
    # df <- df %>% mutate(share_pop = pop_deg_days / sum(pop_deg_days)) %>%
    #   mutate(share_pop = round(share_pop, 4)) %>%
    #   filter(share_pop != 0)
    # if (country == "US") {
    #   df <- df %>% filter(bins > -20 & bins < 30)
    # } else if (country == "MEX") {
    #   df <- df %>% filter(bins > 0 & bins < 40)
    # }

    df <- df %>% filter(pop_deg_days != 0 & !is.na(pop_deg_days)) %>% filter(bins >= t1)
    # df <- df %>% filter(bins >= t1 & bins <= t99)
    response <- df %>% 
      ggplot(aes(x = bins)) +
      geom_line(aes(y = pred, color = key), linewidth = 1) +
      geom_ribbon(aes(ymin = pred_q025, ymax = pred_q975, fill = key), alpha = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed") +

      scale_color_manual(values = setNames(color_palette, unique_levels)) +
      scale_fill_manual(values = setNames(color_palette, unique_levels)) +
      
      labs(x = "Temperature (°C)", y = y_lab_text) +
      facet_wrap(~group_level, scales = "free_y") +
      theme_minimal() +
      theme(
        strip.text = strip_text_element, 
        axis.title.y = axis_title_y_element, 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "inside", 
        legend.position.inside = c(0.7, 0.8), 
        legend.background = element_rect(fill = "transparent", colour = NA), 
        legend.title = element_blank(),
        legend.text = element_text(size = 8, family = 'Helvetica'),
        legend.key.size = unit(0.5, 'lines'), 
        legend.key.spacing.y = unit(0.1, 'lines'),
        plot.margin = margin(t = 0, r = 5, b = 0, l = 10)
      ) + 
      coord_cartesian(ylim = c(0, NA))
    
    if (si_folder == "log_rate") {
      response <- response + scale_y_continuous(labels = scales::percent_format())
    }
    
    # exposure_df <- pwd[pop_deg_days != 0, .(pop_deg_days = mean(pop_deg_days, na.rm = TRUE)), by = .(bins, decade)]
    # exposure_df[, bins := ceiling(as.numeric(as.character(bins)))] 
    # exposure_df[, decade := factor(paste0(decade, "s"), levels = c("1970s", "1980s", "1990s", "2000s", "2010s"))]
    # 
    exposure <- df %>%
      ggplot(aes(x = bins, weight = pop_deg_days, fill = key, color = key)) +
      geom_density(alpha = 0.3, bw = 1) +
      scale_color_manual(values = setNames(color_palette, unique_levels)) +
      scale_fill_manual(values = setNames(color_palette, unique_levels)) +
      labs(x = "Temperature (°C)", y = " \n ") +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.title.y = axis_title_y_element,
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(t = 0, r = 5, b = 0, l = 20)
      )
    
    
    # final_plot <- response/exposure + patchwork::plot_layout(heights = c(1, 0.2))
    
    final_plot <- plot_grid(response, exposure, ncol = 1, rel_heights = c(1, 0.2))

    return(final_plot)
    })

  # Arrange all facets in a single plot
  ncol <- length(plot_list)
  final_plot <- plot_grid(plotlist = plot_list, ncol = ncol, align = 'hv')
  
  # Print final plot
  print(final_plot)
  
  return(final_plot)
}
plot_annual_deaths <- function(data) {
  
  # data <- all_deaths %>% filter(country == "US" & model != "Pooled")
  # data <- us$annual_deaths
  # data <- deaths
  # data <- us$deaths
  
  data <- data %>% 
    # filter(key %in% c("pooled", "1970", "1980", "1990", "2000", "2010", '2020', '10%', '50%', '90%')) %>%
    mutate(group_level = case_when(
      model == 'Pooled' ~ "Pooled",
      model == 'temp' ~ "Climate",
      model == 'income' ~ "Income",
      model == 'baseline' ~ "Baseline Mortality Rate",
      model == 'Young' ~ "Age", 
      model == 'Adult' ~ "Age", 
      model == 'Elderly' ~ "Age", 
      model == 'decade' ~ "Decade",
      .default = as.character(model)
    ))
  
  pooled_result_deaths <- data %>% 
    filter(group_level == "Pooled" | group_level == "Age" | group_level == "Decade")  %>% 
    mutate(key = paste(model, key))
  
  cont_result_deaths <- data %>% 
    filter(group_level != "Age" & group_level != "Decade" & group_level != "Pooled") %>% 
    group_by(model) %>% 
    mutate(percentile = as.numeric(sub("%", "", key))) %>% 
    mutate(percentile = ntile(percentile, n = 3)) %>% 
    mutate(key = case_when(
      is.na(percentile) ~ key,  
      percentile == 1 ~ "Low", 
      percentile == 2 ~ "Medium", 
      percentile == 3 ~ "High"
    )) %>% 
    ungroup() %>% 
    mutate(key = paste(key, model)) %>% 
    group_by(side, key) %>% 
    summarise(group_level = first(group_level),
              country = first(country),
              model = first(model), 
              across(contains("deaths"), ~sum(.x, na.rm = TRUE)),  
              # annual_total_pop = sum(annual_total_pop, na.rm = TRUE),
              # annual_total_death = sum(annual_total_death, na.rm = TRUE),
              .groups = "drop") %>% 
    select(side, total_deaths, deaths_q025, deaths, deaths_q975, key, model, country, group_level)
  
  factor_order <- c("Pooled pooled", 
                    "Young pooled", "Adult pooled", "Elderly pooled", 
                    "Low temp", "Medium temp", "High temp",
                    "Low income", "Medium income", "High income", 
                    "Low baseline", "Medium baseline", "High baseline", 
                    "decade 1970", "decade 1980", "decade 1990", "decade 2000", "decade 2010", "decade 2020")
  
  labels <- c("Pooled", 
              "Young", "Adult", "Elderly", 
              "Low", "Medium", "High",
              "Low", "Medium", "High", 
              "Low", "Medium", "High", 
              "1970s", "1980s", "1990s", "2000s", "2010s", "2020s")
  
  data_list_df <- rbind(pooled_result_deaths, cont_result_deaths) %>% 
    mutate(key = factor(key, levels = factor_order, labels = labels)) %>% 
    mutate(group_level = factor(group_level, levels = c("Pooled", "Age", "Climate", "Income", "Baseline Mortality Rate", "Decade"))) %>% 
    arrange(key) %>% 
    group_split(group_level) 
  
  plot_list <- lapply(data_list_df, function(df){
    # rev(met.brewer('Homer1'))
    # df <- data_list_df[[1]]
    group_level <- c(unique(df$group_level))
    country <- unique(df$country)
    
    if (group_level == "Age") {
      axis_title_y_element <- element_text(size = 10, family = 'Helvetica')
    } else {
      axis_title_y_element <- element_blank()
    }
    
    # if (country == "US"){
    #   breaks = c(0, 0.02, 0.04, 0.08)
    # } else if (country == "EU") {
    #   breaks = c(0, 0.10, 0.20)
    # } else if (country == "MEX") {
    #   breaks = c(0, 0.04, 0.08, 0.12)
    # }
    
    df %>% 
      ggplot(aes(x = key, y = deaths / total_deaths, fill = side)) +
      geom_bar(stat = "identity", position = "stack", width = 0.3) +
      scale_fill_manual(
        values = c("cold" = "#D3D3D3", "hot" = "#FFCC99"),
        labels = c("Cold" = "cold", "Hot" = "hot")) +
      theme_minimal() +
      scale_y_continuous(
        # breaks = breaks,
        labels = scales::percent_format()
      ) +
      ylab('  \n ') +
      theme(
        strip.text.y = element_blank(),
        axis.text.x = element_text(size = 8, angle = 45, family = "Helvetica", face = "bold"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(color="black"),
        axis.title.y = axis_title_y_element,
        panel.grid = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 5, r = 5, b = 20, l = 10)
      )
  })
  
  # ncol
  ncol <- length(plot_list)
  
  final_plot <- plot_grid(plotlist = plot_list, ncol = ncol, align = 'h')
  
  # Print final plot
  print(final_plot)
  
  return(final_plot)
}

############################################################################
# MAIN FIGIRE 2
############################################################################
plot1 <- plot_grid(
  plot_response(combined_response %>% filter(country == "US") %>% filter(model != "Pooled"), si_folder = NULL),
  plot_annual_deaths(combined_deaths %>% filter(country == "US") %>% filter(model != "Pooled")),
  ncol = 1,
  labels = c("(US)"), 
  label_x = 0, 
  label_y = 1, 
  label_size = 10,
  rel_heights = c(1, 0.5)
)

plot2 <- plot_grid(
  plot_response(combined_response %>% filter(country == "EU") %>% filter(model != "Pooled"), si_folder = NULL),
  plot_annual_deaths(combined_deaths %>% filter(country == "EU") %>% filter(model != "Pooled")),
  ncol = 1,
  labels = c("(EU)"), 
  label_x = 0, 
  label_y = 1, 
  label_size = 10,
  rel_heights = c(1, 0.5)
)

plot3 <- plot_grid(
  plot_response(combined_response %>% filter(country == "MEX") %>% filter(model != "Pooled"), si_folder = NULL),
  plot_annual_deaths(combined_deaths %>% filter(country == "MEX") %>% filter(model != "Pooled")),
  ncol = 1,
  labels = c("(MEX)"), 
  label_x = 0, 
  label_y = 1, 
  label_size = 10,
  rel_heights = c(1, 0.5)
)

fig_2 <- plot_grid(
  plot1, 
  plot2, 
  plot3, 
  ncol = 1, 
  nrow = 3
)

fig_2_dir <- file.path("fig/combined/new_approach", paste0("fig_2.pdf"))
ggsave(plot = fig_2, filename = fig_2_dir, width = 12, height = 12, units = "in", dpi = 300)

############################################################################
# MAIN FIGIRE 2: LOG RATE
############################################################################
# plot_annual_deaths_log_rate <- function(data) {
  # data <- us$annual_deaths
  data <- data %>% 
    mutate(group_level = case_when(
      model_name == 'Pooled' ~ "Pooled",
      model_name == 'temp' ~ "Climate",
      model_name == 'income' ~ "Income",
      model_name == 'baseline' ~ "Baseline Mortality Rate",
      model_name == 'Young' ~ "Age", 
      model_name == 'Adult' ~ "Age", 
      model_name == 'Elderly' ~ "Age", 
      model_name == 'decade' ~ "Decade",
      .default = as.character(model_level)
    ))
  
  pooled_result_deaths <- data %>% 
    filter(group_level == "Pooled" | group_level == "Age" | group_level == "Decade") %>% 
    mutate(model_level = paste(model_name, model_level)) %>% 
    select(year, side, model_level, group_level, country, starts_with("deaths"))
  
  cont_result_deaths <- data %>% 
    filter(group_level != "Age" & group_level != "Decade" & group_level != "Pooled") %>% 
    group_by(model_name) %>% 
    mutate(percentile = as.numeric(model_level)) %>% 
    mutate(percentile = ntile(percentile, n = 3)) %>% 
    mutate(model_level = case_when(
      is.na(percentile) ~ model_level,  
      percentile == 1 ~ "Low", 
      percentile == 2 ~ "Medium", 
      percentile == 3 ~ "High"
    )) %>% 
    ungroup() %>% 
    mutate(model_level = paste(model_level, model_name)) %>% 
    group_by(year, side, model_level) %>% 
    summarise(group_level = first(group_level),
              country = first(country),
              across(starts_with("deaths"), ~mean(.x, na.rm = TRUE)), 
              .groups = "drop")
  
  combined_data <- rbind(pooled_result_deaths, cont_result_deaths) %>% 
    group_by(side, model_level) %>% 
    summarise(
      deaths = mean(deaths, na.rm = T),    
      group_level = first(group_level), 
      country = first(country), 
      .groups = "drop"
    ) %>% 
    arrange(group_level)
  
  factor_order <- c("Pooled pooled", 
                    "Young pooled", "Adult pooled", "Elderly pooled", 
                    "Low temp", "Medium temp", "High temp",
                    "Low income", "Medium income", "High income", 
                    "Low baseline", "Medium baseline", "High baseline", 
                    "decade 1970", "decade 1980", "decade 1990", "decade 2000", "decade 2010")
  
  labels <- c("Pooled", 
              "Young", "Adult", "Elderly", 
              "Low", "Medium", "High",
              "Low", "Medium", "High", 
              "Low", "Medium", "High", 
              "1970s", "1980s", "1990s", "2000s", "2010s")
  
  data_list_df <- combined_data %>% 
    mutate(model_level = factor(model_level, levels = factor_order, labels = labels)) %>% 
    mutate(group_level = factor(group_level, levels = c("Pooled", "Age", "Climate", "Income", "Baseline Mortality Rate", "Decade"))) %>% 
    arrange(model_level) %>% 
    group_split(group_level) 
  
  plot_list <- lapply(data_list_df, function(df){
    # rev(met.brewer('Homer1'))
    # df <- data_list_df[[1]]
    group_level <- c(unique(df$group_level))
    country <- unique(df$country)
    
    if (group_level == "Age") {
      axis_title_y_element <- element_text(size = 10, family = 'Helvetica')
    } else {
      axis_title_y_element <- element_blank()
    }
    
    # if (country == "US"){
    #   breaks = c(0, 0.02, 0.04, 0.08)
    # } else if (country == "EU") {
    #   breaks = c(0, 0.10, 0.20)
    # } else if (country == "MEX") {
    #   breaks = c(0, 0.04, 0.08, 0.12)
    # }
    
    df %>% 
      ggplot(aes(x = model_level, y = deaths, fill = side)) +
      geom_bar(stat = "identity", position = "stack", width = 0.3) +
      scale_fill_manual(
        values = c("cold" = "#D3D3D3", "hot" = "#FFCC99"),
        labels = c("Cold" = "cold", "Hot" = "hot")) +
      theme_minimal() +
      scale_y_continuous(
        # breaks = breaks,
        labels = scales::percent_format()
      ) +
      ylab('  \n ') +
      theme(
        strip.text.y = element_blank(),
        axis.text.x = element_text(size = 8, family = "Helvetica", face = "bold"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(color="black"),
        axis.title.y = axis_title_y_element,
        panel.grid = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 5, r = 5, b = 20, l = 10)
      )
  })
  
  # ncol
  ncol <- length(plot_list)
  
  final_plot <- plot_grid(plotlist = plot_list, ncol = ncol, align = 'h')
  
  # Print final plot
  print(final_plot)
  
  return(final_plot)
}
us_log_rate <- read_country_results(country_name = "US", si_folder = "log_rate")
eu_log_rate <- read_country_results(country_name = "EU", si_folder = "log_rate")
mex_log_rate <- read_country_results(country_name = "MEX", si_folder = "log_rate")

combined_response_log_rate <- rbind(us_log_rate$response, eu_log_rate$response, mex_log_rate$response)
combined_death_log_rate <- rbind(us_log_rate$deaths, eu_log_rate$deaths, mex_log_rate$deaths)

plot_pooled_fig(combined_response_log_rate, combined_death_log_rate, si_folder = "log_rate")

log_rate_only_response <- plot_grid(
  plot_response(combined_response_log_rate %>% filter(country == "US") %>% filter(model != "Pooled"), si_folder = "log_rate"), 
  plot_response(combined_response_log_rate %>% filter(country == "EU") %>% filter(model != "Pooled"), si_folder = "log_rate"),
  plot_response(combined_response_log_rate %>% filter(country == "MEX") %>% filter(model != "Pooled"), si_folder = "log_rate"), 
  ncol = 1, 
  labels = c("(US)", "(EU)", "(MEX)")
)

ggsave(plot = log_rate_only_response, 
       filename = file.path("fig/combined/new_approach", paste0("fig_2_log_rate.pdf")), 
       width = 12, 
       height = 12, 
       units = "in", 
       dpi = 300)

# plot1 <- plot_grid(
#   plot_response(combined_response_log_rate %>% filter(country == "US") %>% filter(model != "Pooled"), si_folder = "log_rate"),
#   plot_annual_deaths(combined_death_log_rate %>% filter(country == "US") %>% filter(model != "Pooled")),
#   ncol = 1,
#   labels = c("(US)"), 
#   label_x = 0, 
#   label_y = 1, 
#   label_size = 10,
#   rel_heights = c(1, 0.5)
# )
# 
# plot2 <- plot_grid(
#   plot_response(combined_response_log_rate %>% filter(country == "EU") %>% filter(model != "Pooled"), si_folder = "log_rate"),
#   plot_annual_deaths(combined_death_log_rate %>% filter(country == "EU") %>% filter(model != "Pooled")),
#   ncol = 1,
#   labels = c("(EU)"), 
#   label_x = 0, 
#   label_y = 1, 
#   label_size = 10,
#   rel_heights = c(1, 0.5)
# )
# 
# plot3 <- plot_grid(
#   plot_response(combined_response_log_rate %>% filter(country == "MEX") %>% filter(model != "Pooled"), si_folder = "log_rate"),
#   plot_annual_deaths(combined_death_log_rate %>% filter(country == "MEX") %>% filter(model != "Pooled")),
#   ncol = 1,
#   labels = c("(MEX)"), 
#   label_x = 0, 
#   label_y = 1, 
#   label_size = 10,
#   rel_heights = c(1, 0.5)
# )
# 
# fig_2 <- plot_grid(
#   plot1, 
#   plot2, 
#   plot3, 
#   ncol = 1, 
#   nrow = 3
# )
# 
# fig_2_dir <- file.path("fig/combined/new_approach", paste0("fig_2_log_rate.pdf"))
# ggsave(plot = fig_2, filename = fig_2_dir, width = 12, height = 12, units = "in", dpi = 300)

############################################################################
# MAIN FIGIRE 2: AGE STANDARDIZED RATE (US and MEX)
############################################################################
us_age_standardized <- read_country_results(country_name = "US", si_folder = "age_standardized")
mex_age_standardized <- read_country_results(country_name = "MEX", si_folder = "age_standardized")

combined_response_age_standardized <- rbind(us_age_standardized$response, mex_age_standardized$response)
combined_death_age_standardized <- rbind(us_age_standardized$deaths, mex_age_standardized$deaths)

plot1 <- plot_grid(
  plot_response(combined_response_age_standardized %>% filter(country == "US") %>% filter(model != "Pooled" & model != "Adult" & model != "Elderly" & model != "Young"), si_folder = "age_standardized"),
  plot_annual_deaths(combined_death_age_standardized %>% filter(country == "US") %>%filter(model != "Pooled" & model != "Adult" & model != "Elderly" & model != "Young")),
  ncol = 1,
  labels = c("(US)"), 
  label_x = 0, 
  label_y = 1, 
  label_size = 10,
  rel_heights = c(1, 0.5)
)

plot2 <- plot_grid(
  plot_response(combined_response_age_standardized %>% filter(country == "MEX") %>% filter(model != "Pooled" & model != "Adult" & model != "Elderly" & model != "Young"), si_folder = "age_standardized"),
  plot_annual_deaths(combined_death_age_standardized %>% filter(country == "MEX") %>%filter(model != "Pooled" & model != "Adult" & model != "Elderly" & model != "Young")),
  ncol = 1,
  labels = c("(MEX)"), 
  label_x = 0, 
  label_y = 1, 
  label_size = 10,
  rel_heights = c(1, 0.5)
)

age_std_rate_plot <- plot_grid(
  plot1, 
  plot2,
  ncol = 1, 
  nrow = 2
)

ggsave(plot = age_std_rate_plot, 
       filename = file.path("fig/combined/new_approach", paste0("fig_2_age_std_us_mex.pdf")),
       width = 12, 
       height = 12, 
       units = "in", 
       dpi = 300)

###################################################################################
#   Fig 1: Supplementary
###################################################################################
us_same_year <- read_country_results(country_name = "US", si_folder = "same_year")
eu_same_year <- read_country_results(country_name = "EU", si_folder = "same_year")
mex_same_year <- read_country_results(country_name = "MEX", si_folder = "same_year")

combined_response_same_year <- rbind(us_same_year$response, eu_same_year$response, mex_same_year$response)
combined_death_same_year <- rbind(us_same_year$annual_deaths, eu_same_year$annual_deaths, mex_same_year$annual_deaths)

fig_1_same_year <- plot_pooled_fig(combined_response_same_year, combined_death_same_year, si_folder = "same_year")
fig_1_dir_same_year <- file.path("fig/combined", paste0("fig_1_same_year.pdf"))
ggsave(plot = fig_1_same_year, filename = fig_1_dir_same_year, width = 10, height = 5, units = "in", dpi = 300)


###################################################################################
# Fig 3: Percentage Attributable Deaths by Per Bins 
###################################################################################

plot_deaths_at_bins <- function(){
  country_names <- list("US", "EU", "MEX")
  
  df_list <- lapply(country_names, function(country_name){
    # country_name <- "US"
    processed_country_dir <- paste("processed", country_name, sep = "/")
    processed_files <- list.files(processed_country_dir)
    group_denominator <- T
    
    # bring in deaths by bins
    death_by_bins_files <- processed_files[grepl("deaths_by_bins", processed_files) & 
                                             (grepl("young", processed_files) |
                                                grepl("adult", processed_files) |
                                                grepl("elderly", processed_files) | 
                                                grepl("pooled", processed_files))]
    
    death_by_bins <- lapply(death_by_bins_files, function(file){readRDS(file.path(processed_country_dir, file))})
    death_by_bins <- do.call(rbind, death_by_bins)
    
    # total deaths
    total_deaths_files <- processed_files[(grepl("deaths_by_side_year", processed_files) & !grepl(pattern = "iteration", x = processed_files)) & 
                                            (grepl("young", processed_files) |
                                               grepl("adult", processed_files) |
                                               grepl("elderly", processed_files) | 
                                               grepl("pooled", processed_files))]
    total_deaths <- lapply(total_deaths_files, function(file){readRDS(file.path(processed_country_dir, file))})
    total_deaths <- do.call(rbind, total_deaths)
    total_deaths <- distinct(total_deaths[, .(year, annual_total_death, model_name)])
    
    df <- merge(death_by_bins, total_deaths, by = c("year", "model_name"), all.x = T)
    
    if (group_denominator) {
      df <- df %>% 
        filter(deaths != 0 & deaths > 0) %>%
        mutate(deaths = deaths / annual_total_death) %>% 
        group_by(bins, model_name) %>% 
        summarise(deaths = mean(deaths, na.rm = TRUE), .groups = "drop") %>% 
        ungroup()
    } else { # denominator is overall total deaths
      df <- df %>% group_by(year, bins) %>% mutate(annual_total_death = sum(annual_total_death, na.rm = T)) %>% ungroup() %>% 
        filter(deaths != 0 & deaths > 0) %>%
        mutate(deaths = deaths / annual_total_death) %>% 
        group_by(bins, model_name) %>% 
        summarise(deaths = mean(deaths, na.rm = TRUE), .groups = "drop") %>% 
        ungroup()
    }
    
    df$country_name <- country_name
    return(df)
  })
  df <- do.call(rbind, df_list)
  
  df <- df %>%
    group_by(model_name, country_name) %>%
    mutate(deaths = ifelse(bins < -20, sum(deaths[bins < -20]), deaths)) %>%
    filter(bins > -21) %>% 
    ungroup()
  
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
  
  color_palette <- c('#0072B2', '#CC79A7','#009E73', '#F0E442')
  # color_palette <- rev(met.brewer("Tam", 3))
  
  plot <- df %>% 
    mutate(model_name = factor(model_name, levels = c("Pooled","Elderly", "Adult", "Young"))) %>% 
    arrange(model_name, country_name) %>% 
    mutate(country_name = factor(country_name, levels = c("US", "EU", "MEX"))) %>%
    ggplot(aes(x = bins, y = deaths, fill = model_name)) + 
    geom_bar(stat = 'identity', position = 'stack', width = 1, alpha = 1) + 
    scale_y_continuous(
      labels = scales::percent_format()
    ) + 
    scale_fill_manual(values = color_palette, name = '') + 
    theme_minimal() +
    xlab("Temperature (°C)") + 
    ylab('Percent of total annual deaths') +
    facet_grid(rows = vars(model_name), cols = vars(country_name), scales = "free_y") +
    theme(
      strip.text.x = element_text(size = 10, face = "bold"),
      strip.text.y = element_blank(),
      panel.grid.major.x = element_line(color = "gray80", linetype = "dotted", linewidth = 0.5), 
      panel.grid.minor.x = element_blank(), 
      panel.grid.major.y = element_line(color = "gray90", linetype = "solid", linewidth = 0.5), 
      panel.grid.minor.y = element_blank(),
      panel.spacing = unit(1.5, "lines")
    )

  return(plot)
}
fig_3 <- plot_deaths_at_bins()
fig_3_dir <- file.path("fig/combined", paste0("fig_3.pdf"))
ggsave(plot = fig_3, filename = fig_3_dir, width = 10, height = 6, units = "in", dpi = 300)

###################################################################################
# Fig 3: Supplementary decade fixed 
###################################################################################
country_name = "US"
processed_country_dir <- paste("processed", country_name, sep = "/")
model_name <- "decade"
result_decade <- readRDS(file.path(processed_country_dir, paste0( country_name, "_", "decade", "_results.rds" )))
result_decade_fixed1970 <- readRDS(file.path(processed_country_dir, paste0( country_name, "_", "decade_fixed1970", "_results.rds" )))
deaths_decade <- result_decade$deaths
deaths_decade[, model := "Observed" ]
deaths_decade_fixed1970 <- result_decade_fixed1970$deaths
deaths_decade_fixed1970[, model := "Counterfactual"]

theme_clean <- function() {
  theme_classic() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.title.x = element_text(margin = margin(t = 10)),
      plot.title = element_text(color = "black", size = 12, hjust = 0.5, face = "bold"),
      axis.line = element_line(color = 'black', linewidth = 0.35),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black")
    )
}

decade_fixed1970 <- rbind(deaths_decade, deaths_decade_fixed1970) %>% 
  mutate( model = factor(model, levels = c("Observed", "Counterfactual"))) %>% 
  ggplot(aes(x = model, y = deaths / total_deaths, fill = side)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  facet_wrap(~ key, nrow = 1) +
  scale_fill_manual(
    values = c("cold" = "#D3D3D3", "hot" = "#FFCC99")) +
  theme_minimal() +
  scale_y_continuous(
    labels = scales::percent_format()
  ) +
  ylab('  \n ') +
  theme(
    strip.text.x = element_text(size = 10, family = "Helvetica"),
    axis.text.x = element_text(size = 8, family = "Helvetica", angle = 45),
    axis.title.x = element_blank(),
    axis.ticks.x = element_line(color="black"),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 20, l = 10)
  )

ggsave(
  plot = decade_fixed1970, 
  filename = file.path("fig/combined/new_approach", paste0("decade_fixed1970.pdf")), 
  width = 10, 
  height = 8, 
  units = "in", 
  dpi = 300
)

# ############################################################################
# #   Attributable deaths by each percentile
# ############################################################################
# plot_temp_income <- function(country_name, df) {
#   # country_name = "US"
#   # df <- combined_death %>% filter(country == "US")
#   # Create output directory if it doesn't exist
#   fig_country_dir <- file.path("fig", country_name)
#   if (!dir.exists(fig_country_dir)) dir.create(fig_country_dir, recursive = TRUE)
#   
#   # Define custom theme
#   theme_clean <- function() {
#     theme_classic() +
#       theme(
#         axis.text = element_text(size = 10),
#         axis.title = element_text(size = 10),
#         axis.title.y = element_text(margin = margin(r = 10)),
#         axis.title.x = element_text(margin = margin(t = 10)),
#         plot.title = element_text(color = "black", size = 12, hjust = 0.5, face = "bold"),
#         axis.line = element_line(color = 'black', linewidth = 0.35),
#         axis.ticks = element_line(color = "black", linewidth = 0.35),
#         axis.text.x = element_text(color = "black"),
#         axis.text.y = element_text(color = "black")
#       )
#   }
#   
#   # Define color palette
#   okabe_colors <- c('#D55E00','#0072B2','#009E73','#CC79A7','#000000',
#                     '#E69F00','#F0E442','#4DA3E5','#555555','#D9D9D9')
#   
#   # Process and plot temperature data
#   df_temp <- df %>%
#     filter(model_name == "temp") %>%
#     group_by(model_level, side) %>%
#     summarize(deaths = median(deaths / annual_total_death, na.rm = TRUE), .groups = "drop") %>%
#     mutate(model_level = as.numeric(model_level) + 0.5,
#            side = factor(side, levels = c('hot', 'cold'), labels = c('Heat', 'Cold')))
#   
#   plot_temp <- ggplot(df_temp, aes(x = model_level, y = deaths, fill = side)) +
#     geom_bar(stat = 'identity', position = 'stack') +
#     scale_y_continuous(labels = scales::percent_format()) +
#     scale_fill_manual(values = okabe_colors, name = 'Death due to ...') +
#     theme_clean() +
#     xlab('Percentile of annual mean\ntemperature distribution') +
#     ylab('Percent of total annual deaths')
#   
#   print(plot_temp)
#   
#   temp_out_path <- file.path(fig_country_dir, paste0(country_name, "_temp_percentile_bar_plot.png"))
#   ggsave(plot = plot_temp, filename = temp_out_path, width = 20, height = 8, units = "in", dpi = 300)
#   
#   # Process and plot income data
#   df_income <- df %>%
#     filter(model_name == "income") %>%
#     filter(annual_total_death != 0) %>% 
#     group_by(model_level, side) %>%
#     summarize(deaths = median(deaths / annual_total_death, na.rm = TRUE), .groups = "drop") %>%
#     mutate(model_level = as.numeric(model_level) + 0.5,
#            side = factor(side, levels = c('hot', 'cold'), labels = c('Heat', 'Cold'))) %>% arrange(model_level)
#   
#   plot_income <- ggplot(df_income, aes(x = model_level, y = deaths, fill = side)) +
#     geom_bar(stat = 'identity', position = 'stack') +
#     scale_y_continuous(labels = scales::percent_format()) +
#     scale_fill_manual(values = okabe_colors, name = 'Death due to ...') +
#     theme_clean() +
#     xlab('Percentile of annual mean\nincome distribution') +  # Fixed label for income
#     ylab('Percent of total annual deaths')
#   
#   print(plot_income)
#   
#   income_out_path <- file.path(fig_country_dir, paste0(country_name, "_income_percentile_bar_plot.png"))
#   ggsave(plot = plot_income, filename = income_out_path, width = 20, height = 8, units = "in", dpi = 300)
# }
# 
# combined_response %>% filter(model_level == "97" & model_name == "income" & country == "US") %>% 
#   filter(bins > -20 & bins<40) %>% 
#   ggplot(aes(x = bins)) +
#   geom_line(aes(y = pred, color = model_level), linewidth = 1) +
#   geom_ribbon(aes(ymin = pred_q025, ymax = pred_q975, fill = model_level), alpha = 0.2)
# 
# plot_temp_income("US", us$annual_deaths)
# plot_temp_income("EU", eu$annual_deaths)
# plot_temp_income("MEX", mex$annual_deaths)