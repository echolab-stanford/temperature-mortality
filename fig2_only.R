## fig2_only.R  — generates fig/combined/fig_2.pdf (heterogeneity panels only)
suppressPackageStartupMessages(
  pacman::p_load(fastverse, tidyverse, MetBrewer, cowplot, patchwork, here)
)
source(here::here("config", "countries.R"))

# ── data loading ──────────────────────────────────────────────────────────────
read_country_results <- function(country_name, si_folder = NULL) {
  if (!country_name %in% names(country_config))
    stop("'", country_name, "' not in config/countries.R")
  dir <- here::here("processed", country_name)
  if (!is.null(si_folder)) dir <- file.path(dir, si_folder)
  models <- country_config[[country_name]]$models_to_run
  res <- Filter(Negate(is.null), lapply(setNames(models, models), function(m) {
    p <- file.path(dir, paste0(country_name, "_", m, "_results.rds"))
    if (!file.exists(p)) { warning("Missing: ", p); return(NULL) }
    readRDS(p)
  }))
  resp   <- rbindlist(lapply(res, `[[`, "response"), fill = TRUE)
  deaths <- rbindlist(lapply(res, `[[`, "deaths"),   fill = TRUE)
  resp[,   country := country_name]
  deaths[, country := country_name]
  list(response = resp, deaths = deaths)
}

fig123_countries <- names(Filter(function(x) all(c(1,2,3) %in% x$figs), country_config))
message("Loading: ", paste(fig123_countries, collapse = ", "))
all_results       <- lapply(setNames(fig123_countries, fig123_countries), read_country_results)
combined_response <- rbindlist(lapply(all_results, `[[`, "response"), fill = TRUE)
combined_deaths   <- rbindlist(lapply(all_results, `[[`, "deaths"),   fill = TRUE)

# ── plot_response ─────────────────────────────────────────────────────────────
plot_response <- function(response_df, si_folder = NULL) {
  response_df <- response_df %>%
    filter(!is.na(pop_deg_days)) %>%
    filter(key %in% c("pooled","1970","1980","1990","2000","2010","2020",
                      "10%","50%","90%")) %>%
    mutate(group_level = case_when(
      model == "Pooled"   ~ "Pooled",
      model == "temp"     ~ "Climate",
      model == "income"   ~ "Income",
      model == "baseline" ~ "Baseline Mortality Rate",
      model %in% c("Young","Adult","Elderly") ~ "Age",
      model == "decade"   ~ "Decade",
      .default = as.character(model)
    )) %>%
    mutate(key = paste(model, key))

  country_name <- unique(response_df$country)
  response_df <- response_df %>%
    mutate(across(c(pred_q025, pred, pred_q975), ~ .x * 1e5))
  y_lab_text <- "Change in deaths\nper 100,000 people"

  response_df_list <- response_df %>%
    mutate(key = case_when(
      key == "temp 10%"     ~ "Temp 10th",   key == "temp 50%"     ~ "Temp 50th",
      key == "temp 90%"     ~ "Temp 90th",   key == "income 10%"   ~ "Income 10th",
      key == "income 50%"   ~ "Income 50th", key == "income 90%"   ~ "Income 90th",
      key == "baseline 10%" ~ "Baseline 10th", key == "baseline 50%" ~ "Baseline 50th",
      key == "baseline 90%" ~ "Baseline 90th",
      key == "Adult pooled"   ~ "Adult",   key == "Elderly pooled" ~ "Elderly",
      key == "Young pooled"   ~ "Young",
      key == "decade 1970" ~ "1970s", key == "decade 1980" ~ "1980s",
      key == "decade 1990" ~ "1990s", key == "decade 2000" ~ "2000s",
      key == "decade 2010" ~ "2010s", key == "decade 2020" ~ "2020s",
      .default = as.character(key)
    )) %>%
    mutate(group_level = factor(group_level,
      levels = c("Pooled","Age","Climate","Income","Baseline Mortality Rate","Decade"))) %>%
    group_split(group_level)

  plot_list <- lapply(response_df_list, function(df) {
    gl      <- unique(df$group_level)
    country <- unique(df$country)
    if (gl == "Age")
      df <- df %>% mutate(key = factor(key, levels = c("Young","Adult","Elderly"))) %>% arrange(bins, key)
    else if (gl == "Decade")
      df <- df %>% mutate(key = factor(key, levels = c("1970s","1980s","1990s","2000s","2010s","2020s"))) %>% arrange(key)

    strip_el <- if (country == "US") element_text(size = 10, face = "bold") else element_blank()
    ul <- unique(df$key); n <- length(ul)
    pal <- rev(c("#a62f00","#df7700","#32b2da","#6ad5e8","#fff179")[seq_len(n)])

    df <- df %>% filter(pop_deg_days != 0, !is.na(pop_deg_days), bins >= t1, bins <= t99)

    resp <- df %>%
      ggplot(aes(x = bins)) +
      geom_ribbon(aes(ymin = pred_q025, ymax = pred_q975, fill = key), alpha = 0.2) +
      geom_line(aes(y = pred, color = key), linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_color_manual(values = setNames(pal, ul)) +
      scale_fill_manual(values  = setNames(pal, ul)) +
      labs(y = y_lab_text) +
      facet_wrap(~group_level, scales = "free_y") +
      coord_cartesian(ylim = c(0, NA)) +
      theme_minimal() +
      theme(
        strip.text      = strip_el,
        axis.title.y    = element_blank(),
        axis.title.x    = element_blank(),
        axis.text.x     = element_blank(),
        axis.line.x     = element_blank(),
        axis.ticks.x    = element_blank(),
        panel.grid      = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.7, 0.8),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.title    = element_blank(),
        legend.text     = element_text(size = 8),
        legend.key.size = unit(0.5, "lines"),
        legend.key.spacing.y = unit(0.1, "lines"),
        plot.margin     = margin(0, 5, 0, 10)
      )

    expo <- df %>%
      ggplot(aes(x = bins, weight = pop_deg_days, fill = key, color = key)) +
      geom_density(alpha = 0.3, bw = 1) +
      scale_color_manual(values = setNames(pal, ul)) +
      scale_fill_manual(values  = setNames(pal, ul)) +
      theme_minimal() +
      theme(legend.position = "none", axis.title = element_blank(),
            axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            panel.grid = element_blank(), plot.margin = margin(0, 5, 0, 20))

    plot_grid(resp, expo, ncol = 1, rel_heights = c(1, 0.2))
  })

  plot_grid(plotlist = plot_list, ncol = length(plot_list), align = "hv")
}

# ── plot_annual_deaths ────────────────────────────────────────────────────────
plot_annual_deaths <- function(data) {
  data <- data %>%
    mutate(group_level = case_when(
      model == "Pooled"   ~ "Pooled",
      model == "temp"     ~ "Climate",
      model == "income"   ~ "Income",
      model == "baseline" ~ "Baseline Mortality Rate",
      model %in% c("Young","Adult","Elderly") ~ "Age",
      model == "decade"   ~ "Decade",
      .default = as.character(model)
    ))

  pooled_d <- data %>%
    filter(group_level %in% c("Pooled","Age","Decade")) %>%
    mutate(key = paste(model, key))

  cont_d <- data %>%
    filter(!group_level %in% c("Age","Decade","Pooled")) %>%
    group_by(model) %>%
    mutate(pct = as.numeric(sub("%","", key)),
           tile = ntile(pct, 3),
           key  = case_when(is.na(tile) ~ key, tile == 1 ~ "Low",
                            tile == 2 ~ "Medium", tile == 3 ~ "High")) %>%
    ungroup() %>%
    mutate(key = paste(key, model)) %>%
    group_by(side, key) %>%
    summarise(group_level = first(group_level), country = first(country),
              model = first(model),
              across(contains("deaths"), ~ sum(.x, na.rm = TRUE)),
              .groups = "drop") %>%
    select(side, total_deaths, deaths_q025, deaths, deaths_q975,
           key, model, country, group_level)

  factor_order <- c("Pooled pooled",
                    "Young pooled","Adult pooled","Elderly pooled",
                    "Low temp","Medium temp","High temp",
                    "Low income","Medium income","High income",
                    "Low baseline","Medium baseline","High baseline",
                    "decade 1970","decade 1980","decade 1990",
                    "decade 2000","decade 2010","decade 2020")
  labels <- c("Pooled","Young","Adult","Elderly",
              "Low","Medium","High","Low","Medium","High",
              "Low","Medium","High",
              "1970s","1980s","1990s","2000s","2010s","2020s")

  data_list <- rbind(pooled_d, cont_d) %>%
    mutate(key = factor(key, levels = factor_order, labels = labels),
           group_level = factor(group_level,
             levels = c("Pooled","Age","Climate","Income",
                        "Baseline Mortality Rate","Decade"))) %>%
    arrange(key) %>%
    group_split(group_level)

  plot_list <- lapply(data_list, function(df) {
    gl      <- unique(df$group_level)
    y_el    <- if (gl == "Age") element_text(size = 10) else element_blank()
    df %>%
      ggplot(aes(x = key, y = deaths / total_deaths, fill = side)) +
      geom_bar(stat = "identity", position = "stack", width = 0.3) +
      scale_fill_manual(values = c("cold" = "#D3D3D3", "hot" = "#FFCC99")) +
      scale_y_continuous(labels = scales::percent_format()) +
      theme_minimal() +
      ylab("  \n ") +
      theme(strip.text.y  = element_blank(),
            axis.text.x   = element_text(size = 8, angle = 45, face = "bold"),
            axis.title.x  = element_blank(),
            axis.ticks.x  = element_line(color = "black"),
            axis.title.y  = y_el,
            panel.grid    = element_blank(),
            legend.position = "none",
            plot.margin   = margin(5, 5, 20, 10))
  })
  plot_grid(plotlist = plot_list, ncol = length(plot_list), align = "h")
}

# ── assemble fig 2 ────────────────────────────────────────────────────────────
make_panel <- function(cname, label) {
  plot_grid(
    plot_response(combined_response %>% filter(country == cname, model != "Pooled")),
    plot_annual_deaths(combined_deaths %>% filter(country == cname, model != "Pooled")),
    ncol = 1, labels = label, label_x = 0, label_y = 1, label_size = 10,
    rel_heights = c(1, 0.5)
  )
}

message("Building panels...")
p_us  <- make_panel("US",  "(US)")
p_eu  <- make_panel("EU",  "(EU)")
p_mex <- make_panel("MEX", "(MEX)")
p_bra <- make_panel("BRA", "(BRA)")

fig_2 <- plot_grid(p_us, p_eu, p_mex, p_bra, ncol = 1, nrow = 4)

out <- here::here("fig", "combined", "fig_2.pdf")
ggsave(plot = fig_2, filename = out, width = 12, height = 20, units = "in", dpi = 300)
message("Saved: ", out)
