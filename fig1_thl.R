## fig1_thl.R — figure 1 (pooled response + attributable deaths) for THL_1997_2025
suppressPackageStartupMessages(
  pacman::p_load(fastverse, tidyverse, MetBrewer, cowplot, patchwork, here)
)
source(here::here("config", "countries.R"))

# ── Load THL_1997_2025 Pooled results ─────────────────────────────────────────
rds_path <- here::here("processed", "THL", "THL_1997_2025", "THL_Pooled_results.rds")
if (!file.exists(rds_path)) stop("Results not found: ", rds_path)
thl_res <- readRDS(rds_path)

fig1_response <- as.data.table(thl_res$response)[, country := "THL_1997_2025"]
fig1_deaths   <- as.data.table(thl_res$deaths)[,   country := "THL_1997_2025"]

# ── plot_pooled_fig ───────────────────────────────────────────────────────────
# Adapted from attribution_plot.R; factor levels and palette are dynamic.
plot_pooled_fig <- function(pooled_dataframe, deaths_dataframe, si_folder = NULL) {

  pooled_response_df <- pooled_dataframe %>% filter(model == "Pooled")
  pooled_deaths_df   <- deaths_dataframe %>% filter(model == "Pooled")

  country_order  <- c("US","EU","MEX","BRA","LA","THL_1997_2025","THL_2005_2025")
  present_levels <- intersect(country_order, unique(pooled_response_df$country))
  pooled_response_df <- pooled_response_df %>%
    mutate(country = factor(country, levels = present_levels))

  if (!is.null(si_folder) && si_folder == "log_rate") {
    pooled_response_df <- pooled_response_df %>%
      group_by(country) %>%
      mutate(pred_q025 = expm1(pred_q025), pred = expm1(pred), pred_q975 = expm1(pred_q975)) %>%
      ungroup() %>%
      mutate(country = factor(country, levels = present_levels))
    y_lab_text <- "Percent change in deaths"
  } else {
    pooled_response_df <- pooled_response_df %>%
      group_by(country) %>%
      mutate(pred_q025 = pred_q025 * 1e5, pred = pred * 1e5, pred_q975 = pred_q975 * 1e5) %>%
      ungroup() %>%
      mutate(country = factor(country, levels = present_levels))
    y_lab_text <- "Change in deaths \n per 100,000 people"
  }

  unique_levels <- levels(pooled_response_df$country)
  color_palette <- c("US" = "#009988", "EU" = "#EE7733", "MEX" = "#0072B2",
                     "BRA" = "#CC79A7", "LA" = "#AA4499",
                     "THL_1997_2025" = "#E69F00", "THL_2005_2025" = "#56B4E9")

  theme_clean <- function() {
    theme_classic() +
      theme(
        axis.text        = element_text(size = 10),
        axis.title       = element_text(size = 10),
        axis.title.y     = element_text(margin = margin(r = 10)),
        axis.title.x     = element_text(margin = margin(t = 10)),
        plot.title       = element_text(color = "black", size = 12, hjust = 0.5, face = "bold"),
        axis.line        = element_line(color = "black", linewidth = 0.35),
        axis.ticks       = element_line(colour = "black", linewidth = 0.35),
        axis.text.x      = element_text(color = "black"),
        axis.text.y      = element_text(color = "black"),
        strip.text       = element_text(size = 10, face = "bold")
      )
  }

  response <- pooled_response_df %>%
    filter(pop_deg_days != 0 & !is.na(pop_deg_days), bins >= t1, bins <= t99) %>%
    ggplot(aes(x = bins)) +
    geom_ribbon(aes(ymin = pred_q025, ymax = pred_q975, fill = country), alpha = 0.2) +
    geom_line(aes(y = pred, color = country), linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Temperature (°C)", y = y_lab_text) +
    scale_color_manual(values = color_palette) +
    scale_fill_manual(values = color_palette) +
    theme_clean() +
    theme(
      axis.text.x  = element_blank(), axis.title.x = element_blank(),
      axis.line.x  = element_blank(), axis.ticks.x = element_blank(),
      panel.grid.major.x = element_line(color = "gray80", linetype = "dotted", linewidth = 0.5),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
      legend.background  = element_rect(fill = "transparent", colour = NA),
      legend.position    = "inside", legend.position.inside = c(0.2, 0.7),
      legend.title       = element_blank()
    )

  exposure <- pooled_response_df %>%
    filter(pop_deg_days != 0 & !is.na(pop_deg_days), bins >= t1, bins <= t99) %>%
    ggplot(aes(x = bins, weight = pop_deg_days, fill = country, color = country)) +
    geom_density(alpha = 0.3, bw = 1) +
    scale_color_manual(values = setNames(color_palette, unique_levels)) +
    scale_fill_manual(values  = setNames(color_palette, unique_levels)) +
    labs(x = "Temperature (°C)", y = " \n ") +
    theme_minimal() +
    theme(legend.position = "none", axis.title = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid = element_blank(), plot.margin = margin(0, 5, 0, 20))

  final_plot <- plot_grid(response, exposure, ncol = 1, rel_heights = c(2, 0.5), align = "vh")

  breaks <- c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14)
  deaths <- pooled_deaths_df %>%
    mutate(annual_share_total_mortality = deaths / total_deaths,
           country = factor(country, levels = present_levels)) %>%
    ggplot(aes(x = country, y = annual_share_total_mortality, fill = side)) +
    geom_bar(stat = "identity", position = "stack", width = 0.9, alpha = 1) +
    scale_fill_manual(values = c("cold" = "#D3D3D3", "hot" = "#FFCC99"),
                      labels = c("cold" = "Cold", "hot" = "Hot"),
                      name   = "Mortality due to ...") +
    scale_y_continuous(breaks = breaks, labels = scales::percent_format()) +
    ylab("Percent of total annual deaths") +
    theme_clean() +
    theme(axis.text.x = element_text(size = 10, angle = 45), axis.title.x = element_blank())

  plot_grid(final_plot, deaths, ncol = 2)
}

# ── Save ──────────────────────────────────────────────────────────────────────
fig_1 <- plot_pooled_fig(fig1_response, fig1_deaths, si_folder = NULL)

out_dir <- here::here("fig", "THL", "THL_1997_2025")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_path <- file.path(out_dir, "fig_1_thl_1997_2025.pdf")
ggsave(plot = fig_1, filename = out_path, width = 10, height = 5, units = "in", dpi = 300)
message("Saved: ", out_path)
