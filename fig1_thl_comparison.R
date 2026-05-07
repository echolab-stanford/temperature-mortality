## fig1_thl_comparison.R
##
## Plots pooled temperature-response curves for seven THL specifications
## on the same axes so they can be compared visually:
##
##   1. THL_1997_2025  — crude mortality rate, 913 AMCs, 1997-2025 (winsorized)
##   2. THL_2005_2025  — age-std rate,   920 AMCs, 2005-2025 (full sample)
##   3. THL_2005_no_agg — age-std rate,  725 AMCs, 2005-2025 (no agg-only AMCs)
##   4. THL age-std    — age-std rate,   920 AMCs, 2014-2025
##   5. THL age-std    — age-std rate,   920 AMCs, 2005-2012
##   6. THL crude str-FE — crude rate, 913 AMCs, 1997-2025, str-break FEs
##   7. THL age-std str-FE — age-std rate, 920 AMCs, 2015-2025, str-break FEs
##
## RDS sources:
##   processed/THL/THL_1997_2025/                    THL_Pooled_results.rds
##   processed/THL/THL_2005_2025/age_standardized/   THL_Pooled_results.rds
##   processed/THL/THL_2005_2025/age_std_no_agg/     THL_Pooled_results.rds
##   processed/THL/THL_2005_2025/age_std_2014_2025/  THL_Pooled_results.rds
##   processed/THL/THL_2005_2025/age_std_2005_2012/  THL_Pooled_results.rds
##   processed/THL/THL_1997_2025/crude_str_fe/        THL_Pooled_results.rds
##   processed/THL/THL_2005_2025/age_std_2015_str_fe/ THL_Pooled_results.rds
##
## Output: fig/THL/THL_comparison/fig1_thl_comparison.pdf

suppressPackageStartupMessages(
  pacman::p_load(fastverse, tidyverse, cowplot, here)
)

# ── Helper: load one Pooled_results.rds and label it ─────────────────────────
load_pooled <- function(rds_path, label) {
  if (!file.exists(rds_path)) stop("Results not found: ", rds_path)
  res <- readRDS(rds_path)
  resp   <- as.data.table(res$response)[model == "Pooled"][, country := label]
  deaths <- as.data.table(res$deaths)[model   == "Pooled"][, country := label]
  list(response = resp, deaths = deaths)
}

# ── Load seven specifications ─────────────────────────────────────────────────
specs <- list(
  load_pooled(
    here::here("processed", "THL", "THL_1997_2025", "THL_Pooled_results.rds"),
    "Crude (1997-2025)"
  ),
  load_pooled(
    here::here("processed", "THL", "THL_2005_2025", "age_standardized", "THL_Pooled_results.rds"),
    "Age-std (2005-2025)"
  ),
  load_pooled(
    here::here("processed", "THL", "THL_2005_2025", "age_std_no_agg", "THL_Pooled_results.rds"),
    "Age-std, no-agg (2005-2025)"
  ),
  load_pooled(
    here::here("processed", "THL", "THL_2005_2025", "age_std_2014_2025", "THL_Pooled_results.rds"),
    "Age-std (2014-2025)"
  ),
  load_pooled(
    here::here("processed", "THL", "THL_2005_2025", "age_std_2005_2012", "THL_Pooled_results.rds"),
    "Age-std (2005-2012)"
  ),
  load_pooled(
    here::here("processed", "THL", "THL_1997_2025", "crude_str_fe", "THL_Pooled_results.rds"),
    "Crude str-FE (1997-2025)"
  ),
  load_pooled(
    here::here("processed", "THL", "THL_2005_2025", "age_std_2015_str_fe", "THL_Pooled_results.rds"),
    "Age-std str-FE (2005-2025)"
  )
)

fig1_response <- rbindlist(lapply(specs, `[[`, "response"), fill = TRUE)
fig1_deaths   <- rbindlist(lapply(specs, `[[`, "deaths"),   fill = TRUE)

# ── Interpolate response to 0.1°C resolution ──────────────────────────────────
# attribution.R evaluates the polynomial at integer bins only (1°C steps).
# For Thailand the plotted range is ~21-34°C = only 13-14 points, which
# creates visible straight-line segments in geom_line/geom_ribbon.
# Cubic spline interpolation between the integer-bin predictions restores
# the smooth appearance of the underlying polynomial without re-running the model.
smooth_response <- function(dt, step = 0.1) {
  dt[pop_deg_days > 0 & !is.na(pop_deg_days) & bins >= t1 & bins <= t99][
    , {
      fine <- seq(min(bins), max(bins), by = step)
      list(
        bins         = fine,
        pred         = spline(bins, pred,      xout = fine, method = "natural")$y,
        pred_q025    = spline(bins, pred_q025, xout = fine, method = "natural")$y,
        pred_q975    = spline(bins, pred_q975, xout = fine, method = "natural")$y,
        pop_deg_days = approx(bins, pop_deg_days, xout = fine, rule = 2)$y,
        t1           = t1[1],
        t99          = t99[1],
        model        = model[1]
      )
    },
    by = country
  ]
}

fig1_smooth <- smooth_response(fig1_response)

# ── Color palette & factor order ─────────────────────────────────────────────
label_order  <- c(
  "Crude (1997-2025)", "Age-std (2005-2025)", "Age-std, no-agg (2005-2025)",
  "Age-std (2014-2025)", "Age-std (2005-2012)",
  "Crude str-FE (1997-2025)", "Age-std str-FE (2005-2025)"
)
color_palette <- c(
  "Crude (1997-2025)"            = "#E69F00",
  "Age-std (2005-2025)"          = "#56B4E9",
  "Age-std, no-agg (2005-2025)"  = "#009E73",
  "Age-std (2014-2025)"          = "#CC79A7",
  "Age-std (2005-2012)"          = "#D55E00",
  "Crude str-FE (1997-2025)"     = "#7b1e1ab3",
  "Age-std str-FE (2005-2025)"   = "#03db20"
)

fig1_smooth[,  country := factor(country, levels = label_order)]
fig1_response[, country := factor(country, levels = label_order)]
fig1_deaths[,   country := factor(country, levels = label_order)]

# ── Theme ─────────────────────────────────────────────────────────────────────
theme_clean <- function() {
  theme_classic() +
    theme(
      axis.text       = element_text(size = 10),
      axis.title      = element_text(size = 10),
      axis.title.y    = element_text(margin = margin(r = 10)),
      axis.title.x    = element_text(margin = margin(t = 10)),
      axis.line       = element_line(color = "black", linewidth = 0.35),
      axis.ticks      = element_line(colour = "black", linewidth = 0.35),
      axis.text.x     = element_text(color = "black"),
      axis.text.y     = element_text(color = "black"),
      strip.text      = element_text(size = 10, face = "bold")
    )
}

# ── Response curve panel (uses interpolated data) ─────────────────────────────
response_plot <- fig1_smooth %>%
  mutate(
    pred_q025 = pred_q025 * 1e5,
    pred      = pred      * 1e5,
    pred_q975 = pred_q975 * 1e5
  ) %>%
  ggplot(aes(x = bins)) +
  geom_ribbon(aes(ymin = pred_q025, ymax = pred_q975, fill = country), alpha = 0.2) +
  geom_line(aes(y = pred, color = country), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values  = color_palette) +
  labs(x = "Temperature (°C)", y = "Change in deaths\nper 100,000 people") +
  theme_clean() +
  theme(
    axis.text.x  = element_blank(), axis.title.x = element_blank(),
    axis.line.x  = element_blank(), axis.ticks.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
    legend.background  = element_rect(fill = "transparent", colour = NA),
    legend.position    = "inside", legend.position.inside = c(0.25, 0.78),
    legend.title       = element_blank(),
    legend.text        = element_text(size = 9)
  )

# ── Exposure density panel (uses original 1°C bins) ──────────────────────────
exposure_plot <- fig1_response %>%
  filter(pop_deg_days > 0 & !is.na(pop_deg_days), bins >= t1, bins <= t99) %>%
  ggplot(aes(x = bins, weight = pop_deg_days, fill = country, color = country)) +
  geom_density(alpha = 0.3, bw = 1) +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values  = color_palette) +
  labs(x = "Temperature (°C)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title      = element_blank(),
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    panel.grid      = element_blank(),
    plot.margin     = margin(0, 5, 0, 20)
  )

final_response <- plot_grid(response_plot, exposure_plot, ncol = 1,
                            rel_heights = c(2, 0.5), align = "vh")

# ── Attributable deaths bar panel ─────────────────────────────────────────────
breaks <- c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14)

deaths_plot <- fig1_deaths %>%
  mutate(
    annual_share = deaths / total_deaths,
    country = factor(country, levels = label_order)
  ) %>%
  ggplot(aes(x = country, y = annual_share, fill = side)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6, alpha = 1) +
  scale_fill_manual(
    values = c("cold" = "#D3D3D3", "hot" = "#FFCC99"),
    labels = c("cold" = "Cold", "hot" = "Hot"),
    name   = "Mortality due to ..."
  ) +
  scale_y_continuous(breaks = breaks, labels = scales::percent_format()) +
  ylab("Percent of total annual deaths") +
  theme_clean() +
  theme(
    axis.text.x  = element_text(size = 9, angle = 30, hjust = 1),
    axis.title.x = element_blank()
  )

# ── Assemble & save ───────────────────────────────────────────────────────────
fig_1 <- plot_grid(final_response, deaths_plot, ncol = 2, rel_widths = c(1.6, 1))

out_dir <- here::here("fig", "THL", "THL_comparison")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_path <- file.path(out_dir, "fig1_thl_comparison.pdf")
ggsave(plot = fig_1, filename = out_path, width = 12, height = 5, units = "in", dpi = 300)
message("Saved: ", out_path)
