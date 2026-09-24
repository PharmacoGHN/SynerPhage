# Figures: interaction map, single-pair effect ladder, posterior forest and
# density, and their conversion to interactive plotly widgets.

box::use(
  ggplot2,
  plotly[config, ggplotly],
  rlang[.data],
  scales[label_percent, squish],
)
box::use(
  app/logic/postprocessing[mcmc_node],
)

#' @export
cat_colours <- c(Synergy = "#1F7A5A", Additivity = "#9AA5AD",
                 Antagonism = "#B5443A", Inconclusive = "#D6B24B")

theme_phage <- function() {
  ggplot2$theme_minimal(base_size = 13) +
    ggplot2$theme(panel.grid.minor = ggplot2$element_blank(),
                  panel.grid.major = ggplot2$element_line(colour = "grey92"),
                  plot.title = ggplot2$element_text(face = "bold"),
                  plot.subtitle = ggplot2$element_text(colour = "grey35"),
                  strip.text = ggplot2$element_text(face = "bold"),
                  legend.key.height = ggplot2$unit(1, "lines"))
}

fmt_pct <- function(x) ifelse(is.na(x), "-", sprintf("%.1f%%", 100 * x))

# ggplotly keeps the ggplot geometry and adds hover, drag-to-zoom, pan and
# double-click reset. The `text` aesthetic carries the tooltip, which ggplot
# itself does not know about, hence the suppressed warning.
#' @export
as_interactive <- function(p, tooltip = "text") {
  suppressWarnings(
    ggplotly(p, tooltip = tooltip) |>
      config(displaylogo = FALSE,
             modeBarButtonsToRemove = c("lasso2d", "select2d",
                                        "hoverClosestCartesian", "hoverCompareCartesian"))
  )
}

# Strain order on the y axis: alphabetical, or ranked by the value on display
# so the strongest interactions gather at one end of the map.
#' @export
strain_levels <- function(df, value_var, sort_by) {
  if (identical(sort_by, "effect")) {
    m <- tapply(df[[value_var]], df$bacteria, mean, na.rm = TRUE)
    names(sort(m, na.last = FALSE))
  } else {
    rev(sort(unique(df$bacteria)))
  }
}

# Divergent fill centred on zero by construction, so the colour never suggests
# an interaction the number does not support. A probability is sequential
# instead: it has no meaningful midpoint at zero.
#' @export
matrix_plot <- function(df, fill_var, fill_lab, title, sort_by = "name") {
  df$bacteria <- factor(df$bacteria, levels = strain_levels(df, fill_var, sort_by))
  call_var <- if ("bayes_category" %in% names(df)) "bayes_category" else "point_category"
  df$text <- paste0(
    "<b>", df$bacteria, "  |  ", df$phage_combination, "</b>",
    "<br>observed: ", fmt_pct(df$observed),
    "<br>Bliss expected: ", fmt_pct(df$bliss_expected_point),
    "<br>", fill_lab, ": ", fmt_pct(df[[fill_var]]),
    "<br>call: ", df[[call_var]]
  )

  fill_scale <- if (identical(fill_var, "p_synergy")) {
    ggplot2$scale_fill_gradientn(colours = c("white", "#BFD8CF", cat_colours[["Synergy"]]),
                                 limits = c(0, 1), labels = label_percent(accuracy = 1),
                                 name = fill_lab)
  } else {
    lim <- max(0.05, ceiling(max(abs(df[[fill_var]]), na.rm = TRUE) * 20) / 20)
    ggplot2$scale_fill_gradient2(low = cat_colours[["Antagonism"]], mid = "white",
                                 high = cat_colours[["Synergy"]], midpoint = 0,
                                 limits = c(-lim, lim), oob = squish,
                                 labels = label_percent(accuracy = 1), name = fill_lab)
  }

  ggplot2$ggplot(df, ggplot2$aes(x = .data$phage_combination, y = .data$bacteria,
                                 text = .data$text)) +
    ggplot2$geom_point(ggplot2$aes(size = .data$observed, fill = .data[[fill_var]]),
                       shape = 21, colour = "grey35", stroke = 0.3) +
    fill_scale +
    ggplot2$scale_size_continuous(name = "Observed effect",
                                  limits = c(0, 1), range = c(1.5, 9),
                                  labels = label_percent(accuracy = 1)) +
    ggplot2$labs(x = NULL, y = NULL, title = title) +
    theme_phage() +
    ggplot2$theme(axis.text.x = ggplot2$element_text(angle = 45, hjust = 1))
}

# Single strain / single combination: show the ladder of effects instead of a
# one-bubble matrix.
#' @export
pair_plot <- function(summary, res) {
  comps <- strsplit(res$components, " + ", fixed = TRUE)[[1]]
  singles <- summary[summary$bacteria == res$bacteria & summary$treatment %in% comps, ]
  singles <- singles[match(comps, singles$treatment), ]

  df <- data.frame(
    label = c(singles$treatment, "Bliss expectation", "Observed combination"),
    value = c(singles$mean_effect, res$bliss_expected_point, res$observed),
    kind = c(rep("Phage alone", nrow(singles)), "Expected", "Observed"),
    stringsAsFactors = FALSE
  )
  df$label <- factor(df$label, levels = rev(df$label))
  df$text <- paste0("<b>", df$label, "</b><br>", df$kind, ": ", fmt_pct(df$value))

  ggplot2$ggplot(df, ggplot2$aes(x = .data$value, y = .data$label, colour = .data$kind,
                                 text = .data$text)) +
    ggplot2$geom_segment(ggplot2$aes(x = 0, xend = .data$value, yend = .data$label),
                         linewidth = 1.1) +
    ggplot2$geom_point(size = 5) +
    # Nudged above the line: plotly ignores hjust/vjust, but a position is
    # applied by ggplot before the conversion, so it holds in both renderers.
    ggplot2$geom_text(ggplot2$aes(label = label_percent(accuracy = 0.1)(.data$value)),
                      position = ggplot2$position_nudge(y = 0.3), size = 4,
                      show.legend = FALSE) +
    # Headroom for the value labels, without drawing a tick beyond 100%.
    ggplot2$scale_x_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25),
                               labels = label_percent(accuracy = 1)) +
    ggplot2$scale_colour_manual(values = c("Phage alone" = "#6C8EA0", Expected = "#9AA5AD",
                                           Observed = "#12626B"), name = NULL) +
    ggplot2$labs(x = "Effect (fraction of growth suppressed)", y = NULL,
                 title = sprintf("%s : %s", res$bacteria, res$components),
                 subtitle = sprintf("observed - expected = %+.1f points (ROPE +/- %.1f) : %s",
                                    100 * res$delta_observed, 100 * res$rope,
                                    res$point_category)) +
    theme_phage()
}

#' @export
posterior_forest <- function(post, rope, ncol = 3, sort_by = "name") {
  post$bacteria <- factor(post$bacteria,
                          levels = strain_levels(post, "synergy_delta", sort_by))
  post$text <- paste0(
    "<b>", post$bacteria, "  |  ", post$phage_combination, "</b>",
    "<br>synergy effect: ", fmt_pct(post$synergy_delta),
    "<br>95% CrI: ", fmt_pct(post$synergy_delta_lower), " to ",
    fmt_pct(post$synergy_delta_upper),
    "<br>P(synergy): ", fmt_pct(post$p_synergy),
    "<br>call: ", post$bayes_category
  )

  ggplot2$ggplot(post, ggplot2$aes(x = .data$synergy_delta, y = .data$bacteria,
                                   colour = .data$bayes_category, text = .data$text)) +
    ggplot2$annotate("rect", xmin = -rope, xmax = rope, ymin = -Inf, ymax = Inf,
                     fill = "grey80", alpha = 0.45) +
    ggplot2$geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    ggplot2$geom_linerange(ggplot2$aes(xmin = .data$synergy_delta_lower,
                                       xmax = .data$synergy_delta_upper)) +
    ggplot2$geom_point(size = 2.4) +
    ggplot2$facet_wrap(~phage_combination, ncol = ncol) +
    ggplot2$scale_colour_manual(values = cat_colours, name = "Call", drop = FALSE) +
    ggplot2$scale_x_continuous(labels = label_percent(accuracy = 1)) +
    ggplot2$labs(x = "Synergy effect (posterior median and 95% CrI)", y = NULL,
                 title = "Posterior synergy effect",
                 subtitle = "Shaded band: region of practical equivalence") +
    theme_phage()
}

#' @export
posterior_density <- function(draws, combo, strain_index, strain, rope) {
  x <- mcmc_node(as.matrix(draws[[combo]]), "synergy_effect")[, strain_index]
  ggplot2$ggplot(data.frame(x = x), ggplot2$aes(x = .data$x)) +
    ggplot2$annotate("rect", xmin = -rope, xmax = rope, ymin = -Inf, ymax = Inf,
                     fill = "grey80", alpha = 0.45) +
    ggplot2$geom_density(fill = "#12626B", colour = "#12626B", alpha = 0.25, linewidth = 0.8) +
    ggplot2$geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    ggplot2$scale_x_continuous(labels = label_percent(accuracy = 1)) +
    ggplot2$labs(x = "Synergy effect", y = "Posterior density",
                 title = sprintf("%s : %s", strain, combo),
                 subtitle = "Shaded band: region of practical equivalence") +
    theme_phage()
}
