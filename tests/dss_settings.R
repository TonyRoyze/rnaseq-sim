source("scripts/de_analysis/method_settings.R")

core <- make_method_settings("core")
quick <- make_method_settings("quick")
stopifnot("DSS" %in% vapply(core, `[[`, character(1), "method"))
stopifnot("DSS" %in% vapply(quick, `[[`, character(1), "method"))
stopifnot(!"dispersion_trend-trend-on" %in% vapply(
  core, `[[`, character(1), "setting_id"
))
stopifnot(!"dispersion_trend-trend-on" %in% vapply(
  quick, `[[`, character(1), "setting_id"
))
